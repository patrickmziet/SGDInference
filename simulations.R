rm(list=ls())
source("implicit_sgd.R")
source("exact_inference.R")

get_coverage <- function(ci, theta) {
    ## Determines for each theta whether CI contains theta
    cover = c()
    for(i in 1:nrow(ci)) {
        cover = c(cover, theta[i] <= ci[i, 2] & theta[i] >= ci[i, 1])
    }
    return(as.logical(cover))
}


print_cistat <- function(procedure, ci_stat, j) {
    print(sprintf("%s mean coverage: %.3f -- CI length(x 10-2): %.3f", procedure,
                  mean(ci_stat$cover/j), 100 * ci_stat$len / j))
}

print_gammastar <- function(gamma_star, sgd_stat, j) {
    print(sprintf("Gamma star  = %.3f/ sgd mean coverage(-problematic): %.3f", gamma_star, 
                  mean(sgd_stat$cover[sgd_stat$cover/j > 0.5])/j))
}

plot_coverage <- function(sgd_stat, glm_stat, j=NA, title) {
    ## if nreps param given, must divide coverage & estimate stats by it 
    if(!is.na(j)) {
        sgd_stat$cover = sgd_stat$cover / j
        glm_stat$cover = glm_stat$cover / j
        sgd_stat$len   = sgd_stat$len / j
        glm_stat$len   = glm_stat$len / j
    } 
    
    ## compute plot ylim bounds
    y.upper = max(c(sgd_stat$cover, glm_stat$cover))
    y.lower = min(c(sgd_stat$cover, glm_stat$cover))
    delta   = 0.5 * (y.upper - y.lower)
    
    ## plot
    plot(glm_stat$cover, col='red', pch=18,
         ylab="coverage", xlab="theta_star idx", cex.main=0.75,
         ylim=c(0, 1), main=title)
    abline(h=0.95, col='blue')
    abline(h=mean(glm_stat$cover), lty=3, col='red')
    points(sgd_stat$cover, pch=18)
    abline(h=mean(sgd_stat$cover), lty=3)
    legend(x="bottomright", legend=c("sgd","glm"), col=c("black","red"), pch=18)
}

plot_params <- function(par1, par2, par_true) {
    ## compute plot ylim bounds
    ## plot
    plot(par_true, par_true, type="l", lty=3, cex.main=0.75,col="blue")
    
    points(par_true, par1$est, pch=4, col="black")
    points(par_true, par2$est, pch=3, col="red")
    legend(x="bottomright", legend=c(par1$name, par2$name), col=c("black","red"), pch=c(4,3))
}

plot_summary <- function(gamma_star, sgd_stat, glm_stat, j) {
    ## output strings
    gamma_star.str = sprintf("Gamma* = %.3f", gamma_star)
    ## if nreps param given, must divide coverage & estimate stats by it 
    if(!is.na(j)) {
        sgd_stat$cover = sgd_stat$cover / j
        glm_stat$cover = glm_stat$cover / j
        sgd_stat$len   = sgd_stat$len / j
        glm_stat$len   = glm_stat$len / j
        sgd_stat$l2    = sgd_stat$l2 / j
        glm_stat$l2    = glm_stat$l2 / j
    } 
    sgd_ci.str = sprintf("sgd mean coverage: %.3f \n-- CI length(x 10-2): %.3f\n",
                         mean(sgd_stat$cover), 100 * sgd_stat$len)
    glm_ci.str = sprintf("glm mean coverage: %.3f \n-- CI length(x 10-2): %.3f\n",
                         mean(glm_stat$cover), 100 * glm_stat$len)
    sgd_est.str  = sprintf("sgd mean l2/p(x 10-2): %.3f", 100 * sgd_stat$l2)
    glm_est.str = sprintf("glm mean l2/p(x 10-2): %.3f",  100 * glm_stat$l2)
    nreps.str    = sprintf("num experiments = %d", j)
    
    ## ## plot empty plot
    ## plot(10, ylim=c(0,1), 
    ##      yaxt="n", xaxt="n")
    ## legend(x="bottomleft", bty="n", cex=0.75,
    ##        legend=c(gamma_star.str,
    ##                 sgd_ci.str,
    ##                 glm_ci.str,
    ##                 sgd_est.str,
    ##                 glm_est.str,
    ##                 nreps.str))
}

#### Hyperameters
#### N, p, theta_star
#### - gamma_star
#### - npass (exact_inf) + npass (at gamma_star selection)
#### - Var(theta_n - theta_star)  : what is theta_star?
## Good examples:
## Poisson, toeplitz, classic

check_bias <- function(p, N) {
    
    sample_data <- function(sigma_sqrt=NA) {
        fam <<- "binomial"  ## binomial for logistic
        gen_data(p=p, N=N, sigma_noise=1,
                 model_name = fam, 
                 sigma_x = "id", true_param="usc",
                 sigma_sqrt = sigma_sqrt)
    }
    
    thetas = matrix(0, nrow=0, ncol=p)
    d = sample_data()
    out = exact_inf(d)
    se = head(out$se_sgd, 1)
    for(j in 1:10) {
        thetas = matrix(0, nrow=0, ncol=p)
        for(i in 1:1e3) {
            theta = as.numeric(isgd_avg(d, npass=j))
            thetas = rbind(thetas, theta)
        }   
        se_empirical = apply(thetas, 2, sd)                                                              
        
    }
    
    for(i in 1:100) {
        d = sample_data()
        theta = initialize_sgd(d)
        thetas = rbind(thetas, theta)
    }
    avg = colMeans(thetas)
    plot(d$theta_star, d$theta_star, type="l", lty=3, col="red")
    points(d$theta_star, avg, pch=20)
    
    ci = matrix(0, nrow=p, ncol=2)
    ci[, 1] = lower.ci
    ci[, 2] = upper.ci
}


diagnostics <- function() {
    load("last_run.rda")
    d = last_run$data
    p = ncol(d$X)
    N = length(d$Y)
    sgd_control = last_run$sgd_control
    ## init_control = last_run$init_control
    theta0 = as.numeric(sgd_control$theta0)
    gamma_star = sgd_control$gamma
    
    lam = eigen(fisher(d))$values
    print(sprintf("lmin=%.2f, g>%.2f, lmax=%.2f", min(lam), 1/min(lam), max(lam)))
    theta0[1]
    d$theta_star[1]
    
    ## Improve?
    print("Init control")
    last_run$init_control$npass = as.integer(last_run$init_control$npass)
    print(sprintf("Init control: npass=%d, lr=%.2f",  
                  last_run$init_control$npass, last_run$init_control$lr))
    
    theta0 = as.numeric(aisgd(d, gamma=20, npass=100))
    theta0[1]
    d$theta_star[1]
    
    example_sim(p=100, sigma_x="ill_cond", 
                init_control=list(theta0=rep(0, p), lr=20, npass=100))
    
    par(mfrow=c(1, 2))
    theta0 = naive_init(d)
    gamma_star=20
    theta_hat = isgd_0(theta0, data=d, gamma=gamma_star, npass=10, use_permutation = F)
    theta_hat[1]  + c(-1, 1) * 2 * sqrt(gamma_star/N)
    d$theta_star[1]
    
    plot_params(list(est=theta0, name="theta0"), list(est=d$theta_star, name="true"), d$theta_star)
    plot_params(list(est=theta_hat, name="isgd"), list(est=d$theta_star, name="true"), d$theta_star)
    
    gamma_star=10
    inf = exact_inf(d, sgd_control=list(gamma=gamma_star, theta0=theta0))
    ## CI length
    print("CI for param 1")
    print(inf$ci[1, ])
    print("True value")
    d$theta_star[1]
    print("CI for param 2")
    print(inf$ci[2, ])
    print("True value")
    d$theta_star[2]
    get_coverage(inf$ci, d$theta_star)
}

breakdown <- function() {
    p = 50
    N = 1e4
    sample_data <- function(sigma_sqrt=NA) {
        fam <<- "binomial"  ## binomial for logistic
        gen_data(p=p, N=N, sigma_noise=1,
                 model_name = fam, 
                 sigma_x = "ill_cond", true_param="classic",
                 sigma_sqrt = sigma_sqrt)
    }

    thetas = matrix(0, nrow=0, ncol=p)
    d = sample_data()
    F = fisher(d)
    lam = eigen(F)$values
    lmin = min(lam)
    print(sprintf("Minimum eigenvalue %f", min(lam)))
#### exact_inf breakdown
    ## 1. gamma_star 
    bounds = Lmin_bounds(d)
    print(sprintf("Bounds for gamma: lower = %.2f, upper = %.2f", 
                  bounds$lower, bounds$upper))
    gamma_seq = seq(.5 / bounds$upper, 2 / bounds$lower, length.out=kGammas)
    MVc = par_misgd(d, gamma_seq, init_ctl=list(npass=50, lr=0.3))
    x = MVc$gamma_x
    y = MVc$trace_y
    ## Study SGD Iterates.
    make_plots <- function(D, title) {
        j = which.max(colMeans(D))
        units = c(j, setdiff(1:p, j))
        for(i in units) {
            lo = loess(D[, i] ~ x)
            vi = predict(lo)
            if(i==j) { plot(x, vi, type="l", col=i, main=title) }
            lines(x, vi, col=i)
            m = which.max(vi)
            text(x[m], vi[m] + 100, labels=c(i), col=i)
        }
    }
    par(mfrow=c(2,2))
    make_plots(MVc$Bias, "bias")
    make_plots(MVc$Var, "Variance")
    plot(x, y, type="l", col="magenta")
    lines(x, x* p/2, lty=3)
    lx = function(lam_est) {
        sapply(x, function(g) sum(g^2 * lam_est / (2 * g * lam_est - 1)))
    }
    ## obj = function(par) mean((lx(exp(par)) - y)^2)
    ## out = optim(par=rep(1e-5, p), fn=obj, method="BFGS")
    ## print(exp(out$par))
    lines(x, lx(lam), col="magenta", lty=2)
    lines(x, x*p/2 + .25 * sum(1/lam), lty=3, col="red", lwd=2)
    
    par(mfrow=c(1, 2))
    best_gamma(x, y, lmin, p)
    best_gamma_v2(x, y, lmin, p)
    
#### Condition number
    C = max(lam) / min(lam) ## assume known
    warning("assuming condition number is known")
    
    gstar = best_gamma_v2(x, y, lmin, p)
    print(sprintf("Minimum gamma = %.5f", gstar$value))
    ## gamma_star = calculate_gammaStar(d)
    ## 2. init 
    theta0 = MVc$theta0
    ## 3. theta_hat estimate
    
    theta_hat = isgd0(theta0, d, 
                      gamma=gstar$value, 
                      use_permutation = T)
    psi = get_model_psi(d, theta_hat)
    ## V = psi * gstar * rep(1, p)
    V = psi * gstar$value * gstar$inflation * rep(1, p)
    
    V_star = V / N
    upper.ci = theta_hat + 1.96 * sqrt(V_star)
    lower.ci = theta_hat - 1.96 * sqrt(V_star)
    
    plot(1:p, lower.ci, col="red", type="l")
    lines(1:p, upper.ci, col="blue", type="l")
    points(1:p, d$theta_star, pch=20)

    ci = matrix(0, nrow=p, ncol=2)
    ci[, 1] = lower.ci
    ci[, 2] = upper.ci
    get_coverage(ci, d$theta_star)
    mean(get_coverage(ci, d$theta_star))
}

conditioning_exp <- function() {
    p = 9
    N = 1e3
    X = rmvnorm(N, mean=rep(0, p), sigma=diag(c(rep(1, p-1), 1e6)))
    theta_star = (-1)^seq(1, p) * exp(-.5 * seq(1, p))
    theta_star
    data = gen_data_star(theta_star, X, model_name = "gaussian")
    theta0 = rep(0, p)
    out = isgd_R(data, theta0=theta0, gamma = 10)
    par(mfrow=c(sqrt(p), sqrt(p)))
    for(j in 1:p) {
        plot(out[j, ], type="l", col="blue")
        abline(h=theta_star[j], lty=3, col="red", lwd=2)
    }
    out = isgd_RC(data, gamma=10)
    for(j in 1:p) {
        plot(out[j, ], type="l", col="blue")
        abline(h=theta_star[j], lty=3, col="red", lwd=2)
    }
}


example_sim <- function(p=10, N=1e4, nreps=1e2, 
                        model="gaussian", sigma_x="id", true_param="classic",
                        sgd_control=list(), init_control=list(),
                        verbose=TRUE) {
    print(sprintf("Model = %s, sigma_x = %s, params = %s", model, sigma_x, true_param))
    ## output statistics
    sgd_stat = list(cover=rep(0, p), len=0, l2=0)
    glm_stat = list(cover=rep(0,p), len=0, l2=0) 
    
    ## sample data function
    sample_data <- function(sigma_sqrt=NA) {
        gen_data(p=p, N=N, sigma_noise=1,
                 model_name = model, 
                 sigma_x = sigma_x, true_param=true_param,
                 sigma_sqrt = sigma_sqrt)
    }
    
    ## first compute sigma sqrt
    data = sample_data()
    S.sqrt = data$Sigma_X_sqrt
    
    ## create split screen
    close.screen(all.screens = TRUE)
    par(bg="white")
    split.screen(c(2,2))
    screen(1)
    par(mar=rep(2, 4))
    
    ## calculate or load gamma star
    if(length(sgd_control)==0) {
        if(length(init_control)==0) {
            init_control = default_init(p)
        }
        sgd_control = calculate_gammaStar(data, init_control = init_control)
        print(sprintf("> gamma* = %.2f", sgd_control$gamma))
    }

    ## ci stat update
    stat_update <- function(ci_stat, out, dat) {
        ci_stat$cover = ci_stat$cover + get_coverage(out$ci, dat$theta_star)
        ci_stat$len   = ci_stat$len   + mean(out$ci[ ,2] - out$ci[ ,1])
        ci_stat$l2    = ci_stat$l2    + L2.p(dat$theta_star, out$est)
        return(ci_stat)
    }
    ## timing
    t0 = proc.time()[3]
    ## Interesting method. Not for now...
    for(j in 1:nreps) {
        data = sample_data(sigma_sqrt = S.sqrt)
        ##
        sgd_out = exact_inf(data, sgd_control)
        glm_out = glm_inf(data)
        
        ## update ci stat
        sgd_stat = stat_update(sgd_stat, sgd_out, data)
        glm_stat = stat_update(glm_stat, glm_out, data)
        
        if(verbose) {
            t1 = proc.time()[3]
            if(t1 - t0 > 5) {
                ## Save instance before moving on...  
                last_run = list(data=data, sgd_control=sgd_control, init_control=init_control)
                ## ==> Done saving.
                save(last_run, file="last_run.rda")
                
                print(sprintf("j = %d out of %d", j, nreps))
                print(sgd_stat$cover/j)
                print_cistat("sgd", sgd_stat, j)
                print_gammastar(sgd_control$gamma, sgd_stat, j)
                print_cistat("glm", glm_stat, j)
                t0 = t1
                
                screen(2)
                par(mar=rep(2, 4))
                pstr = ifelse(p==500, "5e2", "1e2")
                title = sprintf("%s, %s, %s, 1e%d", sigma_x, true_param, 
                                pstr, log(N, base=10))
                plot_coverage(sgd_stat=sgd_stat, glm_stat=glm_stat, j, title)
                
                screen(3)
                par(mar=rep(2, 4))
                plot_params(par1=list(est=sgd_out$est, name="sgd"), 
                            par2=list(est=glm_out$est, name="glm"),
                            par_true=data$theta_star)
                
                screen(4)
                par(mar=rep(2, 4))
                plot_params(par1=list(est=sgd_control$theta0, name="theta0"), 
                            par2=list(est=glm_out$est, name="glm"),
                            par_true=data$theta_star)
            }
        }
    }
    
    close.screen(all.screens=TRUE)
    
    ## output coverage data
    return(list(sgd_cover=sgd_stat$cover/nreps, sgd_cover_len=sgd_stat$len/nreps, sgd_l2=sgd_stat$l2/nreps,
                glm_cover=glm_stat$cover/nreps, glm_cover_len=glm_stat$len/nreps, glm_l2=glm_stat$l2/nreps,
                nreps=nreps, gamma_star=sgd_control$gamma))
}

parallel_sim <- function(p=10, N=1e4, nreps=1e2,
                         model="gaussian", sigma_x="id", rho=0.15, true_param="classic",
                         sgd_control=list(), init_control=list(),
                         print_diagnostic=TRUE) {
    ##:ess-bp-start::conditional@:##
    browser(expr={TRUE})##:ess-bp-end:##
    print(sprintf("Model = %s, sigma_x = %s, rho = %f, params = %s, p = %d, N = %d, nreps = %d", 
                  model, sigma_x, rho, true_param, p, N, nreps)) 
    ## output statistics
    sgd_stat = list(cover=rep(0, p), len=0, l2=0)
    glm_stat = list(cover=rep(0,p), len=0, l2=0) 

    ## sample data function
    sample_data <- function(sigma_sqrt=NA) {
        gen_data(p=p, N=N, sigma_noise=1,
                 model_name = model, 
                 sigma_x = sigma_x, rho = rho, 
                 true_param=true_param, sigma_sqrt = sigma_sqrt)
    }
    
    ## first compute sigma sqrt
    data = sample_data()
    S.sqrt = data$Sigma_X_sqrt

    if (print_diagnostic) {
        ## create split screen
        close.screen(all.screens = TRUE)
        par(bg="white")
        split.screen(c(2,2))
        screen(1)
        par(mar=rep(2,4))
    }

    ## calculate or load gamma star
    if(length(sgd_control)==0) {
        ## only rewrite inputs ontop of default init_control
        init_control_def = default_init(p)
        for(n in names(init_control)) init_control_def[[n]] = init_control[[n]]
        init_control = init_control_def
        
        sgd_control = calculate_gammaStar_wrapper(data, init_control=init_control)
        ##print(sprintf("> gamma* = %.2f", sgd_control$gamma))
    } else if (!is.null(sgd_control$use_gamma_half)) {
        ## only rewrite inputs ontop of default init_control
        init_control_def = default_init(p)
        for(n in names(init_control)) init_control_def[[n]] = init_control[[n]]
        init_control = init_control_def
        
        sgd_control = calculate_gammaStar_wrapper(data, init_control=init_control, which.method)
        sgd_control$use_gamma_half = TRUE
    }

    g <- function(n, u, p=NULL, type="naive") {
        ## Threshold function
        switch(type,
               "naive"     = n^u,
               "enhanced1" = n^u + exp(1 / (1 - 2 * u) - 1) - 2,
               "enhanced2" = n^u - (1 - 4 * u) / (1 - 2 * u),
               "dBIC"      = sqrt((n^(1/n) - 1) * (n - p))
               stop("unknown type"))
    }
    
    ## ci stat update
    stat_update <- function(out, dat) {
        ci_stat       = list(cover=rep(0, p), len=0, l2=0,
                             J=rep(FALSE, p), Jhat=rep(FALSE, p),
                             tp=0, fp=0, tn=0, fn=0, prec=0, rec=0, f1=0)
        ci_stat$cover = get_coverage(out$ci, dat$theta_star) # is true theta in CIs
        ci_stat$len   = mean(out$ci[ ,2] - out$ci[ ,1]) # average CI length
        ci_stat$l2    = L2.p(dat$theta_star, out$est)   # L2 norm ||truth - estimate||_2
        ci_stat$J     = dat$theta_star != 0             # True model
        ci_stat$Jhat  = out$pivots > g(n = N, u = 1/3)  # Estimated model
        ci_stat$tp    = sum(ci_stat$J & ci_stat$Jhat)   # True positives
        ci_stat$fp    = sum(!ci_stat$J & ci_stat$Jhat)  # False positives
        ci_stat$tn    = sum(!ci_stat$J & !ci_stat$Jhat) # True negatives
        ci_stat$fn    = sum(ci_stat$J & !ci_stat$Jhat)  # False negatives
        ci_stat$prec  = ci_stat$tp / sum(ci_stat$Jhat)  # Precision
        ci_stat$rec   = ci_stat$tp / sum(ci_stat$J)     # Recall
        ci_stat$f1    = 2 / (1 / ci_stat$prec + 1 / ci_stat$rec) # F1 score
        return(ci_stat)
    }
    
    ## compute one confidence interval
    single_ci <- function(i) {
        i = 1
        data = sample_data(sigma_sqrt = S.sqrt)
        sgd_out = exact_inf(data, sgd_control)
        glm_out = glm_inf(data)
        return(list(sgd=stat_update(sgd_out, data), 
                    glm=stat_update(glm_out, data)))
    }
    
    ## convert mclapply output to confidence interval stats	
    convert_par <- function(mult_ci) {
        sgd_stat = list(cover=rep(0, p), len=0, l2=0, tp=0, fp=0, tn=0, fn=0, prec=0, f1=0)
        glm_stat = list(cover=rep(0, p), len=0, l2=0, tp=0, fp=0, tn=0, fn=0, prec=0, f1=0)
        for(i in 1:nreps) {
            sgd_stat$cover = sgd_stat$cover + mult_ci[[i]]$sgd$cover
            glm_stat$cover = glm_stat$cover + mult_ci[[i]]$glm$cover
            sgd_stat$len   = sgd_stat$len   + mult_ci[[i]]$sgd$len
            glm_stat$len   = glm_stat$len   + mult_ci[[i]]$glm$len
            sgd_stat$l2    = sgd_stat$l2    + mult_ci[[i]]$sgd$l2
            glm_stat$l2    = glm_stat$l2    + mult_ci[[i]]$glm$l2
            sgd_stat$tp    = sgd_stat$tp    + mult_ci[[i]]$sgd$tp
            glm_stat$tp    = glm_stat$tp    + mult_ci[[i]]$glm$tp
            sgd_stat$fp    = sgd_stat$fp    + mult_ci[[i]]$sgd$fp
            glm_stat$fp    = glm_stat$fp    + mult_ci[[i]]$glm$fp
            sgd_stat$tn    = sgd_stat$tn    + mult_ci[[i]]$sgd$tn
            glm_stat$tn    = glm_stat$tn    + mult_ci[[i]]$glm$tn
            sgd_stat$fn    = sgd_stat$fn    + mult_ci[[i]]$sgd$fn
            glm_stat$fn    = glm_stat$fn    + mult_ci[[i]]$glm$fn
            sgd_stat$prec  = sgd_stat$prec  + mult_ci[[i]]$sgd$prec
            glm_stat$prec  = glm_stat$prec  + mult_ci[[i]]$glm$prec
            sgd_stat$f1    = sgd_stat$f1    + mult_ci[[i]]$sgd$f1
            glm_stat$f1    = glm_stat$f1    + mult_ci[[i]]$glm$f1
        }
        return(list(sgd_stat=sgd_stat,
                    glm_stat=glm_stat))
    }
    
    ## run parallel
    print("> Computing confidence intervals in parallel..")
    mult_ci  = mclapply(1:nreps, function(i) single_ci(i))
    
    ## convert to sgd_stat, glm_stat
    out_ci   = convert_par(mult_ci)
    sgd_stat = out_ci$sgd_stat
    glm_stat = out_ci$glm_stat 
    sgd_out = exact_inf(data, sgd_control)
    glm_out = glm_inf(data)

    if (print_diagnostic) {
        ##print summary output
        cat('\n')
        print_cistat("sgd", sgd_stat, nreps)
        cat('\n')
        print_gammastar(sgd_control$gamma, sgd_stat, nreps)
        cat('\n')
        print_cistat("glm", glm_stat, nreps)
        cat('\n')

        ##plot output
        screen(2)                                                                                   
        par(mar=rep(2, 4)) 
        ##pstr = ifelse(p==500, "5e2", "1e2")
        pstr = formatC(p, format='e', digits=0)
        title = sprintf("%s, %s, %s, 1e%d", sigma_x, true_param, 
                        pstr, log(N, base=10))
        plot_coverage(sgd_stat=sgd_stat, glm_stat=glm_stat, nreps, title)
        
        screen(3)
        par(mar=rep(2, 4)) 
        plot_params(par1=list(est=sgd_out$est, name="sgd"), 
                    par2=list(est=glm_out$est, name="glm"),
                    par_true=data$theta_star)
        
        screen(4)
        par(mar=rep(2, 4)) 
        plot_params(par1=list(est=sgd_control$theta0, name="theta0"), 
                    par2=list(est=glm_out$est, name="glm"),
                    par_true=data$theta_star)
        
        close.screen(all=TRUE)
    }
    
    ## output coverage data
    return(list(sgd_cover=sgd_stat$cover/nreps, sgd_cover_avg=mean(sgd_stat$cover/nreps),
                sgd_cover_len=sgd_stat$len/nreps, sgd_l2=sgd_stat$l2/nreps,
                sgd_cover_under50=Reduce('|', sgd_stat$cover/nreps < 0.5),
                glm_cover=glm_stat$cover/nreps, glm_cover_avg=mean(glm_stat$cover/nreps), 
                glm_cover_len=glm_stat$len/nreps, glm_l2=glm_stat$l2/nreps,
                glm_cover_under50=Reduce('|', glm_stat$cover/nreps < 0.5),
                nreps=nreps, gamma_star=sgd_control$gamma))
}

#### Code to be run in R markdown scripts to generate SGD Inference reports
run_report <- function(model, num_experiments=20,
                       init_control=list(), sgd_control=list()) {
    ## read in parameter values and numbering
    report.params = read.csv("report_params.csv", as.is=TRUE)
    model_number = ifelse(model=="gaussian", 3, 4)
    
    ## print section header
    cat(sprintf("## %d %s\n", model_number, model))
    
    for(i in 1:nrow(report.params)) {
        ## parameters
        tracking_number = report.params$tracking_number[i]
        var             = report.params$var[i]
        true_param      = report.params$true_param[i]
        p               = report.params$p[i]
        N               = report.params$N[i]
        
        ## section title
        cat(sprintf('#### %d.%d (Sigma, theta_star, p, N) = (%s, %s, %d, 1e%d)\n',
                    model_number, tracking_number, var, true_param, p, log(N, base=10)))
        
        ## run confidence interval experiment
        out = parallel_sim(p=p, N=N, nreps=num_experiments, 
                           model=model, sigma_x=var, true_param=true_param,
                           sgd_control=sgd_control, init_control=init_control,
                           print_diagnostic=TRUE)
        
        cat('\n\n')
    }
}

colnames.exp = c("model", "var", "true_param", "p", "N",
                 "sgd_cover_avg", "sgd_cover_len", "sgd_l2", "sgd_cover_under50",
                 "glm_cover_avg", "glm_cover_len", "glm_l2", "glm_cover_under50",
                 "nreps", "gamma_star")

write_data <- function(dat, out.fpath) {
    dat = matrix(dat, nrow=1, dimnames=list(c(), colnames.exp))
    write.table(dat, file=out.fpath, sep=";", append=TRUE,
                row.names=FALSE, col.names=!file.exists(out.fpath))
}

#### To run multiple confidence interval experiments, prints data to csv out.
#### No reports generated, just the numbers.
run_experiments <- function(model, num_experiments=20,
                            init_control=list(), sgd_control=list(),
                            params.fpath, out.fpath) {
    ## read in parameter values
    exp.params = read.csv(params.fpath, as.is=TRUE)

    for(i in 1:nrow(exp.params)) {
        ## parameters
        var             = exp.params$var[i]
        true_param      = exp.params$true_param[i]
        p               = exp.params$p[i]
        N               = exp.params$N[i]

        ## run confidence interval experiment
        out = parallel_sim(p=p, N=N, nreps=num_experiments,
                           model=model, sigma_x=var, true_param=true_param,
                           sgd_control=sgd_control, init_control=init_control,
                           print_diagnostic=FALSE)

        ## pull sgd_cover and glm_cover from out
        out$sgd_cover <- NULL
        out$glm_cover <- NULL

        ## append output
        write_data(c(model, var, true_param, p, N, out), out.fpath)
    }
}

