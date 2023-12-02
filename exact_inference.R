## Code for generating exact confidence intervals for glm's using ISGD
## Chooses "best gamma", ie smallest gamma to give exact CI's
## Jerry Chee, Panos Toulis

##to test, uncomment below and run
##dat = gen_data()
##SGD_CI = exact_inf(dat)
##GLM_CI = confint(glm(dat$Y ~ dat$X))[-1,]

##===================
##MVc  = multiple_isgd(data, seq(0.01, 10, length.out=25))
kReps = 100
kGammas = 30
## runs R ISGD chains to calc c_t values, non serial version

default_init <- function(p) {
    return(list(theta0=rep(0, p), npass=as.integer(10 * sqrt(p)), lr=1,
                gamma.method="lmin", bound.method="lmin"))
}

get_model_psi <- function(data, theta) {
    ## Function to compute the dispersion parameter (?)
    psi = NA
    if (data$model == "gaussian") {
        psi = as.numeric(var(data$Y - data$X %*% matrix(theta, ncol=1)))
    } else if (data$model == "binomial") {
        psi = 1
    } else if (data$model == "poisson") {
        psi = 1
    } else {
        stop("psi not specified")
    }
    return(psi)
}

naive_init <- function(data) {
    p = ncol(data$X)
    theta0 = c()
    for(j in 1:p) {
        fit = glm(data$Y ~ data$X[, j] + 0, family=data$model)
        theta0 = c(theta0, coef(fit)[1])
    }
    return(as.numeric(theta0))
}

initialize_sgd <- function(data, init_control) {
    ## if init_control=list(), does averaged isgd init,
    ## otherwise if pass arg use_true or use_naive, will perform that init
    
    ## TODO: Think about better initialization?
    ## TODO: Averaged implicit SGD?
    ## TODO: How about we choose gamma, npass based on 
    ## removing bias to some pre-defined theta_star?
    ##   
    ## binomial, equicor pretty bad with isgd, isgd_const, isgd_avg, isgd_adagrad
    ## return(data$theta_star)
    ## return(isgd(data, gamma=5, npass=100))
    if(!is.null(init_control$use_true)) {
        warning("Using true theta_star values")
        return(data$theta_star)
    }
    if(!is.null(init_control$use_naive)) {
        print(">> Using naive initialization.")
        return(naive_init(data))
    }
    p = ncol(data$X)
    with(init_control, aisgd_0(theta0, data=data, gamma=lr, npass=npass))
}

best_gamma_v2 <- function(gamma_x, trace_y, lmin_low, p) {
    k = length(gamma_x)
    stopifnot(length(gamma_x)==length(trace_y))
    M = matrix(0, nrow=0, ncol=4)
    colnames(M) <- c("i", "j", "good.slope", "is.normal")
    for(i in seq(1, k-1)) {
        for(j in seq(i+1, k)) {
            if(j-i > 4) {
                x1 = gamma_x[i:j]
                y1 = trace_y[i:j]
                fit = lm(y1 ~ x1)
                ci = as.numeric(confint(fit)[2, ])
                is.hit = (p/2 >= ci[1] & p/2 <= ci[2])
                is.normal = shapiro.test(fit$res)$p.value > 0.1
                M = rbind(M, c(i, j, is.hit, is.normal))
            }
        }
    }
    M = as.data.frame(M)
    M0 = subset(M, good.slope==1 & is.normal==1)
    len = apply(M, 1, function(row) row[2] - row[1])
    M$len = len
    hasWarning = FALSE
    if(max(len) < 5 | nrow(M0) == 0) {
        hasWarning = TRUE
        print("[WARNING]: There is no good learning rate for this problem.")
        if(mean(diff(diff(trace_y)) < 0) > 0.6) {
            print("[..] Trace plot has concave shape. Increase #samples?")
        } else {
            print("[WARNING] Please check bias/variance plots for more info.")
        }
    }
    ## M0 = subset(M, good.slope==1 & is.normal==1)
    ## ids = as.numeric(numeric(rownames(M0)))
    m = which.max(M$len * M$good.slope * M$is.normal / sqrt(M$i))
    from = M$i[m]
    to = M$j[m]
    
    plot(gamma_x, trace_y, col="blue", type="l", ylim=c(0, max(trace_y)))
    points(gamma_x[from:to], rep(0.01, to - from + 1), 
           lwd=3, col="magenta", pch=15)
    if(hasWarning) {
        text(gamma_x[to + 5], trace_y[from], "Bad fit.", col="red")
        return(max(gamma_x))
    } else {
        lines(gamma_x[from:to], p * gamma_x[from:to]/2, lwd=2, 
              col="magenta")
        ## return
        med = as.integer((from + to)/2)
        gstar = gamma_x[med]
        abline(v=gstar, lty=3, col="red")
        return(gstar)
    }
}

par_misgd <- function(data, gamma_ls, init_control, R=kReps) {
    N = nrow(data$X)
    p = ncol(data$X)
    ng = length(gamma_ls)
    ## TODO: Should we run init in every thread?
    theta0 = initialize_sgd(data, init_control=init_control)
    
    print(sprintf(">> Running Parallel Boostraps..%d gammas x %d replications", ng, R))
    
    
    single_isgd <- function(gamma_t) {
        out = mclapply(1:R, function(r) {
            isgd_0(theta0, data, gamma=gamma_t, use_permutation = T)
        })
        thetas = matrix(0, nrow=R, ncol=p)
                                        # transform to matrix
        for(i in 1:R) {
            thetas[i,] = out[[i]]
        }
                                        # thetas = R x p matrix of SGD estimates.
                                        # Var(θn)
        tr_t = sum(N * diag(cov(thetas)))
                                        # (E θn - θ*)^2
        m_t = colMeans(thetas)
        c1 = N * (m_t  - data$theta_star)^2
        th_star = matrix(data$theta_star, nrow=R, ncol=p, byrow=T)
                                        # E(θn - θ*)^2
        c2 = N * diag(t(thetas - th_star) %*% (thetas - th_star) / R)  # pxp
        
        return(list(trace=tr_t, bias=c1, variance=c2))
    }
    
    out = mclapply(gamma_ls, function(gamma) {
        single_isgd(gamma)
    })
    
    Tr = c()
    B = matrix(0, nrow=0, ncol=p) # bias
    V = matrix(0, nrow=0, ncol=p) # variance
    for(i in 1:ng) {
        Tr =c(Tr, out[[i]]$trace)
        B = rbind(B, out[[i]]$bias)
        V = rbind(V, out[[i]]$variance)
    }
    return(list(gamma_x=gamma_ls, trace_y=Tr, Bias=B, Var=V, theta0=theta0))
}

Lmin_bounds <- function(data, init_control) {
                                        # Bounds for the minimum eigenvalue of Fisher information.
    if(init_control$bound.method == "lmin") {
        lam = eigen(fisher(data))$values
        lmin = min(lam)
        warning("Using true min eigenvalue bounds")
        return(list(lower=.5 * lmin, upper=2 * lmin))
    } else if( init_control$bound.method == "yamamoto") {
        print("estimating eigenvalue bounds w/Yamamoto")
        p = ncol(data$X)
        N = nrow(data$X)
        thetas = matrix(0, nrow=0, ncol=p)
        for(i in 1:100) {
            theta = initialize_sgd(data, init_control)
            thetas = rbind(thetas, theta)
        }
        theta_hat = colMeans(thetas)
        trA = 0
        for(i in 1:N) {
            xi = data$X[i,]
            yi = data$Y[i]
            gradi = (yi - data$glm_link(sum(xi * theta_hat))) * xi
            trA = trA + sum(gradi^2) / N
        }
        
                                        # Upper bound.
        ubound = trA / p
        
                                        # Lower bound
        J1 = N * cov(thetas)
        trJ1 = sum(diag(J1))  # trace(J^-1)
        trJ2 = sum(diag(J1 %*% J1)) # trace(J^-2)
                                        # Laguerre bound
        a = p * trJ2 / trJ1^2 - 1
        b = sqrt((p -1) * a)
        lbound = (1 / trJ1) * p / (1 + b)
        
        return(list(lower=lbound, upper=ubound))  
    }
    else {
        print("[WARNING]: did not give valid Lmin_bounds(init_control$bound.method) arg")
    }
}

calculate_gammaStar <- function(data, init_control) {
    print("> Calculating gamma*..")
                                        # TODO: speed this up?
    bounds = Lmin_bounds(data, init_control)
    print(sprintf("> Bounds for gamma: lower = %.2f, upper = %.2f", 
                  bounds$lower, bounds$upper))
    gamma_seq = seq(.5 / bounds$upper, 2 / bounds$lower, length.out=kGammas)
    ## Run MVc
    MVc = par_misgd(data, gamma_seq, init_control = init_control)
    gstar = best_gamma_v2(MVc$gamma_x, MVc$trace_y, bounds$lower, p = ncol(data$X))
    list(gamma=gstar, theta0=MVc$theta0)
}

## 01/30/19
## trying different best_gamma calculation methods
## using just knowledge of lmin
calc_gammaStar_lmin <- function(data, init_control) {
    lam = eigen(fisher(data))$values
    lmin = min(lam)
    theta0 = initialize_sgd(data, init_control=init_control)
    return(list(gamma=(1/lmin), theta0=theta0))
}

## 01/30/19
## using just eigenvalue bound
calc_gammaStar_bd <- function(data, init_control) {
    bounds = Lmin_bounds(data, init_control)
    theta0 = initialize_sgd(data, init_control=init_control)
    return(list(gamma=(2/bounds$lower), theta0=theta0))
}

## 2022 Hwanwoo added this method
### using Inverse Power Iteration
calc_gammaStar_invpower <- function(data, init_control) {
    print("select gammaStar with inverPower method")
### Learning the estimate of the theta_star when estimating the fisher matrix
    theta0 = initialize_sgd(data, init_control=init_control)
    ## Gamma(stepsize) still needs to be tuned
    theta_hat =  isgd_0(theta0, data, gamma=0.01, npass=1, use_permutation = FALSE)
    ##print(sum((data$theta_star - theta_hat)**2))
    data$theta_star = theta_hat
    
    F_star = fisher(data)
    
    ##if (is_not_square_numeric_matrix(data))
    ##  stop("\n'power_method()' requires a square numeric matrix")
    
    set.seed(12345)
    q = rnorm(n = ncol(F_star))
    q = q/(sqrt(sum(q*q)))
    q_prev = q
    
    eps = 1e-3
    error = max(abs(q))
    while(error > eps){
        z = solve(F_star, q, tol = 1e-50)
        q_prev = q
        q = z/sqrt(sum(z*z))
        ##print(sum(abs(q+q_prev)))
        if(sum(abs(q+q_prev)) < eps) break
        error = max(abs(q-q_prev))
    }
    lambda_min = t(q)%*%F_star%*%q
    theta0 = initialize_sgd(data, init_control=init_control)
    return(list(gamma=(1/lambda_min[1,1]), theta0=theta0))
}

## 01/30/19 wrapper for different calc gammaStar methods
calculate_gammaStar_wrapper <- function(data, init_control) {
    
    if (init_control$gamma.method == "heuristic") {
        warning("Using heurisitic best gamma selection")
        return(calculate_gammaStar(data, init_control))
    } else if(init_control$gamma.method == "bound") {
        warning("Using bound best gamma selection")
        return(calc_gammaStar_bd(data, init_control))
    } else if (init_control$gamma.method == "lmin") {
        warning("Using 1/lmin for best gamma selection")
        return(calc_gammaStar_lmin(data, init_control))
    } else if (init_control$gamma.method == "ipower") {
        warning("Using inverse iteration for best gamma selection")
        return(calc_gammaStar_invpower(data, init_control))
    } else {
        stop("calc gammaStar method not yet implemented")	
    }
}

exact_inf <- function(data, sgd_control) {
    ## Function to compute CIs as "Plus/Minus" the learning rate
    ## Main method of exact inference.
    ## sgd_control = list(gamma, theta0)

    N = nrow(data$X)
    p = ncol(data$X)
    
    gamma_star = sgd_control$gamma
    ## 1. Estimate theta.
    theta_hat = isgd_0(sgd_control$theta0, data, gamma=gamma_star, 
                       npass=1, use_permutation = FALSE)
    psi = get_model_psi(data, theta_hat)
    
    ## 2. Variance estimate
    if(!is.null(sgd_control$use_gamma_half)) {
        warning("Using gamma_star / 2 in variance")
        gamma_star = gamma_star / 2
    }
    V = psi * gamma_star * rep(1, p) / N
    ## 3. theta +- sqrt(gamma/N)
    upper.ci = theta_hat + 1.96 * sqrt(V)
    lower.ci = theta_hat - 1.96 * sqrt(V)
    
    ## 4. Wald pivots 
    pivots = abs(theta_hat) / sqrt(V)

    return(list(est=theta_hat, ci=cbind(lower.ci, upper.ci), pivots=pivots))
}

glm_inf <- function(data) {
    if(data$model == "gaussian") {
        glm.fit = glm(data$Y ~ 0 + data$X)
    } else if (data$model == "binomial") {
        glm.fit = glm(data$Y ~ 0 + data$X, family=binomial)
    } else if (data$model == "poisson") {
        glm.fit = glm(data$Y ~ 0 + data$X, family=poisson)
    } else {
        stop("glm_inf(): data$model_name not gaussian or logistic, need to implement")
    }
    sum.coef = summary(glm.fit)$coef
    est = sum.coef[,1]
    upper.ci = est + 1.96 * sum.coef[,2]
    lower.ci = est - 1.96 * sum.coef[,2]
    pivots = abs(est) / sum.coef[,2]
    return(list(est=est, ci=cbind(lower.ci,upper.ci), pivots=pivots))
}
