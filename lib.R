expit <- function(x) sapply(x, function(i) { if(i > 60) return(1); if(i < -60) return(0); return(exp(i)/(1 + exp(i)))})

par_rmvnorm <- function(N=1e5,p=1e3) {
  Sigma_X = toeplitz(0.5^abs(1-1:p))
  draw_X <- function(i) {
    i=1
    return(rmvnorm(1e3, mean=rep(0,p), sigma=Sigma_X))
  }
  X = mclapply(1:100, function(i) draw_X(i))
  return(X)
}

du_expit <- function(x) sapply(x, function(i) { if(i > 60) return(1); if(i < -60) return(0); return(exp(-i)/((1 + exp(-i))^2))})
L2 = function(x, y) sqrt(sum((x-y)^2))
L2.p = function(x, y) {
  sqrt(sum((x-y)^2)) / length(x) 
}

fast_sqrtm <- function(M) {
  p = ncol(M)
  I = diag(p)
  U = matrix(1, nrow=p, ncol=p)
  # Check for equicorr matrix.
  equicor_obj <- function(par) {
    Mhat = par[1] * diag(p) + par[2]* U
    sum((Mhat - M)^2)
  }
  out = optim(par=c(0, 0), fn = equicor_obj, 
              lower = c(-100, -100), upper=c(100, 100),
              method="L-BFGS-B")
  ## special case
  if(all(M == I)) {
    out$par = c(1, 0)
    out$value = 0
  }
  if(out$value < 1e-6) {
    ## This is indeed equicor
    c1 = out$par[1]
    c2 = out$par[2]
    a = sqrt(c1)
    all_b = polyroot(c(-c2 * a, 2 * c1, a * p))
    b = Re(all_b[1])
    if(all(eigen(a * I + all_b[2] * U)$values > 0)) {
      b = all_b[2]
    }
    return(a * I + b * U)
  }
  
}

## Compact code for implicit SGD.
##  example:
##    d = gen_data(model_name="logistic", N=1e4)
##    out = implicit_sgd(d, verbose=T)
##
gen_data_star <- function(theta_star, X, model_name="gaussian", sigma_noise=1) {
  # unload experimental params here. TODO: extend to match our current experimental setup.
  # theta_star = (-1)^seq(1, p) * thetaStar_coeff * exp(-0.75 * seq(1, p))
  ##
  N = nrow(X)
  p = ncol(X)
  pred = X %*% theta_star  # predictor.
  # 
  glm_link = NA # glm link function. 
  Y = NA
  #
  if(model_name == "gaussian") {
    glm_link = function(x) x
    du_glm_link = function(x) 1
    noise = rnorm(N, sd=sigma_noise)
    Y = glm_link(pred) + noise
  } else if(model_name == "poisson") {
    glm_link = function(x) exp(x)
    du_glm_link = function(x) exp(x)
    Y = rpois(N, lambda=glm_link(pred))
  } else if(model_name=="binomial") {
    glm_link = function(x) expit(x)
    du_glm_link = function(x) du_expit(x)
    Y = rbinom(N, size=1, prob=glm_link(pred))
  } else {
    stop("model not implemented.")
  }
  
  return(list(model=model_name, X=X, Y=Y, glm_link=glm_link, 
              du_glm_link=du_glm_link, theta_star=theta_star))
}

gen_data <- function(model_name="gaussian", N=1000, p=20, 
                     sigma_x="id", rho=.15, theta_coeff=1,
                     sigma_noise=1, true_param="classic",
                     sigma_sqrt=NA) {
  # unload experimental params here. TODO: extend to match our current experimental setup.
  
  # classic theta star
  if(true_param=="classic") {
    theta_star = (-1)^seq(1, p) * 2 * exp(-0.7 * seq(1, p))
  } else if(true_param=="inc") {
    # increasing theta star 
    theta_star = seq(-3, 3, length.out=p)
  } else if(true_param=="usc") {
    theta_star = seq(0, 1, length.out=p)
  } else if(true_param=="sparse") {
    theta_star = (-1)^seq(1, p) * 2 * exp(-0.7 * seq(1, p))
    i = tail(1:p, 0.8 * p)
    theta_star[i] <- 0
  } else if(true_param=="single") {
    theta_star <- rep(0, p)
    theta_star[1] = theta_coeff
  } else {
    stop("implement additional theta_star definitions")
  }
  Sigma_X = NA
  if(sigma_x=="id") {
    Sigma_X = diag(p)
  } else if(sigma_x=="equicor") {
    Sigma_X = (1 - rho) * diag(p) + rho * matrix(1, nrow=p, ncol=p) 
  } else if (sigma_x=="ill_cond") {
    lam = seq(0.5, 100, length.out=p)
    # Q = qr.Q(qr(u %*% t(u)))
    Sigma_X = diag(lam)
  } else if (sigma_x=="toeplitz") {
    Sigma_X = toeplitz(0.5^seq(0, p-1))
  } else if (sigma_x=="Xdiffvar") {
    #lam = seq(0.1, rho, length.out=p)
    lam = seq(rho, 1, length.out=p)
    # Q = qr.Q(qr(u %*% t(u)))
    Sigma_X = diag(lam)
  } else {
    # TODO: implement different Sigma definitions here.
    stop("implement additional Sigma_X definition")
  } 
  if(nrow(Sigma_X) != ncol(Sigma_X)) {
    stop("Need to feed square matrix here.")
  }
  ##
  if (any(is.na(sigma_sqrt))) {
    S = Sigma_X
    if(sigma_x != "id") {
      S =  sqrtm(Sigma_X)
    }
  } else {
    S = sigma_sqrt
  }
  Z = matrix(rnorm(N * p), ncol=p)
  X = t(S %*% t(Z))
  # rmvnorm(N, mean=rep(0, p), sigma = Sigma_X) #covariates
  pred = X %*% theta_star  # predictor.
  # 
  glm_link = NA # glm link function. 
  Y = NA
  #
  if(model_name == "gaussian") {
    glm_link = function(x) x
    du_glm_link = function(x) 1
    noise = rnorm(N, sd=sigma_noise)
    Y = glm_link(pred) + noise
  } else if(model_name == "poisson") {
    ## TODO: Poisson Y should not be crazy.
    glm_link = function(x) exp(x)
    du_glm_link = function(x) exp(x)
    Y = rpois(N, lambda=glm_link(pred))
  } else if(model_name=="binomial") {
    glm_link = function(x) expit(x)
    du_glm_link = function(x) du_expit(x)
    Y = rbinom(N, size=1, prob=glm_link(pred))
  } else {
    stop("model not implemented.")
  }
  
  return(list(model=model_name, Sigma_X = Sigma_X, X=X, Y=Y, glm_link=glm_link, 
              du_glm_link=du_glm_link, theta_star=theta_star, Sigma_X_sqrt=S,
              cov=sigma_x, true_param=true_param))
}

fisher <- function(data) {
  ## Calculates Fisher information matrix.
  n = nrow(data$X)
  p = ncol(data$X)
  Fi = matrix(0, nrow=p, ncol=p)
  for(i in 1:n) {
    xi = data$X[i, ]
    yi = data$Y[i]
    etai = sum(xi * data$theta_star)
    gradi = matrix((yi - data$glm_link(etai)) * xi, ncol=1)
    Fi = Fi + gradi %*% t(gradi)
  }
  return(Fi / n)
}

bound_over <- function(kappa, alpha=0.05) {
  ## upper bound for over-coverage based on condition number kappa and (1-alpha) nominal
  z = qnorm((1 - alpha/2))
  return(alpha - 2 * pnorm(-z * sqrt(2 - 1/kappa)))
}

bound_under <- function(kappa, alpha=0.05) {
  ## lower bound for over-coverage based on condition number and (1-alpha) nominal
  z = qnorm(1 - alpha/2)
  return(alpha - 2 * pnorm(-z * sqrt(kappa)))
}


