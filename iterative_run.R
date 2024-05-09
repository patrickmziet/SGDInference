## Run model selection and estimation iteratively
library("colorspace")
library("progress")
## Load functions
source("simulations.R")
source("exact_inference.R")
source("iterative_functions.R")

## Hyperameters
model_name <- "binomial"
true_param <- "hht"
nn <- 200
pp <- 10
ss <- 5
iB <- pp * 2 # initial batch size
B <- 1 # batch size
lvl <- 0.95 # confidence set

set.seed(1)
res <- run_sim(model_name = model_name,
               true_param = true_param,
               nn = nn,
               pp = pp,
               ss = ss,
               iB = iB,
               B = B,
               lvl = lvl,
               reest = FALSE)

print_ci(res$mdls)


## Test out oracle properties
nsim <- 300
pb <- progress_bar$new(total = nsim)
coefs <- matrix(0, ncol = pp, nrow = nsim)
ses <- matrix(0, ncol = pp, nrow = nsim)
for (i in seq.int(nsim)) {
    set.seed(i)
    rr <- run_sim(model_name = model_name,
                  true_param = true_param,
                  nn = nn,
                  pp = pp,
                  ss = ss,
                  iB = iB,
                  B = B,
                  lvl = lvl,
                  reest = FALSE,
                  verbose = FALSE,
                  prog = FALSE)
    ## Get model
    Jhat <- with(rr, Jhat[nrow(Jhat), ])
    mm <- paste("X", (1:pp)[Jhat], sep = "")
    ## Re-estimate
    dd <- with(rr, data.frame(data$Y, data$X))
    attr(dd, "formula") <- paste0(c("y ~ -1", mm), collapse = " + ")
    names(dd) <- c("y", paste("X", 1:pp, sep = ""))
    mod <- glm(attr(dd, "formula"),
               family = binomial(link = "logit"),
               data = dd,
               method = "brglmFit")
    coefs[i, Jhat] <- summary(mod)$coefficients[,1]
    X <- as.matrix(dd[, -c(1)])
    XJ <- X[, 1:ss] 
    eta <- X %*% rr$data$theta_star
    phi <- 1
    mu <- exp(-eta) / (1 + exp(-eta))
    W <- diag(as.vector(mu * (1 - mu)))
    FI <- t(XJ) %*% W %*% XJ / phi
    FIinv <- solve(FI)            
    ses[i, c(rep(TRUE, ss), rep(FALSE, pp - ss))] <- sqrt(diag(FIinv))
    pb$tick()
}


coefs
ses

## Plot
## A function to estimate and plot mixtures of continuous densities and point masses at zero
makePlot <- function(est, ref0, ylim = NULL, xlim = NULL, main = "", col = NULL, border = NULL,
                     xlab = "", ylab = "", brks = 20, cl2 = "gray") {
    inds <- ref0
    prop0 <- mean(inds)
    out <- hist(est[!inds], breaks = brks, plot = FALSE)
    out$density <- out$density * (1 - prop0)
    if (is.null(ylim)) {
        ylim <- range(c(prop0, out$density))
    }
    if (is.null(xlim)) {
        xlim <- range(out$breaks)
    }
    plot(c(0, 0),
         type = "l",
         col = "white",
         xlim = xlim,
         ylim = ylim,
         main = main,
         xlab = xlab,
         ylab = ylab)

    ## Plot the histogram for values where inds is FALSE
    plot(out, freq = FALSE, add = TRUE, col = col, border = border)

    if (prop0 > 0) {
        ## Histogram for values where inds is TRUE
        out2 <- hist(est[inds], plot = FALSE, breaks = 4)
        out2$density <- out2$density * prop0

        ## Overlay the histogram for values where inds is TRUE
        plot(out2, freq = FALSE, add = TRUE, col = cl2, border = "black")
    }

    ## Overlay standard normal
    xx <- seq(xlims[1], xlims[2], length.out = 100)
    lines(x = xx, y = dnorm(x = xx) * (1 - prop0))
}


## wald_stats <- with(res, t(apply(coefs, 1, FUN = function(x) x - theta_star)) / ses)[ ,1:ss]
## wald_stats <- with(res, (coefs - 1) / ses)[,1:ss]
wald_stats <- ((coefs - 1) / ses)[,1:ss]

colsT <- rainbow_hcl(1, alpha = 0.7)
colsF <- rainbow_hcl(1, alpha = 1)

par(mfrow = c(3, 2))
for (j in seq.int(ss)) {
    xlims <- range(wald_stats[, j]) + c(-1, 1)
    makePlot(wald_stats[, j],
             ref0 = rep(FALSE, nrow(wald_stats)),
             ylim = c(0, 1),
             xlim = xlims,
             col = colsT[1],
             border = colsF[1],
             main = paste0("Variable ", j))
}
