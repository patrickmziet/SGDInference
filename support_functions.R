make_file_name <- function(x, sim_type) {
    paste0(sim_type,
           "_", "family", x$family,
           "_", "sigma_x", x$sigma_x,
           "_", "true_param", x$true_param, 
           "_", "n", x$n,
           "_", "p", x$p,
           "_", "s", x$s,
           "_", "hard", x$fixed,
           ".rda")
}

get_gamma <- function() {
    cnt <- 0
    isNegative <- TRUE
    while (isNegative) {
        set.seed(i + (cnt - 1) * nsim)
        data <- gen_data(model_name = as.character(sett$family),
                         sigma_x = as.character(sett$sigma_x),
                         N = sett$n,
                         p = sett$p,
                         s = sett$s,
                         true_param = as.character(sett$true_param))

        init_control <- list()
        init_control$gamma.method <- "ipower"
        init_control_def <- default_init(sett$p)
        for(ii in names(init_control)) init_control_def[[ii]] <- init_control[[ii]]
        init_control <- init_control_def
        sgd_control <- calculate_gammaStar_wrapper(data, init_control)
        isNegative <- sgd_control$gamma < 0
    }
    return(list(sgd_control = sgd_control, data = data))
}

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

