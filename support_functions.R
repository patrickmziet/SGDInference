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



get_gamma <- function(i, sett, imth = "ipower") {
    cnt <- 0
    isNegative <- TRUE
    while (isNegative) {
        cnt <- cnt + 1
        set.seed(i + (cnt - 1) * nsim)
        data <- gen_data(model_name = as.character(sett$family),
                         sigma_x = as.character(sett$sigma_x),
                         N = sett$n,
                         p = sett$p,
                         s = sett$s,
                         true_param = as.character(sett$true_param))

        init_control <- list()
        init_control$gamma.method <- imth
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

## Functions from glmodsel
get_file <- function(sim_axis, family, link, n, p, s, rho, base_name) {
    ss <- paste0(paste(base_name,
                       sim_axis,
                       paste0("family-", family),
                       paste0("link-", link),
                       paste0("n", n),
                       paste0("p", p),
                       paste0("s", s),
                       paste0("rho", sub("\\.", "", format(as.numeric(rho), nsmall = 2))), sep = "-"),
                 ".rda")
    ## Remove first occurence of "-" in the above string
    sub("-", "", ss, fixed = TRUE)
}

hline <- function(y = 0, color = "black", dash = "dot") {
    list(
        type = "line",
        x0 = 0,
        x1 = 1,
        xref = "paper",
        y0 = y,
        y1 = y,
        line = list(color = color, dash = dash),
        opacity = 0.2
    )
}

get_object <- function(img, object) {
    if (file.exists(img)) {
        load(img)
        get(object)
    } else {
        stop("Setting is not available")
    }
}

make_file <- function(sim_axis, family, link, n, p, s, rho=0, base_name="summary_sgd_") {
    ff1 <- get_file(sim_axis = sim_axis,
                    family = family,
                    link = link,
                    n = n,
                    p = p,
                    s = s,
                    rho = rho,
                    base_name = base_name)
    file.path(outputs_path, ff1)
}

make_file_name_msglm <- function(ll, base_name, fig = FALSE) {
    ftype <- ".rda"
    if (isTRUE(fig)) ftype <- ".pdf"
    paste0(base_name,
           "_","family", ll$family,
           "_","link", ll$link,
           "_", "n", ll$n,
           "_","p", ll$p,
           "_","s", ll$s,
           "_","rho", ll$rho,
           "_", "nu", ll$nu,
           ftype)
}


conf_set <- function(res, lev) {
    n_reps <- ncol(res)
    active_sets(res) |>
    table() |>
    sort(decreasing = TRUE) |>
    (\(x) cumsum(c(x)/n_reps))() |>
    (\(x) x[x <= lev])()
}

active_sets <- function(mat) {
    if (is.null(rn <- rownames(mat))) {
        rn <- 1:ncol(mat)
    }
    apply(mat, 2, function(x) paste0(rn[which(x)], collapse = ","))
}

make_row_tex <- function(vec) {
    row_tex <- paste(vec, collapse = " & ")
    row_tex <- paste0(row_tex, " \\\\ \n")
    return(row_tex)
}
