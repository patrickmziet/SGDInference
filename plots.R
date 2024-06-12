library("colorspace")
library("dplyr")
library("ggplot2")
library("patchwork")
source("support_functions.R")
source("iterative_functions.R")
## Define colours
base_colours <- colorspace::qualitative_hcl(6)
## Define paths
experiment_path <- "~/repos/SGDInference"
outputs_path <- file.path(experiment_path, "outputs")

## 1. Convergence to the true model
load(file = file.path(outputs_path, "cns_sgd_settings.rda"))

settings_df <- as.data.frame(settings) |>
mutate(sparsity = ifelse(s/p > 5/40, "dense", "sparse")) |>
mutate(method = ifelse(fixed, "hard-BIC", "soft-BIC")) |>
mutate(method_label = paste0(sigma_x, "-", method))

## Group data you want to plot together
grouped_settings <- settings_df |>
group_by(family, p, s, sparsity) |>
summarize()

## Group figures you want together, in this case we will have two in each figure group
figure_groups <- grouped_settings |>
group_by(family, sparsity) |>
summarize()

metrics <- c("f1", "nonzeros")
save_fig <- TRUE
## Loop through figure groups
for (k in seq.int(nrow(figure_groups))){
    ff <- as.character(figure_groups$family[k])
    ss <- figure_groups$sparsity[k]
    cat("Figures for", ff, ss, "\n")
    ## Filter for these
    sim_settings <- settings_df |> filter(family == ff & sparsity == ss)
    results <- data.frame()
    for (i in seq.int(nrow(sim_settings))) {
        sett <- sim_settings[i,]
        sett_nm <- make_file_name(sett, "cns_sgd")
        load(file = file.path(outputs_path, sett_nm))
        rr <- as.data.frame(sett[rep(1, length(stats$sim_vals)),])
        rr$n <- stats$sim_vals
        rr$nonzeros <- sapply(stats$mdls, function(model) {
            mean(rowSums(model))
        })
        true_mod <- with(stats, c(rep(TRUE, sett$s), rep(FALSE, sett$p - sett$s)))
        rr$f1 <- sapply(stats$mdls, function(model) {
            TP <- rowSums(t(apply(X = model, MARGIN = 1, FUN = function(x) x & true_mod))) # True positives
            FP <- rowSums(t(apply(X = model, MARGIN = 1, FUN = function(x) x & !true_mod))) # False positives
            FN <- rowSums(t(apply(X = model, MARGIN = 1, FUN = function(x) !x & true_mod))) # False negatives
            precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
            recall <- ifelse(TP + FN == 0, 0, TP / (TP + FN))            
            f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall))
            mean(f1)
        })
        rownames(rr) <- NULL
        results <- rbind(results, rr)
    }

    num_metrics <- length(unique(sim_settings$p))
    plots <- as.list(rep(NA, 2 * num_metrics)) # 2 for f1 and nonzeros
    names(plots) <- paste0(metrics, c(rep("-1", num_metrics), rep("-2", num_metrics)))
    two_settings <- sim_settings |> group_by(p, s) |> summarize()
    for (j in seq.int(nrow(two_settings))) {
        current_data <- results |> filter(p == two_settings$p[j])
        for (met in metrics) {
            met_nm <- paste0(met, ifelse(j == 1, "-1", "-2"))
            plots[[met_nm]] <- ggplot(current_data,
                                      aes_string(x = "log(n)",
                                                 y = met,
                                                 color = "method_label")) +
                geom_point() +
                geom_line(size = 1) +
                xlab("log(n)") +
                ylab(met) +
                labs(color = "Method") +
                ## scale_color_manual(values = base_colours) +
                theme(legend.position = "bottom") +
                theme_bw() +
                theme(
                    axis.title = element_text(size = 14),
                    axis.text = element_text(size = 12),
                    ) 
            
        }
    }
    plots[["f1-1"]] <- plots[["f1-1"]] + geom_hline(aes(yintercept = 1))
    plots[["f1-2"]] <- plots[["f1-2"]] + geom_hline(aes(yintercept = 1))

    plots[["nonzeros-1"]] <- plots[["nonzeros-1"]] + geom_hline(aes(yintercept = two_settings$s[1]))
    plots[["nonzeros-2"]] <- plots[["nonzeros-2"]] + geom_hline(aes(yintercept = two_settings$s[2]))

    fig <- (plots[["f1-1"]] + plots[["f1-2"]]) / (plots[["nonzeros-1"]] + plots[["nonzeros-2"]]) &
        theme(legend.position = "bottom",
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 14))
    fig <- fig + plot_layout(guides = "collect")
    fig
    if (save_fig) {
        fig_file <- file.path(outputs_path,
                              paste0("sgd_",
                                     figure_groups$family[k],
                                     "_", figure_groups$sparsity[k], ".pdf"))
        pdf(fig_file, width = 10, height = 5)
        print(fig)
        dev.off()
    } else {
        print(fig)
    }
}

## Convergence of simulation settings from MSGLM chapter
## base_col <- colorspace::qualitative_hcl(16, l = 50)
## Plotting performance metrics
## Methods to use
methodsSelection <- c("ms_lasso(relax = TRUE, family = family, link = link)",
                      "ms_thres(w_type = 'bic', g_type = 'enhanced2', fixed = FALSE, family = family, link = link, fit_method = 'brglmFit')",
                      "ms_thres(w_type = 'bic', g_type = 'enhanced2', fixed = TRUE, family = family, link = link, fit_method = 'brglmFit')",
                      "sgd_soft",
                      "sgd_hard")

gg_base_col <- c("#00BFC4", "#619CFF", "#F564E3", "#1AFF1A", "#4B0092")
methodNms <- c("rlasso", "BIC soft", "BIC hard", "BIC SGD soft", "BIC SGD hard")
color_mapping <- setNames(gg_base_col, methodsSelection)

## Adjust settings here to get the plots you want
ff <- "poisson" ## Family: binomial, poisson
ll <- "log"    ## Link function: logit, log
ww <- "dense"    ## dense or sparse setting

if (ww == "sparse") {
    two_settings <- list(
        list(pp = 100 , ss = 5), ## First setting
        list(pp = 40, ss = 5)  ## Second setting
    )
} else if (ww == "dense") {
    two_settings <- list(
        list(pp = 10 , ss = 5), ## First setting
        list(pp = 40, ss = 25)  ## Second setting
    )
}

for (st in seq.int(two_settings)) {
    print(st)
    file1 <- with(two_settings,
                  make_file("Observations",
                       ff,
                       ll,
                       2500,
                       two_settings[[st]]$pp,
                       two_settings[[st]]$ss)
                  )
    summaries <- get_object(file1, "summaries") |>
    mutate(id = 1:n(),
           method = factor(selector_call))## , colours = base_col[method])
    metrics <- unique(summaries$metric)
    if (st == 1) {
        plots <- as.list(rep(NA, length(two_settings) * length(metrics)))
        names(plots) <- paste0(metrics, c(rep("-1", length(metrics)), rep("-2", length(metrics))))
    }
    for (met in metrics) {
        met_nm <- paste0(met, ifelse(st == 1, "-1", "-2"))
        current_data <- summaries |> subset(metric == met) |> subset(method %in% methodsSelection)
        current_data$method <- factor(current_data$method, levels = methodsSelection)
        plots[[met_nm]] <- ggplot(current_data, aes(x = log(n), y = mean, color = method)) +
            geom_point() +
            geom_line(size = 1) +
            xlab("log(n)") +
            ylab(met) +
            labs(color = "Method") +
            scale_color_manual(values = color_mapping,
                               labels = c("rlasso",
                                          "BIC soft",
                                          "BIC hard",
                                          "BIC SGD soft",
                                          "BIC SGD hard")) + 
            theme_bw() +
            theme(
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 12),
                ## legend.position = "none",
                )
    }
}

plots[["f1-1"]] <- plots[["f1-1"]] + geom_hline(aes(yintercept = 1))
plots[["f1-2"]] <- plots[["f1-2"]] + geom_hline(aes(yintercept = 1))

plots[["nonzeros-1"]] <- plots[["nonzeros-1"]] + geom_hline(aes(yintercept = two_settings[[1]]$ss))
plots[["nonzeros-2"]] <- plots[["nonzeros-2"]] + geom_hline(aes(yintercept = two_settings[[2]]$ss))

fig <- (plots[["f1-1"]] + plots[["f1-2"]]) / (plots[["nonzeros-1"]] + plots[["nonzeros-2"]]) &
     theme(legend.position = "bottom",
           legend.text = element_text(size = 12),
           legend.title = element_text(size = 14))
fig <- fig + plot_layout(guides = "collect")
fig

fig_file <- paste0("sgd_comp_",ff,"_",ww,".pdf")
fig_file <- file.path(outputs_path, fig_file)
save_fig <- TRUE
if (save_fig) {
    pdf(fig_file, width = 10, height = 5)
    print(fig)
    dev.off()
} else {
    print(fig)
}


## 2. Oracle property
load(file = file.path(outputs_path, "oracle_sgd_settings.rda"))
## Plots for some settings
colsT <- rainbow_hcl(1, alpha = 0.7)
colsF <- rainbow_hcl(1, alpha = 1)

settings_df <- as.data.frame(settings) |>
mutate(sparsity = ifelse(s/p > 5/40, "dense", "sparse")) |>
mutate(method = ifelse(fixed, "hard-BIC", "soft-BIC")) |>
mutate(method_label = paste0(sigma_x, "-", method))

## Group data you want to plot together
grouped_settings <- settings_df |>
group_by(family, fixed, n, p, s) |>
summarize()

## Single plots
sigma_xs <- c("id", "toeplitz", "equicor")
save_fig <- TRUE
for (k in seq.int(nrow(grouped_settings))) {
    gs <- grouped_settings[k,]
    ff <- gs$family
    hh <- gs$fixed
    pp <- gs$p
    ss <- gs$s
    nn <- gs$n

    wald_stats <- list()
    for (i in seq_along(sigma_xs)) {
        sett <- settings_df |> filter(family == ff
                                      & n == nn
                                      & p == pp
                                      & s == ss
                                      & sigma_x == sigma_xs[i]
                                      & fixed == hh)

        sett_nm <- make_file_name(sett, "oracle_sgd")
        load(file = file.path(outputs_path, sett_nm))
        wald_stats[[sigma_xs[i]]] <- with(stats, (
            t(apply(X = coefs,
                    MARGIN = 1,
                    FUN = function(x) x - theta_star)) / ses)[,1:sett$s]
            )
    }

    fig_file <- file.path(outputs_path,
                          paste0("oracle_sgd",
                                 "_family", ff,
                                 "_fixed", hh,
                                 "_n", nn,
                                 "_p", pp,
                                 "_s", ss,
                                 ".pdf"))

    if (save_fig) pdf(fig_file, width = 10, height = 10/sqrt(2))
    n_mthds <- length(wald_stats)
    par(mfcol = c(n_mthds, 6), mar = c(2,2,2,2), oma = c(3,1,1,1))
    colsT <- rainbow_hcl(n_mthds, alpha = 0.7)
    colsF <- rainbow_hcl(n_mthds, alpha = 1)
    mthds_full <- sigma_xs
    for (k in seq.int(wald_stats)) {
        plot(0, 0, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
        text(0, 0, mthds_full[k], cex = 1.5)
    }
    nsim <- nrow(stats$coefs)
    for (j in seq.int(sett$s)) {
        xlims <- range(sapply(X = wald_stats, FUN = function(x) range(x[,j]))) + c(-1, 1)
        mn <- paste0("Variable ", j)
        makePlot(wald_stats[["id"]][,j],
                 ref0 = rep(FALSE, nsim),
                 ylim = c(0, 1),
                 xlim = xlims,
                 col = colsT[1],
                 border = colsF[1],
                 main = mn)
        makePlot(wald_stats[["toeplitz"]][,j],
                 ref0 = rep(FALSE, nsim),
                 ylim = c(0, 1),
                 xlim = xlims,
                 col = colsT[2],
                 border = colsF[2],
                 main = "")
        makePlot(wald_stats[["equicor"]][,j],
                 ref0 = rep(FALSE, nsim),
                 ylim = c(0, 1),
                 xlim = xlims,
                 col = colsT[3],
                 border = colsF[3],
                 main = "")
    }
    mtext("Distribution of non-zero coefficients",
          side = 1,
          outer = TRUE,
          line = 1)
    if (save_fig) dev.off()
}


## Tables for other settings
table_groups <- settings_df |>
group_by(family, fixed, p, s) |>
summarise()

save_fig <- TRUE
table_groups_ks_tests <- list()
sigma_xs <- c("id", "toeplitz", "equicor")
for (k in seq.int(nrow(table_groups))) {
    gs <- table_groups[k,]
    ff <- gs$family
    hh <- gs$fixed
    pp <- gs$p
    ss <- gs$s
    
    gs_nm <- paste0(
        "family", ff, "_", 
        "fixed", hh, "_", 
        "p", pp, "_", 
        "s", ss
    )

    n_vals <- unique(unlist(as.vector(settings_df |>
                                      subset(family == ff & fixed == hh & p == pp & s == ss)
                                      |> select(n))))
    wald_stats <- list()
    ks_tests <- list()
    for (l in seq.int(length(n_vals))) {
        nn <- n_vals[l]
        for (i in seq_along(sigma_xs)) {
            sett <- settings_df |> filter(family == ff
                                          & n == nn
                                          & p == pp
                                          & s == ss
                                          & sigma_x == sigma_xs[i]
                                          & fixed == hh)
            sett_nm <- make_file_name(sett, "oracle_sgd")
            print(sett_nm)
            load(file = file.path(outputs_path, sett_nm))
            it_nm <- paste0(sigma_xs[i], "-", nn)
            wald_stats[[it_nm]] <- with(stats, (
                t(apply(X = coefs,
                        MARGIN = 1,
                        FUN = function(x) x - theta_star)) / ses)[,1:sett$s]
                )
            ks_tests[[it_nm]] <- apply(X = wald_stats[[it_nm]],
                                       MARGIN = 2,
                                       FUN = function(x) {
                                           obj <- ks.test(x, "pnorm", mean = mean(x), sd = sd(x))
                                           list(statistic = obj$statistic, p.value = obj$p.value)
                                       })
        }

    }
    table_groups_ks_tests[[gs_nm]] <- ks_tests

    ## Set up plot
    fig_file <- file.path(outputs_path,
                          paste0("oracle_sgd_grouped",
                                 "_family", ff,
                                 "_fixed", hh,
                                 "_p", pp,
                                 "_s", ss,
                                 ".pdf"))

    n_mthds <- length(sigma_xs)
    ss <- 5
    if (save_fig) {
        if (ss > 5) {
            pdf(fig_file, width = 14, height = 14/sqrt(2), paper = "a4r")
            par(mfcol = c(length(n_vals) * n_mthds, ss + 1), mar = c(1,1,1,1), oma = c(3,1,1,1))
        } else {
            pdf(fig_file, width = 10, height = 12/sqrt(2))
            par(mfcol = c(length(n_vals) * n_mthds, ss + 1), mar = c(2,2,2,2), oma = c(3,1,1,1))
        }
    }
    ## if (save_fig) pdf(fig_file, width = 10, height = 10/sqrt(2))
    ## par(mfcol = c(length(n_vals) * n_mthds, ss + 1), mar = c(2,2,2,2), oma = c(3,1,1,1))

    ## par(mfcol = c(n_mthds, 6), mar = c(2,2,2,2), oma = c(3,1,1,1))

    colsT <- rainbow_hcl(n_mthds, alpha = 0.7)
    colsF <- rainbow_hcl(n_mthds, alpha = 1)

    mthds_full <- paste0(rep(sigma_xs, length(n_vals)), "-", rep(n_vals, each = 3))
    for (k in seq.int(wald_stats)) {
        plot(0, 0, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
        text(0, 0, mthds_full[k], cex = 1.5)
    }
    nsim <- nrow(wald_stats[[1]])
    for (j in seq.int(ss)) {
        xlims <- range(sapply(X = wald_stats, FUN = function(x) range(x[,j]))) + c(-1, 1)
        for (ii in seq_along(n_vals)) {
            mn <- paste0("Variable ", j)
            if (ii > 1) mn <- ""
            nn <- n_vals[ii]
            makePlot(wald_stats[[paste0("id-",nn)]][,j],
                     ref0 = rep(FALSE, nsim),
                     ylim = c(0, 1),
                     xlim = xlims,
                     col = colsT[1],
                     border = colsF[1],
                     main = mn)
            makePlot(wald_stats[[paste0("toeplitz-", nn)]][,j],
                     ref0 = rep(FALSE, nsim),
                     ylim = c(0, 1),
                     xlim = xlims,
                     col = colsT[2],
                     border = colsF[2],
                     main = "")
            makePlot(wald_stats[[paste0("equicor-", nn)]][,j],
                     ref0 = rep(FALSE, nsim),
                     ylim = c(0, 1),
                     xlim = xlims,
                     col = colsT[3],
                     border = colsF[3],
                     main = "")
        }
    }
    mtext("Distribution of non-zero coefficients",
          side = 1,
          outer = TRUE,
          line = 1)
    if (save_fig) dev.off()
}


format_ks_tests <- function(ks_tests) {
    ## Initialize the output table
    output <- data.frame(
        variable = seq.int(length(ks_tests[[1]])),
        stringsAsFactors = FALSE
    )
    
    ## Define the row and column headers
    column_headers <- c("Variable", names(ks_tests))
    
    ## Iterate over each combination of parameters
    for (key in names(ks_tests)) {
        ## Extract the test results
        test_results <- ks_tests[[key]]
        
        ## Create a vector to store formatted results for this key
        formatted_results <- sapply(test_results, function(test) {
            stat <- round(test$statistic, 3)
            p_val <- round(test$p.value, 3)
            if (p_val < 0.001) {
                p_val_str <- "0.001>"
            } else {
                p_val_str <- formatC(p_val, format = "f", digits = 3)
            }
            paste0(formatC(stat, format = "f", digits = 3), " (", p_val_str, ")")            
            ## paste0(formatC(stat, format = "f", digits = 3),
            ##        " (", formatC(p_val, format = "e", digits = 2), ")")
        })
        
        ## Add the formatted results to the output table
        output[[key]] <- formatted_results
    }

    ## Extract and sort the column names dynamically
    col_names <- names(ks_tests)

    ## Define the desired prefix order
    prefix_order <- c("id", "toeplitz", "equicor")

    ## Function to extract the prefix and numeric suffix
    extract_prefix_suffix <- function(name) {
        parts <- unlist(strsplit(name, "-"))
        prefix <- parts[1]
        suffix <- as.numeric(parts[2])
        return(c(prefix, suffix))
    }

    ## Extract prefixes and suffixes
    prefix_suffix <- t(sapply(col_names, extract_prefix_suffix))

    ## Create a data frame for sorting
    sort_df <- data.frame(
        name = col_names,
        prefix = factor(prefix_suffix[,1], levels = prefix_order),
        suffix = as.numeric(prefix_suffix[,2]),
        stringsAsFactors = FALSE
    )

    ## Sort the data frame by prefix and suffix
    sorted_df <- sort_df[order(sort_df$prefix, sort_df$suffix), ]

    ## Extract the sorted column names
    sorted_col_names <- sorted_df$name

    ## Prepend 'variable' to the sorted column names
    desired_order <- c("Variable", sorted_col_names)

    ## Set the column names and reorder columns
    colnames(output) <- c("Variable", col_names)
    output <- output[, desired_order, drop = FALSE]
    
    ## ## Print the table
    print(output, row.names = FALSE)
}

for (i in seq_along(names(table_groups_ks_tests))) {
    nm <- names(table_groups_ks_tests)[i]
    print(nm)
    format_ks_tests(table_groups_ks_tests[[nm]])
}

## Oracle settings for comparison
load(file = file.path(outputs_path, "oracle_comp_sgd_settings.rda"))
settings_df <- as.data.frame(settings) |> mutate(kappa = s/p)
settings_df <- settings_df |> group_by(kappa, family)
kappa_groups <- settings_df |> summarise()
save_fig <- TRUE
for (k in seq.int(nrow(kappa_groups))) {
    ff <- as.character(kappa_groups$family[k])
    kk <- kappa_groups$kappa[k]
    cat("Figures for", ff, kk, "\n")
    ## Filter for these
    sim_settings <- settings_df |> filter(family == ff & kappa == kk)

    wald_stats <- list()
    for (i in seq.int(nrow(sim_settings))) {
        sett <- sim_settings[i, ]
        sett_nm <- make_file_name_msglm(sett, "oracle_sgd_comp")
        load(file = file.path(outputs_path, sett_nm))

        head(stats$fits[["sgd_soft"]][,1:8])
        soft_nm <- with(sett, paste0(family, "-soft", "-", n))
        hard_nm <- with(sett, paste0(family, "-hard", "-", n))
        wald_stats[[soft_nm]] <- with(stats,
                                      t(apply(X = stats$fits[["sgd_soft"]],
                                              MARGIN = 1,
                                              FUN = function(x) x - theta_star)) / stats$fits[["ses"]])[,1:sett$s]

        wald_stats[[hard_nm]] <- with(stats,
                                      t(apply(X = stats$fits[["sgd_hard"]],
                                              MARGIN = 1,
                                              FUN = function(x) x - theta_star)) / stats$fits[["ses"]])[,1:sett$s]
    }

    fig_file <- file.path(outputs_path,
                          with(sim_settings[1,],
                               paste0("oracle_sgd_comp_grouped",
                                      "_family", family,
                                      "_p", p,
                                      "_s", s,
                                      ".pdf")))

    ss <- 5      # Only show first 5 variables
    n_mthds <- 2 # soft and hard
    n_vals <- sim_settings$n
    if (save_fig) pdf(fig_file, width = 10, height = 12/sqrt(2))
    par(mfcol = c(length(n_vals) * n_mthds, ss + 1), mar = c(2,2,2,2), oma = c(3,1,1,1))
    colsT <- rainbow_hcl(n_mthds, alpha = 0.7)
    colsF <- rainbow_hcl(n_mthds, alpha = 1)

    soft_prefix <- paste0(ff, "-soft")
    hard_prefix <- paste0(ff, "-hard")
    mthds_full <- paste0(rep(c("soft", "hard"), times = length(n_vals)),
                         "-",
                         rep(n_vals, each = n_mthds))

    for (jj in seq.int(wald_stats)) {
        plot(0, 0, type = "n", xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n")
        text(0, 0, mthds_full[jj], cex = 1.25)
    }
    nsim <- nrow(wald_stats[[1]])
    for (j in seq.int(ss)) {
        xlims <- range(sapply(X = wald_stats, FUN = function(x) range(x[,j]))) + c(-1, 1)
        for (ii in seq_along(n_vals)) {
            mn <- paste0("Variable ", j)
            if (ii > 1) {
                mn <- ""
                if (min(xlims) < -20)  xlims <- c(-7, 7)
            }
            nn <- n_vals[ii]
            makePlot(wald_stats[[paste0(soft_prefix, "-", n_vals[ii])]][,j],
                     ref0 = rep(FALSE, nsim),
                     ylim = c(0, 1),
                     xlim = xlims,
                     col = colsT[1],
                     border = colsF[1],
                     main = mn)
            makePlot(wald_stats[[paste0(hard_prefix, "-", n_vals[ii])]][,j],
                     ref0 = rep(FALSE, nsim),
                     ylim = c(0, 1),
                     xlim = xlims,
                     col = colsT[2],
                     border = colsF[2],
                     main = "")
        }
    }
    mtext("Distribution of non-zero coefficients",
          side = 1,
          outer = TRUE,
          line = 1)
    if (save_fig) dev.off()
}

## 3. Coverage of these confidence sets for one-pass SGD
## Using lmin or inversepower
tuning <- "lmin" # "ipower" or "lmin"
fn_prefix <- "cov_sgd"
if (tuning == "lmin") fn_prefix <- paste0(fn_prefix, "_", tuning)
load(file = file.path(outputs_path, paste(fn_prefix, "settings.rda", sep = "_")))

settings_df <- as.data.frame(settings) |>
mutate(sparsity = ifelse(s/p > 5/40, "dense", "sparse")) |>
mutate(sp = s / p) |>
mutate(method = ifelse(fixed, "hard-BIC", "soft-BIC")) |>
mutate(method_label = paste0(sigma_x, "-", method))

grouped_settings <- settings_df |>
group_by(family, fixed, sigma_x, n, sparsity, p, s, sp, true_param) |>
summarize()
sp_vec <- unique(grouped_settings$sp)
table_groups <- settings_df |> group_by(family, fixed) |> summarise()

for (t in seq.int(nrow(table_groups))) {
    table_sett <- table_groups[t, ]
    print(table_sett)
    table_df <- grouped_settings |> filter(family == table_sett$family, fixed == table_sett$fixed)
    sigma_xs <- unique(table_df$sigma_x)
    sparsity_print <- paste0(c(0, 0, sp_vec))
    sparsity_print[1:2] <- c("Covariance", "$n$")
    cat(make_row_tex(sparsity_print))
    for (i in seq_along(sigma_xs)) {
        cat("\\midrule\n")
        sigma_x_df <- table_df |> filter(sigma_x == sigma_xs[i])
        n_vals <- unique(sigma_x_df$n)
        for (j in seq_along(n_vals)) {
            row_print <- paste0(c(0, paste("$", n_vals[j], "$", sep = ""), 0, 0, 0, 0))
            row_print[1] <- ifelse(j == 1, as.character(sigma_xs[i]), "")
            row_df <- sigma_x_df |> filter(n == n_vals[j])
            for (k in seq_along(sp_vec)) {
                if (sp_vec[k] %in% row_df$sp) {
                    sett <-  row_df |> filter(sp == sp_vec[k])
                    sett_nm <- make_file_name(sett, fn_prefix)
                    load(file = file.path(outputs_path, sett_nm))
                    true_mod <- paste0(seq.int(sett$s), collapse = ",")
                    res <- lapply(X = stats$mdls, FUN = function(x) {
                        cs <- conf_set(t(x), 0.95)
                        list(is_in = true_mod %in% names(cs), card = length(cs))
                    })
                    cvrg <- mean(sapply(X = res, FUN = function(x) x$is_in))
                    av_card <- mean(sapply(X = res, FUN = function(x) x$card))
                    row_print[2 + k] <- paste0(round(100 * cvrg, 2), "\\%", " (", round(av_card,2), ")")
                } else {
                    row_print[2 + k] <- "--"
                }
            }
            cat(make_row_tex(row_print))
        }
    }
}



## 4. Confidence sets on the fly as an ad-hoc recommendation
load(file = file.path(outputs_path, "conf_sets.rda"))
## Hard thresholding
print_ci(mdls_hard)
## Soft thresholding
print_ci(mdls_soft)


