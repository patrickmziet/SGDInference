library("colorspace")
library("dplyr")
library("ggplot2")
library("patchwork")
source("support_functions.R")

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


## 2. Oracle property
load(file = file.path(outputs_path, "oracle_sgd_settings.rda"))
## Plots for some settings
colsT <- rainbow_hcl(1, alpha = 0.7)
colsF <- rainbow_hcl(1, alpha = 1)
ls()
settings

settings_df <- as.data.frame(settings) |>
mutate(sparsity = ifelse(s/p > 5/40, "dense", "sparse")) |>
mutate(method = ifelse(fixed, "hard-BIC", "soft-BIC")) |>
mutate(method_label = paste0(sigma_x, "-", method))

## Group data you want to plot together
grouped_settings <- settings_df |>
group_by(sparsity, family, fixed, n, p, s, sigma_x) |>
summarize()

print(grouped_settings, n=42)

## Group figures you want together, in this case we will have two in each figure group
figure_groups <- grouped_settings |>
group_by(family, sparsity) |>
summarize()

save_fig <- TRUE
ff <- "binomial"
hh <- TRUE
pp <- 40
ss <- 5
nn <- 500
sigma_xs <- c("id", "toeplitz", "equicor")
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
    text(0, 0, mthds_full[k])
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


## Tables for other settings

## 3. Coverage of these confidence sets for one-pass SGD

## 4. Confidence sets on the fly as an ad-hoc recommendation



## load(file = file.path(outputs_path, sett_nm))
## ## TODO: You are here, need to fix the differencing with matrix and vector.

## wald_stats <- with(stats, ((coefs - theta_star) / ses)[,1:sett$s])

## wald_stats <- with(stats, (t(apply(X = coefs, MARGIN = 1, FUN = function(x) x - theta_star)) / ses)[,1:sett$s])

## stats$ses
## aa <- apply(X=stats$coefs, MARGIN = 1, FUN = function(x) x - stats$theta_star)
## apply(X=stats$coefs, MARGIN = 2, FUN = function(x) sum(x))


## par(mfrow = c(3, 2))
## for (j in seq.int(5)) {
##     xlims <- range(wald_stats[, j]) + c(-1, 1)
##     makePlot(wald_stats[, j],
##              ref0 = rep(FALSE, nrow(wald_stats)),
##              ylim = c(0, 1),
##              xlim = xlims,
##              col = colsT[1],
##              border = colsF[1],
##              main = paste0("Variable ", j))
## }

