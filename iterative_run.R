## Run model selection and estimation iteratively

## Load functions
source("simulations.R")
source("exact_inference.R")
source("iterative_functions.R")

## Hyperameters
nn <- 5000
pp <- 10
ss <- 5
B <- 50 # batch size
lvl <- 0.95 # confidence set

## Run simulation
## Generate data
data <- gen_data(model_name = "binomial",
                 N = nn,
                 p = pp,
                 s = ss,
                 true_param = "hht")
str(data)
sq <- seq.int(nn)
chnks <- split(sq, ceiling(seq_along(sq) / B))

init_control <- list()
init_control$gamma.method <- "ipower"
init_control_def <- default_init(pp)
for(n in names(init_control)) init_control_def[[n]] <- init_control[[n]]
init_control <- init_control_def

pivots <- matrix(0, nrow = length(chnks), ncol = pp)
Jhat <- matrix(FALSE, nrow = length(chnks), ncol = pp)
## mdls <- matrix(".", nrow = length(chnks), ncol = 1)
mdls <- NA
data_cp <- data

for (i in seq.int(chnks)) {
    print(paste0("Batch ", i, "\n"))
    ## reference data
    ref <- seq(1, max(chnks[[i]]))
    data_cp$X <- data$X[ref,]
    data_cp$Y <- data$Y[ref]

    N <- nrow(data_cp$X)
    p <- ncol(data_cp$X)

    ## Compute gammastar
    sgd_control <- calculate_gammaStar_wrapper(data_cp, init_control=init_control)
    gamma_star <- sgd_control$gamma

    ## Very inefficient
    theta_hat <- isgd_0(sgd_control$theta0, data_cp, gamma=gamma_star, 
                        npass=1, use_permutation = FALSE)
    psi <- get_model_psi(data_cp, theta_hat)

    V <- psi * gamma_star * rep(1, p) / N
    ## Compute pivots
    pivots[i,] <- abs(theta_hat) / sqrt(V)
    thresholds <- optimise_threshold(betas = theta_hat,
                                       stderrs = sqrt(V),
                                       n = N,
                                       p = p,
                                       fixed=TRUE) 
    Jhat[i, ] <- pivots[i, ] > thresholds  # Estimated model
    mdls[i] <- paste0("(",paste0((1:p)[Jhat[i, ]], collapse = ","),")")
    ## Print diagnostics as iterations go on
    print(mdls[i])
    print_ci(mdls)
}


print_ci(mm=mdls)


oo <- order(prps, decreasing = TRUE)
print(prps[oo])

str(prps)

table(mdls)
