## Run model selection and estimation iteratively

## Load functions
source("simulations.R")
source("exact_inference.R")
source("iterative_functions.R")

## Hyperameters
model_name <- "binomial"
true_param <- "hht"
nn <- 2000
pp <- 10
ss <- 5
B <- 50 # batch size
lvl <- 0.95 # confidence set


res <- run_sim(model_name = model_name,
               true_param = true_param,
               nn = nn,
               pp = pp,
               ss = ss,
               B = B,
               lvl = lvl)
str(res)
print_ci(res$mdls)


