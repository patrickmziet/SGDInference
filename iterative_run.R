## Run model selection and estimation iteratively
## Load functions
source("simulations.R")
source("exact_inference.R")
source("iterative_functions.R")

## Hyperameters
model_name <- "binomial"
true_param <- "hht"
nn <- 2000
pp <- 40
ss <- 25
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
               lvl = lvl)

print_ci(res$mdls)

str(res)



