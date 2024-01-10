source("simulations.R")
N = 250
p = 40
s = 5
init_control = default_init(p)
init_control$gamma.method = "ipower"
out = parallel_sim(s, p, N, nreps=20, true_param = "hht", model="binomial", sigma_x="id", init_control=init_control)



## Check oracle distribution of estimates, will need to restimate models after selection


## Perform model selection with each new iteration.
## What needs to be done:
## 1. Perform SGD iteratively: will need to compute gamma^* iteratively as well.
## 2. Compute pivots
## 3. Compute and store model
## 4. And repeat...
## 5. Whilst doing this build a confidence set
