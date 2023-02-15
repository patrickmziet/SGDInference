# "Plus/minus the learning rate": Easy and Scalable Statisitcal Inference with SGD.
To appear in AISTATS'23.

A short example on how to run simulations to evaluate our confidence intervals vs MLE. 
```
source("simulations.R")
N = 1e4
p = 100
init_control = default_init(p)
init_control$gamma.method = "ipow"
out = parallel_sim(p, N, nreps=100, model="gaussian", sigma_x="id", init_control=init_control)
```
- `p` controls the dimension.
- `N` controls the number of samples.
- Set `model` to {`gaussian`, `binomial`, or `poisson`}.
- Set covariance matrix `sigma_x` to values {`id`, `equicor`, `toeplitz`, `ill_cond`}.
- To choose gamma selection method, set `init_control$gamma.method` to one of values {`heuristic`, `bound`, `lmin`, `ipower`}.
