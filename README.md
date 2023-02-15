# "Plus/minus the learning rate": Easy and Scalable Statisitcal Inference with SGD.
To appear in AISTATS'23.

A short example on how to run simulations to evaluate our confidence intervals vs MLE as implemented by the glm() function in R. 
```
source("simulations.R")
N = 1e4
p = 100
init_control = default_init(p)
init_control$gamma.method = "ipower"
out = parallel_sim(p, N, nreps=100, model="gaussian", sigma_x="id", init_control=init_control)
```
Output contains element-wise coverage, average coverage, and average interval length for our method and MLE.
- `p` controls the dimension.
- `N` controls the number of samples.
- `nreps` controls number of separate confidence intervals to generate, with newly generated data for each confidence interval.
- Set `model` to {`gaussian`, `binomial`, or `poisson`}.
- Set covariance matrix `sigma_x` to values {`id`, `equicor`, `toeplitz`, `ill_cond`}.
- To choose gamma selection method, set `init_control$gamma.method` to one of values {`heuristic`, `bound`, `lmin`, `ipower`}.
