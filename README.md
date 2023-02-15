# SGDInference
```
source("simulations.R")
N = 1e4
p = 100
init_control = default_init(p)
init_control$gamma.method = "ipow"
out = parallel_sim(p, N, nreps=100, init_control=init_control)
```
