source("simulations.R")
N = 1e4
p = 10
init_control = default_init(p)
init_control$gamma.method = "ipower"
out = parallel_sim(p, N, nreps=10, model="gaussian", sigma_x="id", init_control=init_control)
