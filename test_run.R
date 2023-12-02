source("simulations.R")
N = 1e4
p = 10
s = 5
init_control = default_init(p)
init_control$gamma.method = "ipower"
out = parallel_sim(s, p, N, nreps=10, true_param = "hht", model="binomial", sigma_x="id", init_control=init_control)

