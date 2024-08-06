
# Compare multinom and poisson for likelihood

set.seed(73)

sim_data <- rpois(10, lambda = 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3))
sim_data2 <- rpois(10, lambda = 15000 * c(0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.15, 0.3))
sim_data3 <- rpois(10, lambda = 15000 * c(0.05, 0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.35))
sim_data4 <- rpois(10, lambda = 15000 * c(0.025, 0.05, 0.05, 0.075, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3))
sim_data5 <- rpois(10, lambda = 15000 * c(0.05, 0.05, 0.05, 0.05, 0.075, 0.1, 0.1, 0.1, 0.125, 0.3))

sim_data6 <- rmultinom(n = 1, size = 15000, prob = c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3))

# sim_data (probs the same between data and likelihood)
dmultinom(x = sim_data, prob = c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3))
# 1.916216e-20
prod(dpois(x = sim_data, lambda = 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 5.004115e-23

# sim_data2 changed size of mid-size clusters
dmultinom(x = sim_data2, prob = c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3))
# 4.136072e-189
prod(dpois(x = sim_data2, lambda = 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 1.08595e-191

# sim_data3 changed size of large clusters
dmultinom(x = sim_data3, prob = c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3))
# 3.997718e-132
prod(dpois(x = sim_data3, lambda = 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 1.054437e-134

# sim_data4 changed size of small clusters
dmultinom(x = sim_data4, prob = c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3))
# 2.981278e-108
prod(dpois(x = sim_data4, lambda = 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 9.375368e-111

# sim_data5 changed size of mid-size clusters
dmultinom(x = sim_data5, prob = c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3))
# 3.884976e-68
prod(dpois(x = sim_data5, lambda = 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 1.247864e-70

# test whether simulating with multinom makes a difference (does not seem so)
dmultinom(x = sim_data6, prob = c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3))
# 8.885988e-20
prod(dpois(x = sim_data6, lambda = 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 2.894461e-22

library(philentropy)
JSD(x = rbind(sim_data, 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 2.327953 
JSD(x = rbind(sim_data2, 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 148.9348 
JSD(x = rbind(sim_data3, 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 102.3734 
JSD(x = rbind(sim_data4, 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 77.56872 
JSD(x = rbind(sim_data5, 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 42.08756 
JSD(x = rbind(t(sim_data6), 15000 * c(0.05, 0.05, 0.05, 0.05, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3)))
# 1.663402

# for all, jsd, multinom, and Poisson, the it's 1, 5, 4, 3, 2 (most likely to least likely). So they don't seem to make much difference.
# multinom has the advantage of producing an actually constant population size
# JSD is slightly better for the multinomial draw --> maybe that means that multinomial is better?