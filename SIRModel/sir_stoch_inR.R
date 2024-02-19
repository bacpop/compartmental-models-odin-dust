stoch_sir_equations <- function(variables, parameters) {
  with(as.list(c(variables, parameters)), {
    n_SI <- rbinom(1,S, 1-exp(-beta * I / (S+I+R)))
    n_IR <- rbinom(1,I, 1 - exp(-gamma))
    dS <- S - n_SI
    dI <-  I + n_SI - n_IR
    dR <-  R + n_IR
    return((c(dS, dI, dR)))
  })
}

parameters_values <- c(
  beta  = 0.25, # infectious contact rate
  gamma = 0.1    # recovery rate
)

initial_values <- c(
  S = 2000,  # number of susceptibles at time = 0
  I =   10,  # number of infectious at time = 0
  R =   0   # number of recovered at time = 0
)

time_values <- seq(0, 150)

stoch_sir_values <- matrix(data = NA, nrow = time_values[length(time_values)]+1, ncol = 4)
colnames(stoch_sir_values) <- c("time", "S", "I", "R")
stoch_sir_values[1,] <- c(0, 2000, 10, 0)
for (i in 1:150) {
  stoch_sir_values[i+1,] <- c(i, stoch_sir_equations(stoch_sir_values[i,c(2,3,4)],parameters_values))

}

cols <- c(S = "#000000", I = "#E69F00", R = "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #8 colorblind friendly colors
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(stoch_sir_values[, 1], stoch_sir_values[, -1], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = cols, lty = 1)
legend("topright", lwd = 1, col = cols, legend = c("S", "I", "R"), bty = "n")
# plot densities
n <- 2000
k <- seq(0, 20, by = 1)
plot (k, dbinom(k, n, 1-exp(-0.25 * 10 / (2010)), log = FALSE), type = "l", ylab = "density")
#lines(k, (dbinom(k, n, 1-exp(-0.25 * 10 / (2010)))), col = "red", lwd = 2)

n <- 10
k <- seq(0, 10, by = 1)
plot (k, dbinom(k, n, 1 - exp(-0.1), log = FALSE), type = "l", ylab = "density")
