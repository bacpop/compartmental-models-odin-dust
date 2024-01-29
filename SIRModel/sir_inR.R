##https://rpubs.com/choisy/sir
#install.packages("deSolve")
library(deSolve) # using the "ode" function
sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dS <- -beta * I * S / (S + I + R)
    dI <-  beta * I * S / (S + I + R) - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}

parameters_values <- c(
  beta  = 0.05, # infectious contact rate
  gamma = 0.1    # recovery rate
)

initial_values <- c(
  S = 2000,  # number of susceptibles at time = 0
  I =   10,  # number of infectious at time = 0
  R =   0   # number of recovered at time = 0
)

time_values <- seq(0, 150)

#install.packages("deSolve")
library(deSolve)
sir_values <- ode(
  y = initial_values,
  times = time_values,
  func = sir_equations,
  parms = parameters_values
)


cols <- c(S = "#000000", I = "#E69F00", R = "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #8 colorblind friendly colors
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(sir_values[, 1], sir_values[, -1], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = cols, lty = 1)
legend("topright", lwd = 1, col = cols, legend = c("S", "I", "R"), bty = "n")


cols <- c(I = "#E69F00", R = "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #8 colorblind friendly colors
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(sir_values[, 1], sir_values[, c(3,4)], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = cols, lty = 1)
legend("topright", lwd = 1, col = cols, legend = c("I", "R"), bty = "n")

#install.packages("phaseR")
library(phaseR)

sir_equations_2d <- function(time, variables, parameters) {
  S <- variables[1]
  I <- variables[2]
  beta <- parameters[1]
  gamma <- parameters[2]
  dS <- - beta * I * S /10
  dI <-  beta * I * S /10 - gamma * I
  return(list(c(dS, dI)))
}


init <- c(999, 1)
parms <- c(beta = 0.002, gamma = 0.1)
steps <- seq(0, 150)

sir_values_2d <- ode(
  y = init,
  times = steps,
  func = sir_equations_2d,
  parms = parms
)

par(mfrow = c(1, 1)) # reset plotting parameters
plot(sir_values_2d[,2], sir_values_2d[,3],
     type = "l", xlab = "S", ylab="I",
     col = "magenta", main = "SI phase plane for SIR model")

plot(sir_values_2d[, c(2, 3)], type = "l", col = "magenta",
     xlim = c(0, 1000), ylim=c(0, 200), xlab = "S", ylab = "I")
flow <- flowField(sir_equations_2d, xlim = c(0, 1000), ylim = c(0, 200),
                  parameters = parms)


