## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(S) <- S - n_SI + n_IS
update(I) <- I + n_SI - n_IS

## Individual probabilities of transition:
p_SI <- 1 - exp(-beta * I / N * dt) # S to I
p_IS <- 1 - exp(-delta * dt) # I to S

## Draws from binomial distributions for numbers changing between
## compartments:
n_IS <- rbinom(I, p_IS)
n_SI <- rbinom(S, p_SI)

## Total population size
N <- S + I

## Initial states:
initial(S) <- S_ini
initial(I) <- I_ini

## User defined parameters - default in parentheses:
S_ini <- user(1000)
I_ini <- user(10)
beta <- user(0.2)
delta <- user(0.1)