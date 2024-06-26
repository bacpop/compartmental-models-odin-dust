## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(S) <- S - n_SE
update(E) <- E + n_SE - n_EI
update(I) <- I + n_EI - n_IR
update(R) <- R + n_IR

## Individual probabilities of transition:
p_SE <- 1 - exp(-beta * I / N * dt) # S to E
p_EI <- 1 - exp(-delta * dt) # E to I
p_IR <- 1 - exp(-gamma * dt) # I to R

## Draws from binomial distributions for numbers changing between
## compartments:
n_IR <- rbinom(I, p_IR)
n_EI <- rbinom(E, p_EI)
n_SE <- rbinom(S, p_SE)

## Total population size
N <- S + E + I + R

## Initial states:
initial(S) <- S_ini
initial(E) <- E_ini
initial(I) <- I_ini
initial(R) <- 0

## User defined parameters - default in parentheses:
S_ini <- user(2000)
E_ini <- user(5)
I_ini <- user(10)
beta <- user(0.2)
delta <- user(0.1)
gamma <- user(0.1)