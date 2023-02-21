# example for a Lotka Volterra Model in the literature: https://doi.org/10.1016/B978-008045405-4.00676-5


## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(R) <- R + n_R - n_RC
update(C) <- C - n_C + n_RC2

## Total population size
#N <- S + I + R

## Individual probabilities of transition:
p_R <- 1 - exp(-alpha * dt)
p_RC <- 1 - exp(-beta * C/(R+C)* dt) # prey killed by pred
p_RC2 <- 1 - exp(-beta * delta * C/(R+C)* dt) # pred benefitting from prey
p_C <- 1 - exp(-gamma * dt) # I to R

## Draws from binomial distributions for numbers changing between
## compartments:
n_R <- rbinom(R, p_R) # birth process for prey
n_RC <- rbinom(R, p_RC) # pred eating prey
n_RC2 <- rbinom(R, p_RC2) # pred benefitting from eating prey
n_C <- rbinom(C, p_C) # death process for pred




## Initial states:
initial(R) <- R_ini
initial(C) <- C_ini

## User defined parameters - default in parentheses:
R_ini <- user(200)
C_ini <- user(50)
alpha <- user(0.2)
beta <- user(0.15)
gamma <- user(0.1)
delta <- user(0.25)