# example for a Lotka Volterra Model in the literature: https://doi.org/10.1016/B978-008045405-4.00676-5


## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(R) <- R + n_R - n_RR - n_RC
update(C) <- C + n_C - n_CC - n_CR

## Total population size
#N <- S + I + R

## Individual probabilities of transition:
p_R <- 1 - exp(-r_R) # growth part of logistic growth
p_RR <- 1 - exp(-r_R / K1 * R ) # decrease of growth rate when population size approaches capacity
p_RC <- 1 - exp(-r_R / K1 * alpha_RC * C/(R+C)) #* dt) # competition between species R and C

p_C <- 1 - exp(-r_C) # growth part of logistic growth
p_CC <- 1 - exp(-r_C / K2 * C ) # decrease of growth rate when population size approaches capacity
p_CR <- 1 - exp(-r_R / K2 * alpha_CR * R/(R+C)) #* dt) # competition between species R and C

## Draws from binomial distributions for numbers changing between
## compartments:
n_R <- rbinom(R, p_R) 
n_RR <- rbinom(R, p_RR)
n_RC <- rbinom(R, p_RC) 

n_C <- rbinom(C, p_C) 
n_CC <- rbinom(C, p_CC)
n_CR <- rbinom(C, p_CR)


## Initial states:
initial(R) <- R_ini
initial(C) <- C_ini

## User defined parameters - default in parentheses:
R_ini <- user(200)
C_ini <- user(50)
r_R <- user(0.2)
K1 <- user(500)
r_C <- user(0.15)
K2 <- user(500)
alpha_RC <- user(.2)
alpha_CR <- user(.3)