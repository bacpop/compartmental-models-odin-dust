# example for a Lotka Volterra Model in the literature: https://doi.org/10.1016/B978-008045405-4.00676-5


## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(R) <- R + n_R - n_RR - n_CR - n_C2R
update(C) <- C + n_C - n_CC - n_RC - n_C2C
update(C2) <- C2 + n_C2 - n_C2C2 - n_RC2 - n_CC2

## Total population size
#N <- S + I + R

## Individual probabilities of transition:
p_R <- 1 - exp(-r_R * dt) # growth part of logistic growth
p_RR <- 1 - exp(-r_R / K1 * R * dt) # decrease of growth rate when population size approaches capacity
p_CR <- 1 - exp(-r_R / K1 * alpha_CR * C * dt) # competition between species R and C
p_C2R <- 1 - exp(-r_R / K1 * alpha_C2R * C2 * dt) # competition between species R and C2

p_C <- 1 - exp(-r_C * dt) # growth part of logistic growth
p_CC <- 1 - exp(-r_C / K2 * C * dt) # decrease of growth rate when population size approaches capacity
p_RC <- 1 - exp(-r_C / K2 * alpha_RC * R * dt) # competition between species R and C
p_C2C <- 1 - exp(-r_C / K2 * alpha_C2C * C2 * dt) # competition between species C and C2

p_C2 <- 1 - exp(-r_C2 * dt) # growth part of logistic growth
p_C2C2 <- 1 - exp(-r_C2 / K3 * C2 * dt) # decrease of growth rate when population size approaches capacity
p_RC2 <- 1 - exp(-r_C2 / K3 * alpha_RC2 * R * dt) # competition between species R and C2
p_CC2 <- 1 - exp(-r_C2 / K3 * alpha_CC2 * C * dt) # competition between species C and C2

### should I divide C and R in p_RC and p_CR, respectively, by (R+C)? I am not sure.
### because for SI term that made sense, did it not?

## Draws from binomial distributions for numbers changing between
## compartments:
n_R <- rbinom(R, p_R) 
n_RR <- rbinom(R, p_RR)
n_CR <- rbinom(R, p_CR) 
n_C2R <- rbinom(R, p_C2R) 

n_C <- rbinom(C, p_C) 
n_CC <- rbinom(C, p_CC)
n_RC <- rbinom(C, p_RC)
n_C2C <- rbinom(C, p_C2C)

n_C2 <- rbinom(C2, p_C2) 
n_C2C2 <- rbinom(C2, p_C2C2)
n_RC2 <- rbinom(C2, p_RC2)
n_CC2 <- rbinom(C2, p_CC2)

## Initial states:
initial(R) <- R_ini
initial(C) <- C_ini
initial(C2) <- C2_ini

## User defined parameters - default in parentheses:
R_ini <- user(200)
C_ini <- user(50)
C2_ini <- user(50)
r_R <- user(0.2)
K1 <- user(500)
r_C <- user(0.15)
K2 <- user(500)
r_C2 <- user(0.15)
K3 <- user(500)

dim(alpha) <- c(3,3)
alpha[,] <- user()
alpha_CR <- alpha[2,1]
alpha_RC <- alpha[1,2]
alpha_C2R <- alpha[3,1]
alpha_C2C <- alpha[3,2]
alpha_RC2 <- alpha[1,3]
alpha_CC2 <- alpha[2,3]