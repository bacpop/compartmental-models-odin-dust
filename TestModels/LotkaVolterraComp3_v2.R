# example for a Lotka Volterra Model in the literature: https://doi.org/10.1016/B978-008045405-4.00676-5


## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(C[]) <- C[i] + n_Ci[i] - sum(n_CjCi[i,]) ### something prob not right with the indexes yet
#update(C) <- C + n_C - n_CC - n_RC - n_C2C
#update(C2) <- C2 + n_C2 - n_C2C2 - n_RC2 - n_CC2

## Total population size
#N <- S + I + R

## Individual probabilities of transition:
p_Ci[] <- 1 - exp(-r_C[i] * dt) # growth part of logistic growth
p_CjCi[,]<- 1 - exp(-r_C[i] / K[i] * alpha[j,i] * C[j] * dt) # S to I
#### BE CAREFUL! THIS SHOULD BE J in alpha and C but did not work!


### should I divide C and R in p_RC and p_CR, respectively, by (R+C)? I am not sure.
### because for SI term that made sense, did it not?

## Draws from binomial distributions for numbers changing between
## compartments:
n_Ci[] <- rbinom(C[i], p_Ci[i])
n_CjCi[,] <- rbinom(C[i], sum(p_CjCi[,i])) ### something prob not right with the indexes yet


## Initial states:
initial(C[]) <- C_ini[i]

## User defined parameters - default in parentheses:

C_ini[] <- user()
K[] <- user()
r_C[] <- user()
alpha[,] <- user()

#dimensions:
dim(C) <- 3
dim(C_ini) <- 3
dim(r_C) <- 3
dim(K) <- 3
dim(alpha) <- c(3,3)
dim(p_Ci) <- 3
dim(p_CjCi) <- c(3,3)
dim(n_Ci) <- 3
dim(n_CjCi) <- c(3,3)