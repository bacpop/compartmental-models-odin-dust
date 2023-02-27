# example for a Lotka Volterra Model in the literature: https://doi.org/10.1016/B978-008045405-4.00676-5


## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(C[]) <- C[i] - n_death[i] + n_birth[i]
#update(C) <- C + n_C - n_CC - n_RC - n_C2C
#update(C2) <- C2 + n_C2 - n_C2C2 - n_RC2 - n_CC2

## Total population size
#N <- S + I + R

## Individual probabilities of transition:
alpha_pop[,] <- alpha[i, j] * C[j] #calculates the current influences of species, based on competition matrix and current population sizes
# note to myself: this is a dot product, two for-loops. This is not standard matrix * vector multiplication (otherwise the result would be a vector, not a matrix)

p_birth[] <- 1 - exp(-r_C[i] * dt) # probabilities for birth process
p_death[] <- 1 - exp(-r_C[i] / K[i] * sum(alpha_pop[i,]) * dt) # probabilities for death processes, including the death process that is part of the logistic growth

## Draws from binomial distributions for births and deaths in each population (=compartment):
n_birth[] <- rbinom(C[i],p_birth[i])
n_death[] <- rbinom(C[i],p_death[i]) 


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
dim(alpha_pop) <- c(3,3)
dim(p_birth) <- 3
dim(p_death) <- 3
dim(n_birth) <- 3
dim(n_death) <- 3