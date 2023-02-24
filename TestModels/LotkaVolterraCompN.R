# example for a Lotka Volterra Model in the literature: https://doi.org/10.1016/B978-008045405-4.00676-5


## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
update(C[]) <- C[i] + n_Ci[i] - sum(n_CjCi[i,]) ### something prob not right with the indexes yet

## Total population size
#N <- S + I + R

## Individual probabilities of transition:
p_Ci[] <- 1 - exp(-r_C[i] * dt) # growth part of logistic growth
p_CjCi[,]<- 1 - exp(-r_C[i] / K[i] * sum(alpha_int[i,]) * dt) 


alpha_int[,] <- alpha[i, j] * C[j]

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
species_no <- user()
dim(C) <- species_no
dim(C_ini) <- species_no
dim(r_C) <- species_no
dim(K) <- species_no
dim(alpha) <- c(species_no,species_no)
dim(alpha_int) <- c(species_no,species_no)
dim(p_Ci) <- species_no
dim(p_CjCi) <- c(species_no,species_no)
dim(n_Ci) <- species_no
dim(n_CjCi) <- c(species_no,species_no)