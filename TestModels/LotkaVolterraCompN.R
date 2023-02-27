# example for a Lotka Volterra Model in the literature: https://doi.org/10.1016/B978-008045405-4.00676-5

## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
# C is the vector of model compartments. Its dimension is determined by the input of the user (species_no)
update(C[]) <- C[i] - n_death[i] + n_birth[i]

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
species_no <- user() # number of species

dim(p_birth) <- species_no
dim(p_death) <- species_no
dim(n_death) <- species_no
dim(n_birth) <- species_no
dim(alpha) <- c(species_no,species_no)
dim(alpha_pop) <- c(species_no, species_no)

dim(C) <- species_no
dim(C_ini) <- species_no
dim(r_C) <- species_no
dim(K) <- species_no



### alternative version with poisson distribution:

#p_birth[] <- r_C[i] * C[i]
#p_death[] <- r_C[i]/K[i] * C[i] * sum(alpha_pop[i,])

#n_death[] <- rpois(p_death[i] * dt)
#n_birth[] <- rpois(p_birth[i] * dt)
