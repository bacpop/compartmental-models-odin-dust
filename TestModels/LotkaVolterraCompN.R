# example for a Lotka Volterra Model in the literature: https://doi.org/10.1016/B978-008045405-4.00676-5

## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equations for transitions between compartments:
# C is the vector of model compartments. Its dimension is determined by the input of the user (species_no)
#update(C[]) <- C[i] - n_change[i]
update(C[]) <- C[i] + n_birth[i] - n_death[i] + n_help[i]

## Individual probabilities of transition:
alpha_pos_pop[,] <- alpha_pos[i, j] * C[j] #calculates the current influences of species, based on competition matrix and current population sizes
# note to myself: this is a dot product, two for-loops. This is not standard matrix * vector multiplication (otherwise the result would be a vector, not a matrix)
alpha_neg_pop[,] <- alpha_neg[i,j] * C[j]



p_birth[] <- r_C[i] * C[i] # linear birth rate (for population size far away from K)
#p_death[] <- r_C[i]/(K[i]*K[i]) * C[i] * sum(alpha_neg_pop[i,])
p_death[] <- 1 - exp(-r_C[i] / K[i] * sum(alpha_neg_pop[i,]) * dt) # negative influence of other species and logistic part of birth rate
p_help[] <- r_C[i]/(K[i]*K[i]) * C[i] * sum(alpha_pos_pop[i,]) # positive influence of other species 

###### version 1: separating birth and death processes:
#n_death[] <- rpois(p_death[i] * dt)
#n_death[] <- rbinom(C[i],p_death[i]) # this is binomial because for the deaths we draw from the population and not from an infinite pool
#n_birth[] <- rpois(p_birth[i] * dt)

##### version 2: now joint draws for births and deaths
#p_change[] <- p_birth[i] - p_death[i]
#n_change[] <- rpois(p_change[i] * dt)
#dim(p_change) <- species_no
#dim(n_change) <- species_no

###### version 3: separating birth and death processes and positive influence:
#n_death[] <- rpois(p_death[i] * dt)
n_death[] <- rbinom(C[i],p_death[i]) # this is binomial because for the deaths we draw from the population and not from an infinite pool
n_birth[] <- rpois(p_birth[i] * dt)
n_help[] <- rpois(p_help[i] * dt)


## Initial states:
initial(C[]) <- C_ini[i]

## User defined parameters - default in parentheses:

C_ini[] <- user()
K[] <- user()
r_C[] <- user()
alpha_pos[,] <- user()
alpha_neg[,] <- user()

#dimensions:
species_no <- user() # number of species

dim(p_birth) <- species_no
dim(p_death) <- species_no
dim(p_help) <- species_no
dim(n_help) <- species_no

dim(n_birth) <- species_no
dim(n_death) <- species_no
dim(alpha_pos) <- c(species_no,species_no)
dim(alpha_pos_pop) <- c(species_no, species_no)
dim(alpha_neg) <- c(species_no, species_no)
dim(alpha_neg_pop) <- c(species_no, species_no)

dim(C) <- species_no
dim(C_ini) <- species_no
dim(r_C) <- species_no
dim(K) <- species_no



### alternative version with binomial distribution:

#p_birth[] <- 1 - exp(-r_C[i] * dt) # probabilities for birth process
#p_death[] <- 1 - exp(-r_C[i] / K[i] * sum(alpha_pop[i,]) * dt) # probabilities for death processes, including the death process that is part of the logistic growth


## Draws from binomial distributions for births and deaths in each population (=compartment):
#n_birth[] <- rbinom(C[i],p_birth[i])
#n_death[] <- rbinom(C[i],p_death[i]) 

#p_birth[] <- r_C[i] * C[i]
#p_death[] <- r_C[i]/K[i] * C[i] * sum(alpha_pop[i,])

#n_death[] <- rpois(p_death[i] * dt)
#n_birth[] <- rpois(p_birth[i] * dt)
