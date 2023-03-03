## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equation for population (assume constant size here, at first):

# number of individuals with B allele is determined by drawing from Bin(n, p)
# with n=pop_size and p = B/pop_size (relative abundance of B alleles in the population)
update(B) <- rbinom(pop_size, B/pop_size) 
update(A) <- pop_size - B # all individuals that don't carry the B allele must carry the A allele


## Initial states:
initial(B) <- pop_size * B_ini # determinstic, user-based start value
initial(A) <- pop_size - (pop_size * B_ini) # determinstic, user-based start value

## User defined parameters - default in parentheses:
B_ini <- user() # initial frequency of B alleles
pop_size <- user() # population size
gen <- user() # number of generations to be simulated