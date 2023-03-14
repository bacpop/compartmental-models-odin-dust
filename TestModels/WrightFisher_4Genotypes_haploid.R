## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equation for population (assume constant size here, at first):

# number of individuals with Genotype A is determined by drawing from Bin(n, p)
# with n=pop_size and p = A/pop_size (relative abundance of Genotype A in the population)

p1 <-  (A/pop_size)
y1 <- rbinom(pop_size, p1)
p_tot <- 1- p1
p2 <- (B/pop_size)
y2 <- rbinom( pop_size-y1, (p2/p_tot))
p_tot2 <- p_tot - B/pop_size
p3 <- ((C/pop_size)/p_tot2)
y3 <- rbinom(pop_size-y1-y2, p3)

### currently issue with random draws
### by re-scaling p it can sometimes become (slightly) bigger than 1!
### not sure why this is not an issue for them:
### https://github.com/mrc-ide/dust/blob/master/inst/include/dust/random/binomial.hpp


#n_inf[2:4] <- rbinom(pop_size - sum(n_inf[1:(i - 1)]), p_inf[i])

update(A) <- y1
update(B) <- y2
update(C) <- y3
update(D) <- pop_size-y1-y2-y3

## Initial states:
initial(A) <- pop_size * A_ini # deterministic, user-based start value
initial(B) <- pop_size * B_ini # deterministic, user-based start value
initial(C) <- pop_size * C_ini # deterministic, user-based start value
initial(D) <- pop_size * D_ini # deterministic, user-based start value

## User defined parameters - default in parentheses:
A_ini <- user(.25) # initial frequency of Genotype A
B_ini <- user(.25) # initial frequency of Genotype B
C_ini <- user(.25) # initial frequency of Genotype C
D_ini <- user(.25) # initial frequency of Genotype D
pop_size <- user() # population size
