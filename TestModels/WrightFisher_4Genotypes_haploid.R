## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equation for population (assume constant size here, at first):

# number of individuals with Genotype A is determined by drawing from Bin(n, p)
# with n=pop_size and p = A/pop_size (relative abundance of Genotype A in the population)

p1 <-  A/pop_size
p2 <- B/pop_size
p3 <- C/pop_size
p4 <- D/pop_size
p_tot <- p1 + p2 + p3 + p4
y1 <- rbinom(pop_size, p1/p_tot)
p_tot1 <- p_tot - p1
y2 <- rbinom(pop_size-y1, (p2/(p2+p3+p4)))
p_tot2 <- p_tot - p1 - p2
y3 <- rbinom(pop_size-y1-y2, p3/(p3+p4))
y4 <- pop_size - y1 - y2 - y3
### currently issue with random draws
### by re-scaling p it can sometimes become (slightly) bigger than 1!
### not sure why this is not an issue for them:
### https://github.com/mrc-ide/dust/blob/master/inst/include/dust/random/multinomial.hpp


#n_inf[2:4] <- rbinom(pop_size - sum(n_inf[1:(i - 1)]), p_inf[i])

update(A) <- y1
update(B) <- y2
update(C) <- y3
update(D) <- y4
#update(p1) <- y1/pop_size
#update(p2) <- y2/pop_size
#update(p3) <- y3/pop_size
#update(p4) <- y4/pop_size

## Initial states:
initial(A) <- pop_size * A_ini # deterministic, user-based start value
initial(B) <- pop_size * B_ini # deterministic, user-based start value
initial(C) <- pop_size * C_ini # deterministic, user-based start value
initial(D) <- pop_size * D_ini # deterministic, user-based start value
#initial(p1) <- A_ini
#initial(p2) <- B_ini
#initial(p3) <- C_ini
#initial(p4) <- D_ini


## User defined parameters - default in parentheses:
A_ini <- user(.25) # initial frequency of Genotype A
B_ini <- user(.25) # initial frequency of Genotype B
C_ini <- user(.25) # initial frequency of Genotype C
D_ini <- user(.25) # initial frequency of Genotype D
pop_size <- user() # population size
