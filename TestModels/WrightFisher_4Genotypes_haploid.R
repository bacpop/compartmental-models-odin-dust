## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equation for population (assume constant size here, at first):


p1 <-  A/pop_size
p2 <- B/pop_size
p3 <- C/pop_size
p4 <- D/pop_size


# number of individuals with Genotype A is determined by drawing from Bin(n, p)
# with n=pop_size and p = A/pop_size (relative abundance of Genotype A in the population)
y1 <- if (p1/(p1+p2+p3+p4) < 1)
  rbinom(pop_size, p1/(p1+p2+p3+p4)) else pop_size

# The n and p for the remaining binomial draws have to be adapted
# because we only want to draw (n-x), i.e. those individuals that are not yet determined to be of Genotype A (and then B, C)
# and we need to adapt the probabilities accordingly to keep the population size the same

y2 <- if (p2/(p2+p3+p4) < 1)
  rbinom(pop_size - y1, p2/(p2+p3+p4)) else pop_size - y1

y3 <- if (p3/(p3+p4) < 1)
  rbinom(pop_size -y1 -y2, p3/(p3+p4)) else pop_size -y1 - y2

y4 <- pop_size - y1 - y2 - y3

update(A) <- y1
update(B) <- y2
update(C) <- y3
update(D) <- y4

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
