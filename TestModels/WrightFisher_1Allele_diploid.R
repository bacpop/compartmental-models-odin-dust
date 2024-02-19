## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equation for population (assume constant size here, at first):

# number of individuals with AA and AB are determined by drawing from Multinom(n, [p1,p2])
# with n=pop_size and p1 = (p + 1/2q) + (1/2 * (p + 1/2 q)) and p2 = (r + 1/2 q) + (1/2 * (p+q+r)) + (p + 1/2 q)
# the number of BB individuals is determined by pop_size - number of AA - number of AB

############## 
### The code does not throw any errors but looking at the development of the population (BB always seems to survive), there is most likely an issue with the random draws.
### I will not work this out because the diploid model is not interesting for us at the moment. 

probAA_AA <- (AA / pop_size + 1/2 * (AB / pop_size)) * AA/pop_size
probAB_AA <- (1/2 * (AA / pop_size + 1/2 * (AB / pop_size))) * AB/pop_size
#nextGen[1] <- rbinom(AA, probAA_AA)  
#nextGen[2] <- rbinom(AB, probAB_AA)
#nextGen[1] <- rbinom(pop_size, probAA_AA)  
#nextGen[2] <- rbinom(pop_size, probAB_AA)
nextGen[1] <- rbinom(pop_size, probAA_AA + probAB_AA)
nextGen[2] <- 0

probAA_AB <- (BB / pop_size + 1/2 * (AB / pop_size)) * AA/pop_size
probAB_AB <- (1/2 * (AA + AB + BB) / pop_size)* AB/pop_size
probBB_AB <- (AA / pop_size + 1/2 * (AB / pop_size))* BB/pop_size
#nextGen[3] <- rbinom(AA - nextGen[1], probAA_AB) + rbinom(AB - nextGen[2], probAB_AB)+ rbinom(BB, probBB_AB)
#nextGen[3] <- rbinom(pop_size , probAA_AB)
#nextGen[4] <- rbinom(pop_size , probAB_AB)
#nextGen[5] <- rbinom(pop_size , probBB_AB)
nextGen[3] <- rbinom(pop_size -nextGen[1], probAA_AB+probAB_AB+probBB_AB)
nextGen[4] <- 0
nextGen[5] <- 0

#p[1] <- .2
#p[2] <- .8
#y[1] <- rbinom(pop_size, p[1])
#y[2] <- rbinom(pop_size - y[1], p[2])

#dim(p) <- 2
#dim(y) <- 2


#multiProb[1] <- 0.2
#multiProb[2] <- 0.8
#nextGen[] <- rmultinom(10,multiProb)


#nextGen[] <- rmultinom(pop_size, multiProb)

update(AA) <- nextGen[1] + nextGen[2]
update(AB) <- nextGen[3] + nextGen[4] + nextGen[5]
update(BB) <- pop_size - AA - AB # all individuals that are not AA or AB must be BB


## Initial states:
initial(AA) <- Pop_ini[1] # deterministic, user-based start value
initial(AB) <- Pop_ini[2] # deterministic, user-based start value
initial(BB) <- Pop_ini[3]




## User defined parameters - default in parentheses:
dim(Pop_ini) <- 3
Pop_ini[] <- user() 
dim(nextGen) <- 5
#dim(multiProb) <- 2
#dim(nextGen) <- 2
pop_size <- user(100) # population size
#gen <- user() # number of generations to be simulated