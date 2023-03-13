## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

## Core equation for population (assume constant size here, at first):

# number of individuals with AA and AB are determined by drawing from Multinom(n, [p1,p2])
# with n=pop_size and p1 = (p + 1/2q) + (1/2 * (p + 1/2 q)) and p2 = (r + 1/2 q) + (1/2 * (p+q+r)) + (p + 1/2 q)
# the number of BB individuals is determined by pop_size - number of AA - number of AB

p[1] <- .2
p[2] <- .8


update(x[]) <- y[i]
dim(x) <- 2
y[] <- rmultinom(pop_size, p) 
### We checked with Rich and rmultinom for vector is currently not supported for odin.dust
### instead, I will implement multiple binomial draws
initial(x[]) <- initx[i]
initx[] <- user()
dim(initx) <- 2
dim(p) <- 2
dim(y) <- 2


#multiProb[1] <- 0.2
#multiProb[2] <- 0.8
#nextGen[] <- rmultinom(10,multiProb)


#nextGen[] <- rmultinom(pop_size, multiProb)

update(AA) <- y[1]
update(AB) <- y[2]
update(BB) <- pop_size - AA - AB # all individuals that are not AA or AB must be BB


## Initial states:
initial(AA) <- Pop_ini[1] # deterministic, user-based start value
initial(AB) <- Pop_ini[2] # deterministic, user-based start value
initial(BB) <- Pop_ini[3]




## User defined parameters - default in parentheses:
dim(Pop_ini) <- 3
Pop_ini[] <- user() 

#dim(multiProb) <- 2
#dim(nextGen) <- 2
pop_size <- user(100) # population size
#gen <- user() # number of generations to be simulated