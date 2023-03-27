## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

# calculate probabilities based on relative sizes of groups and their fitness (does not necessarily sum to 1):
pop_size <- sum(Pop[1:species_no])
probs[] <- (1 + GenotypeFitness[i]) * Pop[i]/pop_size

# number of individuals with Genotype 1 is determined by drawing from Bin(n, p)
# with n=pop_size and p = probs[1]/sum(probs) (relative abundance of Genotype 1 in the population)
# the probabilities need to be normalised because of the fitness term 
#y[1] <- if (probs[1]/sum(probs[1:species_no]) < 1)
#      rbinom(pop_size, probs[1]/sum(probs[1:species_no])) else pop_size

# The n and p for the remaining binomial draws have to be adapted
# because we only want to draw (n-x), i.e. those individuals that are not yet determined to be of Genotype 1 (and then 2, 3 and so on)
# and we need to adapt the probabilities accordingly to keep the population size the same

#y[2:(species_no)] <- if (probs[i]/sum(probs[i:species_no]) < 1)
#  rbinom(pop_size - sum(y[1:(i-1)]), probs[i]/sum(probs[i:species_no])) else pop_size - sum(y[1:(i-1)])

y[] <- if (probs[i]/sum(probs[1:species_no]) < 1)
        rbinom(pop_size, probs[i]/sum(probs[1:species_no])) else pop_size

#y[2:(species_no)] <- if (probs[i]/sum(probs[1:species_no]) < 1)
#    rbinom(pop_size, probs[i]/sum(probs[1:species_no])) else pop_size
  
  
# pop_size is now not strictly constant
# but the way the binomial draws work, it is at least very unlikely that it changes
# could I just keep Pop_size the same and not change p?
# but I do think that those are not multinomial draws. What kind of distribution would I get from that? maybe test that tomorrow!

## Core equation for population (assume constant size here):
update(Pop[]) <- y[i]
update(Pop_size) <- sum(y[1:species_no])

## Initial states:
#Pop_ini_norm[] <- as.integer(Pop_size * Pop_ini[i]/sum(Pop_ini[1:species_no])) # normalizing input vector (then it does not matter whether input consists of absolute or relative values)
# maybe a check that Pop_ini sums to 1 or pop_size would be nice. Not sure how to implement such a check though.

#as.integer(pop_size * Pop_ini_norm[i])
initial(Pop[]) <-  Pop_ini[i] # deterministic, user-based start value
initial(Pop_size) <- sum(Pop_ini[1:species_no])

FitnessMatrix[,] <- (GeneFitness[i]*Genotypes[i,j]) # calculate the fitness for genes in each genotype
GenotypeFitness[] <- sum(FitnessMatrix[,i]) #calculate overall fitness for all genotype (sum over genes*fitness)

## User defined parameters - default in parentheses:
#pop_size <- user() # population size (number of individuals in the whole population, constant)
species_no <- user() #number of species / strains / genotypes in the population
gene_no <- user() # number of genes in the data set
Pop_ini[] <- user() # initial frequency of Genotypes
GeneFitness[] <- user() # fitness vector for different genes
Genotypes[,] <- user() # each column is a genotype, giving the information which genes are present in that genotype and which are not

dim(Pop_ini) <- species_no
#dim(Pop_ini_norm) <- species_no
dim(FitnessMatrix) <- c(gene_no,species_no) # a matrix that stores the presence and the fitness for each gene and genotype
dim(Genotypes) <- c(gene_no, species_no) # we have in each column the genes (present/not present, i.e. 1/0) of one genotype
dim(Pop) <- species_no
dim(y) <- species_no
dim(probs) <- species_no
dim(GeneFitness) <- gene_no
dim(GenotypeFitness) <- species_no