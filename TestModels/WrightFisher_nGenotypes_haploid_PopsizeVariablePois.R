## Definition of the time-step and output as "time"
dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt

# calculate probabilities based on relative sizes of groups and their fitness (does not necessarily sum to 1):
Pop_size <- sum(Pop[1:species_no])

# frequency dependent selection (i.e. the fitness of the genotypes are not constant but based on how frequent the genes/loci they contain are in the population)
# based on Corander et al. (2017)
# frequency of each gene at current time:
gene_freq[,] <-  Genotypes[i,j] * Pop[j]
freq[] <- sum(gene_freq[i,1:species_no]) / Pop_size

# overall deviation of loci for genomes
pi_freq[,] <- Genotypes[i,j] * (eq[i] - freq[i])
pi_genotypes[] <- sum(pi_freq[1:gene_no,i])

# Genotype specific probability to produce offspring
# those are the individuals' probabilities multiplied by the number of individual that have this genotype
#probs[] <- (1 + GenotypeFitness[i])^pi_genotypes[i] * Pop[i]/Pop_size
#probs[] <- (1 + GenotypeFitness[i])^pi_genotypes[i] * Pop[i]
probs[] <- (1 + sigma)^pi_genotypes[i] * Pop[i]


# Okay, my current interpretation is:
# probs stores the probability for a genotype to produce offspring. That has to be relative to the genotype frequency.
# I want to normalise this probability because I want that the sum of the probabilities of all genotypes to have offspring is 1.
# Why? I think because I want X~Pois(lambda), Y~Pois(mu) independent -> X+Y~Pois(lambda+mu) to hold true and make sense.

# The lambda of the poisson distribution then consists of the normalised probability and a factor that describes how close we are to the capacity and a factor for the population size.
# You could simplify this by cancelling out the Pop_size (which makes sense because we will loose less accuracy but it will make it a bit less easy to interpret.)

#old version that is a bit easier to interpret:
#y[] <- if (probs[i]/sum(probs[1:species_no]) < 1) #not strictly necessary for poisson but probs>1 should be avoided anyway
#  rpois((capacity/Pop_size)*(probs[i]/sum(probs[1:species_no]))*Pop_size) else rpois((capacity/Pop_size)*(1)*Pop_size)


y[] <- if (probs[i]/sum(probs[1:species_no]) < 1) #not strictly necessary for poisson but probs>1 should be avoided anyway
  rpois(capacity * (probs[i] / sum(probs[1:species_no]))) else rpois(capacity * 1)


## Core equation for population (assume constant size here):
update(Pop[]) <- y[i] 
#update(probs2[]) <- eq[i]

initial(Pop[]) <- Pop_ini[i] # deterministic, user-based start value
#initial(probs2[]) <- (1 + 0.02) * Pop_ini[i]

#FitnessMatrix[,] <- (GeneFitness[i]*Genotypes[i,j]) # calculate the fitness for genes in each genotype
#GenotypeFitness[] <- sum(FitnessMatrix[1:gene_no,i]) #calculate overall fitness for all genotype (sum over genes*fitness)


#calculate equilibrium frequencies of genes
gene_eq[,] <- Genotypes[i,j] * Pop_eq[j]
eq[] <- sum(gene_eq[i,1:species_no]) / sum(Pop_eq[1:species_no])

## User defined parameters - default in parentheses:
#pop_size <- user() # population size (number of individuals in the whole population, constant)
species_no <- user() #number of species / strains / genotypes in the population
gene_no <- user() # number of genes in the data set
Pop_ini[] <- user() # initial frequency of Genotypes
Pop_eq[] <- user()
capacity <- user()
sigma <- user()
#GeneFitness[] <- user() # fitness vector for different genes
Genotypes[,] <- user() # each column is a genotype, giving the information which genes are present in that genotype and which are not

#dim(probs2) <- gene_no
dim(gene_freq) <- c(gene_no, species_no)
dim(freq) <- gene_no
dim(Pop_ini) <- species_no
dim(Pop_eq) <- species_no
dim(gene_eq) <- c(gene_no, species_no) #frequency of genes at equilibrium
dim(eq) <- gene_no
dim(pi_freq) <- c(gene_no, species_no)
dim(pi_genotypes) <- species_no
#dim(FitnessMatrix) <- c(gene_no,species_no) # a matrix that stores the presence and the fitness for each gene and genotype
dim(Genotypes) <- c(gene_no, species_no) # we have in each column the genes (present/not present, i.e. 1/0) of one genotype
dim(Pop) <- species_no
dim(y) <- species_no
dim(probs) <- species_no
#dim(GeneFitness) <- gene_no
#dim(GenotypeFitness) <- species_no