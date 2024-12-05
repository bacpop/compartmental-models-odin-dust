# model definition derived from Corander et al. (2017)

## Definition of the time-step and output as "time"
initial(time) <- 0
update(time) <- (step + 1) * dt

# calculate probabilities based on relative sizes of groups and their fitness (does not necessarily sum to 1):
Pop_size <- sum(Pop[1:species_no,1:sero_no])

# frequency dependent selection (i.e. the fitness of the genotypes are not constant but based on how frequent the genes/loci they contain are in the population)
# based on Corander et al. (2017)
# frequency of each gene at current time:
gene_freq[,] <-  Genotypes[i,j] * sum(Pop[j,1:sero_no])
freq[] <- sum(gene_freq[i,1:species_no]) / Pop_size

# overall deviation of loci for genomes
# idea 1: make a Boolean vector for strongly selected genes
#delta_bool[] <- if ((delta[i] <= prop_f * gene_no)) 1 else 0
#pi_f_freq[,] <-  Genotypes[i,j] * (eq[i] - freq[i]) * delta_bool[i]
#pi_f_genotypes[] <- sum(pi_f_freq[1:gene_no,i])

#pi_w_freq[,] <-  Genotypes[i,j] * (eq[i] - freq[i]) * (1 - delta_bool[i])
#pi_w_genotypes[] <- sum(pi_w_freq[1:gene_no,i])

pi_f_freq[,] <- if ((delta[i] <= prop_f * gene_no)) Genotypes[i,j] * (eq[i] - freq[i]) else 0
pi_f_genotypes[] <- sum(pi_f_freq[1:gene_no,i])

pi_w_freq[,] <-  if ((delta[i] <= prop_f * gene_no)) 0 else Genotypes[i,j] * (eq[i] - freq[i])
pi_w_genotypes[] <- sum(pi_w_freq[1:gene_no,i])

# Genotype specific probability to produce offspring
# those are the individuals' probabilities multiplied by the number of individual that have this genotype
#probs[,] <- ((1 + exp(sigma_f) + exp(sigma_w))^pi_f_genotypes[i] * (1 + exp(sigma_w))^pi_w_genotypes[i]) * Pop[i,j] * (1- (as.integer(time >= vacc_time) * vaccTypes[j] * v))
#probs[,] <- ((1 + exp(sigma_f + sigma_w))^pi_f_genotypes[i] * (1 + exp(sigma_w))^pi_w_genotypes[i]) * Pop[i,j] * (1- (as.integer(time >= vacc_time) * vaccTypes[j] * v))
# this (above) is a nice idea but it does not work with the fitting
probs[,] <- ((1 + exp(sigma_f))^pi_f_genotypes[i] * (1 + exp(sigma_w))^pi_w_genotypes[i]) * Pop[i,j] * (1- (as.integer(time >= vacc_time) * vaccTypes[j] * v))
#probs[,] <- (1 + exp(sigma_f))^pi_f_genotypes[i] * Pop[i,j] * (1- (as.integer(time >= vacc_time) * vaccTypes[j] * v))



# The lambda of the poisson distribution then consists of the normalised probability and a factor that describes how close we are to the capacity and a factor for the population size.
# You could simplify this by cancelling out the Pop_size (which makes sense because we will loose less accuracy but it will make it a bit less easy to interpret.)

y[,] <- if ((probs[i,j]/sum(probs[1:species_no,1:sero_no])) < 1) #not strictly necessary for poisson but probs>1 should be avoided anyway
  rpois(capacity * (probs[i,j] / sum(probs[1:species_no,1:sero_no])) * (1-exp(m)) ) else rpois(capacity * 1 *(1-exp(m)) )


# m is the migration rate
# fitness of individuals in the community is reduced by this rate
# determining migration number:
mig_num <- rbinom(capacity, exp(m))
Pop_mig[,] <- rbinom(mig_num, migVec[i,j])
# this is a very simple implementation of migration. It does not care what the rates of the different genotypes are and just makes uniform, random draws from all existing genotypes.

## Core equation for population (assume constant size here):
update(Pop[,]) <- y[i,j] + Pop_mig[i,j]
#update(probs2[]) <- eq[i]
update(Pop_tot[]) <- sum(y[i,]) + sum(Pop_mig[i,])

initial(Pop[,]) <- Pop_ini[i,j] # deterministic, user-based start value
initial(Pop_tot[]) <- sum(Pop_ini[i,])
#initial(probs2[]) <- (1 + 0.02) * Pop_ini[i]

#calculate equilibrium frequencies of genes
gene_eq[,] <- Genotypes[i,j] * Pop_eq[j]
eq[] <- sum(gene_eq[i,1:species_no]) / sum(Pop_eq[1:species_no])

## User defined parameters - default in parentheses:
#pop_size <- user() # population size (number of individuals in the whole population, constant)
dt <- user(1)
species_no <- user() #number of species / strains / genotypes in the population
sero_no <- user() # number of serotypes
gene_no <- user() # number of genes in the data set
Pop_ini[,] <- user() # initial frequency of Genotypes
Pop_eq[] <- user()
capacity <- user()
sigma_f <- user()
sigma_w <- user()
prop_f <- user()
delta[] <- user()
m <- user() # migration rate
#GeneFitness[] <- user() # fitness vector for different genes
Genotypes[,] <- user() # each column is a genotype, giving the information which genes are present in that genotype and which are not
vaccTypes[] <- user() # Boolean vector (0/1) whether genotype is affected by vaccine (1) or not (0)
v <- user() # effect size of vaccine on vaccine genotypes
vacc_time <- user() # time when the vaccination happens / starts to have an effect
migVec[,] <- user()

#dim(probs2) <- gene_no
dim(gene_freq) <- c(gene_no, species_no)
dim(freq) <- gene_no
dim(Pop_ini) <- c(species_no,sero_no)
dim(Pop_eq) <- species_no
dim(Pop_mig) <- c(species_no,sero_no)
dim(Pop_tot) <- species_no
dim(gene_eq) <- c(gene_no, species_no) #frequency of genes at equilibrium
dim(eq) <- gene_no
dim(pi_f_freq) <- c(gene_no, species_no)
dim(pi_f_genotypes) <- species_no
dim(pi_w_freq) <- c(gene_no, species_no)
dim(pi_w_genotypes) <- species_no
dim(delta) <- gene_no
#dim(delta_bool) <- gene_no
#dim(FitnessMatrix) <- c(gene_no,species_no) # a matrix that stores the presence and the fitness for each gene and genotype
dim(Genotypes) <- c(gene_no, species_no) # we have in each column the genes (present/not present, i.e. 1/0) of one genotype
dim(Pop) <- c(species_no,sero_no)
dim(y) <- c(species_no,sero_no)
dim(probs) <- c(species_no,sero_no)
dim(vaccTypes) <- sero_no
dim(migVec) <- c(species_no,sero_no)
#dim(GeneFitness) <- gene_no
#dim(GenotypeFitness) <- species_no