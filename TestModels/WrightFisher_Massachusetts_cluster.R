library(odin.dust)
WF_nG_h_vP <- odin.dust::odin_dust("WrightFisher_nGenotypes_haploid_PopsizeVariablePois.R")

# reading in the cluster produced by PopPUNK
clusters <- read.csv("../Data/refined_modelfitk3_clusters.csv")
no_clusters <- max(clusters[,2]) # number of clusters in dataset

# reading in the gene presence absence matrix produced by ggCaller
gene_presence_absence <- read.csv("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/gene_presence_absence.csv", header=FALSE)

# converting the gene presence absence matrix into a boolean df (0 = gene not present, 1 = gene present)
convert_to_bool <- function(x){
  if (x=="") 0 else 1
}
bool_gene_presence_absence <- gene_presence_absence
bool_gene_presence_absence[2:nrow(bool_gene_presence_absence),4:ncol(bool_gene_presence_absence)] <- apply(gene_presence_absence[2:nrow(gene_presence_absence),4:ncol(gene_presence_absence)],c(1,2), convert_to_bool)

# calculate frequency of genes to only keep those which appear in 5-95% of the genomes
#gene_freq <- rep(0, nrow(bool_gene_presence_absence)-1)
#for (i in 1:length(gene_freq)) {
#  gene_freq[i] <- sum(as.integer(bool_gene_presence_absence[i+1,4:ncol(bool_gene_presence_absence)]))
#}
# gene_freq <- gene_freq / (length(bool_gene_presence_absence[2,])-3)

sum_as_int <- function(x){
  sum(as.integer(x))
}

gene_freq <- rep(0, nrow(bool_gene_presence_absence)-1)
gf_bool_gene_presence_absence <- bool_gene_presence_absence[2:nrow(bool_gene_presence_absence),4:ncol(bool_gene_presence_absence)]

gene_freq <- apply(gf_bool_gene_presence_absence,1, sum_as_int)
gene_freq <- as.vector(gene_freq / (length(bool_gene_presence_absence[2,])-3))



# create a dataframe that only contains the genes that appear in 5-95% of the genomes
filtered_bool_gene_presence_absence <- data.frame(matrix(nrow = sum(gene_freq <= 0.95 & gene_freq >= 0.05), ncol = length(bool_gene_presence_absence[1,])))
counter <- 1
for (i in 1:length(gene_freq)){
  if (0.05 <= gene_freq[i] & gene_freq[i] <= 0.95 ){
    filtered_bool_gene_presence_absence[counter,] <- bool_gene_presence_absence[i+1,]
    counter <- counter + 1
  }
}
colnames(filtered_bool_gene_presence_absence) <- bool_gene_presence_absence[1,]
rownames(filtered_bool_gene_presence_absence) <- 1:nrow(filtered_bool_gene_presence_absence)

### make consensus genome for clusters
# attempt 1: always let majority decide (if more than 50% in the cluster don't have gene then 0, else 1)
# should be easy to do with median

cluster_gene_presence_absence <- data.frame(matrix(nrow = nrow(filtered_bool_gene_presence_absence), ncol = no_clusters+3))
colnames(cluster_gene_presence_absence)[1:3] <- colnames(filtered_bool_gene_presence_absence)[1:3]
cluster_gene_presence_absence[1:3] <- filtered_bool_gene_presence_absence[1:3]
colnames(cluster_gene_presence_absence)[4:ncol(cluster_gene_presence_absence)] <- 1:no_clusters

cons_genomes <- function(x){
  as.double(median(as.integer(x)))
}

for (i in 1:no_clusters){
  curr_cluster <- clusters[which(clusters[,"Cluster"]==i),1] # select all genomes in cluster i
  curr_genomes <- as.matrix(filtered_bool_gene_presence_absence[,curr_cluster])
  cluster_gene_presence_absence[,i+3] <- apply(curr_genomes,1,cons_genomes)
}
matrix_cluster_gene_presence_absence <- as.matrix((cluster_gene_presence_absence[,4:ncol(cluster_gene_presence_absence)]))

# well, I'd say that many, many genes are still quite diversely present within the same cluster
# that suggest that taking one consensus genome per cluster might not work

# attempt 2: record relative frequencies. Probably not advantageous.
# attempt 3: keep track of all existing variants. This can be done later

#calculate the frequency of the gene clusters, well actually absolute numbers atm
cluster_freq <- rep(0,no_clusters)
for (i in 1:no_clusters){
  cluster_freq[i] <- length(which(clusters[,"Cluster"]==i))
}

#need the information of when probes were taken
accNo_to_filename <- read.delim("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/filereport_read_run_PRJEB2632_tsv.txt")
accNo_to_filename <- accNo_to_filename[,c(1,8)]
accNo_to_filename[,3] <- matrix(unlist(strsplit(accNo_to_filename[,2],"/")), ncol=6, byrow = TRUE)[,6]
accNo_to_filename[,3] <- matrix(unlist(strsplit(accNo_to_filename[,3],"[.]")), ncol=2, byrow = TRUE)[,1]
accNo_to_filename <- accNo_to_filename[,c(1,3)]
colnames(accNo_to_filename) <- c(colnames(accNo_to_filename)[1], "filenames")
accNoToFilename <- c()
accNoToFilename[as.character(accNo_to_filename[,2])] <- accNo_to_filename[,1]

library(readxl)
Croucher_seqYears <- read_excel("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/Croucher_41588_2013_BFng2625_MOESM28_ESM.xlsx")
Croucher_seqYears <- Croucher_seqYears[,c(1,5)]
sequenceYear <- c()
for (i in 1:nrow(Croucher_seqYears)){
  sequenceYear[as.character(Croucher_seqYears[i,1])]<- Croucher_seqYears[i,2] 
}

# add year to the clusters data set
clusters$seqYear <- rep(0, nrow(clusters))
for (i in 1:nrow(clusters)){
  clusters$seqYear[i] <- sequenceYear[accNoToFilename[clusters$Taxon[i]]]
}

#calculate the frequency of the gene clusters and year
cluster_freq_1 <- rep(0,no_clusters)
cluster_freq_2 <- rep(0,no_clusters)
cluster_freq_3 <- rep(0,no_clusters)
for (i in 1:no_clusters){
  cluster_freq_1[i] <- length(which(clusters[which(clusters[,"Cluster"]==i),]$seqYear == 2001))
  cluster_freq_2[i] <- length(which(clusters[which(clusters[,"Cluster"]==i),]$seqYear == 2004))
  cluster_freq_3[i] <- length(which(clusters[which(clusters[,"Cluster"]==i),]$seqYear == 2007))
}

# make a data frame that contains this information
library(dplyr)
fitting_cluster_freq_df <- data.frame("year" = c(2001, 2004, 2007), rbind(cluster_freq_1, cluster_freq_2, cluster_freq_3))
names(fitting_cluster_freq_df) <- c("year", as.character(1:62))

#plot gene freqs
data_gene_freq_1_m <-  matrix_cluster_gene_presence_absence * cluster_freq_1
#data_gene_freq_1 <- rep(0, length(data_gene_freq_1_m[,1]))
data_gene_freq_1 <- apply(data_gene_freq_1_m[,1:62], 1, sum)

data_gene_freq_2_m <-  matrix_cluster_gene_presence_absence * cluster_freq_2
data_gene_freq_2 <- apply(data_gene_freq_2_m[,1:62], 1, sum)

data_gene_freq_3_m <-  matrix_cluster_gene_presence_absence * cluster_freq_3
data_gene_freq_3 <- apply(data_gene_freq_3_m[,1:62], 1, sum)

#par(mfrow=c(1,3))
#plot(data_gene_freq_1 / nrow(subset(clusters,clusters$seqYear==2001)), data_gene_freq_2 / nrow(subset(clusters,clusters$seqYear==2004)))
#plot(data_gene_freq_1 / nrow(subset(clusters,clusters$seqYear==2001)), data_gene_freq_3 / nrow(subset(clusters,clusters$seqYear==2007)))
#plot(data_gene_freq_2 / nrow(subset(clusters,clusters$seqYear==2004)), data_gene_freq_3 / nrow(subset(clusters,clusters$seqYear==2007)))
# very strong correlation between gene frequencies before and after vaccination
# especially second plot looks a lot like figure 2 d, the first, in Corander et al. (though they only did pre vs. post vacc)
#par(mfrow=c(1,1))
#plot((data_gene_freq_1 + data_gene_freq_2) / (nrow(subset(clusters,clusters$seqYear==2001)) + nrow(subset(clusters,clusters$seqYear==2004))), data_gene_freq_3 / nrow(subset(clusters,clusters$seqYear==2007)))

# correlation of clusters
#par(mfrow=c(1,3))
#plot(cluster_freq_1 / sum(cluster_freq_1), cluster_freq_2 / sum(cluster_freq_2))
#plot(cluster_freq_1 / sum(cluster_freq_1), cluster_freq_3 / sum(cluster_freq_3))
#plot(cluster_freq_2 / sum(cluster_freq_2), cluster_freq_3 / sum(cluster_freq_3))
# the second plot (2001 vs. 2007) looks a lot like figure 2 c, the first, in Corander et al.

#par(mfrow=c(1,1))
#plot((cluster_freq_1 + cluster_freq_2)/ (sum(cluster_freq_1) + sum(cluster_freq_2)), cluster_freq_3 / sum(cluster_freq_3))

# read in information on vaccine types
#library(readxl)
vaccine_types <- read_excel("../Data/Corander_suppData3.xlsx")
vaccine_types <- vaccine_types[,c(2,3,4,6)]
vaccine_types_mass <- vaccine_types[which(vaccine_types$Population == "Massachusetts"),]
isVT <- c()
for (i in 1:nrow(vaccine_types_mass)){
  isVT[vaccine_types_mass$`Accession Code`[i]] <- if (vaccine_types_mass$`Vaccine Type`[i] == "VT") {0} else {1}
}
# I want to create a simple vector instead of a dictionary because this will be easier to use when modelling
# and I also want a consensus for the clusters...
isVT_vec <- rep(0, no_clusters)
for (i in 1:no_clusters){
  curr_cluster <- subset(clusters, clusters$Cluster == i)[,1] # select all genomes in cluster i
  curr_vacc <- accNoToFilename[curr_cluster]
  isVT_vec[i] <- median(isVT[curr_vacc])
}

# make a grouped bar plot similar to Croucher et al.
#library(ggplot2)

#cluster_names <- rep(1:no_clusters, 3)
#cluster_names <- sort(cluster_names)
#cluster_names <- as.character(cluster_names)
#seqTimes <- rep(c("2001", "2004", "2007"), no_clusters)
#cluster_freq_all_times <- rep(0, no_clusters *3)
#for (i in 1:no_clusters) {
#  cluster_freq_all_times[3*(i-1) + 1 ] <- cluster_freq_1[i]/sum(cluster_freq_1)
#  cluster_freq_all_times[3 * (i-1) + 2] <- cluster_freq_2[i]/sum(cluster_freq_2)
#  cluster_freq_all_times[3 * (i-1) + 3] <- cluster_freq_3[i]/sum(cluster_freq_3)
#}
#cluster_freq_df <- data.frame(cluster_names, seqTimes, cluster_freq_all_times)

#ggplot(cluster_freq_df, aes(fill=seqTimes, y=cluster_freq_all_times, x=cluster_names)) +   geom_bar(position="dodge", stat="identity") +
#  scale_fill_manual("legend", values = c("2001" = "#E69F00", "2004" = "#56B4E9", "2007" = "#009E73")) +
#  scale_x_continuous(breaks = 1:no_clusters)

# I want to add low and high levels of selection to the model.
# for that I need to calculate the beta statistics
# first, calculate pre/peri and post vacc frequencies of genes:
pre_peri_vacc_gene_freq <- (data_gene_freq_1 + data_gene_freq_2) / (nrow(subset(clusters,clusters$seqYear==2001)) + nrow(subset(clusters,clusters$seqYear==2004)))
post_vacc_gene_freq <- data_gene_freq_3 / nrow(subset(clusters,clusters$seqYear==2007))

# calculate delta statistic (refer to Corander et al. for more info)
delta_data <- (post_vacc_gene_freq - pre_peri_vacc_gene_freq) ^ 2 / (1 - pre_peri_vacc_gene_freq * (1 - pre_peri_vacc_gene_freq))
delta_ranking <- rank(delta_data)

# implement sampling from pre vaccination population instead of directly using the 2001 population as input
# would it be better to do the sampling in the modelling frame work? maybe.
start_pop <- as.vector(as.double(rmultinom(1,50000,prob = cluster_freq_1/sum(cluster_freq_1))))

### Define model parameters according to the datasets
params_n_vP <- list(dt = 1/36, species_no = no_clusters,  gene_no = nrow(cluster_gene_presence_absence), Pop_ini = start_pop, Pop_eq = start_pop, capacity = sum(start_pop), Genotypes = matrix_cluster_gene_presence_absence, sigma_f = 0.14, sigma_w = 0.002, prop_f = 0.25, delta = delta_ranking, m = 0.03, vaccTypes = isVT_vec, v = 0.07, vacc_time = 100) 
### Running the model:
WFmodel_nG_h_vP <- WF_nG_h_vP$new(pars = params_n_vP,
                                  time = 1,
                                  n_particles = 10L,
                                  n_threads = 4L,
                                  seed = 1L)

### fitting
#install.packages("drat") # -- if you don't have drat installed
drat:::add("ncov-ic")
#install.packages("mcstate")
library(mcstate)

# process data with particle filter:
dt <- 1/36 # we assume that the generation time of Strep. pneumo is 1 month
# we have data from 2001, 2004, 2007, so we want 3 (years) * 12 (months) = 36 updates in-between

# adding additional times to df, filled with NAs
vecNA <- rep(NA, 62)
pfd_fitting_cluster_freq_df <- fitting_cluster_freq_df
pfd_fitting_cluster_freq_df[4,] <- pfd_fitting_cluster_freq_df[2,]
pfd_fitting_cluster_freq_df[5,] <- c(2005, vecNA)
pfd_fitting_cluster_freq_df[6,] <- c(2006, vecNA)
pfd_fitting_cluster_freq_df[7,] <- pfd_fitting_cluster_freq_df[3,]
pfd_fitting_cluster_freq_df[2,] <- c(2002, vecNA)
pfd_fitting_cluster_freq_df[3,] <- c(2003, vecNA)

pfd_simpleInts_fitting_cluster_freq_df <- pfd_fitting_cluster_freq_df
pfd_simpleInts_fitting_cluster_freq_df$year <- c(1,2,3,4,5,6,7) # this is a test
# I have all years from 2001-2007 now, but the some are NAs obviously
# and I renamed the years to 1-7, to not cause weird effects of the model running 2000 years before already

# No! This is unnecessary because it means more work for the fitting algorithm where we actually do not have data

#instead, I want to excluded 2001 data from fitting, since that is just the start (and the start population is sampled from that)

peripost_fitting_cluster_freq_df <- fitting_cluster_freq_df[2:3,]

mass_data <- mcstate::particle_filter_data(data = peripost_fitting_cluster_freq_df,
                                           time = "year",
                                           rate = 1 / dt,
                                           initial_time = 2001)

# heavily based on this: https://mrc-ide.github.io/mcstate/articles/sir_models.html
# trying to generalise the comparison function so that we can later fit to all strain data
# log-likelihood of Poisson count
ll_pois <- function(obs, model) {
  exp_noise <- 1e6
  if (is.na(obs)) {
    # Creates vector of zeros in ll with same length, if no data
    ll_obs <- numeric(length(model))
  } else {
    lambda <- model + rexp(n = length(model), rate = exp_noise)
    #ll_obs <- dpois(x = obs, lambda = lambda, log = TRUE) # this is the simplest prior possible
    # it would be good to think about a binomial, neg binomial or beta binomial version
    # ll_obs <- dbinom() 
    # for neg binom: maybe n is the capacity of the system? the model value is the mu, and the y is the data?
    ll_obs <- dnbinom(x = obs, size = capacity, mu = lambda, log = TRUE)
  }
  ll_obs
}

combined_compare <- function(state, observed, pars = NULL) {
  result <- 0
  for (i in 1:62){
    #print(observed[as.character(ind)])
    #print(state[1+i, , drop = TRUE])
    result <- result + ll_pois(observed[[as.character(i)]], state[1+i, , drop = TRUE])
  }
  result
}

### Inferring parameters

n_particles <- 100
# n_particles <- 100 nicer but takes about 35 minutes atm

filter <- mcstate::particle_filter$new(data = mass_data,
                                       model = WF_nG_h_vP,
                                       n_particles = n_particles,
                                       compare = combined_compare,
                                       seed = 1L)

#original:
params_n_vP$dt <- 1/36
#params_n_vP$dt <- 1/12
filter$run(save_history = TRUE, pars = params_n_vP)
#I have an issue with the filter since it only outputs three time point. Is that an issue? I think so.
# Or are those just the three time points that we have data to compare to?
# but I thought that was the reason to set dt <- 1/36, so that I have shorter time steps in the model than in the data set...?

# Using MCMC to infer parameters
pmcmc_sigma_f <- mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0)
pmcmc_sigma_w <- mcstate::pmcmc_parameter("sigma_w", 0.05, min = 0)
pmcmc_prop_f <- mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1)
pmcmc_m <- mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1)
pmcmc_v <- mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)
species_no <- no_clusters
gene_no <- nrow(cluster_gene_presence_absence)
Pop_ini <- cluster_freq_1
Pop_eq <- cluster_freq_1
Genotypes <-matrix_cluster_gene_presence_absence
capacity <- sum(cluster_freq_1)
delta <- delta_ranking
vaccTypes <- isVT_vec
vacc_time <- 100
dt <- 1/36


complex_params <- c(Pop_ini, Pop_eq, Genotypes, capacity, delta, vaccTypes, species_no, gene_no, vacc_time, dt, 0.15, 0.05, 0.25, 0.03, 0.05)
#set_names <- function(x, nms) {
#  names(x) <- nms
#  x
#}
#transform <- set_names(lapply(complex_params, make_transform), complex_params)
make_transform <- function(p) {
  list(Pop_ini = p[1:62],
       Pop_eq = p[63 : 124],
       Genotypes = matrix(p[125 : (125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1)], nrow = gene_no, ncol = no_clusters),
       capacity = p[(125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) +1],
       delta = p[((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) : (((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no -1)],
       vaccTypes = p[(((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no) : (((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no + 1 + no_clusters -2)],
       species_no = p[(((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no + 1 + no_clusters -1)],
       gene_no = p[(((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no + 1 + no_clusters)],
       vacc_time = p[(((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no + 1 + no_clusters) + 1],
       dt = p[(((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no + 1 + no_clusters) + 2],
       sigma_f = p[(((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no + 1 + no_clusters) + 3],
       sigma_w = p[(((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no + 1 + no_clusters) + 4],
       prop_f = p[(((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no + 1 + no_clusters) + 5],
       m = p[(((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no + 1 + no_clusters) + 6],
       v = p[(((125 + (nrow(matrix_cluster_gene_presence_absence) * ncol(matrix_cluster_gene_presence_absence))-1) + 2) + gene_no + 1 + no_clusters) + 7])
}

transform <- function(x) {
  make_transform(complex_params)}
proposal_matrix <- diag(0.1, 5) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
# here, all parameters are proposed independently. 
# think about this, this might not actually be true
#mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0), mcstate::pmcmc_parameter("sigma_w", 0.05, min = 0), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, transform)
#= make_transform(c(Pop_ini, Pop_eq, Genotypes, capacity, delta, vaccTypes, species_no, gene_no, vacc_time)))
#mcmc_pars$names()
#mcmc_pars$model(mcmc_pars$initial())
mcmc_pars$initial()
# read this: https://mrc-ide.github.io/mcstate/reference/pmcmc_parameters.html
# it explains how to not fit all parameters but just the ones I want
# non-scalar parameters have to be transformed for this.

n_steps <- 50
n_burnin <- 20

control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE, 
  save_trajectories = TRUE,
  progress = TRUE)
pmcmc_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)

# plot particle filter
plot_particle_filter <- function(history, true_history, times, obs_end = NULL) {
  if (is.null(obs_end)) {
    obs_end <- max(times)
  }
  
  
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  cols <- rainbow(62)
  
  matplot(times, t(history[2, , -1]), type = "l",
          xlab = "Time", ylab = "Number of individuals",
          col = cols[1], lty = 1, ylim = range(history[2:5, , -1]))
  for (species in 3:length(history[,1,1])) {
    matlines(times, t(history[species, , -1]), col = cols[species], lty = 1)
  }  
  for (species in 2:length(history[,1,1])) {
    matpoints(times[1:obs_end], (true_history[1:2,species , -1]), pch = 19, col = cols[species-1])
  }
  legend("left", lwd = 1, col = cols, legend = 1:62, bty = "n")
}
true_history <- peripost_fitting_cluster_freq_df[]
#par(mfrow = c(1,1))
#plot_particle_filter(pmcmc_run$trajectories$state, true_history, 1:2)

processed_chains <- mcstate::pmcmc_thin(pmcmc_run, burnin = n_burnin, thin = 2)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
parameter_mean_hpd

library(coda)
mcmc1 <- coda::as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars))

#summary(mcmc1)
#plot(mcmc1)
#coda::effectiveSize(mcmc1) #effective samle size
#1 - coda::rejectionRate(mcmc1) # acceptance rate of MCMC proposals

proposal_matrix <- cov(pmcmc_run$pars)
mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0), mcstate::pmcmc_parameter("sigma_w", 0.05, min = 0), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, transform)
proposal_matrix

control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE,
  save_trajectories = TRUE,
  progress = TRUE,
  n_chains = 2)
pmcmc_tuned_run <- mcstate::pmcmc(mcmc_pars, filter, control = control)

# runs in 4 minutes per chain right now
mcmc2 <- coda::as.mcmc(cbind(
  pmcmc_tuned_run$probabilities, pmcmc_tuned_run$pars))

true_history <- peripost_fitting_cluster_freq_df[]
par(mfrow = c(1,1))
plot_particle_filter(pmcmc_tuned_run$trajectories$state, true_history, 1:2)

summary(mcmc2)
plot(mcmc2)
coda::effectiveSize(mcmc2)
1 - coda::rejectionRate(mcmc2)
