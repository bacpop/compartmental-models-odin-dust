### R script for fitting my different model versions on the cluster
# run using 
# Rscript FitModelsOnCluster.R "ggCaller" "PopPUNK"
# pre-processed data files should be in the same folder

### Loading packages
# install.packages("drat") # -- if you don't have drat installed
# drat:::add("ncov-ic")
# install.packages("odin.dust")
library(odin.dust)
#install.packages("mcstate")
library(mcstate)
#install.packages("mcstate")
library(mcstate)
library(coda)

## command line arguments
args <- commandArgs(trailingOnly = TRUE)

#stop script if no arguments
if(length(args)==0){
  print("Please let me know which version of the model you want to run!")
  print("As the first argument, specify ggCaller or COGtriangles.")
  print("As the second argument, specify PopPUNK or manualSeqClusters.")
  stop("Requires command line argument.")
}


# read in model from file
WF <- odin.dust::odin_dust("NFDS_Model_logisticYshift.R")

# likelihood for fitting:
ll_pois <- function(obs, model) {
  exp_noise <- 1e6
  
  if (is.na(obs)) {
    # Creates vector of zeros in ll with same length, if no data
    ll_obs <- numeric(length(model))
  } else {
    lambda <- model + rexp(n = length(model), rate = exp_noise)
    ll_obs <- dpois(x = obs, lambda = lambda, log = TRUE)
  }
  ll_obs
}

combined_compare <- function(state, observed, pars = NULL) {
  result <- 0
  #data_size <- sum(mass_cluster_freq_1)
  #model_size = 15000
  data_size <- sum(unlist(observed))
  model_size = sum(unlist(state))
  
  for (i in 1:mass_clusters){
    result <- result + ll_pois(observed[[as.character(i)]], state[1+i, , drop = TRUE]/model_size * data_size)
  }
  result
}

if(args[1] == "ggCaller" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("PopPUNK_clusters.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "ggCPP_intermed_gene_presence_absence_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  model_start_pop <- readRDS(file = "PP_model_start_pop.rds")
  delta_ranking <- readRDS(file = "ggC_inv_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "PP_mass_cluster_freq_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "PP_mass_cluster_freq_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "PP_mass_cluster_freq_3.rds")
  mass_VT <- readRDS(file = "PP_mass_VT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
  output_filename <- "ggCaller_PopPUNK_logYshift"
} else if(args[1] == "COGtriangles" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("PopPUNK_clusters.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "PP_intermed_gene_presence_absence_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  model_start_pop <- readRDS(file = "PP_model_start_pop.rds")
  delta_ranking <- readRDS(file = "inv_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "PP_mass_cluster_freq_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "PP_mass_cluster_freq_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "PP_mass_cluster_freq_3.rds")
  mass_VT <- readRDS(file = "PP_mass_VT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
  output_filename <- "COGtriangles_PopPUNK_logYshift"
} else if(args[1] == "ggCaller" & args[2] == "manualSeqClusters"){
  seq_clusters <- readRDS("Mass_Samples_accCodes.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "ggC_intermed_gene_presence_absence_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  model_start_pop <- readRDS(file = "model_start_pop.rds")
  delta_ranking <- readRDS(file = "ggC_inv_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "mass_cluster_freq_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "mass_cluster_freq_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "mass_cluster_freq_3.rds")
  mass_VT <- readRDS(file = "mass_VT.rds")
  mass_clusters <- length(unique(seq_clusters$SequenceCluster))
  avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
  output_filename <- "ggCaller_manSeqClusters_logYshift"
} else if(args[1] == "COGtriangles" & args[2] == "manualSeqClusters"){
  seq_clusters <- readRDS("Mass_Samples_accCodes.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "intermed_gene_presence_absence_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  model_start_pop <- readRDS(file = "model_start_pop.rds")
  delta_ranking <- readRDS(file = "inv_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "mass_cluster_freq_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "mass_cluster_freq_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "mass_cluster_freq_3.rds")
  mass_VT <- readRDS(file = "mass_VT.rds")
  mass_clusters <- length(unique(seq_clusters$SequenceCluster))
  avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
  output_filename <- "COGtriangles_manSeqClusters_logYshift"
}


# process data with particle filter:
dt <- 1/36 # we assume that the generation time of Strep. pneumo is 1 month
# we have data from 2001, 2004, 2007, so we want 3 (years) * 12 (months) = 36 updates in-between

peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(mass_cluster_freq_2, mass_cluster_freq_3))
names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))

fitting_mass_data <- mcstate::particle_filter_data(data = peripost_mass_cluster_freq,
                                                   time = "year",
                                                   rate = 1 / dt,
                                                   initial_time = 0)

det_filter <- particle_deterministic$new(data = fitting_mass_data,
                                         model = WF,
                                         compare = combined_compare)

# Using MCMC to infer parameters
pmcmc_L <- mcstate::pmcmc_parameter("L", 0.15, min = 0, max = 1)
pmcmc_Y <- mcstate::pmcmc_parameter("K", 0.02, min = 0, max = 1)
pmcmc_x0 <- mcstate::pmcmc_parameter("x0", 0.2, min = 0, max = 1)
pmcmc_m <- mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1)
pmcmc_v <- mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)
species_no <- mass_clusters
no_clusters <- mass_clusters
gene_no <- nrow(intermed_gene_presence_absence_consensus_matrix)

Pop_ini <- model_start_pop
Pop_eq <- model_start_pop
Genotypes <- intermed_gene_presence_absence_consensus_matrix

capacity <- sum(model_start_pop)
delta <- delta_ranking
vaccTypes <- mass_VT
vacc_time <- 0
dt <- 1/36
migVec <- avg_cluster_freq

complex_params <- c(Pop_ini, Pop_eq, Genotypes, capacity, delta, vaccTypes, species_no, gene_no, vacc_time, dt, migVec)

make_transform <- function(p) {
  function(theta){
    c(list(Pop_ini = p[1:mass_clusters],
           Pop_eq = p[(mass_clusters +1) : (mass_clusters + mass_clusters)],
           Genotypes = matrix(p[(mass_clusters + mass_clusters + 1): ((mass_clusters + mass_clusters + 1) + (gene_no * species_no) - 1)], nrow = gene_no, ncol = species_no),
           capacity = p[((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 1],
           delta = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2) : (((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no -1)],
           vaccTypes = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) : ((((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters -1)],
           species_no = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters],
           gene_no = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters + 1],
           vacc_time = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters + 2],
           dt = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters + 3],
           migVec = p[((((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters + 4):((((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters + 4 + species_no - 1)]), as.list(theta))
  }
}

transform <- function(x) {
  make_transform(complex_params)}
proposal_matrix <- diag(0.1, 5) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
# here, all parameters are proposed independently. 
# think about this, this might not actually be true
#mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
#mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0.075, max = 1), mcstate::pmcmc_parameter("sigma_w", 0.05, min = 0, max = 0.0749), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
#= make_transform(c(Pop_ini, Pop_eq, Genotypes, capacity, delta, vaccTypes, species_no, gene_no, vacc_time)))
#mcmc_pars$names()
#mcmc_pars$model(mcmc_pars$initial())

# read this: https://mrc-ide.github.io/mcstate/reference/pmcmc_parameters.html
# it explains how to not fit all parameters but just the ones I want
# non-scalar parameters have to be transformed for this.

mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("L", 0.1432, min = 0, max = 1), mcstate::pmcmc_parameter("Y", 0.02, min = 0, max = 1), mcstate::pmcmc_parameter("x0", 0.2, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 0.2), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 0.5)), proposal_matrix, make_transform(complex_params))
mcmc_pars$initial()

det_filter <- particle_deterministic$new(data = fitting_mass_data,
                                         model = WF,
                                         compare = combined_compare)

n_steps <- 1000
n_burnin <- 0

control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE, 
  save_trajectories = TRUE,
  progress = TRUE,
  adaptive_proposal = TRUE,
  n_chains = 4)
det_pmcmc_run <- mcstate::pmcmc(mcmc_pars, det_filter, control = control)
processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run, burnin = 250, thin = 1)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
print(parameter_mean_hpd)

det_mcmc1 <- coda::as.mcmc(cbind(det_pmcmc_run$probabilities, det_pmcmc_run$pars))
pdf(file = paste(output_filename,"det_mcmc1.pdf",sep = "_"),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 12)
plot(det_mcmc1)
dev.off()
print("det_mcmc_1 log likelihood")
processed_chains$probabilities[nrow(processed_chains$probabilities),2]
print("det_mcmc_1 mean log likelihood")
mean(processed_chains$probabilities[,2])
det_proposal_matrix <- cov(processed_chains$pars)
#det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0.075, max = 0.22), mcstate::pmcmc_parameter("sigma_w", 0.05, min = 0.000001, max = 0.0749), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 0.2), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 0.5)), det_proposal_matrix, make_transform(complex_params))
det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("L", parameter_mean_hpd[1], min = 0, max = 1), mcstate::pmcmc_parameter("Y", parameter_mean_hpd[2], min = 0, max = 1), mcstate::pmcmc_parameter("x0", parameter_mean_hpd[3], min = 0, max = 1), mcstate::pmcmc_parameter("m", parameter_mean_hpd[4], min = 0, max = 0.2), mcstate::pmcmc_parameter("v", parameter_mean_hpd[5], min = 0, max = 0.5)), det_proposal_matrix, make_transform(complex_params))

det_filter <- particle_deterministic$new(data = fitting_mass_data,
                                         model = WF,
                                         compare = combined_compare)

n_steps <- 20000
n_burnin <- 0


control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE, 
  save_trajectories = TRUE,
  progress = TRUE,
  adaptive_proposal = TRUE,
  n_chains = 4)
det_pmcmc_run2 <- mcstate::pmcmc(det_mcmc_pars, det_filter, control = control)
par(mfrow = c(1,1))

det_mcmc2 <- coda::as.mcmc(cbind(det_pmcmc_run2$probabilities, det_pmcmc_run2$pars))

pdf(file = paste(output_filename,"det_mcmc2.pdf",sep = "_"),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 12)
plot(det_mcmc2)
dev.off()

processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run2, burnin = 1000, thin = 1)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
parameter_mean_hpd
print("det_mcmc_2 log likelihood")
processed_chains$probabilities[nrow(processed_chains$probabilities),2]
print("det_mcmc_2 mean log likelihood")
mean(processed_chains$probabilities[,2])
