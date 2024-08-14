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
WF <- odin.dust::odin_dust("NFDS_Model_VT_2vacc.R")

# likelihood for fitting:
ll_pois <<- function(obs, model) {
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
  #data_size <- sum(unlist(observed))
  data_size <- sum(unlist(observed[as.character(1:(length(unlist(observed))-4))]))
  model_size = sum(unlist(state[-1, , drop = TRUE]))
  exp_noise <- 1e6
  
  for (i in 1:(length(unlist(observed))-4)){ 
    state_name <- paste("sum_clust", i, sep = "")
    if (is.na(observed[[as.character(i)]])) {
      #Creates vector of zeros in ll with same length, if no data
      ll_obs <- numeric(length( state[state_name, , drop = TRUE]))
    } else {
      lambda <-  state[state_name, , drop = TRUE]/model_size * data_size + rexp(n = length( state[state_name, , drop = TRUE]/model_size * data_size), rate = exp_noise)
      ll_obs <- dpois(x = observed[[as.character(i)]], lambda = lambda, log = TRUE)
    }
    
    result <- result + ll_obs
  }
  result
}

if(args[1] == "Nepal" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("Nepal_PP.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "Nepal_ggCaller_intermed_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  model_start_pop <- readRDS(file = "Nepal_model_start_pop.rds")
  delta_ranking <- readRDS(file = "Nepal_delta_ranking.rds")
  #mass_cluster_freq_1 <- readRDS(file = "Nepal_cluster_freqs_1.rds")
  #mass_cluster_freq_2 <- readRDS(file = "Nepal_cluster_freqs_2.rds")
  #mass_cluster_freq_3 <- readRDS(file = "Nepal_cluster_freqs_3.rds")
  #mass_cluster_freq_4 <- readRDS(file = "Nepal_cluster_freqs_4.rds")
  #mass_cluster_freq_5 <- readRDS(file = "Nepal_cluster_freqs_5.rds")
  mass_cluster_freq_6 <- readRDS(file = "Nepal_cluster_freqs_6.rds")
  mass_cluster_freq_7 <- readRDS(file = "Nepal_cluster_freqs_7.rds")
  mass_cluster_freq_8 <- readRDS(file = "Nepal_cluster_freqs_8.rds")
  mass_cluster_freq_9 <- readRDS(file = "Nepal_cluster_freqs_9.rds")
  mass_cluster_freq_10 <- readRDS(file = "Nepal_cluster_freqs_10.rds")
  mass_cluster_freq_11 <- readRDS(file = "Nepal_cluster_freqs_11.rds")
  mass_cluster_freq_12 <- readRDS(file = "Nepal_cluster_freqs_12.rds")
  mass_cluster_freq_13 <- readRDS(file = "Nepal_cluster_freqs_13.rds")
  mass_cluster_freq_14 <- readRDS(file = "Nepal_cluster_freqs_14.rds")
  mass_VT <- readRDS(file = "Nepal_VT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
  output_filename <- "Nepal_VT_ggCaller_PopPUNK"
  
  dt <- 1/12
  #peripost_mass_cluster_freq <- data.frame("year" = 1:13, rbind(mass_cluster_freq_2, mass_cluster_freq_3, mass_cluster_freq_4, mass_cluster_freq_5, mass_cluster_freq_6, mass_cluster_freq_7, mass_cluster_freq_8, mass_cluster_freq_9,mass_cluster_freq_10, mass_cluster_freq_11, mass_cluster_freq_12, mass_cluster_freq_13,mass_cluster_freq_14))
  # now using 2009 as the start population (Nepal_cluster_freqs_5.rds)
  
  ceil_mass_NVT <- rep(0, mass_clusters)
  for (i in 1:mass_clusters){
    ceil_mass_NVT[i] <- ceiling(mean(as.integer(seq_clusters[seq_clusters$Cluster == unique(seq_clusters$Cluster)[i],"VT"]=="NVT")))
  }
  #ceil_mass_NVT[is.nan(ceil_mass_NVT)] <- 0
  mean_mass_VT_2001 <- rep(0, mass_clusters)
  for (i in 1:mass_clusters){
    mean_mass_VT_2001[i] <- (mean(as.integer(seq_clusters[seq_clusters$SeqYear == 2001 & seq_clusters$Cluster == unique(seq_clusters$Cluster)[i],"VT"]=="VT")))
  }
  mean_mass_VT_2001[is.nan(mean_mass_VT_2001)] <- 0
  
  NVT_mig <- rep(1/mass_clusters, mass_clusters) * (ceil_mass_NVT)
  VT_mig <- rep(1/mass_clusters, mass_clusters) * (1-ceil_mass_NVT)
  avg_cluster_freq <- data.frame(as.matrix(cbind(NVT_mig,VT_mig)))
  model_start_pop <- matrix(as.double(c(model_start_pop * (1-mean_mass_VT_2001),model_start_pop * mean_mass_VT_2001)), byrow = FALSE, nrow = mass_clusters, ncol = 2)
  peripost_mass_cluster_freq <- data.frame("year" = 1:9, rbind(mass_cluster_freq_6, mass_cluster_freq_7, mass_cluster_freq_8, mass_cluster_freq_9,mass_cluster_freq_10, mass_cluster_freq_11, mass_cluster_freq_12, mass_cluster_freq_13,mass_cluster_freq_14))
  
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  
  vacc_time <- 6
} else if(args[1] == "Navajo" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("Navajo_PP.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "Navajo_ggCaller_intermed_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  model_start_pop <- readRDS(file = "Navajo_model_start_pop.rds")
  delta_ranking <- readRDS(file = "Navajo_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "Navajo_cluster_freqs_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "Navajo_cluster_freqs_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "Navajo_cluster_freqs_3.rds")
  mass_cluster_freq_4 <- readRDS(file = "Navajo_cluster_freqs_4.rds")
  mass_cluster_freq_5 <- readRDS(file = "Navajo_cluster_freqs_5.rds")
  mass_cluster_freq_6 <- readRDS(file = "Navajo_cluster_freqs_6.rds")
  mass_cluster_freq_7 <- readRDS(file = "Navajo_cluster_freqs_7.rds")
  mass_cluster_freq_8 <- readRDS(file = "Navajo_cluster_freqs_8.rds")
  mass_cluster_freq_9 <- readRDS(file = "Navajo_cluster_freqs_9.rds")
  mass_cluster_freq_10 <- readRDS(file = "Navajo_cluster_freqs_10.rds")
  mass_cluster_freq_11 <- readRDS(file = "Navajo_cluster_freqs_11.rds")
  mass_cluster_freq_12 <- readRDS(file = "Navajo_cluster_freqs_12.rds")
  mass_cluster_freq_13 <- readRDS(file = "Navajo_cluster_freqs_13.rds")
  mass_cluster_freq_14 <- readRDS(file = "Navajo_cluster_freqs_14.rds")
  mass_cluster_freq_15 <- readRDS(file = "Navajo_cluster_freqs_15.rds")
  mass_VT <- readRDS(file = "Navajo_VT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
  output_filename <- "Navajo_VT_ggCaller_PopPUNK"
  
  dt <- 1/12
  peripost_mass_cluster_freq <- data.frame("year" = 1:14, rbind(mass_cluster_freq_2,mass_cluster_freq_3,mass_cluster_freq_4,mass_cluster_freq_5,mass_cluster_freq_6, mass_cluster_freq_7, mass_cluster_freq_8, mass_cluster_freq_9,mass_cluster_freq_10, mass_cluster_freq_11, mass_cluster_freq_12, mass_cluster_freq_13,mass_cluster_freq_14,mass_cluster_freq_15))
  
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  
  vacc_time1 <- 2 # trying vacc time =2 instead of 5 as before
  vacc_time2 <- 12
  
  ceil_mass_NVT1 <- rep(0, mass_clusters)
  for (i in 1:mass_clusters){
    ceil_mass_NVT1[i] <- ceiling(mean(as.integer(seq_clusters[seq_clusters$Cluster == unique(seq_clusters$Cluster)[i],"VT"]=="NVT")))
  }
  ceil_mass_NVT2 <- rep(0, mass_clusters)
  for (i in 1:mass_clusters){
    ceil_mass_NVT2[i] <- ceiling(mean(as.integer(seq_clusters[seq_clusters$Cluster == unique(seq_clusters$Cluster)[i],"VT_PCV13"]=="NVT")))
  }
  #ceil_mass_NVT[is.nan(ceil_mass_NVT)] <- 0
  mean_mass_VT1_start <- rep(0, mass_clusters)
  for (i in 1:mass_clusters){
    mean_mass_VT1_start[i] <- (mean(as.integer(seq_clusters[seq_clusters$SeqYear == 1999 & seq_clusters$Cluster == unique(seq_clusters$Cluster)[i],"VT"]=="VT")))
  }
  mean_mass_VT1_start[is.nan(mean_mass_VT1_start)] <- 0
  mean_mass_VT2_start <- rep(0, mass_clusters)
  for (i in 1:mass_clusters){
    mean_mass_VT2_start[i] <- (mean(as.integer(seq_clusters[seq_clusters$SeqYear == 1999 & seq_clusters$Cluster == unique(seq_clusters$Cluster)[i],"VT_PCV13"]=="VT")))
  }
  mean_mass_VT2_start[is.nan(mean_mass_VT2_start)] <- 0
  mean_mass_VT2_start <- mean_mass_VT2_start - mean_mass_VT1_start
  
  NVT_mig <- rep(1/mass_clusters, mass_clusters) * (ceil_mass_NVT2)
  VT_mig <- rep(1/mass_clusters, mass_clusters) * (1-ceil_mass_NVT1) 
  VT2_mig <- rep(1/mass_clusters, mass_clusters) * (1-ceil_mass_NVT2) * ceil_mass_NVT1
  avg_cluster_freq <- data.frame(as.matrix(cbind(NVT_mig,VT_mig,VT2_mig)))
  model_start_pop <- matrix(as.double(c(model_start_pop * (1-(mean_mass_VT1_start + mean_mass_VT2_start)),model_start_pop * mean_mass_VT1_start, model_start_pop * mean_mass_VT2_start)), byrow = FALSE, nrow = mass_clusters, ncol = 3)
} else if(args[1] == "UK" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("UK_PP.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "UK_ggCaller_intermed_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  model_start_pop <- readRDS(file = "UK_model_start_pop.rds")
  delta_ranking <- readRDS(file = "UK_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "UK_cluster_freqs_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "UK_cluster_freqs_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "UK_cluster_freqs_3.rds")
  mass_cluster_freq_4 <- readRDS(file = "UK_cluster_freqs_4.rds")
  mass_cluster_freq_5 <- readRDS(file = "UK_cluster_freqs_5.rds")
  mass_cluster_freq_6 <- readRDS(file = "UK_cluster_freqs_6.rds")
  mass_cluster_freq_7 <- readRDS(file = "UK_cluster_freqs_7.rds")
  mass_VT <- readRDS(file = "UK_VT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
  output_filename <- "UK_VT_ggCaller_PopPUNK"
  
  dt <- 1/12
  peripost_mass_cluster_freq <- data.frame("year" = 1:6, rbind(mass_cluster_freq_2,mass_cluster_freq_3,mass_cluster_freq_4,mass_cluster_freq_5,mass_cluster_freq_6, mass_cluster_freq_7))
  
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  
  vacc_time1 <- 0
  vacc_time2 <- 4
  
  ceil_mass_NVT1 <- rep(0, mass_clusters)
  for (i in 1:mass_clusters){
    ceil_mass_NVT1[i] <- ceiling(mean(as.integer(seq_clusters[seq_clusters$Cluster == unique(seq_clusters$Cluster)[i],"VT_PCV7"]=="NVT")))
  }
  ceil_mass_NVT2 <- rep(0, mass_clusters)
  for (i in 1:mass_clusters){
    ceil_mass_NVT2[i] <- ceiling(mean(as.integer(seq_clusters[seq_clusters$Cluster == unique(seq_clusters$Cluster)[i],"VT_PCV13"]=="NVT")))
  }
  #ceil_mass_NVT[is.nan(ceil_mass_NVT)] <- 0
  mean_mass_VT1_start <- rep(0, mass_clusters)
  for (i in 1:mass_clusters){
    mean_mass_VT1_start[i] <- (mean(as.integer(seq_clusters[seq_clusters$Time == 2006 & seq_clusters$Cluster == unique(seq_clusters$Cluster)[i],"VT_PCV7"]=="VT")))
  }
  mean_mass_VT1_start[is.nan(mean_mass_VT1_start)] <- 0
  mean_mass_VT2_start <- rep(0, mass_clusters)
  for (i in 1:mass_clusters){
    mean_mass_VT2_start[i] <- (mean(as.integer(seq_clusters[seq_clusters$Time == 2006 & seq_clusters$Cluster == unique(seq_clusters$Cluster)[i],"VT_PCV13"]=="VT")))
  }
  mean_mass_VT2_start[is.nan(mean_mass_VT2_start)] <- 0
  mean_mass_VT2_start <- mean_mass_VT2_start - mean_mass_VT1_start
  
  NVT_mig <- rep(1/mass_clusters, mass_clusters) * (ceil_mass_NVT2)
  VT_mig <- rep(1/mass_clusters, mass_clusters) * (1-ceil_mass_NVT1) 
  VT2_mig <- rep(1/mass_clusters, mass_clusters) * (1-ceil_mass_NVT2) * ceil_mass_NVT1
  avg_cluster_freq <- data.frame(as.matrix(cbind(NVT_mig,VT_mig,VT2_mig)))
  model_start_pop <- matrix(as.double(c(model_start_pop * (1-(mean_mass_VT1_start + mean_mass_VT2_start)),model_start_pop * mean_mass_VT1_start, model_start_pop * mean_mass_VT2_start)), byrow = FALSE, nrow = mass_clusters, ncol = 3)
}

threads_total <- 1
if(length(args)>=3){
  print(paste("Setting the number of total threads to ", args[3]))
  threads_total <- as.integer(args[3])
}
worker_nodes <- 1
if(length(args)>=4){
  print(paste("Setting the number of workers to ", args[4]))
  worker_nodes <- as.integer(args[4])
}
stoch_run <- FALSE
if(length(args)>=5 & args[5]=="stoch"){
  stoch_run <- TRUE
}

fitting_mass_data <- mcstate::particle_filter_data(data = peripost_mass_cluster_freq,
                                                   time = "year",
                                                   rate = 1 / dt,
                                                   initial_time = 0)

det_filter <- particle_deterministic$new(data = fitting_mass_data,
                                         model = WF,
                                         compare = combined_compare)


# process data with particle filter:
#dt <- 1/36 # we assume that the generation time of Strep. pneumo is 1 month
# we have data from 2001, 2004, 2007, so we want 3 (years) * 12 (months) = 36 updates in-between

#peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(mass_cluster_freq_2, mass_cluster_freq_3))
#names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))

#fitting_mass_data <- mcstate::particle_filter_data(data = peripost_mass_cluster_freq,
#                                                   time = "year",
#                                                   rate = 1 / dt,
#                                                   initial_time = 0)

#det_filter <- particle_deterministic$new(data = fitting_mass_data,
#                                         model = WF,
#                                         compare = combined_compare)

# Using MCMC to infer parameters
pmcmc_sigma_f <- mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0, max = 1)
#pmcmc_sigma_w <- 0
pmcmc_sigma_w <- -1000
pmcmc_prop_f <- mcstate::pmcmc_parameter("prop_f", 0.2, min = 0, max = 1)
pmcmc_m <- mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1)
pmcmc_v <- mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)
species_no <- mass_clusters
no_clusters <- mass_clusters
gene_no <- nrow(intermed_gene_presence_absence_consensus_matrix)

Pop_ini <- data.frame(model_start_pop)
Pop_eq <- rowSums(model_start_pop)
Genotypes <- intermed_gene_presence_absence_consensus_matrix
#Genotypes_vec <- as.double(c(intermed_gene_presence_absence_consensus_matrix))

capacity <- sum(model_start_pop)
delta <- delta_ranking
vT1 = c(0,1,0)
vT2 = c(0,0,1)
migVec <- data.frame(avg_cluster_freq)

complex_params = list(species_no = species_no, Pop_ini = Pop_ini, Pop_eq = Pop_eq, Genotypes = intermed_gene_presence_absence_consensus[-1,-1], capacity = capacity, delta = delta, vT1 = vT1,  vT2 = vT2, gene_no = gene_no, vacc_time1 = vacc_time1, vacc_time2 = vacc_time2, dt = dt, sigma_w = pmcmc_sigma_w, migVec = (migVec))
#complex_params <- c(Pop_ini, Pop_eq, Genotypes, capacity, delta, species_no, gene_no, vacc_time, dt, migVec,vT)

make_transform <- function(m) {
  function(theta) {
    as_double_mtx <- function(x){
      sapply(x,as.double)
    }
    c(lapply(m, as_double_mtx), as.list(theta))
  }
}

# as_double_mtx nice idea but can't find the function in parallelisation

take_list <- function(x){
  print(x$Genotype)
}

transform <- function() {
  make_transform(complex_params)}

#make_transform <- function(p) {
#  function(theta){
#    as_double_mtx <- function(x){
#      sapply(x,as.double)
#    }
#    c(list(species_no = species_no, Pop_ini = Pop_ini, Pop_eq = Pop_eq, Genotypes = as.matrix(sapply(intermed_gene_presence_absence_consensus[-1,-1], as.double), nrow = 1770, ncol = 55), capacity = capacity, delta = delta, vT = vT, gene_no = gene_no, vacc_time = vacc_time, dt = dt, sigma_w = pmcmc_sigma_w, migVec = as.matrix(migVec)), as.list(theta))
#  }
#}

#transformed_params <<- make_transform(complex_params)



index <- function(info) {
  list(run = c(sum_clust = info$index$Pop_tot),
       state = c(Pop = info$index$Pop))
}

#WF_model <- WF$new(pars = list(), time = 0, n_particles = 1L)
#index(WF$info())


proposal_matrix <- diag(0.1,4) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
# here, all parameters are proposed independently. 
# think about this, this might not actually be true
#mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
#mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.1432, min = 0, max = 1), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 01), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", 0.125, min = 0, max = 1), mcstate::pmcmc_parameter("m", -4, min = -1000, max = 0), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
mcmc_pars$initial()
#mcmc_pars$model(mcmc_pars$initial())

#WF$public_methods$has_openmp()
det_filter <- particle_deterministic$new(data = fitting_mass_data,
                                         model = WF,
                                         index = index,
                                         compare = combined_compare)

n_steps <- 5
n_burnin <- 0


control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE, 
  save_trajectories = TRUE,
  progress = TRUE,
  adaptive_proposal = TRUE,
  n_chains = 1)
det_pmcmc_run <- mcstate::pmcmc(mcmc_pars, det_filter, control = control)

n_steps <- 1000
n_burnin <- 0


control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE, 
  save_trajectories = TRUE,
  progress = TRUE,
  adaptive_proposal = TRUE,
  n_chains =4, n_workers = 4,
  n_threads_total = 4)

#n_chains = 8, n_workers = 8,
#n_threads_total = 8

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
print("det_mcmc_1 final log likelihood")
processed_chains$probabilities[nrow(processed_chains$probabilities),2]
print("det_mcmc_1 mean log likelihood")
mean(processed_chains$probabilities[,2])
det_proposal_matrix <- cov(processed_chains$pars)

det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", parameter_mean_hpd[1], min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", parameter_mean_hpd[2], min = 0, max = 1),mcstate::pmcmc_parameter("m", parameter_mean_hpd[3], min = -1000, max = 0), mcstate::pmcmc_parameter("v", parameter_mean_hpd[4], min = 0, max = 1)), det_proposal_matrix, make_transform(complex_params))

det_filter <- particle_deterministic$new(data = fitting_mass_data,
                                         model = WF,
                                         index = index,
                                         compare = combined_compare)

n_steps <- 5
n_burnin <- 0


control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE, 
  save_trajectories = TRUE,
  progress = TRUE,
  adaptive_proposal = TRUE,
  n_chains = 1)
det_pmcmc_run <- mcstate::pmcmc(det_mcmc_pars, det_filter, control = control)

n_steps <- 20000
n_burnin <- 0


control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE, 
  save_trajectories = TRUE,
  progress = TRUE,
  adaptive_proposal = TRUE,
  n_chains = 4, n_workers = 4, n_threads_total = 4)
det_pmcmc_run2 <- mcstate::pmcmc(det_mcmc_pars, det_filter, control = control)
processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run2, burnin = 2000, thin = 1)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
print(parameter_mean_hpd)
par(mfrow = c(1,1))

det_mcmc2 <- coda::as.mcmc(cbind(det_pmcmc_run2$probabilities, det_pmcmc_run2$pars))
pdf(file = paste(output_filename,"det_mcmc2.pdf",sep = "_"),   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 12)
plot(det_mcmc2)
dev.off()
print("det_mcmc_2 final log likelihood")
processed_chains$probabilities[nrow(processed_chains$probabilities),2]
print("det_mcmc_2 mean log likelihood")
mean(processed_chains$probabilities[,2])

saveRDS(det_pmcmc_run2, paste(output_filename, "_det_pmcmc_run2.rds", sep = ""))

if(stoch_run == TRUE){
  det_proposal_matrix <- cov(processed_chains$pars)
  #det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0.075, max = 0.22), mcstate::pmcmc_parameter("sigma_w", 0.05, min = 0.000001, max = 0.0749), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 0.2), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 0.5)), det_proposal_matrix, make_transform(complex_params))
  det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", parameter_mean_hpd[1], min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", parameter_mean_hpd[2], min = 0, max = 1),mcstate::pmcmc_parameter("m", parameter_mean_hpd[3], min = -1000, max = 0), mcstate::pmcmc_parameter("v", parameter_mean_hpd[4], min = 0, max = 1)), det_proposal_matrix, make_transform(complex_params))
  
  
  filter <- mcstate::particle_filter$new(data = fitting_mass_data,
                                         model = WF,
                                         n_particles = 6,
                                         index = index,
                                         compare = combined_compare,
                                         n_threads = 8)
  
  n_steps <- 2
  n_burnin <- 0
  
  control <- mcstate::pmcmc_control(n_steps, n_chains = 1, n_workers = 1,save_state = TRUE,
                                    save_trajectories = TRUE,
                                    progress = TRUE,
                                    n_threads_total = 1)
  pmcmc_run <- mcstate::pmcmc(det_mcmc_pars, filter, control = control)
  
  filter <- mcstate::particle_filter$new(data = fitting_mass_data,
                                         model = WF,
                                         n_particles = 96,
                                         index = index,
                                         compare = combined_compare,
                                         n_threads = 8)
  n_steps <- 2000
  n_burnin <- 0
  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE, 
    save_trajectories = TRUE,
    progress = TRUE, 
    n_chains = 4, n_workers = worker_nodes, 
    n_threads_total = threads_total)
  #control <- mcstate::pmcmc_control(
  #  n_steps,
  #  save_state = TRUE, 
  #  save_trajectories = TRUE,
  #  progress = TRUE, 
  #  n_chains = 2)
  
  stoch_pmcmc_run2 <- mcstate::pmcmc(det_mcmc_pars, filter, control = control)
  par(mfrow = c(1,1))
  
  stoch_mcmc2 <- coda::as.mcmc(cbind(stoch_pmcmc_run2$probabilities, stoch_pmcmc_run2$pars))
  
  
  pdf(file = paste(output_filename,"stoch_mcmc2.pdf",sep = "_"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 12)
  plot(stoch_mcmc2)
  dev.off()
  
  processed_chains <- mcstate::pmcmc_thin(stoch_pmcmc_run2, burnin = 200, thin = 1)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  print(parameter_mean_hpd)
  print("stoch_mcmc_2 final log likelihood")
  print(processed_chains$probabilities[nrow(processed_chains$probabilities),2])
  print("stoch_mcmc_2 mean log likelihood")
  print(mean(processed_chains$probabilities[,2]))
  
  saveRDS(stoch_pmcmc_run2, paste(output_filename, "_stoch_pmcmc_run2.rds", sep = ""))
}