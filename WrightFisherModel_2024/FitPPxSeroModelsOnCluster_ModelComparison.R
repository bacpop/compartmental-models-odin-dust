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
WF <- odin.dust::odin_dust("NFDS_Model_PPxSero.R")

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

#combined_compare <- function(state, observed, pars = NULL) {
#  result <- 0
#data_size <- sum(unlist(observed))
#  data_size <- sum(unlist(observed[as.character(1:(length(unlist(observed))-4))]))
#  model_size = sum(unlist(state[-1, , drop = TRUE]))
#  exp_noise <- 1e6

#  for (i in 1:(length(unlist(observed))-4)){ 
#    state_name <- paste("sum_clust", i, sep = "")
#    if (is.na(observed[[as.character(i)]])) {
#      #Creates vector of zeros in ll with same length, if no data
#      ll_obs <- numeric(length( state[state_name, , drop = TRUE]))
#    } else {
#      lambda <-  state[state_name, , drop = TRUE]/model_size * data_size + rexp(n = length( state[state_name, , drop = TRUE]/model_size * data_size), rate = exp_noise)
#      ll_obs <- dpois(x = observed[[as.character(i)]], lambda = lambda, log = TRUE)
#    }

#    result <- result + ll_obs
#  }
#  result
#}

combined_compare <- function(state, observed, pars = NULL) {
  result <- 0
  #data_size <- sum(unlist(observed))
  data_size <- sum(unlist(observed[as.character(1:(length(unlist(observed))-4))]))
  model_size = sum(unlist(state[-1, , drop = TRUE]))
  exp_noise <- 1e6
  data_vals <- unlist(observed[as.character(1:(length(unlist(observed))-4))])
  #model_vals <- state[-1, , drop = TRUE]
  model_vals <- rep(0, length(unlist(observed))-4)
  data_missing <- FALSE
  for (i in 1:(length(unlist(observed))-4)){ 
    state_name <- paste("sum_clust", i, sep = "")
    model_vals[i] <- state[state_name, , drop = TRUE]
    if (is.na(observed[[as.character(i)]])) {
      #Creates vector of zeros in ll with same length, if no data
      #ll_obs <- numeric(length( state[state_name, , drop = TRUE]))
      data_missing <- TRUE
    } 
  }
  models_vals_err <- model_vals + rexp(n = length(model_vals), rate = exp_noise)
  if(data_missing){
    ll_obs <- 0
  }
  else{
    ll_obs <- dmultinom(x = (data_vals), prob = models_vals_err/model_size, log = TRUE)   
  }
  result <- ll_obs
  #for (i in 1:(length(unlist(observed))-4)){ 
  #  state_name <- paste("sum_clust", i, sep = "")
  #  if (is.na(observed[[as.character(i)]])) {
  #    #Creates vector of zeros in ll with same length, if no data
  #    ll_obs <- numeric(length( state[state_name, , drop = TRUE]))
  #  } else {
  #lambda <-  state[state_name, , drop = TRUE]/model_size * data_size + rexp(n = length( state[state_name, , drop = TRUE]/model_size * data_size), rate = exp_noise)
  #ll_obs <- dpois(x = observed[[as.character(i)]], lambda = lambda, log = TRUE)
  #    ll_obs <- dmultinom(x = (data_vals), prob = model_vals/model_size, log = TRUE)
  #  }
  
  #  result <- result + ll_obs
  #}
  result
}
# data pre-processing in NFDS_COGtriangles_ManualSeqClust.Rmd
if(args[1] == "ggCaller" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("PopPUNK_clusters.rds")
  sero_no = length(unique(seq_clusters$Serotype))
  intermed_gene_presence_absence_consensus <- readRDS(file = "ggCPP_intermed_gene_presence_absence_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  #model_start_pop <- readRDS("PP_mass_cluster_freq_1_sero.rds")
  #model_start_pop <- model_start_pop / 133 * 15708
  model_start_pop <- readRDS("PPsero_startpop6.rds") 
  #model_start_pop <- readRDS("PP_mass_cluster_freq_1_sero.rds") # try using data directly
  #model_start_pop <- readRDS(file = "PPsero_startpop4.rds")
  #model_start_pop <- readRDS(file = "PPsero_startpop.rds")
  delta_ranking <- readRDS(file = "ggC_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "PP_mass_cluster_freq_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "PP_mass_cluster_freq_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "PP_mass_cluster_freq_3.rds")
  #mass_VT <- readRDS(file = "SeroVT.rds")
  mass_VT <- readRDS(file = "SeroVT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- readRDS(file = "PPsero_mig.rds")
  dt <- 1/36
  peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(mass_cluster_freq_2, mass_cluster_freq_3))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  vacc_time <- 0
  output_filename <- "PPxSero_ggCaller_PopPUNK"
} else if(args[1] == "Nepal" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("Nepal_PP.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "Nepal_ggCaller_intermed_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  
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
  mass_cluster_freq_15 <- readRDS(file = "Nepal_cluster_freqs_15.rds")
  mass_cluster_freq_16 <- readRDS(file = "Nepal_cluster_freqs_16.rds")
  mass_cluster_freq_17 <- readRDS(file = "Nepal_cluster_freqs_17.rds")
  
  mass_clusters <- length(unique(seq_clusters$Cluster))
  sero_no = length(unique(seq_clusters$Serotype))
  
  model_start_pop <- readRDS(file = "Nepal_PPsero_startpop.rds")
  
  mass_VT <- readRDS(file = "Nepal_SeroVT.rds")
  mass_clusters <- length(unique(seq_clusters$GPSC))
  avg_cluster_freq <- readRDS(file = "Nepal_PPsero_mig.rds")
  dt <- 1/12
  peripost_mass_cluster_freq <- data.frame("year" = 1:10, rbind(mass_cluster_freq_8,mass_cluster_freq_9, mass_cluster_freq_10, mass_cluster_freq_11,mass_cluster_freq_12, mass_cluster_freq_13, mass_cluster_freq_14,mass_cluster_freq_15,mass_cluster_freq_16,mass_cluster_freq_17))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  vacc_time <- 6
  #cov_matrix <- readRDS("Nepal_cov_matrix.rds")
  output_filename <- "Nepal_PPxSero_ggCaller_PopPUNK"
} else if(args[1] == "Navajo" & args[2] == "PopPUNK"){
  # this is the test what happens when I only focus on the first vacc introduction in Navajo
  # for this, I will ignore data from 2011 and 2012
  seq_clusters <- readRDS("Navajo_PP.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "Navajo_ggCaller_intermed_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  
  delta_ranking <- readRDS(file = "Navajo_delta_ranking.rds") # this is as before, so does include data from 2011/12
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
  #mass_cluster_freq_14 <- readRDS(file = "Navajo_cluster_freqs_14.rds")
  #mass_cluster_freq_15 <- readRDS(file = "Navajo_cluster_freqs_15.rds")
  
  mass_clusters <- length(unique(seq_clusters$Cluster))
  sero_no = length(unique(seq_clusters$Serotype))
  
  #model_start_pop <- readRDS(file = "Navajo_PPsero_startpop.rds")
  model_start_pop <- readRDS(file = "Navajo_PPsero_startpop.rds")
  #mass_VT <- readRDS(file = "Navajo_SeroVT.rds")
  #PCV13_VTs <- rep(0,sero_no)
  #names(PCV13_VTs) <- unique(seq_clusters$Serotype)
  #PCV13_VTs[intersect(names(PCV13_VTs), c("1", "3", "4","5","6A","6B", "7F", "9V", "14", "18C", "19A", "19F", "23F"))] <- 1
  #vaccTypes2 <- unname(PCV13_VTs)
  
  # add 6A to PCV7 because there is strong cross-immunity btw PVC7 and 6A (4.Croucher, N. J. et al. Population genomics of post-vaccine changes in pneumococcal epidemiology. Nat. Genet. 45, 656–663 (2013).)
  #PCV7_VTs <- rep(0,sero_no)
  #names(PCV7_VTs) <- unique(seq_clusters$Serotype)
  #PCV7_VTs[intersect(PCV7_VTs,c("4", "6A","6B", "9V", "14", "18C", "19F", "23F"))] <- 1
  #vaccTypes1 <- unname(PCV7_VTs)
  mass_VT <- readRDS("Navajo_SeroVT.rds")
  #vaccTypes2 <- readRDS("Navajo_SeroVT_PCV13.rds")
  
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- readRDS(file = "Navajo_PPsero_mig.rds")
  dt <- 1/12
  peripost_mass_cluster_freq <- data.frame("year" = 1:12, rbind(mass_cluster_freq_2,mass_cluster_freq_3,mass_cluster_freq_4,mass_cluster_freq_5,mass_cluster_freq_6, mass_cluster_freq_7,mass_cluster_freq_8,mass_cluster_freq_9, mass_cluster_freq_10, mass_cluster_freq_11,mass_cluster_freq_12, mass_cluster_freq_13))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  vacc_time <- 2
  #vacc_time2 <- 12
  output_filename <- "Navajo_PPxSero_ggCaller_PopPUNK"
}

# number of parameters to fit in model
# 2 == null model, 3 == homogeneous, 4 = partial NFDS, 5 == heterogeneous
params_total <- 4
if(length(args)>=3){
  print(paste("Fitting the model with parameter number ", args[3]))
  params_total <- as.integer(args[3])
}

threads_total <- 1
if(length(args)>=4){
  print(paste("Setting the number of total threads to ", args[4]))
  threads_total <- as.integer(args[4])
}
worker_nodes <- 1
if(length(args)>=5){
  print(paste("Setting the number of workers to ", args[5]))
  worker_nodes <- as.integer(args[5])
}
stoch_run <- FALSE
if(length(args)>=6 & args[6]=="stoch"){
  stoch_run <- TRUE
}

# process data with particle filter:
#dt <- 1/36 # we assume that the generation time of Strep. pneumo is 1 month
# we have data from 2001, 2004, 2007, so we want 3 (years) * 12 (months) = 36 updates in-between

#peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(mass_cluster_freq_2, mass_cluster_freq_3))
#names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))

fitting_mass_data <- mcstate::particle_filter_data(data = peripost_mass_cluster_freq,
                                                   time = "year",
                                                   rate = 1 / dt,
                                                   initial_time = 0)

det_filter <- particle_deterministic$new(data = fitting_mass_data,
                                         model = WF,
                                         compare = combined_compare)

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

capacity <- sum(model_start_pop)
delta <- delta_ranking
vaccTypes <- mass_VT

migVec <- data.frame(avg_cluster_freq)

#complex_params <- c(Pop_ini, Pop_eq, Genotypes, capacity, delta, vaccTypes, species_no, gene_no, vacc_time, dt, migVec, pmcmc_sigma_w, sero_no)

### have to re-write the whole make_transform function because of the matrices. 
#dim(Pop_ini)

#make_transform <- function(p) {
#  function(theta){
#    c(list(Pop_ini = matrix(p[1:(nrow(Pop_ini) * ncol(Pop_ini))], nrow = mass_clusters, ncol = sero_no),
#           Pop_eq = p[((nrow(Pop_ini) * ncol(Pop_ini)) +1) : ((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq))],
#           Genotypes = matrix(p[(((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq))+ 1): (((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes))], nrow = gene_no, ncol = species_no),
#           capacity = p[(((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 1],
#           delta = p[((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2) : ((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no -1)],
#           vaccTypes = p[((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no) : (((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no) + sero_no -1)],
#           species_no = p[(((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no) + sero_no )],
#           gene_no = p[(((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no) + sero_no)+ 1],
#           vacc_time = p[(((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no) + sero_no) + 2],
#           dt = p[(((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no) + sero_no) + 3],
#           migVec = matrix(p[((((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no) + sero_no) + 4):((((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no) + sero_no) + 4 + (nrow(migVec) * ncol(migVec)) -1)], nrow = mass_clusters, ncol = sero_no),
#           sigma_w = p[((((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no) + sero_no) + 4 + (nrow(migVec) * ncol(migVec)))],
#           sero_no = p[((((((nrow(Pop_ini) * ncol(Pop_ini)) + length(Pop_eq)) + nrow(Genotypes) * ncol(Genotypes)) + 2 + gene_no) + sero_no) + 4 + (nrow(migVec) * ncol(migVec)) +1)]), as.list(theta))
#  }
#}
##################################

if(params_total == 2){

  complex_params = list(species_no = species_no, Pop_ini = Pop_ini, Pop_eq = Pop_eq, Genotypes = intermed_gene_presence_absence_consensus[-1,-1], capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, sigma_w = pmcmc_sigma_w, migVec = (migVec), sero_no = sero_no, sigma_f = -1000, prop_f = 1)
  #complex_params <- c(Pop_ini, Pop_eq, Genotypes, capacity, delta, species_no, gene_no, vacc_time, dt, migVec,vT)
  
  make_transform <- function(m) {
    function(theta) {
      as_double_mtx <- function(x){
        sapply(x,as.double)
      }
      c(lapply(m, as_double_mtx), as.list(theta))
    }
  }
  
  
  #proposal_matrix[1,1] <- exp(0.1)
  #proposal_matrix[3,3] <- exp(0.1)
  # here, all parameters are proposed independently. 
  # think about this, this might not actually be true
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0, max = 1), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  #= make_transform(c(Pop_ini, Pop_eq, Genotypes, capacity, delta, vaccTypes, species_no, gene_no, vacc_time)))
  #mcmc_pars$names()
  #mcmc_pars$model(mcmc_pars$initial())
  # read this: https://mrc-ide.github.io/mcstate/reference/pmcmc_parameters.html
  # it explains how to not fit all parameters but just the ones I want
  # non-scalar parameters have to be transformed for this.
  
  index <- function(info) {
    list(run = c(sum_clust = info$index$Pop_tot),
         state = c(Pop = info$index$Pop))
  }
  
  #WF_model <- WF$new(pars = list(), time = 0, n_particles = 1L)
  #index(WF$info())
  
  
  #proposal_matrix <- diag(0.1,4) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
  # here, all parameters are proposed independently. 
  # think about this, this might not actually be true
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.1432, min = 0, max = 1), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 01), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", 0.125, min = 0, max = 1), mcstate::pmcmc_parameter("m", -4, min = -1000, max = -2), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
 # mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", 0.125, min = 0, max = 1), mcstate::pmcmc_parameter("m", -4, min = -1000, max = 0), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  proposal_matrix <- diag(c(exp(1), 0.1))
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("m", runif(n=1, min=-10, max=0), min = -1000, max = 0), mcstate::pmcmc_parameter("v", runif(n=1, min=0, max=1), min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  
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
  pdf(file = paste(output_filename,"Null_det_mcmc1.pdf",sep = "_"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 12)
  plot(det_mcmc1)
  dev.off()
  print("det_mcmc_1 final log likelihood")
  processed_chains$probabilities[nrow(processed_chains$probabilities),2]
  print("det_mcmc_1 mean log likelihood")
  mean(processed_chains$probabilities[,2])
  det_proposal_matrix <- cov(processed_chains$pars)
  
  det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("m", parameter_mean_hpd[1], min = -1000, max = 0), mcstate::pmcmc_parameter("v", parameter_mean_hpd[2], min = 0, max = 1)), det_proposal_matrix, make_transform(complex_params))
  
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
  pdf(file = paste(output_filename,"Null_det_mcmc2.pdf",sep = "_"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 12)
  plot(det_mcmc2)
  dev.off()
  print("det_mcmc_2 final log likelihood")
  processed_chains$probabilities[nrow(processed_chains$probabilities),2]
  print("det_mcmc_2 mean log likelihood")
  mean(processed_chains$probabilities[,2])
  
  saveRDS(det_pmcmc_run2, paste(output_filename, "_Null_det_pmcmc_run2.rds", sep = ""))
} else if(params_total == 3){
  
  complex_params = list(species_no = species_no, Pop_ini = Pop_ini, Pop_eq = Pop_eq, Genotypes = intermed_gene_presence_absence_consensus[-1,-1], capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, sigma_w = pmcmc_sigma_w, migVec = (migVec), sero_no = sero_no, prop_f = 1)
  #complex_params <- c(Pop_ini, Pop_eq, Genotypes, capacity, delta, species_no, gene_no, vacc_time, dt, migVec,vT)
  
  make_transform <- function(m) {
    function(theta) {
      as_double_mtx <- function(x){
        sapply(x,as.double)
      }
      c(lapply(m, as_double_mtx), as.list(theta))
    }
  }
  
  
  #proposal_matrix <- diag(0.1,3) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
  #proposal_matrix[1,1] <- exp(0.1)
  #proposal_matrix[3,3] <- exp(0.1)
  # here, all parameters are proposed independently. 
  # think about this, this might not actually be true
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0, max = 1), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  #= make_transform(c(Pop_ini, Pop_eq, Genotypes, capacity, delta, vaccTypes, species_no, gene_no, vacc_time)))
  #mcmc_pars$names()
  #mcmc_pars$model(mcmc_pars$initial())
  # read this: https://mrc-ide.github.io/mcstate/reference/pmcmc_parameters.html
  # it explains how to not fit all parameters but just the ones I want
  # non-scalar parameters have to be transformed for this.
  
  index <- function(info) {
    list(run = c(sum_clust = info$index$Pop_tot),
         state = c(Pop = info$index$Pop))
  }
  
  #WF_model <- WF$new(pars = list(), time = 0, n_particles = 1L)
  #index(WF$info())
  
  
  #proposal_matrix <- diag(0.1,4) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
  # here, all parameters are proposed independently. 
  # think about this, this might not actually be true
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.1432, min = 0, max = 1), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 01), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", 0.125, min = 0, max = 1), mcstate::pmcmc_parameter("m", -4, min = -1000, max = -2), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("m", -4, min = -1000, max = 0), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  proposal_matrix <- diag(c(exp(1), exp(1), 0.1))
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", runif(n=1, min=-10, max=0), min = -1000, max = 0), mcstate::pmcmc_parameter("m", runif(n=1, min=-10, max=0), min = -1000, max = 0), mcstate::pmcmc_parameter("v", runif(n=1, min=0, max=1), min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  
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
  pdf(file = paste(output_filename,"3param_det_mcmc1.pdf",sep = "_"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 12)
  plot(det_mcmc1)
  dev.off()
  print("det_mcmc_1 final log likelihood")
  processed_chains$probabilities[nrow(processed_chains$probabilities),2]
  print("det_mcmc_1 mean log likelihood")
  mean(processed_chains$probabilities[,2])
  det_proposal_matrix <- cov(processed_chains$pars)
  
  det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", parameter_mean_hpd[1], min = -1000, max = 0), mcstate::pmcmc_parameter("m", parameter_mean_hpd[2], min = -1000, max = 0), mcstate::pmcmc_parameter("v", parameter_mean_hpd[3], min = 0, max = 1)), det_proposal_matrix, make_transform(complex_params))
  
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
  pdf(file = paste(output_filename,"3param_det_mcmc2.pdf",sep = "_"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 12)
  plot(det_mcmc2)
  dev.off()
  print("det_mcmc_2 final log likelihood")
  processed_chains$probabilities[nrow(processed_chains$probabilities),2]
  print("det_mcmc_2 mean log likelihood")
  mean(processed_chains$probabilities[,2])
  
  saveRDS(det_pmcmc_run2, paste(output_filename, "_3param_det_pmcmc_run2.rds", sep = ""))
} else if(params_total == 4){
  complex_params = list(species_no = species_no, Pop_ini = Pop_ini, Pop_eq = Pop_eq, Genotypes = intermed_gene_presence_absence_consensus[-1,-1], capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, sigma_w = pmcmc_sigma_w, migVec = (migVec), sero_no = sero_no)
  #complex_params <- c(Pop_ini, Pop_eq, Genotypes, capacity, delta, species_no, gene_no, vacc_time, dt, migVec,vT)
  
  make_transform <- function(m) {
    function(theta) {
      as_double_mtx <- function(x){
        sapply(x,as.double)
      }
      c(lapply(m, as_double_mtx), as.list(theta))
    }
  }
  
  
  proposal_matrix <- diag(0.1,4) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
  #proposal_matrix[1,1] <- exp(0.1)
  #proposal_matrix[3,3] <- exp(0.1)
  # here, all parameters are proposed independently. 
  # think about this, this might not actually be true
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0, max = 1), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  #= make_transform(c(Pop_ini, Pop_eq, Genotypes, capacity, delta, vaccTypes, species_no, gene_no, vacc_time)))
  #mcmc_pars$names()
  #mcmc_pars$model(mcmc_pars$initial())
  # read this: https://mrc-ide.github.io/mcstate/reference/pmcmc_parameters.html
  # it explains how to not fit all parameters but just the ones I want
  # non-scalar parameters have to be transformed for this.
  
  index <- function(info) {
    list(run = c(sum_clust = info$index$Pop_tot),
         state = c(Pop = info$index$Pop))
  }
  
  #WF_model <- WF$new(pars = list(), time = 0, n_particles = 1L)
  #index(WF$info())
  
  
  #proposal_matrix <- diag(0.1,4) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
  # here, all parameters are proposed independently. 
  # think about this, this might not actually be true
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.1432, min = 0, max = 1), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 01), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", 0.125, min = 0, max = 1), mcstate::pmcmc_parameter("m", -4, min = -1000, max = -2), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", 0.125, min = 0, max = 1), mcstate::pmcmc_parameter("m", -4, min = -1000, max = 0), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  proposal_matrix <- diag(c(exp(1), 0.1, exp(1), 0.1))
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", runif(n=1, min=-10, max=0), min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", runif(n=1, min=0, max=1), min = 0, max = 1), mcstate::pmcmc_parameter("m", runif(n=1, min=-10, max=0), min = -1000, max = 0), mcstate::pmcmc_parameter("v", runif(n=1, min=0, max=1), min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  
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
  #n_steps <- 500
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
  #processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run, burnin = 100, thin = 1)
  processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run, burnin = 250, thin = 1)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  print(parameter_mean_hpd)
  
  det_mcmc1 <- coda::as.mcmc(cbind(det_pmcmc_run$probabilities, det_pmcmc_run$pars))
  pdf(file = paste(output_filename,"4param_det_mcmc1.pdf",sep = "_"),   # The directory you want to save the file in
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
  #n_steps <- 5000
  n_burnin <- 0
  
  
  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE, 
    save_trajectories = TRUE,
    progress = TRUE,
    adaptive_proposal = TRUE,
    n_chains = 4, n_workers = 4, n_threads_total = 4)
  det_pmcmc_run2 <- mcstate::pmcmc(det_mcmc_pars, det_filter, control = control)
  #processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run2, burnin = 1000, thin = 1)
  processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run2, burnin = 2000, thin = 1)
  parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
  print(parameter_mean_hpd)
  par(mfrow = c(1,1))
  
  det_mcmc2 <- coda::as.mcmc(cbind(det_pmcmc_run2$probabilities, det_pmcmc_run2$pars))
  pdf(file = paste(output_filename,"_4param_det_mcmc2.pdf",sep = "_"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 12)
  plot(det_mcmc2)
  dev.off()
  print("det_mcmc_2 final log likelihood")
  processed_chains$probabilities[nrow(processed_chains$probabilities),2]
  print("det_mcmc_2 mean log likelihood")
  mean(processed_chains$probabilities[,2])
  
  saveRDS(det_pmcmc_run2, paste(output_filename, "_4param_det_pmcmc_run2.rds", sep = ""))
  
} else if(params_total == 5){
  WF <- odin.dust::odin_dust("NFDS_Model_PPxSero_5param.R")
  complex_params = list(species_no = species_no, Pop_ini = Pop_ini, Pop_eq = Pop_eq, Genotypes = intermed_gene_presence_absence_consensus[-1,-1], capacity = capacity, delta = delta, vaccTypes = vaccTypes, gene_no = gene_no, vacc_time = vacc_time, dt = dt, migVec = (migVec), sero_no = sero_no)
  #complex_params <- c(Pop_ini, Pop_eq, Genotypes, capacity, delta, species_no, gene_no, vacc_time, dt, migVec,vT)
  
  make_transform <- function(m) {
    function(theta) {
      as_double_mtx <- function(x){
        sapply(x,as.double)
      }
      c(lapply(m, as_double_mtx), as.list(theta))
    }
  }
  
  
  proposal_matrix <- diag(0.1,5) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
  #proposal_matrix[1,1] <- exp(0.1)
  #proposal_matrix[3,3] <- exp(0.1)
  # here, all parameters are proposed independently. 
  # think about this, this might not actually be true
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0, max = 1), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  #= make_transform(c(Pop_ini, Pop_eq, Genotypes, capacity, delta, vaccTypes, species_no, gene_no, vacc_time)))
  #mcmc_pars$names()
  #mcmc_pars$model(mcmc_pars$initial())
  # read this: https://mrc-ide.github.io/mcstate/reference/pmcmc_parameters.html
  # it explains how to not fit all parameters but just the ones I want
  # non-scalar parameters have to be transformed for this.
  
  index <- function(info) {
    list(run = c(sum_clust = info$index$Pop_tot),
         state = c(Pop = info$index$Pop))
  }
  
  #WF_model <- WF$new(pars = list(), time = 0, n_particles = 1L)
  #index(WF$info())
  
  
  #proposal_matrix <- diag(0.1,4) # the proposal matrix defines the covariance-variance matrix for a mult normal dist
  # here, all parameters are proposed independently. 
  # think about this, this might not actually be true
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(pmcmc_sigma_f, pmcmc_sigma_w, pmcmc_prop_f, pmcmc_m, pmcmc_v), proposal_matrix, transform)
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.1432, min = 0, max = 1), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 01), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", 0.125, min = 0, max = 1), mcstate::pmcmc_parameter("m", -4, min = -1000, max = -2), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", 0.125, min = 0, max = 1), mcstate::pmcmc_parameter("sigma_w", -0.597837, min = -1000, max = 0), mcstate::pmcmc_parameter("m", -4, min = -1000, max = 0), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  proposal_matrix <- diag(c(exp(1), 0.1, exp(1), exp(1), 0.1))
  #mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", runif(n=1, min=-3.4, max=0), min = -3.5, max = 0), mcstate::pmcmc_parameter("prop_f", runif(n=1, min=0, max=1), min = 0, max = 1), mcstate::pmcmc_parameter("sigma_w", runif(n=1, min=-1000, max=-3.6), min = -1000, max = -3.5), mcstate::pmcmc_parameter("m", runif(n=1, min=-10, max=0), min = -1000, max = 0), mcstate::pmcmc_parameter("v", runif(n=1, min=0, max=1), min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", runif(n=1, min=-10, max=0), min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", runif(n=1, min=0, max=1), min = 0, max = 1), mcstate::pmcmc_parameter("sigma_w", runif(n=1, min=-1000, max=0), min = -1000, max = -3.5), mcstate::pmcmc_parameter("m", runif(n=1, min=-10, max=0), min = -1000, max = 0), mcstate::pmcmc_parameter("v", runif(n=1, min=0, max=1), min = 0, max = 1)), proposal_matrix, make_transform(complex_params))
  
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
  pdf(file = paste(output_filename,"5param_det_mcmc1.pdf",sep = "_"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 12)
  plot(det_mcmc1)
  dev.off()
  print("det_mcmc_1 final log likelihood")
  print(processed_chains$probabilities[nrow(processed_chains$probabilities),2])
  print("det_mcmc_1 mean log likelihood")
  mean(processed_chains$probabilities[,2])
  print(det_proposal_matrix <- cov(processed_chains$pars))
  
  det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", parameter_mean_hpd[1], min = -1000, max = 0), mcstate::pmcmc_parameter("prop_f", parameter_mean_hpd[2], min = 0, max = 1), mcstate::pmcmc_parameter("sigma_w", parameter_mean_hpd[3], min = -1000, max = 0),mcstate::pmcmc_parameter("m", parameter_mean_hpd[4], min = -1000, max = 0), mcstate::pmcmc_parameter("v", parameter_mean_hpd[5], min = 0, max = 1)), det_proposal_matrix, make_transform(complex_params))
  
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
  pdf(file = paste(output_filename,"_5param_det_mcmc2.pdf",sep = "_"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 12)
  plot(det_mcmc2)
  dev.off()
  print("det_mcmc_2 final log likelihood")
  processed_chains$probabilities[nrow(processed_chains$probabilities),2]
  print("det_mcmc_2 mean log likelihood")
  mean(processed_chains$probabilities[,2])
  
  saveRDS(det_pmcmc_run2, paste(output_filename, "_5param_det_pmcmc_run2.rds", sep = ""))
  
} else{
  print(paste("Number of parameters chosen not okay:", params_total))
}




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



