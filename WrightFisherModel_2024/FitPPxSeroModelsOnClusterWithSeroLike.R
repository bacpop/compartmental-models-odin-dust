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
WF <- odin.dust::odin_dust("NFDS_Model_PPxSeroWithSero.R")

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
  data_size <- sum(unlist(observed[as.character(1:(length(unlist(observed))-4))]))/2
  model_size = sum(unlist(state[-1, , drop = TRUE]))/2
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
# data pre-processing in NFDS_COGtriangles_ManualSeqClust.Rmd
if(args[1] == "ggCaller" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("PopPUNK_clusters.rds")
  sero_no = length(unique(seq_clusters$Serotype))
  intermed_gene_presence_absence_consensus <- readRDS(file = "ggCPP_intermed_gene_presence_absence_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  #model_start_pop <- readRDS("PP_mass_cluster_freq_1_sero.rds")
  #model_start_pop <- model_start_pop / 133 * 15708
  model_start_pop <- readRDS("PPsero_startpop5.rds") # this is drawing from the start population and very low frequency for all other occurring combinations
  #model_start_pop <- readRDS("PP_mass_cluster_freq_1_sero.rds") # try using data directly
  #model_start_pop <- readRDS(file = "PPsero_startpop4.rds")
  #model_start_pop <- readRDS(file = "PPsero_startpop.rds")
  delta_ranking <- readRDS(file = "ggC_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "PP_mass_cluster_freq_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "PP_mass_cluster_freq_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "PP_mass_cluster_freq_3.rds")
  mass_VT <- readRDS(file = "SeroVT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- readRDS(file = "PPsero_mig.rds")
  dt <- 1/36
  
  vacc_time <- 0
  PP_mass_sero_freq_1 <- readRDS(file = "PP_mass_sero_freq_1.rds")
  PP_mass_sero_freq_2 <- readRDS(file = "PP_mass_sero_freq_2.rds")
  PP_mass_sero_freq_3 <- readRDS(file = "PP_mass_sero_freq_3.rds")
  
  peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(c(mass_cluster_freq_2, PP_mass_sero_freq_2), c(mass_cluster_freq_3,PP_mass_sero_freq_3)))
  names(peripost_mass_cluster_freq) <- c("year", as.character(c(1:mass_clusters, (mass_clusters+1):(mass_clusters+sero_no))))
  output_filename <- "PPxSeroWithSero_ggCaller_PopPUNK"
} else if(args[1] == "COGtriangles" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("PopPUNK_clusters.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "PP_intermed_gene_presence_absence_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  model_start_pop <- readRDS(file = "PP_model_start_pop.rds")
  delta_ranking <- readRDS(file = "delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "PP_mass_cluster_freq_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "PP_mass_cluster_freq_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "PP_mass_cluster_freq_3.rds")
  mass_VT <- readRDS(file = "PP_mass_VT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
  dt <- 1/36
  peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(mass_cluster_freq_2, mass_cluster_freq_3))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  vacc_time <- 0
  output_filename <- "4param_COGtriangles_PopPUNK"
} else if(args[1] == "ggCaller" & args[2] == "manualSeqClusters"){
  seq_clusters <- readRDS("Mass_Samples_accCodes.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "ggC_intermed_gene_presence_absence_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  model_start_pop <- readRDS(file = "model_start_pop.rds")
  delta_ranking <- readRDS(file = "ggC_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "mass_cluster_freq_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "mass_cluster_freq_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "mass_cluster_freq_3.rds")
  mass_VT <- readRDS(file = "mass_VT.rds")
  mass_clusters <- length(unique(seq_clusters$SequenceCluster))
  avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
  dt <- 1/36
  peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(mass_cluster_freq_2, mass_cluster_freq_3))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  vacc_time <- 0
  output_filename <- "4param_ggCaller_manSeqClusters"
} else if(args[1] == "COGtriangles" & args[2] == "manualSeqClusters"){
  seq_clusters <- readRDS("Mass_Samples_accCodes.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "intermed_gene_presence_absence_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  model_start_pop <- readRDS(file = "model_start_pop.rds")
  delta_ranking <- readRDS(file = "delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "mass_cluster_freq_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "mass_cluster_freq_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "mass_cluster_freq_3.rds")
  mass_VT <- readRDS(file = "mass_VT.rds")
  mass_clusters <- length(unique(seq_clusters$SequenceCluster))
  avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
  dt <- 1/36
  peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(mass_cluster_freq_2, mass_cluster_freq_3))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  vacc_time <- 0
  output_filename <- "4param_COGtriangles_manSeqClusters"
} else if(args[1] == "Navajo" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("Navajo_PP.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "Navajo_ggCaller_intermed_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  
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
  
  mass_clusters <- length(unique(seq_clusters$Cluster))
  sero_no = length(unique(seq_clusters$Serotype))
  model_start_pop <- readRDS(file = "Navajo_PPsero_startpop3.rds")
  #model_start_pop <- readRDS(file = "Navajo_PPsero_startpop.rds")
  
  mass_VT <- readRDS(file = "Navajo_SeroVT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- readRDS(file = "Navajo_PPsero_mig.rds")
  dt <- 1/12
  peripost_mass_cluster_freq <- data.frame("year" = 1:12, rbind(mass_cluster_freq_4,mass_cluster_freq_5,mass_cluster_freq_6, mass_cluster_freq_7,mass_cluster_freq_8,mass_cluster_freq_9, mass_cluster_freq_10, mass_cluster_freq_11,mass_cluster_freq_12, mass_cluster_freq_13, mass_cluster_freq_14,mass_cluster_freq_15))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  vacc_time <- 2
  output_filename <- "Navajo_PPxSero_ggCaller_PopPUNK"
} else if(args[1] == "UK" & args[2] == "PopPUNK"){
  seq_clusters <- readRDS("UK_PP.rds")
  intermed_gene_presence_absence_consensus <- readRDS(file = "UK_ggCaller_intermed_consensus.rds")
  intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
  mass_clusters <- length(unique(seq_clusters$Cluster))
  delta_ranking <- readRDS(file = "UK_delta_ranking.rds")
  mass_cluster_freq_1 <- readRDS(file = "UK_cluster_freqs_1.rds")
  mass_cluster_freq_2 <- readRDS(file = "UK_cluster_freqs_2.rds")
  mass_cluster_freq_3 <- readRDS(file = "UK_cluster_freqs_3.rds")
  mass_cluster_freq_4 <- readRDS(file = "UK_cluster_freqs_4.rds")
  mass_cluster_freq_5 <- readRDS(file = "UK_cluster_freqs_5.rds")
  mass_cluster_freq_6 <- readRDS(file = "UK_cluster_freqs_6.rds")
  mass_cluster_freq_7 <- readRDS(file = "UK_cluster_freqs_7.rds")
  mass_VT <- readRDS(file = "UK_SeroVT.rds")
  
  avg_cluster_freq <- readRDS(file = "UK_PPsero_mig.rds")
  output_filename <- "UK_PPxSero_ggCaller_PopPUNK"
  
  dt <- 1/12
  peripost_mass_cluster_freq <- data.frame("year" = 1:6, rbind(mass_cluster_freq_2,mass_cluster_freq_3,mass_cluster_freq_4,mass_cluster_freq_5,mass_cluster_freq_6, mass_cluster_freq_7))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  vacc_time <- 0
  sero_no = length(unique(seq_clusters$Serotype))
  #model_start_pop <- readRDS(file = "UK_PPsero_startpop.rds")
  model_start_pop <- readRDS(file = "UK_PPsero_startpop3.rds")
}  else if(args[1] == "Nepal" & args[2] == "PopPUNK"){
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
  
  mass_clusters <- length(unique(seq_clusters$Cluster))
  sero_no = length(unique(seq_clusters$Serotype))
  
  model_start_pop <- readRDS(file = "Nepal_PPsero_startpop.rds")
  
  mass_VT <- readRDS(file = "Nepal_SeroVT.rds")
  mass_clusters <- length(unique(seq_clusters$Cluster))
  avg_cluster_freq <- readRDS(file = "Nepal_PPsero_mig.rds")
  dt <- 1/12
  peripost_mass_cluster_freq <- data.frame("year" = 1:9, rbind(mass_cluster_freq_8,mass_cluster_freq_9, mass_cluster_freq_10, mass_cluster_freq_11,mass_cluster_freq_12, mass_cluster_freq_13, mass_cluster_freq_14,mass_cluster_freq_15,mass_cluster_freq_16))
  names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
  vacc_time <- 6
  output_filename <- "Nepal_PPxSero_ggCaller_PopPUNK"
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
  list(run = c(sum_clust = c(info$index$Pop_tot, info$index$Pop_sero)),
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



