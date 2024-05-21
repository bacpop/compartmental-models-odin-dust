### Loading packages
# install.packages("drat") # -- if you don't have drat installed
# drat:::add("ncov-ic")
# install.packages("odin.dust")
library(odin.dust)


# import data
seq_clusters <- readRDS("PopPUNK_clusters.rds")
intermed_gene_presence_absence_consensus <- readRDS(file = "ggCPP_intermed_gene_presence_absence_consensus.rds")
intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
model_start_pop <- readRDS(file = "PP_model_start_pop.rds")
#delta_ranking <- readRDS(file = "ggC_delta_ranking.rds")
mass_cluster_freq_1 <- readRDS(file = "PP_mass_cluster_freq_1.rds")
mass_cluster_freq_2 <- readRDS(file = "PP_mass_cluster_freq_2.rds")
mass_cluster_freq_3 <- readRDS(file = "PP_mass_cluster_freq_3.rds")
mass_VT <- readRDS(file = "PP_mass_VT.rds")
mass_clusters <- length(unique(seq_clusters$Cluster))
avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
output_filename <- "ggCaller_PopPUNK"

FindGenes_ggCPP_params <- list(dt = 1/36, species_no = mass_clusters,  gene_no = nrow(intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(model_start_pop), Pop_eq = as.double(model_start_pop), capacity = sum(model_start_pop), Genotypes = intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.3090376, sigma_w = 0, prop_f = 1, m = 0.03104461, migVec = avg_cluster_freq, vaccTypes = mass_VT, v = 0.15977862, vacc_time = 0)


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
  data_size <- sum(observed)
  model_size = sum(state)
  
  for (i in 1:length(observed)){
    result <- result + ll_pois(observed[i], state[i]/model_size * data_size)
  }
  result
}

### load model
WF <- odin.dust::odin_dust("NFDS_Model_FindGenes.R")

### try using genetic / evolutionary algorithms for finding best genes
#install.packages("GA", repos = 'https://cran.ma.imperial.ac.uk/')
library(GA)



decode2 <- function(x)
{ 
  x <- round(x)         
  return(x)
}



fitting_closure_max_decode <- function(all_other_params, data1, data2){
  null_fit_dfoptim_fl <- function(fit_params){
    fit_params <- decode2(fit_params)
    rnd_vect_full <- fit_params
    
    all_other_params$delta_bool = rnd_vect_full
    WFmodel_ggCPP <- WF$new(pars = all_other_params,
                            time = 1,
                            n_particles = 10L,
                            n_threads = 4L,
                            seed = 1L)
    #n_particles <- 10L
    #n_times <- 73
    #x <- array(NA, dim = c(WFmodel_ggCPP$info()$len, n_particles, n_times))
    
    #for (t in seq_len(n_times)) {
    #  x[ , , t] <- WFmodel_ggCPP$run(t)
    #}
    #time <- x[1, 1, ]
    #x <- x[-1, , ]
    simMeanggCPP2 <- rowMeans(WFmodel_ggCPP$run(36)[-1,])
    simMeanggCPP3 <- rowMeans(WFmodel_ggCPP$run(72)[-1,])
    combined_compare(simMeanggCPP2,data1) + combined_compare(simMeanggCPP3,data2) 
    #- combined_compare(x[,1,37],data1) - combined_compare(x[,1,73],data2) 
  }
}

gene_number <-  nrow(intermed_gene_presence_absence_consensus)-1

ga_fit_FindGenes_ggCPP_dec <- fitting_closure_max_decode(FindGenes_ggCPP_params, mass_cluster_freq_2, mass_cluster_freq_3)
gann <- ga(type = "real-valued", fitness = ga_fit_FindGenes_ggCPP_dec, lower = rep(0, gene_number), upper = rep(1,gene_number), 
           seed = 123, elitism = 40, maxiter = 20, popSize = 80, run = 30)
summary(gann)
as.vector(t(apply(gann@solution, 1, decode2)))
sum(as.vector(t(apply(gann@solution, 1, decode2))))/gene_number
saveRDS(gann,"gann.rds")
