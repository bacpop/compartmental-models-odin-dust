# idea randomly "mutate" 0-1-vector that I can map to real gene vector (based on common absence and presence, i.e. treat genes that always occur together as one gene)

#### run this part 

### Likelihood
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

#setwd("/nfs/research/jlees/leonie/WF_fitting_2024/run3")
### load model
WF <- odin.dust::odin_dust("NFDS_Model_FindGenes.R")

###########################################

### try binary ga

fitting_closure_max_binary <- function(all_other_params, data1, data2){
  null_fit_dfoptim_fl <- function(fit_params){
    rnd_vect_full <- fit_params
    
    all_other_params$delta_bool = rnd_vect_full
    WFmodel_ggCPP <- WF$new(pars = all_other_params,
                            time = 1,
                            n_particles = 1L,
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
    #simMeanggCPP2 <- rowMeans(WFmodel_ggCPP$run(36)[-1,])
    #simMeanggCPP3 <- rowMeans(WFmodel_ggCPP$run(72)[-1,])
    #combined_compare(simMeanggCPP2,data1) + combined_compare(simMeanggCPP3,data2) 
    sim1 <- WFmodel_ggCPP$run(37)
    sim2 <- WFmodel_ggCPP$run(73)
    combined_compare(sim1,data1) + combined_compare(sim2,data2)
    #- combined_compare(x[,1,37],data1) - combined_compare(x[,1,73],data2) 
  }
}


library(parallel)
#install.packages("doParallel",repos = "http://cran.us.r-project.org")
library(doParallel)
#install.packages("doSNOW",repos = "http://cran.us.r-project.org")
library(doSNOW)
#install.packages("GA",repos = "http://cran.us.r-project.org")
library(GA)

seq_clusters <- readRDS("PopPUNK_clusters.rds")
intermed_gene_presence_absence_consensus <- readRDS(file = "ggCPP_intermed_gene_presence_absence_consensus.rds")
intermed_gene_presence_absence_consensus_matrix <- sapply(intermed_gene_presence_absence_consensus[-1,-1],as.double)
model_start_pop <- readRDS(file = "PP_model_start_pop.rds")
delta_ranking <- readRDS(file = "ggC_delta_ranking.rds")
mass_cluster_freq_1 <- readRDS(file = "PP_mass_cluster_freq_1.rds")
mass_cluster_freq_2 <- readRDS(file = "PP_mass_cluster_freq_2.rds")
mass_cluster_freq_3 <- readRDS(file = "PP_mass_cluster_freq_3.rds")
mass_VT <- readRDS(file = "PP_mass_VT_mean.rds")
mass_VT <- readRDS(file = "PP_mass_VT.rds")
mass_clusters <- length(unique(seq_clusters$Cluster))
avg_cluster_freq <- rep(1/mass_clusters, mass_clusters)
output_filename <- "4param_ggCaller_PopPUNK"
# process data with particle filter:
dt <- 1/36 # we assume that the generation time of Strep. pneumo is 1 month
# we have data from 2001, 2004, 2007, so we want 3 (years) * 12 (months) = 36 updates in-between

peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(mass_cluster_freq_2, mass_cluster_freq_3))
names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))
vacc_time <- 0



FindGenes_ggCPP_params <- list(dt = 1/36, species_no = mass_clusters,  gene_no = nrow(intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(model_start_pop), Pop_eq = as.double(model_start_pop), capacity = sum(model_start_pop), Genotypes = intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.3090376, sigma_w = 0, prop_f = 1, m = 0.03104461, migVec = avg_cluster_freq, vaccTypes = mass_VT, v = 0.15977862, vacc_time = 0)

ga_fit_FindGenes_ggCPP_bin <- fitting_closure_max_binary(FindGenes_ggCPP_params, mass_cluster_freq_2, mass_cluster_freq_3)

#gann <- ga(type = "binary", nBits = 1770, fitness = ga_fit_FindGenes_ggCPP_bin, lower = rep(0, 1770), upper = rep(1,1770), 
#           elitism = 50, maxiter = 30, popSize = 3000, run = 10, pcrossover = 0.8, pmutation = 0.3, crossover = gabin_spCrossover, mutation = gabin_raMutation, parallel = 16)
# better seed: 123  Mean = -685.3895 | Best = -453.6379

gann <- ga(type = "binary", nBits = 1770, fitness = ga_fit_FindGenes_ggCPP_bin, lower = rep(0, 1770), upper = rep(1,1770), 
           elitism = 1000, maxiter = 100, popSize = 100000, run = 10, pcrossover = 0.8, pmutation = 0.1, crossover = gabin_spCrossover, mutation = gabin_raMutation, parallel = 48)

#plot(gann)
print(summary(gann))
print(sum(gann@solution))
#[1] 905

saveRDS(gann, "gann.rds")


pdf(file = "gann_plot.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 12)
plot(gann)
dev.off()
