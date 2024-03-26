# idea randomly "mutate" 0-1-vector that I can map to real gene vector (based on common absence and presence, i.e. treat genes that always occur together as one gene)

### find genes that always occur together
which(ggCPP_intermed_gene_presence_absence_consensus[2,]==1)

ggCPP_intermed_gene_presence_absence_consensus[-1,-1] <-  sapply(ggCPP_intermed_gene_presence_absence_consensus[-1,-1],as.integer)
grouped_genes_df <- data.frame(matrix(0,nrow = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1,ncol = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1))

for (i in 1:(nrow(ggCPP_intermed_gene_presence_absence_consensus)-3)){
  #local_rowSum <- apply(abs(sapply(data.frame(matrix(rep(ggCPP_intermed_gene_presence_absence_consensus[(i+1),-1],nrow(ggCPP_intermed_gene_presence_absence_consensus)-i),nrow = nrow(ggCPP_intermed_gene_presence_absence_consensus)-i, byrow = TRUE)), as.integer) - sapply(ggCPP_intermed_gene_presence_absence_consensus[(i+2):(nrow(ggCPP_intermed_gene_presence_absence_consensus)),-1], as.integer)), MARGIN = 1, max)
  local_vec <- as.integer(ggCPP_intermed_gene_presence_absence_consensus[(i+1),-1])
  local_diff <- t(apply(sapply(ggCPP_intermed_gene_presence_absence_consensus[(i+2):(nrow(ggCPP_intermed_gene_presence_absence_consensus)),-1], as.integer), 1, function(x) x - local_vec))
  local_rowSum <- apply(abs(local_diff), MARGIN = 1, max)
  grouped_genes_df[i,(i+1):(nrow(ggCPP_intermed_gene_presence_absence_consensus)-1)] <- 1 - local_rowSum
}
local_vec <- as.integer(ggCPP_intermed_gene_presence_absence_consensus[nrow(ggCPP_intermed_gene_presence_absence_consensus)-1,-1])
local_diff <- as.integer(ggCPP_intermed_gene_presence_absence_consensus[(nrow(ggCPP_intermed_gene_presence_absence_consensus)),-1]) - local_vec
local_rowSum <- max(abs(local_diff))
grouped_genes_df[nrow(ggCPP_intermed_gene_presence_absence_consensus)-1,nrow(ggCPP_intermed_gene_presence_absence_consensus)-1] <- 1 - local_rowSum

diag(grouped_genes_df) <- 0
sum(grouped_genes_df)
# [1] 1876
# I found 1876 gene pairs.

grouped_ind <- which(matrix(grepl("1", unlist(grouped_genes_df)), dim(grouped_genes_df)), arr.ind = TRUE)
# these are the indexes

group_gene_cl <- rep(0,length(unique(grouped_ind[,2])))
names_group_gene_cl <- rep(0,length(unique(grouped_ind[,2])))
j <- 0
for (i in 1:nrow(grouped_ind)){
  if(!(grouped_ind[i,2] %in% names_group_gene_cl)){
    j <- j + 1
    group_gene_cl[j] <- grouped_ind[i,1]
    names_group_gene_cl[j] <- grouped_ind[i,2]
  }
}
names(group_gene_cl) <- names_group_gene_cl
length((group_gene_cl))
# 520 genes that are in a cluster with some other gene(s)
length(unique(group_gene_cl))
# 160 group gene clusters
# unname

# So of 1774 genes, I can remove 520 when I am mutating the vector because they should always have the same value as their "gene group leader"

# create vector of length 1174 - 520 with randomly distributed 0s and 1s
rnd_vect <- rbinom(nrow(ggCPP_intermed_gene_presence_absence_consensus)-1 -length((group_gene_cl)) ,1,0.1)
# genes that are not part of any gene group (except gene group leaders)
unique_genes <- setdiff(1:(nrow(ggCPP_intermed_gene_presence_absence_consensus)-1), names(group_gene_cl))
rnd_vect_full <- rep(0, (nrow(ggCPP_intermed_gene_presence_absence_consensus)-1))
rnd_vect_full[unique_genes] <- rnd_vect
rnd_vect_full[as.integer(names(group_gene_cl))] <- rnd_vect_full[unname(group_gene_cl)]
# create random vector of length of all unique genes and the gene group leaders
# then expand to full length (1774 genes)
# then set all the members of gene groups to the value of their "group leader"

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

### load model
WF <- odin.dust::odin_dust("NFDS_Model_FindGenes.R")

# run ggCaller PopPUNK 3-param model multiple times for the different delta_bool vectors
FindNFDSgenes <- function(repeats = 100, frac = 0.1){
  best_vec <- rep(0, nrow(ggCPP_intermed_gene_presence_absence_consensus)-1)
  best_like <- -100000
  for (i in 1:repeats) {
    # create vector of length 1174 - 520 with randomly distributed 0s and 1s
    rnd_vect <- rbinom(nrow(ggCPP_intermed_gene_presence_absence_consensus)-1 -length((group_gene_cl)) ,1,frac)
    # genes that are not part of any gene group (except gene group leaders)
    unique_genes <- setdiff(1:(nrow(ggCPP_intermed_gene_presence_absence_consensus)-1), names(group_gene_cl))
    rnd_vect_full <- rep(0, (nrow(ggCPP_intermed_gene_presence_absence_consensus)-1))
    rnd_vect_full[unique_genes] <- rnd_vect
    rnd_vect_full[as.integer(names(group_gene_cl))] <- rnd_vect_full[unname(group_gene_cl)]
    
    
    ### ggCaller PopPUNK 3-param
    delta_bool <- rnd_vect_full
    params_ggCPP <- list(dt = 1/36, species_no = PP_mass_clusters,  gene_no = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(PP_model_start_pop), Pop_eq = as.double(PP_model_start_pop), capacity = sum(PP_model_start_pop), Genotypes = ggCPP_intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.01931485, sigma_w = 0, prop_f = 1, delta_bool = delta_bool, m = 0.03104461, migVec = PP_avg_cluster_freq, vaccTypes = PP_mass_VT, v = 0.15977862, vacc_time = 0)
    WFmodel_ggCPP <- WF$new(pars = params_ggCPP,
                            time = 1,
                            n_particles = 10L,
                            n_threads = 4L,
                            seed = 1L)
    simMeanggCPP2 <- rowMeans(WFmodel_ggCPP$run(36)[-1,])
    simMeanggCPP3 <- rowMeans(WFmodel_ggCPP$run(72)[-1,])
    
    # ggC PopPUNK
    local_like <- combined_compare(simMeanggCPP2,PP_mass_cluster_freq_2) + combined_compare(simMeanggCPP3,PP_mass_cluster_freq_3)
    if(local_like > best_like){
      best_like <- local_like
      best_vec <- rnd_vect_full
    }
    #print(local_like)
  }
  print(best_like)
  best_vec
}

best_100_perc_vect <- FindNFDSgenes(repeats = 100, frac = 1)
# -316.8468
best_90_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.9)
# -309.0821
best_80_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.8)
# -311.1709
best_70_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.7)
# -310.0685
best_60_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.6)
# -313.4829
best_50_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.5)
# -320.779
best_40_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.4)
# -327.5759
best_35_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.35)
# -327.4471
best_30_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.3)
#  -339.4924
best_25_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.25)
# -343.4345
best_20_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.2)
# -348.0353
best_15_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.15)
# -358.3352
best_10_perc_vec <- FindNFDSgenes(repeats = 100, frac = 0.15)
# -354.3937
best_5_perc_vec <- FindNFDSgenes(repeats = 100, frac = 0.05)
# -376.391


### more repeats (doesn't seem to bring a qualitative improvement)
best_90_perc_vect <- FindNFDSgenes(repeats = 500, frac = 0.9)
# -308.6588
best_80_perc_vect <- FindNFDSgenes(repeats = 500, frac = 0.8)
# -310.1969
best_70_perc_vect <- FindNFDSgenes(repeats = 500, frac = 0.7)
# -308.3332
best_60_perc_vect <- FindNFDSgenes(repeats = 500, frac = 0.6)
# -310.633