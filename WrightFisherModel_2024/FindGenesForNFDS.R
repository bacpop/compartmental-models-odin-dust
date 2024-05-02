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
  
  best_vec_df <- data.frame(matrix(0, nrow = repeats, ncol = nrow(ggCPP_intermed_gene_presence_absence_consensus)))
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
    #params_ggCPP <- list(dt = 1/36, species_no = PP_mass_clusters,  gene_no = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(PP_model_start_pop), Pop_eq = as.double(PP_model_start_pop), capacity = sum(PP_model_start_pop), Genotypes = ggCPP_intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.01931485, sigma_w = 0, prop_f = 1, delta_bool = delta_bool, m = 0.03104461, migVec = PP_avg_cluster_freq, vaccTypes = PP_mass_VT, v = 0.15977862, vacc_time = 0)
    #params_ggCPP <- list(dt = 1/36, species_no = PP_mass_clusters,  gene_no = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(PP_model_start_pop), Pop_eq = as.double(PP_model_start_pop), capacity = sum(PP_model_start_pop), Genotypes = ggCPP_intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.0386297, sigma_w = 0, prop_f = 1, delta_bool = delta_bool, m = 0.03104461, migVec = PP_avg_cluster_freq, vaccTypes = PP_mass_VT, v = 0.15977862, vacc_time = 0)
    #params_ggCPP <- list(dt = 1/36, species_no = PP_mass_clusters,  gene_no = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(PP_model_start_pop), Pop_eq = as.double(PP_model_start_pop), capacity = sum(PP_model_start_pop), Genotypes = ggCPP_intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.0772594, sigma_w = 0, prop_f = 1, delta_bool = delta_bool, m = 0.03104461, migVec = PP_avg_cluster_freq, vaccTypes = PP_mass_VT, v = 0.15977862, vacc_time = 0)
    #params_ggCPP <- list(dt = 1/36, species_no = PP_mass_clusters,  gene_no = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(PP_model_start_pop), Pop_eq = as.double(PP_model_start_pop), capacity = sum(PP_model_start_pop), Genotypes = ggCPP_intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.1545188, sigma_w = 0, prop_f = 1, delta_bool = delta_bool, m = 0.03104461, migVec = PP_avg_cluster_freq, vaccTypes = PP_mass_VT, v = 0.15977862, vacc_time = 0)
    params_ggCPP <- list(dt = 1/36, species_no = PP_mass_clusters,  gene_no = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(PP_model_start_pop), Pop_eq = as.double(PP_model_start_pop), capacity = sum(PP_model_start_pop), Genotypes = ggCPP_intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.3090376, sigma_w = 0, prop_f = 1, delta_bool = delta_bool, m = 0.03104461, migVec = PP_avg_cluster_freq, vaccTypes = PP_mass_VT, v = 0.15977862, vacc_time = 0)
    #params_ggCPP <- list(dt = 1/36, species_no = PP_mass_clusters,  gene_no = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(PP_model_start_pop), Pop_eq = as.double(PP_model_start_pop), capacity = sum(PP_model_start_pop), Genotypes = ggCPP_intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.6180752, sigma_w = 0, prop_f = 1, delta_bool = delta_bool, m = 0.03104461, migVec = PP_avg_cluster_freq, vaccTypes = PP_mass_VT, v = 0.15977862, vacc_time = 0)
    #params_ggCPP <- list(dt = 1/36, species_no = PP_mass_clusters,  gene_no = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(PP_model_start_pop), Pop_eq = as.double(PP_model_start_pop), capacity = sum(PP_model_start_pop), Genotypes = ggCPP_intermed_gene_presence_absence_consensus_matrix, sigma_f = 1.081632, sigma_w = 0, prop_f = 1, delta_bool = delta_bool, m = 0.03104461, migVec = PP_avg_cluster_freq, vaccTypes = PP_mass_VT, v = 0.15977862, vacc_time = 0)
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
    best_vec_df[i,1] <- local_like
    best_vec_df[i,-1] <- rnd_vect_full
    #print(local_like)
  }
  print(best_like)
  #best_vec
  best_vec_df
}

# first values with fitted parameters from 3-param model
# second values from same parameter set but sigma_f twice as large (0.0386297)
# third: doubled sigma_f again (0.07725940)
# fourth: doubled sigma_f again (0.1545188)
# fifth: doubled sigma_f again (0.3090376)
# sixth: doubled sigma_f again (0.6180752)
best_100_perc_vect <- FindNFDSgenes(repeats = 100, frac = 1)
# -316.8468
# -336.6054
# -377.6104
# -485.8535
# -1177.413
# -1710.235
best_90_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.9)
# -309.0821
# -325.0176
# -355.2654
# -382.8643
# -1247.047
best_80_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.8)
# -311.1709
# -318.3955
# -343.2562
# -366.5913
# -1180.972
# -1114.129
best_70_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.7)
# -310.0685
# -312.7839
# -336.6817
# -360.9399
# -985.3221
# -1084.834
best_60_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.6)
# -313.4829
# -304.5797
# -329.0323
# -357.7246
# -486.8475
# -1053.929
best_50_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.5)
# -320.779
# -303.6713
# -322.5524
# -347.8738
# -400.3433
# -724.7756
best_40_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.4)
# -327.5759
# -303.1324
# -313.384
# -331.0309
# -359.3175
# -602.6219
best_35_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.35)
# -327.4471
# -307.2633
# -308.0027
# -324.1988
# -338.6673
# -477.3307
best_30_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.3)
#  -339.4924
# -308.8178
# -306.5803
# -319.527
# -343.9595
# -407.3833
best_25_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.25)
# -343.4345
# -313.2131
# -304.1711
# -307.3205
# -338.3931
# -371.6217
best_20_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.2)
# -348.0353
# -327.1084
# -308.3151
# -312.0316
# -323.8821
# -344.4585
best_15_perc_vect <- FindNFDSgenes(repeats = 100, frac = 0.15)
# -358.3352
# -337.041
# -304.6377
# -303.4799
# -299.9564
# -323.5808
best_10_perc_vec <- FindNFDSgenes(repeats = 100, frac = 0.10)
# -354.3937
# -349.9868
# -329.6862
# -306.7359
# -288.273
# -304.1316
best_5_perc_vec <- FindNFDSgenes(repeats = 100, frac = 0.05)
# -376.391
# -371.9671
# -349.9106
# -319.9057
# -312.8814
# -306.8226
# -301.8122 (with sigma_f = 1.081632 = 0.6180752 * 1.75)

### more repeats (doesn't seem to bring a qualitative improvement)
best_90_perc_vect <- FindNFDSgenes(repeats = 500, frac = 0.9)
# -308.6588
best_80_perc_vect <- FindNFDSgenes(repeats = 500, frac = 0.8)
# -310.1969
best_70_perc_vect <- FindNFDSgenes(repeats = 500, frac = 0.7)
# -308.3332
best_60_perc_vect <- FindNFDSgenes(repeats = 500, frac = 0.6)
# -310.633


# changed the FindNFDSgenes function slightly to produce data frame
best_10_perc_vec_df <- FindNFDSgenes(repeats = 1000, frac = 0.10)

# 10% best likelihood
which(best_10_perc_vec_df$X1>(max(best_10_perc_vec_df$X1) + (min(best_10_perc_vec_df$X1) - max(best_10_perc_vec_df$X1))/10))
best_10_perc_vec_df$X1[which(best_10_perc_vec_df$X1>(max(best_10_perc_vec_df$X1) + (min(best_10_perc_vec_df$X1) - max(best_10_perc_vec_df$X1))/10))]
colMeans(best_10_perc_vec_df[which(best_10_perc_vec_df$X1>(max(best_10_perc_vec_df$X1) + (min(best_10_perc_vec_df$X1) - max(best_10_perc_vec_df$X1))/10)),])
which(colMeans(best_10_perc_vec_df[which(best_10_perc_vec_df$X1>(max(best_10_perc_vec_df$X1) + (min(best_10_perc_vec_df$X1) - max(best_10_perc_vec_df$X1))/10)),])>0.2)
# not even one is present in more than 20% of the vectors. looks like they are just randomly turned on/off
plot(colMeans(best_10_perc_vec_df[which(best_10_perc_vec_df$X1>(max(best_10_perc_vec_df$X1) + (min(best_10_perc_vec_df$X1) - max(best_10_perc_vec_df$X1))/10)),-1]))

which(colMeans(best_10_perc_vec_df[which(best_10_perc_vec_df$X1>(max(best_10_perc_vec_df$X1) + (min(best_10_perc_vec_df$X1) - max(best_10_perc_vec_df$X1))/10)),-1])>0.2)

opt_genes <- rep(0, (nrow(ggCPP_intermed_gene_presence_absence_consensus)-1))
opt_genes[which(colMeans(best_10_perc_vec_df[which(best_10_perc_vec_df$X1>(max(best_10_perc_vec_df$X1) + (min(best_10_perc_vec_df$X1) - max(best_10_perc_vec_df$X1))/10)),-1])>0.2)] <- 1

delta_bool <- opt_genes
params_ggCPP <- list(dt = 1/36, species_no = PP_mass_clusters,  gene_no = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(PP_model_start_pop), Pop_eq = as.double(PP_model_start_pop), capacity = sum(PP_model_start_pop), Genotypes = ggCPP_intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.3090376, sigma_w = 0, prop_f = 1, delta_bool = delta_bool, m = 0.03104461, migVec = PP_avg_cluster_freq, vaccTypes = PP_mass_VT, v = 0.15977862, vacc_time = 0)
WFmodel_ggCPP <- WF$new(pars = params_ggCPP,
                        time = 1,
                        n_particles = 10L,
                        n_threads = 4L,
                        seed = 1L)
simMeanggCPP2 <- rowMeans(WFmodel_ggCPP$run(36)[-1,])
simMeanggCPP3 <- rowMeans(WFmodel_ggCPP$run(72)[-1,])

# ggC PopPUNK
local_like <- combined_compare(simMeanggCPP2,PP_mass_cluster_freq_2) + combined_compare(simMeanggCPP3,PP_mass_cluster_freq_3)
local_like

#saveRDS(opt_genes,"opt_genes_v1.rds")
opt_genes1 <- readRDS("opt_genes_v1.rds")
#saveRDS(opt_genes,"opt_genes_v2.rds")
opt_genes2 <- readRDS("opt_genes_v2.rds")
#saveRDS(opt_genes,"opt_genes_v3.rds")
opt_genes3 <- readRDS("opt_genes_v3.rds")
#saveRDS(opt_genes,"opt_genes_v4.rds")
opt_genes4 <- readRDS("opt_genes_v4.rds")
#saveRDS(opt_genes,"opt_genes_v5.rds")
opt_genes5 <- readRDS("opt_genes_v5.rds")

plot(opt_genes4[order(ggC_delta_data3)])

tail(opt_genes4[ggC_delta_ranking])

which((opt_genes1 + opt_genes2 + opt_genes)>0)

plot(opt_genes1+opt_genes2+opt_genes3+opt_genes4+opt_genes5)

combined_opt <- opt_genes1+opt_genes2+opt_genes3+opt_genes4+opt_genes5
new_best_vec <- rep(0, (nrow(ggCPP_intermed_gene_presence_absence_consensus)-1))
new_best_vec[which(combined_opt>1)] <- 1
# [1] -266.0292
plot(new_best_vec[order(ggC_delta_data3)])
# does not seem to correlate at all with the gene frequencies

NFDSgenes_df <- data.frame(matrix(0, nrow = 8,ncol = 15))
colnames(NFDSgenes_df)[1] <- "sigma_f"
NFDSgenes_df[2:7,1] <- c(0.01931485, 0.0386297, 0.07725940, 0.1545188, 0.3090376, 0.6180752)
NFDSgenes_df[1,-1] <- c(1,0.9,0.8,0.7,0.6,0.5,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.05)
rownames(NFDSgenes_df)[1] <- "frac"
NFDSgenes_df[-1,2] <- c(-316.8468, -336.6054, -377.6104, -485.8535, -1633.545, -1710.235, NA)
NFDSgenes_df[-1,3] <- c(-309.0821, -325.0176, -355.2654, -382.8643, -1177.413, -1247.047, NA)
NFDSgenes_df[-1,4] <- c(-311.1709, -318.3955, -343.2562, -366.5913, -1180.972, -1114.129, NA)
NFDSgenes_df[-1,5] <- c(-310.0685, -312.7839, -336.6817, -360.9399, -985.3221, -1084.834, NA)
NFDSgenes_df[-1,6] <- c(-313.4829, -304.5797, -329.0323, -357.7246, -486.8475, -1053.929, NA)
NFDSgenes_df[-1,7] <- c(-320.779, -303.6713, -322.5524, -347.8738, -400.3433, -724.7756, NA)
NFDSgenes_df[-1,8] <- c(-327.5759, -303.1324, -313.384, -331.0309, -359.3175, -602.6219, NA)
NFDSgenes_df[-1,9] <- c(-327.4471, -307.2633, -308.0027, -324.1988, -338.6673, -477.3307, NA)
NFDSgenes_df[-1,10] <- c(-339.4924, -308.8178, -306.5803, -319.527, -343.9595, -407.3833, NA)
NFDSgenes_df[-1,11] <- c(-343.4345, -313.2131, -304.1711, -307.3205, -338.3931, -371.6217, NA)
NFDSgenes_df[-1,12] <- c(-348.0353, -327.1084, -308.3151, -312.0316, -323.8821, -344.4585, NA)
NFDSgenes_df[-1,13] <- c(-358.3352, -337.041, -304.6377, -303.4799, -299.9564, -323.5808, NA)
NFDSgenes_df[-1,14] <- c(-354.3937, -349.9868, -329.6862, -306.7359, -288.273, -304.1316, NA)
NFDSgenes_df[-1,15] <- c(-376.391, -371.9671, -349.9106, -319.9057, -312.8814, -306.8226, NA)

matplot(NFDSgenes_df$sigma_f[2:7], NFDSgenes_df[2:7,-1], type = "l", lty = 1, 
        col = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#b15928", "#6a3d9a", "#cab2d6", "#fdbf6f", "#fb9a99", "#e31a1c"),
        ylab = "Likelihood", main = "Likelihood for different fractions of genes under NFDS and various sigma_f values",
        ylim = c(-600,-250), log="x", xlab = "sigma_f values",xaxt='n')
axis(side = 1, at = NFDSgenes_df$sigma_f[2:7], labels = NFDSgenes_df$sigma_f[2:7])
legend("bottomleft", legend = NFDSgenes_df[1,-1],
       col = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#b15928", "#6a3d9a", "#cab2d6", "#fdbf6f", "#fb9a99", "#e31a1c"),
       lty = 1)





### Try finding optimal gene set using simulated annealing
library(optimization)
start_vect <- rbinom(nrow(ggCPP_intermed_gene_presence_absence_consensus)-1 -length((group_gene_cl)) ,1, 0.1)

fitting_closure <- function(all_other_params, data1, data2){
  null_fit_dfoptim <- function(fit_params){
    unique_genes <- setdiff(1:(nrow(ggCPP_intermed_gene_presence_absence_consensus)-1), names(group_gene_cl))
    rnd_vect_full <- rep(0, (nrow(ggCPP_intermed_gene_presence_absence_consensus)-1))
    rnd_vect_full[unique_genes] <- fit_params
    rnd_vect_full[as.integer(names(group_gene_cl))] <- rnd_vect_full[unname(group_gene_cl)]
    
    params$delta_bool = rnd_vect_full
    WFmodel_ggCPP <- WF$new(pars = params,
                                time = 1,
                                n_particles = 10L,
                                n_threads = 4L,
                                seed = 1L)
    n_particles <- 10L
    n_times <- 73
    x <- array(NA, dim = c(WFmodel_ggCPP$info()$len, n_particles, n_times))
    
    for (t in seq_len(n_times)) {
      x[ , , t] <- WFmodel_ggCPP$run(t)
    }
    time <- x[1, 1, ]
    x <- x[-1, , ]
    - combined_compare(x[,1,37],data1) - combined_compare(x[,1,73],data2) 
  }
}


FindGenes_ggCPP_params <- list(dt = 1/36, species_no = PP_mass_clusters,  gene_no = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1, Pop_ini = as.double(PP_model_start_pop), Pop_eq = as.double(PP_model_start_pop), capacity = sum(PP_model_start_pop), Genotypes = ggCPP_intermed_gene_presence_absence_consensus_matrix, sigma_f = 0.3090376, sigma_w = 0, prop_f = 1, m = 0.03104461, migVec = PP_avg_cluster_freq, vaccTypes = PP_mass_VT, v = 0.15977862, vacc_time = 0)

fit_FindGenes_ggCPP <- fitting_closure(FindGenes_ggCPP_params, PP_mass_cluster_freq_2, PP_mass_cluster_freq_3)

FindGenes_sa <- optim_sa(fun = fit_FindGenes_ggCPP, start = start_vect, maximization = TRUE, trace = TRUE,lower=rep(0,length(start_vect)), upper=rep(1,length(start_vect)))
plot(FindGenes_sa$par)
length(which(FindGenes_sa$par>0.6))
#[1] 130
length(which(FindGenes_sa$par>0.6))/length(FindGenes_sa$par)
#[1] 0.1036683

# that is not really a surprise because I was using the sigma_f that is optimal at 10% NFDS genes, according to my previous investigation
# but good to corroborate that I guess
# and the plot really is quite interesting because while the values are not just 0 or 1, most genes cluster at 0.2 and lower or 0.8 and higher

simAnn_results <- data.frame(matrix(0, nrow = nrow(ggCPP_intermed_gene_presence_absence_consensus)-1 -length((group_gene_cl)), ncol = 10))
simAnn_values <- rep(0, 10)

nextfun <- function(x) rbinom(nrow(ggCPP_intermed_gene_presence_absence_consensus)-1 -length((group_gene_cl)) ,1,runif(1))

for (i in 1:10) {
  FindGenes_sa2 <- optim(fn=fit_FindGenes_ggCPP, par=start_vect, gr=nextfun, method="SANN", 
                         control=list(maxit=1000,fnscale=1,trace=10))
  simAnn_results[,i] <- FindGenes_sa2$par
  simAnn_values[i] <- FindGenes_sa2$value
}

#FindGenes_sa2 <- optim(fn=fit_FindGenes_ggCPP, par=start_vect, gr=nextfun, method="SANN", 
#                       control=list(maxit=10,fnscale=1,trace=10))
plot(rowSums(simAnn_results))
plot(simAnn_results[,1])
points(simAnn_results[,2], col = "red")
colMeans(simAnn_results)

new_start_vec <- rep(0,length(start_vect))
new_start_vec[which(rowSums(simAnn_results)>=2)] <- 1
fit_FindGenes_ggCPP(new_start_vec)
FindGenes_sa3 <- optim(fn=fit_FindGenes_ggCPP, par=new_start_vec, gr=nextfun, method="SANN", 
                       control=list(maxit=1000,fnscale=1,trace=10))
plot(FindGenes_sa3$par + new_start_vec)
