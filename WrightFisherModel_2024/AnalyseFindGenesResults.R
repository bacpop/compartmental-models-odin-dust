library(readr)
output_ggCmanSeqCluster_FindGenes <- read_csv("~/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC/2024_03_12/output_ggCmanSeqCluster_FindGenes.txt")

View(output_ggCmanSeqCluster_FindGenes)
colnames(output_ggCmanSeqCluster_FindGenes) <- output_ggCmanSeqCluster_FindGenes[10,1]
output_ggCmanSeqCluster_FindGenes <- output_ggCmanSeqCluster_FindGenes[-(1:10),]

likelihoods_FindGenes_df <- data.frame(matrix(0, ncol = 4, nrow = 50))
#output_ggCmanSeqCluster_FindGenes[9,1]
colnames(likelihoods_FindGenes_df) <- c("ggCmanSeqCl","COGmanSeqCl","ggCPP","COGPP")

for (i in 1:49) {
  likelihoods_FindGenes_df[i,1] <- strsplit(output_ggCmanSeqCluster_FindGenes[9 * i,1]$`[1] "ggCaller_manSeqClusters"`, "\\[1\\]")[[1]][2]
}

#sort.list(likelihoods_FindGenes_df$ggCmanSeqCl) # should be decreasing=TRUE but it does not seem to recognise the minus sign
sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]


# Find the genes that are in these categories:
intermed_gene_presence_absence_consensus <- readRDS(file = "ggC_intermed_gene_presence_absence_consensus.rds")
# compute boolean gene vectors
n_groups <- ceiling((nrow(intermed_gene_presence_absence_consensus)-1)/25)
find_genes_df <- data.frame(matrix(0,nrow = nrow(intermed_gene_presence_absence_consensus)-1, ncol = 50))
for (i in 1:25) {
  find_genes_df[,i] <- rep(c(rep(0,24),1),(n_groups + 4))[i:(nrow(intermed_gene_presence_absence_consensus)-2 + i)]
  find_genes_df[,i+25] <- c(rep(0,(nrow(intermed_gene_presence_absence_consensus) - 1-n_groups)),rep(1,n_groups),rep(0,(nrow(intermed_gene_presence_absence_consensus) -n_groups)))[(1 + (i-1) * n_groups) : (nrow(intermed_gene_presence_absence_consensus)-1 + (i-1) * n_groups)]
}

rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]])
sum(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]])==1) # 571
sum(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]])==2) # 69

new_fitting_vec <- rep(0, nrow(intermed_gene_presence_absence_consensus)-1)
new_fitting_vec <- as.integer(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]])>=2)
saveRDS(new_fitting_vec, "ggCmanSeqCl_FindGenesResults1.rds")

# Fit for ggCmanSeq 30 44 46 21  5 49 47 22  1 43 (all genes that appear at least twice)

#sigma_f           m           v 
#0.176375633 0.009003682 0.061537212 
#[1] "det_mcmc_2 log likelihood"
#log_likelihood 
#-260.6447 
#[1] "det_mcmc_2 mean log likelihood"
#[1] -260.1266
### That is a relatively good, but not outstanding, fit for the 3-param model.
# there was one fit with a log-likelihood of -253.9187 (likelihoods_FindGenes_df$ggCmanSeqCl[30])

new_fitting_vec2 <- as.integer(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCmanSeqCl)[1:10]])>=1)
saveRDS(new_fitting_vec2, "ggCmanSeqCl_FindGenesResults2.rds")

#    sigma_f          m          v 
#0.06377337 0.02299571 0.18618946 
#[1] "det_mcmc_2 log likelihood"
#log_likelihood 
#-261.6746 
#[1] "det_mcmc_2 mean log likelihood"
#[1] -258.4692

### ggCPP FindGenes version
output_ggCPP_FindGenes <- read_csv("~/Documents/PhD_Project/Code/1st_project/WF_plots_postTAC/2024_03_12/output_ggCPP_FindGenes.txt")
View(output_ggCPP_FindGenes)
colnames(output_ggCPP_FindGenes) <- output_ggCPP_FindGenes[10,1]
output_ggCPP_FindGenes <- output_ggCPP_FindGenes[-(1:10),]
for (i in 1:49) {
  likelihoods_FindGenes_df[i,3] <- strsplit(output_ggCPP_FindGenes[9 * i,1]$`[1] "ggCaller_PopPUNK"`, "\\[1\\]")[[1]][2]
}
#sort.list(likelihoods_FindGenes_df$ggCmanSeqCl) # should be decreasing=TRUE but it does not seem to recognise the minus sign
sort.list(likelihoods_FindGenes_df$ggCPP)[1:10]
#  5 44 46 21 12 38 30  2 19 28

rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCPP)[1:10]])
sum(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCPP)[1:10]])==1) # 568
sum(rowSums(find_genes_df[,sort.list(likelihoods_FindGenes_df$ggCPP)[1:10]])==2) # 71
