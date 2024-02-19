# reading in the cluster produced by PopPUNK
database_k3_clusters <- read.csv("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/refined_modelfitk3_clusters.csv")
no_clusters <- max(database_k3_clusters[,2]) # number of clusters in dataset

# reading in the gene presence absence matrix produced by ggCaller
gene_presence_absence <- read.csv("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/gene_presence_absence.csv", header=FALSE)

# converting the gene presence absence matrix into a boolean df (0 = gene not present, 1 = gene present)
convert_to_bool <- function(x){
  if (x=="") 0 else 1
}

bool_gene_presence_absence <- gene_presence_absence

for (i in 4:length(bool_gene_presence_absence[1,])) {

  bool_gene_presence_absence[,i] <- unlist(lapply(bool_gene_presence_absence[,i], FUN =  convert_to_bool))
}
bool_gene_presence_absence[1,] <- gene_presence_absence[1,]
  
# calculate frequency of genes to only keep those which appear in 5-95% of the genomes
gene_freq <- rep(0, nrow(bool_gene_presence_absence)-1)
for (i in 1:length(gene_freq)) {
  gene_freq[i] <- sum(as.integer(bool_gene_presence_absence[i+1,4:ncol(bool_gene_presence_absence)]))
}
gene_freq <- gene_freq / (length(bool_gene_presence_absence[2,])-3)

# create a dataframe that only contains the genes that appear in 5-95% of the genomes
filtered_bool_gene_presence_absence <- data.frame(matrix(nrow = 0, ncol = length(bool_gene_presence_absence[1,])))

for (i in 1:length(gene_freq)){
  if (0.05 < gene_freq[i] & gene_freq[i] < 0.95 ){
    filtered_bool_gene_presence_absence[nrow(filtered_bool_gene_presence_absence)+1,] <- bool_gene_presence_absence[i+1,]
  }
}
colnames(filtered_bool_gene_presence_absence) <- bool_gene_presence_absence[1,]

### make consensus genome for clusters
# attempt 1: always let majority decide (if more than 50% in the cluster don't have gene then 0, else 1)
# should be easy to do with median
# would be nice to have the genome names as column names though. yupp.

cluster_gene_presence_absence <- data.frame(matrix(nrow = nrow(filtered_bool_gene_presence_absence), ncol = no_clusters+3))
colnames(cluster_gene_presence_absence)[1:3] <- colnames(filtered_bool_gene_presence_absence)[1:3]
cluster_gene_presence_absence[1:3] <- filtered_bool_gene_presence_absence[1:3]
colnames(cluster_gene_presence_absence)[4:ncol(cluster_gene_presence_absence)] <- 1:no_clusters
for (i in 1:no_clusters){
  curr_cluster <- subset(database_k3_clusters, database_k3_clusters$Cluster == i)[,1] # select all genomes in cluster i
  curr_genomes <- as.matrix(filtered_bool_gene_presence_absence[,curr_cluster])
  for (j in 1:nrow(filtered_bool_gene_presence_absence)){
    cluster_gene_presence_absence[j,i+3] <- as.integer(median(as.integer(curr_genomes[j,])))
  }
}

# attempt 2: record relative frequencies. Probably not advantageous.
# attempt 3: keep track of all existing variants. This can be done later


#need the information of when probes where taken
accNo_to_filename <- read.delim("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/filereport_read_run_PRJEB2632_tsv.txt")
accNo_to_filename <- accNo_to_filename[,c(1,8)]
accNo_to_filename[3] <- unlist(strsplit(accNo_to_filename[,2],"/"))[6]
accNo_to_filename[3] <- unlist(strsplit(accNo_to_filename[,3],"[.]"))[1]
accNo_to_filename <- accNo_to_filename[,c(1,3)]
colnames(accNo_to_filename) <- c(colnames(accNo_to_filename)[1], "filenames")

Croucher_seqYears <- read_excel("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/Croucher_41588_2013_BFng2625_MOESM28_ESM.xlsx")
Croucher_seqYears <- Croucher_seqYears[,c(1,5)]

# read in information on vaccine types
library(readxl)
vaccine_types <- read_excel("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/Corander_suppData3.xlsx")
vaccine_types <- vaccine_types[,c(2,3,4,6)]
vaccine_types_mass <- subset(vaccine_types, vaccine_types$Population == "Massachusetts")
isVT <- c()
for (i in 1:nrow(vaccine_types_mass)){
  isVT[vaccine_types_mass$`Accession Code`[i]] <- (vaccine_types_mass$`Vaccine Type`[i] == "VT")
}
