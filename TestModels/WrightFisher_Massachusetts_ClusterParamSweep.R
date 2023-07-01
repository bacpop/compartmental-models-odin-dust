# use the following command on the cluster to load r installation
#module load r/4.2.2
# then R to start R

#install.packages("drat", repos='http://cran.us.r-project.org') # -- if you don't have drat installed
#drat:::add("ncov-ic")
#install.packages("odin")
#install.packages("hrbrthemes", repos='http://cran.us.r-project.org')
#install.packages("gridExtra", repos='http://cran.us.r-project.org')
#install.packages("readxl", repos='http://cran.us.r-project.org')
#install.packages("dplyr", repos='http://cran.us.r-project.org')
#install.packages(c("hrbrthemes","gridExtra","readxl","dplyr"), repos='http://cran.us.r-project.org')
library(odin.dust)
library(hrbrthemes)
library(gridExtra)
library(readxl)
library(dplyr)
library(mcstate)

path_to_scripts <- "/nfs/research/jlees/leonie/WF_parameterSweep/"
path_to_data <- "/nfs/research/jlees/leonie/WF_parameterSweep/Data/"

local_path <- "/Users/llorenz/Documents/PhD_Project/Code/1st_project/odin-dust-examples/TestModels/"
#local_path_to_script <- paste(local_path,"WrightFisher_nGenotypes_haploid_PopsizeVariablePois.R", sep="")
#WF_nG_h_vP <- odin.dust::odin_dust(local_path_to_script)

path_to_rscript <- paste(path_to_scripts,"WrightFisher_nGenotypes_haploid_PopsizeVariablePois.R", sep="")
WF_nG_h_vP <- odin.dust::odin_dust(path_to_rscript)

# reading in the cluster produced by PopPUNK
clusters <- read.csv(paste(path_to_data,"refined_modelfitk3_clusters.csv", sep=""))
no_clusters <- max(clusters[,2]) # number of clusters in dataset

# reading in the gene presence absence matrix produced by ggCaller
gene_presence_absence <- read.csv(paste(path_to_data,"gene_presence_absence.csv", sep=""), header=FALSE)

# converting the gene presence absence matrix into a boolean df (0 = gene not present, 1 = gene present)
convert_to_bool <- function(x){
  if (x=="") 0 else 1
}
bool_gene_presence_absence <- gene_presence_absence
bool_gene_presence_absence[2:nrow(bool_gene_presence_absence),4:ncol(bool_gene_presence_absence)] <- apply(gene_presence_absence[2:nrow(gene_presence_absence),4:ncol(gene_presence_absence)],c(1,2), convert_to_bool)

# calculate frequency of genes to only keep those which appear in 5-95% of the genomes
#gene_freq <- rep(0, nrow(bool_gene_presence_absence)-1)
#for (i in 1:length(gene_freq)) {
#  gene_freq[i] <- sum(as.integer(bool_gene_presence_absence[i+1,4:ncol(bool_gene_presence_absence)]))
#}
# gene_freq <- gene_freq / (length(bool_gene_presence_absence[2,])-3)

sum_as_int <- function(x){
  sum(as.integer(x))
}

gene_freq <- rep(0, nrow(bool_gene_presence_absence)-1)
gf_bool_gene_presence_absence <- bool_gene_presence_absence[2:nrow(bool_gene_presence_absence),4:ncol(bool_gene_presence_absence)]

gene_freq <- apply(gf_bool_gene_presence_absence,1, sum_as_int)
gene_freq <- as.vector(gene_freq / (length(bool_gene_presence_absence[2,])-3))

# create a dataframe that only contains the genes that appear in 5-95% of the genomes
filtered_bool_gene_presence_absence <- data.frame(matrix(nrow = sum(gene_freq <= 0.95 & gene_freq >= 0.05), ncol = length(bool_gene_presence_absence[1,])))
counter <- 1
for (i in 1:length(gene_freq)){
  if (0.05 <= gene_freq[i] & gene_freq[i] <= 0.95 ){
    filtered_bool_gene_presence_absence[counter,] <- bool_gene_presence_absence[i+1,]
    counter <- counter + 1
  }
}
colnames(filtered_bool_gene_presence_absence) <- bool_gene_presence_absence[1,]
rownames(filtered_bool_gene_presence_absence) <- 1:nrow(filtered_bool_gene_presence_absence)

### make consensus genome for clusters
# attempt 1: always let majority decide (if more than 50% in the cluster don't have gene then 0, else 1)
# should be easy to do with median

cluster_gene_presence_absence <- data.frame(matrix(nrow = nrow(filtered_bool_gene_presence_absence), ncol = no_clusters+3))
colnames(cluster_gene_presence_absence)[1:3] <- colnames(filtered_bool_gene_presence_absence)[1:3]
cluster_gene_presence_absence[1:3] <- filtered_bool_gene_presence_absence[1:3]
colnames(cluster_gene_presence_absence)[4:ncol(cluster_gene_presence_absence)] <- 1:no_clusters

cons_genomes <- function(x){
  as.double(median(as.integer(x)))
}

for (i in 1:no_clusters){
  curr_cluster <- clusters[which(clusters[,"Cluster"]==i),1] # select all genomes in cluster i
  curr_genomes <- as.matrix(filtered_bool_gene_presence_absence[,curr_cluster])
  cluster_gene_presence_absence[,i+3] <- apply(curr_genomes,1,cons_genomes)
}
matrix_cluster_gene_presence_absence <- as.matrix((cluster_gene_presence_absence[,4:ncol(cluster_gene_presence_absence)]))

# well, I'd say that many, many genes are still quite diversely present within the same cluster
# that suggest that taking one consensus genome per cluster might not work

# attempt 2: record relative frequencies. Probably not advantageous.
# attempt 3: keep track of all existing variants. This can be done later

#calculate the frequency of the gene clusters, well actually absolute numbers atm
cluster_freq <- rep(0,no_clusters)
for (i in 1:no_clusters){
  cluster_freq[i] <- length(which(clusters[,"Cluster"]==i))
}

#need the information of when probes were taken
accNo_to_filename <- read.delim(paste(path_to_data,"filereport_read_run_PRJEB2632_tsv.txt", sep=""))
accNo_to_filename <- accNo_to_filename[,c(1,8)]
accNo_to_filename[,3] <- matrix(unlist(strsplit(accNo_to_filename[,2],"/")), ncol=6, byrow = TRUE)[,6]
accNo_to_filename[,3] <- matrix(unlist(strsplit(accNo_to_filename[,3],"[.]")), ncol=2, byrow = TRUE)[,1]
accNo_to_filename <- accNo_to_filename[,c(1,3)]
colnames(accNo_to_filename) <- c(colnames(accNo_to_filename)[1], "filenames")
accNoToFilename <- c()
accNoToFilename[as.character(accNo_to_filename[,2])] <- accNo_to_filename[,1]


Croucher_seqYears <- read_excel(paste(path_to_data,"Croucher_41588_2013_BFng2625_MOESM28_ESM.xlsx", sep=""))
Croucher_seqYears <- Croucher_seqYears[,c(1,5)]
sequenceYear <- c()
for (i in 1:nrow(Croucher_seqYears)){
  sequenceYear[as.character(Croucher_seqYears[i,1])]<- Croucher_seqYears[i,2] 
}

# add year to the clusters data set
clusters$seqYear <- rep(0, nrow(clusters))
for (i in 1:nrow(clusters)){
  clusters$seqYear[i] <- sequenceYear[accNoToFilename[clusters$Taxon[i]]]
}

#calculate the frequency of the gene clusters and year
cluster_freq_1 <- rep(0,no_clusters)
cluster_freq_2 <- rep(0,no_clusters)
cluster_freq_3 <- rep(0,no_clusters)
for (i in 1:no_clusters){
  cluster_freq_1[i] <- length(which(clusters[which(clusters[,"Cluster"]==i),]$seqYear == 2001))
  cluster_freq_2[i] <- length(which(clusters[which(clusters[,"Cluster"]==i),]$seqYear == 2004))
  cluster_freq_3[i] <- length(which(clusters[which(clusters[,"Cluster"]==i),]$seqYear == 2007))
}

# make a data frame that contains this information

fitting_cluster_freq_df <- data.frame("year" = c(2001, 2004, 2007), rbind(cluster_freq_1, cluster_freq_2, cluster_freq_3))
names(fitting_cluster_freq_df) <- c("year", as.character(1:62))

#plot gene freqs
data_gene_freq_1_m <-  matrix_cluster_gene_presence_absence * cluster_freq_1
data_gene_freq_1 <- apply(data_gene_freq_1_m[,1:62], 1, sum)

data_gene_freq_2_m <-  matrix_cluster_gene_presence_absence * cluster_freq_2
data_gene_freq_2 <- apply(data_gene_freq_2_m[,1:62], 1, sum)

data_gene_freq_3_m <-  matrix_cluster_gene_presence_absence * cluster_freq_3
data_gene_freq_3 <- apply(data_gene_freq_3_m[,1:62], 1, sum)



vaccine_types <- read_excel(paste(path_to_data,"Corander_suppData3.xlsx", sep=""))
vaccine_types <- vaccine_types[,c(2,3,4,6)]
vaccine_types_mass <- vaccine_types[which(vaccine_types$Population == "Massachusetts"),]
isVT <- c()
for (i in 1:nrow(vaccine_types_mass)){
  isVT[vaccine_types_mass$`Accession Code`[i]] <- if (vaccine_types_mass$`Vaccine Type`[i] == "VT") {0} else {1}
}
# I want to create a simple vector instead of a dictionary because this will be easier to use when modelling
# and I also want a consensus for the clusters...
isVT_vec <- rep(0, no_clusters)
for (i in 1:no_clusters){
  curr_cluster <- subset(clusters, clusters$Cluster == i)[,1] # select all genomes in cluster i
  curr_vacc <- accNoToFilename[curr_cluster]
  isVT_vec[i] <- median(isVT[curr_vacc])
}


# I want to add low and high levels of selection to the model.
# for that I need to calculate the beta statistics
# first, calculate pre/peri and post vacc frequencies of genes:
pre_peri_vacc_gene_freq <- (data_gene_freq_1 + data_gene_freq_2) / (nrow(subset(clusters,clusters$seqYear==2001)) + nrow(subset(clusters,clusters$seqYear==2004)))
post_vacc_gene_freq <- data_gene_freq_3 / nrow(subset(clusters,clusters$seqYear==2007))

# calculate delta statistic (refer to Corander et al. for more info)
delta_data <- (post_vacc_gene_freq - pre_peri_vacc_gene_freq) ^ 2 / (1 - pre_peri_vacc_gene_freq * (1 - pre_peri_vacc_gene_freq))
delta_ranking <- rank(delta_data)

# implement sampling from pre vaccination population instead of directly using the 2001 population as input
# would it be better to do the sampling in the modelling frame work? maybe.
start_pop <- as.vector(as.double(rmultinom(1,50000,prob = cluster_freq_1/sum(cluster_freq_1))))

### Define model parameters according to the datasets
dt <- 1/36
capacity <- sum(start_pop)

params_n_vP <- list(dt = 1/36, 
                    species_no = no_clusters,  
                    gene_no = nrow(cluster_gene_presence_absence), 
                    Pop_ini = start_pop, 
                    Pop_eq = start_pop, 
                    capacity = sum(start_pop), 
                    Genotypes = matrix_cluster_gene_presence_absence, 
                    sigma_f = 0.14, 
                    sigma_w = 0.002, 
                    prop_f = 0.25, 
                    delta = delta_ranking, 
                    m = 0.03, 
                    vaccTypes = isVT_vec, 
                    v = 0.07,
                    vacc_time = 100) 
### Running the model:
#WFmodel_nG_h_vP <- WF_nG_h_vP$new(pars = params_n_vP,
#                                  time = 1,
#                                  n_particles = 10L,
#                                  n_threads = 4L,
#                                  seed = 1L)
peripost_fitting_cluster_freq_df <- fitting_cluster_freq_df[2:3,]
big_population_data <- peripost_fitting_cluster_freq_df
big_population_data[,2:63] <- big_population_data[,2:63] * 208

big_mass_data <- mcstate::particle_filter_data(data = big_population_data,
                                               time = "year",
                                               rate = 1 / dt,
                                               initial_time = 2001)


ll_pois <- function(obs, model) {
  exp_noise <- 1e6
  if (is.na(obs)) {
    # Creates vector of zeros in ll with same length, if no data
    ll_obs <- numeric(length(model))
  } else {
    lambda <- model + rexp(n = length(model), rate = exp_noise)
    #ll_obs <- dpois(x = obs, lambda = lambda, log = TRUE) # this is the simplest prior possible
    # it would be good to think about a binomial, neg binomial or beta binomial version
    # ll_obs <- dbinom() 
    # for neg binom: maybe n is the capacity of the system? the model value is the mu, and the y is the data?
    ll_obs <- dnbinom(x = obs, size = capacity, mu = lambda, log = TRUE)
  }
  ll_obs
}

combined_compare <- function(state, observed, pars = NULL) {
  result <- 0
  for (i in 1:62){
    #print(observed[as.character(ind)])
    #print(state[1+i, , drop = TRUE])
    result <- result + ll_pois(observed[[as.character(i)]], state[1+i, , drop = TRUE])
  }
  result
}

n_particles <- 100


big_filter <- mcstate::particle_filter$new(data = big_mass_data,
                                           model = WF_nG_h_vP,
                                           n_particles = n_particles,
                                           compare = combined_compare,
                                           seed = 1L)



#sw_sigma_f <- c(0.1)
#sw_sigma_w <- c(0.001)
#sw_prop_f <- c(0.1)
#sw_m <- c(0.01)
#sw_v <- c(0.01)


sw_sigma_f <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
sw_sigma_w <- c(0.001, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
sw_prop_f <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35)
sw_m <- c(0, 0.025, 0.05, 0.75, 0.1, 0.125, 0.15)
sw_v <- c(0.005, 0.01, 0.25, 0.05, 0.1)


counter <- 1

likelihood_vec <- rep(0, length(sw_sigma_f) * length(sw_sigma_w) * length(sw_prop_f) * length(sw_m) * length(sw_v))
param_values <- matrix(rep(0, 5* length(sw_sigma_f) * length(sw_sigma_w) * length(sw_prop_f) * length(sw_m) * length(sw_v)), nrow = length(sw_sigma_f) * length(sw_sigma_w) * length(sw_prop_f) * length(sw_m) * length(sw_v), ncol = 5) 


for (swsf in sw_sigma_f) {
  for (swsw in sw_sigma_w) {
    for (swpf in sw_prop_f) {
      for (swm in sw_m) {
        for (swv in sw_v) {
          params_n_vP <- list(dt = 1/36, 
                              species_no = no_clusters,  
                              gene_no = nrow(cluster_gene_presence_absence), 
                              Pop_ini = start_pop, 
                              Pop_eq = start_pop, 
                              capacity = sum(start_pop), 
                              Genotypes = matrix_cluster_gene_presence_absence, 
                              delta = delta_ranking,
                              vaccTypes = isVT_vec, 
                              sigma_f = swsf, 
                              sigma_w = swsw, 
                              prop_f = swpf, 
                              m = swm, 
                              v = swv,
                              vacc_time = 100)
          
          dt <- 1/36
          n_particles <- 10L
          WFmodel_nG_h_vP <- WF_nG_h_vP$new(pars = params_n_vP,
                                            time = 1,
                                            n_particles = 10L,
                                            n_threads = 4L,
                                            seed = 1L)
          n_times <- 72
          x <- array(NA, dim = c(WFmodel_nG_h_vP$info()$len, n_particles, n_times))
          
          for (t in seq_len(n_times)) {
            x[ , , t] <- WFmodel_nG_h_vP$run(t)
          }
          x <- x[,,c(1,36,72)]
          time <- x[1, 1, ]
          x <- x[-1, , ]
          x_mean <- apply(x, c(1,3), mean)
          
          
          likelihood_vec[counter] <- big_filter$run(save_history = TRUE, pars = params_n_vP)
          param_values[counter,] <- c(swsf, swsw, swpf, swm, swv)
          
          counter <- counter + 1
        }
      }
    }
  }
}     
#which.max(likelihood_vec)

#which(likelihood_vec > -51000)
param_comb_best_fits <- param_values[which(likelihood_vec > -51000),]
print("These are the best parameter combinations (better than -51000):\\")
print(param_comb_best_fits)
print("this is the likelihood of the best combination:\\")
print(max(likelihood_vec))
print("This is the best parameter combination:\\")
print(param_values[which.max(likelihood_vec),])