library(odin.dust)
library(hrbrthemes)
library(gridExtra)
library(readxl)
library(dplyr)
library(ggplot2)

WF_nG_h_vP <- odin.dust::odin_dust("WrightFisher_nGenotypes_haploid_PopsizeVariablePois.R")

# reading in the cluster produced by PopPUNK
clusters <- read.csv("../Data/refined_modelfitk3_clusters.csv")
no_clusters <- max(clusters[,2]) # number of clusters in dataset

# reading in the gene presence absence matrix produced by ggCaller
gene_presence_absence <- read.csv("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/gene_presence_absence.csv", header=FALSE)

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
accNo_to_filename <- read.delim("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/filereport_read_run_PRJEB2632_tsv.txt")
accNo_to_filename <- accNo_to_filename[,c(1,8)]
accNo_to_filename[,3] <- matrix(unlist(strsplit(accNo_to_filename[,2],"/")), ncol=6, byrow = TRUE)[,6]
accNo_to_filename[,3] <- matrix(unlist(strsplit(accNo_to_filename[,3],"[.]")), ncol=2, byrow = TRUE)[,1]
accNo_to_filename <- accNo_to_filename[,c(1,3)]
colnames(accNo_to_filename) <- c(colnames(accNo_to_filename)[1], "filenames")
accNoToFilename <- c()
accNoToFilename[as.character(accNo_to_filename[,2])] <- accNo_to_filename[,1]


Croucher_seqYears <- read_excel("~/Documents/PhD_Project/Code/1st_project/odin-dust-examples/Data/Croucher_41588_2013_BFng2625_MOESM28_ESM.xlsx")
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
#data_gene_freq_1 <- rep(0, length(data_gene_freq_1_m[,1]))
data_gene_freq_1 <- apply(data_gene_freq_1_m[,1:62], 1, sum)

data_gene_freq_2_m <-  matrix_cluster_gene_presence_absence * cluster_freq_2
data_gene_freq_2 <- apply(data_gene_freq_2_m[,1:62], 1, sum)

data_gene_freq_3_m <-  matrix_cluster_gene_presence_absence * cluster_freq_3
data_gene_freq_3 <- apply(data_gene_freq_3_m[,1:62], 1, sum)

#par(mfrow=c(1,3))
#plot(data_gene_freq_1 / nrow(subset(clusters,clusters$seqYear==2001)), data_gene_freq_2 / nrow(subset(clusters,clusters$seqYear==2004)))
#plot(data_gene_freq_1 / nrow(subset(clusters,clusters$seqYear==2001)), data_gene_freq_3 / nrow(subset(clusters,clusters$seqYear==2007)))
#plot(data_gene_freq_2 / nrow(subset(clusters,clusters$seqYear==2004)), data_gene_freq_3 / nrow(subset(clusters,clusters$seqYear==2007)))
# very strong correlation between gene frequencies before and after vaccination
# especially second plot looks a lot like figure 2 d, the first, in Corander et al. (though they only did pre vs. post vacc)
#par(mfrow=c(1,1))
#plot((data_gene_freq_1 + data_gene_freq_2) / (nrow(subset(clusters,clusters$seqYear==2001)) + nrow(subset(clusters,clusters$seqYear==2004))), data_gene_freq_3 / nrow(subset(clusters,clusters$seqYear==2007)))

# correlation of clusters
#par(mfrow=c(1,3))
#plot(cluster_freq_1 / sum(cluster_freq_1), cluster_freq_2 / sum(cluster_freq_2))
#plot(cluster_freq_1 / sum(cluster_freq_1), cluster_freq_3 / sum(cluster_freq_3))
#plot(cluster_freq_2 / sum(cluster_freq_2), cluster_freq_3 / sum(cluster_freq_3))
# the second plot (2001 vs. 2007) looks a lot like figure 2 c, the first, in Corander et al.

#par(mfrow=c(1,1))
#plot((cluster_freq_1 + cluster_freq_2)/ (sum(cluster_freq_1) + sum(cluster_freq_2)), cluster_freq_3 / sum(cluster_freq_3))

# read in information on vaccine types
#library(readxl)
vaccine_types <- read_excel("../Data/Corander_suppData3.xlsx")
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

# make a grouped bar plot similar to Croucher et al.
#library(ggplot2)

#cluster_names <- rep(1:no_clusters, 3)
#cluster_names <- sort(cluster_names)
#cluster_names <- as.character(cluster_names)
#seqTimes <- rep(c("2001", "2004", "2007"), no_clusters)
#cluster_freq_all_times <- rep(0, no_clusters *3)
#for (i in 1:no_clusters) {
#  cluster_freq_all_times[3*(i-1) + 1 ] <- cluster_freq_1[i]/sum(cluster_freq_1)
#  cluster_freq_all_times[3 * (i-1) + 2] <- cluster_freq_2[i]/sum(cluster_freq_2)
#  cluster_freq_all_times[3 * (i-1) + 3] <- cluster_freq_3[i]/sum(cluster_freq_3)
#}
#cluster_freq_df <- data.frame(cluster_names, seqTimes, cluster_freq_all_times)

#ggplot(cluster_freq_df, aes(fill=seqTimes, y=cluster_freq_all_times, x=cluster_names)) +   geom_bar(position="dodge", stat="identity") +
#  scale_fill_manual("legend", values = c("2001" = "#E69F00", "2004" = "#56B4E9", "2007" = "#009E73")) +
#  scale_x_continuous(breaks = 1:no_clusters)

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
#start_time <- Sys.time()
start_pop <- as.vector(as.double(rmultinom(1,50000,prob = cluster_freq_1/sum(cluster_freq_1))))
#end_time <- Sys.time()
#end_time - start_time
# Time difference of 0.001049042 secs

#start_time <- Sys.time()
#sample_vec <- sum(cluster_freq_1)
#ind_sum <- 1
#for (i in 1:no_clusters) {
#  if(cluster_freq_1[i] != 0){
#    sample_vec[ind_sum : (ind_sum + cluster_freq_1[i] - 1)] <- rep(i, cluster_freq_1[i])
#    ind_sum <- ind_sum + cluster_freq_1[i]
#  }
#}
#sampled_vec <- sample(sample_vec, 50000, replace = TRUE)
#start_pop <- rep(0, no_clusters)
#j <- 1
#counter <- 0
#for (i in sampled_vec) {
#  start_pop[i] <- start_pop[i] + 1
#}
#end_time <- Sys.time()
#end_time - start_time
# Time difference of 0.03029203 secs
#Yes, I could make the sampling faster by using apply but I don't think that it can be faster than the multinomial
# so I want to stick to the multinomial sampling 
# because the sampling frequencies are virtually the same
# (see plots below)
#matplot(1:no_clusters, start_pop, type = "p", pch = 1, col = "#56B4E9")
#matpoints(1:no_clusters, start_pop_old, type = "p", pch = 2, col = "#E69F00")

#plot frequencies instead of absolute numbers
#matplot(1:no_clusters, start_pop/sum(start_pop), type = "p", pch = 1, col = "#56B4E9")
#matpoints(1:no_clusters, start_pop_old/sum(start_pop_old), type = "p", pch = 2, col = "#E69F00")
#matpoints(1:no_clusters, cluster_freq_1/sum(cluster_freq_1), type = "p", pch = 3, col = "black")


### Define model parameters according to the datasets
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

sw_sigma_f <- c(0.1, 0.2, 0.3)
sw_sigma_w <- c(0.001, 0.005, 0.01)
sw_prop_f <- c(0.1, 0.2, 0.3)
sw_m <- c(0.01, 0.05, 0.1)
sw_v <- c(0.01, 0.05, 0.1)

sw_sigma_f <- c(0.1500)
sw_sigma_w <- c(0.0075)
sw_prop_f <- c(0.1000)
sw_m <- c(0.0250)
sw_v <- c(0.0050)

#corander params
sw_sigma_f <- c(0.1363)
sw_sigma_w <- c(0.0023)
sw_prop_f <- c(0.2483)
sw_m <- c(0.0044)
sw_v <- c(0.0812)

counter <- 1

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
          
          plot_title <- paste(paste(paste(paste(paste("sigma_f", as.character(swsf)), paste(", sigma_w", as.character(swsw))), paste(", prop_f", as.character(swpf))), paste(", m", as.character(swm))), paste(", v", as.character(swv)))
          
          # Create data
          lollipop_data_2001 <- data.frame(
            x=1:no_clusters, 
            model2001=x_mean[,1] / sum(x_mean[,1]),
            data2001=as.numeric(cluster_freq_1 / sum(cluster_freq_1))
          )
          lollipop_data_2004 <- data.frame(
            x=1:no_clusters, 
            model2004=x_mean[,2] / sum(x_mean[,2]),
            data2004=as.numeric(big_population_data[1,2:63]/sum(big_population_data[1,2:63]))
          )
          lollipop_data_2007 <- data.frame(
            x=1:no_clusters, 
            model2007=x_mean[,3] / sum(x_mean[,3]),
            data2007=as.numeric(big_population_data[2,2:63]/sum(big_population_data[2,2:63]))
          )
          # Change baseline
          lollipop_plot_2001 <- ggplot(lollipop_data_2001) +
            geom_segment( aes(x=x, xend=x, y=model2001, yend=data2001), color="grey") +
            geom_point( aes(x=x, y=model2001, color="Model"), size=3 ) +
            geom_point( aes(x=x, y=data2001, color="Data"), size=3 ) +
            scale_color_manual(values = c("#E69F00", "#56B4E9"),
                               guide  = guide_legend(), 
                               name   = "Group") +
            coord_flip()+
            theme_ipsum() +
            theme(legend.position = "none") +
            ggtitle(plot_title) +
            xlab("") +
            ylab("Value of Y") +
            ylim(0, max(max(lollipop_data_2001$model2001),max(lollipop_data_2001$data2001),max(lollipop_data_2004$model2004),max(lollipop_data_2004$data2004),max(lollipop_data_2007$model2007),max(lollipop_data_2007$data2007)))
          lollipop_plot_2004 <- ggplot(lollipop_data_2004) +
            geom_segment( aes(x=x, xend=x, y=model2004, yend=data2004), color="grey") +
            geom_point( aes(x=x, y=model2004, color="Model"), size=3 ) +
            geom_point( aes(x=x, y=data2004, color="Data"), size=3 ) +
            scale_color_manual(values = c("#E69F00", "#56B4E9"),
                               guide  = guide_legend(), 
                               name   = "Group") +
            coord_flip()+
            theme_ipsum() +
            theme(legend.position = "none") +
            xlab("") +
            ylab("Value of Y") +
            ylim(0, max(max(lollipop_data_2001$model2001),max(lollipop_data_2001$data2001),max(lollipop_data_2004$model2004),max(lollipop_data_2004$data2004),max(lollipop_data_2007$model2007),max(lollipop_data_2007$data2007)))          
          lollipop_plot_2007 <- ggplot(lollipop_data_2007) +
            geom_segment( aes(x=x, xend=x, y=model2007, yend=data2007), color="grey") +
            geom_point( aes(x=x, y=model2007, color="Model"), size=3 ) +
            geom_point( aes(x=x, y=data2007, color="Data"), size=3 ) +
            scale_color_manual(values = c("#E69F00", "#56B4E9"),
                               guide  = guide_legend(), 
                               name   = "Group") +
            coord_flip()+
            theme_ipsum() +
            theme(legend.position = c(.8,.8)) +
            xlab("") +
            ylab("Value of Y") + 
            ylim(0, max(max(lollipop_data_2001$model2001),max(lollipop_data_2001$data2001),max(lollipop_data_2004$model2004),max(lollipop_data_2004$data2004),max(lollipop_data_2007$model2007),max(lollipop_data_2007$data2007)))          
        
          #save plot as a png
          g <- arrangeGrob(lollipop_plot_2001,lollipop_plot_2004, lollipop_plot_2007, ncol =3) #generates g
          ggsave(file=paste(paste("plot",as.character(counter), sep=""), ".png", sep=""), g, width = 16, height = 16) #saves g
          counter <- counter + 1
        }
      }
    }
    
  }
}

# calculate likelihood instead of plots
big_population_data <- peripost_fitting_cluster_freq_df
big_population_data[,2:63] <- big_population_data[,2:63] * 208

big_mass_data <- mcstate::particle_filter_data(data = big_population_data,
                                               time = "year",
                                               rate = 1 / dt,
                                               initial_time = 2001)

big_filter <- mcstate::particle_filter$new(data = big_mass_data,
                                           model = WF_nG_h_vP,
                                           n_particles = n_particles,
                                           compare = combined_compare,
                                           seed = 1L)


counter <- 1

likelihood_vec <- rep(0, length(sw_sigma_f) * length(sw_sigma_w) * length(sw_prop_f) * length(sw_m) * length(sw_v))

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
          
          
          counter <- counter + 1
        }
      }
    }
  }
}     
likelihood_vec
max(likelihood_vec)
which.max(likelihood_vec)

which(likelihood_vec > -51000)
# answer is 232 

param_values <- matrix(rep(0, 5* length(sw_sigma_f) * length(sw_sigma_w) * length(sw_prop_f) * length(sw_m) * length(sw_v)), nrow = length(sw_sigma_f) * length(sw_sigma_w) * length(sw_prop_f) * length(sw_m) * length(sw_v), ncol = 5) 
counter <- 1
for (swsf in sw_sigma_f) {
  for (swsw in sw_sigma_w) {
    for (swpf in sw_prop_f) {
      for (swm in sw_m) {
        for (swv in sw_v) {
            param_values[counter,] <- c(swsf, swsw, swpf, swm, swv)
            counter <- counter + 1
        }
      }
    }
  }
}

param_comb_best_fits <- param_values[which(likelihood_vec > -51000),]
param_comb_best_fits
# [,1]  [,2] [,3] [,4] [,5]
#[1,]  0.1 0.005  0.2 0.05 0.01
#[2,]  0.2 0.001  0.2 0.05 0.01
#[3,]  0.2 0.005  0.2 0.05 0.01
#[4,]  0.2 0.010  0.2 0.05 0.01
#[5,]  0.3 0.001  0.2 0.10 0.01
#[6,]  0.3 0.005  0.3 0.10 0.01
#[7,]  0.3 0.010  0.2 0.10 0.01
#[8,]  0.3 0.010  0.3 0.10 0.01
# okay, interesting: vaccination influence is always 0.01 in these (lowest tested setting)
# migration is better when higher (0.05, 0.1) (not surprising, can compensate issues in model)
# prop f is better when higher (0.2, 0.3)
# sigma_w is better when higher (0.005, 0.01)
# sigma_f is better when higher (0.2, 0.3)




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
which.max(likelihood_vec)

which(likelihood_vec > -50000 & likelihood_vec < 0)
param_comb_best_fits <- param_values[which(likelihood_vec > -50000 & likelihood_vec < 0),]
param_comb_best_fits
### best parameter values:
# (I did interrupt the run before it was done though)
#      [,1]   [,2] [,3]  [,4]  [,5]
#[1,] 0.10 0.0075 0.10 0.025 0.005
#[2,] 0.10 0.0100 0.10 0.025 0.005
#[3,] 0.10 0.0150 0.20 0.025 0.005
#[4,] 0.10 0.0150 0.20 0.050 0.005
#[5,] 0.10 0.0150 0.25 0.050 0.005
#[6,] 0.10 0.0150 0.30 0.050 0.005
#[7,] 0.10 0.0175 0.15 0.050 0.005
#[8,] 0.10 0.0175 0.15 0.050 0.010
#[9,] 0.10 0.0175 0.25 0.050 0.005
#[10,] 0.15 0.0050 0.10 0.025 0.005
#[11,] 0.15 0.0075 0.10 0.025 0.005
#[12,] 0.15 0.0075 0.20 0.050 0.005
#[13,] 0.15 0.0075 0.25 0.050 0.005
#[14,] 0.15 0.0100 0.10 0.025 0.005
#[15,] 0.15 0.0125 0.20 0.050 0.005
#[16,] 0.15 0.0150 0.15 0.050 0.005
#[17,] 0.15 0.0150 0.20 0.050 0.005
#[18,] 0.15 0.0200 0.25 0.050 0.005
#[19,] 0.20 0.0010 0.20 0.050 0.005
#[20,] 0.20 0.0050 0.10 0.025 0.005

#sigma f has only been tested for 0.1 and 0.15 so far I think


sw_sigma_f <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
sw_sigma_w <- c(0.001, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
sw_prop_f <- c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35)
sw_m <- c(0, 0.025, 0.05, 0.75, 0.1, 0.125, 0.15)
sw_v <- c(0.005, 0.01, 0.25, 0.05, 0.1)


#param_values[which.max(Filter(smaller_0, likelihood_vec)),]
#[1] 0.1500 0.0075 0.1000 0.0250 0.0050