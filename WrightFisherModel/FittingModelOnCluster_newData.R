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
#install.packages("coda", repos='http://cran.us.r-project.org')
library(odin.dust)
library(hrbrthemes)
library(ggplot2)
library(gridExtra)
library(readxl)
library(dplyr)
library(mcstate)
library(coda)

path_to_scripts <- "/nfs/research/jlees/leonie/WF_fitting/"
path_to_data <- "/nfs/research/jlees/leonie/WF_fitting/Data/"
path_to_results <- "/nfs/research/jlees/leonie/WF_fitting/Plots/"

#local_path <- "/Users/llorenz/Documents/PhD_Project/Code/1st_project/odin-dust-examples/TestModels/"
#local_path_to_script <- paste(local_path,"WrightFisher_nGenotypes_haploid_PopsizeVariablePois.R", sep="")
#WF_nG_h_vP <- odin.dust::odin_dust(local_path_to_script)

#read in Nick Croucher's data set
mass_data <- read.delim(paste(path_to_data,"mass.input", sep=""))

###
# compute gene frequencies
sum_as_int <- function(x){
  sum(as.integer(x))
}

gene_freq <- rep(0, ncol(mass_data)-5)
gene_freq_mass_data <- mass_data[1:nrow(mass_data),6:ncol(mass_data)]

gene_freq <- apply(gene_freq_mass_data,2, sum_as_int)
gene_freq <- gene_freq / (length(gene_freq_mass_data[,1]))

#create sub-data set for the intermediate frequency genes (5-95% across all sequences)
intermed_mass_data <- data.frame(matrix(ncol = sum(gene_freq <= 0.95 & gene_freq >= 0.05) + 5, nrow = nrow(mass_data)))
intermed_mass_data[,1:5] <- mass_data[,1:5]
colnames(intermed_mass_data)[1:5] <- colnames(mass_data)[1:5]
counter <- 1
for (i in 1:length(gene_freq)){
  if (0.05 <= gene_freq[i] & gene_freq[i] <= 0.95 ){
    intermed_mass_data[,counter+5] <- mass_data[,i+5]
    colnames(intermed_mass_data)[counter+5] <- colnames(mass_data)[i+5]
    counter <- counter + 1
  }
}

### 
# compute consensus genomes for sequence clusters
# and I will transform the data frame, so that the genes are the rows and the columns are the clusters
mass_clusters <- length(unique(mass_data$SC))
mass_consensus_presence_absence <- data.frame(matrix(nrow = sum(gene_freq <= 0.95 & gene_freq >= 0.05), ncol = mass_clusters+1))

mass_consensus_presence_absence[,1] <- colnames(intermed_mass_data)[-(1:5)]
colnames(mass_consensus_presence_absence) <- c('genes', sort(unique(mass_data$SC)))


mean_mass_cluster_freq <- mass_consensus_presence_absence # this is just to assess how valid it is to make consensus genomes

cons_genomes <- function(x){
  as.double(median(as.integer(x)))
}

cons_genomes_mean <- function(x){
  as.double(mean(as.integer(x)))
}

for (i in 0:(mass_clusters-1)){
  curr_genomes <- intermed_mass_data[which(intermed_mass_data$SC==i),-(1:5)] # select all genomes in cluster i
  #cluster_gene_presence_absence[j,i+3] <- as.double(median(as.integer(curr_genomes[j,])))
  mass_consensus_presence_absence[,i+2] <- apply(curr_genomes,2,cons_genomes)
  
  mean_mass_cluster_freq[,i+2] <- apply(curr_genomes,2,cons_genomes_mean)
}
###
#calculate the frequency of the gene clusters and year
mass_cluster_freq_1 <- rep(0,mass_clusters)
mass_cluster_freq_2 <- rep(0,mass_clusters)
mass_cluster_freq_3 <- rep(0,mass_clusters)
for (i in 1:mass_clusters){
  mass_cluster_freq_1[i] <- length(which(mass_data[which(mass_data$SC==i-1),]$Time==0))
  mass_cluster_freq_2[i] <- length(which(mass_data[which(mass_data$SC==i-1),]$Time==36))
  mass_cluster_freq_3[i] <- length(which(mass_data[which(mass_data$SC==i-1),]$Time==72))
}

###
# calculate Vaccine Type consensus for clusters
mass_VT <- rep(0, mass_clusters)
mass_VT_mean <- rep(0, mass_clusters)

for (i in 1:mass_clusters){
  curr_VTs <- subset(mass_data,mass_data$SC == i-1)$VT # select all genomes in cluster i-1
  mass_VT[i] <- ceiling(median(curr_VTs))
  mass_VT_mean[i] <- mean(curr_VTs)
  #print(mass_VT_mean[i]) #calculate mean as check for consensus
}

###
# calculate the beta statistics for including low and high levels of selection

# calculate gene frequencies first, separate for three times
mass_gene_freq_0 <- apply(intermed_mass_data[which(intermed_mass_data$Time==0),][,-(1:5)], 2, sum)
mass_gene_freq_36 <- apply(intermed_mass_data[which(intermed_mass_data$Time==36),][,-(1:5)], 2, sum)
mass_gene_freq_72 <- apply(intermed_mass_data[which(intermed_mass_data$Time==72),][,-(1:5)], 2, sum)

# first, calculate pre/peri and post vacc frequencies of genes:
pre_peri_vacc_gene_freq <- (mass_gene_freq_0 + mass_gene_freq_36) / (nrow(subset(mass_data,mass_data$Time==0)) + nrow(subset(mass_data,mass_data$Time==36)))
pre_vacc_gene_freq <- (mass_gene_freq_0) / nrow(subset(mass_data,mass_data$Time==0))
post_vacc_gene_freq <- mass_gene_freq_72 / nrow(subset(mass_data,mass_data$Time==72))
peri_post <- (mass_gene_freq_36 + mass_gene_freq_72) / (nrow(subset(mass_data,mass_data$Time==36)) + nrow(subset(mass_data,mass_data$Time==72)))

# calculate delta statistic (refer to Corander et al. for more info)
#delta_data <- (post_vacc_gene_freq - pre_peri_vacc_gene_freq) ^ 2 / (1 - pre_peri_vacc_gene_freq * (1 - pre_peri_vacc_gene_freq))
#delta_ranking <- rank(delta_data)

delta_data3 <- (peri_post - pre_vacc_gene_freq) ^ 2 / (1 - pre_vacc_gene_freq * (1 - pre_vacc_gene_freq))
delta_ranking3 <- rank(delta_data3)
delta_ranking <- delta_ranking3

###
### create initial population that is based on the 2001 data set but not an exact sampling from it
# but a Poisson process
expand_factor <- 15000 / sum(mass_cluster_freq_1)
exp_noise <- 10
model_start_pop <- (sapply((mass_cluster_freq_1 + rexp(n = length(mass_cluster_freq_1), rate = exp_noise)) * expand_factor, rpois, n=1))

### 
#loading the model from the odin file
path_to_rscript <- paste(path_to_scripts,"WrightFisher_newData_nGenotypes_haploid_PopsizeVariablePois.R", sep="")
WF_nG_h_vP <- odin.dust::odin_dust(path_to_rscript)

avg_cluster_freq <- (mass_cluster_freq_1 + mass_cluster_freq_2 + mass_cluster_freq_3)/(sum(mass_cluster_freq_1)+sum(mass_cluster_freq_2)+sum(mass_cluster_freq_3))

params_n_vP <- list(dt = 1/36, species_no = mass_clusters,  gene_no = nrow(mass_consensus_presence_absence), Pop_ini = as.double(model_start_pop), Pop_eq = as.double(model_start_pop), capacity = sum(model_start_pop), Genotypes = as.matrix(mass_consensus_presence_absence[,-1]), sigma_f = 0.14, sigma_w = 0.002, prop_f = 0.25, delta = delta_ranking, m = 0.03, migVec = avg_cluster_freq, vaccTypes = mass_VT, v = 0.1, vacc_time = 0)

WFmodel_nG_h_vP <- WF_nG_h_vP$new(pars = params_n_vP,
                                  time = 1,
                                  n_particles = 10L,
                                  n_threads = 4L,
                                  seed = 1L)

###
# Function for creating lollipop plots
make_lollipop_plot <- function(data1, data2, data3, model1, model2, model3){
  
  # Create data
  lollipop_data_2001_nVT <- data.frame(
    x=which(mass_VT==0),
    model2001=model1[which(mass_VT==0)] / sum(model1),
    data2001=as.numeric(data1[which(mass_VT==0)] / sum(data1))
  )
  lollipop_data_2001_VT <- data.frame(
    x=which(mass_VT==1),
    model2001=model1[which(mass_VT==1)] / sum(model1),
    data2001=as.numeric(data1[which(mass_VT==1)] / sum(data1))
  )
  lollipop_data_2004_nVT <- data.frame(
    x=which(mass_VT==0),
    model2004=model2[which(mass_VT==0)] / sum(model2),
    data2004=as.numeric(data2[which(mass_VT==0)] / sum(data2))
  )
  lollipop_data_2004_VT <- data.frame(
    x=which(mass_VT==1),
    model2004=model2[which(mass_VT==1)] / sum(model2),
    data2004=as.numeric(data2[which(mass_VT==1)] / sum(data2))
  )
  lollipop_data_2007_nVT <- data.frame(
    x=which(mass_VT==0),
    model2007=model3[which(mass_VT==0)] / sum(model3),
    data2007=as.numeric(data3[which(mass_VT==0)] / sum(data3))
  )
  lollipop_data_2007_VT <- data.frame(
    x=which(mass_VT==1),
    model2007=model3[which(mass_VT==1)] / sum(model3),
    data2007=as.numeric(data3[which(mass_VT==1)] / sum(data3))
  )
  # Change baseline
  lollipop_plot_2001_nVT <- ggplot(lollipop_data_2001_nVT) +
    geom_segment( aes(x=x, xend=x, y=model2001, yend=data2001), color="grey") +
    geom_point( aes(x=x, y=model2001, color="Model non-VT"), size=3 ) +
    geom_point( aes(x=x, y=data2001, color="Data non_VT"), size=3 ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = "none") +
    ggtitle("2001 Non-Vaccine Types") +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_2001_nVT$model2001),max(lollipop_data_2001_nVT$data2001),max(lollipop_data_2004_VT$model2004),max(lollipop_data_2004_VT$data2004),max(lollipop_data_2007_nVT$model2007),max(lollipop_data_2007_nVT$data2007)))
  lollipop_plot_2001_VT <- ggplot(lollipop_data_2001_VT) +
    geom_segment( aes(x=x, xend=x, y=model2001, yend=data2001), color="grey") +
    geom_point( aes(x=x, y=model2001, color="Model non-VT"), size=3 ) +
    geom_point( aes(x=x, y=data2001, color="Data non_VT"), size=3 ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = "none") +
    ggtitle("2001 Vaccine Types") +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_2001_nVT$model2001),max(lollipop_data_2001_nVT$data2001),max(lollipop_data_2004_VT$model2004),max(lollipop_data_2004_VT$data2004),max(lollipop_data_2007_nVT$model2007),max(lollipop_data_2007_nVT$data2007)))
  lollipop_plot_2004_nVT <- ggplot(lollipop_data_2004_nVT) +
    geom_segment( aes(x=x, xend=x, y=model2004, yend=data2004), color="grey") +
    geom_point( aes(x=x, y=model2004, color="Model"), size=3 ) +
    geom_point( aes(x=x, y=data2004, color="Data"), size=3 ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = "none") +
    ggtitle("2004 Non-Vaccine Types") +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_2001_nVT$model2001),max(lollipop_data_2001_nVT$data2001),max(lollipop_data_2004_VT$model2004),max(lollipop_data_2004_VT$data2004),max(lollipop_data_2007_nVT$model2007),max(lollipop_data_2007_nVT$data2007)))          
  lollipop_plot_2004_VT <- ggplot(lollipop_data_2004_VT) +
    geom_segment( aes(x=x, xend=x, y=model2004, yend=data2004), color="grey") +
    geom_point( aes(x=x, y=model2004, color="Model"), size=3 ) +
    geom_point( aes(x=x, y=data2004, color="Data"), size=3 ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = "none") +
    ggtitle("2004 Vaccine Types") +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_2001_nVT$model2001),max(lollipop_data_2001_nVT$data2001),max(lollipop_data_2004_VT$model2004),max(lollipop_data_2004_VT$data2004),max(lollipop_data_2007_nVT$model2007),max(lollipop_data_2007_nVT$data2007))) 
  lollipop_plot_2007_nVT <- ggplot(lollipop_data_2007_nVT) +
    geom_segment( aes(x=x, xend=x, y=model2007, yend=data2007), color="grey") +
    geom_point( aes(x=x, y=model2007, color="Model"), size=3 ) +
    geom_point( aes(x=x, y=data2007, color="Data"), size=3 ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       guide  = guide_legend(), 
                       name   = "Group") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = "none") +
    ggtitle("2007 Non-Vaccine Types") +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_2001_nVT$model2001),max(lollipop_data_2001_nVT$data2001),max(lollipop_data_2004_VT$model2004),max(lollipop_data_2004_VT$data2004),max(lollipop_data_2007_nVT$model2007),max(lollipop_data_2007_nVT$data2007)))          
  lollipop_plot_2007_VT <- ggplot(lollipop_data_2007_VT) +
    geom_segment( aes(x=x, xend=x, y=model2007, yend=data2007), color="grey") +
    geom_point( aes(x=x, y=model2007, color="Model"), size=3 ) +
    geom_point( aes(x=x, y=data2007, color="Data"), size=3 ) +
    scale_color_manual(values = c("#E69F00", "#56B4E9"),
                       guide  = guide_legend(), 
                       name   = "Legend") +
    coord_flip()+
    #theme_ipsum() +
    theme(legend.position = c(.8,.8),legend.text = element_text(size = 20),legend.title = element_text(size = 20)) +
    ggtitle("2007 Vaccine Types") +
    ylab("Frequency") +
    xlab("Clusters") +
    theme(axis.title  = element_text(size = 20), axis.text = element_text(size = 20), plot.title = element_text(size = 25,hjust = 0.5))  +
    ylim(0, max(max(lollipop_data_2001_nVT$model2001),max(lollipop_data_2001_nVT$data2001),max(lollipop_data_2004_VT$model2004),max(lollipop_data_2004_VT$data2004),max(lollipop_data_2007_nVT$model2007),max(lollipop_data_2007_nVT$data2007)))
  
  grid.arrange(lollipop_plot_2001_nVT + theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm")),lollipop_plot_2004_nVT+ theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm")), lollipop_plot_2007_nVT + theme(plot.margin = unit(c(.5,0.5,1,0.5), "cm")),lollipop_plot_2001_VT + theme(plot.margin = unit(c(.5,0.5,.5,0.5), "cm")), lollipop_plot_2004_VT + theme(plot.margin = unit(c(.5,0.5,.5,0.5), "cm")), lollipop_plot_2007_VT + theme(plot.margin = unit(c(.5,0.5,.5,0.5), "cm")), ncol = 3, nrow=2)
}

########
#FITTING

###
# process data with particle filter:
dt <- 1/36 # we assume that the generation time of Strep. pneumo is 1 month
# we have data from 2001, 2004, 2007, so we want 3 (years) * 12 (months) = 36 updates in-between

peripost_mass_cluster_freq <- data.frame("year" = c(1, 2), rbind(mass_cluster_freq_2, mass_cluster_freq_3))
names(peripost_mass_cluster_freq) <- c("year", as.character(1:mass_clusters))


# (until now I had 2004 and 2007 there and then in the model that resulted in times 72144 and 72252. does not make sense)
# when changing this back, also change initial_time back to 2001

fitting_mass_data <- mcstate::particle_filter_data(data = peripost_mass_cluster_freq,
                                                   time = "year",
                                                   rate = 1 / dt,
                                                   initial_time = 0)

### 
# Likelihood definition

# log-likelihood of with Poisson distribution
ll_pois <- function(obs, model) {
  exp_noise <- 1e6
  
  if (is.na(obs)) {
    # Creates vector of zeros in ll with same length, if no data
    ll_obs <- numeric(length(model))
  } else {
    lambda <- model + rexp(n = length(model), rate = exp_noise)
    ll_obs <- dpois(x = obs, lambda = lambda, log = TRUE)
    #print(model)
  }
  ll_obs
}

combined_compare <- function(state, observed, pars = NULL) {
  result <- 0
  #data_size <- sum(observed[as.character(1:62)])
  # cannot get this to work for now
  # instead I will use this:
  data_size <- sum(mass_cluster_freq_1)
  #model_size <- sum(state[, , drop = TRUE])
  model_size = 15000
  
  for (i in 1:mass_clusters){
    #browser()
    #print(observed[[as.character(i)]])
    #print(state[1+i, , drop = TRUE]/model_size * data_size)
    #print(state[1+i, , drop = TRUE])
    result <- result + ll_pois(observed[[as.character(i)]], state[1+i, , drop = TRUE]/model_size * data_size)
    #test_likelihood <- ll_pois(observed[[as.character(i)]], state[1+i, , drop = TRUE]/model_size * data_size)
    #print(result)
  }
  #print(result)
  result
}

###
#setting up deterministic particle filter
det_filter <- particle_deterministic$new(data = fitting_mass_data,
                                         model = WF_nG_h_vP,
                                         compare = combined_compare)
###
# Defining and transforming parameters
pmcmc_sigma_f <- mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0, max = 1)
pmcmc_sigma_w <- mcstate::pmcmc_parameter("sigma_w", 0.05, min = 0, max = 1)
pmcmc_prop_f <- mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1)
pmcmc_m <- mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1)
pmcmc_v <- mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)
species_no <- mass_clusters
no_clusters <- mass_clusters
gene_no <- nrow(as.matrix(mass_consensus_presence_absence[,-1]))
Pop_ini <- model_start_pop
Pop_eq <- model_start_pop
Genotypes <- as.matrix(mass_consensus_presence_absence[,-1])
capacity <- sum(model_start_pop)
delta <- delta_ranking
vaccTypes <- mass_VT
vacc_time <- 0
dt <- 1/36
migVec <- avg_cluster_freq

complex_params <- c(Pop_ini, Pop_eq, Genotypes, capacity, delta, vaccTypes, species_no, gene_no, vacc_time, dt, migVec)

make_transform <- function(p) {
  function(theta){
    c(list(Pop_ini = p[1:mass_clusters],
           Pop_eq = p[(mass_clusters +1) : (mass_clusters + mass_clusters)],
           Genotypes = matrix(p[(mass_clusters + mass_clusters + 1): ((mass_clusters + mass_clusters + 1) + (gene_no * species_no) - 1)], nrow = gene_no, ncol = species_no),
           capacity = p[((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 1],
           delta = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2) : (((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no -1)],
           vaccTypes = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) : ((((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters -1)],
           species_no = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters],
           gene_no = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters + 1],
           vacc_time = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters + 2],
           dt = p[(((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters + 3],
           migVec = p[((((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters + 4):((((2 * mass_clusters + 1) + (gene_no * species_no) - 1) + 2 + gene_no) + no_clusters + 4 + species_no - 1)]), as.list(theta))
  }
  
}

transform <- function(x) {
  make_transform(complex_params)}
proposal_matrix <- diag(0.1, 5) # the proposal matrix defines the covariance-variance matrix for a mult normal dist

mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.15, min = 0.075, max = 1), mcstate::pmcmc_parameter("sigma_w", 0.05, min = 0, max = 0.0749), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 1), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 1)), proposal_matrix, make_transform(complex_params))

mcmc_pars$initial()

### 
# DETERMINISTIC FITTING 1

#mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.1432, min = 0.075, max = 0.22), mcstate::pmcmc_parameter("sigma_w", 0.0011, min = 0, max = 0.0749), mcstate::pmcmc_parameter("prop_f", 0.25, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.03, min = 0, max = 0.2), mcstate::pmcmc_parameter("v", 0.05, min = 0, max = 0.5)), proposal_matrix, make_transform(complex_params))
#det_filter <- particle_deterministic$new(data = fitting_mass_data,
#                                         model = WF_nG_h_vP,
#                                         compare = combined_compare)
#n_steps <- 10000
#n_burnin <- 0

#control <- mcstate::pmcmc_control(
#  n_steps,
#  save_state = TRUE, 
#  save_trajectories = TRUE,
#  progress = TRUE,
#  adaptive_proposal = TRUE,
#  n_chains = 4)
#det_pmcmc_run <- mcstate::pmcmc(mcmc_pars, det_filter, control = control)
#par(mfrow = c(1,1))

#det_mcmc1 <- coda::as.mcmc(cbind(det_pmcmc_run$probabilities, det_pmcmc_run$pars))

# save mcmc diagnostic plots
# %03d guarantees that there are multiple files for the two diagnostic plots
#png(file=paste(path_to_results, "detRun1_diagnostics%03d.png", sep=""),width=1500, height=1300)
#plot(det_mcmc1)
#dev.off()

# remove burn-in, estimate param values
#processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run, burnin = 500, thin = 1)
#parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
#print("Parameters of Deterministic Fit 1")
#print(parameter_mean_hpd)
### I will need to check whether this printed stuff will go into the output file of the codon cluster or whether this will get lost somehow

#use lollipop plotting function
#png(file=paste(path_to_results, "detRun1_lollipop.png", sep=""),width=2000, height=1700)
#make_lollipop_plot(data1 = mass_cluster_freq_1, data2 = mass_cluster_freq_2, data3 = mass_cluster_freq_3, model1 = processed_chains$trajectories$state[-1,1,1], model2 = processed_chains$trajectories$state[-1,1,2], model3 = processed_chains$trajectories$state[-1,1,3])
#dev.off()

###
#DETERMINISTIC FIT 2
#det_proposal_matrix <- cov(processed_chains$pars)
#det_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", parameter_mean_hpd[1], min = 0.075, max = 0.22), mcstate::pmcmc_parameter("sigma_w", parameter_mean_hpd[2], min = 0.000001, max = 0.0749), mcstate::pmcmc_parameter("prop_f", parameter_mean_hpd[3], min = 0, max = 1), mcstate::pmcmc_parameter("m", parameter_mean_hpd[4], min = 0, max = 0.2), mcstate::pmcmc_parameter("v", parameter_mean_hpd[5], min = 0, max = 0.5)), det_proposal_matrix, make_transform(complex_params))

#n_steps <- 50000
#n_burnin <- 0


#control <- mcstate::pmcmc_control(
#  n_steps,
#  save_state = TRUE, 
#  save_trajectories = TRUE,
#  progress = TRUE,
#  adaptive_proposal = TRUE,
#  n_chains = 4)
#det_pmcmc_run2 <- mcstate::pmcmc(det_mcmc_pars, det_filter, control = control)
#det_mcmc2 <- coda::as.mcmc(cbind(det_pmcmc_run2$probabilities, det_pmcmc_run2$pars))
#plotting diagnostics and saving the plots
#png(file=paste(path_to_results, "detRun2_diagnostics%03d.png", sep=""),width=1500, height=1300)
#plot(det_mcmc2)
#dev.off()

# remove burn-in, estimate param values
#processed_chains <- mcstate::pmcmc_thin(det_pmcmc_run2, burnin = 1000, thin = 1)
#parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
#print("Parameters of Deterministic Fit 2")
#print(parameter_mean_hpd)

#use lollipop plotting function
#png(file=paste(path_to_results, "detRun2_lollipop.png", sep=""),width=2000, height=1700)
#make_lollipop_plot(data1 = mass_cluster_freq_1, data2 = mass_cluster_freq_2, data3 = mass_cluster_freq_3, model1 = processed_chains$trajectories$state[-1,1,1], model2 = processed_chains$trajectories$state[-1,1,2], model3 = processed_chains$trajectories$state[-1,1,3])
#dev.off()


#print squared error:
#print("Squared Error for Deterministic Fit 2:")
#print("Squared Error for Time Point 1:")
#print(sum((processed_chains$trajectories$state[-1,1,1]/sum(processed_chains$trajectories$state[-1,1,1]) - mass_cluster_freq_1/sum(mass_cluster_freq_1))^2))
#print("Squared Error for Time Point 2:")
#print(sum((processed_chains$trajectories$state[-1,1,2]/sum(processed_chains$trajectories$state[-1,1,2]) - mass_cluster_freq_2/sum(mass_cluster_freq_2))^2))
#print("Squared Error for Time Point 3:")
#print(sum((processed_chains$trajectories$state[-1,1,3]/sum(processed_chains$trajectories$state[-1,1,3]) - mass_cluster_freq_3/sum(mass_cluster_freq_3))^2))

###
# 1st STOCHASTIC FIT
#stoch_proposal_matrix <- cov(processed_chains$pars)
#stoch_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", parameter_mean_hpd[1], min = 0.075, max = 0.22), mcstate::pmcmc_parameter("sigma_w", parameter_mean_hpd[2], min = 0.000001, max = 0.0749), mcstate::pmcmc_parameter("prop_f", parameter_mean_hpd[3], min = 0, max = 1), mcstate::pmcmc_parameter("m", parameter_mean_hpd[4], min = 0, max = 0.2), mcstate::pmcmc_parameter("v", parameter_mean_hpd[5], min = 0, max = 0.5)), stoch_proposal_matrix, make_transform(complex_params))
stoch_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", 0.145663676, min = 0.075, max = 0.22), mcstate::pmcmc_parameter("sigma_w", 0.001221254, min = 0.000001, max = 0.0749), mcstate::pmcmc_parameter("prop_f", 0.058899571, min = 0, max = 1), mcstate::pmcmc_parameter("m", 0.149041718, min = 0, max = 0.2), mcstate::pmcmc_parameter("v", 0.087650161, min = 0, max = 0.5)), proposal_matrix, make_transform(complex_params))
n_particles <- 100
filter <- mcstate::particle_filter$new(data = fitting_mass_data,
                                       model = WF_nG_h_vP,
                                       n_particles = n_particles,
                                       compare = combined_compare,
                                       seed = 1L)


n_steps <- 1000
n_burnin <- 0
control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE, 
  save_trajectories = TRUE,
  progress = TRUE,
  n_chains = 3)
pmcmc_run <- mcstate::pmcmc(stoch_mcmc_pars, filter, control = control)

stoch_mcmc1 <- coda::as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars))

#plotting diagnostics and saving the plots
png(file=paste(path_to_results, "stochRun1_diagnostics%03d.png", sep=""),width=1500, height=1300)
plot(stoch_mcmc1)
dev.off()

processed_chains <- mcstate::pmcmc_thin(pmcmc_run, burnin = 100, thin = 1)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
print("Parameters of Stochastic Fit 1")
print(parameter_mean_hpd)

#use lollipop plotting function
png(file=paste(path_to_results, "stochRun1_lollipop.png", sep=""),width=2000, height=1700)
make_lollipop_plot(data1 = mass_cluster_freq_1, data2 = mass_cluster_freq_2, data3 = mass_cluster_freq_3, model1 = processed_chains$trajectories$state[-1,1,1], model2 = processed_chains$trajectories$state[-1,1,2], model3 = processed_chains$trajectories$state[-1,1,3])
dev.off()

###
# 2nd STOCHASTIC FIT
stoch_proposal_matrix <- cov(processed_chains$pars)
stoch_mcmc_pars <- mcstate::pmcmc_parameters$new(list(mcstate::pmcmc_parameter("sigma_f", parameter_mean_hpd[1], min = 0.075, max = 0.22), mcstate::pmcmc_parameter("sigma_w", parameter_mean_hpd[2], min = 0.000001, max = 0.0749), mcstate::pmcmc_parameter("prop_f", parameter_mean_hpd[3], min = 0, max = 1), mcstate::pmcmc_parameter("m", parameter_mean_hpd[4], min = 0, max = 0.2), mcstate::pmcmc_parameter("v", parameter_mean_hpd[5], min = 0, max = 0.5)), stoch_proposal_matrix, make_transform(complex_params))
n_particles <- 100
filter <- mcstate::particle_filter$new(data = fitting_mass_data,
                                       model = WF_nG_h_vP,
                                       n_particles = n_particles,
                                       compare = combined_compare,
                                       seed = 1L)


mcmc_pars <- stoch_mcmc_pars

n_steps <- 5000
n_burnin <- 0
control <- mcstate::pmcmc_control(
  n_steps,
  save_state = TRUE, 
  save_trajectories = TRUE,
  progress = TRUE,
  n_chains = 3)
pmcmc_run2 <- mcstate::pmcmc(stoch_mcmc_pars, filter, control = control)


stoch_mcmc2 <- coda::as.mcmc(cbind(pmcmc_run2$probabilities, pmcmc_run2$pars))


png(file=paste(path_to_results, "stochRun2_diagnostics%03d.png", sep=""),width=1500, height=1300)
plot(stoch_mcmc2)
dev.off()

processed_chains <- mcstate::pmcmc_thin(pmcmc_run2, burnin = 100, thin = 1)
parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
print("Parameters of Stochastic Fit 2")
print(parameter_mean_hpd)

#use lollipop plotting function
png(file=paste(path_to_results, "stochRun2_lollipop.png", sep=""),width=2000, height=1700)
make_lollipop_plot(data1 = mass_cluster_freq_1, data2 = mass_cluster_freq_2, data3 = mass_cluster_freq_3, model1 = processed_chains$trajectories$state[-1,1,1], model2 = processed_chains$trajectories$state[-1,1,2], model3 = processed_chains$trajectories$state[-1,1,3])
dev.off()


#print squared error:
print("Squared Error for Stochastic Fit 2:")
print("Squared Error for Time Point 1:")
print(sum((processed_chains$trajectories$state[-1,1,1]/sum(processed_chains$trajectories$state[-1,1,1]) - mass_cluster_freq_1/sum(mass_cluster_freq_1))^2))
print("Squared Error for Time Point 2:")
print(sum((processed_chains$trajectories$state[-1,1,2]/sum(processed_chains$trajectories$state[-1,1,2]) - mass_cluster_freq_2/sum(mass_cluster_freq_2))^2))
print("Squared Error for Time Point 3:")
print(sum((processed_chains$trajectories$state[-1,1,3]/sum(processed_chains$trajectories$state[-1,1,3]) - mass_cluster_freq_3/sum(mass_cluster_freq_3))^2))