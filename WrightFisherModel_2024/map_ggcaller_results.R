mmseq_results_MassNavajo <- read.delim("~/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/mmseq_results_MassNavajo.m8", header=FALSE)

nrow(mmseq_results_MassNavajo)

gene_Mass_Navajo_dict <- c()
curr_gene <- ""
curr_val <- 0
for (i in 1:nrow(mmseq_results_MassNavajo)) {
  if(curr_gene == mmseq_results_MassNavajo$V1[i]){
    if(curr_val < mmseq_results_MassNavajo$V3[i]){
      curr_val <- mmseq_results_MassNavajo$V3[i]
      gene_Mass_Navajo_dict[curr_gene] <- mmseq_results_MassNavajo$V2[i]
    }
  }
  else{
    curr_gene <- mmseq_results_MassNavajo$V1[i]
    curr_val <- mmseq_results_MassNavajo$V3[i]
    gene_Mass_Navajo_dict[curr_gene] <- mmseq_results_MassNavajo$V2[i]
  }
}

mmseq_results_NavajoMass <- read.delim("~/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/mmseq_results_NavajoMass.m8", header=FALSE)

gene_Navajo_Mass_dict <- c()
curr_gene <- ""
curr_val <- 0
for (i in 1:nrow(mmseq_results_NavajoMass)) {
  if(curr_gene == mmseq_results_NavajoMass$V1[i]){
    if(curr_val < mmseq_results_NavajoMass$V3[i]){
      curr_val <- mmseq_results_NavajoMass$V3[i]
      gene_Navajo_Mass_dict[curr_gene] <- mmseq_results_NavajoMass$V2[i]
    }
  }
  else{
    curr_gene <- mmseq_results_NavajoMass$V1[i]
    curr_val <- mmseq_results_NavajoMass$V3[i]
    gene_Navajo_Mass_dict[curr_gene] <- mmseq_results_NavajoMass$V2[i]
  }
}

### check whether Mass to Navajo mapping is the same as Navajo to Mass mapping
mismatch_example_navajo <- ""
mismatch_example_mass <- ""
mismatch_count <- 0
for (i in 1:length(gene_Navajo_Mass_dict)) {
  navajo <- names(gene_Navajo_Mass_dict)[i]
  mass <- gene_Navajo_Mass_dict[i]
  if(gene_Mass_Navajo_dict[mass] != navajo){
    mismatch_count <- mismatch_count +1
    mismatch_example_navajo <- navajo
    mismatch_example_mass <- mass
  }
}

mismatch_count2 <- 0
for (i in 1:length(gene_Mass_Navajo_dict)) {
  mass <- names(gene_Mass_Navajo_dict)[i]
  navajo <- gene_Mass_Navajo_dict[i]
  if(gene_Navajo_Mass_dict[navajo] != mass){
    mismatch_count2 <- mismatch_count2 +1
  }
}

# do greedy version
mmseq_results_MassNavajo_sorted <- mmseq_results_MassNavajo[order(mmseq_results_MassNavajo$V3, decreasing = TRUE),]

gene_Mass_Navajo_dict_greedy <- rep("", length(unique(mmseq_results_MassNavajo_sorted$V1)))
names(gene_Mass_Navajo_dict_greedy) <- unique(mmseq_results_MassNavajo_sorted$V1)
val_Mass_Navajo_dict_greedy <- rep(0, length(unique(mmseq_results_MassNavajo_sorted$V1)))
names(val_Mass_Navajo_dict_greedy) <- unique(mmseq_results_MassNavajo_sorted$V1)
removed_genes <- c()
curr_gene <- ""
#curr_val <- 0
for (i in 1:nrow(mmseq_results_MassNavajo_sorted)) {
  curr_gene <- mmseq_results_MassNavajo_sorted$V1[i]
  if(val_Mass_Navajo_dict_greedy[curr_gene] < mmseq_results_MassNavajo_sorted$V3[i]){
    if(! is.element(mmseq_results_MassNavajo_sorted$V2[i], removed_genes)){
      val_Mass_Navajo_dict_greedy[curr_gene] <- mmseq_results_MassNavajo_sorted$V3[i]
      gene_Mass_Navajo_dict_greedy[curr_gene] <- mmseq_results_MassNavajo_sorted$V2[i]
      removed_genes <- c(removed_genes, mmseq_results_MassNavajo_sorted$V2[i])
    }
  }  
}

# do greedy version Navajo Mass
mmseq_results_NavajoMass_sorted <- mmseq_results_NavajoMass[order(mmseq_results_NavajoMass$V3, decreasing = TRUE),]

gene_Navajo_Mass_dict_greedy <- rep("", length(unique(mmseq_results_NavajoMass_sorted$V1)))
names(gene_Navajo_Mass_dict_greedy) <- unique(mmseq_results_NavajoMass_sorted$V1)
val_Navajo_Mass_dict_greedy <- rep(0, length(unique(mmseq_results_NavajoMass_sorted$V1)))
names(val_Navajo_Mass_dict_greedy) <- unique(mmseq_results_NavajoMass_sorted$V1)
removed_genes <- c()
curr_gene <- ""
#curr_val <- 0
for (i in 1:nrow(mmseq_results_NavajoMass_sorted)) {
  curr_gene <- mmseq_results_NavajoMass_sorted$V1[i]
  if(gene_Navajo_Mass_dict_greedy[curr_gene] < mmseq_results_NavajoMass_sorted$V3[i]){
    if(! is.element(mmseq_results_NavajoMass_sorted$V2[i], removed_genes)){
      val_Mass_Navajo_dict_greedy[curr_gene] <- mmseq_results_NavajoMass_sorted$V3[i]
      gene_Navajo_Mass_dict_greedy[curr_gene] <- mmseq_results_NavajoMass_sorted$V2[i]
      removed_genes <- c(removed_genes, mmseq_results_NavajoMass_sorted$V2[i])
    }
  }  
}

### check whether greedy Mass to Navajo mapping is more similar
mismatch_example_navajo <- ""
mismatch_example_mass <- ""
mismatch_count <- 0
for (i in 1:length(gene_Navajo_Mass_dict_greedy)) {
  navajo <- names(gene_Navajo_Mass_dict_greedy)[i]
  mass <- gene_Navajo_Mass_dict_greedy[i]
  if(gene_Mass_Navajo_dict_greedy[mass] != navajo | ! is.element(mass, names(gene_Mass_Navajo_dict_greedy))){
    mismatch_count <- mismatch_count + 1
    mismatch_example_navajo <- navajo
    mismatch_example_mass <- mass
  }
}

# new attempt: maximize the 1-to-1 mapping by taking the maximum of the sum of map_Nav_Mass and map_Mass_Nav
length(unique(mmseq_results_MassNavajo$V1))
length(unique(mmseq_results_NavajoMass$V1))

Mass_ggCaller <- read.csv(paste(path_to_data, "Massachusetts_ggcaller/run_withFuncAnn/ggCaller_output/gene_presence_absence.csv", sep = ""), header=FALSE)
Navajo_ggCaller <- read.csv(paste(path_to_data, "StrepPneumo_NavajoNew/ggCaller_output/gene_presence_absence.csv", sep = ""), header=FALSE)

Mass_to_Ind_dict <- 1:(nrow(Mass_ggCaller)-1)
names(Mass_to_Ind_dict) <- Mass_ggCaller$V1[-1]

Nav_to_Ind_dict <- 1:(nrow(Navajo_ggCaller)-1)
names(Nav_to_Ind_dict) <-  Navajo_ggCaller$V1[-1]


map_MassNav_mtx <- matrix(NA, nrow = (nrow(Mass_ggCaller)-1), ncol = (nrow(Navajo_ggCaller)-1))
map_NavMass_mtx_tr <- matrix(NA, nrow = (nrow(Mass_ggCaller)-1), ncol = (nrow(Navajo_ggCaller)-1))

for (i in 1:nrow(mmseq_results_MassNavajo)) {
  map_MassNav_mtx[Mass_to_Ind_dict[mmseq_results_MassNavajo[i,1]], Nav_to_Ind_dict[mmseq_results_MassNavajo[i,2]]] <- mmseq_results_MassNavajo[i,3]
}

for(i in 1:nrow(mmseq_results_NavajoMass)){
  map_NavMass_mtx_tr[Mass_to_Ind_dict[mmseq_results_NavajoMass[i,2]], Nav_to_Ind_dict[mmseq_results_NavajoMass[i,1]]] <- mmseq_results_NavajoMass[i,3]
}

map_sum <- map_MassNav_mtx + map_NavMass_mtx_tr

whichmax_NA <- function(row){
  row_whichmax <- which.max(row)
  if(length(row_whichmax)==0){
    row_whichmax <- NA
  }
  row_whichmax
}

max_map_dict <- names(Nav_to_Ind_dict)[unlist(apply(map_sum, 1, whichmax_NA))]
names(max_map_dict) <- Mass_ggCaller$V1[-1]

max_map_dict2 <- names(Mass_to_Ind_dict)[apply(map_sum, 2, whichmax_NA)]
names(max_map_dict2) <- Navajo_ggCaller$V1[-1]

Mass_intermed_gene_presence_absence <- readRDS("ggC_intermed_gene_presence_absence.rds")
Mass_intermed_genes <- Mass_intermed_gene_presence_absence$V1[-1]

Navajo_ggCaller_intermed <- readRDS("Navajo_ggCaller_intermed.rds")
Nav_intermed_genes <- Navajo_ggCaller_intermed$Gene[-1]

length(max_map_dict[Mass_intermed_genes])
length(unique(max_map_dict[Mass_intermed_genes]))

hit <- 0
miss <- 0
miss_type1 <- 0
mass_val2 <- ""
for (mass_val in Mass_intermed_gene_presence_absence$V1[-1]) {
  nav_val <- max_map_dict[mass_val]
  mass_val2 <- max_map_dict2[nav_val]
  #print(mass_val)
  #print(mass_val2)
  if(!is.na(mass_val2) & mass_val2 == mass_val){
    hit <- hit + 1
  }
  else if(is.na(mass_val2)){
    miss_type1 <- miss_type1 + 1
  }
  else{
    miss <- miss +1
  }
}
hit
#[1] 1167
miss
#[1] 587
miss_type1
#[1] 16