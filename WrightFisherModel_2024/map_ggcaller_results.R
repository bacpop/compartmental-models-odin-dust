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


### try mapping annotation instead
mass_gene_presence_absence <- read.csv("~/Documents/PhD_Project/Data/Massachusetts_ggcaller/run_withFuncAnn/ggCaller_output/gene_presence_absence.csv")
navajo_gene_presence_absence <- read.csv("~/Documents/PhD_Project/Data/StrepPneumo_NavajoNew/ggCaller_output/gene_presence_absence.csv")

mass_gene_anno_dict <- mass_gene_presence_absence$Gene
names(mass_gene_anno_dict) <- mass_gene_presence_absence$Annotation

navajo_anno_gene_dict <- navajo_gene_presence_absence$Annotation
names(navajo_anno_gene_dict) <- navajo_gene_presence_absence$Gene

mass_gene_anno_dict <- c()
for (i in 1:nrow(mass_gene_presence_absence)) {
  anno <- strsplit(mass_gene_presence_absence[i,3],"SPARC1\\_")[[1]][length(strsplit(mass_gene_presence_absence[i,3],"SPARC1\\_")[[1]])]
  if(! is.na(anno)){
    mass_gene_anno_dict[anno] <- mass_gene_presence_absence[i,1]
  }
}
# nrow(mass_gene_presence_absence) 5256
# length(mass_gene_anno_dict) 3549

mass_anno_gene_dict <- c()
for (i in 1:nrow(mass_gene_presence_absence)) {
  anno <- strsplit(mass_gene_presence_absence[i,3],"SPARC1\\_")[[1]][length(strsplit(mass_gene_presence_absence[i,3],"SPARC1\\_")[[1]])]
  if(! is.na(anno)){
    mass_anno_gene_dict[mass_gene_presence_absence[i,1]] <- anno 
  }
}
# nrow(mass_gene_presence_absence) 5256
# length(mass_anno_gene_dict) 4938
# 0.9394977 have annotation
# there are 3549 unique annotations

navajo_anno_gene_dict <- c()
for(i in 1:nrow(navajo_gene_presence_absence)){
  anno <- strsplit(navajo_gene_presence_absence[i,3],"SPARC1\\_")[[1]][length(strsplit(navajo_gene_presence_absence[i,3],"SPARC1\\_")[[1]])]
  if(! is.na(anno)){
    navajo_anno_gene_dict[mass_gene_presence_absence[i,1]] <- anno
  }
}
# nrow(navajo_gene_presence_absence) 5712
# length(navajo_anno_gene_dict) 5024
# 0.8795518 have annotation
# there are 3511 unique annotations

Navajo_matching_Mass <- data.frame(matrix(NA, nrow = length(navajo_anno_gene_dict), ncol = 2))
colnames(Navajo_matching_Mass) <- c("Navajo", "MatchingMass")
Navajo_matching_Mass$Navajo <- unname(navajo_anno_gene_dict)
missing_match_count <- 0
for (i in 1:nrow(Navajo_matching_Mass)) {
  matching_mass <- mass_gene_anno_dict[navajo_anno_gene_dict[i]]
  Navajo_matching_Mass$MatchingMass[i] <- matching_mass
  if(is.na(matching_mass)){
    missing_match_count <- missing_match_count +1
  }
}
# did not find a match for 222 sequences.
# 4802/5712 = 0.8406863 have mapped annotation to Mass

# plan:
# create dict with names=genes and values=anno for each dataset
# translate ggCaller presence-absence matrix into anno-space
# take translated matrix and translate it into other dataset
# problem: annotations that only exist in one dataset. Please think about this.

# 05.12.2024
mmseq_results_FindNavajoInMass <- read.delim("~/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/run3/bestResultNavajoInMass.m8", header=FALSE)

mmseq_results_FindMassInNavajo <- read.delim("~/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/run3_FindMassInNavajo/bestResultsMassInNavajo.m8", header=FALSE)

NavajoInMass_dict <- mmseq_results_FindNavajoInMass$V2
names(NavajoInMass_dict) <- mmseq_results_FindNavajoInMass$V1

MassInNavajo_dict <- mmseq_results_FindMassInNavajo$V2
names(MassInNavajo_dict) <- mmseq_results_FindMassInNavajo$V1

match_count <- 0
no_match_count <- 0
for (i in 1:length(names(NavajoInMass_dict))){
  name <- names(NavajoInMass_dict)[i]
  val <- NavajoInMass_dict[name]
  return_val <- MassInNavajo_dict[val]
  if(name != return_val){
    if(mmseq_results_FindNavajoInMass$V3[i]>0){
      no_match_count <- no_match_count + 1 
    }
  }
  else{
    if(mmseq_results_FindNavajoInMass$V3[i]>0){
    match_count <- match_count + 1
    }
  }
}
# no quality check
#> no_match_count (0.0)
#[1] 1343
#> match_count (0.0)
#[1] 3872

#> no_match_count (0.90)
#[1] 463
#> match_count (0.90)
#[1] 3103

#> no_match_count (0.95)
#[1] 204
#> match_count (0.95)
#[1] 1382

#> no_match_count (0.99)
#[1] 64
#> match_count (0.99)
#[1] 252

# this is already much more promising.
# I did the bestResult filtering / search with mmseqs2, which makes the mapping back and forth already much more consistent
# I need to think about how many matches I really expect, considering the singleton genes in both datasets

# 06.12.2024
# after a short conversation with Sam: filter representatives first (by 5-95% intermed freq genes), then do mapping with mmseqs and then compare forward/backward search
# things that do no match at all (prob not under NFDS)
# try 98% amino acids sequence identity
# rep is longest member btw

library(stringr)

Mass_gene_cluster_rep_str <- paste(readLines("~/Documents/PhD_Project/Data/Mapping_ggCaller/Massachusetts/pangenome_reference.fa"), collapse="\n")
Mass_gene_cluster_rep_str_split <- strsplit(Mass_gene_cluster_rep_str, ">")
Mass_gene_cluster_rep_dict <- c()
for (i in 2:length(Mass_gene_cluster_rep_str_split[[1]])) {
  local_split <- str_split_fixed(Mass_gene_cluster_rep_str_split[[1]][i],"\n",2)
  #cluster_name <- paste(">",local_split[1,1], sep = "")
  cluster_name <- local_split[1,1]
  cluster_seq <- local_split[1,2]
  Mass_gene_cluster_rep_dict[cluster_name] <- cluster_seq
}


Mass_ggC_intermed_gene_names <- readRDS(file = "Mass_ggC_intermed_gene_names.rds")
print_str <- ""
for (gene_name in Mass_ggC_intermed_gene_names) {
  print_str <- paste(print_str, ">", gene_name, "\n",Mass_gene_cluster_rep_dict[gene_name], sep = "")
}
fileConn<-file("~/Documents/PhD_Project/Data/Mapping_ggCaller/Massachusetts/pangenome_reference_filtered.fa")
writeLines(print_str, fileConn)
close(fileConn)

# repeat this for UK and Navajo
# and then run mmseqs again

filter_cluster_rep_file <- function(cluster_rep_file, intermed_gene_names, filtered_file){
  # read in fasta file with cluster representatives as string
  gene_cluster_rep_str_split <- strsplit(cluster_rep_file, ">")
  gene_cluster_rep_dict <- c()
  for (i in 2:length(gene_cluster_rep_str_split[[1]])) {
    local_split <- str_split_fixed(gene_cluster_rep_str_split[[1]][i],"\n",2)
    #cluster_name <- paste(">",local_split[1,1], sep = "")
    cluster_name <- local_split[1,1]
    cluster_seq <- local_split[1,2]
    gene_cluster_rep_dict[cluster_name] <- cluster_seq
  }
  # add intermediate gene names and sequences to file string
  print_str <- ""
  for (gene_name in intermed_gene_names) {
    print_str <- paste(print_str, ">", gene_name, "\n",gene_cluster_rep_dict[gene_name], sep = "")
  }
  # write output file
  writeLines(print_str, filtered_file)
  close(filtered_file)
}

### UK
UK_ggC_intermed_gene_names <- readRDS("UK_ggC_intermed_gene_names.rds")
UK_gene_cluster_rep_str <- paste(readLines("~/Documents/PhD_Project/Data/Mapping_ggCaller/UK/pangenome_reference.fa"), collapse="\n")
UK_filtered_reps <- file("~/Documents/PhD_Project/Data/Mapping_ggCaller/UK/pangenome_reference_filtered.fa")
filter_cluster_rep_file(UK_gene_cluster_rep_str, UK_ggC_intermed_gene_names, UK_filtered_reps)

### Navajo
Navajo_gene_cluster_rep_str <- paste(readLines("~/Documents/PhD_Project/Data/Mapping_ggCaller/Navajo/pangenome_reference.fa"), collapse="\n")
Navajo_ggC_intermed_gene_names <- readRDS("Navajo_ggC_intermed_gene_names.rds")
Navajo_filtered_reps <- file("~/Documents/PhD_Project/Data/Mapping_ggCaller/Navajo/pangenome_reference_filtered.fa")
filter_cluster_rep_file(Navajo_gene_cluster_rep_str, Navajo_ggC_intermed_gene_names, Navajo_filtered_reps)

# 11.12.2024
mmseq_results_FindUKInMass <- read.delim("~/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/FindIntermedUKinIntermedMass/bestResultUKinMass.m8", header=FALSE)

mmseq_results_FindMassInUK <- read.delim("~/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/FindIntermedMassInIntermedUK/bestResultMassInUK.m8", header=FALSE)

#UKInMass_dict_val <- 

UKInMass_dict <- mmseq_results_FindUKInMass$V2
names(UKInMass_dict) <- mmseq_results_FindUKInMass$V1

MassInUK_dict <- mmseq_results_FindMassInUK$V2
names(MassInUK_dict) <- mmseq_results_FindMassInUK$V1

UKMass_seq_id_vec <- mmseq_results_FindUKInMass$V3
MassUK_seq_id_vec <- mmseq_results_FindMassInUK$V3

recip_matching <- function(AinB_dict, BinA_dict, seq_id_vec, seq_identity = 0.95){
  match_count <- 0
  no_match_count <- 0
  match_dict <- c()
  
  for (i in 1:length(names(AinB_dict))){
    name <- names(AinB_dict)[i]
    val <- AinB_dict[name]
    return_val <- BinA_dict[val]
    if(is.na(return_val) | name != return_val){
      if(seq_id_vec[i]>=seq_identity){
        no_match_count <- no_match_count + 1 
      }
    }
    else{
      if(seq_id_vec[i]>=seq_identity){
        match_count <- match_count + 1
        match_dict[name] <- val
      }
    }
  }
  print(paste("recip matches", match_count))
  print(paste("non-recip matches", no_match_count))
  match_dict
}

#> match_count
#[1] 397
#> no_match_count
#[1] 78

Mass_ggC_intermed_gene_freqs <- readRDS("Mass_ggC_intermed_gene_freqs.rds")
names(Mass_ggC_intermed_gene_freqs) <- Mass_ggC_intermed_gene_names

UK_ggC_intermed_gene_freqs <- readRDS("UK_ggC_intermed_gene_freqs.rds")
names(UK_ggC_intermed_gene_freqs) <- UK_ggC_intermed_gene_names

match_UKMass_95 <- recip_matching(UKInMass_dict, MassInUK_dict, UKMass_seq_id_vec, 0.95)
match_UKMass_90 <- recip_matching(UKInMass_dict, MassInUK_dict, UKMass_seq_id_vec, 0.90)
match_UKMass_99 <- recip_matching(UKInMass_dict, MassInUK_dict, UKMass_seq_id_vec, 0.99)

plot(UK_ggC_intermed_gene_freqs[names(match_UKMass_90)], Mass_ggC_intermed_gene_freqs[match_UKMass_90[names(match_UKMass_90)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="Intermed. Gene Frequencies, 90% sequence identity")
abline(0,1)
plot(UK_ggC_intermed_gene_freqs[names(match_UKMass_95)], Mass_ggC_intermed_gene_freqs[match_UKMass_95[names(match_UKMass_95)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="Intermed. Gene Frequencies, 95% sequence identity")
abline(0,1)
plot(UK_ggC_intermed_gene_freqs[names(match_UKMass_99)], Mass_ggC_intermed_gene_freqs[match_UKMass_99[names(match_UKMass_99)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="Intermed. Gene Frequencies, 99% sequence identity")
abline(0,1)

# same with unfiltered matches
mmseq_results_FindUKInMass_unfiltered <- read.delim("~/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/FindAllUKinAllMass/bestResultUKinMass.m8", header=FALSE)
mmseq_results_FindMassInUK_unfiltered <- read.delim("~/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/FindAllMassInAllUK/bestResultMassInUK.m8", header=FALSE)

UKInMassUnfiltered_dict <- mmseq_results_FindUKInMass_unfiltered$V2
names(UKInMassUnfiltered_dict) <- mmseq_results_FindUKInMass_unfiltered$V1

MassInUKUnfiltered_dict <- mmseq_results_FindMassInUK_unfiltered$V2
names(MassInUKUnfiltered_dict) <- mmseq_results_FindMassInUK_unfiltered$V1

UKMassUnfiltered_seq_id_vec <- mmseq_results_FindUKInMass_unfiltered$V3
MassUKUnfiltered_seq_id_vec <- mmseq_results_FindMassInUK_unfiltered$V3

# make a plot that shows how frequent genes are with recip matches and those that do not have a recip match
Mass_ggC_all_gene_freqs_dict <- readRDS("Mass_ggC_all_gene_freqs.rds")
names(Mass_ggC_all_gene_freqs_dict) <- readRDS("Mass_ggC_all_gene_names.rds")

UK_ggC_all_gene_freqs_dict <- readRDS("UK_ggC_all_gene_freqs.rds")
names(UK_ggC_all_gene_freqs_dict) <- readRDS("UK_ggC_all_gene_names.rds")

match_UKMassUnfiltered_90 <- recip_matching(UKInMassUnfiltered_dict, MassInUKUnfiltered_dict, UKMassUnfiltered_seq_id_vec, 0.90)
match_UKMassUnfiltered_95 <- recip_matching(UKInMassUnfiltered_dict, MassInUKUnfiltered_dict, UKMassUnfiltered_seq_id_vec, 0.95)
match_UKMassUnfiltered_99 <- recip_matching(UKInMassUnfiltered_dict, MassInUKUnfiltered_dict, UKMassUnfiltered_seq_id_vec, 0.99)
match_MassUKUnfiltered_90 <- recip_matching(MassInUKUnfiltered_dict, UKInMassUnfiltered_dict, MassUKUnfiltered_seq_id_vec, 0.90)
match_MassUKUnfiltered_95 <- recip_matching(MassInUKUnfiltered_dict, UKInMassUnfiltered_dict, MassUKUnfiltered_seq_id_vec, 0.95)
match_MassUKUnfiltered_99 <- recip_matching(MassInUKUnfiltered_dict, UKInMassUnfiltered_dict, MassUKUnfiltered_seq_id_vec, 0.99)
par(pty="s")
plot(UK_ggC_all_gene_freqs_dict[names(match_UKMassUnfiltered_90)], Mass_ggC_all_gene_freqs_dict[match_UKMassUnfiltered_90[names(match_UKMassUnfiltered_90)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="All Gene Frequencies, 90% sequence identity")
abline(0,1)
plot(UK_ggC_all_gene_freqs_dict[names(match_UKMassUnfiltered_95)], Mass_ggC_all_gene_freqs_dict[match_UKMassUnfiltered_95[names(match_UKMassUnfiltered_95)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="All Gene Frequencies, 95% sequence identity")
abline(0,1)
plot(UK_ggC_all_gene_freqs_dict[names(match_UKMassUnfiltered_99)], Mass_ggC_all_gene_freqs_dict[match_UKMassUnfiltered_99[names(match_UKMassUnfiltered_99)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="All Gene Frequencies, 99% sequence identity")
abline(0,1)

### 13.12.2024
# colour these plots in by whether genes are predicted to be under NFDS or not
Mass_delta_ranking <- readRDS(file = "ggC_delta_ranking.rds")
UK_delta_ranking <- readRDS(file = "UK_delta_ranking.rds")
# UK genes under NFDS 4-param model (Dec 2024): 0.1432
# Mass genes under NFDS 4-param model (Dec 2024): 0.2797
# gene under NFDS if delta[i] <= prop_f * gene_no
Mass_delta_ranking_names <- Mass_delta_ranking
names(Mass_delta_ranking_names) <- Mass_ggC_intermed_gene_names
Mass_intermed_gene_underNFDS <- rep(0, length(Mass_delta_ranking_names))
names(Mass_intermed_gene_underNFDS) <- names(Mass_delta_ranking_names)
# new GPSC fit suggests Mass 0.2784
Mass_intermed_gene_underNFDS[which((Mass_delta_ranking_names <= 0.2784 * length(Mass_delta_ranking_names)))] <- 1
#if(Mass_delta_ranking <= 0.2797 * length(Mass_delta_ranking)) 1 else 0
UK_delta_ranking_names <- UK_delta_ranking
names(UK_delta_ranking_names) <- UK_ggC_intermed_gene_names
UK_intermed_gene_underNFDS <- rep(0, length(UK_delta_ranking_names))
names(UK_intermed_gene_underNFDS) <- names(UK_delta_ranking_names)
# and UK GPSC fit: 0.2095
UK_intermed_gene_underNFDS[which((UK_delta_ranking_names <= 0.2095 * length(UK_delta_ranking_names)))] <- 1
# expand this vector to all genes, NFDS or not
# Mass
Mass_underNFDS <- rep(0, length(Mass_ggC_all_gene_freqs_dict))
names(Mass_underNFDS) <- readRDS("Mass_ggC_all_gene_names.rds")
Mass_underNFDS[names(which(Mass_intermed_gene_underNFDS==1))] <- 1
# UK
UK_underNFDS <- rep(0, length(UK_ggC_all_gene_freqs_dict))
names(UK_underNFDS) <- readRDS("UK_ggC_all_gene_names.rds")
UK_underNFDS[names(which(UK_intermed_gene_underNFDS==1))] <- 1
# create color vector
col_clb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #8 colorblind friendly colors

colours_UKMassUnfiltered_90 <- rep("grey", length(match_UKMassUnfiltered_90))
names(colours_UKMassUnfiltered_90) <- names(match_UKMassUnfiltered_90)
colours_UKMassUnfiltered_90[names(which(UK_underNFDS[names(match_UKMassUnfiltered_90)]==1))] <- col_clb[2]
colours_UKMassUnfiltered_90[match_MassUKUnfiltered_90[names(which(Mass_underNFDS[match_UKMassUnfiltered_90]==1))]] <- col_clb[3]
colours_UKMassUnfiltered_90[intersect(names(which(UK_underNFDS[names(match_UKMassUnfiltered_90)]==1)),match_MassUKUnfiltered_90[names(which(Mass_underNFDS[match_UKMassUnfiltered_90]==1))])] <- col_clb[8]

colours_UKMassUnfiltered_95 <- rep("grey", length(match_UKMassUnfiltered_95))
names(colours_UKMassUnfiltered_95) <- names(match_UKMassUnfiltered_95)
colours_UKMassUnfiltered_95[names(which(UK_underNFDS[names(match_UKMassUnfiltered_95)]==1))] <- col_clb[2]
colours_UKMassUnfiltered_95[match_MassUKUnfiltered_95[names(which(Mass_underNFDS[match_UKMassUnfiltered_95]==1))]] <- col_clb[3]
colours_UKMassUnfiltered_95[intersect(names(which(UK_underNFDS[names(match_UKMassUnfiltered_95)]==1)),match_MassUKUnfiltered_95[names(which(Mass_underNFDS[match_UKMassUnfiltered_95]==1))])] <- col_clb[8]

colours_UKMassUnfiltered_99 <- rep("grey", length(match_UKMassUnfiltered_99))
names(colours_UKMassUnfiltered_99) <- names(match_UKMassUnfiltered_99)
colours_UKMassUnfiltered_99[names(which(UK_underNFDS[names(match_UKMassUnfiltered_99)]==1))] <- col_clb[2]
colours_UKMassUnfiltered_99[match_MassUKUnfiltered_99[names(which(Mass_underNFDS[match_UKMassUnfiltered_99]==1))]] <- col_clb[3]
colours_UKMassUnfiltered_99[intersect(names(which(UK_underNFDS[names(match_UKMassUnfiltered_99)]==1)),match_MassUKUnfiltered_99[names(which(Mass_underNFDS[match_UKMassUnfiltered_99]==1))])] <- col_clb[8]

par(pty="s")
plot(UK_ggC_all_gene_freqs_dict[names(match_UKMassUnfiltered_90)], Mass_ggC_all_gene_freqs_dict[match_UKMassUnfiltered_90[names(match_UKMassUnfiltered_90)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="All Gene Frequencies, 90% sequence identity", col = colours_UKMassUnfiltered_90, pch = 19)
abline(0,1)
plot(UK_ggC_all_gene_freqs_dict[names(match_UKMassUnfiltered_95)], Mass_ggC_all_gene_freqs_dict[match_UKMassUnfiltered_95[names(match_UKMassUnfiltered_95)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="All Gene Frequencies, 95% sequence identity", col = colours_UKMassUnfiltered_95, pch = 19)
abline(0,1)
plot(UK_ggC_all_gene_freqs_dict[names(match_UKMassUnfiltered_99)], Mass_ggC_all_gene_freqs_dict[match_UKMassUnfiltered_99[names(match_UKMassUnfiltered_99)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="All Gene Frequencies, 99% sequence identity", col = colours_UKMassUnfiltered_99, pch = 19)
abline(0,1)

# 14.01.2025
# map ggCaller results to COGs
mmseq_results_FindMassInMassCOGs_unfiltered <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/FindAllMassInMassCOGs/bestResultMassInMassCOGs.m8", header=FALSE)
mmseq_results_FindMassCOGsInMass_unfiltered <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/FindAllMassCOGsInMass/bestResultMassCOGsInMass.m8", header=FALSE)

MassInMassCOGs_unfiltered_dict <- mmseq_results_FindMassInMassCOGs_unfiltered$V2
names(MassInMassCOGs_unfiltered_dict) <- mmseq_results_FindMassInMassCOGs_unfiltered$V1

FindMassCOGsInMass_unfiltered_dict <- mmseq_results_FindMassCOGsInMass_unfiltered$V2
names(FindMassCOGsInMass_unfiltered_dict) <- mmseq_results_FindMassCOGsInMass_unfiltered$V1

MassInMassCOGs_unfiltered_seq_id_vec <- mmseq_results_FindMassInMassCOGs_unfiltered$V3
MassCOGsInMass_unfiltered_seq_id_vec <- mmseq_results_FindMassCOGsInMass_unfiltered$V3

# make a plot that shows how frequent genes are with recip matches and those that do not have a recip match
Mass_ggC_all_gene_freqs_dict <- readRDS("Mass_ggC_all_gene_freqs.rds")
names(Mass_ggC_all_gene_freqs_dict) <- readRDS("Mass_ggC_all_gene_names.rds")

Mass_cog_all_gene_freqs_dict <- readRDS("Mass_cog_all_gene_freqs.rds")
names(Mass_cog_all_gene_freqs_dict) <- readRDS("Mass_cog_all_gene_names.rds")

match_MassMassCOGsUnfiltered_90 <- recip_matching(MassInMassCOGs_unfiltered_dict, FindMassCOGsInMass_unfiltered_dict, MassInMassCOGs_unfiltered_seq_id_vec, 0.90)
match_MassMassCOGsUnfiltered_95 <- recip_matching(MassInMassCOGs_unfiltered_dict, FindMassCOGsInMass_unfiltered_dict, MassInMassCOGs_unfiltered_seq_id_vec, 0.95)
match_MassMassCOGsUnfiltered_99 <- recip_matching(MassInMassCOGs_unfiltered_dict, FindMassCOGsInMass_unfiltered_dict, MassInMassCOGs_unfiltered_seq_id_vec, 0.99)
match_MassCOGsMassUnfiltered_90 <- recip_matching(FindMassCOGsInMass_unfiltered_dict, MassInMassCOGs_unfiltered_dict, MassCOGsInMass_unfiltered_seq_id_vec, 0.90)
match_MassCOGsMassUnfiltered_95 <- recip_matching(FindMassCOGsInMass_unfiltered_dict, MassInMassCOGs_unfiltered_dict, MassCOGsInMass_unfiltered_seq_id_vec, 0.95)
match_MassCOGsMassUnfiltered_99 <- recip_matching(FindMassCOGsInMass_unfiltered_dict, MassInMassCOGs_unfiltered_dict, MassCOGsInMass_unfiltered_seq_id_vec, 0.99)
par(pty="s")
plot(Mass_ggC_all_gene_freqs_dict[names(match_MassMassCOGsUnfiltered_90)], Mass_cog_all_gene_freqs_dict[match_MassMassCOGsUnfiltered_90[names(match_MassMassCOGsUnfiltered_90)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity")
abline(0,1)
plot(Mass_ggC_all_gene_freqs_dict[names(match_MassMassCOGsUnfiltered_95)], Mass_cog_all_gene_freqs_dict[match_MassMassCOGsUnfiltered_95[names(match_MassMassCOGsUnfiltered_95)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 95% sequence identity")
abline(0,1)
plot(Mass_ggC_all_gene_freqs_dict[names(match_MassMassCOGsUnfiltered_99)], Mass_cog_all_gene_freqs_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 99% sequence identity")
abline(0,1)

# colur by genes that are under NFDS in Mass ggCaller
colours_MassUnfiltered_95 <- rep("grey", length(match_MassMassCOGsUnfiltered_95))
names(colours_MassUnfiltered_95) <- names(match_MassMassCOGsUnfiltered_95)
colours_MassUnfiltered_95[names(which(Mass_underNFDS[names(match_MassMassCOGsUnfiltered_95)]==1))] <- col_clb[3]

plot(Mass_ggC_all_gene_freqs_dict[names(match_MassMassCOGsUnfiltered_95)], Mass_cog_all_gene_freqs_dict[match_MassMassCOGsUnfiltered_95[names(match_MassMassCOGsUnfiltered_95)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 95% sequence identity", col = colours_MassUnfiltered_95, pch = 19)
abline(0,1)

# e.g.
# group_386 
# 0.8947368 
# from ggCaller mapped to COG
# CLS03441 
# 0 

MassCOGs_seqlengths <- read.csv("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/Mass_COGs/sequence_lengths.csv", header=TRUE)
MassCOGs_seqlengths_dict <- MassCOGs_seqlengths$Length
names(MassCOGs_seqlengths_dict) <- MassCOGs_seqlengths$GeneCluster

plot(MassCOGs_seqlengths_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]], Mass_ggC_all_gene_freqs_dict[names(match_MassMassCOGsUnfiltered_99)] - Mass_cog_all_gene_freqs_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 99% sequence identity")

MassggC_seqlengths <- read.csv("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/Massachusetts/sequence_lengths_AllMass.csv", header=TRUE)
MassggC_seqlengths_dict <- MassggC_seqlengths$Length
names(MassggC_seqlengths_dict) <- MassggC_seqlengths$GeneCluster

plot(MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)], MassCOGs_seqlengths_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]], xlab = "Mass ggCaller gene length", ylab = "Mass COG gene length",main="All Gene Frequencies, 99% sequence identity")

plot(MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)] - MassCOGs_seqlengths_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]], Mass_ggC_all_gene_freqs_dict[names(match_MassMassCOGsUnfiltered_99)] - Mass_cog_all_gene_freqs_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]], xlab = "Mass ggCaller gene length - Mass COG gene length", ylab = "Mass ggCaller gene freq - Mass COG gene freq",main="All Gene Frequencies, 99% sequence identity")

plot(MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)], col = "#56B4E980", pch = 19)
points(MassCOGs_seqlengths_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]], col = "#E69F0080", pch = 19)

plot(sort(MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)]), col = "#56B4E980", pch = 19)
points((MassCOGs_seqlengths_dict[match_MassMassCOGsUnfiltered_99[names(sort(MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)]))]]), col = "#E69F0080", pch = 19)

plot(abs(MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)] - MassCOGs_seqlengths_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]])/MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)])

# maybe filter matches by length difference?
plot(Mass_ggC_all_gene_freqs_dict[names(match_MassMassCOGsUnfiltered_99)], Mass_cog_all_gene_freqs_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 99% sequence identity")
points(Mass_ggC_all_gene_freqs_dict[names(which(abs(MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)] - MassCOGs_seqlengths_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]])/MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)] > 0.2))], Mass_cog_all_gene_freqs_dict[match_MassMassCOGsUnfiltered_99[names(which(abs(MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)] - MassCOGs_seqlengths_dict[match_MassMassCOGsUnfiltered_99[names(match_MassMassCOGsUnfiltered_99)]])/MassggC_seqlengths_dict[names(match_MassMassCOGsUnfiltered_99)] > 0.2))]], col = "red")
# hm, does not look very successful.
# maybe better change what mmseqs considers a good fit
# also, check whether you have the same issues btw ggCaller Mass and ggCaller UK (or whether this is just an issue with how I selected the COG reps)

# 15.01.2025
# using reciprocal match function of mmseqs
mmseq_results_RecipMassggCMassCOGs_unfiltered <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/RecipBestHit_ggCMassCOGMass/MassCOGggC_recipBestHit", header=TRUE)
autRecipMatch_MassCOGggC <- mmseq_results_RecipMassggCMassCOGs_unfiltered$target
names(autRecipMatch_MassCOGggC) <- mmseq_results_RecipMassggCMassCOGs_unfiltered$query

plot(sort(abs(mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen - mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen)/mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen))
plot(mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen, mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen)
abline(0,1)

par(pty="s")
plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity")
abline(0,1)
plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.9)]], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.9)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity", pch = 19)
abline(0,1)
plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.95)]], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.95)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 95% sequence identity", pch = 19)
abline(0,1)
plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.99)]], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.99)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity", pch = 19)
abline(0,1)

#library("viridis")

heat_col_vec <- heat.colors(150)
heat_col_vec_plt <- heat_col_vec[round((abs(mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen - mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen)/mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen)*100)]

plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity", col = heat_col_vec_plt, pch = 19)
abline(0,1)
plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.9)]], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.9)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity", col = heat_col_vec_plt[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.9)])
abline(0,1)
plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.95)]], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.95)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity")
abline(0,1)
plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.99)]], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.99)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity")
abline(0,1)

plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query[which(abs(mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen - mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen)/mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen < 0.01)]], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target[which(abs(mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen - mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen)/mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen < 0.01)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity", pch = 19)
abline(0,1)
# hm, no filtering for length does not seem to solve the problem

autRecipMatch_MassCOGggC[names(which(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query] > 0.95))]
# find genes that have freq of >0.95 in ggC and <0.05 in COG:
which (Mass_cog_all_gene_freqs_dict[autRecipMatch_MassCOGggC[names(which(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query] > 0.95))]] < 0.05 )
# seven
# e.g. mmseq_results_RecipMassggCMassCOGs_unfiltered[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$target == "CLS03500"),]
#   target fident alnlen qlen tlen    evalue
# CLS01535  0.964    252  255  255 2.425e-41
# CLS339898  0.978    138 1362  153 1.517e-19
# CLS06224  0.888    108  111  111 9.059e-13 (query is group_33)
# CLS03742   0.98    150  153  153 2.214e-21
# CLS04709      1    177 1095  180 5.454e-29
# CLS02134  0.999   4827 4845 4830      0
# CLS03500  0.897    234  234  234 1.103e-38
# 
# top left corner
which (Mass_cog_all_gene_freqs_dict[autRecipMatch_MassCOGggC[names(which(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query] < 0.1))]] > 0.95 )
# two
#   target fident alnlen qlen tlen    evalue
# CLS03675  0.879    174  177  177 1.911e-26
# CLS03660  0.756    111  114  111 1.61e-11

plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$evalue < 1.1e-40)]], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target[which(mmseq_results_RecipMassggCMassCOGs_unfiltered$evalue < 1.1e-40)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity", pch = 19)
abline(0,1)
#abline(v=0.95)
#abline(h=0.2)
names(which(Mass_cog_all_gene_freqs_dict[autRecipMatch_MassCOGggC[names(which(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query] > 0.95))]] < 0.2 ))
# 13 (only looking at those that I didn't with the seven)
#                                               query   target fident alnlen qlen tlen    evalue
#                                     ERS070118_02504 CLS04644      1    267  270  270 2.792e-48
#                                     ERS069964_02041 CLS04346  0.982    174  177  177 6.881e-29
#                                           group_120 CLS04219   0.98    153  180  213 3.286e-25 (query is group_120)
# ERS070168_01801~~~ERS044041_01743~~~ERS069947_01849 CLS00638  0.921    192  309  198 1.137e-29
#                                            group_69 CLS03161      1    243  309  246 3.259e-44 (query is group_69)
#                                                 ... CLS01052      1    234  396  240 1.612e-42

which (Mass_cog_all_gene_freqs_dict[autRecipMatch_MassCOGggC[names(which(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query] < 0.4))]] > 0.9 )
#                                                                 query   target fident alnlen qlen tlen     evalue
# ERS044061_00151~~~ERS044046_02303~~~ERS043891_02198~~~ERS070194_02271 CLS00178  0.986    900  903  903 1.416e-173
#                                     ERS070176_02395~~~ERS044033_02244 CLS01049  0.994    516  519  846 1.676e-102
#                                                       ERS070206_00208 CLS02387  0.998   1677 1842 1680      0
#                             ERS043950_00161~~~ERS043873_02260~~~(...) CLS02387  0.998   1677 1914 1680      0
#                   ERS044058_01722~~~ERS070061_01847~~~ERS070031_01836 CLS01448      1   1101 1104 1152 2.779e-206
#                               ERS043904_00576~~~ERS044064_00746 (...) CLS00675  0.953   1035 1266 1038 4.073e-197
#                               ERS070030_00907~~~ERS069998_00570 (...) CLS04068  0.916    432  666  429 1.419e-70
#                               ERS044147_02707~~~ERS043873_01127 (...) CLS02742  0.955    951 1173  954 7.086e-178

plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query[intersect(intersect(which(mmseq_results_RecipMassggCMassCOGs_unfiltered$alnlen/mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen > 0.9), which(mmseq_results_RecipMassggCMassCOGs_unfiltered$alnlen/mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen > 0.9)), which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.95))]], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target[intersect(intersect(which(mmseq_results_RecipMassggCMassCOGs_unfiltered$alnlen/mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen > 0.9), which(mmseq_results_RecipMassggCMassCOGs_unfiltered$alnlen/mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen > 0.9)), which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.95))]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity", pch = 19)
abline(0,1)

plot(Mass_ggC_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$query[intersect(intersect(which(mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen/mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen >= 0.8), which(mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen/mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen >= 0.8)), which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.95))]], Mass_cog_all_gene_freqs_dict[mmseq_results_RecipMassggCMassCOGs_unfiltered$target[intersect(intersect(which(mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen/mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen >= 0.8), which(mmseq_results_RecipMassggCMassCOGs_unfiltered$tlen/mmseq_results_RecipMassggCMassCOGs_unfiltered$qlen >= 0.8)), which(mmseq_results_RecipMassggCMassCOGs_unfiltered$fident >= 0.95))]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 90% sequence identity", pch = 19)
abline(0,1)

# new try: take all matches and filter them for at least 80% sequence length of one another and at least 95% quality match. Afterwards try to find the best reciprocal match.
AllMatches_MassggCinCOG <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/AllMatches_ggCMassinCOGMass/MassggCinCOG_AllMatches", header=TRUE)
AllMatches_MassCOGinggC <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/AllMatches_COGMassInggCMass/MassCOGinggC_AllMatches", header=TRUE)
# filter for fident>=0.95
AllMatches_MassggCinCOG_hq <- AllMatches_MassggCinCOG[which(AllMatches_MassggCinCOG$fident >= 0.95),]
AllMatches_MassCOGinggC_hq <- AllMatches_MassCOGinggC[which(AllMatches_MassCOGinggC$fident >= 0.95),]
# filter for sequence length matches
AllMatches_MassggCinCOG_hq_seqlen <- AllMatches_MassggCinCOG_hq[intersect(which(AllMatches_MassggCinCOG_hq$qlen/AllMatches_MassggCinCOG_hq$tlen >= 0.8), which(AllMatches_MassggCinCOG_hq$tlen/AllMatches_MassggCinCOG_hq$qlen >= 0.8)),]
AllMatches_MassCOGinggC_hq_seqlen <- AllMatches_MassCOGinggC_hq[intersect(which(AllMatches_MassCOGinggC_hq$qlen/AllMatches_MassCOGinggC_hq$tlen >= 0.8), which(AllMatches_MassCOGinggC_hq$tlen/AllMatches_MassCOGinggC_hq$qlen >= 0.8)),]
# filter for alignment length
# the same sequence seems to be aligned multiple times, just for shorter alignments
AllMatches_MassggCinCOG_hq_seqlen_alnlen <- AllMatches_MassggCinCOG_hq_seqlen[intersect(which(AllMatches_MassggCinCOG_hq_seqlen$alnlen/AllMatches_MassggCinCOG_hq_seqlen$qlen >= 0.8), which(AllMatches_MassggCinCOG_hq_seqlen$alnlen/AllMatches_MassggCinCOG_hq_seqlen$tlen >= 0.8)),]
AllMatches_MassCOGinggC_hq_seqlen_alnlen <- AllMatches_MassCOGinggC_hq_seqlen[intersect(which(AllMatches_MassCOGinggC_hq_seqlen$alnlen/AllMatches_MassCOGinggC_hq_seqlen$qlen >= 0.8), which(AllMatches_MassCOGinggC_hq_seqlen$alnlen/AllMatches_MassCOGinggC_hq_seqlen$tlen >= 0.8)),]

# investigate the correlation between fident and the difference in gene frequence
# filtered
plot(abs(Mass_ggC_all_gene_freqs_dict[AllMatches_MassggCinCOG_hq_seqlen_alnlen$query] - Mass_cog_all_gene_freqs_dict[AllMatches_MassggCinCOG_hq_seqlen_alnlen$target]),AllMatches_MassggCinCOG_hq_seqlen_alnlen$fident, col = "#00000030")
# unfiltered
plot(abs(Mass_ggC_all_gene_freqs_dict[AllMatches_MassggCinCOG$query] - Mass_cog_all_gene_freqs_dict[AllMatches_MassggCinCOG$target]),AllMatches_MassggCinCOG$fident, col = "#00000030")
# there is definitively some trend but it is not very clear at all.


# now, check whether forward and backward match is the same
AllMatches_MassggCinCOG_dict <- AllMatches_MassggCinCOG_hq_seqlen_alnlen$target
names(AllMatches_MassggCinCOG_dict) <- AllMatches_MassggCinCOG_hq_seqlen_alnlen$query

AllMatches_MassCOGinggC_dict <- AllMatches_MassCOGinggC_hq_seqlen_alnlen$target
names(AllMatches_MassCOGinggC_dict) <- AllMatches_MassCOGinggC_hq_seqlen_alnlen$query

fake_seq_id <- rep(0.96, length(names(AllMatches_MassggCinCOG_dict)))
names(fake_seq_id) <- names(AllMatches_MassggCinCOG_dict)
recip_matches_MassggCCOG <- recip_matching(AllMatches_MassggCinCOG_dict, AllMatches_MassCOGinggC_dict, fake_seq_id, 0.95)
recip_matches_MassCOGggC <- recip_matching(AllMatches_MassCOGinggC_dict, AllMatches_MassggCinCOG_dict, fake_seq_id, 0.95)

par(pty="s")
plot(Mass_ggC_all_gene_freqs_dict[names(recip_matches_MassggCCOG)], Mass_cog_all_gene_freqs_dict[recip_matches_MassggCCOG[names(recip_matches_MassggCCOG)]], xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 95% sequence identity", pch = 19)
abline(0,1)
#basically looks like before
# the alignments are longer, i.e. better, but usually the same clusters are still matched

# now look at non-recip matches
# I have 2171 rows in AllMatches_MassggCinCOG_hq_seqlen_alnlen, 2171 unique query sequences, 2136 unique target sequences
# I have 2278 rows in AllMatches_MassCOGinggC_hq_seqlen_alnlen, 2134 unique query seqs (which are the target seqs in AllMatches_MassggCinCOG_hq_seqlen_alnlen), 2168 unique target sequences (which are the query seqs in AllMatches_MassggCinCOG_hq_seqlen_alnlen)
# that matches almost perfectly.
# So I think I want to cluster based on that.

# I will just add the gene freqs together based on the mapping and then see how the correlation is
# will do it one direction first
Mass_ggC_all_gene_freqs_dict_with_matches <- Mass_ggC_all_gene_freqs_dict[unique(AllMatches_MassggCinCOG_hq_seqlen_alnlen$query)]
Mass_ggC_all_gene_freqs_dict_with_matches_matchfreqs <- rep(0, length(unique(AllMatches_MassggCinCOG_hq_seqlen_alnlen$query)))
names(Mass_ggC_all_gene_freqs_dict_with_matches_matchfreqs) <- unique(AllMatches_MassggCinCOG_hq_seqlen_alnlen$query)

match_vec <- c()
once <- c()
more_than_once <- c()

for (i in 1:length(AllMatches_MassggCinCOG_hq_seqlen_alnlen$query)) {
  query_name <- AllMatches_MassggCinCOG_hq_seqlen_alnlen$query[i]
  target_name <- AllMatches_MassggCinCOG_hq_seqlen_alnlen$target[i]
  match_name <- paste(query_name, target_name)
  if(!is.element(match_name,match_vec)){
    Mass_ggC_all_gene_freqs_dict_with_matches_matchfreqs[query_name] <- Mass_ggC_all_gene_freqs_dict_with_matches_matchfreqs[query_name] + Mass_cog_all_gene_freqs_dict[target_name]
    match_vec <- append(match_vec, match_name)
    if(is.element(query_name,once)){
      more_than_once <- append(query_name, more_than_once)
    }
    if(is.element(target_name,once)){
      more_than_once <- append(target_name, more_than_once)
    }
    once <- append(query_name, once)
    once <- append(target_name, once)
  }
  
}

length(which(abs(Mass_ggC_all_gene_freqs_dict_with_matches - Mass_ggC_all_gene_freqs_dict_with_matches_matchfreqs)>0.05))
# 106
length(which(abs(Mass_ggC_all_gene_freqs_dict_with_matches - Mass_ggC_all_gene_freqs_dict_with_matches_matchfreqs)<=0.05))
# 2065
# so most of them are actually within +-5%
# most of them ~ >95%
# this might be good enough?
par(pty="s")
plot(Mass_ggC_all_gene_freqs_dict_with_matches, Mass_ggC_all_gene_freqs_dict_with_matches_matchfreqs, pch = 19, col = "#00000030")
# overplotting!

# try mapping back and forth to form new clusters
AllGenesMass_COGggC <- c(unique(AllMatches_MassggCinCOG_hq_seqlen_alnlen$query))
removed_items <- c()
i_help <- 1
global_gene_clusters <- list()
global_gene_clusters_ggC <- list()
global_gene_clusters_COG <- list()

for (i in 1:length(AllGenesMass_COGggC)) {
  curr_gene <- AllGenesMass_COGggC[i]
  if(!is.element(curr_gene, removed_items)){
    local_genes <- c()
    local_matches <- AllMatches_MassggCinCOG_hq_seqlen_alnlen$target[which(AllMatches_MassggCinCOG_hq_seqlen_alnlen$query == curr_gene)]
    #local_re_matches <- AllMatches_MassggCinCOG_hq_seqlen_alnlen$query[which(AllMatches_MassggCinCOG_hq_seqlen_alnlen$target == local_matches)]
    local_re_matches <- c()
    for (j in 1:length(local_matches)) {
      match_gene <-local_matches[j]
    #  local_genes <- c(local_genes, match_gene)
      local_re_matches <- c(local_re_matches,AllMatches_MassggCinCOG_hq_seqlen_alnlen$query[which(AllMatches_MassggCinCOG_hq_seqlen_alnlen$target == match_gene)])
    }
    local_genes <- c(curr_gene, local_matches, local_re_matches)
    local_genes <- unique(local_genes)
    global_gene_clusters[[i_help]] <- local_genes
    global_gene_clusters_ggC[[i_help]] <- unique(c(curr_gene, local_re_matches))
    global_gene_clusters_COG[[i_help]] <- unique(c(local_matches))
    removed_items <- append(removed_items, c(curr_gene, local_re_matches))
    i_help <- i_help +1
  }
}

# print cluster sizes
cluster_sizes <- rep(0, length(global_gene_clusters))
for (i in 1:length(global_gene_clusters)) {
  cluster_sizes[i] <- length(global_gene_clusters[[i]])
}
# okay, now I have the global clusters
# now I can properly sum them and re-plot the correlation plot

#compute gene freqs
global_gene_clusters_ggC_freqs <- rep(0, length(cluster_sizes))
global_gene_clusters_COG_freqs <- rep(0, length(cluster_sizes))
for (i in 1:length(cluster_sizes)) {
  global_gene_clusters_ggC_freqs[i] <- min(1, sum(Mass_ggC_all_gene_freqs_dict[global_gene_clusters_ggC[[i]]])) # there were a couple >1 which should be split up paralogs?
  global_gene_clusters_COG_freqs[i] <- min(1,sum(Mass_cog_all_gene_freqs_dict[global_gene_clusters_COG[[i]]]))
}


#plot them!
plot(global_gene_clusters_ggC_freqs, global_gene_clusters_COG_freqs, xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 95% sequence identity")
abline(0,1)

plot(global_gene_clusters_ggC_freqs, global_gene_clusters_COG_freqs, pch = 19, col = "#00000030")
length(which(abs(global_gene_clusters_ggC_freqs - global_gene_clusters_COG_freqs)>0.05))
# 85
length(which(abs(global_gene_clusters_ggC_freqs - global_gene_clusters_COG_freqs)<=0.05))
# 2050
# again, more than 95% are less than 0.05 different.
# I think that's fine.


### find best filtering values
# tradeoff btw high quality matches vs having genes that do not match at all
# 1) 0.95, 0.8, 0.8
# 2) 0.9, 0.75, 0.75
# 3) 0.95, 0.25
# Let's compare this across datasets!
filter_matches <- function(AllMatches_data, fident_filter = 0.95, len_filter = 0.8){
  AllMatches_data_hq <- AllMatches_data[which(AllMatches_data$fident >= fident_filter),] 
  # keep only matches that have a squence identity of at least 0.95
  AllMatches_data_hq_seqlen <- AllMatches_data_hq[intersect(which(AllMatches_data_hq$qlen/AllMatches_data_hq$tlen >= len_filter), which(AllMatches_data_hq$tlen/AllMatches_data_hq$qlen >= len_filter)),]
  # keep only matches that are within 80% of each others sequence length
  AllMatches_data_hq_seqlen_alnlen <- AllMatches_data_hq_seqlen[intersect(which(AllMatches_data_hq_seqlen$alnlen/AllMatches_data_hq_seqlen$qlen >= len_filter), which(AllMatches_data_hq_seqlen$alnlen/AllMatches_data_hq_seqlen$tlen >= len_filter)),]
  # keep only matches of which at least 80% of the sequences are matched (this removes the high-quality but really short matches)
  AllMatches_data_hq_seqlen_alnlen
}

# filter matches just by best (which is the first match)
filter_matches_v2 <- function(AllMatches_data){
  Best1on1Matches <- as.data.frame(AllMatches_data[1,])
  old_query <- AllMatches_data$query[1]
  old_target <- AllMatches_data$target[1]
  new_idx <- 1
  for (i in 1:nrow(AllMatches_data)) {
    curr_query <- AllMatches_data$query[i]
    curr_target <- AllMatches_data$target[i]
    if(curr_query != old_query | curr_target != old_target){
      new_idx <- new_idx + 1
      Best1on1Matches[new_idx,] <- AllMatches_data[i,]
    }
    old_query <- AllMatches_data$query[i]
    old_target <- AllMatches_data$target[i]
  }
  Best1on1Matches
}

# Mass vs UK
AllMatches_MassInUK <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/AllMatches_MassInUK/MassInUK_AllMatches", header=TRUE)
AllMatches_UKinMass <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/AllMatches_UKinMass/UKinMass_AllMatches", header=TRUE)
# new with cov-mode 1
AllMatches_MassInUK_covmode1 <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/MMseqs2_results/AllMatches_MassInUK_covmode1/MassInUK_AllMatches", header=TRUE)
AllMatches_MassInUK_covmode1_fident70 <- AllMatches_MassInUK_covmode1[which(AllMatches_MassInUK_covmode1$fident >= 0.7),]
AllMatches_MassInUK_covmode1_Best1on1 <- filter_matches_v2(AllMatches_MassInUK_covmode1_fident70)

# filter them:
AllMatches_MassInUK_filtered <- filter_matches(AllMatches_MassInUK)
AllMatches_UKinMass_filtered <- filter_matches(AllMatches_UKinMass)

AllMatches_MassInUK_filtered <- filter_matches(AllMatches_MassInUK, fident_filter = 0.95, len_filter = 0.2)
AllMatches_UKinMass_filtered <- filter_matches(AllMatches_UKinMass, fident_filter = 0.95, len_filter = 0.2)

create_global_clusters <- function(AllMatches_data){
  AllGenesMass <- c(unique(AllMatches_data$query))
  removed_items <- c()
  i_help <- 1
  global_gene_clusters_both <- list()
  global_gene_clusters_a <- list()
  global_gene_clusters_b <- list()
  
  for (i in 1:length(AllGenesMass)) {
    curr_gene <- AllGenesMass[i]
    if(!is.element(curr_gene, removed_items)){
      local_genes <- c()
      local_matches <- AllMatches_data$target[which(AllMatches_data$query == curr_gene)]
      #local_re_matches <- AllMatches_MassggCinCOG_hq_seqlen_alnlen$query[which(AllMatches_MassggCinCOG_hq_seqlen_alnlen$target == local_matches)]
      local_re_matches <- c()
      for (j in 1:length(local_matches)) {
        match_gene <-local_matches[j]
        #  local_genes <- c(local_genes, match_gene)
        local_re_matches <- c(local_re_matches,AllMatches_data$query[which(AllMatches_data$target == match_gene)])
      }
      local_genes <- c(curr_gene, local_matches, local_re_matches)
      local_genes <- unique(local_genes)
      global_gene_clusters_both[[i_help]] <- local_genes
      global_gene_clusters_a[[i_help]] <- unique(c(curr_gene, local_re_matches))
      global_gene_clusters_b[[i_help]] <- unique(c(local_matches))
      removed_items <- append(removed_items, c(curr_gene, local_re_matches))
      i_help <- i_help +1
    }
  }
  return(list(global_gene_clusters_both, global_gene_clusters_a, global_gene_clusters_b))
}

# try this for new approach
#AllMatches_MassInUK_covmode1_Best1on1
global_clusters_return_Mass_UK <- create_global_clusters(AllMatches_MassInUK_covmode1_Best1on1)
# create cluster for Mass-UK
global_clusters_return_Mass_UK <- create_global_clusters(AllMatches_MassInUK_filtered)
global_gene_clusters_Mass_UK <- global_clusters_return_Mass_UK[[1]]
global_gene_clusters_Mass <- global_clusters_return_Mass_UK[[2]]
global_gene_clusters_UK <- global_clusters_return_Mass_UK[[3]]

#compute gene freqs
global_gene_clusters_Mass_freqs <- rep(0, length(global_gene_clusters_Mass_UK))
global_gene_clusters_UK_freqs <- rep(0, length(global_gene_clusters_Mass_UK))
for (i in 1:length(global_gene_clusters_Mass_UK)) {
  global_gene_clusters_Mass_freqs[i] <- min(1, sum(Mass_ggC_all_gene_freqs_dict[global_gene_clusters_Mass[[i]]])) # there were a couple >1 which should be split up paralogs?
  global_gene_clusters_UK_freqs[i] <- min(1,sum(UK_ggC_all_gene_freqs_dict[global_gene_clusters_UK[[i]]]))
}
plot(global_gene_clusters_Mass_freqs, global_gene_clusters_UK_freqs, pch = 19, col = "#00000030")
length(which(abs(global_gene_clusters_Mass_freqs - global_gene_clusters_UK_freqs)>0.05))
# 760
length(which(abs(global_gene_clusters_Mass_freqs - global_gene_clusters_UK_freqs)<=0.05))
# 2381
# >75% are within 0.05 of one another

# colour them in by NFDS
global_gene_clusters_Mass_NFDS <- rep(0, length(global_gene_clusters_Mass_UK))
global_gene_clusters_UK_NFDS <- rep(0, length(global_gene_clusters_Mass_UK))
for (i in 1:length(global_gene_clusters_Mass_UK)) {
  global_gene_clusters_Mass_NFDS[i] <- min(1, sum(Mass_underNFDS[global_gene_clusters_Mass[[i]]]))
  global_gene_clusters_UK_NFDS[i] <- min(1, sum(UK_underNFDS[global_gene_clusters_UK[[i]]]))
}

colours_global_gene_UKMass_95 <- rep("grey", length(global_gene_clusters_Mass_freqs))

colours_global_gene_UKMass_95[which(global_gene_clusters_UK_NFDS==1)] <- col_clb[2]
colours_global_gene_UKMass_95[which(global_gene_clusters_Mass_NFDS==1)] <- col_clb[3]
colours_global_gene_UKMass_95[intersect(which(global_gene_clusters_UK_NFDS==1),which(global_gene_clusters_Mass_NFDS==1))] <- col_clb[8]
par(pty="s")
plot(global_gene_clusters_Mass_freqs, global_gene_clusters_UK_freqs, pch = 19, col = colours_global_gene_UKMass_95)

length(which(colours_global_gene_UKMass_95 == col_clb[2]))
# 161
length(which(colours_global_gene_UKMass_95 == col_clb[3]))
# 218
length(which(colours_global_gene_UKMass_95 == col_clb[8]))
# 53

# save global cluster assignment in dictionary
Mass_ggC_intermed_gene_names <- readRDS(file = "Mass_ggC_intermed_gene_names.rds")
local_to_global_Mass <- rep(NA, length(Mass_ggC_intermed_gene_names))
names(local_to_global_Mass) <- Mass_ggC_intermed_gene_names
for (i in 1:length(global_gene_clusters_Mass_UK)) {
  for (mass_gene in global_gene_clusters_Mass[[i]]) {
    #print(mass_gene)
    if(is.element(mass_gene, Mass_ggC_intermed_gene_names)){
      local_to_global_Mass[mass_gene] <- i
    }
  }
}

length(which(!is.na(local_to_global_Mass)))
#[1] 1034
length(which(is.na(local_to_global_Mass)))
#[1] 736

### investigate hits that are very frequent in ggCaller and very rare in COGs
intersect(which(global_gene_clusters_ggC_freqs > 0.99), which(global_gene_clusters_COG_freqs < 0.05))
# [1]   27 1444
global_gene_clusters_ggC[27]
# "ERS044043_02758~~~ERS044119_02218"
View(AllMatches_MassggCinCOG_hq_seqlen_alnlen)
# match in COGs: CLS01535

# new approach: 24.01.2025
# sub-cluster clusters and then take the best reciprocal hits
BestHits_Mass_ggC_gCOG <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/Map_SubClusterReps/Mass_BestRecipHit_SubclusterReps/MassCOGggC_recipBestHit", header=TRUE)
BestHits_Mass_ggC_gCOG_filterfident <- BestHits_Mass_ggC_gCOG[which(BestHits_Mass_ggC_gCOG$fident >= 0.7),]

# try this with easy search. 
AllHits_Mass_ggC_gCOG <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/Map_SubClusterReps/Mass_easysearch_SubclusterReps/MassCOGggC_AllHits", header=TRUE)

AllHits_Mass_ggC_gCOG_filterfident <- AllHits_Mass_ggC_gCOG[which(AllHits_Mass_ggC_gCOG$fident >= 0.7),]
BestHits_Mass_ggC_gCOG_filterfident_len1 <- AllHits_Mass_ggC_gCOG_filterfident[which(AllHits_Mass_ggC_gCOG_filterfident$qlen/AllHits_Mass_ggC_gCOG_filterfident$tlen  >= 0.5),]
BestHits_Mass_ggC_gCOG_filterfident_len <- BestHits_Mass_ggC_gCOG_filterfident_len1[which(BestHits_Mass_ggC_gCOG_filterfident_len1$tlen/BestHits_Mass_ggC_gCOG_filterfident_len1$qlen  >= 0.5),]
# only keep the first hit for each query (i.e. the best)
BestHits_Mass_ggC_gCOG_filterfident <- AllHits_Mass_ggC_gCOG_filterfident
# or
#BestHits_Mass_ggC_gCOG_filterfident <- BestHits_Mass_ggC_gCOG_filterfident_len
  

besthit_Mass_ggC_gCOG_dict <- BestHits_Mass_ggC_gCOG_filterfident$target
names(besthit_Mass_ggC_gCOG_dict) <- BestHits_Mass_ggC_gCOG_filterfident$query
besthit_Mass_gCOG_ggC_dict <- BestHits_Mass_ggC_gCOG_filterfident$query
names(besthit_Mass_gCOG_ggC_dict) <- BestHits_Mass_ggC_gCOG_filterfident$target
# 6000 of 6516 have fident >= 0.7
# create dictionary that maps each sequence to the cluster it belongs to
# ggCaller
SeqsClust_Mass_ggC <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/Map_SubClusterReps/Mass_BestRecipHit_SubclusterReps/Mass_ggCaller_SeqsAndClusters.tsv", header=FALSE)
SeqsClust_Mass_ggC_dict <- SeqsClust_Mass_ggC$V2
names(SeqsClust_Mass_ggC_dict) <- SeqsClust_Mass_ggC$V1
names_which <- function(search_val, target_vec){
  names(which(target_vec == search_val))
}
Clust_Mass_ggC_list <- lapply(SeqsClust_Mass_ggC$V2, names_which, SeqsClust_Mass_ggC_dict)
names(Clust_Mass_ggC_list) <- SeqsClust_Mass_ggC$V2
# gCOG
SeqsClust_Mass_gCOG <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/Map_SubClusterReps/Mass_BestRecipHit_SubclusterReps/Mass_gCOG_SeqsAndClusters.tsv", header=FALSE)
SeqsClust_Mass_gCOG_dict <- SeqsClust_Mass_gCOG$V2
names(SeqsClust_Mass_gCOG_dict) <- SeqsClust_Mass_gCOG$V1
Clust_Mass_gCOG_list <- lapply(SeqsClust_Mass_gCOG$V2, names_which, SeqsClust_Mass_gCOG_dict)
names(Clust_Mass_gCOG_list) <- SeqsClust_Mass_gCOG$V2
# iterate over best hits
#

CreateConsensusClusters <- function(local_clust_a_list, local_clust_b_list, seq_clust_dict_a, seq_clust_dict_b, besthit_a_to_b_dict, besthit_b_to_a_dict){
  # create empty list of consensus gene clusters
  consensus_gene_clusters_both <- list()
  consensus_gene_clusters_a <- list()
  consensus_gene_clusters_b <- list()
  # save which clusters have already been clustered
  Clust_a_used <- rep(FALSE, length(local_clust_a_list))
  names(Clust_a_used) <- names(local_clust_a_list)
  #Clust_b_used <- rep(FALSE, length(local_clust_b_list))
  #names(Clust_b_used) <- names(local_clust_b_list)
  #
  #local_genes <- list()
  #local_clusters <- c()
  i_help <- 1 # keeps track of number of consensus clusters
  for (cluster_id in names(local_clust_a_list)) {
    if(!Clust_a_used[cluster_id]){
      local_clusters <- "a" # saves info that this cluster is from dataset a
      names(local_clusters) <- cluster_id # saves cluster name
      local_genes <- list(cluster_id = "") # creates empty entry for all genes in the cluster
      i <- 1
      while(i <= length(local_clusters)){
        curr_clust <- names(local_clusters)[i]
        if(local_clusters[i] == "a"){
          curr_genes <- c(setdiff(c(local_clust_a_list[[curr_clust]]), local_genes[[curr_clust]])) # genes that are in the cluster but have not been looked at yet
          #print(curr_genes)
          for (j in 1:length(curr_genes)) {
            now_gene <- curr_genes[j]
            #print(now_gene)
            if(!is.na(besthit_a_to_b_dict[now_gene])){
              new_gene <- besthit_a_to_b_dict[now_gene]
              #print(new_gene)
              new_clust <- seq_clust_dict_b[new_gene]
              #print(new_clust)
              local_clusters[new_clust] <- "b"
              local_genes[[new_clust]] <- c(local_genes[[new_clust]], new_gene)
            }
          }
        }
        else{ # is from dataset b
          curr_genes <- setdiff(local_clust_b_list[[curr_clust]], local_genes[[curr_clust]]) # genes that are in the cluster but have not been looked at yet
          for (j in curr_genes) {
            new_gene <- curr_genes[j]
            if(!is.na(besthit_b_to_a_dict[now_gene])){
              new_clust <- seq_clust_dict_a[new_gene]
              local_clusters[new_clust] <- "a"
              local_genes[[new_clust]] <- c(local_genes[[new_clust]], new_gene)
            }
          }
        }
        i <- i + 1
      }
      Clust_a_used[names(local_clusters[which(local_clusters == "a")])] <- TRUE # makes sure to not iterate over clusters later that have already been assigned to global clusters
      consensus_gene_clusters_both[[i_help]] <- names(local_clusters) # saves all clusters (from a and b) to consensus cluster
      consensus_gene_clusters_a[[i_help]] <- names(local_clusters[which(local_clusters == "a")])
      consensus_gene_clusters_b[[i_help]] <- names(local_clusters[which(local_clusters == "b")])
      i_help <- i_help + 1
    }
  }
  return(list(consensus_gene_clusters_both, consensus_gene_clusters_a, consensus_gene_clusters_b))
}

consensus_clusters_Mass_ggC_gCOG_return <- CreateConsensusClusters(local_clust_a_list = Clust_Mass_ggC_list, local_clust_b_list = Clust_Mass_gCOG_list, seq_clust_dict_a = SeqsClust_Mass_ggC_dict, seq_clust_dict_b = SeqsClust_Mass_gCOG_dict, besthit_a_to_b_dict = besthit_Mass_ggC_gCOG_dict, besthit_b_to_a_dict = besthit_Mass_gCOG_ggC_dict)

consensus_clusters_Mass_ggC_gCOG_both <- consensus_clusters_Mass_ggC_gCOG_return[[1]]
consensus_clusters_Mass_ggC <- consensus_clusters_Mass_ggC_gCOG_return[[2]]
consensus_clusters_Mass_gCOG <- consensus_clusters_Mass_ggC_gCOG_return[[3]]

# there are 1776 without any match
# 2789 with at least two clusters in one
# filter for real matches
consensus_clusters_Mass_ggC_gCOG_both_match <- consensus_clusters_Mass_ggC_gCOG_both[which(lengths(consensus_clusters_Mass_ggC_gCOG_both)>1)]
consensus_clusters_Mass_ggC_match <- consensus_clusters_Mass_ggC[which(lengths(consensus_clusters_Mass_ggC_gCOG_both)>1)]
consensus_clusters_Mass_gCOG_match <- consensus_clusters_Mass_gCOG[which(lengths(consensus_clusters_Mass_ggC_gCOG_both)>1)]

# replace "-" by "~~~"
replace_minus <- function(orig_str){
  gsub("-","~~~",orig_str, fixed = TRUE)
}

#compute gene freqs
consensus_gene_clusters_ggC_freqs <- rep(0, length(consensus_clusters_Mass_ggC_gCOG_both_match))
consensus_gene_clusters_COG_freqs <- rep(0, length(consensus_clusters_Mass_ggC_gCOG_both_match))
for (i in 1:length(consensus_clusters_Mass_ggC_gCOG_both_match)) {
  consensus_gene_clusters_ggC_freqs[i] <- sum(Mass_ggC_all_gene_freqs_dict[sapply(consensus_clusters_Mass_ggC_match[[i]], replace_minus)])
  consensus_gene_clusters_COG_freqs[i] <- sum(Mass_cog_all_gene_freqs_dict[consensus_clusters_Mass_gCOG_match[[i]]])
}


#plot them!
plot(consensus_gene_clusters_ggC_freqs, consensus_gene_clusters_COG_freqs, xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 95% sequence identity")
abline(0,1)

plot(consensus_gene_clusters_ggC_freqs, consensus_gene_clusters_COG_freqs, pch = 19, col = "#00000030")
length(which(abs(consensus_gene_clusters_ggC_freqs - consensus_gene_clusters_COG_freqs)>0.05))
# 434
length(which(abs(consensus_gene_clusters_ggC_freqs - consensus_gene_clusters_COG_freqs)<=0.05))
# 2355

CreateConsensusClusters_multipleMatches <- function(local_clust_a_list, local_clust_b_list, seq_clust_dict_a, seq_clust_dict_b, match_df){
  # create empty list of consensus gene clusters
  consensus_gene_clusters_both <- list()
  consensus_gene_clusters_a <- list()
  consensus_gene_clusters_b <- list()
  # save which clusters have already been clustered
  Clust_a_used <- rep(FALSE, length(local_clust_a_list))
  names(Clust_a_used) <- names(local_clust_a_list)
  #Clust_b_used <- rep(FALSE, length(local_clust_b_list))
  #names(Clust_b_used) <- names(local_clust_b_list)
  #
  #local_genes <- list()
  #local_clusters <- c()
  i_help <- 1 # keeps track of number of consensus clusters
  for (cluster_id in names(local_clust_a_list)) {
    if(!Clust_a_used[cluster_id]){
      local_clusters <- "a" # saves info that this cluster is from dataset a
      names(local_clusters) <- cluster_id # saves cluster name
      local_genes <- list(cluster_id = "") # creates empty entry for all genes in the cluster
      i <- 1
      while(i <= length(local_clusters)){
        curr_clust <- names(local_clusters)[i]
        if(local_clusters[i] == "a"){
          curr_genes <- c(setdiff(c(local_clust_a_list[[curr_clust]]), local_genes[[curr_clust]])) # genes that are in the cluster but have not been looked at yet
          #print(curr_genes)
          for (j in 1:length(curr_genes)) {
            now_gene <- curr_genes[j]
            #print(now_gene)
            if(length(which(match_df$query == now_gene)) != 0){
              new_genes <- match_df$target[which(match_df$query == now_gene)]
              for (k in 1:length(new_genes)) {
                new_gene <- new_genes[k]
                #print(new_gene)
                new_clust <- seq_clust_dict_b[new_gene]
                #print(new_clust)
                local_clusters[new_clust] <- "b"
                local_genes[[new_clust]] <- c(local_genes[[new_clust]], new_gene)
              }
            }
          }
        }
        else{ # is from dataset b
          curr_genes <- setdiff(local_clust_b_list[[curr_clust]], local_genes[[curr_clust]]) # genes that are in the cluster but have not been looked at yet
          for (j in curr_genes) {
            now_gene <- curr_genes[j]
            if(length(which(match_df$target == now_gene)) != 0){
              new_genes <- match_df$query[which(match_df$target == now_gene)]
              for (k in 1:length(new_genes)) {
                new_gene <- new_genes[k]
                new_clust <- seq_clust_dict_a[new_gene]
                local_clusters[new_clust] <- "a"
                local_genes[[new_clust]] <- c(local_genes[[new_clust]], new_gene)
              }
            }
          }
        }
        i <- i + 1
      }
      Clust_a_used[names(local_clusters[which(local_clusters == "a")])] <- TRUE # makes sure to not iterate over clusters later that have already been assigned to global clusters
      consensus_gene_clusters_both[[i_help]] <- names(local_clusters) # saves all clusters (from a and b) to consensus cluster
      consensus_gene_clusters_a[[i_help]] <- names(local_clusters[which(local_clusters == "a")])
      consensus_gene_clusters_b[[i_help]] <- names(local_clusters[which(local_clusters == "b")])
      i_help <- i_help + 1
    }
  }
  return(list(consensus_gene_clusters_both, consensus_gene_clusters_a, consensus_gene_clusters_b))
}
consensus_clusters_Mass_ggC_gCOG_return <- CreateConsensusClusters_multipleMatches(local_clust_a_list = Clust_Mass_ggC_list, local_clust_b_list = Clust_Mass_gCOG_list, seq_clust_dict_a = SeqsClust_Mass_ggC_dict, seq_clust_dict_b = SeqsClust_Mass_gCOG_dict, match_df = BestHits_Mass_ggC_gCOG_filterfident)

consensus_clusters_Mass_ggC_gCOG_both <- consensus_clusters_Mass_ggC_gCOG_return[[1]]
consensus_clusters_Mass_ggC <- consensus_clusters_Mass_ggC_gCOG_return[[2]]
consensus_clusters_Mass_gCOG <- consensus_clusters_Mass_ggC_gCOG_return[[3]]

consensus_clusters_Mass_ggC_gCOG_both_match <- consensus_clusters_Mass_ggC_gCOG_both[which(lengths(consensus_clusters_Mass_ggC_gCOG_both)>1)]
consensus_clusters_Mass_ggC_match <- consensus_clusters_Mass_ggC[which(lengths(consensus_clusters_Mass_ggC_gCOG_both)>1)]
consensus_clusters_Mass_gCOG_match <- consensus_clusters_Mass_gCOG[which(lengths(consensus_clusters_Mass_ggC_gCOG_both)>1)]

#compute gene freqs
consensus_gene_clusters_ggC_freqs <- rep(0, length(consensus_clusters_Mass_ggC_gCOG_both_match))
consensus_gene_clusters_COG_freqs <- rep(0, length(consensus_clusters_Mass_ggC_gCOG_both_match))
for (i in 1:length(consensus_clusters_Mass_ggC_gCOG_both_match)) {
  consensus_gene_clusters_ggC_freqs[i] <- sum(Mass_ggC_all_gene_freqs_dict[sapply(consensus_clusters_Mass_ggC_match[[i]], replace_minus)])
  consensus_gene_clusters_COG_freqs[i] <- sum(Mass_cog_all_gene_freqs_dict[consensus_clusters_Mass_gCOG_match[[i]]])
}


#plot them!
plot(consensus_gene_clusters_ggC_freqs, consensus_gene_clusters_COG_freqs, xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 95% sequence identity")
abline(0,1)

plot(consensus_gene_clusters_ggC_freqs, consensus_gene_clusters_COG_freqs, pch = 19, col = "#00000030")
length(which(abs(consensus_gene_clusters_ggC_freqs - consensus_gene_clusters_COG_freqs)>0.05))
# 332
length(which(abs(consensus_gene_clusters_ggC_freqs - consensus_gene_clusters_COG_freqs)<=0.05))
# 1409

# analyse issue based on previous example, ERS044043_02758-ERS044119_02218 vs CLS01535
Clust_Mass_ggC_list[["ERS044043_02758-ERS044119_02218"]]
#[1] "6925_1#49_1651" "6930_8#9_8492"  "7622_2#29_602" 
# all hits for these three are filtered out by the 50% length requirement

# I would like to not filter by length hits. But if I do not filter by it, the data frame has more than 600,000 rows. That's a bit much.
# Instead, I will check which clusters all of the entries belong to and then filter by unique cluster matches.
# translate AllHits to clusters
AllHits_Mass_ggC_gCOG <- read.delim("/Users/llorenz/Documents/PhD_Project/Data/Mapping_ggCaller/Map_SubClusterReps/Mass_easysearch_SubclusterReps/MassCOGggC_AllHits", header=TRUE)

AllHits_Mass_ggC_gCOG_filterfident <- AllHits_Mass_ggC_gCOG[which(AllHits_Mass_ggC_gCOG$fident >= 0.7),]
#AllHits_Mass_ggC_gCOG_filterfident_len1 <- AllHits_Mass_ggC_gCOG_filterfident[which(AllHits_Mass_ggC_gCOG_filterfident$qlen/AllHits_Mass_ggC_gCOG_filterfident$tlen  >= 0.7),]
#AllHits_Mass_ggC_gCOG_filterfident_len <- AllHits_Mass_ggC_gCOG_filterfident_len1[which(AllHits_Mass_ggC_gCOG_filterfident_len1$tlen/AllHits_Mass_ggC_gCOG_filterfident_len1$qlen  >= 0.7),]
#AllHits_Mass_ggC_gCOG_filterfident_len1 <- AllHits_Mass_ggC_gCOG[which(AllHits_Mass_ggC_gCOG_filterfident$alnlen/AllHits_Mass_ggC_gCOG_filterfident$tlen  >= 0.8),]
#AllHits_Mass_ggC_gCOG_filterfident_len <- AllHits_Mass_ggC_gCOG[which(AllHits_Mass_ggC_gCOG_filterfident_len1$alnlen/AllHits_Mass_ggC_gCOG_filterfident_len1$qlen  >= 0.8),]
# no genes within 70% of each other's length!! (qlen/tlen >= 0.7 and tlen/qlen>=0.7)

#AllHits_Mass_ggC_gCOG_filterfident <- AllHits_Mass_ggC_gCOG_filterfident_len
AllHits_Mass_ggC_gCOG_filterfident$query_clust <- SeqsClust_Mass_ggC_dict[AllHits_Mass_ggC_gCOG_filterfident$query]
AllHits_Mass_ggC_gCOG_filterfident$target_clust <- SeqsClust_Mass_gCOG_dict[AllHits_Mass_ggC_gCOG_filterfident$target]
AllHits_Mass_ggC_gCOG_filterfident$combined_clust <- paste(AllHits_Mass_ggC_gCOG_filterfident$query_clust, AllHits_Mass_ggC_gCOG_filterfident$target_clust)
length(AllHits_Mass_ggC_gCOG_filterfident$combined_clust)
#[1] 611133
# filtered by length 50% 20050
length(unique(AllHits_Mass_ggC_gCOG_filterfident$combined_clust))
#[1] 5843
# filtered by length 50% 1743
AllHits_Mass_ggC_gCOG_filterfident_unique <- AllHits_Mass_ggC_gCOG_filterfident[!duplicated(AllHits_Mass_ggC_gCOG_filterfident$combined_clust), ]

CreateConsensusClusters_bycluster <- function(match_df){
  # create empty list of consensus gene clusters
  consensus_gene_clusters_both <- list()
  consensus_gene_clusters_a <- list()
  consensus_gene_clusters_b <- list()
  # save which clusters have already been clustered
  Clust_a_used <- rep(FALSE, nrow(match_df))
  names(Clust_a_used) <- match_df$query_clust
  #Clust_b_used <- rep(FALSE, length(local_clust_b_list))
  #names(Clust_b_used) <- names(local_clust_b_list)
  #
  #local_genes <- list()
  #local_clusters <- c()
  i_help <- 1 # keeps track of number of consensus clusters
  for (cluster_id in match_df$query_clust) {
    if(!Clust_a_used[cluster_id]){
      local_clusters <- "a" # saves info that this cluster is from dataset a
      names(local_clusters) <- cluster_id # saves cluster name
      #local_genes <- list(cluster_id = "") # creates empty entry for all genes in the cluster
      i <- 1
      while(i <= length(local_clusters)){
        curr_clust <- names(local_clusters)[i]
        if(local_clusters[i] == "a"){
          #curr_genes <- c(setdiff(c(local_clust_a_list[[curr_clust]]), local_genes[[curr_clust]])) # genes that are in the cluster but have not been looked at yet
          #print(curr_genes)
          if(length(which(match_df$query_clust == curr_clust)) != 0){
            new_clusters <- match_df$target_clust[which(match_df$query_clust == curr_clust)]
            for (k in 1:length(new_clusters)) {
              new_clust <- new_clusters[k]
              #print(new_gene)
              #new_clust <- seq_clust_dict_b[new_gene]
              #print(new_clust)
              local_clusters[new_clust] <- "b"
              #local_genes[[new_clust]] <- c(local_genes[[new_clust]], new_gene)
            }
          }
        }
        else{ # is from dataset b
          if(length(which(match_df$target_clust == curr_clust)) != 0){
            new_clusters <- match_df$query_clust[which(match_df$target_clust == curr_clust)]
            for (k in 1:length(new_clusters)) {
              new_clust <- new_clusters[k]
              #print(new_gene)
              #new_clust <- seq_clust_dict_b[new_gene]
              #print(new_clust)
              local_clusters[new_clust] <- "a"
              #local_genes[[new_clust]] <- c(local_genes[[new_clust]], new_gene)
            }
          }
        }
        i <- i + 1
      }
      Clust_a_used[names(local_clusters[which(local_clusters == "a")])] <- TRUE # makes sure to not iterate over clusters later that have already been assigned to global clusters
      #print(local_clusters)
      consensus_gene_clusters_both[[i_help]] <- names(local_clusters) # saves all clusters (from a and b) to consensus cluster
      consensus_gene_clusters_a[[i_help]] <- names(local_clusters[which(local_clusters == "a")])
      consensus_gene_clusters_b[[i_help]] <- names(local_clusters[which(local_clusters == "b")])
      i_help <- i_help + 1
    }
  }
  return(list(consensus_gene_clusters_both, consensus_gene_clusters_a, consensus_gene_clusters_b))
}

consensus_clusters_Mass_ggC_gCOG_return <- CreateConsensusClusters_bycluster(match_df = AllHits_Mass_ggC_gCOG_filterfident_unique)

consensus_clusters_Mass_ggC_gCOG_both <- consensus_clusters_Mass_ggC_gCOG_return[[1]]
consensus_clusters_Mass_ggC <- consensus_clusters_Mass_ggC_gCOG_return[[2]]
consensus_clusters_Mass_gCOG <- consensus_clusters_Mass_ggC_gCOG_return[[3]]

consensus_clusters_Mass_ggC_gCOG_both_match <- consensus_clusters_Mass_ggC_gCOG_both[which(lengths(consensus_clusters_Mass_ggC_gCOG_both)>1)]
consensus_clusters_Mass_ggC_match <- consensus_clusters_Mass_ggC[which(lengths(consensus_clusters_Mass_ggC_gCOG_both)>1)]
consensus_clusters_Mass_gCOG_match <- consensus_clusters_Mass_gCOG[which(lengths(consensus_clusters_Mass_ggC_gCOG_both)>1)]

#compute gene freqs
consensus_gene_clusters_ggC_freqs <- rep(0, length(consensus_clusters_Mass_ggC_gCOG_both_match))
consensus_gene_clusters_COG_freqs <- rep(0, length(consensus_clusters_Mass_ggC_gCOG_both_match))
for (i in 1:length(consensus_clusters_Mass_ggC_gCOG_both_match)) {
  consensus_gene_clusters_ggC_freqs[i] <- min(1,sum(Mass_ggC_all_gene_freqs_dict[sapply(consensus_clusters_Mass_ggC_match[[i]], replace_minus)]))
  consensus_gene_clusters_COG_freqs[i] <- min(1,sum(Mass_cog_all_gene_freqs_dict[consensus_clusters_Mass_gCOG_match[[i]]]))
}


#plot them!
par(pty="s")
plot(consensus_gene_clusters_ggC_freqs, consensus_gene_clusters_COG_freqs, xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 95% sequence identity", ylim = c(0,1), xlim = c(0,1))
abline(0,1)
plot(consensus_gene_clusters_ggC_freqs, consensus_gene_clusters_COG_freqs, pch = 19, col = "#00000030")
plot(consensus_gene_clusters_ggC_freqs, consensus_gene_clusters_COG_freqs, pch = 19, col = "#00000030",ylim = c(0,1), xlim = c(0,1))
abline(0,1)
length(which(abs(consensus_gene_clusters_ggC_freqs - consensus_gene_clusters_COG_freqs)>0.05))
# 105
length(which(abs(consensus_gene_clusters_ggC_freqs - consensus_gene_clusters_COG_freqs)<=0.05))
# 1378

#intersect(which(consensus_gene_clusters_ggC_freqs > 0.95), which(consensus_gene_clusters_COG_freqs < 0.05))

# new way of calculating gene_frews
Mass_ggCaller_gene_presence_absence_2001 <- readRDS(file = "ggCaller_gene_presence_absence_2001.rds")
Mass_COG_gene_presence_absence_2001 <- readRDS(file = "COG_Mass_gene_presence_absence_2001.rds")
# create a list dictionary with cluster names and presence-absence vectors as content
Mass_ggCaller_gene_presence_absence_2001_listdict <- list()
for (i in 2:nrow(Mass_ggCaller_gene_presence_absence_2001)) {
  Mass_ggCaller_gene_presence_absence_2001_listdict[[Mass_ggCaller_gene_presence_absence_2001[i,1]]] <- as.integer(Mass_ggCaller_gene_presence_absence_2001[i,-1])
  
}
Mass_COG_gene_presence_absence_2001_listdict <- list()
for (i in 2:nrow(Mass_COG_gene_presence_absence_2001)) {
  Mass_COG_gene_presence_absence_2001_listdict[[Mass_COG_gene_presence_absence_2001[i,1]]] <- as.integer(Mass_COG_gene_presence_absence_2001[i,-1])
  
}

consensus_gene_clusters_ggC_freqs_v2 <- rep(0, length(consensus_clusters_Mass_ggC_gCOG_both_match))
consensus_gene_clusters_COG_freqs_v2 <- rep(0, length(consensus_clusters_Mass_ggC_gCOG_both_match))
for (i in 1:length(consensus_clusters_Mass_ggC_gCOG_both_match)) {
  cluster_names_ggC <- sapply(consensus_clusters_Mass_ggC_match[[i]], replace_minus)
  cluster_presence_vec <- rep(0, length(Mass_ggCaller_gene_presence_absence_2001_listdict[[1]]))
  for (j in cluster_names_ggC) {
    cluster_presence_vec <- (cluster_presence_vec | Mass_ggCaller_gene_presence_absence_2001_listdict[[j]])
  }
  cluster_presence_vec <- as.integer(cluster_presence_vec)
  consensus_gene_clusters_ggC_freqs_v2[i] <- sum(cluster_presence_vec)/length(cluster_presence_vec)
  cluster_names_COG <- consensus_clusters_Mass_gCOG_match[[i]]
  cluster_presence_vec_COG <- rep(0, length(Mass_COG_gene_presence_absence_2001_listdict[[1]]))
  for (k in cluster_names_COG) {
    cluster_presence_vec_COG <- cluster_presence_vec_COG | Mass_COG_gene_presence_absence_2001_listdict[[k]]
  }
  cluster_presence_vec_COG <- as.integer(cluster_presence_vec_COG)
  consensus_gene_clusters_COG_freqs_v2[i] <- sum(cluster_presence_vec_COG)/length(cluster_presence_vec_COG)
}

#plot them!
par(pty="s")
plot(consensus_gene_clusters_ggC_freqs_v2, consensus_gene_clusters_COG_freqs_v2, xlab = "Mass ggCaller gene frequencies", ylab = "Mass COG gene frequencies",main="All Gene Frequencies, 95% sequence identity", ylim = c(0,1), xlim = c(0,1))
abline(0,1)
plot(consensus_gene_clusters_ggC_freqs_v2, consensus_gene_clusters_COG_freqs_v2, pch = 19, col = "#00000030")
plot(consensus_gene_clusters_ggC_freqs_v2, consensus_gene_clusters_COG_freqs_v2, pch = 19, col = "#00000030",ylim = c(0,1), xlim = c(0,1))
abline(0,1)
length(which(abs(consensus_gene_clusters_ggC_freqs_v2 - consensus_gene_clusters_COG_freqs_v2)>0.05))
# 58
length(which(abs(consensus_gene_clusters_ggC_freqs_v2 - consensus_gene_clusters_COG_freqs_v2)<=0.05))
# 806

plot(consensus_gene_clusters_ggC_freqs_v2, consensus_gene_clusters_ggC_freqs, pch = 19, col = "#00000030")
plot(consensus_gene_clusters_COG_freqs_v2, consensus_gene_clusters_COG_freqs, pch = 19, col = "#00000030")

length(which(lengths(consensus_clusters_Mass_ggC_gCOG_both_match)>2))
# 372


hist(table(cut(consensus_gene_clusters_ggC_freqs_v2, breaks = c(0,0.2,0.4,0.6,),include.lowest = T)),10, plot=F)

#intersect(which(consensus_gene_clusters_ggC_freqs_v2 > 0.9),which(consensus_gene_clusters_COG_freqs_v2 < 0.1))
#[1]  191  461  637 1053 1313 1631 2086
# > consensus_clusters_Mass_ggC_match[[461]]
#[1] "ERS043853_02646-ERS044068_02410" "ERS044131_01836"                
#> Mass_ggC_all_gene_freqs_dict["ERS043853_02646~~~ERS044068_02410"]
#ERS043853_02646~~~ERS044068_02410 
#0.04511278 
#> Mass_ggC_all_gene_freqs_dict["ERS044131_01836" ]
#ERS044131_01836 
#1 
Clust_Mass_ggC_list[["ERS043853_02646-ERS044068_02410"]]
# [1] "6925_1#66_4548"    "7553_4#39_10923"   "7622_2#12_-156897" "7622_3#51_-173523"
AllHits_Mass_ggC_gCOG_filterfident$qlen[which(AllHits_Mass_ggC_gCOG_filterfident$query=="6925_1#66_4548")]
#52
AllHits_Mass_ggC_gCOG_filterfident$qlen[which(AllHits_Mass_ggC_gCOG_filterfident$query=="7553_4#39_10923")]
# 79
AllHits_Mass_ggC_gCOG_filterfident$qlen[which(AllHits_Mass_ggC_gCOG_filterfident$query=="7622_2#12_-156897")]
# 54
AllHits_Mass_ggC_gCOG_filterfident$qlen[which(AllHits_Mass_ggC_gCOG_filterfident$query=="7622_3#51_-173523")]
# 43
Clust_Mass_ggC_list[["ERS044131_01836"]]
AllHits_Mass_ggC_gCOG_filterfident$qlen[which(AllHits_Mass_ggC_gCOG_filterfident$query=="6925_1#49_-138")]
#56
AllHits_Mass_ggC_gCOG_filterfident$qlen[which(AllHits_Mass_ggC_gCOG_filterfident$query=="7001_1#9_-103094")]
# filtered out?
AllHits_Mass_ggC_gCOG_filterfident$qlen[which(AllHits_Mass_ggC_gCOG_filterfident$query=="7553_4#44_7068")]
# 55
AllHits_Mass_ggC_gCOG_filterfident$qlen[which(AllHits_Mass_ggC_gCOG_filterfident$query=="7622_2#25_-161860")]
# filtered out?
AllHits_Mass_ggC_gCOG_filterfident$qlen[which(AllHits_Mass_ggC_gCOG_filterfident$query=="7622_4#19_-191370")]
# 56