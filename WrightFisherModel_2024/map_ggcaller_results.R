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
    if(name != return_val){
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
Mass_intermed_gene_underNFDS[which((Mass_delta_ranking_names <= 0.2797 * length(Mass_delta_ranking_names)))] <- 1
#if(Mass_delta_ranking <= 0.2797 * length(Mass_delta_ranking)) 1 else 0
UK_delta_ranking_names <- UK_delta_ranking
names(UK_delta_ranking_names) <- UK_ggC_intermed_gene_names
UK_intermed_gene_underNFDS <- rep(0, length(UK_delta_ranking_names))
names(UK_intermed_gene_underNFDS) <- names(UK_delta_ranking_names)
UK_intermed_gene_underNFDS[which((UK_delta_ranking_names <= 0.1432 * length(UK_delta_ranking_names)))] <- 1
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
colours_UKMassUnfiltered_90 <- rep("grey", length(match_UKMassUnfiltered_90))
names(colours_UKMassUnfiltered_90) <- names(match_UKMassUnfiltered_90)
colours_UKMassUnfiltered_90[names(which(UK_underNFDS[names(match_UKMassUnfiltered_90)]==1))] <- "blue"
colours_UKMassUnfiltered_90[match_MassUKUnfiltered_90[names(which(Mass_underNFDS[match_UKMassUnfiltered_90]==1))]] <- "red"
colours_UKMassUnfiltered_90[intersect(names(which(UK_underNFDS[names(match_UKMassUnfiltered_90)]==1)),match_MassUKUnfiltered_90[names(which(Mass_underNFDS[match_UKMassUnfiltered_90]==1))])] <- "purple"

colours_UKMassUnfiltered_95 <- rep("grey", length(match_UKMassUnfiltered_95))
names(colours_UKMassUnfiltered_95) <- names(match_UKMassUnfiltered_95)
colours_UKMassUnfiltered_95[names(which(UK_underNFDS[names(match_UKMassUnfiltered_95)]==1))] <- "blue"
colours_UKMassUnfiltered_95[match_MassUKUnfiltered_95[names(which(Mass_underNFDS[match_UKMassUnfiltered_95]==1))]] <- "red"
colours_UKMassUnfiltered_95[intersect(names(which(UK_underNFDS[names(match_UKMassUnfiltered_95)]==1)),match_MassUKUnfiltered_95[names(which(Mass_underNFDS[match_UKMassUnfiltered_95]==1))])] <- "purple"

colours_UKMassUnfiltered_99 <- rep("grey", length(match_UKMassUnfiltered_99))
names(colours_UKMassUnfiltered_99) <- names(match_UKMassUnfiltered_99)
colours_UKMassUnfiltered_99[names(which(UK_underNFDS[names(match_UKMassUnfiltered_99)]==1))] <- "blue"
colours_UKMassUnfiltered_99[match_MassUKUnfiltered_99[names(which(Mass_underNFDS[match_UKMassUnfiltered_99]==1))]] <- "red"
colours_UKMassUnfiltered_99[intersect(names(which(UK_underNFDS[names(match_UKMassUnfiltered_99)]==1)),match_MassUKUnfiltered_99[names(which(Mass_underNFDS[match_UKMassUnfiltered_99]==1))])] <- "purple"

par(pty="s")
plot(UK_ggC_all_gene_freqs_dict[names(match_UKMassUnfiltered_90)], Mass_ggC_all_gene_freqs_dict[match_UKMassUnfiltered_90[names(match_UKMassUnfiltered_90)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="All Gene Frequencies, 90% sequence identity", col = colours_UKMassUnfiltered_90, pch = 19)
abline(0,1)
plot(UK_ggC_all_gene_freqs_dict[names(match_UKMassUnfiltered_95)], Mass_ggC_all_gene_freqs_dict[match_UKMassUnfiltered_95[names(match_UKMassUnfiltered_95)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="All Gene Frequencies, 95% sequence identity", col = colours_UKMassUnfiltered_95, pch = 19)
abline(0,1)
plot(UK_ggC_all_gene_freqs_dict[names(match_UKMassUnfiltered_99)], Mass_ggC_all_gene_freqs_dict[match_UKMassUnfiltered_99[names(match_UKMassUnfiltered_99)]], xlab = "UK gene frequencies", ylab = "Mass gene frequencies",main="All Gene Frequencies, 99% sequence identity", col = colours_UKMassUnfiltered_99, pch = 19)
abline(0,1)
