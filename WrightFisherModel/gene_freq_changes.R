gene_freq_df_1 <- data.frame(cbind(mass_gene_freq_0/length(which(mass_data$Time==0)), mass_gene_freq_36/length(which(mass_data$Time==36)), mass_gene_freq_72/length(which(mass_data$Time==72))))
colnames(gene_freq_df) <- 1:3

gene_freq_df <- data.frame(x = c(rep(1,length(mass_gene_freq_0)), rep(2,length(mass_gene_freq_0)),rep(3,length(mass_gene_freq_0))),
                           y = c(mass_gene_freq_0, mass_gene_freq_36, mass_gene_freq_72),
                           z = rep(names(mass_gene_freq_0),3))
gene_freq_df_test <- rbind(gene_freq_df[1:10,], gene_freq_df[(length(mass_gene_freq_0)+1):(length(mass_gene_freq_0)+10),],gene_freq_df[(2*length(mass_gene_freq_0)+1):(2*length(mass_gene_freq_0)+10),])


ggplot(data = gene_freq_df_test, aes(x=x, y=y)) + geom_line(aes(colour=z))

par(mfrow=c(1,1))
plot(0,0,xlim = c(1,3),ylim = c(0,1),type = "n")
cl <- rainbow(5)
for (i in 1:length(mass_gene_freq_0)){
  points(1:3,gene_freq_df_1[i,],col = "black",type = 'b')
}



#plot(1:3, gene_freq_df[1,], type = "b")
#for (i in 2:nrow(gene_freq_df)){
#  lines(gene_freq_df[i,], pch=20)
#}


#ggplot(gene_freq_df, aes(y=cluster_sizes, x = as.integer(row.names(df_cluster_sizes)))) +
#  geom_bar(position="dodge", stat="identity", colour="black", width = 1,  fill="#009E73")
