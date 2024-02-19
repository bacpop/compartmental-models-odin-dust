test_variety <- intermed_mass_data[,-(1:5)]

always_together <- matrix(1, nrow = (ncol(test_variety)), ncol = (ncol(test_variety)))

for (j in 1:nrow(intermed_mass_data)){
  curr_vec <- which(test_variety[j,]==1)
  for (i in 1:length(curr_vec)) {
    always_together[curr_vec[i],] <- always_together[curr_vec[i],] * as.numeric((test_variety)[j,])
  }
}
rowSums(always_together)
which(always_together[1,]==1)
which(always_together[2,]==1)
which(always_together[3,]==1)

df_always_together <- as.data.frame(always_together)
cp_always_together <- always_together

for (i in 1:ncol(always_together)){
  for (j in (i):nrow(always_together)){
    if(always_together[j,i]!=always_together[i,j]){
      cp_always_together[i,j] <- 0
      cp_always_together[j,i] <- 0
    }
    cp_always_together[i,i] <- 0
  }
}
rowSums(cp_always_together)
colSums(cp_always_together)
which(cp_always_together[10,]==1)
which(cp_always_together[20,]==1)
which(cp_always_together[3,]==1)

unidentified <- c()
for (i in 1:nrow(cp_always_together)){
  if(!(i %in% unidentified)){
    unidentified <- append(unidentified, which(cp_always_together[i,]==1))
  }
}
unidentified <- unique(unidentified)
length(unidentified)
# [1] 352
# this means that 352 genes could be removed from the calculation because they are always present together with some other gene(s)
# this means that they cannot be identified as causal or just carried along
# I mean this is just 352 of 1062 but still (it's a third)
# could be helpful to make the functional analysis more feasible
sort(unidentified)
