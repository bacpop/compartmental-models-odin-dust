### this is just a test script that I used to have a look at the properties of multiple poisson and binomial draws
# maybe go back to this if model behavior of poisson or multinomial model is weird?
### will test visually if rbinom(n, p) multiple draws are multinomially distributed



pop_size <- 400
probs <- rep(0.1,10)
#y <- probs/pop_size
species_no <- 10

# here, I correct n and p for each draw based on the draws before and the probabilities of the other components
maybe_mult_v1<- function(repeats, pop_size, probs){
  species_no <- length(probs)
  maybe_multi <- matrix(0, species_no, repeats) # empty matrix to store the results of the random draws
  y <- rep(0, species_no) # empty vector for random draws
  
  for (j in 1:repeats) {
    y[1] <- if (probs[1]/sum(probs[1:species_no]) < 1)
      rbinom(1,pop_size, probs[1]/sum(probs[1:species_no])) else pop_size
    for (k in 2:10) {
      y[k] <- if (probs[k]/sum(probs[k:species_no]) < 1)
        rbinom(1,pop_size - sum(y[1:(k-1)]), probs[k]/sum(probs[k:species_no])) else pop_size - sum(y[1:(k-1)])
    }
    maybe_multi[,j] <- y
  }
  maybe_multi
}

# here, I do not correct n or p. I just normalise p so that the sum is 1. 
maybe_mult_v2<- function(repeats, pop_size, probs){
  species_no <- length(probs)
  maybe_multi <- matrix(0, species_no, repeats) # empty matrix to store the results of the random draws
  y <- rep(0, species_no) # empty vector for random draws
  
  for (j in 1:repeats) {
    for (k in 1:species_no) {
      y[k] <- if (probs[k]/sum(probs[1:species_no]) < 1)
      rbinom(1,pop_size, probs[k]/sum(probs[1:species_no])) else pop_size
    }
    maybe_multi[,j] <- y
  }
  maybe_multi
}


maybe_multi_version1 <- maybe_mult_v1(100000,400,rep(0.1,10))
maybe_multi_version2 <- maybe_mult_v2(100000,400,rep(0.1,10))

# plot function to visually compare the components to binomial distribution
plot_binomial <- function(multinom_matrix, index, repeats, pop_size){
  maybe_mult_slice <- multinom_matrix[index,]
  rel_maybe_mult <- rep(0, pop_size)
  for (i in 1:repeats){
    rel_maybe_mult[maybe_mult_slice[i]] <- rel_maybe_mult[maybe_mult_slice[i]] + 1
  }
  rel_maybe_mult <- rel_maybe_mult/repeats
  
  plot(1:pop_size, rel_maybe_mult)
  lines(dbinom(0:pop_size, pop_size, .1), col = "red")
}

# both look binomially distributed, I think
plot_binomial(maybe_multi_version1, 10, 100000,400)
plot_binomial(maybe_multi_version2, 10, 100000,400)

# for version 1, the sum of all draws is always 400, i.e. the pop_size is constant
plot(1:100000,colSums(maybe_multi_version1))
# for version 2, the pop_size varies. 
plot(1:100000,(colSums(maybe_multi_version2)))
# But it looks normally distributed. Let's check that:
rel_colSums <- rep(0, 1000)
abs_colSums <- colSums(maybe_multi_version2)
for (i in 1:100000){
  rel_colSums[abs_colSums[i]] <- rel_colSums[abs_colSums[i]] + 1
}
rel_colSums <- rel_colSums/100000
plot(1:1000,rel_colSums)
lines(dnorm(0:1000, mean=400, sd=sqrt(400)), col = "red")
# also the shapiro wilk test suggests that the population size with version 2 is normally distributed
shapiro.test(rel_colSums)



#########################
# Let's compare all that to the behaviour of poisson draws
poisson_population <- function(repeats, pop_size, probs, capacity){
  probs_var <- probs
  Pop_size <- pop_size
  pop_size_over_time <- rep(0,repeats)
  species_no <- length(probs)
  poisson_matrix <- matrix(0, species_no, repeats) # empty matrix to store the results of the random draws
  y <- rep(0, species_no) # empty vector for random draws
  for (j in 1:repeats) {
    for (k in 1:10) {
      y[k] <- if ((capacity/Pop_size)*probs[k]/sum(probs[1:species_no]) < 1)
        rpois(1,(capacity/Pop_size)*(probs[k]/sum(probs[1:species_no]))*Pop_size) else Pop_size 
    }
    Pop_size <- sum(y)
    pop_size_over_time[j] <- Pop_size
    probs_var <- y/Pop_size
    poisson_matrix[,j] <- y
  }
  return_list <- list(poisson_matrix, pop_size_over_time)
  return_list
}

poisson_pop_return <- poisson_population(100000, 400, rep(0.1,10), 400)
poisson_data <- poisson_pop_return[[1]]
pop_size_vector <- poisson_pop_return[[2]]
plot(1:100000, pop_size_vector)

poisson_data_slice <- poisson_data[1,]
rel_poisson <- rep(0, 1500)
for (i in 1:100000){
  rel_poisson[poisson_data_slice[i]] <- rel_poisson[poisson_data_slice[i]] + 1
}
rel_poisson <- rel_poisson/100000

plot(1:1500, rel_poisson)
# the population size is HIGHLY variable with the poisson draws. Maybe it get better when we have logistic growth?

poisson_population_v2 <- function(repeats, pop_size, probs){
  #Pop_size <- pop_size
  probs_var <- probs
  species_no <- length(probs)
  poisson_matrix <- matrix(0, species_no, repeats) # empty matrix to store the results of the random draws
  y <- rep(0, species_no) # empty vector for random draws
  for (j in 1:repeats) {
    for (k in 1:10) {
      y[k] <- if (probs_var[k]/sum(probs_var[1:species_no]) < 1)
        rpois(1,(probs_var[k]/sum(probs_var[1:species_no]))*pop_size) else pop_size 
    }
    #probs_var <- y/pop_size
    #print(sum(probs_var))
    poisson_matrix[,j] <- y
  }
  poisson_matrix
}

poisson_data_v2 <- poisson_population_v2(100000, 400, rep(0.1,10))
poisson_data_v2_slice <- poisson_data_v2[1,]
rel_poisson_v2 <- rep(0, 400)
for (i in 1:100000){
  rel_poisson_v2[poisson_data_v2_slice[i]] <- rel_poisson_v2[poisson_data_v2_slice[i]] + 1
}
rel_poisson_v2 <- rel_poisson_v2/100000

plot(1:400, rel_poisson_v2)
lines(dbinom(0:pop_size, pop_size, .1), col = "red")



# for poisson version 1, the population size varies but looks normally, when including the capacity
plot(1:100000,colSums(poisson_data))
# for poisson version 2, we again have something that looks normally distributed. 
plot(1:100000,(colSums(poisson_data_v2)))

rel_colSums_pois <- rep(0, 1000)
abs_colSums_pois <- colSums(poisson_data)
for (i in 1:100000){
  rel_colSums_pois[abs_colSums_pois[i]] <- rel_colSums_pois[abs_colSums_pois[i]] + 1
}
rel_colSums_pois <- rel_colSums_pois/100000
plot(1:1000,rel_colSums_pois)
lines(dnorm(0:1000, mean=400, sd=sqrt(400)), col = "red")
# also the shapiro wilk test suggests that the population size with version 2 is normally distributed
shapiro.test(rel_colSums_pois)



rel_colSums_pois2 <- rep(0, 1000)
abs_colSums_pois2 <- colSums(poisson_data_v2)
for (i in 1:100000){
  rel_colSums_pois2[abs_colSums_pois2[i]] <- rel_colSums_pois2[abs_colSums_pois2[i]] + 1
}
rel_colSums_pois2 <- rel_colSums_pois2/100000
plot(1:1000,rel_colSums_pois2)
lines(dnorm(0:1000, mean=400, sd=sqrt(400)), col = "red")
# also the shapiro wilk test suggests that the population size with version 2 is normally distributed
shapiro.test(rel_colSums_pois)






real_multi <- rep(0,10)
for (j in 1:100000) {
  real_multi <- real_multi + rmultinom(1,pop_size, probs)
}

plot(real_multi/100000)
y

plot(rmultinom(1,400,probs))

plot(0:400, dbinom(0:400, 400, .1))

binom_no <- rbinom(10000,400,.1)
rel_binom_no <- rep(0,400)
for (i in 1:10000){
  rel_binom_no[binom_no[i]] <- rel_binom_no[binom_no[i]] + 1
}
rel_binom_no <- rel_binom_no/10000
plot(1:400,rel_binom_no)
lines(dbinom(0:400, 400, .1), col = "red")


mbinom_no <- rmultinom(10000,400,probs)
mbinom_no_dim1 <-  mbinom_no[10,]
rel_mbinom_no <- rep(0,400)
for (i in 1:10000){
  rel_mbinom_no[mbinom_no_dim1[i]] <- rel_mbinom_no[mbinom_no_dim1[i]] + 1
}
rel_mbinom_no <- rel_mbinom_no/10000
plot(1:400,rel_mbinom_no)
lines(dbinom(0:400, 400, .1), col = "red")

