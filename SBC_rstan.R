library(rethinking)
library(phytools)
library(phangorn)
source("chop_tree.R")

## Simulate a random phylogenetic tree (pure birth)
N_tips <- 25

# Vectors to hold rank statistics for A matrix
A_1_2_rank <- c()
A_2_1_rank <- c()


# Compile stan models
sim_mod <- stan_model(file="DPM_sim.stan")
fit_mod <- stan_model(file="DPM_fit.stan")
###################################################
#### Loop over simulations ########################
N_sims <- 1000

for ( s in 1:N_sims ) {
  
  tree <- pbtree(n = N_tips)
  
  tree_data <- chop_tree(tree)
  
  #  Specify number of traits
  J <- 2
  
  tree_data$J <- 2
  tree_data$resp_type <- c(1,1)

  # Draw samples from generative model
  sim <- sampling(sim_mod, data=tree_data, iter=1, chains=1, algorithm="Fixed_param")
  
  post_sim <- extract.samples(sim)
  
  A_sim <- post_sim$A[1,,]
  
  # Add simulated outcome data
  tree_data$X <- post_sim$X[1,,]
  
  n_warmup <- 100
  n_samps <- 350
  n_chains <- 8
  
  fit <- sampling(fit_mod, data=tree_data, warmup=n_warmup, iter=n_samps, chains=8, cores=8, control=list(adapt_delta=0.9), init=0.1)
  
  post_fit <- extract.samples(fit)
  A_fit <- post_fit$A
  
  A_1_2_rank[s] <- rank( c(A_sim[1,2], A_fit[,1,2]) )[1]
  A_2_1_rank[s] <- rank( c(A_sim[2,1], A_fit[,2,1]) )[1]

} # end loop over simulations

# Save results
write.csv(data.frame(A_1_2_rank, A_2_1_rank), file="SBC_ranks.csv")

#ranks <- read.csv("SBC_ranks.csv")
#off_diag_ranks <- cbind(ranks$A_1_2_rank, ranks$A_2_1_rank)
#hist(off_diag_ranks/2000, breaks=100, main="", xlab="Rank (prop max)", col="skyblue", border=NA)

#hist(off_diag_ranks[,1]/2001, breaks=30, main="", xlab="Rank (prop max)", col="skyblue", border=NA)

