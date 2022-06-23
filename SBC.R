library(cmdstanr)
library(phytools)
library(phangorn)
library(posterior)
source("chop_tree.R")

## Simulate a random phylogenetic tree (pure birth)
N_tips <- 25

# Vectors to hold rank statistics for A matrix
A_1_2_rank <- c()
A_2_1_rank <- c()

# Compile stan models
sim_mod <- cmdstan_model(stan_file="DPM_sim.stan")
fit_mod <- cmdstan_model(stan_file="DPM_fit.stan")
###################################################
#### Loop over simulations ########################
N_sims <- 250

for ( s in 1:N_sims ) {
  
  tree <- pbtree(n = N_tips)
  
  tree_data <- chop_tree(tree)
  
  #  Specify number of traits
  J <- 2
  
  tree_data$J <- 2
  tree_data$resp_type <- c(1,1)

# Draw samples from generative model
sim <- sim_mod$sample(
  data=tree_data,
  iter_warmup = 0,
  iter_sampling = 1,
  fixed_param = T,
  chains = 1
)

post_sim <- posterior::as_draws_rvars(sim$draws())

A_sim <- draws_of(post_sim$A)[1,,]

# Add simulated outcome data
tree_data$X <- draws_of(post_sim$X)[1,,]

n_samps <- 250
n_chains <- 8

fit <- fit_mod$sample(
  data=tree_data,
  iter_warmup = n_chains,
  iter_sampling = 250,
  init = 0,
  chains = n_chains,
  parallel_chains = 8,
  adapt_delta = 0.9
)

post_fit <- posterior::as_draws_rvars(fit$draws())
A_fit <- draws_of(post_fit$A)

A_1_2_rank[s] <- rank( c(A_sim[1,2], A_fit[,1,2]) )[1]
A_2_1_rank[s] <- rank( c(A_sim[2,1], A_fit[,2,1]) )[1]

} # end loop over simulations

# Save results
write.csv(data.frame(A_1_2_rank, A_2_1_rank), file="SBC_ranks.csv")
