
# A simple application of the diffusion model for a more efficient simulation
optimal_single_diffusion_seed <- function(SBM_list, q = 0.2){
  SBM_tmp = SBM_list[[1]]
  K = nrow(SBM_tmp$P)
  B = length(SBM_list)

  # List of variance of the estimates
  average_variance = numeric(K)
  sd_variance = numeric(K)

  for(k_treat in seq(K)){

    variance_boot_vec = numeric(B)
    for(b in seq(B)){
      SBM_b = SBM_list[[b]]
      P = SBM_b$P
      Z = SBM_b$Z
      X_reg = P[k_treat, Z]

      obs_var = q*X_reg*(1 - q*X_reg)
      variance = sum(X_reg**2)**(-2)*sum(X_reg**2*obs_var)
      variance_boot_vec[b] = variance
    }
    average_variance[k_treat] = mean(variance_boot_vec)
    sd_variance[k_treat] = sd(variance_boot_vec)
  }
  opt_allocation = which.min(average_variance)
  return(list('average_variance' = average_variance, 'sd_variance' = sd_variance, 'opt_allocation' = opt_allocation))
}



# example for diffusion estimation over graphs

set.seed(1)

# Generate ARD
P = matrix(0.05, 5, 5)
diag(P) = seq(0.1,0.5, length.out = 5)
Z = rep(seq(5), each = 100)
traits <- Z # cluster membership are traits directly


SBM_draw <- generateSBM(P = P, Z = Z)
G = SBM_draw$G
n = nrow(G)
X <- computeARD(traits, G)

# Diffusion simulation:

q_oracle = 0.05 #

q_working_est = 0.2
SBM_est = SBM_estimate_ARD(traits, X, K = 5)
SBM_list <- SBM_uncertainty_ARD(traits, X, K = 5, L = 1000)

optimal_diffusion = optimal_single_diffusion_seed(SBM_list)

opt_K = optimal_diffusion$opt_allocation



# Simulation using optimal assignments and random assignments
n_sims = 100

random_assignment_q_est = numeric(n_sims)
random_assignment_q_est_partial_data = numeric(n_sims)
degree_opt_assignment_q_est = numeric(n_sims)
degree_opt_assignment_q_est_partial_data = numeric(n_sims)
model_opt_assignment_q_est = numeric(n_sims)
model_opt_assignment_q_est_partial_data = numeric(n_sims)


# indicates no intercept
fmla <- formula(Y ~ 0 + treated_neighbor)


for(sim in seq(n_sims)){
  cat(paste('Sim:', sim, '/', n_sims), end='\r')
  A.random = numeric(n)
  A.degree.opt = numeric(n)
  A.model.opt = numeric(n)

  A.random[sample(1:n, size = 1)] = 1
  max_degree = max(Matrix::colSums(G))
  A.degree.opt[sample(which(Matrix::colSums(G) == max_degree), size = 1)] = 1
  A.model.opt[sample(which(SBM_est$Z == opt_K), size = 1)] = 1

  # Simulate data using the true graph
  Y.random <- simSimpleContagion(A.random,G, T_steps = 1, q = q_oracle)
  Y.degree.opt <- simSimpleContagion(A.degree.opt,G, T_steps = 1, q = q_oracle)
  Y.model.opt <- simSimpleContagion(A.model.opt,G, T_steps = 1, q = q_oracle)


  # full data regressions:
  data.full.random = data.frame(Y = Y.random, A = A.random)
  data.full.degree.opt = data.frame(Y = Y.degree.opt, A = A.degree.opt)
  data.full.model.opt = data.frame(Y = Y.model.opt, A = A.model.opt)

  data.full.random$treated_neighbor <- treated_neighbors(data.full.random, G)
  data.full.degree.opt$treated_neighbor <- treated_neighbors(data.full.degree.opt, G)
  data.full.model.opt$treated_neighbor <- treated_neighbors(data.full.model.opt, G)

  q.est.full.random = lm(fmla, data = data.full.random)$coefficients
  q.est.full.degree.opt = lm(fmla, data = data.full.degree.opt)$coefficients
  q.est.full.model.opt = lm(fmla, data = data.full.model.opt)$coefficients


  random_assignment_q_est[sim] = q.est.full.random
  degree_opt_assignment_q_est[sim] = q.est.full.degree.opt
  model_opt_assignment_q_est[sim] = q.est.full.model.opt


  # partial data regressions:
  data.partial.random = data.frame(Y = Y.random, A = A.random)
  data.partial.degree.opt = data.frame(Y = Y.degree.opt, A = A.degree.opt)
  data.partial.model.opt = data.frame(Y = Y.model.opt, A = A.model.opt)

  sbm_lm_random <- SBM_lm(fmla, data.partial.random, SBM_est, network_feature_list = list('treated_neighbor'= treated_neighbors), B = 200, verbose = F)
  sbm_lm_degree_opt <- SBM_lm(fmla, data.partial.degree.opt, SBM_est, network_feature_list = list('treated_neighbor'= treated_neighbors), B = 200, verbose = F)
  sbm_lm_model_opt <- SBM_lm(fmla, data.partial.model.opt, SBM_est, network_feature_list = list('treated_neighbor'= treated_neighbors), B = 200, verbose = F)

  q.est.partial.random = sbm_lm_random$model$coefficients
  q.est.degree.opt = sbm_lm_degree_opt$model$coefficients
  q.est.model.opt = sbm_lm_model_opt$model$coefficients

  random_assignment_q_est_partial_data[sim] = q.est.partial.random
  degree_opt_assignment_q_est_partial_data[sim] = q.est.degree.opt
  model_opt_assignment_q_est_partial_data[sim] = q.est.model.opt

}

random_assignment_q_est
degree_opt_assignment_q_est
model_opt_assignment_q_est
random_assignment_q_est_partial_data
degree_opt_assignment_q_est_partial_data
model_opt_assignment_q_est_partial_data
rmse_vec = c(sqrt(mean((random_assignment_q_est - q_oracle)^2)),
             sqrt(mean((degree_opt_assignment_q_est - q_oracle)^2)),
             sqrt(mean((model_opt_assignment_q_est - q_oracle)^2)),
             sqrt(mean((random_assignment_q_est_partial_data - q_oracle)^2)),
             sqrt(mean((degree_opt_assignment_q_est_partial_data - q_oracle)^2)),
             sqrt(mean((model_opt_assignment_q_est_partial_data - q_oracle)^2))
             )
bias_vec = c(mean((random_assignment_q_est - q_oracle)),
             mean((degree_opt_assignment_q_est - q_oracle)),
             mean((model_opt_assignment_q_est - q_oracle)),
             mean((random_assignment_q_est_partial_data - q_oracle)),
             mean((degree_opt_assignment_q_est_partial_data - q_oracle)),
             mean((model_opt_assignment_q_est_partial_data - q_oracle))
)


