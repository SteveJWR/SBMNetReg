

#### Create the graph model

set.seed(1)
K = 5
n = 100
P = matrix(0.05, K, K)
diag(P) = 0.3
Z = rep(seq(5), each = n)
traits <- Z

SBM_draw <- generateSBM(P = P, Z = Z)
G = SBM_draw$G
ARD <- computeARD(traits, G)


data = data.frame('X' = rnorm(n*K), 'A' = rbinom(K*n, 1, c(rep(0.1, 2*n), rep(0.9, 3*n))))
data_wide <- network_features(data, G, feat_list = list('degree'= degree, 'frac.treated' =  frac_treated_neighbors))



# Simulate the linear model
coef = c(1,2,3,4)
X_fmla = formula(~ X +  A + frac.treated)
sim_lm <- simLinearModel(data_wide,X_fmla, coef, sd = 1)

data_obs <- data
data_obs$Y = sim_lm$y


# Estimate the SBM
traits <- Z
ARD_norm = ARD/rowSums(ARD)

Z.hat = clusterARDHclust(ARD_norm, K,alg = 'complete')
SBM = SBM_estimate_ARD(traits, ARD, K)

fmla = formula(Y ~ X +  A + frac.treated)

SBM_model = SBM_lm(fmla, data_obs, SBM, network_feature_list = list('degree'= degree, 'frac.treated' =  frac_treated_neighbors), B = 1000)






# Simulate a Logistic Model
# fraction of treated is unknown, we need to simulate
coef = c(-2,2,3,4) # negative intercept here
X_fmla = formula(~ X +  A + frac.treated)
sim_log <- simLogisticModel(data_wide,X_fmla, coef)


data_obs <- data
data_obs$Y = sim_log$y


# Estimate the SBM
traits <- Z
ARD_norm = ARD/rowSums(ARD)

Z.hat = clusterARDHclust(ARD_norm, K,alg = 'complete')
SBM = SBM_estimate_ARD(traits, ARD, K)

data_obs$Y = sim_log$y

fmla = formula(Y ~ X +  A + frac.treated)

SBM_glm_model = SBM_glm(fmla,
                        data_obs,
                        SBM = SBM,
                        family=quasibinomial(),
                        network_feature_list = list(degree = degree,
                                                    frac.treated = frac_treated_neighbors),
                        B = 200
                        )

SBM_glm_model$model






