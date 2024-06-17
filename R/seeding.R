


#' Optimal Seeding for Contagion
#'
#' @param SBM Model Parameters for the SBM
#' @param n_seeds Number of total seeds
#' @param num_sim Number of Simulations to conduct for each seed combination
#' @param contagion_type Simple or Complex
#' @param complex_mean Complex Contagion Thresholds mean
#' @param complex_sd Complex Contagion Thresholds sd
#' @param ... Additional arguments
#'
#' @return Dataframe with results of the optimal seeding
#' @export
ContagionOptimalSeeding <- function(SBM,n_seeds = 2, num_sim = 2000, contagion_type = 'simple', complex_mean = NULL, complex_sd = NULL, ...){
  K = nrow(SBM$P) # number of blocks in the SBM
  n = length(SBM$Z) # number of nodes

  grid <- expand.grid(rep(list(1:K), n_seeds))
  n_clusters = as.numeric(table(SBM$Z))

  n_comb = nrow(grid)
  avg_vec = c()
  sd_vec = c()
  comb_list = c()

  results_row = 1
  for(i in seq(n_comb)){
    cat(paste('Seed Combination', i, '/', n_comb), end = '\r') # print progress
    cluster_seeds = as.numeric(grid[i,])
    if(identical(cluster_seeds, sort(cluster_seeds))){
      treat_counts = numeric(K)
      for(j in seq(length(cluster_seeds))){
        treat_counts[cluster_seeds[j]] = treat_counts[cluster_seeds[j]] + 1
      }
      treat_levels = treat_counts / n_clusters

      total_adopted = numeric(num_sim)
      for(sim in seq(num_sim)){

        g.sim <- generateSBM(P = SBM$P,Z = SBM$Z)
        G = g.sim$G

        # treatment levels
        A = saturationRandomizationTreatment(SBM$Z, treat_levels)
        if(contagion_type == 'simple'){
          adopted = simSimpleContagion(A, G, ...)
        } else if(contagion_type == 'complex'){
          thresholds = BeamanAdoptionThreshold(n, mean = complex_mean, sd = complex_sd)
          adopted = simComplexContagion(A, G, thresholds, ...)
        }

        total_adopted[sim] = sum(adopted)

      }
      avg_adopted = mean(total_adopted)
      se_adopted = sd(total_adopted) / sqrt(length(total_adopted))

      avg_vec[results_row] = avg_adopted
      sd_vec[results_row] = se_adopted
      comb_list[results_row] = paste(cluster_seeds, collapse = ", ")

      results_row = results_row + 1
    }
  }
  results = data.frame('avg_adopted' = avg_vec, 'se_adopted' = sd_vec, 'cluster_seeds' = comb_list, 'rank' = rank(-avg_vec))
  return(results)
}


ContagionOptimalSeedingFullGraph <- function(G,n_seeds = 2, num_sim = 2000, contagion_type = 'simple', complex_mean = NULL, complex_sd = NULL, ...){

  n = length(SBM$G) # number of nodes

  grid <- expand.grid(rep(list(1:n), n_seeds))
  n_clusters = as.numeric(table(SBM$Z))

  n_comb = nrow(grid)
  avg_vec = c()
  sd_vec = c()
  comb_list = c()

  results_row = 1
  for(i in seq(n_comb)){
    cat(paste('Seed Combination', i, '/', n_comb), end = '\r') # print progress
    cluster_seeds = as.numeric(grid[i,])

    # treatment seeds
    A = rep(0,n)
    A[cluster_seeds] = 1
    for(sim in seq(num_sim)){
      if(contagion_type == 'simple'){
        adopted = simSimpleContagion(A, G, ...)
      } else if(contagion_type == 'complex'){
        thresholds = BeamanAdoptionThreshold(n, mean = complex_mean, sd = complex_sd)
        adopted = simComplexContagion(A, G, thresholds, ...)
      }

      total_adopted[sim] = sum(adopted)

    }
    avg_adopted = mean(total_adopted)
    se_adopted = sd(total_adopted) / sqrt(length(total_adopted))

    avg_vec[results_row] = avg_adopted
    sd_vec[results_row] = se_adopted
    comb_list[results_row] = paste(cluster_seeds, collapse = ", ")

    results_row = results_row + 1
  }
  results = data.frame('avg_adopted' = avg_vec, 'se_adopted' = sd_vec, 'cluster_seeds' = comb_list, 'rank' = rank(-avg_vec))
  return(results)
}



