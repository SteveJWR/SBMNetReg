# create simple SBM
K = 2
Z = c(rep(1,5), rep(2,5))
G = matrix(1, length(Z),length(Z))
diag(G) = 0

test_that("Subgraph estimation works", {
  expect_equal(SBM_estimate_subgraph(Z,G)$P, matrix(c(1,1,1,1),2,2))
})


# TODO: Add ARD Estimator
# test_that("ARD estimation works", {
#   expect_equal(SBM_estimate_ARD(traits,X,Z)$P, matrix(c(1,1,1,1),2,2))
# })





