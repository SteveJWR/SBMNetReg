test_that("network features work", {
  G = matrix(c(0,1,1,1,0,0,1,0,0), nrow=3, ncol=3)
  data = data.frame( 'X' = c(1,2,3), 'A' = c(0,1,1))
  deg = degree(data = NULL, G = G)
  t_neigh = treated_neighbors(data, G, treatment_column = 'A')
  frac_t_neigh = frac_treated_neighbors(data, G, treatment_column = 'A')
  power2_treated_neighbors = power_treated_neighbors(data, G, treatment_column = 'A', power = 2)

  expect_equal(deg, c(2,1,1))
  expect_equal(t_neigh, c(2,0,0))
  expect_equal(frac_t_neigh, c(1,0,0))
  expect_equal(power2_treated_neighbors, c(0,2,2))
})
