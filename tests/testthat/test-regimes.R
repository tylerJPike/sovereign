test_that("Regime learning functions", {

  # simple time series
  AA = c(1:100) + rnorm(100)
  BB = c(1:100) + rnorm(100)
  CC = AA + BB + rnorm(100)
  date = seq.Date(from = as.Date('2000-01-01'), by = 'month', length.out = 100)
  Data = data.frame(date = date, AA, BB, CC)
  Data = dplyr::mutate(Data, reg = dplyr::if_else(AA > median(AA), 1, 0))

 # run ml
 regimes.rf =
  regimes(
    data = Data,
    method = 'rf'
  )

 regime.kmeans =
   regimes(
     data = Data,
     method = 'kmeans'
   )

 regime.em =
   regimes(
     data = Data,
     method = 'EM'
   )

 regime.bp =
   regimes(
     data = Data[, c(1,2)],
     method = 'BP'
   )

 expect_true(is.data.frame(regimes.rf))
 expect_true(is.data.frame(regime.kmeans))
 expect_true(is.data.frame(regime.em))
 expect_true(is.data.frame(regime.bp))

})
