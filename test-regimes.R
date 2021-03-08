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
  learn_regimes(
    data = Data,
    engine = 'rf'
  )

 regime.kmeans =
   learn_regimes(
     data = Data,
     engine = 'kmeans'
   )

 regime.em =
   learn_regimes(
     data = Data,
     engine = 'EM'
   )

 expect_true(is.data.frame(regimes.rf))
 expect_true(is.data.frame(regime.kmeans))
 expect_true(is.data.frame(regime.em))

})
