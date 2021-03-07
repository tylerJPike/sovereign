test_that("threshold VAR workflow", {

  # simple time series
  AA = c(1:100) + rnorm(100)
  BB = c(1:100) + rnorm(100)
  CC = AA + BB + rnorm(100)
  date = seq.Date(from = as.Date('2000-01-01'), by = 'month', length.out = 100)
  Data = data.frame(date = date, AA, BB, CC)
  Data = dplyr::mutate(Data, reg = dplyr::if_else(AA > median(AA), 1, 0))

  # estimate VAR
  tvar =
    threshold_VAR(
      data = Data,
      regime = 'reg',
      p = 1,
      horizon = 10,
      freq = 'month')

  expect_true(is.list(tvar))
  expect_true(is.list(tvar$model))
  expect_true(is.list(tvar$forecasts))
  expect_true(is.list(tvar$residuals))

  # estimate IRF
  irf =
    threshold_var_irf(
      tvar,
      bootstraps.num = 10,
      CI = c(0.05,0.95))

  expect_true(is.data.frame(irf[[1]]))

  # estimate forecast error variance decomposition
  fevd =
    threshold_var_fevd(
      tvar,
      horizon = 10)

  expect_true(is.list(fevd))

})
