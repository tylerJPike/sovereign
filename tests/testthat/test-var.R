test_that("VAR workflow", {

  # simple time series
  AA = c(1:100) + rnorm(100)
  BB = c(1:100) + rnorm(100)
  CC = AA + BB + rnorm(100)
  date = seq.Date(from = as.Date('2000-01-01'), by = 'month', length.out = 100)
  Data = data.frame(date = date, AA, BB, CC)

  # estimate VAR
  var =
    VAR(
      data = Data,
      p = 1,
      horizon = 10,
      freq = 'month')

  expect_true(is.list(var))
  expect_true(is.list(var$model))
  expect_true(is.list(var$forecasts))
  expect_true(is.list(var$residuals))

  # estimate IRF
  irf =
    var_irf(
      var,
      bootstraps.num = 10,
      CI = c(0.05,0.95))

  expect_true(is.list(irf))
  expect_true(is.data.frame(irf$ci.lower))
  expect_true(is.data.frame(irf$ci.upper))

  # estimate forecast error variance decomposition
  fevd =
    var_fevd(
      var,
      horizon = 10)

  expect_true(is.data.frame(fevd))

})
