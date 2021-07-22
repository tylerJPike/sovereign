test_that("VAR workflow", {

  # simple time series
  AA = c(1:100) + rnorm(100)
  BB = c(1:100) + rnorm(100)
  CC = AA + BB + rnorm(100)
  date = seq.Date(from = as.Date('2015-01-01'), by = 'month', length.out = 100)
  Data = data.frame(date = date, AA, BB, CC)


  # estimate VAR (without lag selection)
  var =
    VAR(
      data = Data,
      p = 2,
      horizon = 10,
      freq = 'month',
      type = 'both')

  expect_true(is.list(var))
  expect_true(is.list(var$model))
  expect_true(is.list(var$forecasts))
  expect_true(is.list(var$residuals))

  # estimate VAR (with lag selection)
  var =
    VAR(
      data = Data,
      horizon = 10,
      freq = 'month',
      lag.ic = 'BIC',
      lag.max = 3)

  expect_true(is.list(var))
  expect_true(is.list(var$model))
  expect_true(is.list(var$forecasts))
  expect_true(is.list(var$residuals))

  # plot forecasts
  plot.forecast =
    plot_forecast(var$forecasts[[1]])

  expect_true(is.list(plot.forecast))

  plot.errors =
    plot_error(var$residuals[[1]])

  expect_true(is.list(plot.errors))

  # estimate IRF
  irf =
    var_irf(
      var,
      bootstrap.num = 10,
      CI = c(0.05,0.95))

  expect_true(is.data.frame(irf))

  # plot IRF
  plot.irf = plot_irf(irf)

  expect_true(is.list(plot.irf))

  # estimate forecast error variance decomposition
  fevd =
    var_fevd(
      var,
      horizon = 10,
      scale = TRUE)

  expect_true(is.data.frame(fevd))

  # plot fevd
  plot.fevd = plot_fevd(fevd)

  expect_true(is.list(plot.fevd))

  # estimate hd
  hd =
    var_hd(var)

  expect_true(is.data.frame(hd))

  # plot hd
  plot.hd = plot_hd(hd)

  expect_true(is.list(plot.hd))

  # test covid correction
  var.corrected = covid_volatility_correction(var)

  expect_true(is.list(var.corrected))
  expect_true(is.list(var.corrected$model))
  expect_true(is.list(var.corrected$forecasts))
  expect_true(is.list(var.corrected$residuals))
  expect_true(is.vector(var.corrected$correction.factors))

})
