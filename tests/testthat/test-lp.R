test_that("Local projection workflow", {

  # simple time series
  AA = c(1:100) + rnorm(100)
  BB = c(1:100) + rnorm(100)
  CC = AA + BB + rnorm(100)
  date = seq.Date(from = as.Date('2000-01-01'), by = 'month', length.out = 100)
  Data = data.frame(date = date, AA, BB, CC)
  Data = dplyr::mutate(Data, reg = dplyr::if_else(AA > median(AA), 1, 0))

  # local projection forecasts (new workflow)
  lp =
    LP(
      data = Data,
      p = 1,
      horizon = 1,
      freq = 'month')

  expect_true(is.list(lp))
  expect_true(is.list(lp$model))
  expect_true(is.data.frame(lp$forecasts))
  expect_true(is.data.frame(lp$residuals))

  lp =
    LP(
      data = Data,
      p = 1,
      horizon = c(1:10),
      NW = TRUE,
      freq = 'month')

  expect_true(is.list(lp))
  expect_true(is.list(lp$model))
  expect_true(is.list(lp$forecasts))
  expect_true(is.list(lp$residuals))

  # estimate IRF
  irf = lp_irf(lp)
  expect_true(is.data.frame(irf))

  # plot IRF
  plot.all = irf_plot(irf)
  expect_true(is.list(plot.all))

  # multi-regime local projection
  tlp =
    threshold_LP(
      data = Data,
      regime = 'reg',
      p = 1,
      horizon = c(1:10),
      NW = FALSE,
      freq = 'month')

  expect_true(is.list(tlp))

  #tirf = threshold_lp_irf(tlp)
  #expect_true(is.list(irf))

})
