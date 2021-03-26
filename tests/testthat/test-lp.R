test_that("Local projection workflow", {

  # simple time series
  AA = c(1:100) + rnorm(100)
  BB = c(1:100) + rnorm(100)
  CC = AA + BB + rnorm(100)
  date = seq.Date(from = as.Date('2000-01-01'), by = 'month', length.out = 100)
  Data = data.frame(date = date, AA, BB, CC)
  Data = dplyr::mutate(Data, reg = dplyr::if_else(AA > median(AA), 1, 0))

  # local projection forecasts
  lp =
    LP(
      data = Data,
      horizon = 1,
      lag.ic = 'AIC',
      lag.max = 4,
      type =  'none',
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
      type = 'both',
      NW = TRUE,
      NW_prewhite = FALSE,
      NW_lags = 1,
      freq = 'month')

  expect_true(is.list(lp))
  expect_true(is.list(lp$model))
  expect_true(is.list(lp$forecasts))
  expect_true(is.list(lp$residuals))

  # estimate IRF
  irf = lp_irf(lp)
  expect_true(is.data.frame(irf))

  # plot IRF
  plot.all = plot_irf(irf)
  expect_true(is.list(plot.all))

  # multi-regime local projection
  rlp =
    RLP(
      data = Data,
      regime = 'reg',
      p = 1,
      horizon = c(1:10),
      NW = FALSE,
      freq = 'month')


  expect_true(is.list(rlp))

  rlp =
    RLP(
      data = Data,
      p = 1,
      horizon = c(1:10),
      NW = FALSE,
      freq = 'month')


  expect_true(is.list(rlp))


  rlp =
    RLP(
      data = Data,
      regime = 'reg',
      type = 'both',
      horizon = c(1:10),
      lag.ic = 'AIC',
      lag.max = 4,
      freq = 'month')

  expect_true(is.list(rlp))

  # multi-regime horizon irf
  rlp.irf = rlp_irf(rlp)

  expect_true(is.list(rlp.irf))

})
