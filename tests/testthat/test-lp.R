test_that("Local projection workflow", {

  # simple time series
  AA = c(1:100) + rnorm(100)
  BB = c(1:100) + rnorm(100)
  CC = AA + BB + rnorm(100)
  date = seq.Date(from = as.Date('2000-01-01'), by = 'month', length.out = 100)
  Data = data.frame(date = date, AA, BB, CC)
  Data = dplyr::mutate(Data, reg = dplyr::if_else(AA > median(AA), 1, 0))

  # local proejctions
  irfs =
    lp_irf(
      data = Data,
      shock = 'AA',
      target = 'BB',
      horizons = 20,
      lags = 2)

  t.irfs =
    threshold_lp_irf(
      data = dplyr::select(Data, -reg),
      thresholdVar = dplyr::select(Data, reg, date),
      shock = 'AA',
      target = 'BB',
      horizons = 20,
      lags = 2)

  expect_true(is.data.frame(irfs))
  expect_true(is.list(t.irfs))


})
