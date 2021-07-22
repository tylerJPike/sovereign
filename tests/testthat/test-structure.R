test_that("VAR workflow", {

  # simple time series
  AA = c(1:100) + rnorm(100)
  BB = c(1:100) + rnorm(100)
  CC = AA + BB + rnorm(100)
  inst = rbinom(100, size = 1, prob = .5)
  date = seq.Date(from = as.Date('2015-01-01'), by = 'month', length.out = 100)
  Data = data.frame(date, inst, AA, BB, CC)

  # estimate VAR without structure
  var =
    VAR(
      data = Data,
      p = 2,
      horizon = 10,
      freq = 'month',
      type = 'both',
      structure = NA)

  expect_true(is.list(var))
  expect_true(is.list(var$model))
  expect_true(is.list(var$forecasts))
  expect_true(is.list(var$residuals))

  var.irf = IRF(var, bootstrap.num = 10)

  expect_true(is.data.frame(var.irf))

  # estimate VAR with short-run restrictions
  var.short =
    VAR(
      data = Data,
      p = 2,
      horizon = 10,
      freq = 'month',
      type = 'both',
      structure = 'short')

  expect_true(is.list(var.short))
  expect_true(is.list(var.short$model))
  expect_true(is.list(var.short$forecasts))
  expect_true(is.list(var.short$residuals))

  var.irf = IRF(var.short, bootstrap.num = 10)
  var.hd = HD(var.short)
  var.fevd = FEVD(var.short)

  expect_true(is.data.frame(var.irf))
  expect_true(is.data.frame(var.hd))
  expect_true(is.data.frame(var.fevd))

  # estimate VAR with external instrument
  var.iv =
    VAR(
      data = Data,
      p = 2,
      horizon = 10,
      freq = 'month',
      type = 'both',
      structure = 'IV',
      instrument = 'inst',
      instrumented = 'AA')

  expect_true(is.list(var.iv))
  expect_true(is.list(var.iv$model))
  expect_true(is.list(var.iv$forecasts))
  expect_true(is.list(var.iv$residuals))

  var.irf = IRF(var.iv, bootstrap.num = 10)

  expect_true(is.data.frame(var.irf))


  # estimate VAR with external instrument and short-run restrictions
  var.iv_short =
    VAR(
      data = Data,
      p = 2,
      horizon = 10,
      freq = 'month',
      type = 'both',
      structure = 'IV-short',
      instrument = 'inst',
      instrumented = 'AA')

  expect_true(is.list(var.iv_short))
  expect_true(is.list(var.iv_short$model))
  expect_true(is.list(var.iv_short$forecasts))
  expect_true(is.list(var.iv_short$residuals))

  var.irf = IRF(var.iv_short, bootstrap.num = 10)
  var.hd = HD(var.iv_short)
  var.fevd = FEVD(var.iv_short)

  expect_true(is.data.frame(var.irf))
  expect_true(is.data.frame(var.hd))
  expect_true(is.data.frame(var.fevd))

  # estimate RVAR without structure

  Data = dplyr::mutate(Data, reg = dplyr::if_else(AA > median(AA), 1, 0))

  var =
    RVAR(
      data = Data,
      p = 2,
      regime = 'reg',
      horizon = 10,
      freq = 'month',
      type = 'both',
      structure = NA)

  expect_true(is.list(var))
  expect_true(is.list(var$model))
  expect_true(is.list(var$forecasts))
  expect_true(is.list(var$residuals))

  var.irf = IRF(var, bootstrap.num = 10)

  expect_true(is.list(var.irf))

  # estimate VAR with short-run restrictions
  var.short =
    RVAR(
      data = Data,
      p = 2,
      regime = 'reg',
      horizon = 10,
      freq = 'month',
      type = 'both',
      structure = 'short')

  expect_true(is.list(var.short))
  expect_true(is.list(var.short$model))
  expect_true(is.list(var.short$forecasts))
  expect_true(is.list(var.short$residuals))

  var.irf = IRF(var.short, bootstrap.num = 10)
  var.hd = HD(var.short)
  var.fevd = FEVD(var.short)

  expect_true(is.list(var.irf))
  expect_true(is.data.frame(var.hd))
  expect_true(is.list(var.fevd))

  # estimate VAR with external instrument
  var.iv =
    RVAR(
      data = Data,
      p = 2,
      regime = 'reg',
      horizon = 10,
      freq = 'month',
      type = 'both',
      structure = 'IV',
      instrument = 'inst',
      instrumented = 'AA')

  expect_true(is.list(var.iv))
  expect_true(is.list(var.iv$model))
  expect_true(is.list(var.iv$forecasts))
  expect_true(is.list(var.iv$residuals))

  var.irf = IRF(var.iv, bootstrap.num = 10)

  expect_true(is.list(var.irf))


  # estimate VAR with external instrument and short-run restrictions
  var.iv_short =
    RVAR(
      data = Data,
      p = 2,
      regime = 'reg',
      horizon = 10,
      freq = 'month',
      type = 'both',
      structure = 'IV-short',
      instrument = 'inst',
      instrumented = 'AA')

  expect_true(is.list(var.iv_short))
  expect_true(is.list(var.iv_short$model))
  expect_true(is.list(var.iv_short$forecasts))
  expect_true(is.list(var.iv_short$residuals))

  var.irf = IRF(var.iv_short, bootstrap.num = 10)
  var.hd = HD(var.iv_short)
  var.fevd = FEVD(var.iv_short)

  expect_true(is.list(var.irf))
  expect_true(is.data.frame(var.hd))
  expect_true(is.list(var.fevd))


})
