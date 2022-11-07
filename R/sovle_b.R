
#------------------------------------------
# Function to estimate the structural
#  error impact matrix
#
# See Cesa-Bianchi's VAR-Toolbox for
# MATLAB implementation
#
# IV-short implementation based on matlab
# code by Nick von Turkovich and Paco
# Vazquez-Grande
#------------------------------------------
solve_B = function(var, report_iv = TRUE){

  # function variables
  forecast.date = NA

  # set residuals
  if(class(var) == 'VAR'){
    residuals = var$residuals[[1]]
  }else{
    residuals = var$residuals
  }

  if(is.na(var$structure) == TRUE){

    # retrieve reduced-form residuals
    data.residuals = residuals

    # reduced form variance-covariance matrix
    cov.matrix = stats::var(stats::na.omit(dplyr::select(data.residuals, -date, -forecast.date)))

    B = cov.matrix

    return(B)

  }else if(var$structure == 'short'){

    # retrieve reduced-form residuals
    data.residuals = residuals

    # reduced form variance-covariance matrix
    cov.matrix = stats::var(stats::na.omit(dplyr::select(data.residuals, -date, -forecast.date)))

    # take cholesky decomposition
    B = t(chol(cov.matrix))

    return(B)

  }else if(var$structure == 'IV'){

    # retrieve reduced-form residuals
    data.residuals = residuals
    covariates = colnames(dplyr::select(data.residuals, -date, -forecast.date))
    col_to_instrument = var$instrumented

    data.residuals[,covariates] = -1*data.residuals[,covariates]

    # retrieve instrument
    data.instrument = var$instrument
    instrument = colnames(dplyr::select(data.instrument, -date))

    # combine data
    data =
      dplyr::inner_join(data.residuals, data.instrument, by = 'date') %>%
      stats::na.omit()

    # first stage least squares
    model.first_stage =
      stats::lm(col_to_instrument ~.,
         data = data %>%
           dplyr::select(col_to_instrument = dplyr::all_of(col_to_instrument), dplyr::all_of(instrument)) %>%
           stats::na.omit())
    p_hat = model.first_stage$fitted.values

    if(report_iv == TRUE){
      print(summary(model.first_stage))
    }

    # second stage least squares
    # to estimate first column of B
    # (automatically scales first entry to 1)
    instrumented_shocks = covariates %>%
      purrr::map_df(.f = function(X){

        model.second_stage = stats::lm(data[,X] ~ p_hat)
        second_stage_beta = model.second_stage$coefficients[2]
        return(second_stage_beta)

      })

    # scale size of the shock
    #  see Gertler and Karadi (2015) for background
    #  see Cesa-Bianchi's VAR-Toolbox for MATLAB implementation
    X.demean =
      dplyr::select(data, -date, -forecast.date, -dplyr::all_of(instrument)) %>%
      dplyr::mutate_all(function(X){return(X - mean(X, na.rm = T))}) %>%
      as.matrix()

    sigma_b = 1/(nrow(data) - ncol(var$model$coef) + 1) * t(X.demean) %*% X.demean

    s21s11 = instrumented_shocks[2:length(covariates),] %>% as.matrix()
    S11 = sigma_b[1,1] %>% as.numeric()
    S21 = sigma_b[2:nrow(sigma_b),1] %>%as.vector()
    S22 = sigma_b[2:nrow(sigma_b),2:ncol(sigma_b)]

    Q = (s21s11 * S11)  %*% t(s21s11) - (S21 %*% t(s21s11) + as.matrix(s21s11) %*% t(S21)) + S22
    sp = sqrt( S11 - t(S21 - as.matrix(s21s11) %*% S11) %*% as.matrix(solve(Q) %*%  (S21 - s21s11 * S11)) )

    scaled_instrumented_shocks = instrumented_shocks * as.numeric(sp)

    # prepare B matrix
    B = matrix(0, ncol = (length(covariates) - 1), nrow = length(covariates))
    B = cbind(scaled_instrumented_shocks, B)

    return(B)

  }else if(var$structure == 'IV-short'){

    # align instrument and residuals
    valid_dates =
      dplyr::inner_join(
        dplyr::select(stats::na.omit(residuals), date),
        dplyr::select(stats::na.omit(var$instrument), date),
        by = 'date'
      )

    # set residuals matrix
    residuals = residuals %>%
      dplyr::inner_join(valid_dates, by = 'date') %>%
      dplyr::select(-date, -forecast.date) %>%
      as.matrix()

    residuals = -1*residuals

    # set instrument
    instrument = var$instrument %>%
      data.frame() %>%
      dplyr::inner_join(valid_dates, by = 'date') %>%
      dplyr::select(-date)

    instrument.mean = instrument %>%
      dplyr::mutate_all(mean, na.rm = T)

    # number of observations
    n.obs = nrow(residuals)

    # solve for IV implied impact
    pi = solve(n.obs^(-1) * t(residuals) %*% residuals) %*%
      (n.obs^(-1) * t(residuals) %*% as.matrix(instrument - instrument.mean))

    phi_sq =
      (n.obs^(-1) * t(instrument - instrument.mean) %*% residuals) %*%
      solve( n.obs^(-1) * t(residuals) %*% residuals ) %*%
      ( n.obs^(-1) * t(residuals) %*% as.matrix(instrument - instrument.mean) )

    B1 =
      n.obs^(-1) *
      ( t(residuals) %*% as.matrix(instrument - instrument.mean) ) %*%
      (phi_sq)^(-0.5)

    B = matrix(ncol = ncol(residuals), nrow = ncol(residuals))
    B[,1:ncol(instrument)] = B1
    rownames(B) = colnames(B) = colnames(residuals)

    # solve for orthogonalized structural shock
    model.first_stage = stats::lm(instrument[,1] ~ residuals)
    orthogonal_instrument = instrument - model.first_stage$residuals
    orthogonal_instrument = orthogonal_instrument[,] / stats::sd(orthogonal_instrument[,], na.rm = T)

    shock_matrix = matrix(nrow = nrow(residuals), ncol = ncol(residuals))
    shock_matrix[,1] = orthogonal_instrument

    # reduced form variance-covariance matrix
    sigma = stats::var(residuals)

    # solve additional entries
    # with a lower triangular restriction
    order_sequence = c(1:ncol(residuals))

    for(i in order_sequence){

      Y = residuals[,i]
      X = shock_matrix[,c(1:i)]
      model.second_stage = stats::lm(Y~X)

      if(i != utils::tail(order_sequence,1)){

        B[i,] =
          c( model.second_stage$coef[-1],
             stats::sd(model.second_stage$residuals),
             rep(0, length(order_sequence) - ncol(instrument) - i) )

        shock_matrix[,i+1] = model.second_stage$residuals / stats::sd(model.second_stage$residuals)

      }else{

        B[i,] =  model.second_stage$coef[-1]

      }

    }

    # make diagonal entries positive
    shock_ordering = data.frame(residuals) %>%
      dplyr::select(var$instrumented, dplyr::everything()) %>%
      colnames()
    B.sign = diag(sign(diag(B[shock_ordering,])))
    B = B %*% B.sign

    return(B)

  }else{

    stop('structure must be set to either NA, "short", "IV", or "IV-short".')

  }

}
