# Sovereign: State-Dependent Empirical Analysis  

## Version 1.2.0 
(2021-07-23)

- Major updates
    - VARs may now be estimated as Proxy-SVARs
    - IRF functions now support
      - Choice of bootstrapping techniques 
        - standard residual resampling
        - Wild bootstrapping
      - Choice of structural shocks 
        - no structural restrictions, 
        - short-term restrictions (via Cholesky decomposition), 
        - instrumental variable estimation  
        - combination of instrument variable and short-term restrictions  
      - Parallel processing for greater computational efficiency in computing bootstrapped confidence intervals   
      
- Minor updates
  - IRF plots use darker blue for CI
  - Analysis function wrappers implemented for IRFs, FEVDs, and HDs
  - Updated model documentation with academic references and corrected minor errors 
   
- Bug fixes
    - plot_irf() plots all targets by default 
    - test-threvhold_var.R changed to test-regime_var.R
    - Fixed IRF results format listed in documentation

## Version 1.1.0   
(2021-06-01)

- Major updates
    - Include COVID shock correction a la Lenza and Primiceri (2020)
    - Include historical decomposition of VAR shocks
      
- Minor updates
    - Create object classes for LP and VAR output  
    - n.lag no longer requires a 'date' column 
    - Specify structural assumptions in IRFs  
    - Clarify model specification of RVAR  

- Bug fixes
    - Convert errorCondition() to stop()   

## Version 1.0.0 
(2021-04-02)

- Initial release 
