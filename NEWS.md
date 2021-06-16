# Sovereign: State-Dependent Empirical Analysis  

## Version 1.2.0

- Major updates
    - VARs may now be estimated as Proxy-SVARs (with multiple instruments)
    - IRF functions now support
      - The choice between standard and Wild bootstrapping techniques
      - The choice between no structural restrictions, short-term restrictions (via Cholesky decomposition), or instrumental variable restrictions
      - The choice run bootstrapping routines in parallel for greater computational efficiency  
   
- Bug fixes
    - plot_irf() plots all targets by default 
    - test-threvhold_var.R changed to test-threshold_var.R

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
    - Convert errorConditiont() to stop()   

## Version 1.0.0 
(2021-04-02)

- Initial release 
