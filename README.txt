1. comlasso_0.2.zip: R package
-it depends on the Rglpk / biom packages
-user should install manually.

2. Figure
- bsi data is publicly avaiable 
- the package includes the bsi data(./data/bsi_new.RData)
- FigureX.R: script for Figure X 

3. Table
- Table1.R, Table2.R: script for table 1 and 2
	* Note that the functions used in Table1.R , Table2.R are not included in "comlasso" package.
	  The function 'comlassoC'  in ./library/comLassoC.R only fit only  L2 loss function.
	  It depends on Rcpp/RcppArmadillo. The C code is compiled by Rtoods 3.4.
	  Currently the genlasso package is not available in CRAN, so it should be installed manually for the simulation.
- Table3.R: script for table 3 in the real data analysis
