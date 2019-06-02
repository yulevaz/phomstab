require(devtools)
require(roxygen2)

Rcpp::compileAttributes()
system("cp ./R/RcppExports.R.bkp ./R/RcppExports.R")
roxygen2::roxygenize()
devtools::document()
devtools::build()
system("R CMD INSTALL ../phomstab_0.0.0.9000.tar.gz")
