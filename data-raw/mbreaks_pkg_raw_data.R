#' -----------------------------------------------------------------------------
#' US Real interest rate 1961Q1 - 1986Q4 
rint <- as.matrix(read.table('inst/extdata/real.csv', sep = ',', header = F))
colnames(rint) <- c('rint')
usethis::use_data(rint, overwrite = TRUE)
#' -----------------------------------------------------------------------------
#' Inflation series (quarterly) from 1959Q1 - 2018Q4
inflation <- as.matrix(read.table('inst/extdata/inflation.csv', sep = ',', header = F))
colnames(rint) <- c('inflation')
usethis::use_data(inflation, overwrite = TRUE)
#' -----------------------------------------------------------------------------

