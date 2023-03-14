# mbreaks
[![CRAN](http://www.r-pkg.org/badges/version/mbreaks)](https://cran.r-project.org/package=mbreaks)[![Downloads](https://cranlogs.r-pkg.org/badges/mbreaks)](https://cran.r-project.org/package=mbreaks) [![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/mbreaks)](https://cran.r-project.org/package=mbreaks)


## Examples 


### US real interest rate data available in `mbreaks`

```r
library(mbreaks)

y <- mbreaks::rint

z <- matrix(1,nrow(y),1)

x <- matrix(0,0,0)

m <- 1
n <- 1
h <- 10


mbreaks::pslr0(y, m, h, z) 

mbreaks::pslr1(y, n, h, z) 

mbreaks::pslr2(y, m, n, h, z) 

mbreaks::pslr3(y, m, n, h, z) 

mbreaks::pslr4(y, m, n, h, z) 

mbreaks::pslr9(y, m, n, h, z) 

mbreaks::pslr10(y, m, n, h, z) 


```

## References

**Bai, Jushan & Pierre Perron (1998).**
Estimating and Testing Linear Models with Multiple Structural Changes,
_Econometrica_, vol 66, 47-78. 
[https://doi.org/10.2307/2998540](https://doi.org/10.2307/2998540)

**Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020).** 
Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model
_Quantitative Economics_, vol 11, 1019-1057. 
[https://doi.org/10.3982/QE1332](https://doi.org/10.3982/QE1332)

