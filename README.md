# mbreaks
[![CRAN](http://www.r-pkg.org/badges/version/mbreaks)](https://cran.r-project.org/package=mbreaks)[![Downloads](https://cranlogs.r-pkg.org/badges/mbreaks)](https://cran.r-project.org/package=mbreaks) [![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/mbreaks)](https://cran.r-project.org/package=mbreaks)


## Installation


```r
install.packages("mbreaks")
```

## Load Package

Once package has been installed it can be loaded. 
```{r}
library(mbreaks)
```

## Examples 


### US real interest rate data available in `mbreaks`

```r
y <- mbreaks::rint

z <- matrix(1,nrow(y),1)

x <- matrix(0,0,0)

m <- 1
n <- 1
M <- 3
N <-2


trm <- 0.10


st <- proc.time()

mbreaks::pslr0(y, m, trm, z) 

mbreaks::pslr1(y, n, trm, z) 

mbreaks::pslr2(y, m, n, trm, z) 

mbreaks::pslr3(y, m, n, trm, z) 

mbreaks::pslr4(y, m, n, trm, z) 

mbreaks::pslr00(y, M, trm, z)

mbreaks::pslr5(y, N, trm, z)

mbreaks::pslr6(y, m, N, trm, z)

mbreaks::pslr7(y, M, n, trm, z)

mbreaks::pslr8(y, M, N, trm, z)

mbreaks::pslr9(y, m, n, trm, z) 

mbreaks::pslr10(y, m, n, trm, z) 

end <- proc.time() - st
print(end)

```

### US inflation data available in `mbreaks`


## References

**Bai, Jushan & Pierre Perron (1998),**
Estimating and Testing Linear Models with Multiple Structural Changes,
_Econometrica_, vol 66, 47-78. 
[https://doi.org/10.2307/2998540](https://doi.org/10.2307/2998540)

**Bai, J. and Perron, P. (2003),** 
Computation and analysis of multiple structural change models,
_Journal of Applied Econometrics_, 18: 1-22. 
[https://doi.org/10.1002/jae.659](https://doi.org/10.1002/jae.659)

**Bai, J. and Perron, P. (2003),** 
Critical values for multiple structural change tests, 
_The Econometrics Journal_, 6: 72-78. [https://doi.org/10.1111/1368-423X.00102](https://doi.org/10.1111/1368-423X.00102)

**Perron, Pierre, Yohei Yamamoto, and Jing Zhou (2020),** 
Testing Jointly for Structural Changes in the Error Variance and Coefficients of a Linear Regression Model,
_Quantitative Economics_, vol 11, 1019-1057. 
[https://doi.org/10.3982/QE1332](https://doi.org/10.3982/QE1332)

