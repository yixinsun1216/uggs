# uggs

A package that builds on the [bcaboot](https://github.com/bnaras/bcaboot) package to compute bias corrected and accelerated bootstrap confidence limits . 
By implementing parallel computing using the [furrr](https://github.com/DavisVaughan/furrr) package, we are able to make substantial speed improvements.

For further information the bca bootstrapping method, refer to 
* Computer Age Satistical Inference found [here](https://web.stanford.edu/~hastie/CASI/order.html)
* Or the [manual](https://statweb.stanford.edu/~ckirby/brad/papers/2018Automatic-Construction-BCIs.pdf) for `bcaboot`

# Installation

To install, use `devtools::install_github("yixinsun1216/uggs")`

# Example
```
 library(lfe)
 library(uggs)
 
 ## create covariates
 x1 <- rnorm(1000)
 x2 <- rnorm(length(x1))
 
 ## fixed effects
 fe <- factor(sample(20, length(x1), replace=TRUE))
 
 ## effects for fe
 fe_effs <- rnorm(nlevels(fe))
 
 ## creating left hand side y
 u <- rnorm(length(x1))
 y <- 2 * x1 + x2 + fe_effs[fe] + u
 
 # create dataframe to pass into uggs
 df_test <- as.data.frame(cbind(y, x1, x2, fe))
 
 # function that returns parameter of interest, x1
 est_test <- function(df){
 	m <- felm(y ~ x1 + x2 | fe, df)
 	as.numeric(coef(m)["x1"])
 }
 
 x1_boot <- uggs(df_test, 1000, est_test, jcount = 40, jreps = 5)
 
 x1_boot
$`limits`
      bca      std      pct   jacksd     
0.025 1.977667 1.975631 0.029 0.004010308
0.05  1.987274 1.985915 0.054 0.001749787
0.1   1.996758 1.99777  0.102 0.004436592
0.5   2.038755 2.039592 0.486 0.001477231
0.9   2.08201  2.081414 0.902 0.002506374
0.95  2.096626 2.09327  0.954 0.004449781
0.975 2.106012 2.103553 0.979 0.001599195

$`stats`
       theta       sdboot          z0          a    sdjack
est 2.039592 0.0326336866 -0.01754730 0.02673691 0.0315359
jsd 0.000000 0.0007583239  0.03020335 0.00000000 0.0000000

$B.mean
[1] 1000.000000    2.039816

$ustats
    ustat       sdu 
2.0393682 0.1367483 
 ```
