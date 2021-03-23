# Generalized Network Structured Models with Mixed Responses subject to Measurement Error and Misclassification

2021-03-22

Data and code for reproducing the results in the article "Generalized Network Structured Models with Mixed Responses subject to Measurement Error and Misclassification""
by Qihuang Zhang and Grace Y. Yi.



## File Structure

Before implementing the code, need to install R package [GeneErrorMis](https://https://github.com/QihuangZhang/GeneErrorMis):

``` r
# install devtools if necessary
install.packages('devtools')

# install the GeneErrorMis package
devtools::install_github('QihuangZhang/GeneErrorMis')

# load
library(GeneErrorMis)
```


### Data Analysis
* [Analysis of Mice SNP data](https://github.com/QihuangZhang/GEEmix/blob/main/code/DataAnalysis/DataAnalysis_905_IVoptimcompcomp.R)


### Simulations
* [Setting 1](https://github.com/QihuangZhang/GNSM/blob/main/code/Simulation/Simulation2.R): GNSM without measurement error in comparison to RelNet method.
* Setting 2: GNSM with different measurement error processes
  * [Known Parameter](https://github.com/QihuangZhang/GNSM/blob/main/code/Simulation/Simulation3.R)
  * [Internal Validation](https://github.com/QihuangZhang/GNSM/blob/main/code/Simulation/Simulation4.R)
  * [External Validation](https://github.com/QihuangZhang/GNSM/blob/main/code/Simulation/Simulation5.R)

