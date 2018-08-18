# rasteraggregation
Weighted Resampling Algorithms in Rcpp

Download here: https://github.com/danilomalzao/rasteraggregation/archive/master.zip
Than extract this to: "C:/"

## Install the plugin (copying this into R)
```
Rcpp::compileAttributes("C:/rasteraggregation-master/", verbose = TRUE)

install.packages("C:/rasteraggregation-master/", repos = NULL, type = "source")
```

## In Dinamica EGO
### If R is installed:
Just open the following submodel: C:/rasteraggregation-master/submodel/sumAlgorithmSubmodel.ego
Run and you got a working Example.

### IF R is not installed:
Follow: https://csr.ufmg.br/dinamica/dokuwiki/doku.php?id=calculate_r_expression#local_r_installation
To install R in dinamica.

## To run directly in R
```
library(rasteraggregation)
aggregation_resamplingSum("C:/rasteraggregation-master/submodel/pop_density_estimate_2015.tif", "newMap.tif")
```
