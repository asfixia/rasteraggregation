# rasteraggregation
Weighted Resampling Algorithms in Rcpp

Extract this to:
C:/rasteraggregation

## Compile
```
Rcpp::compileAttributes("C:/rasteraggregation/", verbose = TRUE)
```

## Install rasterAgregation
```
install.packages("C:/rasteraggregation/", repos = NULL, type = "source")
```

## In Dinamica EGO
Follow: https://csr.ufmg.br/dinamica/dokuwiki/doku.php?id=calculate_r_expression#local_r_installation
To install R in dinamica, then use the submodel contained in:
C:/rasteraggregation/submodel/sumAlgorithmSubmodel.ego
And run.

## To run directly in R
```
library(rasteraggregation)
aggregation_resamplingSum("./pop_density_estimate_2015.tif", "./menor.tif")
```
