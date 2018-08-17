# rasteraggregation
Weighted Resampling Algorithms in Rcpp

//Compile
Rcpp::compileAttributes("./", verbose = TRUE)

//Install rasterAgregation
install.packages("./", repos = NULL, type = "source")

//Carregar rasteraggregation
library(rasteraggregation)
aggregation_resamplingSum("./pop_density_estimate_2015.tif", "./menor.tif")

