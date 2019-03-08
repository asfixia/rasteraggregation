# rasteraggregation
Weighted Resampling Algorithms in Rcpp

Download here: https://github.com/danilomalzao/rasteraggregation/archive/master.zip
Than extract to: "C:/"

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
aggregation_resamplingSum("C:/rasteraggregation-master/submodel/pop_density_estimate_2015/pop_density_estimate_2015.tif", "C:/rasteraggregation-master/submodel/pop_density_estimate_2015_resampled.tif")
```

# Solving Errors
### Cannot access RWinLib URL
1) Download file: https://github.com/rwinlib/gdal2/archive/v2.2.3.zip
2) Edit by hand the winlibs.R:

Comment the download line

Put the downloaded file path

Comment the unlink command

```
   #download.file(sprintf("https://github.com/rwinlib/gdal2/archive/v%s.zip", VERSION), "lib.zip", quiet = TRUE)
   unzip("C:/Users/administrador.CSR/AppData/Local/Temp/2/gdal2-2.2.3.zip", exdir = "../windows")
   #unlink("lib.zip")
```
