# Demo

## Install CellTools

```R
start_time <- Sys.time()
devtools::install_github("TrakhtenbergLab/CellTools")
end_time <- Sys.time()
end_time-start_time
```
[1] Time difference of 48.56251 secs

## Read reference and query seurat object files into the R

```r
reference <- readRDS(“../vignette_data/reference.rds”) 
query <- readRDS(“../vignette_data/query.rds”)
```
You can download the reference and query data [here](https://basespace.illumina.com/s/L6R0Dqm6r7R6). The Atlas data which generated the reference and query data is also avaible to download from the same link.

## Denoise the Atlas Seurat Object
```r
start_time <- Sys.time()
reference <- DeNoise(ref = reference, lab = "type", coef = 1.5)
end_time <- Sys.time()
end_time-start_time
```
[1] Time difference of 22.7701 secs


## Find the optimal k-weight parameter for KNN
```r
start_time <- Sys.time()
reference <- OptiK(reference, lab = "type", range = c(5,50), dims = 50, perc = 0.2)
end_time <- Sys.time()
end_time-start_time
```
[1] Time difference of 17.58038 mins
```r
reference@misc$CellTools$opti_k
```
[1] 19

## Map query Seurat object to reference Seurat object
```r
start_time <- Sys.time()
result<-MapTo(reference,query,lab="type",k=19)
end_time <- Sys.time()
end_time-start_time
```
[1] Time difference of 15.09316 mins

## Accuracy of mapping query object to reference object

```r
length(which(result$predicted.id == result$type))/ncol(result)
```
[1] 0.9922059


## Optimized Tracing of Query to Reference Data
```r
start_time <- Sys.time()
correlation<-GlobalProbe(reference,lab="type")
end_time <- Sys.time()
end_time-start_time
```
[1] Time difference of 1.795279 mins





