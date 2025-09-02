# PoDA
Post-selection differential abundance analysis of microbiome compositional data

Install and load PoDA
```
devtools::install_github('weihaowang01/PoDA')
library(PoDA)
help(PODA)
```

An example
```
library(PoDA)
library(mvtnorm)
library(foreach)
library(doParallel)
data("example_data") #example data with 500 taxa and 100 samples.

Y=example_data$count # observed abundance matrix
u=example_data$interest # phenotype of interest
confoun=example_data$confoun # two confounders

rej=PODA(Y = Y, u = u, confoun = confoun, rep_time=10, tau = 1, fdrnomial=0.05, fdrmethod="BH", lib_cut=1000, pre_cut=0.2) #PoDA

rej$rej # the index of differential abundant taxa identified by PoDA
rej$rejtaxon.name # the names of differential abundant taxa identified by PoDA
rej$rej_atleast # the names of taxa that are considered differentially abundant taxa in at least one of the multiple information splitting replications.
rej$rej_freq # for all taxa, the number of times each was considered a differential abundant taxa across multiple information splitting replications.
rej$all_freq # for those taxa identified by PoDA as differential abundant taxa, the number of times each was considered a differential abundant taxa across multiple information splitting replications.
```
