# library(mvtnorm)
# library(foreach)
# library(doParallel)
source("~/PODA/R/poda-function/R/poda.R")
example_data=readRDS("~~/PODA/R/example_data.rds")  #example data with 500 taxa and 100 samples.

Y=example_data$count # observed abundance matrix
u=example_data$interest # phenotype of interest
confun=example_data$confun # all phenotyes (1st column is represents the phenotype of interest)

rej=PODA(Y = Y, u = u, confun = confun[,-1] ,rep_time=10,tau = 1,fdrnomial=0.05,fdrmethod="BH",lib_cut=1000, pre_cut=0.2) #PoDA

rej$rej
rej$rejtaxon.name
rej$rej_atleast
head(rej$rej_freq)
head(rej$all_freq)





