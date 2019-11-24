#Start Rserve
library(Rserve)
library(Rcpp)
library(tsAnomaly)
setwd('/')
run.Rserve(6311, args = NULL, config.file = "/Rserv.conf")
