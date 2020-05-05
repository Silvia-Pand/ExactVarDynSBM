rm(list = ls())
source("est_var_dyn_approximate.R")
source("est_var_dyn_exact.R")
load("example_data.RData")

tol = 10^-8
n = dim(YY)[1]                   # number of individuals
TT = dim(YY)[3]                  # number of time occasions
k = 3                          # number latent states
maxit = 5000

cat("\n Approximate algorithm\n")
varApprox = est_var_dyn_approximate(YY,k,start=0,tol=tol,maxit=maxit)

cat("\n Exact algorithm\n")
varExact = est_var_dyn_exact(YY,k,start=0,maxit=maxit)

