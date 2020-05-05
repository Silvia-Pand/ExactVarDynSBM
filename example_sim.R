rm(list = ls())
source("draw_sn_dyn.R")
source("est_var_dyn_approximate.R")
source("est_var_dyn_exact.R")

tol = 10^-8
n = 50                          # number of individuals
TT = 6                          # number of time occasions
k = 2                           # number latent states
maxit = 5000

# transition probs
rho=0.3
Pi0 = matrix(0,k,k)
for(u in 1:k) for(v in 1:k) Pi0[u,v] = rho^abs(v-u)
Pi0 = diag(1/rowSums(Pi0))%*%Pi0

# initial probs
piv0 = rep(1/k, k)


#connectivity parameters
Psi0 = matrix(0.03,k,k) 
grid = 0.3*runif(k,0.5,1.5)
diag(Psi0) = grid

sim = draw_sn_dyn(n,k,TT,piv0,Pi0,Psi0)
Y = sim$Y
Utrue = sim$U+1

YY = Y
for(t in 1:TT){
  diag(YY[,,t]) = NA
}

varExact = est_var_dyn_exact(YY,k,start=0,maxit=maxit)

print(table(varExact$cl,Utrue))


