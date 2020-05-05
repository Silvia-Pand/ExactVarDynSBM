draw_sn_dyn = function(n, k, TT, piv,Pi,Psi){
  
  
  #  draw Markov chains
  U = matrix(0,n,TT)
  if(k>1){
    for(i in 1:n){
      U[i,1] = which(rmultinom(1,1,piv)==1)-1
      for(t in 2:TT) 	U[i,t] = which(rmultinom(1,1,Pi[U[i,t-1]+1,])==1)-1
    }	
  }
  
  # generate the network
  Y = array(0, c(n,n,TT))
  for(t in 1:TT) for(i in 1:(n-1)) for(j in (i+1):n){
    if(k==1) Y[i,j,t] = Y[j,i,t] = 1*(runif(1)<Psi)
    else Y[i,j,t] = Y[j,i,t] = 1*(runif(1)<Psi[U[i,t]+1, U[j,t]+1])
  }
  # ind = order(U)
  # U = U[ind]
  # Y = Y[ind,ind]
  
  out = list(U = U, Y = Y)
}