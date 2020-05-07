est_var_dyn_exact = function(Y,k,start=0,tol=10^-10,maxit=1000,Tau1=NULL,TAU=NULL){
  
  # Exact algorithm for variational estimation of the dynamic SBMs
  #INPUT:
  #Y      = array of dynamic social networks data of dimension nxnxTT
  #k      = number of latent blocks
  #start  = starting values
            # start = 0 -- deterministic assignment to the classes via Kmeans
            # start = 1 -- random assignment to the classes
            # start = 2 -- external input
  #tol    = tollerance level for convergence (optional)
  #maxit  = max number of iterations (optional)
  #Tau1   = intial values for tau(i,q) when start=2
  #TAU    = initial values for tau(t,i,q,q') when start=2
  
  #OUTPUT:
  #piv    = intial probabilities vector
  #Pi     = transition probabilities matrix
  #Psi    = connectivity paramters matrix
  #J      = target function
  #cl     = classification of nodes
  #Tau1   = tau(i,q)
  #TAU    = tau(t,i,q,q')
  #lk     = log-likelihood 
  
  ## preliminaries
  n = dim(Y)[1]
  TT = dim(Y)[3]

  # Starting values
  
  # start = 0 -- deterministic assignment to the classes via Kmeans
  # start = 1 -- random assignment to the classes
  # start = 2 -- external input
  
  # Starting values for Tau1 and TAU
	if(start==0){
		YY = Y
		for(t in 1:TT) diag(YY[,,t]) = 0
		YYY = YY[,,1]
		for(t in 2:TT) YYY = rbind(YYY,YY[,,t])
		cl0 = matrix(kmeans(YYY, k, nstart=100)$cl,n,TT)
		Tau1 = matrix(0, n, k)
		for(i in 1:n) Tau1[i,cl0[i,1]]  = 1
		TAU = array(0,c(n,k,k,TT))
		for(i in 1:n) for(t in 2:TT) TAU[i,,cl0[i,t],t]  = 1
	}else if(start==1){
		Tau1 = matrix(runif(n*k), n, k)
		Tau1 = (1/rowSums(Tau1))*Tau1
		TAU = array(runif(n*k*TT), c(n,k,k,TT))
		TAU[,,,1] = 0
		for(i in 1:n) for(t in 1:TT) TAU[i,,,t] = (1/rowSums(TAU[i,,,t]))*TAU[i,,,t]
    }
    if(start==2) if(is.null(Tau1) | is.null(TAU)) stop("initial value for tau must be given in input")
    
    Tau1c = pmax(Tau1,10^-300)
    TAUc = pmax(TAU,10^-300)
  
# marginal tau and derivatives
	Taum = array(0,c(n,k,TT))
	Taum[,,1] = Tau1
	for(i in 1:n) for(t in 2:TT) Taum[i,,t] = t(TAU[i,,,t])%*%Taum[i,,t-1]
	Taumc = pmax(Taum,10^-300)
	d1Taum = array(0,c(n,k,TT,k))  # with respect the initial tau
	for(i in 1:n){
		d1Taum[i,,1,] = diag(k)
		for(t in 2:TT) d1Taum[i,,t,] = t(TAU[i,,,t])%*%d1Taum[i,,t-1,]
	}
	d2Taum = array(0,c(n,k,TT,k,k,TT))  # with respect the transition tau
  	for(i in 1:n) for(s in 2:TT){
  		for(u in 1:k) d2Taum[i,,s,u,,s] = Taum[i,u,s-1]*diag(k)
		if(s<TT) for(t in (s+1):TT) for(u in 1:k) d2Taum[i,,t,u,,s] = t(TAU[i,,,t])%*%d2Taum[i,,t-1,u,,s]
  	}
 
# starting values of the parameters
	  Psi = matrix(NA, k,k)
    for(u in 1:k) for(v in 1:k){
    	num = den = 0
    	for(i in 1:(n-1)){ 
    		ind = (i+1):n
    		for(t in 1:TT){
    			num = num + Taum[i,u,t] * sum(Taum[ind,v,t] * Y[i,ind,t]) + 
    			                              Taum[i,v,t] * sum(Taum[ind,u,t] * Y[ind,i,t])
    			den = den + Taum[i,u,t] * sum(Taum[ind,v,t]) + Taum[i,v,t] * sum(Taum[ind,u,t])
    		}
    	}
    	Psi[u,v] = num/den
    }
	  Psic = pmax(Psi,10^-300)
	  
	#piv = colMeans(Tau1)
	piv = rep(1/k,k)  
	pivc = pmax(piv,10^-300)
	
	# Pi = matrix(0,k,k)
	# for(u in 1:k){
	# 	for(v in 1:k) for(t in 2:TT) Pi[u,v] = Pi[u,v]+sum(Taum[,u,t-1]*TAU[,u,v,t])
	# 	Pi[u,] = Pi[u,]/sum(Pi[u,])
	# }
	Pi = matrix(1/k,k,k) 
	Pic = pmax(Pi,10^-300)
	
# compute target function
	J = sum(Tau1%*%log(pivc))-sum(Tau1*log(Tau1c))
	for(i in 1:n) for(t in 2:TT) J = J+sum(Taum[i,,t-1]*TAU[i,,,t]*(log(Pic)-log(TAUc[i,,,t])))
	for(i in 1:(n-1)){
		jj = (i+1):n
		for(t in 1:TT) for(u in 1:k) for(v in 1:k){
			J = J + sum(Taum[i,u,t]*Taum[jj,v,t] *
			            (Y[i,jj,t] * log(Psic[u,v]) + (1-Y[i,jj,t])*log(1-Psic[u,v])))
		}
	}
  
# iterate until convergence
	J0 = J; it = 0
	cat("------------|-------------|-------------|-------------|\n")
	cat("      it    |      J      | discrepancy | rel. discr  |\n")
	cat("------------|-------------|-------------|-------------|\n")
	cat(sprintf("%11g",c(it,J)),"\n",sep=" | ")  
	while((abs(J-J0)/(abs(J0))>tol | it==0) & it<maxit){
		J0 = J; it = it+1
    
# E-step
# update Tau1
		for(i in 1:n){
			ltau = log(pivc)
			for(t in 2:TT){
				ltau = ltau+t(d1Taum[i,,t-1,])%*%rowSums(TAU[i,,,t]*(log(Pic)-log(TAUc[i,,,t])))
			}
			for(j in (1:n)[-i]) for(t in 1:TT){
				lPhi = Y[i,j,t] * log(Psic) + (1-Y[i,j,t])*log(1-Psic)
				ltau = ltau + t(d1Taum[i,,t,])%*%lPhi%*%Taum[j,,t]
			}
			ltau = ltau-max(ltau)
			Tau1[i,] = exp(ltau)/sum(exp(ltau))
			Tau1c[i,] = pmax(Tau1[i,],10^-300)			
			Taum[i,,1] = Tau1[i,]
			for(t1 in 2:TT) Taum[i,,t1] = t(TAU[i,,,t1])%*%Taum[i,,t1-1]
			Taumc[i,,] = pmax(Taum[i,,],10^-300)
			for(s1 in 2:TT){
				for(u1 in 1:k) d2Taum[i,,s1,u1,,s1] = Taum[i,u1,s1-1]*diag(k)
				if(s1<TT) for(t1 in (s1+1):TT) for(u1 in 1:k) d2Taum[i,,t1,u1,,s1] = t(TAU[i,,,t1])%*%d2Taum[i,,t1-1,u1,,s1]
			}
		}

# update TAU
		for(i in 1:n) for(s in 2:TT) for(u in 1:k){
			ltau = log(Pic[u,])
			if(s<TT) for(t in (s+1):TT) ltau = ltau+(t(d2Taum[i,,t-1,u,,s])%*%rowSums(TAU[i,,,t]*(log(Pic)-log(TAUc[i,,,t]))))/Taumc[i,u,s-1]
			for(j in (1:n)[-i]){
				for(t in s:TT){
					lPhi = Y[i,j,t] * log(Psic) + (1-Y[i,j,t])*log(1-Psic)
					ltau = ltau +(t(d2Taum[i,,t,u,,s])%*%lPhi%*%Taum[j,,t])/Taumc[i,u,s-1]
				}
			}
			ltau = ltau-max(ltau)
			TAU[i,u,,s] = exp(ltau)/sum(exp(ltau))
			TAUc[i,u,,s] = pmax(TAU[i,u,,s],10^-300)
			for(t1 in s:TT) Taum[i,,t1] = t(TAU[i,,,t1])%*%Taum[i,,t1-1]
			Taumc[i,,s:TT] = pmax(Taum[i,,s:TT],10^-300)
			for(t1 in s:TT) d1Taum[i,,t1,] = t(TAU[i,,,t1])%*%d1Taum[i,,t1-1,]
			for(s1 in s:TT){
				for(u1 in 1:k) d2Taum[i,,s1,u1,,s1] = Taum[i,u1,s1-1]*diag(k)
				if(s1<TT) for(t1 in (s1+1):TT) for(u1 in 1:k) d2Taum[i,,t1,u1,,s1] = t(TAU[i,,,t1])%*%d2Taum[i,,t1-1,u1,,s1]
			}
		}		
		
# M-step
		Psi = matrix(NA, k,k)
		for(u in 1:k) for(v in 1:k){
			num = den = 0
			for(i in 1:(n-1)){
	    		ind = (i+1):n
				for(t in 1:TT){
					num = num + Taum[i,u,t] * sum(Taum[ind,v,t] * Y[i,ind,t]) +
    			                Taum[i,v,t] * sum(Taum[ind,u,t] * Y[ind,i,t])
    				den = den + Taum[i,u,t] * sum(Taum[ind,v,t]) + Taum[i,v,t] * sum(Taum[ind,u,t])
    			}
    		}
    		Psi[u,v] = num/den
		}
		Psic = pmax(Psi,10^-300)
		piv = colMeans(Tau1)
		pivc = pmax(piv,10^-300)
		Pi = matrix(0,k,k)
		for(u in 1:k){
			for(v in 1:k) for(t in 2:TT) Pi[u,v] = Pi[u,v]+sum(Taum[,u,t-1]*TAU[,u,v,t])
			Pi[u,] = Pi[u,]/sum(Pi[u,])
		}
		Pic = pmax(Pi,10^-300)
		
# compute function J
		J = sum(Tau1%*%log(pivc))-sum(Tau1*log(Tau1c))
		for(i in 1:n) for(t in 2:TT) J = J+sum(Taum[i,,t-1]*TAU[i,,,t]*(log(Pic)-log(TAUc[i,,,t])))
		for(i in 1:(n-1)){
			jj = (i+1):n
			for(t in 1:TT) for(u in 1:k) for(v in 1:k){
				J = J + sum(Taum[i,u,t]*Taum[jj,v,t] *
	          	            (Y[i,jj,t] * log(Psic[u,v]) + (1-Y[i,jj,t])*log(1-Psic[u,v])))
			}
		}

# display output
		if(it/1 == floor(it/1)) cat(sprintf("%11g",c(it,J,J-J0,(J-J0)/abs(J0))),"\n",sep=" | ")
  }
  if(it/1>floor(it/1)) cat(sprintf("%11g",c(it,J,J-J0,(J-J0)/abs(J0))),"\n",sep=" | ")
  cat("------------|-------------|-------------|-------------|\n");
  
# final log-likelihood and classification
	cl = apply(Taum, c(1,3), which.max)
	lk = 0
	for(i in 1:n){
		lk = lk+log(pivc[cl[i,1]])
		for(t in 2:TT) lk = lk+log(Pi[cl[i,t-1],cl[i,t]])
	}
	for(i in 1:(n-1)) for(t in 1:TT){
		jj = (i+1):n
		lk = lk+sum(Y[i,jj,t] * log(Psic[cl[i,t],cl[jj,t]]) + (1-Y[i,jj,t])*log(1-Psic[cl[i,t],cl[jj,t]]))
	}

	  # output
  out = list(piv=piv,Pi=Pi,Psi=Psi,J=J,cl=cl,Tau1=Tau1,TAU=TAU,lk=lk)
  return(out)
  
  
}
