# Sampling functions
PROBAB <- function(j,i,YY,TTheta,TTau){
  # Fucntion to calculate probabilities for the ith example for all j
  #Input:
  #'w is the the list of configuration
  #'i is the example we want to computes its probability
  # Output
  #' list of probabilities for all j and the i th example.
  #######################
  YY[[j]]=as(YY[[j]],"sparseMatrix")
  YY[[3]][[j]]=as(YY[[3]][[j]],"sparseMatrix")
  TTheta[[j]] = as(TTheta[[j]],"sparseMatrix")
  TTheta[[3]][[j]] = as(TTheta[[3]][[j]],"sparseMatrix")
  #######################
  
  n=ncol(YY[[j]])
  if(is.null(n)| n==0){

      q=integer()

  }else{
  q=rep(0,n)
  
  mu=lapply(1:length(w), function(j)t(X[[j]])%*%Beta[[j]]+B[[j]]%*%Phi[[j]])
  
  for (k in 1:n) {

    muk = mu[[j]]+A[[j]]%*%TTheta[[j]][,k]
    CH  = as((TTau[[j]][i])^(-1)*Diagonal(ndata[[j]]),"sparseMatrix")
    CH     = Cholesky(CH) 
    q[k]= drop((1-e)*exp(dmvn.sparse(t(YY[[j]][,i]),muk,CH,prec=FALSE)))
  
  }
  q[i]=0
  }
  ############# start: Calculate q0 --> q[n+1] ########### 
  
  invE = (TTau[[j]][i])*t(A[[j]])%*%A[[j]]+Lambda[[j]]*QSigmaTheta[[j]]
  E    = solve(invE)
  W    = (TTau[[j]][i])*t(A[[j]])%*%(YY[[j]][,i]-mu[[j]])
  d    = (TTau[[j]][i])*(t(YY[[j]][,i])%*%(YY[[j]][,i]-2*mu[[j]])+
                             t(t(X[[j]])%*%Beta[[j]])%*%t(X[[j]])%*%Beta[[j]]+
                             2*t(t(X[[j]])%*%Beta[[j]])%*%B[[j]]%*%Phi[[j]]+
                             t(B[[j]]%*%Phi[[j]])%*%B[[j]]%*%Phi[[j]]
  )
  auxlogExp = -0.5*(d- t(W)%*%E%*%W)
  auxlognorm = (0.5)*((nrow(E))*log(6.283185)+logdet(as.matrix(E)))-(0.5)*log( det(6.283185*(TTau[[j]][i])^(-1)*Diagonal(ndata[[j]]) ))-
    (0.5)*((nrow(SigmaTheta[[j]]))*log((6.283185*Lambda[[j]]^(-1)))+logdet(as.matrix(SigmaTheta[[j]]) ) )
  q[n+1]=drop(alpha[[j]]*(1-e)*exp(auxlognorm+auxlogExp))
  ############ End: Calculate q0 --> q[n+1] #########
  
  ########### Calculate shared q---> qs
  
  TTheta[[3]][[j]]= as.matrix(TTheta[[3]][[j]])
  
  ns = ncol(TTheta[[3]][[j]])
  test = is.null(ns) | ns==0
  if(test){
      qs = integer()
  }else{
  qsrex =NULL
  for (k in seq_len(ns)) {
    muk = mu[[j]]+A[[j]]%*%TTheta[[3]][[j]][,k]
    CH  = as((TTau[[j]][i])^(-1)*Diagonal(ndata[[j]]),"sparseMatrix")
    CH     = Cholesky(CH) 
    qsrex=c(qsrex,drop((e)*exp(dmvn.sparse(t(YY[[j]][,i]),muk,CH,prec=FALSE))))
  }
  qs=qsrex 
  }


  ############ Start: Calculate q0 for shared --> qs[j+1] #########
 
  invE = (TTau[[j]][i])*t(A[[j]])%*%A[[j]]+Lambda[[3]]*QSigmaTheta[[3]] 
  E= solve(invE)
  W    = (TTau[[j]][i])*t(A[[j]])%*%(YY[[j]][,i]-mu[[j]])
  d    = (TTau[[j]][i])*(t(YY[[j]][,i])%*%(YY[[j]][,i]-2*mu[[j]])+
                             t(t(X[[j]])%*%Beta[[j]])%*%t(X[[j]])%*%Beta[[j]]+
                             2*t(t(X[[j]])%*%Beta[[j]])%*%B[[j]]%*%Phi[[j]]+
                             t(B[[j]]%*%Phi[[j]])%*%B[[j]]%*%Phi[[j]]
  )
  auxlogExp = -0.5*(d- t(W)%*%E%*%W)
  auxlognorm = (0.5)*((nrow(E))*log(6.283185)+logdet(as.matrix(E)))-(0.5)*log( det(6.283185*(Tau[[j]][i])^(-1)*Diagonal(ndata[[j]]) ))-
    (0.5)*((nrow(SigmaTheta[[3]]))*log((6.283185*Lambda[[3]]^(-1)))+logdet(as.matrix(SigmaTheta[[3]])))
  if(test){
    qs=drop(alpha[[3]]*e*exp(auxlognorm+auxlogExp))
    }else{
      qs[ns+1]=drop(alpha[[3]]*e*exp(auxlognorm+auxlogExp))
    }
  
  
  return(list(q,qs))
}
  
# Funtion to update elements in the shared component. Jump across groups within the shared component
PROBAB0 <- function(j,i,YY,TTheta,TTau){
  # Fucntion to calculate probabilities for the ith example for all j
  #Input:
  #'w is the the list of configuration
  #'i is the example we want to computes its probability
  # Output
  #' list of probabilities for all j and the i th example.
  #
  #######################
  YY[[j]]=as(YY[[j]],"sparseMatrix")
  YY[[3]][[j]]=as(YY[[3]][[j]],"sparseMatrix")
  TTheta[[j]] = as(TTheta[[j]],"sparseMatrix")
  TTheta[[3]][[j]] = as(TTheta[[3]][[j]],"sparseMatrix")
  #######################
  
  qs=list()
  nk=list()
  for (jk in 1:length(w)) {
    
  YY[[3]][[jk]] =as.matrix(YY[[3]][[jk]])
  nk[[jk]]      =ncol(YY[[3]][[jk]])
  n= ncol(YY[[3]][[j]])
  ####
  test = is.null(nk[[jk]])| nk[[jk]]==0 | n==0
  if(test){
    
    qs[[jk]]=integer()

  }else{
    
  TTheta[[3]][[jk]] =as.matrix(TTheta[[3]][[jk]])
  qrex=NULL
  
  mu=lapply(1:length(w), function(j)t(X[[j]])%*%Beta[[j]]+B[[j]]%*%Phi[[j]])
  
  for (k in 1:nk[[jk]]) {
    
    muk = mu[[j]]+A[[j]]%*%TTheta[[3]][[jk]][,k]
    CH =as((TTau[[3]][[j]][i])^(-1)*Diagonal(ndata[[j]]),"sparseMatrix")
    CH     = Cholesky(CH) 
    ###
    if(j==jk & k==i){
    qrex=c(qrex, 0)
    ###
  }else{
    qrex=c(qrex, drop((e)*exp(dmvn.sparse(t(YY[[3]][[j]][,i]),muk,CH,prec=FALSE)))  )
     }
  }
  qs[[jk]]=qrex
  }
}
  ############# start: Calculate q00 --> as q[n+1] ###########
  
  invE = (TTau[[3]][[j]][i])*t(A[[j]])%*%A[[j]]+Lambda[[3]]*QSigmaTheta[[3]]
  E    = solve(invE)
  W    = (TTau[[3]][[j]][i])*t(A[[j]])%*%(YY[[3]][[j]][,i]-mu[[j]])
  d    = (TTau[[3]][[j]][i])*(t(YY[[3]][[j]][,i])%*%(YY[[3]][[j]][,i]-2*mu[[j]])+
                           t(t(X[[j]])%*%Beta[[j]])%*%t(X[[j]])%*%Beta[[j]]+
                           2*t(t(X[[j]])%*%Beta[[j]])%*%B[[j]]%*%Phi[[j]]+
                           t(B[[j]]%*%Phi[[j]])%*%B[[j]]%*%Phi[[j]]
  )
  auxlogExp = -0.5*(d- t(W)%*%E%*%W)
  auxlognorm = (0.5)*((nrow(E))*log(6.283185)+logdet(as.matrix(E)))-(0.5)*log( det(6.283185*(TTau[[3]][[j]][i])^(-1)*Diagonal(ndata[[j]]) ))-
    (0.5)*((nrow(SigmaTheta[[3]]))*log((6.283185*Lambda[[3]]^(-1)))+logdet(as.matrix(SigmaTheta[[3]]) ) )
  qs[[3]]=drop(alpha[[j]]*(e)*exp(auxlognorm+auxlogExp))
  ############ End: Calculate q00 #########
  
  ####################################################
  ########### Calculate specific probability ---> q###
  ####################################################
  # Only visit the corresponding j
  TTheta[[j]] = as.matrix(TTheta[[j]])
  
  n = ncol(TTheta[[j]])
   test = is.null(n) | n==0
  if(test){
  
    q=integer()

  }else{
  qsrex =NULL
  for (k in 1: n) {
    muk = mu[[j]]+A[[j]]%*%TTheta[[j]][,k]
    CH  = as((TTau[[3]][[j]][i])^(-1)*Diagonal(ndata[[j]]),"sparseMatrix")
    CH     = Cholesky(CH) 
    qsrex=c(qsrex,drop((1-e)*exp(dmvn.sparse(t(YY[[3]][[j]][,i]),muk,CH,prec=FALSE)))  )
  }
  q=qsrex 
  }
  
  
  ############ Start: Calculate q0 for specific --> q[n+1] #########
  
  invE = (TTau[[3]][[j]][i])*t(A[[j]])%*%A[[j]]+Lambda[[j]]*QSigmaTheta[[j]]
  E= solve(invE)
  W    = (TTau[[3]][[j]][i])*t(A[[j]])%*%(YY[[3]][[j]][,i]-mu[[j]])
  d    = (TTau[[3]][[j]][i])*(t(YY[[3]][[j]][,i])%*%(YY[[3]][[j]][,i]-2*mu[[j]])+
                           t(t(X[[j]])%*%Beta[[j]])%*%t(X[[j]])%*%Beta[[j]]+
                           2*t(t(X[[j]])%*%Beta[[j]])%*%B[[j]]%*%Phi[[j]]+
                           t(B[[j]]%*%Phi[[j]])%*%B[[j]]%*%Phi[[j]]
  )
  auxlogExp = -0.5*(d- t(W)%*%E%*%W)
  auxlognorm = (0.5)*((nrow(E))*log(6.283185)+logdet(as.matrix(E)))-(0.5)*log( det(6.283185*(TTau[[3]][[j]][i])^(-1)*Diagonal(ndata[[j]]) ))-
    (0.5)*((nrow(SigmaTheta[[j]]))*log((6.283185*Lambda[[j]]^(-1)))+logdet(as.matrix(SigmaTheta[[j]])))
  if(test){
    q = drop(alpha[[j]]*(1-e)*exp(auxlognorm+auxlogExp))
  }else{
    q[n+1]=drop(alpha[[j]]*(1-e)*exp(auxlognorm+auxlogExp))
  }
  
  return(list(q,qs))
}

