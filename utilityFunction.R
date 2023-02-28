plotspatial <- function(data,...){
  # Function to plot spatialpoint data frame using ggplot
  # input: dataFrame with x,y,value
  # output: plot object
   require(sp)
   require(inlabru)
  if(!is.data.frame(data)) stop("not a dataframe")
  if(ncol(data)!=3) stop("Number of columns must be 3 (x,y,z)")
  
  aux1 = list(...)
  aux2 = names(aux1)
  if('xlab' %in% aux2) {
    xlab =aux1$xlab
  }else{xlab=''}
  if('ylab' %in% aux2) {
    ylab =aux1$ylab
  }else{ylab=''}
  colnames(data) = c('x','y','z')
  gridded(data) = ~x+y
  plt = ggplot()+gg(data,aes(fill=z))+
    scale_fill_viridis_c(option = "H", direction = 1,name="")+#scale_fill_gradient(low="yellow", high="red")+
    theme_bw()+theme(legend.title =element_blank())+labs(x=xlab,y=ylab)
  return(plt)
}

AGenerator <- function(dataa,datab,psi=1,c=NULL){
  # Function to compute A (projection matrix)
  # Input:
  #'dataa: dataframe-> data location (locations where data samples are collected)
  #'datab: dataframe-> field location (locations where field were sampled)
  #'psi: smoothness parameter
  #'c : truncation parameter
  #'Output
  #'Matrix A'
  c = round(c)
  if (!is.null(c)){
  if (c < 0 )stop("Conditions: c>=0",'\n')
  }
  if (psi <= 0) stop("Conditions: psi>0",'\n')
  
  if (ncol(dataa)!=2 | ncol(datab)!=2){
    if(ncol(dataa)<2 | ncol(datab)<2 ) {
      stop("Function needs two data columns",'\n')
    }else{
      dataa2=dataa[,c(1,2)]
      datab2=datab[,c(1,2)]
      cat("Using the first two columns",'\n')
    }
  }else{
    dataa2=dataa
    datab2=datab
  }
  # vectorized Euclidean distance function
  Eulean  = function(p, vp){
    # p is a 1 x2 dataframe 
    # vp is an n x 2 dataframe
    aux = sqrt((vp[,1]-p[1,1])^2+(vp[,2]-p[1,2])^2)
    return(aux)
  }

  d= matrix(unlist(lapply(1:nrow(dataa2),function(i)Eulean(data.frame(dataa2[i,]),data.frame(datab2)) )),
                 nrow=nrow(dataa2),ncol = nrow(datab2),byrow = TRUE)
  d= exp(-d/psi)
  if (!is.null(c)) {
  b=t(apply(d,1,order))
  for(i in 1:nrow(d)){
    d[i,b[i,c(1:(ncol(d)-c))]]=0
  }
  }else{
    d=d
  }

  #d=d/rowSums(d)
  d=as(d,"sparseMatrix")
  cat("Nearest neigbour ",c)
  return(d)
}

cMatern <- function(Loc1,Loc2,nu,inc.dist=2, kappa,c=5, UseNN=T) {
  ## Input
  # Loc1: Location dataframe, size (,2)
  # Loc2: Location dataframe size  (,2)
  # nu: parameter nu 
  # kappa: parameter kappa
  # Vectorize == T
  ## Output
  # Covariance matrix .
  dataa2=Loc1
  datab2=Loc2
  Eulean  = function(p, vp){
    # p is a 1 x2 dataframe 
    # vp is an n x 2 dataframe
    aux = sqrt((vp[,1]-p[1,1])^2+(vp[,2]-p[1,2])^2)
    return(aux)
  }
  
  h= inc.dist*matrix(unlist(lapply(1:nrow(dataa2),function(i)Eulean(dataa2[i,],datab2))),
            nrow=nrow(dataa2),ncol = nrow(datab2),byrow = TRUE) # to 
  denseCov=ifelse(h > 0, besselK(h * kappa, nu) * (h * kappa)^nu / 
           (gamma(nu) * 2^(nu - 1)), 1)
  
  if(UseNN){
  # Use nearest neighbor matrix to reduce matrix complexity
  NN <- NearestN(h,round(c))
  aux= computeAD(denseCov,NN)
  H=aux$A
  VV=aux$D
  Qtild <- solve(Diagonal(nrow(H))-H)%*%VV%*%t(solve(Diagonal(nrow(H))-H))
  Qtild=as(Qtild,"sparseMatrix")
  return(Qtild)
  }else{
    denseCov=as(denseCov,"sparseMatrix")
    return(denseCov)
  }
}
  
NearestN <- function(dmat,neigb){
  
  ## Input
  #  dmat: distance matrix
  #  neigb: number of neigbors
  ## Output: List of nearest neigbhors
  
  N=list()
  n=nrow(dmat)
  diag(dmat) <- NA
  
  for(i in 1 : n){
    roworder <-  (dmat[i,] %>% order) # order the row according to dist.
    roworderid <- roworder<i           # determin j < i
    roworder <-  roworder[roworderid]  # Select index of j <i
    
    if(i<=neigb){
      
      N[[i]] <- roworder
    }else {
      N[[i]] <- roworder[1:neigb]
    }
    
  } 
  return(N)
}


computeAD <- function(C,N){
  
  ## Input
  # C : dense covariance matrix
  # N: a list of Neigbhors
  ## Output
  ## Sparse matrix of A and D.
  
  n=nrow(C)
  A=D=matrix(0,n,n)
  
  for(i in 1:(n-1)) {
    A[i+1,N[[i+1]]] = solve(C[N[[i+1]],N[[i+1]]], C[N[[i+1]],i+1])
    D[i+1,i+1] = C[i+1,i+1] - sum(C[i+1,N[[i+1]]]*A[i+1,N[[i+1]]])
  }
  D[1,1]=C[1,1]
  A[1,]=0
  
  return(list(A=as(A,"sparseMatrix"),
              D=as(D,"sparseMatrix")))
}  


RW_2 <-  function(n,p){
  # Function to generate random walk two precision matrix
  # Input:
  # n is the size of the walk
  # p is the variance
  # Output:
  # Generates nxn precision matrix
  W <- matrix(0,n,n)
  W[1,1] <- 1
  W[2,2] <- 1
  W[3:n,1] <- -c(1:(n-2))
  for (i in 3:n) {
    W[i,2:(length((i-1):1)+1)] <- (i-1):1
  }
  Wt <- W%*%diag(rep(p,n))%*%t(W)
  Wq <- chol2inv(chol(Wt))
  return(Wq)
}


# Non-Linear effect
nonlinear =  function(z,B,PHI,ci=FALSE){
  # Function to plot non-linear function
  if (!ci){
    x=1:length(z)
    z= scale(z,scale =T)
    zhat=   B %*% rowMeans(PHI)
    zhat= scale(zhat,scale=T)
    dat = data.frame(z,zhat)
    p = ggplot(dat)+
      geom_line(aes(x=x,y=z))+geom_point(aes(x=x,y=z))+
      geom_line(aes(x=x,y=zhat),color="red",linetype="dashed")+
      geom_line(aes(x=x,y=zhat),color="red",linetype="dashed",size=1)+
      labs(x="Index",y="Effect")+
      theme_bw()
  }else{
    
    x=1:length(z)
    z= scale(z,scale =T)
    zhat=  rowMeans(B %*% PHI)
    VarEff= B %*%cov(t(PHI)) %*%t(B)
    VarEff = sqrt(diag(VarEff))
    zup =zhat+VarEff *2
    zlw =zhat-VarEff *2
    zhat= scale(zhat,scale=T)
    zup= scale(zup,scale=T)
    zlw= scale(zlw,scale=T)
    dat = data.frame(z,zhat,zup,zlw)
    p = ggplot(dat)+
      geom_line(aes(x=x,y=z))+geom_point(aes(x=x,y=z))+
      geom_line(aes(x=x,y=zhat),color="red",size=1)+
      geom_line(aes(x=x,y=zup),color="blue",linetype="dashed")+
      geom_line(aes(x=x,y=zlw),color="red",linetype="dashed")+
      labs(x="Index",y="Effect")+
      theme_bw()
    
  }
  
  return(p)
}
sdFun = function(M,r=T){
  # Function to calculate variance of matrix with col
  # r is replication direction
  if (r){
    n= ncol(M)
    s = rowSums(M)
    s2= rowSums(M^2)
    sig = s2-s^2/n
    sig= sig/(n-1)
    sig=sqrt(sig)
  }else{
    n= nrow(M)
    s = colSums(M)
    s2= colSums(M^2)
    sig = s2-s^2/n
    sig= sig/(n-1)
    sig=sqrt(sig)
  }
  
  return(sig)
  
}

#m=function(h){ifelse(h > 0, besselK(h * kappa, nu) * (h * kappa)^nu / 
#(gamma(nu) * 2^(nu - 1)), 1)}