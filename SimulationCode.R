# Code for example 1
################
library(gstat) #
library(ggpubr)
library(sparseMVN)#
library(Matrix)
library(PenCoxFrail)
library(msos)
library(tidyverse)
library(truncnorm)
################
rm(list=ls())
# Load external functions
source("utilityFunction.R")
source("Sampling3.R")
################
#Simulate data
################
# Spatial dormain
xy.xy <- expand.grid((1:150)/150, (1:150)/150)
set.seed(1234456)
# Data location
xy = data.frame(x=sample(xy.xy[,1],150),y=sample(xy.xy[,2],150))
names(xy) <- c('x','y')
## Set parameters ##
sigma_sqr = 2.5^2
phi_k1     = 2.5
phi_k2     = 3.0
phi_k3     = 5
tau_1     = 0.2
tau_2       =0.2
sigma_Sx = sigma_sqr *( (xy[,1]-(max(xy[,1])-min(xy[,1]))/2)^2+
              (xy[,2]-(max(xy[,2])-min(xy[,2]))/2)^2 )
sigma_S = ifelse(sigma_Sx<1,1,sigma_Sx)
Eulean  = function(p, vp){
  # p is a 1 x2 dataframe 
  # vp is an n x 2 dataframe
  aux = sqrt((vp[,1]-p[1,1])^2+(vp[,2]-p[1,2])^2)
  return(aux)
}
SIGMA1=SIGMA2=SIGMA3 = matrix(NA,nrow(xy),nrow(xy))

for (i in 1:nrow(xy)) {
  SIGMA1[i,] =sigma_S[i] * sigma_S*exp(-phi_k1 *( Eulean(xy[i,],xy) ) )
  SIGMA2[i,] =sigma_S[i] * sigma_S*exp(-phi_k2 *( Eulean(xy[i,],xy) ) )
  SIGMA3[i,] =sigma_S[i] * sigma_S*exp(-phi_k3 *( Eulean(xy[i,],xy) ) )
}
mu = matrix(0,nrow = nrow(xy),1)
SIGMA1 = as(SIGMA1,"sparseMatrix")
SIGMA2 = as(SIGMA2,"sparseMatrix")
SIGMA3 = as(SIGMA3,"sparseMatrix")
CH1=Cholesky(SIGMA1)
CH2=Cholesky(SIGMA2)
CH3=Cholesky(SIGMA3)

#  Z1 = t(rmvn.sparse(1,(mu),CH1,prec=FALSE))
#  Z2 = t(rmvn.sparse(1,(mu),CH2,prec=FALSE))
#  Z3 = t(rmvn.sparse(1,(mu),CH3,prec=FALSE))

###
n1=n2=nrow(xy)
Beta1 <- matrix(c(2,-2,1))
Beta2 <- matrix(c(3.5,1,-2.5,0.1))
p1 <- nrow(Beta1)
p2 <- nrow(Beta2)

x11 <- rnorm(n1,0,1)
x12 <- rnorm(n1,0,1)
x13 <- rnorm(n1,0,1)

X1 <- t(as.matrix(cbind(x11,x12,x13)))

x21 <- rnorm(n2,0,2)
x22 <- rnorm(n2,0,2)
x23 <- rnorm(n2,0,2)
x24 <- rnorm(n2,0,2)

X2 <-  t(as.matrix(cbind(x21,x22,x23,x24)))


# Non-linear effect (random walk 1) 
z1=vector("numeric",n1)
z1[1]=rnorm(1,0,0.5)
for (i in 2:n1) {
  z1[i]=rnorm(1,z1[i-1],0.5)
}
z2=vector("numeric",n2)
z2[1]=rnorm(1,0,0.5)
for (i in 2:n2) {
  z2[i]=rnorm(1,z2[i-1],0.5)
}

CH1_response =Cholesky( Diagonal(nrow(xy),1/tau_1))
CH2_response =Cholesky( Diagonal(nrow(xy),1/tau_2))
nsamp = 10
Y1 =Y2= matrix(0,nrow(xy),300)
#Y2=NULL
e = 0.5
e.prop=e
# Simulate response variable
for (i in 1:300) {
  
  Z1 = t(rmvn.sparse(1,(mu),CH1,prec=FALSE))
  Z2 = t(rmvn.sparse(1,(mu),CH2,prec=FALSE))
  Z3 = t(rmvn.sparse(1,(mu),CH3,prec=FALSE))
  sele = rbinom(1,1,e)
  mu1 = t(X1)%*%Beta1 + z1 + (1-sele) * Z1  +  sele * Z3
  sele = rbinom(1,1,e) 
  mu2 = t(X2)%*%Beta2 + z2 + (1-sele) * Z2  +  sele * Z3
  for (j in 1:nrow(mu1)) {
    Y1[j,i] = rnorm(1,mu1[j,1],1/tau_1)
    Y2[j,i] = rnorm(1,mu2[j,1],1/tau_2)
  }
  
}
DATA = list(Y1,Y2, xy, X1, X2,z1,z2)
TrainID = sample(1:150,100)
TrainIDcol = sample(1:300,nsamp)
Y1 = DATA[[1]][TrainID,TrainIDcol]
Y2 = DATA[[2]][TrainID,TrainIDcol]
X1 = DATA[[4]][,TrainID]
X2 = DATA[[5]][,TrainID]
z1 = DATA[[6]][TrainID]
z2 = DATA[[7]][TrainID]

##################
#Preparation for inference
##################
# Extract only the location from the true field
datLoc1 <- DATA[[3]][TrainID,]
datLoc2 <- DATA[[3]][TrainID,]
set.seed(1234)
# Select random field location for the estimation process
#SampleFieldLoc1 =  as.data.frame( cbind(sample(seq(0,1,0.001),170),sample(seq(0,1,0.001),170)))
#SampleFieldLoc2 =  as.data.frame(cbind(sample(seq(0,1,0.001),170),sample(seq(0,1,0.001),170)))

# Select Equidistant spatial field location for the estimation process
samloc =expand.grid(seq(0.1,0.9,length.out=17),seq(0.1,0.9,length.out=10))
SampleFieldLoc1 =  as.data.frame(samloc)
SampleFieldLoc2 =  as.data.frame(samloc)

# Determine projection matrix A
A1 <- AGenerator(datLoc1,SampleFieldLoc1,psi=4,c=20)
A2 <- AGenerator(datLoc2,SampleFieldLoc2,psi=4,c=20)
A1 =A1/matrix(rowSums(A1),nrow=nrow(A1),ncol=ncol(A1))
A1= as(A1,"sparseMatrix")
A2 =A2/matrix(rowSums(A2),nrow=nrow(A2),ncol=ncol(A2))
A2= as(A2,"sparseMatrix")

YY1 = Y1
YY2 = Y2
n1= nrow(YY1)
n2= nrow(YY2)
########### Artificial Set ######## would be rewritting during estimation

Y <- list(YY1,YY2)
X= list(X1,X2)
A=list(A1,A2)
ndata=list(n1,n2)
########################
#Initialization
########################
# Fixed effect

############
Beta1=Matrix(rnorm( nrow(X[1][[1]])))
Beta2=Matrix(rnorm( nrow(X[2][[1]])))
Beta=list(Beta1,Beta2)

# spatial effect
Theta1=Matrix(rnorm( ncol(A[1][[1]])*ncol(YY1)),nrow =ncol(A[1][[1]]))
Theta2=Matrix(rnorm( ncol(A[2][[1]])*ncol(YY1)),nrow =ncol(A[2][[1]] ))
Theta=list(Theta1,Theta2)

# Nonlinear effect
bz1= TrainID
bz2= TrainID
#cubic spline basis matrix
B1 <- bs.design(x=bz1, xl=min(bz1), xr=max(bz1), spline.degree=3, nbasis=8)
B2 <- bs.design(x=bz2, xl=min(bz2), xr=max(bz2), spline.degree=3, nbasis=8)

B=list(B1,B2)

Phi1=Matrix(rnorm( ncol(B[1][[1]])))
Phi2=Matrix(rnorm( ncol(B[2][[1]])))
Phi= list(Phi1,Phi2)

#Prior variance

Lambda=list()
Lambda[1]=rgamma(1,shape = 5,rate = 5)
Lambda[2]=rgamma(1,shape = 5,rate = 5)
Lambda[3]=rgamma(1,shape = 5,rate = 5)

# kappa and nu 
range=list(1.39,0.25,1)
kappa=list(2,2,2,2)
nu=list((range[[1]]*kappa[[1]])^2/8,
        (range[[2]]*kappa[[2]])^2/8,
        (range[[3]]*kappa[[3]])^2/8
)
# Prior covariance matrix
SigmaTheta1=cMatern(SampleFieldLoc1[,1:2],SampleFieldLoc1[,1:2],inc.dist=10,nu=nu[[1]],kappa=kappa[[1]],c=10) # Nearest-Neigb=10
SigmaTheta2=cMatern(SampleFieldLoc2[,1:2],SampleFieldLoc2[,1:2],inc.dist=10,nu=nu[[2]],kappa=kappa[[2]],c=10) # Neares-Neigb=10

# randomly draw shared field locations across the related studies
# to obtain the prior variance for the shared field
#set.seed(123445)
set.seed(34459999)
# Random field location
#Loc.sh = as.data.frame (cbind(sample(seq(0,1,0.001),170),sample(seq(0,1,0.001),170)))
#SigmaTheta.sh=cMatern(Loc.sh,Loc.sh,inc.dist=10,nu=nu[[3]],kappa=kappa[[3]],c=5) # Nearest-Neigb=5
### Equidistant location
# Equidistant field location 
samloc =expand.grid(seq(0.07,0.95,length.out=17),seq(0.07,0.95,length.out=10))
Loc.sh =  as.data.frame(samloc)
SigmaTheta.sh=cMatern(Loc.sh,Loc.sh,inc.dist=10,nu=nu[[3]],kappa=kappa[[3]],c=5) # Nearest-Neigb=5

SigmaTheta=list(SigmaTheta1,SigmaTheta2,SigmaTheta.sh)
QSigmaTheta=list(solve(SigmaTheta1),solve(SigmaTheta2),solve(SigmaTheta.sh))
SampleField = list(SampleFieldLoc1[,1:2],SampleFieldLoc2[,1:2],Loc.sh) # to be used to update kappa and nu

# Initialize probability between share and specific
e=0.1

# indicator variable
w=list()
w[1][[1]]=sample(c(0,1),nsamp,prob = c(0.5,0.5),replace = TRUE)
w[2][[1]]=sample(c(0,1),nsamp,prob = c(0.5,0.5),replace = TRUE)

# likelihood variance for each j,i
Tau1=rgamma(nsamp,shape = 0.5,rate = 2)
Tau2=rgamma(nsamp,shape = 0.5,rate = 2)

Tau=list(Tau1,Tau2)

# Alpha for DP prior
alpha1=alpha2=alpha3=2
alpha=list(alpha1,alpha2,alpha3)
# Fix a and b for alpha prior
a.alpha=b.alpha=list()
a.alpha[[1]]=a.alpha[[2]]=a.alpha[[3]]=2
b.alpha[[1]]=b.alpha[[2]]=b.alpha[[3]]=1

# Prior variance for Beta
QBeta =list()
QBeta[[1]]=Diagonal(length(Beta[[1]]),1/100000)
QBeta[[2]]=Diagonal(length(Beta[[2]]),1/100000)

# Prior variance for Beta
QPhi =list()
tau.phi=c(1/1000,1/1000)
QPhi[[1]]=as(RW_2(ncol(B[[1]]),1),"sparseMatrix")#Diagonal(length(Phi[[1]]),1/10)
QPhi[[2]]=as(RW_2(ncol(B[[2]]),1),"sparseMatrix")#Diagonal(length(Phi[[2]]),1/10)


#Hyper parameter for Tau
a.tau=5 #Assume equal value for all j and i.
b.tau=5
#Hyper parameter for lambda
a.lambda=5
b.lambda=5
#hyper parameter for phi precision
a.tauphi=2
b.tauphi=2
# prior for e
e.pi0=e.p1=0.0025
e.pie=1-e.pi0-e.p1
a.e=(1/3)*2*nsamp
b.e=(1/3)*2*nsamp
# hyper-parameter for kapppa
kappa.mu=1
kappa.sig=2
# hyper-parameter for nu
nu.mu =1
nu.sig=2
#####################
##### Begin inference
#####################
# Storage
# Number of posterior draws
posit = 3   # position of the shared effect
Num= 3000
count=0
burn=floor(Num/2)
BETA1 <- matrix(NA,nrow = nrow(Beta1),ncol=(Num-burn))
BETA2 <- matrix(NA,nrow = nrow(Beta2),ncol=(Num-burn))

PHI1  <- matrix(NA,nrow = nrow(Phi1),ncol=(Num-burn))
PHI2  <- matrix(NA,nrow = nrow(Phi2),ncol=(Num-burn))

THETA1<- array(NA,dim= c(nrow(Theta1),ncol(Theta1),(Num-burn)))
THETA2<- array(NA,dim= c(nrow(Theta1),ncol(Theta1),(Num-burn)))

W1    <- matrix(NA,nrow = nsamp,ncol=(Num-burn))
W2    <- matrix(NA,nrow = nsamp,ncol=(Num-burn))

TAU1  <- matrix(NA,nrow = nsamp,ncol=(Num-burn))
TAU2  <- matrix(NA,nrow = nsamp,ncol=(Num-burn))

EE    <- vector("numeric",Num-burn)
LAMBDA<- matrix(NA,nrow = 3,ncol=(Num-burn))
ALPHA <- matrix(NA,nrow = 3,ncol=(Num-burn))
KAPPA<- matrix(NA,nrow = 3,ncol=(Num-burn))
NU <- matrix(NA,nrow = 3,ncol=(Num-burn))
RESID <- vector("numeric",(Num-burn))
###############

while(count<Num){
  YY=TTheta=TTau=list()
  for (j in 1:length(w)) {
    
    YY[[j]]=Y[[j]][,which(w[[j]]==0)]
    TTheta[[j]]=Theta[[j]][,which(w[[j]]==0)]
    TTau[[j]]  = Tau[[j]][which(w[[j]]==0)]
    
  }
  
  YY[[posit]]=list(Y[[1]][,which(w[[1]]!=0)],Y[[2]][,which(w[[2]]!=0)])
  TTheta[[posit]]=list(Theta[[1]][,which(w[[1]]!=0)],Theta[[2]][,which(w[[2]]!=0)])
  TTau[[posit]]  = list(Tau[[1]][which(w[[1]]!=0)],Tau[[2]][which(w[[2]]!=0)])

  # From specific ---> specific and shared
  ww=list()
  ww[[1]]=ww[[2]]=integer()#=ww[[3]]=integer()
  ID=list(0,0,list(0,0))
  ID2=list(0,0,list(0,0))
  for (j in 1:2) {
    ###########
    TTheta[[j]] = as(TTheta[[j]],"sparseMatrix")
    TTheta[[posit]][[j]] = as(TTheta[[posit]][[j]],"sparseMatrix")
    YY[[posit]][[j]] = as(YY[[posit]][[j]],"sparseMatrix")
    YY[[j]] = as(YY[[j]],"sparseMatrix")
    t=ncol(as.matrix(YY[[j]]))
    if(is.null(t) | t==0 ){next}
    ##########
    for (i in seq_len(t)){
      AUX=PROBAB(j,i,YY,TTheta,TTau)
      n1 = length(AUX[[1]])
      n2 = length(AUX[[2]])
      s=n1+n2
    ###########
      rex.AUX= c(AUX[[1]],AUX[[2]])
      d=data.frame(sn=1:s,
                   rex.AUX,
                   id=c(rep(1,n1),rep(2,n2)),
                   idd=c(seq_len(n1),seq_len(n2)))
      d$id[n1]=d$idd[n1]=99
      d$id[s] =d$idd[s]=999
    ###########
      if(sum(d$rex.AUX)!=0){
        postsamp =sample(1:s,1,prob=(d$rex.AUX)/sum(d$rex.AUX)) 
      }else{
        postsamp =n1
      }
      
      ID[[j]][[i]]=postsamp
      if(d$id[postsamp]==99){
        ww[[j]][[i]] = 0
        invauxE= TTau[[j]][[i]]*t(A[[j]])%*%A[[j]]+Lambda[[j]]*QSigmaTheta[[j]]
        auxCH  = Cholesky(forceSymmetric(solve(invauxE)))
        auxmu  =  (TTau[[j]][[i]])*t(A[[j]])%*%(YY[[j]][,i]-B[[j]]%*%Phi[[j]]-t(X[[j]])%*%Beta[[j]])
        TTheta[[j]][,i] = rmvn.sparse(1,auxmu,auxCH,prec=FALSE)
        ID2[[j]][[i]]=99
      }else if(d$id[postsamp]==999){ # case 2
        ww[[j]][[i]] = 1
        invauxE=TTau[[j]][[i]]*t(A[[j]])%*%A[[j]]+Lambda[[posit]]*QSigmaTheta[[posit]] 
        auxCH  = Cholesky(forceSymmetric(solve(invauxE)))
        auxmu  =  (TTau[[j]][[i]])*t(A[[j]])%*%(YY[[j]][,i]-B[[j]]%*%Phi[[j]]-t(X[[j]])%*%Beta[[j]])
        TTheta[[j]][,i] = rmvn.sparse(1,auxmu,auxCH,prec=FALSE)
        ID2[[j]][[i]]=999
      }else if(d$id[postsamp]==1){ # case 3
        ww[[j]][[i]] = 0
        TTheta[[j]][,i]=TTheta[[j]][,d$idd[postsamp]]
        ID2[[j]][[i]]=d$idd[postsamp]
      }else{                                         # case 4
        ww[[j]][[i]] = 1
        TTheta[[j]][,i]=as(TTheta[[posit]][[j]][,d$idd[postsamp]],"sparseMatrix")
        ID2[[j]][[i]]=d$idd[postsamp]
      }
      
    } # end i
  } #end j
  
  ww0=list(integer(),integer(),integer())
  
  for (j in 1:2) {
    ##########
    TTheta[[posit]][[j]] = as(TTheta[[posit]][[j]],"sparseMatrix")
    YY[[posit]][[j]] = as(YY[[posit]][[j]],"sparseMatrix")
    t=ncol(TTheta[[3]][[j]])
    if(is.null(t) | t==0 ){ next }
    #########
    for (i in seq_len(t)) {
      
      AUX=  PROBAB0(j,i,YY,TTheta,TTau) 
      n0Tojspec = length(AUX[[1]])
      n0TOshareFrm1= length(AUX[[2]][[1]])
      n0TOshareFrm2= length(AUX[[2]][[2]])
      n0TOshareNew = length(AUX[[2]][[3]])
      s=n0Tojspec+n0TOshareFrm1+n0TOshareFrm2+n0TOshareNew
      rex.AUX= c(AUX[[1]],AUX[[2]][[1]],AUX[[2]][[2]],AUX[[2]][[3]])
      ############
      d=data.frame(sn=1:s,
                   rex.AUX,
                   id=c(rep(1,n0Tojspec),rep(2,n0TOshareFrm1),rep(3,n0TOshareFrm2),rep(5,n0TOshareNew)),
                   idd=c(seq_len(n0Tojspec),seq_len(n0TOshareFrm1),seq_len(n0TOshareFrm2),seq_len(n0TOshareNew)))
      d$id[n0Tojspec]=d$idd[n0Tojspec]=99
      d$id[s]=d$idd[s]=999
      ############
      if(sum(d$rex.AUX)!=0){
        postsamp= sample(1:s,1,prob=d$rex.AUX/sum(d$rex.AUX))
      }else{
        postsamp=s
      }
      
      ID[[3]][[j]][[i]]=postsamp
      
      if(d$id[postsamp]==99){
        ww0[[j]][[i]] = 0
        invauxE= TTau[[3]][[j]][[i]]*t(A[[j]])%*%A[[j]]+Lambda[[j]]*QSigmaTheta[[j]]
        auxCH  = Cholesky(forceSymmetric(solve(invauxE)))
        auxmu  =  (TTau[[3]][[j]][[i]])*t(A[[j]])%*%(YY[[3]][[j]][,i]-B[[j]]%*%Phi[[j]]-t(X[[j]])%*%Beta[[j]])
        TTheta[[3]][[j]][,i] = rmvn.sparse(1,auxmu,auxCH,prec=FALSE)
        ID2[[3]][[j]][[i]] = 99
      }else if(d$id[postsamp]==999){ # case 2
        ww0[[j]][[i]] = 1
        invauxE=TTau[[3]][[j]][[i]]*t(A[[j]])%*%A[[j]]+Lambda[[3]]*QSigmaTheta[[3]]
        auxCH  = Cholesky(forceSymmetric(solve(invauxE)))
        auxmu  =  (TTau[[3]][[j]][[i]])*t(A[[j]])%*%(YY[[3]][[j]][,i]-B[[j]]%*%Phi[[j]]-t(X[[j]])%*%Beta[[j]])
        TTheta[[3]][[j]][,i] = rmvn.sparse(1,auxmu,auxCH,prec=FALSE)
        ID2[[3]][[j]][[i]]=999
      }else if(d$id[postsamp]==1){
        ww0[[j]][[i]] = 0
        TTheta[[3]][[j]][,i]=as(TTheta[[j]][,d$idd[postsamp]],"sparseMatrix")
        ID2[[3]][[j]][[i]]=d$idd[postsamp]
      }else if(d$id[postsamp]==2){#postsamp>n0Tojspec & postsamp <= (n0Tojspec+n0TOshareFrm1)){ # case 3
        ww0[[j]][[i]] = 1
        TTheta[[3]][[j]][,i]=as(TTheta[[3]][[1]][,d$idd[postsamp]],"sparseMatrix")
        ID2[[3]][[j]][[i]]=d$idd[postsamp]
      }else if(d$id[postsamp]==3){#postsamp > (n0Tojspec+n0TOshareFrm1) & postsamp <= (n0Tojspec+n0TOshareFrm1+n0TOshareFrm2)){
        ww0[[j]][[i]] = 1
        TTheta[[3]][[j]][,i]=as(TTheta[[3]][[2]][,d$idd[postsamp]],"sparseMatrix")
        ID2[[3]][[j]][[i]]=d$idd[postsamp]
      }
      
      
    }# end for i
    #  } # end if
    
  }# end j
  
  if(length(ww[[1]])==nsamp){
    w[[1]]=c(ww[[1]])
  }else if(length(ww0[[1]])==nsamp){
    w[[1]]=c(ww0[[1]])
  }else{
    w[[1]]=c(ww[[1]],ww0[[1]])
  }
  
  if(length(ww[[2]])==nsamp){
    w[[2]]=c(ww[[2]])
  }else if(length(ww0[[2]])==nsamp){
    w[[2]]=c(ww0[[2]])
  }else{
    w[[2]]=c(ww[[2]],ww0[[2]])
  }
  
  if(ncol(as.matrix(YY[[1]]))==nsamp){
    Y[[1]] = as(YY[[1]],"sparseMatrix")
  }else if (ncol(as.matrix(YY[[3]][[1]]))==nsamp){
    Y[[1]] = as(YY[[3]][[1]],"sparseMatrix")
  }else{
    Y[[1]] = as(cbind(YY[[1]],YY[[3]][[1]]),"sparseMatrix")
  }
  
  if(ncol(as.matrix(YY[[2]]))==nsamp){
    Y[[2]] = as(YY[[2]],"sparseMatrix")
  }else if (ncol(as.matrix(YY[[3]][[2]]))==nsamp){
    Y[[2]] = as(YY[[3]][[2]],"sparseMatrix")
  }else{
    Y[[2]] = as(cbind(YY[[2]],YY[[3]][[2]]),"sparseMatrix")
  }
  
  if(ncol(as.matrix(TTheta[[1]]))==nsamp){
    Theta[[1]] = as(TTheta[[1]],"sparseMatrix")
  }else if (ncol(as.matrix(TTheta[[3]][[1]]))==nsamp){
    Theta[[1]] = as(TTheta[[3]][[1]],"sparseMatrix")
  }else{
    Theta[[1]] = as(cbind(TTheta[[1]],TTheta[[3]][[1]]),"sparseMatrix")
  }
  
  if(ncol(as.matrix(TTheta[[2]]))==nsamp){
    Theta[[2]] = as(TTheta[[2]],"sparseMatrix")
  }else if (ncol(as.matrix(TTheta[[3]][[2]]))==nsamp){
    Theta[[2]] = as(TTheta[[3]][[2]],"sparseMatrix")
  }else{
    Theta[[2]] = as(cbind(TTheta[[2]],TTheta[[3]][[2]]),"sparseMatrix")
  }
  

  if(length(TTau[[1]])==nsamp){
    Tau[[1]]   = TTau[[1]]
  }else if(length(TTau[[3]][[1]])==nsamp){
    Tau[[1]]   = TTau[[3]][[1]]
  }else{
    Tau[[1]]   = c(TTau[[1]],TTau[[3]][[1]])
  }
  
  if(length(TTau[[2]])==nsamp){
    Tau[[2]]   = TTau[[2]]
  }else if(length(TTau[[3]][[2]])==nsamp){
    Tau[[2]]   = TTau[[3]][[2]]
  }else{
    Tau[[2]]   = c(TTau[[2]],TTau[[3]][[2]])
  }

  ID3=list(integer(),integer())
  
  ##########################
  for (j in 1:length(w)) {
    
    sx =sample(1:70,1)
    aaux= unique(Theta[[j]][sx,])
    for (k in 1:length(aaux)) {
      ID3[[j]][which(Theta[[j]][sx,]==aaux[k])]=k
    }
  }
  rex= data.frame(id=rep(1:nsamp,2),spec=rep(c(1,2),c(nsamp,nsamp)),w=c(w[[1]],w[[2]]),
                  id1=c(ID3[[1]],ID3[[2]]))
  ##########################
  
  # Update cluster effects
  mu=lapply(1:length(w), function(j)t(X[[j]])%*%Beta[[j]]+B[[j]]%*%Phi[[j]])
  for (j in 1:length(w)) {
    
    aux = rex[rex$spec==j,] # Extract index j
    aux = aux[aux$w==0,]
    clusters = unique(aux$id1 )
    if(length(clusters)!=0){
      for (k in clusters) {
        #  if(k!=99){   # Update those indices not sampled from H
        clust.mem =aux$id[which(aux$id1==k)]
        aux.star.invE = sum(Tau[[j]][clust.mem])*t(A[[j]])%*%A[[j]]+Lambda[[j]]*QSigmaTheta[[j]]
        aux.star.E    = solve(aux.star.invE)
        aux.star.E    = as(aux.star.E  ,"sparseMatrix")
        Wjh           = Y[[j]][,clust.mem]-matrix(mu[[j]],nrow=nrow(mu[[j]]),ncol =length(clust.mem))
        Wjh           = t(A[[j]])%*%Wjh
        if(length(clust.mem)==1){
          Wjh = Wjh *(Tau[[j]][clust.mem])
        }else{
          Wjh = Wjh %*% diag(Tau[[j]][clust.mem])
        }
        Wjh           = matrix(rowSums(Wjh))
        muk            = aux.star.E %*% Wjh
        CH             = Cholesky(forceSymmetric(aux.star.E))
        Theta[[j]][,clust.mem]=t(rmvn.sparse(length(clust.mem),muk,CH,prec=FALSE))
        
      }# end for k
    } # end if
  } # end for
  aux = rex[rex$w !=0,] # Extract all index j for which w=0
  clusters = unique(aux$id1 )
  
  if(length(clusters)>0){
    
    for (k in clusters) {
      #      if(k!=999){
      aux1=aux[which(aux$id1==k),]
      invE0h=W0h=0
      for (kk in 1:nrow(aux1)) { # Sum over all members ---> [spec][id]
        invE0h  = invE0h+Tau[[aux1$spec[kk] ]][[aux1$id[kk]]] *t(A[[aux1$spec[kk]]])%*%A[[aux1$spec[kk]]]
        W0h     =  W0h+Tau[[aux1$spec[kk] ]][[aux1$id[kk]]] *t(A[[aux1$spec[kk]]]) %*% (Y[[aux1$spec[kk]]][,aux1$id[kk]]-
                                                                                          t(X[[aux1$spec[kk]]])%*%Beta[[aux1$spec[kk]]]-
                                                                                          B[[aux1$spec[kk]]]%*%Phi[[aux1$spec[kk]]] )
      }
      invE0h = invE0h+Lambda[[3]]*QSigmaTheta[[3]]
      E0h    = solve(invE0h)
      muk    = E0h %*% W0h
      CH     = Cholesky(forceSymmetric(E0h))
      Tkpost =t(rmvn.sparse(1,muk,CH,prec=FALSE))
      for (kk in 1:nrow(aux1)) { # Update for every member that share the same cluster
        Theta[[aux1$spec[kk]]][,aux1$id[kk]]= Tkpost
      }
      
      
    }
    
  }
  
  # Update alpha for all spatial priors 
  ns = list(nrow(rex[rex$spec==1&rex$w==0,]), nrow(rex[rex$spec==2&rex$w==0,]),
            nrow(rex[rex$spec==3&rex$w==0,]),nrow(rex[rex$w!=0,]))
  nstar=list()
  nstar[[1]]= length(unique(rex[rex$spec==1&rex$w==0&rex$id1!=99,]$id1))  + length(rex[rex$spec==1&rex$w==0&rex$id1==99,]$id1)
  nstar[[2]]= length(unique(rex[rex$spec==2&rex$w==0&rex$id1!=99,]$id1))  + length(rex[rex$spec==2&rex$w==0&rex$id1==99,]$id1)
  nstar[[3]]= length(unique(rex[rex$w!=0&rex$id1!=999,]$id1)) + length(rex[rex$w!=0&rex$id1==999,]$id1)
  # Update alpha for all specific effects
  for (j in 1:(length(w)+1)) { # add +1 for alpha0
    k= rbeta(1,alpha[[j]]+1,ns[[j]])
    c= (a.alpha[[j]]+nstar[[j]]-1)/(nsamp*(b.alpha[[j]])-log(k))
    pi.alpha= c/(1+c)
    rex.alpha=sample(c(1,2),1,prob = c(pi.alpha,1-pi.alpha))
    if(rex.alpha==1){
      alpha[[j]] =rgamma(1,a.alpha[[j]]+nstar[[j]],b.alpha[[j]]-log(k))
    }else{
      alpha[[j]] =rgamma(1,a.alpha[[j]]+nstar[[j]]-1,b.alpha[[j]]-log(k)) #let a.alpha>1, to avoid error
    }
  }
  
  
  
  # Update Beta and Phi at the same time
  for (j in 1:length(w)) {
    
    invE.beta= X[[j]]%*%t(X[[j]])*sum(Tau[[j]]) + QBeta[[j]]
    E.beta   = solve(invE.beta)
    E.beta   = as(E.beta ,"sparseMatrix")
    invE.phi= t(B[[j]])%*%B[[j]]*sum(Tau[[j]]) + tau.phi[j]*QPhi[[j]]
    E.phi   = solve(invE.phi)
    E.phi   = as(E.phi,"sparseMatrix")
    # Compute Wbeta
    Wbeta=Wphi=0
    for (i in 1:nsamp) {
      Wbeta=Wbeta+ Tau[[j]][[i]]*X[[j]]%*%(Y[[j]][,i]-B[[j]]%*%Phi[[j]]-A[[j]]%*%Theta[[j]][,i])
      Wphi=Wphi + Tau[[j]][[i]]*t(B[[j]])%*%(Y[[j]][,i]-t(X[[j]])%*%Beta[[j]]-A[[j]]%*%Theta[[j]][,i])
    }
    mubeta=E.beta %*%Wbeta
    CH= Cholesky(forceSymmetric(E.beta))
    Beta[[j]]=t(rmvn.sparse(1,mubeta,CH,prec=FALSE))
    muphi=E.phi %*%Wphi
    CH= Cholesky(forceSymmetric(E.phi))
    Phi[[j]]=t(rmvn.sparse(1,muphi,CH,prec=FALSE))
    #Phi[[j]]= Phi[[j]]-mean(Phi[[j]]) # sum-to-zero constriant
  }
  
  # update Tau for all J and i
  
  for (j in 1:length(w)) {
    muk = t(X[[j]])%*%Beta[[j]]+B[[j]]%*%Phi[[j]]
    for (i in 1:nsamp) {
      a.rex= 0.5*(nrow(Y[[j]])+2*a.tau)
      b.rex= (Y[[j]][,i]-muk-A[[j]]%*%Theta[[j]][,i])
      b.rex=0.5*(t(b.rex)%*%b.rex+2*b.tau)
      Tau[[j]][[i]]= rgamma(1,drop(a.rex),drop(b.rex))
    }
    
  } 
  
  ## Update lambda for all the spatial effects
  
  for (j in 1:(length(w)+1)) {
    if(j <=2){  # For specific lambdas
      a.rex= 0.5*(nstar[[j]]*nrow(Theta[[j]])+2*a.lambda)   # shape parameter
      # Calculate shape
      # Determine the whole structure
      rex1     = rex[rex$spec==j&rex$w==0,]                   # Extract structure only for j component and specific effect (w=0)
      lamb.clusta= unique(rex1$id1)                           #Determine the unique clusters. Note that, all (2*nsamp+2)/2 are unique as they were 
      
      if(length(lamb.clusta)>0){                              # and must be accounted for, since they were drawn from H. 
        rex.lambda=0
        for (k in 1:length(lamb.clusta)){                     # Iterate over all clusters
          i= rex1$id[which(rex1$id1==lamb.clusta[k])] 
          
          rex.lambda=rex.lambda + t(Theta[[j]][,i[1]])%*%QSigmaTheta[[j]]%*%Theta[[j]][,i[1]]
          
          
        }
        
      }
      
      rex.lambda = 0.5*(rex.lambda+2*b.lambda)
      Lambda[[j]]= rgamma(1,drop(a.rex),drop(rex.lambda))
      
    }else{
      # Extract structure only for j component
      rex1     = rex[rex$w!=0,]                # Extract j-th structure where w=0
      lamb.clusta= unique(rex1$id1)    
      
      if(length(lamb.clusta)>0){
        rex.lambda=0
        
        for (k in 1:length(lamb.clusta)){
          ii= rex1$id[which(rex1$id1 ==lamb.clusta[k])]             # 
          jj= rex1$spec[which(rex1$id1 ==lamb.clusta[k])] 
          # Allow only one member of a cluster to be used
          rex.lambda=rex.lambda + t(Theta[[jj[1]]][,ii[1]])%*%QSigmaTheta[[3]]%*%Theta[[jj[1]]][,ii[1]]
          
        }
        rex.lambda = 0.5*(rex.lambda+2*b.lambda)
        a.rex= 0.5*(nstar[[3]]*nrow(Theta[[jj[1]]])+2*a.lambda)########## Suspected      
        Lambda[[j]]= rgamma(1,drop(a.rex),drop(rex.lambda))   ###########
      }
    }
  }
  # Update precision for nonlinear effect
  
  #for (j in 1:length(w)) {
  # tau.phi[j] = rgamma(1,0.5*nrow(Phi[[j]])+a.tauphi , drop(0.5*t(Phi[[j]])%*%QPhi[[j]]%*%Phi[[j]]+b.tauphi) )
  #}
  #}
  
  # Update e. with beta prior
  # a.e and b.e,e.pi,e.p0,e.p1 are hyperparameters
  N1 = sum(unlist(w))
  N0 = length(unlist(w))-N1
  #if(N1==length(unlist(w))){
  #   e=1
  # }else if(N0==length(unlist(w))){
  #   e=0
  #  }else{
  e  = rbeta(1,N1+a.e,N0+b.e)
  # }
  
  ### Update kappa and nu
  # Update all kappa and nu
  # and Update all prior variance
  

  Sigmarex0=cMatern(SampleField[[3]],SampleField[[3]],inc.dist=10,nu=nu[[3]],kappa=kappa[[3]],c=5)
  CH0      = Cholesky(Lambda[[3]]^(-1)*Sigmarex0)
  ####
  test = TRUE
  while(test){
  kappa0.new = rtruncnorm(1,a=0,b=Inf, mean = kappa[[3]], sd=0.3)
  nu0.new    = rtruncnorm(1,a=0,b=Inf, mean = nu[[3]],    sd=0.3)
  ####
  Sigmarex0.new=cMatern(SampleField[[3]],SampleField[[3]],inc.dist=10,nu=nu0.new,kappa=kappa0.new,c=5)
  if(det(Sigmarex0.new)>0){
    CH0.new      = Cholesky(Lambda[[3]]^(-1)*Sigmarex0.new)
    test = FALSE
  }
  }
  aux0=0       # set accumulation variable for shared effect
  aux0.new=0
  ###
  
  
  for (j in 1:length(w)) {
    test = TRUE
    while(test){
    kappa.new = rtruncnorm(1,a=0,b=Inf, mean = kappa[[j]], sd=0.3)
    nu.new    = rtruncnorm(1,a=0,b=Inf, mean = nu[[j]],    sd=0.3)
    Sigmarexj.new=cMatern(SampleField[[j]],SampleField[[j]],inc.dist=10,nu=nu.new,kappa=kappa.new,c=5)
    if(det(Sigmarexj.new)>0){
      CH.new      = Cholesky(Lambda[[j]]^(-1)*Sigmarexj.new)
      test = FALSE
    }
    }
    auxj.new    =0
    auxj        =0 # set accumulation variable for independent effects
    ########
    Sigmarexj=cMatern(SampleField[[j]],SampleField[[j]],inc.dist=10,nu=nu[[j]],kappa=kappa[[j]],c=5)
    CHj      = Cholesky(Lambda[[j]]^(-1)*Sigmarexj)

    for (i in 1:nsamp) {
      if (w[[j]][[i]]==0){
        
        auxj=auxj+dmvn.sparse(Theta[[j]][,i],matrix(rep(0,nrow(CHj))),CHj,prec=FALSE)
        auxj.new=auxj.new+dmvn.sparse(Theta[[j]][,i],matrix(rep(0,nrow(CH.new))),CH.new,prec=FALSE)
      }else{
        
        aux0=aux0+dmvn.sparse(Theta[[j]][,i],matrix(rep(0,nrow(CH0))),CH0,prec=FALSE)
        aux0.new=aux0.new+dmvn.sparse(Theta[[j]][,i],matrix(rep(0,nrow(CH0.new))),CH0.new,prec=FALSE)
      }
    } 
    # Update kj,nuj here
    
    xi= auxj.new+dlnorm(kappa.new,kappa.mu,kappa.sig,log = TRUE)+dlnorm(nu.new,nu.mu,nu.sig,log = TRUE)+ 
      log(dtruncnorm(kappa[[j]],a=0,b=Inf, mean = kappa.new, sd=0.3))+ log(dtruncnorm(nu[[j]],a=0,b=Inf, mean = nu.new, sd=0.3))-
      auxj-dlnorm(kappa[[j]],kappa.mu,kappa.sig,log = TRUE)-dlnorm(nu[[j]],nu.mu,nu.sig,log = TRUE)-
      log(dtruncnorm(kappa.new,a=0,b=Inf, mean = kappa[[j]], sd=0.3))-log(dtruncnorm(nu.new,a=0,b=Inf, mean = nu[[j]], sd=0.3))
    xi= min(1,exp(xi))
    U=runif(1)
    if(U<xi){
      kappa[[j]]=kappa.new
      nu[[j]]=nu.new
    }else{
      kappa[[j]]=kappa[[j]]
      nu[[j]]   =nu[[j]]
    }
    # Update Sigma using this updated kappa and nu for specific effect
    SigmaTheta[[j]]= cMatern(SampleField[[j]],SampleField[[j]],inc.dist=10,nu=nu[[j]],kappa=kappa[[j]],c=5)
  }
  # Update kappa[[3]] and nu[[3]]
  
  xi0= aux0.new+dlnorm(kappa0.new,kappa.mu,kappa.sig,log = TRUE)+dlnorm(nu0.new ,nu.mu,nu.sig,log = TRUE)+ 
    log(dtruncnorm(kappa[[3]],a=0,b=Inf, mean = kappa0.new, sd=0.3))+ log(dtruncnorm(nu[[3]],a=0,b=Inf, mean = nu0.new, sd=0.3))-
    aux0-dlnorm(kappa[[3]],kappa.mu,kappa.sig,log = TRUE)-dlnorm(nu[[3]],nu.mu,nu.sig,log = TRUE)-
    log(dtruncnorm(kappa0.new,a=0,b=Inf, mean = kappa[[3]], sd=0.3))-log(dtruncnorm(nu0.new,a=0,b=Inf, mean = nu[[3]], sd=0.3))
  xi0= min(1,exp(xi0))
  U=runif(1)
  if(U<xi0){
    kappa[[3]]=kappa0.new
    nu[[3]]=nu0.new
  }else{
    kappa[[3]]=kappa[[3]]
    nu[[3]]   =nu[[3]]
  }
  # Update Sigma using this updated kappa and nu for shared effect
  SigmaTheta[[3]]= cMatern(SampleField[[3]],SampleField[[3]],inc.dist=10,nu=nu[[3]],kappa=kappa[[3]],c=5)
  
  SigmaTheta=list(SigmaTheta[[1]],SigmaTheta[[2]],SigmaTheta[[3]])
  QSigmaTheta=list(solve(SigmaTheta[[1]]),solve(SigmaTheta[[2]]),solve(SigmaTheta[[3]]))
  
  # Calculate mean sqr error 1/n[sum(Y-E(Y|D))^2]
  residual= 0
  for (j in 1:length(w)) {
    residualaux=lapply(1:nsamp, function(i)as.matrix(Y[[j]][,i]- t(X[[j]])%*%Beta[[j]]-
                                                    B[[j]]%*%Phi[[j]]- A[[j]]%*%Theta[[j]][,i]))
    residual=residual+ mean(unlist(residualaux)^2)
  }
  
  #### Storage 
  
  ###Save
  count=count+1
  #burn=0
  if(count>burn){
    storindx=count-burn
    # Save posterior draws for the marginals posterior distribution.
    BETA1[,storindx] <- Beta[[1]][,1]
    BETA2[,storindx] <- Beta[[2]][,1]
   
    PHI1[,storindx] <- Phi[[1]][,1]
    PHI2[,storindx] <- Phi[[2]][,1]
 
    THETA1[,,storindx] <- Theta[[1]] %>% as.matrix()
    THETA2[,,storindx] <- Theta[[2]] %>% as.matrix()
 
    W1[,storindx] <- w[[1]]
    W2[,storindx] <- w[[2]]
 
    TAU1[,storindx]<- Tau[[1]]
    TAU2[,storindx]<- Tau[[2]]
 
    EE[storindx]   <- e
    LAMBDA[,storindx]<- matrix(c(Lambda[[1]],Lambda[[2]],Lambda[[3]]))
    ALPHA[,storindx] <- matrix(c(alpha[[1]],alpha[[2]],alpha[[3]]))
    KAPPA[,storindx] <- matrix(c(kappa[[1]],kappa[[2]],kappa[[3]]))
    NU[,storindx]    <- matrix(c(nu[[1]],nu[[2]],nu[[3]]))
    RESID[storindx]  <- residual
  }
  cat(sum(w[[1]]),sum(w[[2]]),"\n")
  cat(count,"\n")
  print(Beta[[1]])
}
##########
par(mfrow=c(3,3))
plot(BETA1[1,] ,type="l")
abline(h=mean(BETA1[1,]),col="green")
plot(BETA2[2,] ,type="l")
abline(h=mean(BETA2[2,]),col="green")
plot(PHI1[3,] ,type="l")
abline(h=mean(PHI1[3,]),col="green")
plot(THETA1[1,1,] ,type="l")
abline(h=mean(THETA1[1,1,]),col="green")
plot(THETA2[1,1,] ,type="l")
abline(h=mean(THETA2[1,1,]),col="green")
plot(EE ,type="l")
abline(h=mean(EE),col="green")
plot(KAPPA[1,] ,type="l")
abline(h=mean(KAPPA[1,]),col="green")
plot(NU[1,] ,type="l")
abline(h=mean(NU[1,]),col="green")
plot(LAMBDA[1,] ,type="l")
abline(h=mean(LAMBDA[1,]),col="green")
