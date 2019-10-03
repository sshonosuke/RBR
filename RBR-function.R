library(statmod)
library(MCMCpack)

##-------------------------------------------------##
##  Robust Bayesian regression with Laplace prior  ##
##-------------------------------------------------##
RBR.L=function(Y,X,mcmc=1500,burn=500,gam=0.2){
  n=dim(X)[1]; p=dim(X)[2]

  # function for computing mode (MM-algorithm)
  MM=function(YY,XX,ww,uu,Beta=NA,Int=NA){
    if(is.na(Beta[1])){ Beta=rep(0,p) }
    if(is.na(Int)){ Int=0 }
    sig=1
    th=10^(-6); dd=1; count=1
    while(dd>th & count<100){
      # weight
      Mu=Int+as.vector(XX%*%Beta)
      val=ww*dnorm(YY,Mu,sig)^(gam)
      weight=val/sum(val)*n
      # update Int
      resid=as.vector(YY-XX%*%Beta)
      Int=sum(weight*resid)/sum(weight)
      # update Beta
      new.Beta=ginv(t(XX)%*%diag(weight)%*%XX+diag(sig^2/uu))%*%t(XX)%*%diag(weight)%*%(YY-Int)
      new.Beta=as.vector(new.Beta)
      # update sigma
      resid=as.vector(YY-Int-XX%*%new.Beta)
      new.sig=sqrt((1+gam)*sum(weight*resid^2)/sum(weight))
      # difference
      dd=sum(abs(c(new.Beta,new.sig)-c(Beta,sig)))/sum(abs(c(Beta,sig)))
      Beta=new.Beta; sig=new.sig
      count=count+1
    }
    return(c(Int,Beta,sig))
  }
  
  # initial values 
  Alpha.pos=c()
  Beta.pos=matrix(NA,mcmc,p)
  Sigma.pos=c()
  init=coef(lm(Y~X))
  Alpha=init[1]
  Beta=init[-1]
  Sigma=1
  U.pos=matrix(NA,mcmc,p); U=rep(10,p)
  Lam.pos=c()
  a=1; b=1; Lam=1
  
  # MCMC
  for(k in 1:mcmc){
    w=as.vector(rdirichlet(1,rep(1,n)))
    Para=MM(Y,X,w,U,Beta,Alpha)
    Alpha.pos[k]=Para[1]
    Beta.pos[k,]=Para[-c(1,p+2)]
    Sigma.pos[k]=Para[p+2]
    mu=sqrt(Lam^2/Beta^2)
    U=1/rinvgauss(p,mu,Lam^2)
    U.pos[k,]=U
    Lam=sqrt(rgamma(1,p+a,sum(U)/2+b))
    Lam.pos[k]=Lam
  }
  # summary
  om=1:burn
  Res=list(Alpha.pos[-om],Beta.pos[-om,],Sigma.pos[-om])
  names(Res)=c("Alpha","Beta","Sigma")
  return(Res)
}






##-------------------------------------------------##
##  Robust Bayesian regression with Horseshoe prior  ##
##-------------------------------------------------##
RBR.HS=function(Y,X,mcmc=1500,burn=500,gam=0.2){
  n=dim(X)[1]; p=dim(X)[2]
  
  # function for computing mode (MM-algorithm)
  MM=function(YY,XX,ww,uu,Beta=NA,Int=NA){
    if(is.na(Beta[1])){ Beta=rep(0,p) }
    if(is.na(Int)){ Int=0 }
    sig=1
    th=10^(-6); dd=1; count=1
    while(dd>th & count<100){
      # weight
      Mu=Int+as.vector(XX%*%Beta)
      val=ww*dnorm(YY,Mu,sig)^(gam)
      weight=val/sum(val)*n
      # update Int
      resid=as.vector(YY-XX%*%Beta)
      Int=sum(weight*resid)/sum(weight)
      # update Beta
      new.Beta=ginv(t(XX)%*%diag(weight)%*%XX+diag(sig^2/uu))%*%t(XX)%*%diag(weight)%*%(YY-Int)
      new.Beta=as.vector(new.Beta)
      # update sigma
      resid=as.vector(YY-Int-XX%*%new.Beta)
      new.sig=sqrt((1+gam)*sum(weight*resid^2)/sum(weight))
      # difference
      dd=sum(abs(c(new.Beta,new.sig)-c(Beta,sig)))/sum(abs(c(Beta,sig)))
      Beta=new.Beta; sig=new.sig
      count=count+1
    }
    return(c(Int,Beta,sig))
  }

  # Initial values
  Alpha.pos=c()
  Beta.pos=matrix(NA,mcmc,p)
  Sigma.pos=c()
  init=coef(lm(Y~X))
  Alpha=init[1]
  Beta=init[-1]
  Sigma=1
  U.pos=matrix(NA,mcmc,p); U=rep(1,p)
  Xi.pos=matrix(NA,mcmc,p); Xi=rep(1,p)
  Tau2=1  
  Eta=2   # scale parameter
  cc=1    # hyperparameter in prior for Tau2 
  
  # MCMC
  for(k in 1:mcmc){
    # parameters in regression model
    w=as.vector(rdirichlet(1,rep(1,n)))
    Para=MM(Y,X,w,U,Beta,Alpha)
    Alpha.pos[k]=Para[1]
    Beta.pos[k,]=Para[-c(1,p+2)]
    Sigma.pos[k]=Para[p+2]
    # parameters/latent variables in HS prior
    U=rinvgamma(p,1,Beta^2/2+1/Xi)
    U.pos[k,]=U
    Xi=rinvgamma(p,Eta+1,1/U+1/Tau2)
    Xi.pos[k,]=Xi
    Tau2=rinvgamma(1,cc+p/2+p*Eta,cc+sum(1/Xi))
  }
  # summary
  om=1:burn
  Res=list(Alpha.pos[-om],Beta.pos[-om,],Sigma.pos[-om])
  names(Res)=c("Alpha","Beta","Sigma")
  return(Res)
}








##-------------------------------------------------##
##              Bayesian Lasso                     ##
##-------------------------------------------------##
BL=function(Y,X,mc=2000,burn=500){
  n=dim(X)[1]; p=dim(X)[2]
  
  # initial values
  Int.pos=c(); Int=1
  Beta.pos=matrix(NA,mc,p); Beta=rep(0,p)
  sig.pos=c(); sig=1
  Lam.pos=c(); Lam=1
  U.pos=matrix(NA,mc,p);  U=rep(1,p)
  a=1; b=1    # prior for lambda^2
  cc=1   # prior for sigma^2
  
  # MCMC 
  for(k in 1:mc){
    # Int
    resid=as.vector(Y-X%*%Beta)
    Int=rnorm(1,mean(resid),sig/sqrt(n))
    Int.pos[k]=Int
    # Beta
    Yt=Y-Int
    invA=ginv(t(X)%*%X/sig^2+diag(1/U))
    Beta=mvrnorm(1,invA%*%t(X)%*%Yt/sig^2,invA)
    Beta.pos[k,]=Beta
    # sigma
    resid=as.vector(Y-Int-X%*%Beta)
    sig2=rinvgamma(1,n/2+cc,sum(resid^2)/2+cc)
    sig=sqrt(sig2)
    sig.pos[k]=sig
    # U
    mu=sqrt(Lam^2/Beta^2)
    U=1/rinvgauss(p,mu,Lam^2)
    U.pos[k,]=U
    # Lam
    Lam=sqrt(rgamma(1,p+a,sum(U)/2+b))
    Lam.pos[k]=Lam
  }
  # summary
  om=1:burn
  Res=list(Int.pos[-om],Beta.pos[-om,],sig.pos[-om])
  names(Res)=c("Alpha","Beta","Sigma")
  return(Res)
}

