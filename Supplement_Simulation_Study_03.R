require(stats)
library(xtable)
set.seed(12345)

Greenberg.RR.data = function(nn,pp,cc,beta.1, alpha.00)
{

  X1 = rbinom(nn,1,0.5)
  X2 = rbinom(nn,1,0.5)
  
  chi.1 = cbind(1,X1,X2)
  
  alpha.0 = alpha.00[c(1:3)]
  alpha.1 = alpha.00[c(4:6)]
  
  b1 = chi.1%*%beta.1
  a0 = chi.1%*%alpha.0
  a1 = chi.1%*%alpha.1
  
  H1   = 1/(1+exp(-b1))
  HZ.0 = 1/(1+exp(-a0))
  HZ.1 = 1/(1+exp(-a1))
  
  Y1 = rbinom(nn,1,H1)
  Z0 = rbinom(nn,1,HZ.0)
  Z1 = rbinom(nn,1,HZ.1)
  ZZ = Y1*Z1+(1-Y1)*Z0
  
  TT   = rbinom(nn,1,pp)
  DD   = rbinom(nn,1,cc)
  Y1.0 = TT*Y1+(1-TT)*DD
  
  return(cbind(Y1, Y1.0, ZZ, chi.1))
}


### Conditional Expection of Y[i] given Y0[i] and ZZ[i] and X
CEY = function(Y1.0,ZZ,chi.1,pp,cc, Theta.hat)
{
  KK = ncol(chi.1)
  
  beta.h = Theta.hat[1:KK]
  
  alpha.00h = Theta.hat[(KK+1):(3*KK)]
  alpha.0h  = alpha.00h[1:KK]
  alpha.1h  = alpha.00h[(KK+1):(2*KK)]
  
  b1.h = chi.1%*%beta.h
  a0.h = chi.1%*%alpha.0h
  a1.h = chi.1%*%alpha.1h
  
  H1   = 1/(1+exp(-b1.h))
  HZ.0 = 1/(1+exp(-a0.h))
  HZ.1 = 1/(1+exp(-a1.h))
  
  P1.Y10.ZZ = (pp+(1-pp)*cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1)) + (1-pp)*cc*(1-H1)*(ZZ*HZ.0+(1-ZZ)*(1-HZ.0))
  P0.Y10.ZZ = (1-pp)*(1-cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1))  + (pp+(1-pp)*(1-cc))*(1-H1)*(ZZ*HZ.0+(1-ZZ)*(1-HZ.0))
  
  ff.11.ZZ  = (pp+(1-pp)*cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1))/P1.Y10.ZZ
  ff.10.ZZ  = (1-pp)*(1-cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1))/P0.Y10.ZZ
  
  EE.Y = Y1.0*ff.11.ZZ + (1-Y1.0)*ff.10.ZZ
  
  return(EE.Y)
}

## Three part of Q-function
#
Greenberg.and.Z.reg1 = function(pp,cc,beta.1,E.Y,ZZ,chi.1)
{
  KK = ncol(chi.1)
  b1 = chi.1%*%beta.1
  H1 = 1/(1+exp(-b1))
  
  log.like.GZ = sum(E.Y*log(H1) + (1-E.Y)*log(1-H1))
  
  return(-log.like.GZ)
}

#
Greenberg.and.Z.reg11 = function(pp,cc,alpha.1,E.Y,ZZ,chi.1)
{
  KK = ncol(chi.1)
  a1   = chi.1%*%alpha.1
  HZ.1 = 1/(1+exp(-a1))
  
  log.like.GZ1 = sum(E.Y*ZZ*log(HZ.1) + E.Y*(1-ZZ)*log(1-HZ.1))
  
  return(-log.like.GZ1)
}

#
Greenberg.and.Z.reg10 = function(pp,cc,alpha.0,E.Y,ZZ,chi.1)
{
  KK = ncol(chi.1)
  a0   = chi.1%*%alpha.0
  HZ.0 = 1/(1+exp(-a0))
  
  log.like.GZ0 = sum((1-E.Y)*ZZ*log(HZ.0) + (1-E.Y)*(1-ZZ)*log(1-HZ.0))
  
  return(-log.like.GZ0)
}

#
EM.alogoritmH1 = function(Y1.0,ZZ,chi.1,pp,cc,Theta.hat)
{
  KK = ncol(chi.1)
  KS = 0
  Err= 1
  index.convergence = 1
  
  Est.Theta = rep(-99999,length(Theta.hat))
  
  while(KS <= 100 & Err >= 0.00001)
  {
    E.Y = CEY(Y1.0,ZZ,chi.1,pp,cc, Theta.hat)
    
    beta.hat    = Theta.hat[1:KK]
    alpha.0.hat = Theta.hat[(KK+1):(2*KK)]
    alpha.1.hat = Theta.hat[(2*KK+1):(3*KK)]
    
    beta.1.est  = nlminb(beta.hat,    Greenberg.and.Z.reg1,gr=NULL,  E.Y=E.Y,ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,hessian=TRUE)
    alpha.0.est = nlminb(alpha.0.hat, Greenberg.and.Z.reg10,gr=NULL, E.Y=E.Y,ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,hessian=TRUE)
    alpha.1.est = nlminb(alpha.1.hat, Greenberg.and.Z.reg11,gr=NULL, E.Y=E.Y,ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,hessian=TRUE)
    
    Theta.new.hat = c(beta.1.est$par, alpha.0.est$par, alpha.1.est$par)
    
    Err = sum(abs(Theta.hat - Theta.new.hat))/9
    
    if(Err<=0.00001) {Est.Theta = Theta.hat; index.convergence=0 }
    Theta.hat = Theta.new.hat
    
    KS = KS+1
  }
  return(c(Est.Theta,Err,index.convergence))
}

##
CP.95 = function(est1,est1.se,Theta.t)
{
  NP     = length(Theta.t)
  A.025  = matrix(0,nrow(est1),NP)
  A.975  = matrix(0,nrow(est1),NP)
  A.95CP = matrix(0,nrow(est1),NP)
  
  for(k in 1:nrow(est1))
  {
    A.025[k,] = est1[k,] - 1.96*est1.se[k,]
    A.975[k,] = est1[k,] + 1.96*est1.se[k,]
    
    for(j in 1:NP)
    {
      if(Theta.t[j]>=A.025[k,j] & Theta.t[j]<=A.975[k,j]) {A.95CP[k,j]=1}
    } #end for j
  } ## end for k
  
  return(A.95CP)
}

## Orginal Log likelihood for (Y1.0,ZZ)
Log.like.GRY0Z = function(pp,cc,Y1.0,ZZ,chi.1,Theta.1)
{
  KK = ncol(chi.1)
  
  beta.h    = Theta.1[1:KK]
  alpha.00h = Theta.1[(KK+1):(3*KK)]
  alpha.0h  = alpha.00h[1:KK]
  alpha.1h  = alpha.00h[(KK+1):(2*KK)]
  
  b1.h = chi.1%*%beta.h
  a0.h = chi.1%*%alpha.0h
  a1.h = chi.1%*%alpha.1h
  
  H1   = 1/(1+exp(-b1.h))
  HZ.0 = 1/(1+exp(-a0.h))
  HZ.1 = 1/(1+exp(-a1.h))
  
  P1.Y10.ZZ = (pp+(1-pp)*cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1)) + (1-pp)*cc*(1-H1)*(ZZ*HZ.0+(1-ZZ)*(1-HZ.0))
  P0.Y10.ZZ = (1-pp)*(1-cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1))  + (pp+(1-pp)*(1-cc))*(1-H1)*(ZZ*HZ.0+(1-ZZ)*(1-HZ.0))
  
  log.like.Y0Z = sum(Y1.0*log(P1.Y10.ZZ) + (1-Y1.0)*log(P0.Y10.ZZ))
  
  return(log.like.Y0Z)
}

# Partial derivative of Q in Theta.hat part,  to let Q function be a vectors
Df1.LogLike.Y0Z = function(pp,cc,Y1.0,ZZ,chi.1, Theta.hat)
{
  delta = 0.00005
  k    = length(Theta.hat)
  d1   = diag(rep(delta,k))                        
  df.1 = matrix(0,1,k)                            
  
  for(i in 1:k)                                 
  {       
    LL1 = Log.like.GRY0Z(pp, cc, Y1.0, ZZ, chi.1, Theta.hat+d1[i,])
    LL2 = Log.like.GRY0Z(pp, cc, Y1.0, ZZ, chi.1, Theta.hat-d1[i,])
    
    df.1[1,i] = (LL1-LL2)/(2*delta) 
  }       
  return(t(df.1)) 
}
##
Df2.LogLike.Y0Z = function(pp,cc,Y1.0,ZZ,chi.1, Theta.hat)
{
  delta = 0.00005
  
  k    =length(Theta.hat)                               
  d1   = diag(rep(delta,k))                        
  df.1 = matrix(0,k,k)                            
  
  for(i in 1:k)                                 
  {       
    DF1 = Df1.LogLike.Y0Z(pp, cc, Y1.0, ZZ, chi.1, Theta.hat+d1[i,])
    DF2 = Df1.LogLike.Y0Z(pp, cc, Y1.0, ZZ, chi.1, Theta.hat-d1[i,])
    
    df.1[,i] = (DF1-DF2)/(2*delta)
  }       
  
  return(df.1) 
}

### Log.like.Y0
Log.like.Y0.H0 = function(pp,cc,Y1.0,ZZ,chi.1,beta.H0)
{
  KK = ncol(chi.1)
  
  b1 = chi.1%*%beta.H0
  H1 = 1/(1+exp(-b1))
  
  P.Y0.1 = (pp+(1-pp)*cc)*H1+(1-pp)*cc*(1-H1)
  P.Y0.0 = (1-pp)*(1-cc)*H1+(pp+(1-pp)*(1-cc))*(1-H1)
  
  Log.like.Y0 = sum(Y1.0*log(P.Y0.1)+(1-Y1.0)*log(P.Y0.0))
  
  return(-Log.like.Y0)
}


## Log.like.ZZ
Log.like.ZZ.H0 = function(pp,cc,Y1.0,ZZ,chi.1,alpha.H0)
{
  KK = ncol(chi.1)
  
  a0 = chi.1%*%alpha.H0
  HZ = 1/(1+exp(-a0))
  
  Log.like.ZZ = sum(ZZ*log(HZ)+(1-ZZ)*log(1-HZ))
  
  return(-Log.like.ZZ)
}

##
Estimation.H0 = function(pp,cc,Y1.0,ZZ,chi.1,Theta.H0)
{
  KK = ncol(chi.1)
  
  beta.H0  = Theta.H0[1:KK]
  alpha.H0 = Theta.H0[(KK+1):(2*KK)]
  
  beta.1.H0  = nlminb(beta.H0,  Log.like.Y0.H0,gr=NULL,Y1.0=Y1.0, ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,hessian=TRUE)
  alpha.0.H0 = nlminb(alpha.H0, Log.like.ZZ.H0,gr=NULL, Y1.0=Y1.0,ZZ=ZZ,chi.1=chi.1, pp=pp, cc=cc,hessian=TRUE)
  
  Theta.hat.H0=c(beta.1.H0$par,alpha.0.H0$par)
  
  return(Theta.hat.H0)
}


#####
# Greenberg Part Only
Log.like.GR = function(GZ.data, beta.t)
{
  # Theta.1 = Theta.t
  GZ.data = as.data.frame(GZ.data)
  Y1.0 = GZ.data$Y1.0
  X1 = GZ.data$X1
  X2 = GZ.data$X2
  chi.1   = cbind(1,X1,X2)
  
  b1.h = as.matrix(chi.1)%*%beta.t
  H    = 1/(1+exp(-b1.h))

  P1.Y10 = pp*H + cc*(1-pp)
  P0.Y10 = 1 - pp*H - cc*(1-pp)
    
  log.like = sum(Y1.0*log(P1.Y10) + (1-Y1.0)*log(P0.Y10))
  
  return(log.like)
}

M.GR = function(GZ.data, beta.t)  # It is J_n part
{
  GZ.data = as.data.frame(GZ.data)
  Y1.0 = GZ.data$Y1.0
  X1   = GZ.data$X1
  X2   = GZ.data$X2
  chi.1   = cbind(1,X1,X2)
  
  b1.h = as.matrix(chi.1)%*%beta.t
  H    = c(1/(1+exp(-b1.h)))
  H1   = c(H*(1-H))
  
  AA = pp*H1/((pp*H+(1-pp)*cc)*(1-pp*H-(1-pp)*cc))
  
  SS = chi.1*AA*(Y1.0-pp*H-pp*(1-cc))
  
  MM = t(SS)%*%SS
  
  return(MM)
}


################################
pp = 0.7
cc = 0.5
nn = 2000


beta.t = c(-1.0,-1,1.2)
alpha.t0 = c( -1.0,-0.5,1.5)
alpha.t  = c( -1.0,-0.5,1.5,-1.2, 1.2,1.2)
Theta.t0 = c(beta.t,alpha.t0)
Theta.t  = c(beta.t,alpha.t)

KT = 1000
est1    = matrix(0,KT,length(Theta.t))
est1.se = matrix(0,KT,length(Theta.t))
est1.cp = matrix(0,KT,length(Theta.t))

est0    = matrix(0,KT,2*length(beta.t))
est0.se = matrix(0,KT,2*length(beta.t))
est0.cp = matrix(0,KT,2*length(beta.t))

est.GR  = matrix(0,KT,length(beta.t))
est.GR.se = matrix(0,KT,length(beta.t))
est.GR.cp = matrix(0,KT,length(beta.t))

LR.test = matrix(0,KT,2)

jj  = 1
KTR = 1000+300
KTR.0 = 0


while(jj <= KT & KTR.0 <= (KT + 3000))
{
  KTR.0     = KTR.0 + 1
  cat(jj)
  
  GZ.data = Greenberg.RR.data(nn,pp,cc,beta.t, alpha.t)
  
  Y1.0  = GZ.data[,2]
  ZZ    = GZ.data[,3]
  chi.1 = GZ.data[,4:ncol(GZ.data)]
  
  b1 = chi.1%*%beta.t
  a0 = chi.1%*%alpha.t[c(1:3)]
  a1 = chi.1%*%alpha.t[c(4:6)]
  
  H1   = 1/(1+exp(-b1))
  HZ.0 = 1/(1+exp(-a0))
  HZ.1 = 1/(1+exp(-a1))
  
  mean(H1)
  mean(HZ.0)
  mean(HZ.1)
  
  KK = ncol(chi.1)
  Est.Theta.all = EM.alogoritmH1(Y1.0,ZZ,chi.1,pp,cc, Theta.t)
  
  if(Est.Theta.all[3*KK+1]>0.00001 | Est.Theta.all[3*KK+2]!=0) 
  {
    cat("Error", "\n")  
    next
  }
  
  est1[jj,]=Est.Theta.all[1:(3*KK)]
  
  Theta.H = Est.Theta.all[1:(3*KK)]
  
  AA.11        = Df2.LogLike.Y0Z(pp,cc,Y1.0,ZZ,chi.1, Theta.H)
  est1.se[jj,] = sqrt(diag(solve(-AA.11)))
  
  Theta.H0     = c(beta.t,alpha.t[1:3])
  Est.Theta.H0 = Estimation.H0(pp,cc,Y1.0,ZZ,chi.1,Theta.H0)
  est0[jj,]    = Est.Theta.H0
  
  Theta.H0.1   = c(Est.Theta.H0,Est.Theta.H0[(KK+1):(2*KK)])
  
  LR.test[jj,1]=2*(Log.like.GRY0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.H)-Log.like.GRY0Z(pp,cc,Y1.0,ZZ,chi.1,Theta.H0.1))
  LR.test[jj,2]=1-pchisq(LR.test[jj,1],KK)
  
  Theta.est.GR = optim(beta.t, Log.like.GR, GZ.data = GZ.data, gr = NULL, method="BFGS", control=list(fnscale=-1), hessian=TRUE)
  
  if (Theta.est.GR$convergence != 0){
    cat("Possible error in GR:", "\n")
    next
  }
  
  if (min(abs(diag(Theta.est.GR$hessian)))==0){
    cat("Hessian matrix in GR: diag=0", "\n") 
    next
  }
  
  Inv.Hess.GR = solve(-Theta.est.GR$hessian)
  var.GR  = diag(Inv.Hess.GR)
  
  if (min(var.GR) <= 0){
    cat("min(var.GR) <= 0", "\n")
    next
  }
  
  est.GR[jj,]     = Theta.est.GR$par
  M.theta.GR      = M.GR(GZ.data, Theta.est.GR$par) 
  est.GR.se[jj,]  = sqrt(var.GR)
  
  
  jj = jj + 1
}


est1.cp   = CP.95(est1, est1.se, Theta.t)
est.GR.cp = CP.95(est.GR, est.GR.se, beta.t)

output.1 = cbind(Theta.t,round(apply(est1,2,mean),4),round(apply(est1,2,sd),4),round(apply(est1.se,2,mean),4),round(apply(est1.cp,2,mean),4))
colnames(output.1) = c("True par.","EM-estimaor","SD-simul.","ASE","95%CP")
rownames(output.1) = c("beta_0","beta_1","beta_2","alpha0_0","alpha0_1","alpha0_2","alpha1_0","alpha1_1","alpha1_2")
output.1

output.2 = cbind(round(apply(est0,2,mean),4),round(apply(est0,2,sd),4))
colnames(output.2) = c("EM-estimaor","SD-simul.")
rownames(output.2) = c("beta_0","beta_1","beta_2","alpha0_0","alpha0_1","alpha0_2")
output.2

output.3 = cbind(beta.t,round(apply(est.GR,2,mean),4),round(apply(est.GR,2,sd),4),round(apply(est.GR.se,2,mean),4),round(apply(est.GR.cp,2,mean),4))
colnames(output.3) = c("True par.","ML-estimaor","SD-simul.","ASE","95%CP")
rownames(output.3) = c("beta_0","beta_1","beta_2")
output.3


length(LR.test[LR.test[,2]<=0.05,2])/1000  #Power
pp; cc; nn

