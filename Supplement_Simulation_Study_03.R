############################################################
# Simulation Script (Case 3 / Binary X1, X2): 
# Conditional distributional framework for joint modeling of 
# sensitive attribute and an observed Variable under randomized response designs
#
# This script reproduces a simulation setting in Section 4 
# and the comparison with the Greenberg-only model reported 
# in the manuscript.
#
# Model:
#   P(Y, Z | X) = P(Y | X) P(Z | Y, X)
#   P(Y=1 | X)      = H(beta^T X)
#   P(Z=1 | Y=y, X) = H(alpha_y^T X), y=0,1
# where H(t) is the logistic link.
#
# RRT mechanism (unrelated-question design):
#   With probability p, respondent answers the sensitive question (Y);
#   otherwise answers an unrelated question that yields "Yes" with prob c.
# The observed randomized response is Y* (denoted Y1.0 in the code).
#
# Goals of the simulation:
#   - Evaluate EM estimator: Bias, SD, ASE, 95% CP
#   - Likelihood ratio test (LRT) for H0: alpha0 = alpha1
#   - Compare beta estimates with the Greenberg-only model
#     (which uses only the RRT part and ignores Z).
#
# Script structure:
#   (1) Data generation under the joint model + RRT
#   (2) EM algorithm for the full joint model (H1)
#   (3) Restricted estimation under H0: alpha0 = alpha1
#   (4) Greenberg-only estimation (baseline comparison)
#   (5) Monte Carlo loop (KT replications)
#   (6) Output summary tables + LRT empirical rejection rate
############################################################

require(stats)
library(xtable)
set.seed(12345)

############################################################
# Greenberg.RR.data(nn, pp, cc, beta.1, alpha.00)
#
# Generates one synthetic dataset under the joint model and
# unrelated-question RRT mechanism.
#
# Inputs:
#   nn       : sample size
#   pp (p)   : probability of answering sensitive question
#   cc (c)   : probability of forced "Yes" under unrelated question
#   beta.1   : true beta parameters for P(Y=1|X)
#   alpha.00 : stacked true alpha parameters:
#              alpha.0 = alpha.00[1:3] for P(Z=1|Y=0,X)
#              alpha.1 = alpha.00[4:6] for P(Z=1|Y=1,X)
#
# Data generating steps:
#   1) Generate covariates: X1, X2 ~ Bernoulli(0.5)
#   2) Generate latent Y ~ Bernoulli(H(beta^T X))
#   3) Generate Z conditional on Y using alpha0/alpha1
#   4) Apply RRT to obtain observed Y* (Y1.0)
#
# Output columns:
#   Y1   : true latent sensitive outcome (unobserved in practice)
#   Y1.0 : observed randomized response Y* under RRT
#   ZZ   : observed binary variable Z
#   chi.1: design matrix (Intercept, X1, X2)
############################################################


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
############################################################
# CEY(Y1.0, ZZ, chi.1, pp, cc, Theta.hat)
#
# E-step of EM algorithm:
# Computes posterior expectation E(Y | Y*, Z, X; Theta_hat).
#
# This is the key step that links the observed randomized response (Y*)
# and observed Z to the latent sensitive variable Y.
############################################################


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
############################################################
# Q-function components (M-step objectives)
#
# Given E(Y | data), we maximize three separate weighted logistic
# log-likelihoods corresponding to:
#   - beta   : model for P(Y=1|X)
#   - alpha1 : model for P(Z=1|Y=1,X)
#   - alpha0 : model for P(Z=1|Y=0,X)
#
# Each function returns the negative log-likelihood (for nlminb).
############################################################


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

############################################################
# EM.alogoritmH1(Y1.0, ZZ, chi.1, pp, cc, Theta.hat)
#
# EM algorithm for the full joint model (H1).
#
# Iteration:
#   - E-step: compute E.Y = E(Y | Y*, Z, X; Theta_hat)
#   - M-step: update beta, alpha0, alpha1 via nlminb
#
# Convergence:
#   Err = average absolute parameter change (scaled)
#   stop when Err < 1e-5 or KS > 100
#
# Returns:
#   Estimated parameters (beta, alpha0, alpha1),
#   final Err, and convergence flag (0 indicates success).
############################################################

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

############################################################
# CP.95(est1, est1.se, Theta.t)
#
# Computes empirical 95% coverage probability for each parameter:
#   CI_k = est ± 1.96 * SE
# Coverage = proportion of replications whose CI contains truth.
############################################################

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

############################################################
# Observed log-likelihood for (Y*, Z) under the full joint model
# and numerical derivatives (gradient/Hessian).
#
# Log.like.GRY0Z:
#   Used to compute LRT:
#     LR = 2[ logL(H1) - logL(H0) ]
#
# Df1.LogLike.Y0Z and Df2.LogLike.Y0Z:
#   Numerical gradient and Hessian for logL.
#   The inverse observed information approximates Var(Theta_hat):
#     Var(Theta_hat) ≈ (-Hessian)^{-1}
############################################################

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

############################################################
# Estimation under H0: alpha0 = alpha1
#
# Under H0, Z does not depend on Y after conditioning on X:
#   P(Z=1 | Y=0,X) = P(Z=1 | Y=1,X) = H(alpha^T X)
#
# Estimation.H0 returns MLE under the restricted model, used in LRT.
############################################################

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
############################################################
# Greenberg-only baseline (ignoring Z)
#
# This block implements the model that uses only (Y*, X) under
# Greenberg's unrelated-question design:
#   P(Y* | X; beta)
#
# Log.like.GR:
#   Log-likelihood for beta based only on randomized response Y*.
#
# M.GR:
#   Computes the empirical outer product of scores (J_n / M_n-type term),
#   used for variance-related calculations in the Greenberg-only approach.
#
# This provides a baseline comparison for beta estimation efficiency.
############################################################

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
############################################################
# SIMULATION SETTINGS (one scenario)
#
# pp (p)   : probability of answering sensitive question
# cc (c)   : forced "Yes" probability under unrelated question
# nn       : sample size
#
# True parameters:
#   beta.t   : parameters for P(Y=1|X)
#   alpha.t  : stacked alpha0 and alpha1 for P(Z=1|Y,X)
#
# Monte Carlo:
#   KT       : target number of successful replications
#   KTR.0    : total attempts (allows skipping non-converged runs)
############################################################

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

############################################################
# MONTE CARLO LOOP
#
# For each replication:
#   (1) Generate data: (Y, Y*, Z, X)
#   (2) Fit full joint model via EM (H1)
#       - If EM fails to converge, skip and try next replication
#   (3) Compute ASE from inverse observed information (-Hessian)^{-1}
#   (4) Fit restricted model under H0: alpha0 = alpha1
#   (5) Compute LRT p-value
#   (6) Fit Greenberg-only baseline for beta via optim (BFGS)
#       - skip if convergence/Hessian issues occur
#
# Stored objects:
#   est1, est1.se: EM estimates and ASEs under H1
#   est0         : restricted estimates under H0
#   est.GR       : Greenberg-only beta estimates and ASEs
#   LR.test      : LRT statistic and p-value
############################################################

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

############################################################
# OUTPUT TABLES
#
# output.1 (Full joint model, EM under H1):
#   - True par.   : true parameter values
#   - EM-estimaor : Monte Carlo mean of EM estimates
#   - SD-simul.   : empirical SD across replications
#   - ASE         : mean asymptotic SE (from observed information)
#   - 95%CP       : empirical 95% coverage probability
#
# output.2 (Restricted model under H0):
#   - EM-estimaor, SD-simul. for (beta, alpha) under H0
#
# output.3 (Greenberg-only baseline):
#   - ML-estimaor : Monte Carlo mean of beta MLE from Greenberg-only model
#   - SD/ASE/CP   : corresponding performance metrics
#
# LRT summary:
#   length(LR.test[p<=0.05])/KT gives empirical rejection rate 
#   (power under H1, type I error under H0).
############################################################


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


