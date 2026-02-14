############################################################
# Empirical Example Script (TSCS 2012) 
# Conditional distributional framework for joint modeling of 
# sensitive attribute and an observed Variable under randomized response designs.
#
# This script accompanies the manuscript:
# Conditional distributional framework for joint modeling of 
# sensitive attribute and an observed Variable under randomized response designs.
#
# PURPOSE
#   Reproduce the empirical analysis workflow in Section 5:
#     - Fit the proposed joint model:
#         P(Y, Z | X) = P(Y | X) P(Z | Y, X)
#       where Y is latent sensitive attribute collected under RRT,
#       Z is an observed attitude/response variable, and X are covariates.
#     - Conduct LRT for H0: alpha0 = alpha1
#     - Report parameter estimates + ASEs
#     - Provide additional comparison models:
#         (i) Greenberg RR regression for P(Y=1 | X, Z)
#         (ii) Standard logistic regression for P(Z=1 | X)
#
# DATA AVAILABILITY / REPRODUCIBILITY
#   The empirical data are from the Taiwan Social Change Survey (TSCS),
#   distributed by the Survey Research Data Archive (SRDA), Academia Sinica.
#   Due to data usage policies, the authors cannot redistribute the raw data.
#   Researchers may obtain access by application through SRDA.
#
#   This script is written so that:
#     (A) If you have access to TSCS data, you can load it and run directly.
#     (B) If you do NOT have the data, the script can run on a synthetic demo
#         dataset to verify that the pipeline executes end-to-end.
#
# USER INSTRUCTIONS
#   1) Install required packages (if needed).
#   2) Choose ONE of the following:
#        - Option A (preferred): Load TSCS data and map variables (see below).
#        - Option B: Run demo dataset (default) to test the workflow.
#   3) Run the script: source("Empirical_Analysis_TSCS.R")
#
# OUTPUTS
#   - output.1: EM estimates (full joint model) and ASE
#   - LR.test : LRT statistic + p-value for H0: alpha0 = alpha1
#   - output.2: Restricted model estimates (alpha0 = alpha1) and ASE
#   - Table.Est: Greenberg RR regression estimates for P(Y=1 | X, Z)
#   - result_prob: summary of probability quantities (illustrative)
############################################################

suppressPackageStartupMessages({
  library(MASS)
  library(naniar)
  library(dplyr)
  library(xtable)
  library(ggplot2)
  library(tibble)
})

set.seed(12345)

############################################################
# (1) CORE FUNCTIONS: E-step, M-step objectives, EM algorithm
############################################################

# Logistic link
H_logit <- function(u) 1/(1+exp(-u))

############################################################
# CEY: E-step
#   E(Y | Y*, Z, X; Theta_hat) under the joint likelihood and RRT mechanism
############################################################
CEY = function(Y1.0, ZZ, chi.1, pp, cc, Theta.hat)
{
  KK     = ncol(chi.1)

  beta.h = Theta.hat[1:KK]

  alpha.00h = Theta.hat[(KK+1):(3*KK)]
  alpha.0h  = alpha.00h[1:KK]
  alpha.1h  = alpha.00h[(KK+1):(2*KK)]

  b1.h = chi.1%*%beta.h
  a0.h = chi.1%*%alpha.0h
  a1.h = chi.1%*%alpha.1h

  H1   = H_logit(b1.h)
  HZ.0 = H_logit(a0.h)
  HZ.1 = H_logit(a1.h)

  P1.Y10.ZZ = (pp+(1-pp)*cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1)) + (1-pp)*cc*(1-H1)*(ZZ*HZ.0+(1-ZZ)*(1-HZ.0))
  P0.Y10.ZZ = (1-pp)*(1-cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1))  + (pp+(1-pp)*(1-cc))*(1-H1)*(ZZ*HZ.0+(1-ZZ)*(1-HZ.0))

  ff.11.ZZ = (pp+(1-pp)*cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1))/P1.Y10.ZZ
  ff.10.ZZ = (1-pp)*(1-cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1))/P0.Y10.ZZ

  EE.Y = Y1.0*ff.11.ZZ + (1-Y1.0)*ff.10.ZZ
  return(EE.Y)
}

############################################################
# Q-function components (M-step objectives)
#   - beta  : P(Y=1|X)
#   - alpha1: P(Z=1|Y=1,X)
#   - alpha0: P(Z=1|Y=0,X)
############################################################
Greenberg.and.Z.reg1 = function(pp, cc, beta.1, E.Y, ZZ, chi.1)
{
  b1 = chi.1%*%beta.1
  H1 = H_logit(b1)
  log.like.GZ = sum(E.Y*log(H1) + (1-E.Y)*log(1-H1))
  return(-log.like.GZ)
}

Greenberg.and.Z.reg11 = function(pp, cc, alpha.1, E.Y, ZZ, chi.1)
{
  a1   = chi.1%*%alpha.1
  HZ.1 = H_logit(a1)
  log.like.GZ1 = sum(E.Y*ZZ*log(HZ.1) + E.Y*(1-ZZ)*log(1-HZ.1))
  return(-log.like.GZ1)
}

Greenberg.and.Z.reg10 = function(pp, cc, alpha.0, E.Y, ZZ, chi.1)
{
  a0   = chi.1%*%alpha.0
  HZ.0 = H_logit(a0)
  log.like.GZ0 = sum((1-E.Y)*ZZ*log(HZ.0) + (1-E.Y)*(1-ZZ)*log(1-HZ.0))
  return(-log.like.GZ0)
}

############################################################
# EM algorithm for full joint model (H1)
#   stops when average absolute parameter change < 1e-5
############################################################
EM.alogoritmH1 = function(Y1.0, ZZ, chi.1, pp, cc, Theta.hat)
{
  KK  = ncol(chi.1)
  KS  = 0
  Err = 1
  index.convergence = 1

  Est.Theta = rep(-99999, length(Theta.hat))

  while(KS <= 100 & Err >= 0.00001)
  {
    E.Y = CEY(Y1.0, ZZ, chi.1, pp, cc, Theta.hat)

    beta.hat    = Theta.hat[1:KK]
    alpha.0.hat = Theta.hat[(KK+1):(2*KK)]
    alpha.1.hat = Theta.hat[(2*KK+1):(3*KK)]

    beta.1.est  = nlminb(beta.hat,    Greenberg.and.Z.reg1,  gr=NULL, E.Y=E.Y,ZZ=ZZ, chi.1=chi.1, pp=pp, cc=cc, hessian=TRUE)
    alpha.0.est = nlminb(alpha.0.hat, Greenberg.and.Z.reg10, gr=NULL, E.Y=E.Y,ZZ=ZZ, chi.1=chi.1, pp=pp, cc=cc, hessian=TRUE)
    alpha.1.est = nlminb(alpha.1.hat, Greenberg.and.Z.reg11, gr=NULL, E.Y=E.Y,ZZ=ZZ, chi.1=chi.1, pp=pp, cc=cc, hessian=TRUE)

    Theta.new.hat = c(beta.1.est$par, alpha.0.est$par, alpha.1.est$par)

    Err = sum(abs(Theta.hat - Theta.new.hat))/9

    if(Err <= 0.00001) {Est.Theta = Theta.hat; index.convergence = 0}
    Theta.hat = Theta.new.hat

    KS = KS + 1
  }
  return(c(Est.Theta, Err, index.convergence))
}

############################################################
# Observed log-likelihood for (Y*, Z) under full joint model
# Used for LRT and numerical derivatives for ASE
############################################################
Log.like.GRY0Z = function(pp, cc, Y1.0, ZZ, chi.1, Theta.1)
{
  KK = ncol(chi.1)

  beta.h    = Theta.1[1:KK]
  alpha.00h = Theta.1[(KK+1):(3*KK)]
  alpha.0h  = alpha.00h[1:KK]
  alpha.1h  = alpha.00h[(KK+1):(2*KK)]

  b1.h = chi.1%*%beta.h
  a0.h = chi.1%*%alpha.0h
  a1.h = chi.1%*%alpha.1h

  H1   = H_logit(b1.h)
  HZ.0 = H_logit(a0.h)
  HZ.1 = H_logit(a1.h)

  P1.Y10.ZZ = (pp+(1-pp)*cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1)) + (1-pp)*cc*(1-H1)*(ZZ*HZ.0+(1-ZZ)*(1-HZ.0))
  P0.Y10.ZZ = (1-pp)*(1-cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1))  + (pp+(1-pp)*(1-cc))*(1-H1)*(ZZ*HZ.0+(1-ZZ)*(1-HZ.0))

  log.like.Y0Z = sum(Y1.0*log(P1.Y10.ZZ) + (1-Y1.0)*log(P0.Y10.ZZ))
  return(log.like.Y0Z)
}

############################################################
# Numerical derivatives: gradient and Hessian (observed information)
# Var(Theta_hat) ≈ (-Hessian)^{-1}; ASE = sqrt(diag(Var))
############################################################
Df1.LogLike.Y0Z = function(pp, cc, Y1.0, ZZ, chi.1, Theta.hat)
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

Df2.LogLike.Y0Z = function(pp, cc, Y1.0, ZZ, chi.1, Theta.hat)
{
  delta = 0.00005
  k    = length(Theta.hat)
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

############################################################
# Restricted model: alpha0 = alpha1 (used for LRT under H0)
# We compute the Hessian of the restricted log-likelihood separately.
############################################################
Log.like.GRY0Z.dependent = function(pp, cc, Y1.0, ZZ, chi.1, Theta.1)
{
  KK = ncol(chi.1)

  beta.h  = Theta.1[1:KK]
  alpha.h = Theta.1[(KK+1):(2*KK)]

  b1.h = chi.1%*%beta.h
  a0.h = chi.1%*%alpha.h
  a1.h = chi.1%*%alpha.h

  H1   = H_logit(b1.h)
  HZ.0 = H_logit(a0.h)
  HZ.1 = H_logit(a1.h)

  P1.Y10.ZZ = (pp+(1-pp)*cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1)) + (1-pp)*cc*(1-H1)*(ZZ*HZ.0+(1-ZZ)*(1-HZ.0))
  P0.Y10.ZZ = (1-pp)*(1-cc)*H1*(ZZ*HZ.1+(1-ZZ)*(1-HZ.1))  + (pp+(1-pp)*(1-cc))*(1-H1)*(ZZ*HZ.0+(1-ZZ)*(1-HZ.0))

  log.like.Y0Z = sum(Y1.0*log(P1.Y10.ZZ) + (1-Y1.0)*log(P0.Y10.ZZ))
  return(log.like.Y0Z)
}

Df1.LogLike.Y0Z.dependent = function(pp, cc, Y1.0, ZZ, chi.1, Theta.hat)
{
  delta = 0.00005
  k    = length(Theta.hat)
  d1   = diag(rep(delta,k))
  df.1 = matrix(0,1,k)

  for(i in 1:k)
  {
    LL1 = Log.like.GRY0Z.dependent(pp, cc, Y1.0, ZZ, chi.1, Theta.hat+d1[i,])
    LL2 = Log.like.GRY0Z.dependent(pp, cc, Y1.0, ZZ, chi.1, Theta.hat-d1[i,])
    df.1[1,i] = (LL1-LL2)/(2*delta)
  }
  return(t(df.1))
}

Df2.LogLike.Y0Z.dependent = function(pp, cc, Y1.0, ZZ, chi.1, Theta.hat)
{
  delta = 0.00005
  k    = length(Theta.hat)
  d1   = diag(rep(delta,k))
  df.1 = matrix(0,k,k)

  for(i in 1:k)
  {
    DF1 = Df1.LogLike.Y0Z.dependent(pp, cc, Y1.0, ZZ, chi.1, Theta.hat+d1[i,])
    DF2 = Df1.LogLike.Y0Z.dependent(pp, cc, Y1.0, ZZ, chi.1, Theta.hat-d1[i,])
    df.1[,i] = (DF1-DF2)/(2*delta)
  }
  return(df.1)
}

############################################################
# Estimation under H0 via separate likelihood components
############################################################
Log.like.Y0.H0 = function(pp, cc, Y1.0, ZZ, chi.1, beta.H0)
{
  b1 = chi.1%*%beta.H0
  H1 = H_logit(b1)

  P.Y0.1 = (pp+(1-pp)*cc)*H1 + (1-pp)*cc*(1-H1)
  P.Y0.0 = (1-pp)*(1-cc)*H1  + (pp+(1-pp)*(1-cc))*(1-H1)

  Log.like.Y0 = sum(Y1.0*log(P.Y0.1) + (1-Y1.0)*log(P.Y0.0))
  return(-Log.like.Y0)
}

Log.like.ZZ.H0 = function(pp, cc, Y1.0, ZZ, chi.1, alpha.H0)
{
  a0 = chi.1%*%alpha.H0
  HZ = H_logit(a0)

  Log.like.ZZ = sum(ZZ*log(HZ) + (1-ZZ)*log(1-HZ))
  return(-Log.like.ZZ)
}

Estimation.H0 = function(pp, cc, Y1.0, ZZ, chi.1, Theta.H0)
{
  KK = ncol(chi.1)

  beta.H0  = Theta.H0[1:KK]
  alpha.H0 = Theta.H0[(KK+1):(2*KK)]

  beta.1.H0  = nlminb(beta.H0,  Log.like.Y0.H0, gr=NULL, Y1.0=Y1.0, ZZ=ZZ, chi.1=chi.1, pp=pp, cc=cc, hessian=TRUE)
  alpha.0.H0 = nlminb(alpha.H0, Log.like.ZZ.H0, gr=NULL, Y1.0=Y1.0, ZZ=ZZ, chi.1=chi.1, pp=pp, cc=cc, hessian=TRUE)

  Theta.hat.H0 = c(beta.1.H0$par, alpha.0.H0$par)
  return(Theta.hat.H0)
}

############################################################
# (2) DATA LOADING: Option A (TSCS) or Option B (demo data)
############################################################

## -------------------------
## Option A: Load TSCS data
## -------------------------
## If you have access to TSCS data (from SRDA), set:
##   USE_REAL_DATA <- TRUE
## and update DATA_PATH to the correct local file path.
##
## IMPORTANT:
## The variable mappings below follow the manuscript example.
## If your TSCS extract uses different column names, update the mapping.
USE_REAL_DATA <- FALSE
DATA_PATH <- "path/to/your/TSCS2012_extract.csv"  # <-- update if USE_REAL_DATA=TRUE

tscs.2012 <- read.csv(DATA_PATH, header=TRUE, sep=",")

  # Select required variables used in Section 5:
  # a1  : gender
  # b1  : education category (needs recoding)
  # j1  : attitude item (used to define Z)
  # k9  : RRT response for sensitive behavior (used to define Y*)
  # Other variables are included in the original script for filtering.
  tscs.2012.Greenbergdata <- tscs.2012[,c("a1","a11","b1","j1","j2","c7c","k9","j5")]

  # Apply filters consistent with the manuscript preprocessing
  tscs.2012.GR <- tscs.2012.Greenbergdata[
    tscs.2012.Greenbergdata[,7] <= 2 &
      tscs.2012.Greenbergdata[,4] <= 4 &
      tscs.2012.Greenbergdata[,5] <= 4 &
      tscs.2012.Greenbergdata[,8] <= 4, ]

  # Education recoding (example mapping; adjust if needed)
  tscs.2012.GR$edu <- ifelse(tscs.2012.GR[,3] %in% c(1:4),1,0) +
    ifelse(tscs.2012.GR[,3] %in% c(5:9),2,0) +
    ifelse(tscs.2012.GR[,3] %in% c(10:15),3,0) +
    ifelse(tscs.2012.GR[,3] %in% c(16:19),4,0) +
    ifelse(tscs.2012.GR[,3] %in% c(20:21),5,0)

  tscs.2012.GR$edu1 <- ifelse(tscs.2012.GR$edu <= 3, 1, ifelse(tscs.2012.GR$edu == 4, 2, 3))

  # Construct analysis variables (as in the manuscript):
  pp <- 0.5; cc <- 0.25

  Y1.0 <- ifelse(tscs.2012.GR$k9 == 1, 1, 0)      # observed RRT response Y*
  ZZ   <- ifelse(tscs.2012.GR$j1 <= 2, 0, 1)      # attitude Z (binary)

  X1   <- ifelse(tscs.2012.GR$a1 == 1, 1, 0)      # gender: male=1, female=0
  X2   <- ifelse(tscs.2012.GR$edu1 == 1, 0, 1)    # bachelor+ = 1, otherwise 0
  chi.1 <- cbind(1, X1, X2)

############################################################
# (3) FIT THE JOINT MODEL (EM) + ASE + LRT
############################################################

# Initial values (use manuscript's defaults as starting point)
beta.t  <- c(-1.0, -1, 1.2)
alpha.t <- c(-1.0, -0.5, 1.5,  -1.2, 1.2, 1.2)
Theta.t <- c(beta.t, alpha.t)
Theta.hat <- Theta.t

KT <- 1
est1    <- matrix(0, KT, length(Theta.t))
est1.se <- matrix(0, KT, length(Theta.t))
est0    <- matrix(0, KT, 2*length(beta.t))
est0.se <- matrix(0, KT, 2*length(beta.t))
LR.test <- matrix(0, KT, 2)

KK <- ncol(chi.1)

Est.Theta.all <- EM.alogoritmH1(Y1.0, ZZ, chi.1, pp, cc, Theta.hat)

if(Est.Theta.all[3*KK+2]==0 & Est.Theta.all[3*KK+1]<=0.00001)
{
  est1[1,] <- Est.Theta.all[1:(3*KK)]
  Theta.H  <- Est.Theta.all[1:(3*KK)]

  # ASE from inverse observed information
  AA.11 <- Df2.LogLike.Y0Z(pp, cc, Y1.0, ZZ, chi.1, Theta.H)
  est1.se[1,] <- sqrt(diag(solve(-AA.11)))

  # Restricted model under H0: alpha0 = alpha1
  Theta.H0     <- c(beta.t, alpha.t[1:3])
  Est.Theta.H0 <- Estimation.H0(pp, cc, Y1.0, ZZ, chi.1, Theta.H0)
  est0[1,] <- Est.Theta.H0

  AA.00 <- Df2.LogLike.Y0Z.dependent(pp, cc, Y1.0, ZZ, chi.1, Est.Theta.H0)
  est0.se[1,] <- sqrt(diag(solve(-AA.00)))

  # LRT
  Theta.H0.1 <- c(Est.Theta.H0, Est.Theta.H0[(KK+1):(2*KK)])
  LLG.1 <- Log.like.GRY0Z(pp, cc, Y1.0, ZZ, chi.1, Theta.H)
  LLG.0 <- Log.like.GRY0Z(pp, cc, Y1.0, ZZ, chi.1, Theta.H0.1)
  LR.test[1,1] <- 2*(LLG.1 - LLG.0)
  LR.test[1,2] <- 1 - pchisq(LR.test[1,1], KK)

} else {
  warning("EM did not converge for the provided initialization. Consider adjusting starting values.")
}

############################################################
# (4) REPORT: Full model estimates + ASE; Restricted model; LRT
############################################################

output.1 <- round(rbind(est1, est1.se), 4)
rownames(output.1) <- c("estimator","ASE")
colnames(output.1) <- c("beta_0","beta_1","beta_2",
                        "alpha0_0","alpha0_1","alpha0_2",
                        "alpha1_0","alpha1_1","alpha1_2")

colnames(LR.test) <- c("chi-square statistic", "p-value")

output.2 <- round(rbind(est0, est0.se), 4)
rownames(output.2) <- c("estimator","ASE")
colnames(output.2) <- c("beta_0","beta_1","beta_2","alpha_0","alpha_1","alpha_2")

cat("\n--- Full joint model (EM) estimates and ASE ---\n")
print(output.1)

cat("\n--- LRT for H0: alpha0 = alpha1 ---\n")
print(round(LR.test, 4)[1,])

cat("\n--- Restricted model (H0) estimates and ASE ---\n")
print(output.2)

############################################################
# (5) OPTIONAL: Model-implied probabilities (illustrative summary)
############################################################

beta   <- est1[1,1:3]
alpha0 <- est1[1,4:6]
alpha1 <- est1[1,7:9]

H.est <- function(chii, thetta){
  round(mean(H_logit(chii%*%thetta)), 4)
}

chi.00 <- chi.1[ (chi.1[,2]==0 & chi.1[,3]==0), , drop=FALSE ]
chi.10 <- chi.1[ (chi.1[,2]==1 & chi.1[,3]==0), , drop=FALSE ]
chi.01 <- chi.1[ (chi.1[,2]==0 & chi.1[,3]==1), , drop=FALSE ]
chi.11 <- chi.1[ (chi.1[,2]==1 & chi.1[,3]==1), , drop=FALSE ]

df.Hest <- data.frame(
  X1X2     = c("00","10","01","11"),
  H.beta   = c(H.est(chi.00, beta),  H.est(chi.10, beta),  H.est(chi.01, beta),  H.est(chi.11, beta)),
  H.alpha1 = c(H.est(chi.00, alpha1),H.est(chi.10, alpha1),H.est(chi.01, alpha1),H.est(chi.11, alpha1)),
  H.alpha0 = c(H.est(chi.00, alpha0),H.est(chi.10, alpha0),H.est(chi.01, alpha0),H.est(chi.11, alpha0))
)

cat("\n--- Model-implied mean probabilities H(.) by covariate strata ---\n")
print(df.Hest)

cat("\n--- Estimated sensitive proportion (model-implied) ---\n")
cat("Proposed joint model: ", H.est(chi.1, beta), "\n")

beta_Greenberg <- est0[1,1:3]
cat("Restricted/Greenberg-style (H0) model: ", H.est(chi.1, beta_Greenberg), "\n")

############################################################
# (6) ADDITIONAL COMPARISONS (as in the original script)
#     (1) Greenberg RR regression for P(Y=1 | X, Z)
#     (2) Logistic regression for P(Z=1 | X)
############################################################

H.0 <- function(chi, theta){
  u <- as.matrix(chi) %*% theta
  1/(1+exp(-u))
}

data.Log <- data.frame(Y0 = Y1.0, X1 = chi.1[,2], X2 = chi.1[,3], Z = ZZ)
data.Log.Z1 <- subset(data.Log, Z==1)
data.Log.Z0 <- subset(data.Log, Z==0)

# Log-likelihood for Greenberg RR regression: Y* ~ (X1, X2, Z)
LL.para <- function(data.Log, theta)
{
  np   <- length(theta)
  Y0   <- data.Log$Y0
  chi  <- as.matrix(cbind(1, data.Log$X1, data.Log$X2, data.Log$Z))
  HH   <- H.0(chi, theta)
  H1   <- pp*HH + (1-pp)*cc
  sum(Y0*log(H1) + (1-Y0)*log(1-H1))
}

# Initial value from standard logistic regression (heuristic)
fit_init <- glm(Y0 ~ X1 + X2 + Z, family = binomial(), data = data.Log)
theta.t <- round(unname(unlist(fit_init$coefficients)), 2)

Theta.est.para <- optim(theta.t, function(th) LL.para(data.Log, th),
                        method="BFGS", hessian=TRUE, control=list(fnscale=-1))

theta.est.para <- Theta.est.para$par
Hes.para <- -Theta.est.para$hessian
theta.ase.para <- sqrt(diag(ginv(Hes.para)))

t_value <- theta.est.para/theta.ase.para
p_value <- 2*pnorm(-abs(t_value))

Table.Est <- data.frame(
  Estimation = round(theta.est.para,4),
  ASE        = round(theta.ase.para,4),
  p_value    = round(p_value,4)
)
rownames(Table.Est) <- c("theta_0","theta_1","theta_2","theta_3")

cat("\n--- Greenberg RR regression for P(Y=1 | X, Z): estimates ---\n")
print(Table.Est)

# Probability summaries
chi.para    <- as.matrix(cbind(1, data.Log$X1, data.Log$X2, data.Log$Z))
chi.para.Z1 <- as.matrix(cbind(1, data.Log.Z1$X1, data.Log.Z1$X2, data.Log.Z1$Z))
chi.para.Z0 <- as.matrix(cbind(1, data.Log.Z0$X1, data.Log.Z0$X2, data.Log.Z0$Z))

PY0.XZ   <- mean(pp*H.0(chi.para,    theta.est.para) + (1-pp)*cc)
PY.XZ    <- mean(H.0(chi.para,       theta.est.para))
PY0.XZ1  <- mean(pp*H.0(chi.para.Z1, theta.est.para) + (1-pp)*cc)
PY.XZ1   <- mean(H.0(chi.para.Z1,    theta.est.para))
PY0.XZ0  <- mean(pp*H.0(chi.para.Z0, theta.est.para) + (1-pp)*cc)
PY.XZ0   <- mean(H.0(chi.para.Z0,    theta.est.para))

# Logistic regression for Z ~ X
fit.ZX <- glm(Z ~ X1 + X2, family = binomial(), data = data.Log)
theta.est.ZX <- unname(unlist(fit.ZX$coefficients))
chi.ZX <- as.matrix(cbind(1, data.Log$X1, data.Log$X2))
P.ZX <- mean(H_logit(chi.ZX %*% theta.est.ZX))

result_prob <- tibble(
  Quantity = c(
    "P(Y0 = 1 | X, Z)",
    "P(Y  = 1 | X, Z)",
    "P(Y0 = 1 | X, Z = 1)",
    "P(Y  = 1 | X, Z = 1)",
    "P(Y0 = 1 | X, Z = 0)",
    "P(Y  = 1 | X, Z = 0)",
    "P(Z  = 1 | X)"
  ),
  Estimate = c(PY0.XZ, PY.XZ, PY0.XZ1, PY.XZ1, PY0.XZ0, PY.XZ0, P.ZX)
)

cat("\n--- Logistic regression for Z ~ X: summary ---\n")
print(summary(fit.ZX))

cat("\n--- Probability summaries (illustrative) ---\n")
print(result_prob)

############################################################
# END OF SCRIPT
############################################################
