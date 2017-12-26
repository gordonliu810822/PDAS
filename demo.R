############################################
rm(list=ls())
library(PDAS)
ls("package:PDAS")

##---- small scale Regression problem ----- #
#p: length of signal x
#n: number of samples
#K: number of nonzero elements in x
#ratio: range of value in x (10^ratio)
#sigma: noise standard deviation
#kind: standard Regression setting
#seednum: the seed number
#rho: correlation
#snr: signal-noise-ratio

n = 200; p = 1000; K = 10; ratio = 1; sigma=0.5; kind = 1; seednum = NULL; rho = 0.5;
corr = NULL; snr = 2;
dat = gendata(n,p,K,sigma,ratio,seednum,kind,rho,corr,snr)
X = dat$Xp; y = dat$y;  D = dat$D; ye = dat$ye; A = dat$A; xe = dat$xe;

#alpha  ---  regularization parameter (default: 0)                #
#method ---  nonconvex model                                      #
#tau    ---  concavity parameter                                  #
#mu     ---  stopping parameter  (default: n/log(n)               #
#MaxIt  ---  maximum number of iterations (default: 1)            #
#del    --- noise level (default: 0)                              #
#sel    --- tunning parameter selection method (defaut:vote)      #
#Lmin   --- minimun in Lam (defaut: 1e-3)                         #
#N      --- length of path (defaut: 100)                          #
method = "lasso";
sel="bic"; alpha = 0; mu = n/log(n); del = 0; MaxIt = 1; Lmax = 1; tau = 0; Lmin = 0.01; Nstep = 100;

out=BIC.PDAS(X, y, method, sel, alpha, tau, mu, del,MaxIt, Lmax, Lmin, Nstep);
out2=list(ithist_x=matrix(rep(1/D,dim(out$ithist_x)[2]),ncol=dim(out$ithist_x)[2])*
out$ithist_x,ithist_as = out$ithist_as, ithist_bic = out$ithist_bic);
plot.PDAS(out2,xvar="Index");
bic = out2$ithist_bic;
bh = out2$ithist_x[,which(bic==min(bic))];
