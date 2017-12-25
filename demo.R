############################################
rm(list=ls())
library(PDAS)
ls("package:PDAS")

n = 200; p = 1000; K = 10; ratio = 1; sigma=0.5; kind = 1; seednum = NULL; rho = 0.5;
corr = NULL; snr = 2;
dat = gendata(n,p,K,sigma,ratio,seednum,kind,rho,corr,snr)
X = dat$Xp; y = dat$y;  D = dat$D; ye = dat$ye; A = dat$A; xe = dat$xe;

out=BIC.PDAS(X, y, 'lasso', 'vote', 0, tau= 2.7, 200, 0, 1, 1, 0.01, 100);
out2=list(ithist_x=matrix(rep(1/D,dim(out$ithist_x)[2]),ncol=dim(out$ithist_x)[2])*
out$ithist_x,ithist_as = out$ithist_as, ithist_bic = out$ithist_bic);
plot.PDAS(out2,xvar="Index");
bic = out2$ithist_bic;
bh = out2$ithist_x[,which(bic==min(bic))];


#library(glmnet)
#fit = glmnet(X, y,lambda.min.ratio=0.01)
#fit$beta=(1/D)*fit$beta
#plot(fit)
#cv = cv.glmnet(X,y,lambda.min.ratio=0.01)
#plot(cv)
