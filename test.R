############################################
rm(list=ls())
library(PDAS)
ls("package:PDAS")

#setwd("/Users/nus/GoogleDrive/Work/JianHuang/updas")
#setwd("/Users/nus/GoogleDrive/Work/JianHuang/updas")
setwd("G:\\Data\\GoogleDrive\\Work\\JianHuang\\updas")
Ae = as.numeric(unlist(read.table("Ae.mat",header=F)))
D = as.numeric(unlist(read.table("D.mat",header=F)))
X = matrix(as.numeric(unlist(read.table("X.mat",header=F))),nrow=200);
xe = as.numeric(unlist(read.table("xe.mat",header=F)))
y = as.numeric(unlist(read.table("y.mat",header=F)))
ye = as.numeric(unlist(read.table("ye.mat",header=F)))
Xty = t(X)%*%y;
p = dim(X)[2]
#pdas(arma::mat X, arma::vec y, double al, double lam, std::string method, double tau, double mu, int MaxIt, arma::vec Xty, arma::vec x0, arma::vec d0, arma::uvec Aoo, arma::mat Go)
#List pdas_path(arma::mat X, arma::vec y, std::string method, std::string sel, double al,
#double tau, double mu, double del, int MaxIt, double Lmax, double Lmin, int N)
# outp = pdas_path(X, y, 'l0', 'vote', 0, 0, 37.7478, 0, 1, 1, 1e-8, 100);
#outp = pdas_path(X, y, 'lasso', 'vote', 0, 1, 37.7478, 0, 1, 1, 1e-8, 100);
#out=list(ithist_x=D*outp$ithist_x,ithist_as = outp$ithist_as);
#out$ithist_x = D*outp$ithist_x
n = 200; p = 1000; K = 10; ratio = 1; sigma=0.5; kind = 1; seednum = NULL; rho = 0.5;
corr = NULL; snr = 2;
dat = gendata(n,p,K,sigma,ratio,seednum,kind,rho,corr,snr)
X = dat$Xp; y = dat$y;  D = dat$D; ye = dat$ye; A = dat$A; xe = dat$xe;#Sigma=dat$Sigma; tmp = chol(Sigma)
write.table(X,"X.dat",sep="\t",quote=F, row.names=F,col.names=F)
write.table(y,"y.dat",sep="\t",quote=F, row.names=F,col.names=F)
write.table(D,"D.dat",sep="\t",quote=F, row.names=F,col.names=F)

#library(matlabr)
#options(matlab.path = "G:\\Data\\GoogleDrive\\Work\\JianHuang\\updas\\use all")
#have_matlab()
#code=paste("n = 200; p = 1000; K = 10; ratio = 1; sigma=0.5; kind = 1; seednum = NULL; rho = 0.5;
#corr = NULL; snr = 50;[X,D,y,ye,xe,Ae] = gendata(n,p,K,sigma,ratio,seednum,kind,rho);")
#res = run_matlab_code(code)

#List pdas_path(arma::mat X, arma::vec y, std::string method, std::string sel, double al, double tau,
#double mu, double del, int MaxIt, double Lmax, double Lmin, int N)
out=BIC.PDAS(X, y, 'l0', 'vote', 0, tau= 2.7, 40, 0, 1, 1, 0.01, 100);
out2=list(ithist_x=matrix(rep(1/D,dim(out$ithist_x)[2]),ncol=dim(out$ithist_x)[2])*
out$ithist_x,ithist_as = out$ithist_as, ithist_bic = out$ithist_bic);
plot.PDAS(out2,xvar="Index");
bic = out2$ithist_bic;
bh = out2$ithist_x[,which(bic==min(bic))];


#bic =outp$ithist_bic
#ind = which(bic == min(bic))
#x = outp$ithist_x[,ind]
#x = D*x;
#
#err_l2 = sqrt(sum((x-xe)^2))/sqrt(sum(xe^2))
#err_linf = max(abs(x-xe));
#index = 1:p;
#A = index[x!=0];
#lengthA = length(A);
#
#plot(1:p,xe,pch = 2)
#lines(1:p,x,type = 'p',pch = 4)


library(glmnet)
fit = glmnet(X, y,lambda.min.ratio=0.0001,standardize=T)
fit$beta=D*fit$beta
plot(fit)

#tLL <- fit$nulldev - deviance(fit)
#k <- fit$df
#BIC<-log(n)*k - tLL
#cv = cv.glmnet(X, y,lambda.min.ratio=0.2)
#plot(cv)
#########################################






