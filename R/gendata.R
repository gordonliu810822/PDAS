gendata <- function(n,p,K,sigma,ratio,seednum,kind,rho,corr,snr){
#========================================================================#
# INPUTS:                                                                #
#     n   ---- number of samples                                         #
#     p   ---- signal length                                             #
#     K   ---- number of nonzero elements in the signal                  #
#  ratio  ---- range of value in the signal (= 10^ratio)                 #
#  sigma  ---- noise variance                                            #
# seednum ---- seed number                                               #
#   kind  ---- type of sample matirx                                     #
#          1 for Random Gaussian with auto coorelation (small p)         #
#          2 for correlated Random Gaussian  (large p)                   #
#    rho  ---- corr. coeff. for small p                                  #
#    corr  ---- corr. coeff. for large p (default 0.2)                   #
#    snr  ---- signal to noise ratio                                     #
# OUTPUTS:                                                               #
#     X   ---- normalized sample matrix                                  #
#     D   ---- vetor  makes X normalized                                 #
#     y   ---- data vector with noise                                    #
#    ye   ---- data vector without noise                                 #
#    xe   ---- true signal                                               #
#     A   ---- the support of xe                                         #
# (c) by Yuling Jiao(yulingjiaomath@whu.edu.cn) and Bangti Jin           #
# Created on Mar, 17, 2013                                               #
#========================================================================#

   if (!is.null(seednum)){
      set.seed(seednum)
   }
   xe = numeric(p);
   q =  sample(p);
   A = q[1:K];
   if (ratio !=0){
      vx = ratio*rnorm(K);
      vx = vx - min(vx);
      vx = vx/max(vx)*ratio;
      xe[A] = 10^vx*sign(rnorm(K));
   }
   else if (ratio == 0){
      xe[A] = sign(rnorm(K));
   }

   if (kind == 1){
      X = matrix(rnorm(n*p),nrow=n,ncol=p);
      if (rho!=0){
          Sigma = diag(p);
          for (k in 1:p){
             for ( l in 1:p){
                 Sigma[k,l] = rho^(abs(k-l));
             }
          }
          X = X%*%chol(Sigma);
      }
   }
   else if (kind == 2){
      X = matrix(rnorm(n*p),nrow=n,ncol=p);
      for (j in 2:(p-1)){
            X[,j] = X[,j] + corr*X[,j+1]+ corr*X[,j-1];
      }
   }
   mean.x = apply(X,2,mean);
   Xm = sweep(X,2,mean.x);
   X_l2 = apply(Xm^2,2,sum)
   D = X_l2^0.5
   Xp = sweep(Xm,2,D,"/")

   ye  = X%*%xe;
   noise = rnorm(n);

   if (is.null(snr)){
     y = ye + sigma*noise
   }
   else{
     sigma = sd(ye)/sqrt(snr)#norm(ye)/(sqrt(n*snr)); #% to gerate signal and noise with predifined snr
     y   = ye + sigma*noise;
   }

   obj = list(Sigma= Sigma,X = X, Xp = Xp, D = D, y = y, ye = ye, xe = xe, A = A);
   obj
}