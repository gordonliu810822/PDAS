plot.PDAS=function(out, xvar="norm",label=FALSE,...){
  #xvar=match.arg(xvar)
  #plotCoef(x$beta,lambda=x$lambda,df=x$df,dev=x$dev.ratio,label=label,xvar=xvar,...)
  #par(mfrow = c(1,1))
  plotCoef(out$ithist_x,lambda=out$Lam,df=out$ithist_as,bic = out$ithist_bic, label=label,xvar=xvar,...)
}
