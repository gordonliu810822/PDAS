rm(list=ls())
library(PDAS)
ls("package:PDAS")

#List pdas_geno(std::string stringname, std::string method, std::string sel, double al, double tau,
#double mu,double del, int MaxIt, double Lmax, double Lmin, int Nstep)
setwd("/home/projects/11000369/Work/NFBC")
#setwd("G:\\Data\\GWAS\\NFBC")
stringname = "NFBC_filter_mph10"#"NFBC_filter_mph10_chr16";
#method = "lasso"; tau = 0; Lmin = 0.3;
method = "l0"; tau = 0; Lmin = 1e-3;
#method = "bridge"; tau = 0.5; Lmin = 0.05;
sel="bic";
al = 0; mu = 600; del = 0; MaxIt = 1; Lmax = 1; Nstep = 100;
out = pdas_geno(stringname, method, sel, al, tau, mu, del ,MaxIt, Lmax, Lmin, Nstep);
out$ithist_x = matrix(rep(1/out$x_var,dim(out$ithist_x)[2]),ncol=dim(out$ithist_x)[2]) * out$ithist_x
#plot.PDAS(out,xvar="Index");

#setwd("/home/projects/11000369/Work/NFBC")
write.table(out$ithist_x,paste("NFBC_HDL_x_",method,".txt",sep=""),sep=" ",quote=F, row.names=F,col.names=F);
write.table(out$ithist_bic,paste("NFBC_HDL_bic_",method,".txt",sep=""),sep=" ",quote=F, row.names=F,col.names=F);
write.table(out$x_var,paste("NFBC_HDL_x_var_",method,".txt",sep=""),sep=" ",quote=F, row.names=F,col.names=F);
write.table(out$Lam,paste("NFBC_HDL_Lam_",method,".txt",sep=""),sep=" ",quote=F, row.names=F,col.names=F);
