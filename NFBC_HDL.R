#######################process NFBC HDL chr16(for test)#################
rm(list=ls())
setwd("G:\\Data\\GWAS\\NFBC")

trait = c("CRP","Glucose","Insulin","TC","HDL","LDL","TG", "BMI","SysBP","DiaBP");

cmd = "plink --noweb --bfile NFBC_filter_mph10 --chr 16 --make-bed --out NFBC_filter_mph10_chr16";
system(cmd)
fam =read.table("NFBC_filter_mph10_chr16.fam",header=F)
pheno = read.table("NFBC_filter_mph10.pheno",header=F);

fam[,6] = pheno[,5+2];
write.table(fam,"NFBC_filter_mph10_chr16.fam",sep=" ",quote=F, row.names=F,col.names=F);

#################################################################
setwd("/home/projects/11000369/Work/NFBC")
fam =read.table("NFBC_filter_mph10.fam",header=F)
pheno = read.table("NFBC_filter_mph10.pheno",header=F);
fam[,6] = pheno[,5+2];
write.table(fam[,1:6],"NFBC_filter_mph10.fam",sep=" ",quote=F, row.names=F,col.names=F);

####################################################################
rm(list=ls())
library(PDAS)
ls("package:PDAS")

#List pdas_geno(std::string stringname, std::string method, std::string sel, double al, double tau,
#double mu,double del, int MaxIt, double Lmax, double Lmin, int Nstep)
#setwd("/home/projects/11000369/Work/NFBC")
setwd("G:\\Data\\GWAS\\NFBC")
stringname = "NFBC_filter_mph10_chr16"#"NFBC_filter_mph10_chr16";
method = "l0"; #"l0","bridge","scad","mcp","cl1","lasso"
sel="bic";
al = 0; tau = 0; mu = 200; del = 0; MaxIt = 1; Lmax = 1; Lmin = 0.01; Nstep = 100;
out = pdas_geno(stringname, method, sel, al, tau, mu, del ,MaxIt, Lmax, Lmin, Nstep);
out$ithist_x = matrix(rep(1/out$x_var,dim(out$ithist_x)[2]),ncol=dim(out$ithist_x)[2]) * out$ithist_x
plot.PDAS(out,xvar="Index");

#setwd("/home/projects/11000369/Work/NFBC")
write.table(out$ithist_x,"NFBC_HDL_x.txt",sep=" ",quote=F, row.names=F,col.names=F);
write.table(out$ithist_bic,"NFBC_HDL_bic.txt",sep=" ",quote=F, row.names=F,col.names=F);
write.table(out$x_var,"NFBC_HDL_x_var.txt",sep=" ",quote=F, row.names=F,col.names=F);
write.table(out$Lam,"NFBC_HDL_Lam.txt",sep=" ",quote=F, row.names=F,col.names=F);

#List pdas_path(arma::mat X, arma::vec y, std::string method, std::string sel, double al, double tau,
#double mu, double del, int MaxIt, double Lmax, double Lmin, int N)
#out=BIC.PDAS(X, y, 'l0', 'vote', 0, tau= 2.7, 40, 0, 1, 1, 0.01, 100);
