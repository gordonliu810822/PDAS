#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <armadillo>
#include <iostream>
#include "plinkfun.hpp"
#include "read_genotype_lasso.hpp"
#include "pdas.hpp"

using namespace std;
using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List pdas_geno(std::string stringname, std::string method, std::string sel, double al, double tau, double mu, double del, int MaxIt, double Lmax, double Lmin, int Nstep){

	std::string famfile = stringname;
	famfile += ".fam";
	std::string bimfile = stringname;
	bimfile += ".bim";
	clock_t t0 = clock();

	int P = getLineNum(bimfile);
	int N = getLineNum(famfile);

	//mat X_std(N,P);
	arma::mat tmp(N, P);
	arma::mat* X_std = new arma::mat(tmp.memptr(), N, P, false, false);
	Rcpp::List input = ReadGenotype(stringname, X_std, N, P);
	arma::vec pheno = as<vec>(input["Phenotype"]);
	arma::vec x_var = input["sqrtsum"];
	//centering outcome
	vec Y = pheno;
	double meanY = sum(pheno) / N;
	Y = Y - meanY;

	clock_t t1 = clock();
	cout << "Start fitting " << method << ":" << endl;
	List out = pdas_path(tmp, Y, method, sel, al, tau, mu, del, MaxIt, Lmax, Lmin, Nstep);
	cout << "Finished fitting " << method << " in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
	cout << endl;
	cout << "Total analysis done in " << (clock() - t0)*1.0 / CLOCKS_PER_SEC << " sec." << endl;

	out["x_var"] = x_var;
	out["x"] = tmp;
	return out;
}