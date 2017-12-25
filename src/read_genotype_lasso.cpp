#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <stdio.h>
#include "plinkfun.hpp"
//#include "readPlink.hpp"

using namespace std;
using namespace arma;
using namespace Rcpp;

#define MAX_LEN 20

void ReadFamFile(std::string stringname, CharacterVector FID, CharacterVector IID, NumericVector Phenotype, int N)
{
	FILE *stream;
	int a, b, c;
	char s1[MAX_LEN + 1], s2[MAX_LEN + 1];
	double pheno;

	stream = fopen(stringname.c_str(), "r");
	clock_t t1 = clock();
	//int i = 0;
	/* Put in various data. */
	for (int i = 0; i < N; i++){
		if (i % 100000 == 0 && i != 0){
			cout << i << "-th SNP" << ",";
			cout << "Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}

		fscanf(stream, "%s %s %i %i %i %lf", &s1[0], &s2[0], &a, &b, &c, &pheno);

		FID(i) = s1;
		IID(i) = s2;
		Phenotype(i) = pheno;
	}
}

List ReadGenotype(std::string stringname, arma::mat* X_std, int N, int P) {
	std::string famfile = stringname;
	famfile += ".fam";
	std::string bimfile = stringname;
	bimfile += ".bim";

	clock_t t1 = clock();
	//int N = getLineNum(famfile);
	
	CharacterVector FID(N), IID(N);
	NumericVector Phenotype(N);

	ReadFamFile(famfile, FID, IID, Phenotype, N);

	cout << endl;
	cout << "Start loading genotype: " << endl;

	//int P = getLineNum(bimfile);
	//cout << bimfile << endl;
	//cout << "N: " << N << endl;
	//cout << "P: " << P << endl;
	arma::mat bimmat;
	bimmat.load(bimfile, arma::csv_ascii);
	arma::vec c1 = bimmat.col(0);
	cout << bimmat.size() << endl;
	arma::vec index(22);
	double sum2 = 0;

	for (int i = 1; i <= 22; i++){
		index[i - 1] = sum(c1 == i);
		sum2 += index[i - 1];
		cout << "number of snps in chromsome " << i << " =" << index[i - 1] << endl;
	}
	unsigned* X = new unsigned[N * P];

	readPlink(stringname, N, P, X);
	cout << "Finished loading genotype " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
	arma::Mat<unsigned>* Xdata = new arma::Mat<unsigned>(X, N, P, false, false);
	cout << "Dimension of genotype: " << size(*Xdata) << endl;

	uword size1 = Xdata->n_rows;
	uword size2 = Xdata->n_cols;

	Xdata->replace(3, 0); //replace the missing value indicated by 3
	vec meanX(size2);
	vec sqrtsum(size2);
	//arma::fmat X_std(size1, size2);

	cout << "Start standardizing genotype: " << endl;
	t1 = clock();
	for (int i = 0; i < (int)(size2); i++) { //calculate the mean of the vector and sqrt sum
		meanX[i] = sum(Xdata->col(i))*1.0 / size1;
		vec v_i = conv_to<vec>::from(Xdata->col(i));
		v_i = v_i - meanX[i];
		(X_std->col(i)) = v_i;
		mat pd = v_i.t() * v_i;
		sqrtsum[i] = sqrt(pd.at(0));
		//X_std.col(i) = conv_to<fvec>::from(v_i) / sqrtsum[i];
		(X_std->col(i)) = v_i / sqrtsum[i];
	}
	cout << "Finished standardizing genotype in " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec." << endl;
	//arma::Mat<unsigned> XX = *Xdata;
	List output = List::create(//Rcpp::Named("X_std") = X_std,
		Rcpp::Named("sqrtsum") = sqrtsum,
		Rcpp::Named("Phenotype") = Phenotype);
	return output;
}
