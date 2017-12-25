//
//  plinkfun.hpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef read_genotype_lasso_hpp
#define read_genotype_lasso_hpp

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <string.h>
using namespace std;
using namespace Rcpp;
using namespace arma;

void ReadFamFile(std::string stringname, CharacterVector FID, CharacterVector IID, NumericVector Phenotype, int N);
List ReadGenotype(std::string stringname, arma::mat* X_std, int N, int P);

#endif /* read_genotype_lasso_hpp */
