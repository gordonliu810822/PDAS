//
//  plinkfun.hpp
//  PlinkRead
//
//  Created by DaiMingwei on 16/10/6.
//  Copyright © 2016年 daviddai. All rights reserved.
//

#ifndef pdas_hpp
#define pdas_hpp

#include <RcppArmadillo.h>
#include <Rcpp.h>
//#include <omp.h>
//#include "rcpp_hello.hpp"
using namespace Rcpp;
using namespace arma;


vec getdualA(vec x, vec d, vec pd, uvec A, double lam, double tau, std::string method);
Rcpp::List pdas(arma::mat X, arma::vec y, double al, double lam, std::string method, double tau, double mu, int MaxIt, arma::vec Xty, arma::vec x0, arma::vec d0, arma::uvec Aoo, arma::mat Goo);
Rcpp::List pdas_path(arma::mat X, arma::vec y, std::string method, std::string sel, double al, double tau, double mu, double del, int MaxIt, double Lmax, double Lmin, int N);

#endif /* plinkfun_hpp */
