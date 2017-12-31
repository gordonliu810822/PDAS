#include <RcppArmadillo.h>
#include <Rcpp.h>
//#include <omp.h>
//#include "rcpp_hello.hpp"
using namespace Rcpp;
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

// Jin Liu, Duke-NUS Medical School
// email: jin.liu@duke-nus.edu.sg


vec getdualA(vec x, vec d, vec pd, uvec A, double lam, double tau, char methodn){
    //-----------------------------------------------------------------------//
    // compute the dual variable for nonconvex penalties on the active set   //
    // INPUTS:                                                               //
    //     x  --- primal variable                                            //
    //     d  --- dual primal variable                                       //
    //     pd  --- current (x + d)                                           //
    //     lam  ---  regualrization paramter for sparsity                    //
    //     method ---  nonconvex model                                       //
    //     tau  ---  concavity parameter                                     //
    //     A  --- active set                                                 //
    // OUTPUTS:                                                              //
    //     dA  --- dual variable on A                                        //
    // (c) by Jin Liu (jin.liu@duke-nus.edu.sg)                              //
    // Created on Dec, 17, 2017                                              //
    //-----------------------------------------------------------------------//
    int s = A.n_elem;
    vec dualA;   
	dualA.zeros(s);
    
    vec dA = d(A), pdA = pd(A), xA = x(A);
    double Tstar;
    IntegerVector ii, i0, i1;
    
	switch (methodn){
	case 'A':
		pdA = pd(A);
		dualA = lam*sign(pdA);
		//ii = find(abs(xA) >= lam);
		//dualA(as<uvec>(ii)) == lam*sign(pdA(as<uvec>(ii)));
		break;
	case 'C':
		Tstar = pow(2 * lam*(1 - tau), 1 / (2 - tau));
		ii = find(abs(xA) >= Tstar);
		dualA(as<uvec>(ii)) = lam*tau*pow((abs(xA(as<uvec>(ii)))), tau) / (xA(as<uvec>(ii)));
		break;
	case 'D':
		pdA = pd(A);
		ii = intersect(as<IntegerVector>(wrap(find(abs(pdA)>lam))),
			as<IntegerVector>(wrap(find(abs(pdA)<(tau + 0.5)*lam))));
		dualA(as<uvec>(ii)) = lam*sign(pdA(as<uvec>(ii)));
		break;
	case 'E':
		ii = intersect(as<IntegerVector>(wrap(find(abs(pdA) <= 2 * lam))),
			as<IntegerVector>(wrap(find(abs(pdA)>lam))));
		dualA(as<uvec>(ii)) = lam*sign(pdA(as<uvec>(ii)));
		i0 = find((xA % dA) >= 0);
		i1 = intersect(as<IntegerVector>(wrap(find(abs(pdA)>2 * lam))),
			as<IntegerVector>(wrap(find(abs(pdA)< tau*lam))));
		ii = intersect(i1, i0);
		dualA(as<uvec>(ii)) = (tau*lam*sign(pdA(as<uvec>(ii))) - xA(as<uvec>(ii))) / (tau - 1);
		break;
	case 'F':
		i0 = as<IntegerVector>(wrap(find((xA % dA) >= 0)));
		i1 = intersect(as<IntegerVector>(wrap(find(abs(pdA)>lam))),
			as<IntegerVector>(wrap(find(abs(pdA)< tau*lam))));
		ii = intersect(i1, i0);
		dualA(as<uvec>(ii)) = (lam*tau*sign(pdA(as<uvec>(ii))) - xA(as<uvec>(ii))) / tau;
		break;
	default:
		dualA.zeros(s); //'l0'
	}

    return dualA;
}

// [[Rcpp::export]]
Rcpp::List pdas(arma::mat X, arma::vec y, double al, double lam, char methodn, double tau, double mu, double weight, int MaxIt, arma::vec Xty, arma::vec x0, arma::vec d0, arma::uvec Aoo, arma::mat Goo) {
    //-------------------------------------------------------------------------//
    // Solving nonconvex problem                                               //
    //      1/2||X*x-y||^2 + alpha/2||x||^2 + rho_{lam,tau}(x)                 //
    // by a primal-dual active set algorithm                                   //
    //                                                                         //
    // INPUTS:                                                                 //
    //     X   ---  matrix, normalized                                         //
    //     y   ---  data vector                                                //
    //        al ---  regularization parameter (default: 0)                    //
    //     lam  ---  regualrization paramter for sparsity                      //
    //     method ---  nonconvex model                                         //
    //     tau  ---  concavity parameter                                       //
    //     mu  ---  stopping parameter                                         //
    //     MaxIt  ---  maximum number of iterations (default: 5)               //
    //     x0    ---  initial guess for x (default: 0)                         //
    //     d0    ---  initial guess for d (default: Xty)                       //
    // OUTPUTS:                                                                //
    //      x  ---  solution                                                   //
    //      d  ---  dual                                                       //
    //      G  ---  X(:,A)'*X(:,A) + alpha I                                   //
    //      A  ---  active set                                                 //
    //      s  ---  size of active set                                         //
    //      L  ---  loss - norm(y)^2                                           //
    //      It ---  # of iter for the PDAS stop                                //
    // (c) by Jin Liu (jin.liu@duke-nus.edu.sg)                                //
    // Created on Dec, 17, 2017                                                //
    //-------------------------------------------------------------------------//

	uword p = X.n_cols;
	//uword n = X.n_rows;
    uvec Ao = Aoo;
    
    vec x = x0, d = d0;
    mat Go = Goo, G;
    int so = Ao.n_elem;
     
    double T;
	switch (methodn){
	case 'B':
		T = sqrt(2 * lam);
		break;
	case 'C':
		T = (2 - tau)*pow((2 * (1 - tau)), ((tau - 1) / (2 - tau)))*pow(lam, (1 / (2 - tau)));
		break;
	default:
		T = lam;
	}
  
    // initializing
	//cout << "break1 ..." << weight << endl;
	vec pd = weight*x + (1-weight)*d;
	//cout << "break2 ..." << pd(0) << ";" << pd(1) << endl;
	//cout << "break1 ... " << pd.n_elem << endl;
    uvec A = find(abs(pd) > T);
    int s = A.size();
    uvec pdt = A;//abs(pd) > T;
    uvec tpd;
    double L = 0;
    int It = 0;
    vec rhs, xA, dA;
    mat Xa;
    
    
	while (It  < MaxIt && s < mu){
        It = It + 1;
        //getdualA(vec x, vec d, vec pd, IntegerVector A, double lam, double tau, std::string method)
		dA = getdualA(x, d, pd/weight, A, lam, tau, methodn);
 		vec XtyA = Xty(A);
        rhs = XtyA - dA;
        
		if (s == so && s == Go.n_rows){
			if (s == 0){
				G = Go;
			}
			else {
				IntegerVector Atmp = intersect(as<IntegerVector>(wrap(Ao)),as<IntegerVector>(wrap(A)));
				if (Atmp.size() == A.n_elem){
				//if (Ao.n_elem == A.n_elem && sum(Ao) == sum(A)){
					G = Go;
				}
				else {
					Xa = X.cols(A);
					mat tmp; tmp.eye(s, s);
					G = Xa.t() * Xa + al * tmp;
				}
			}      
        }
        else {
            Xa = X.cols(A);
            mat tmp; tmp.eye(s,s);
            G = Xa.t() * Xa + al * tmp;
        }
        
        x.zeros(p);
        xA = solve(G,rhs);
        x(A) = xA;
        
        L = dot(xA, XtyA+dA);
        d = X.t()*(y-X*x);
		pd = weight*x + (1 - weight)*d;
		//cout << "break3 ..." << pd(0) << ";" << pd(1) << endl;
        Ao = A;
        Go = G;
        so = s;
        A = find(abs(pd) > T);
        s = A.n_elem;
        tpd = pdt;
        pdt = find(abs(pd) > T);
                
        if (A.n_elem > mu){
            break;
        }
        //IntegerVector pdtmp = intersect(as<IntegerVector>(wrap(pdt)),
        //                               as<IntegerVector>(wrap(tpd)));
        if ( tpd.n_elem == pdt.n_elem && sum(tpd) == sum(pdt)){
            break;
        }
		
	}
    if ( It > 0 ){
        s = so;
        A = Ao;
    }
    
	List out = List::create(Rcpp::Named("x") = x,
		Rcpp::Named("d") = d,
		Rcpp::Named("G") = G,
		Rcpp::Named("A") = A,
        Rcpp::Named("s") = s,
        Rcpp::Named("L") = L,
        Rcpp::Named("It") = It);
	return out;
}


// [[Rcpp::export]]
Rcpp::List pdas_path(arma::mat X, arma::vec y, std::string method, std::string sel, double al, double tau, double mu, double del, double weight, int MaxIt, double Lmax, double Lmin, int N){
    //-----------------------------------------------------------------------//
    // minimize 1/2||X*x-y||^2 + alpha/2||x||^2 + rho_{lambda,tau}(x)        //
    // with PDAS algorithm, stopping criterion: discrepency principle        //
    // Lam ={lam_{1},...,lam_{N}}.                                           //
    //                                                                       //
    // INPUTS:                                                               //
    //     X  ---  sampling matrix (normalized)                              //
    //     y  ---  data vector                                               //
    //   opts ---  structure containing                                      //
    //        tau  --- concavity parameter                                   //
    //       alpha --- regularization parameter (defaut: 0)                  //
    //         N   --- length of path (defaut: 100)                          //
    //       Lmin  --- minimun in Lam (defaut: 1e-8)                         //
    //        del  --- noise level                                           //
    //        mu   --- stop if ||x_{lambda_{k}}||> opts.mu (defaut:n/log(n)  //
    //        x0   --- initial value for PDAS (defaut: 0)                    //
    //     MaxIter --- maximum number of iterations in PDAS  (defaut: 1)     //
    //         sel --- tunning parameter selection method (defaut:vote)      //
    // OUTPUTS:                                                              //
    //     x   ---- reconstructed signal                                     //
    //    lam  ---- regu. parameter (by discrepancy principle/bic/vote)      //
    //     id  ---- seletion index of solution                               //
    //  ithist ---- structure on iteration history, containing               //
    //          .x  --- solution path                                        //
    //          .as --- size of active set                                   //
    //          .it --- # of iteration on the path                           //
    //          .res--- residual on the path                                 //
    // (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)                        //
    //         Bangti Jin (bangti.jin@gmail.com)                             //
    // Created on March 17, 2013                                             //
    // ======================================================================//

    uword p = X.n_cols;
    uword n = X.n_rows;
	char methodn;
	if (method.compare("l0") == 0){
		cout << "l_0 model is running ... " << endl;
		methodn = 'B';
	}
	else if (method.compare("bridge") == 0){
		cout << "Bridge model is running ... " << endl;
		methodn = 'C';
	}
	else if (method.compare("scad") == 0){
		cout << "SCAD model is running ... " << endl;
		methodn = 'E';
	}
	else if (method.compare("mcp") == 0){
		cout << "MCP model is running ... " << endl;
		methodn = 'F';
	}
	else if (method.compare("cl1") == 0){
		cout << "Capped-l_1 model is running ... " << endl;
		methodn = 'D';
	}
	else if (method.compare("lasso") == 0){
		cout << "Lasso model is running ... " << endl;
		methodn = 'A';
	}
	else {
		cout << "Undefined method is PDAS! " << endl;
	}
    
    vec Lam1   = exp(linspace(log(Lmax),log(Lmin),N));
	vec Lam = Lam1;//Lam1.subvec(1,Lam1.n_elem-1);
    vec Xty = X.t()*y;
    double linf = norm(Xty,"inf");
    double cnst;

	switch (methodn){
	case 'B':
		cnst = pow(linf, 2) / 2;
		break;
	case 'C':
		cnst = pow(linf / (2 - tau)*pow(2 * (1 - tau), (tau - 1) / (2 - tau)), 2 - tau);
		break;
	default:
		cnst = linf;
	}
    /*if (method.compare("l0") == 0){
        cnst = pow(linf,2)/2;
    }
    else if (method.compare("bridge") == 0){
        cnst = pow(linf/(2-tau)*pow(2*(1-tau),(tau-1)/(2-tau)),2-tau);
    }
    else {
        cnst = linf;
    }*/
	Lam = Lam*cnst*(1-weight);
    
    mat ithist_x = zeros<mat>(p,N);
    uvec ithist_as = zeros<uvec>(N);
    uvec ithist_it = zeros<uvec>(N);
    vec ithist_res = zeros<vec>(N);
    vec ithist_bic = zeros<vec>(N);
    
    double ny = pow(norm(y),2);
    vec d0 = Xty, x0 = zeros<vec>(p);
    uvec A0 = find(x0 == 0), A;
    mat G0, G;
    int s, It, id;
    vec x, d;
    double L, res, lam;

	clock_t t1 = clock();
    for (unsigned int k = 0; k < Lam.n_elem; k ++ ){
        lam = Lam(k);
		//pdas(arma::mat X, arma::vec y, double al, double lam, char methodn, double tau, double mu, double weight, int MaxIt, arma::vec Xty, arma::vec x0, arma::vec d0, arma::uvec Aoo, arma::mat Goo)

		List out = pdas(X, y, al, lam, methodn, tau, mu, weight, MaxIt, Xty, x0, d0, A0, G0);
        vec tmp1 = out["x"]; x = tmp1;
        vec tmp2 = out["d"]; d = tmp2;
        mat tmp3 = out["G"]; G = tmp3;
        uvec tmp4 = out["A"]; A = tmp4;
        int tmp5 = out["s"]; s = tmp5;
        double tmp6 = out["L"]; L = tmp6;
        int tmp7 = out["It"]; It = tmp7;
        // warm start
        x0 = x;
        d0 = d;
        G0 = G;
        A0 = A;
        
        ithist_x.col(k) = x;
        ithist_it(k) = It;
        ithist_as(k) = s;
        res = ny - L;
        ithist_res(k) = res;

		ithist_bic(k) = log(res / n) + s*log(n) / n;//log(res / n) + s*log(n) / n;
		//log(res / n) + s*log(n) / n;
        
		if (k % 5 == 4){
			cout << k + 1 << "-th lambda: df " << s << "; bic " << ithist_bic(k) << "; lambda " << lam;
			cout << "; Elapsed time is " << (clock() - t1)*1.0 / CLOCKS_PER_SEC << " sec" << endl;
		}

        if ( res <= del ){
            cout << "Discrepancy principle is satisfied. " <<endl;
            lam = Lam(k);
            x = ithist_x.col(k);
            id = k;
            break;
        }
        if ( s > mu ){
            cout << "Exceed maximum degree of freedoms. " <<endl;
            lam = Lam(k);
            break;
        }
    }
    
	uvec tmplst = find(ithist_as != 0);
	uvec lst = zeros<uvec>(tmplst.n_elem);
	lst(0) = 0;
	lst.subvec(1, tmplst.n_elem - 1) = tmplst.subvec(0, tmplst.n_elem - 2);
	ithist_x = ithist_x.cols(lst);
	ithist_it = ithist_it(lst);
	ithist_bic = ithist_bic(lst);
	ithist_res = ithist_res(lst);
	Lam = Lam(lst);
	ithist_as = ithist_as(lst);

    /*if ( res > del){
        if (sel.compare("vote") == 0){
            
            ithist_x = ithist_x.cols(find(ithist_as !=0));
            ithist_bic = ithist_bic(find(ithist_as !=0));
            ithist_as = ithist_as(find(ithist_as !=0));
            //ii = find(ithist_as == )
        }
        else {
            ithist_x = ithist_x.cols(find(ithist_as !=0));
            ithist_bic = ithist_bic(find(ithist_as !=0));
            ithist_as = ithist_as(find(ithist_as !=0));
            
            uvec ii = find(ithist_bic == min(ithist_bic));
            id = ii(ii.n_elem-1);
            x  = ithist_x.col(id);
            lam= Lam(id);
        }
    }*/
    List out = List::create(/*Rcpp::Named("x") = x,
                            Rcpp::Named("lam") = lam,
                            Rcpp::Named("Lam") = Lam,*/
							Rcpp::Named("id") = id,
							Rcpp::Named("Lam") = Lam,
                            Rcpp::Named("ithist_x") = ithist_x,
                            Rcpp::Named("ithist_it") = ithist_it,
                            Rcpp::Named("ithist_as") = ithist_as,
                            Rcpp::Named("ithist_res") = ithist_res,
                            Rcpp::Named("ithist_bic") = ithist_bic);
    return out;
}

                                          
                                          

