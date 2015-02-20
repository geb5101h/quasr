// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

double thresh(
	double x, 
	double lam
){
	double out;
	out = (std::abs(x)-lam);
	if(out>0){
		out = out * x/std::abs(x);
		return(out);
	}else{
		return(0.0);
	}
};
// [[Rcpp::export]]
arma::mat sparseprec(
	arma::mat sigma,
	double lam
){
	int d= sigma.n_rows;
	arma::mat omega(d,d);
	arma::mat omega_old(d,d);
	omega.eye();
	//omega = sigma;
	//omega += lam*arma::eye(d,d);
	double s=0;
	double tol=1e5;
	int iter = 1;
	while(tol>1e-4 & iter<10000){
	omega_old=omega;
	for(int i=0; i<d ; ++i){
		for(int j=i; j<d; ++j){
			s=0;
			if(i>0){
			s += .5 * dot(sigma(i,arma::span(0,i-1)), omega(j,arma::span(0,i-1)));
			}
			if(j>0){
				s+= .5*dot(sigma(j,arma::span(0,j-1)), omega(i,arma::span(0,j-1)));
			}
				if(i+1<d){
					s+=.5*dot(sigma(i,arma::span(i+1,d-1)) , omega(j,arma::span(i+1,d-1)));
				}
				if(j+1<d){
	 				s+=.5*dot(sigma(j,arma::span(j+1,d-1)) , omega(i,arma::span(j+1,d-1)));
 				}
			if(i==j) s-=1;
			s/=(-.5*sigma(i,i)-.5*sigma(j,j));
			omega(i,j)=thresh(s,lam);
			omega(j,i)=thresh(s,lam);
		}
	}
	//tol = sqrt(accu(arma::pow(omega-omega_old,2)));
	tol= arma::accu(arma::abs(omega-omega_old))/(d*d);
	iter+=1;
	}
	if(iter == 10000){
		Rcpp::Rcout << "Max iter reached: " << iter << std::endl;	
	}
	return(omega);
};