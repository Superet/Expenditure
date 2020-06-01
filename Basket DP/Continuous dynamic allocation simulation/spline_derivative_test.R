# Test partial derivative of spline 
# Compute the partial derivative f_1(x1,x2)

# 1. Compute f_1(x1, x2nodes) along x2nodes, and then interpolate f_1(x1, x2) 
# 2. Compute f(x1nodes, x2), and then compute f_1(x1, x2) 

# Conclusion: both approaches yield the same

library(Rcpp)
library(RcppGSL)

src <- '
#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]

struct spl{
	gsl_spline *spline; gsl_interp_accel *acc; 
	NumericVector x; NumericVector y; double slope_first; double slope_end; 
			
	double splEval(double x0){
		gsl_set_error_handler_off();
		int nx = x.length();
		double out; 
		if((x0>=x[0] )& (x0<=x[nx-1])){
			out = gsl_spline_eval(spline, x0, acc);
		}else if(x0< x[0]){
			out = y[0] + slope_first * (x0 - x[0]);
		}else{
			out = y[nx-1] + slope_end * (x0 - x[0]);
		}
		return(out);
	}
	
	double spldevEval(double x0){
		int nx = x.length(); 
		double out; 
		if((x0>=x[0] )& (x0<=x[nx-1])){
			out = gsl_spline_eval_deriv(spline, x0, acc);
		}else if(x0< x[0]){
			out = slope_first;
		}else{
			out = slope_end;
		}
		return(out);
	}
	
	void splfree(){
		gsl_interp_accel_free(acc);
		gsl_spline_free (spline);
	}
};

struct spl spl_init(NumericVector x, NumericVector y){
	gsl_set_error_handler_off();
	struct spl out;
	int nx = x.size();
	out.x = x;
	out.y = y;
	out.acc= gsl_interp_accel_alloc();    	
	out.spline = gsl_spline_alloc( gsl_interp_cspline , nx );
	double eps = 0.0001;
	gsl_spline_init(out.spline, x.begin(), y.begin(), nx);	
	out.slope_first= (gsl_spline_eval (out.spline, x[0]+eps, out.acc) - y[0])/eps;
	out.slope_end	 = (y[nx-1] - gsl_spline_eval (out.spline, x[nx-1]-eps, out.acc))/eps;
	return(out);
}

// [[Rcpp::export]]
double myf1(NumericMatrix x, NumericVector y, NumericVector x0){
// This might be wrong;
	NumericVector x1=unique(x(_,0))	, x2 = unique(x(_,1));
	int N1 = x1.length(), N2 = x2.length(); 
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
		// Initiate the c-spline; 
	int nspl = N2;
	struct spl myspl[nspl];
	for(int i=0; i<nspl; i++){
		IntegerVector index	= i*N1 + seq_len(N1) - 1;
		NumericVector tmpy 	= y[index];
		myspl[i] = spl_init(x1, tmpy);
	}
	
	NumericVector y2(N2); 
	for(int i=0; i<N2; i++){
		y2(i) = myspl[i].spldevEval(x0(0)); 
	}
	
	struct spl new_spl; 
	new_spl = spl_init(x2, y2); 
	double out = new_spl.splEval(x0(1));
	new_spl.splfree(); 
	for(int i=0; i<nspl; i++){
		myspl[i].splfree(); 
	}
	return(out);
}

// [[Rcpp::export]]
double myf2(NumericMatrix x, NumericVector y, NumericVector x0){
// This might be wrong;
	NumericVector x1=unique(x(_,0))	, x2 = unique(x(_,1));
	int N1 = x1.length(), N2 = x2.length(); 
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
		// Initiate the c-spline; 
	int nspl = N1;
	struct spl myspl[nspl];
	for(int i=0; i<nspl; i++){
		IntegerVector index	= (seq_len(N2)-1)*N1 + i;
		NumericVector tmpy 	= y[index];
		myspl[i] = spl_init(x2, tmpy);
	}
	
	NumericVector y1(N1); 
	for(int i=0; i<N1; i++){
		y1(i) = myspl[i].splEval(x0(1)); 
	}
	
	struct spl new_spl; 
	new_spl = spl_init(x1, y1); 
	double out = new_spl.spldevEval(x0(0));
	new_spl.splfree(); 
	for(int i=0; i<nspl; i++){
		myspl[i].splfree(); 
	}
	return(out);
}
'
sourceCpp(code = src)

x1 <- c(1:10)
x2 <- seq(2, 8, 1)
x <- cbind(rep(x1, length(x2)), rep(x2, each=length(x1)))
y  <- runif(nrow(x))
x0 <- c(3.5, 5.9)

myf1(x, y, x0)
myf2(x, y, x0)




