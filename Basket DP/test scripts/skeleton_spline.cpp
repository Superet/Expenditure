#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
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
		if(x0>=x[0] & x0<=x[nx-1]){
			out = gsl_spline_eval(spline, x0, acc);
		}else if(x0< x[0]){
			out = y[0] + slope_first * (x0 - x[0]);
		}else{
			out = y[nx-1] + slope_end * (x0 - x[0]);
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
double splineC(NumericVector x, NumericVector y, double xnew){
	gsl_set_error_handler_off();
	struct spl myspl = spl_init(x, y);
	double out = myspl.splEval(xnew);
	return(out);
}