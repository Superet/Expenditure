#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_cdf.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppGSL)]]

int findInterval1(double x, NumericVector breaks);

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

struct spl3d{
	 NumericMatrix x;  NumericVector y; struct spl *spl_vec;  int n; 
	 NumericVector x1;  NumericVector x2;  NumericVector x3;
	
	double splEval(NumericVector xnew){
		gsl_set_error_handler_off();
		int Nx1 = x1.length(), Nx2 = x2.length(), Nx3 = x3.length();
		int x3_index = findInterval1(xnew(2), x3);		//Find which discrete value it is in the 3rd dimension; 
		NumericVector y2(Nx2);
		
		// Interpolate the first dimension conditional the other values; 
		for(int j=0; j<Nx2; j++){
			y2(j) = spl_vec[(x3_index-1)*Nx2+j].splEval(xnew(0));
		}
		
		// Set up another spl object for the 2nd dimension; 
		struct spl spl2 = spl_init(x2, y2);
		double out = spl2.splEval(xnew(1));
		spl2.splfree();
		return(out);		
	}
	
	void splfree(){
		(*spl_vec).splfree();
// 		for(int i=0; i<n; i++){
// 			spl_vec[i].splfree();
// 		}
	}
};

struct spl3d spl3d_init(NumericMatrix x, NumericVector y){
	struct spl3d out; 
	out.x 	= x;
	out.y 	= y;
	NumericVector x1 = unique(x(_,0)), x2 = unique(x(_,1)), x3 = unique(x(_,2));
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
	std::sort(x3.begin(), x3.end()); 
	out.x1 	= x1;
	out.x2 	= x2;
	out.x3 	= x3;
	int Nx1 = x1.length(), Nx2 = x2.length(), Nx3 = x3.length();
	int n 	= Nx2*Nx3;
	out.n	= n;
	struct spl spl_vec[n];
	for(int i=0; i<n; i++){
		IntegerVector index	= i*Nx1 + seq_len(Nx1) - 1;
		NumericVector tmpy 	= y[index];
		spl_vec[i] = spl_init(x1, tmpy);
	}
	out.spl_vec = spl_vec;
	return(out);
}

// [[Rcpp::export]]
int findInterval1(double x, NumericVector breaks) {
	int out;
	out = std::distance(breaks.begin(), std::upper_bound(breaks.begin(), breaks.end(), x));
	return out;
}

// [[Rcpp::export]]
double my_fn(NumericMatrix x, NumericVector y, NumericVector xnew){
	gsl_set_error_handler_off();
	struct spl3d out = spl3d_init(x,y);
	double a = out.splEval(xnew);
// 	out.splfree();
	return(a);
}

// [[Rcpp::export]]
NumericVector my_fn1(NumericMatrix x, NumericVector y, NumericMatrix xnew){
	gsl_set_error_handler_off();
	int n = xnew.nrow();
	NumericVector out1(n);
	for(int i=0; i<n; i++){
		struct spl3d tmp = spl3d_init(x,y);;
		out1(i) = tmp.splEval(xnew(i,_));
		tmp.splfree();	
	}
	return(out1);
}

// Another way ;
struct spl3ds{
	NumericMatrix x; NumericVector y;
	
	double spl3dEval(NumericVector xnew){
		NumericVector x1 = unique(x(_,0)), x2 = unique(x(_,1)), x3 = unique(x(_,2));
		std::sort(x1.begin(), x1.end()); 
		std::sort(x2.begin(), x2.end()); 
		std::sort(x3.begin(), x3.end()); 
		int N1 = x1.length(), N2 = x2.length(), N3 = x3.length();
		int x3_index = findInterval1(xnew(2), x3);		//Find which discrete value it is in the 3rd dimension; 
		
		// Interpolating the 1st value conditional on the rest; 
		NumericVector y2(N2);
		for(int i=0; i<N2; i++){
			IntegerVector index = (x3_index -1)*N1*N2 + i*N1 + seq_len(N1) - 1;
			NumericVector tmpy = y[index];
			struct spl myspl = spl_init(x1, tmpy);
			y2(i) = myspl.splEval(xnew(0));
			myspl.splfree();
		}
		
		// Interpolate the 2nd value ; 
		struct spl spl2 = spl_init(x2, y2);
		double out = spl2.splEval(xnew(1));
		spl2.splfree();
		return(out);
	}
};

// [[Rcpp::export]]
double my_fnS(NumericMatrix x, NumericVector y, NumericVector xnew){
	gsl_set_error_handler_off();
	struct spl3ds out = {x, y};
	double a = out.spl3dEval(xnew);
	return(a);
}

// [[Rcpp::export]]
NumericVector my_fnS1(NumericMatrix x, NumericVector y, NumericMatrix xnew){
	gsl_set_error_handler_off();
	int n = xnew.nrow();
	NumericVector out1(n);
	struct spl3ds out = {x, y};
	for(int i=0; i<n; i++){
		out1(i) = out.spl3dEval(xnew(i,_));
	}
	return(out1);
}

