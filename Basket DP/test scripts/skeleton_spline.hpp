#infndef __UTILITIES__
#define __UTILITIES__

struct spl{
	gsl_spline *spline; gsl_interp_accel *acc; 
	NumericVector x; NumericVector y; double slope_first; double slope_end; 
}

struct spl spl_init(NumericVector x, NumericVector y);

double splineC(NumericVector x, NumericVector y, double xnew);

#endif	//__UTILITIES__