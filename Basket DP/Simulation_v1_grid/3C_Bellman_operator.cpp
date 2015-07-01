#include <Rcpp.h>
#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_interp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]

struct my_f_params { double I_current; double lambda; double tau1; double tau3;
					double beta; 
					double q;
					const gsl_interp *myint; gsl_interp_accel *acc; 
					NumericVector x;NumericVector y;};

double my_valuefn( double cs, void *p){
    struct my_f_params *params = (struct my_f_params *)p;	
    double lambda				= (params->lambda);
	double tau1					= (params->tau1);
	double tau3					= (params->tau3);
	double beta					= (params->beta);
	double q					= (params->q);
	double I_current			= (params->I_current);
	const gsl_interp *myint    	= (params->myint);
    gsl_interp_accel *acc      	= (params->acc);
    NumericVector x           	= (params->x);
    NumericVector y            	= (params->y);
    int k = x.length();
    	
	double I_next=0, ev = 0, out = 0;
	I_next = I_current + q -cs;
	if(I_next>=x(0) & I_next<=x(k-1)){
		ev = gsl_interp_eval(myint, x.begin(), y.begin(), I_next, acc);
	}else if(I_next < x(0)){
		ev = y(0);
// 		ev = y(0) + (y(1)-y(0))/(x(1)-x(0)) * (I_next - x(0));
	}else {
		ev = y(k-1);
// 		ev = y(k-1) + (y(k-2)-y(k-1))/(x(k-2)-x(k-1)) * (I_next - x(k-1));
	}
	out = lambda*log(cs+.01) - tau1*I_next - tau3*q + beta*ev;
	return(-out);
}

// [[Rcpp::export]]
List Solve_c_fnC(NumericVector I_grid, NumericVector x, NumericVector V, NumericVector param, 
				double beta, NumericVector Q_grid){
	const double lambda = param(0), tau1 = param(1), tau3 = param(2);
	int ns = I_grid.length(), nQ = Q_grid.length();
	int k = x.length();
	double Imax = max(I_grid);
	NumericMatrix cstar(ns,nQ);
	NumericMatrix value(ns,nQ);
	NumericVector policy_k(ns), value_max(ns);
	
	// Linear interpolation set up;
	gsl_interp_accel *accP	= gsl_interp_accel_alloc() ;
	gsl_interp *linP 		= gsl_interp_alloc( gsl_interp_linear , k );
	gsl_interp_init( linP, x.begin(), V.begin(), k);
	
	struct my_f_params params = {I_grid(0), lambda, tau1, tau3, beta, Q_grid(0), linP, accP, x, V};
		
	// Optimizer setup;
	int status, sat;
	int iter=0,max_iter=200;
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	gsl_set_error_handler_off(); //TURN OFF GSL ERROR HANDLER (GSL JUST PRINT ERROR MESSAGE AND KILL THE PROGRAM IF FLAG IN ON)
	
	// GSL objective function ;
	gsl_function F; 
	F.function 	= &my_valuefn;
	F.params	= &params;
	
	// Create optimizer;
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	double low = x(0), up = x(k-1);
	double m = 0.2, v= 0;
			
	// Solve for c for every state in a loop;
	for(int i=0; i < ns; i++){
		for(int j=0; j < nQ; j++){
			params.I_current	= I_grid(i);
			params.q			= Q_grid(j);
			if(I_grid(i) + Q_grid(j)<=0){
				cstar(i,j) = 0;
				value(i,j) = -my_valuefn(0, &params); 
			}else{
				// Re-initiate minimizer and parameters;
				low			= 0;
				up			= I_grid(i) + Q_grid(j);
				iter		= 0;
				int m_iter	= 1;
				
				do{
					m = m_iter*0.1*(I_grid(i) + Q_grid(j));
					sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
					m_iter++;
				}
				while(sat == GSL_EINVAL && m_iter < 10);

				if(sat == GSL_EINVAL){
					cstar(i,j) =  up;
					value(i,j) =  -my_valuefn(cstar(i,j), &params);	 
				}else{
					//Optimization loop;
					do {
					iter++;
					status = gsl_min_fminimizer_iterate(s);
					m      = gsl_min_fminimizer_x_minimum(s);
					low    = gsl_min_fminimizer_x_lower(s);
					up     = gsl_min_fminimizer_x_upper(s);
					v	   = gsl_min_fminimizer_f_minimum(s);

					status = gsl_min_test_interval(low,up,0.0000001,0.0);
					}	
					while (status == GSL_CONTINUE && iter < max_iter);
					cstar(i,j) = m;
					value(i,j) = -v;
				}			
			}
			value_max(i)	= max(value(i,_));	
			for(int j= 0; j <nQ; j++){
				if(value(i,j)==value_max(i)) { policy_k(i) = j+1;}
			}
		}
	}
	gsl_min_fminimizer_free(s);
	gsl_interp_free(linP);
    gsl_interp_accel_free (accP);

	return Rcpp::List::create(_["policy"] = List::create(_["k"] = policy_k, _["c"] = cstar), 
							  _["value"]=value_max);
}

// [[Rcpp::export]]
List Bellman_operatorC(List DP_list){
	NumericVector I_grid = DP_list["state"], V = DP_list["value_fn"], Q_grid = DP_list["Q_grid"];
	double	beta = DP_list["beta"];
	List param_list = DP_list["param_list"];
	NumericVector param(3);
	param(0) = param_list["lambda"];
	param(1) = param_list["tau1"];
	param(2) = param_list["tau3"];
	List out;
	out  = Solve_c_fnC(I_grid, I_grid, V, param, beta, Q_grid);
	return(out);
}
