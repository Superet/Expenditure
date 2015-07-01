#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_spline.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]

//------------------//;
// Basic functions 	//;
//------------------//;
// [[Rcpp::export]]
int findInterval1(double x, NumericVector breaks) {
	int out;
	out = std::distance(breaks.begin(), std::upper_bound(breaks.begin(), breaks.end(), x));
	return out;
}

int my_rmultinorm(NumericVector probs) {
    int k = probs.size();
    IntegerVector ans(k);
    rmultinom(1, probs.begin(), k, ans.begin());
    IntegerVector index = seq_len(k);
    IntegerVector out = index[ans==1];
    return(out[0]);
}

// [[Rcpp::export]]
List param_assign(NumericVector param){
	double lambda = param[0], tau = param[2], mu = param[3], sigma = exp(param[4]); 
	return Rcpp::List::create(_["lambda"] = lambda, _["tau"] = tau, _["mu"] = mu, _["sigma"] = sigma);
}

//------------------//;
// Define structure //;
//------------------//;
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

struct my_f_params {NumericVector param_vec; NumericVector state_vec; 
					double beta; Function Q_fn; Function omega_fn; double lnzero; 
					NumericVector Inc_nodes; NumericMatrix Inc_weight; 
					NumericVector lnZ_nodes; NumericVector lnZ_weight; 
					NumericVector x1; NumericVector x2; NumericVector x3; NumericVector x4; struct spl *spl2d;};

struct mymin{
// s is the NM minimizer, the function F is negative of value function. The parameters need to be converted when computing value function; 
// s1d is brent one-dimensional minimizer when y = 0, the function F1d is negative value function; 
	gsl_multimin_fminimizer *s; 	gsl_multimin_function F; gsl_min_fminimizer *s1d; gsl_function F1d; 
	gsl_multimin_function orig_F; 	double bound_eps; 
	double NM_sizetol; 				double NM_startsize;	int NM_iter; 
	double brent_tol; 				int brent_iter; 
		
	NumericVector inner_multimin(gsl_vector *x){
		int iter = 0, status, n=F.n; 
		double size; 
		gsl_vector *ss; 
		ss	= gsl_vector_alloc(F.n); 
		
		/* Set initial step sizes to 1 */
		gsl_vector_set_all (ss, NM_startsize);
		gsl_multimin_fminimizer_set (s, &F, x, ss);
		iter = 0; 
		status = 0; 
		do{
		  iter++;
		  status = gsl_multimin_fminimizer_iterate(s);
		  if (status) break;
		  size = gsl_multimin_fminimizer_size (s);
		  status = gsl_multimin_test_size (size, NM_sizetol);
		}
		while (status == GSL_CONTINUE && iter < NM_iter);
		
		// Return the minimizer and optimal value; 
		gsl_vector_free(ss);  
		NumericVector out(F.n+1); 
		for(int i=0; i<n; i++){
			out(i) = gsl_vector_get(s->x, i); 
		}
		out(n) = s->fval; 
		return(out); 
	}
	
	NumericVector opt1d(){
		// When y = 0; 
		struct my_f_params *params 	= (struct my_f_params *)orig_F.params;
		NumericVector s1			= (params->state_vec); 
		NumericVector x1			= (params->x1);			
		double lnzero				= (params->lnzero);
		double 	lower_c 			= std::max(0.0, s1(0)+0.0-x1[x1.length()-1]), 
				upper_c 			= s1(0);
		double m = 0, v= 0, low = lower_c, up = upper_c;
		int sat, status;
		
		int iter	= 0;
		// Try initial values close to the boundaries;
		m		= upper_c - brent_tol;
		sat = gsl_min_fminimizer_set ( s1d, &F1d, m, low, up );
		if(sat == GSL_EINVAL){
			m	= lower_c + brent_tol;
			sat = gsl_min_fminimizer_set ( s1d, &F1d, m, low, up );
		}
		if(sat == GSL_EINVAL){
			m 	=  up;
			v 	=  F1d.function(up, F1d.params);	 
		}else{
			//Optimization loop;
			do {
				iter++;
				status = gsl_min_fminimizer_iterate(s1d);
				m      = gsl_min_fminimizer_x_minimum(s1d);
				low    = gsl_min_fminimizer_x_lower(s1d);
				up     = gsl_min_fminimizer_x_upper(s1d);
				v	   = gsl_min_fminimizer_f_minimum(s1d);

				status = gsl_min_test_interval(low,up, brent_tol, 0.0);
			}	
			while (status == GSL_CONTINUE && iter < brent_iter);
		}
		NumericVector out(F.n+1); 
		out(0) = lnzero; 
		out(1) = m; 
		out(2) = v; 
		return(out); 	
	}	
			
	NumericVector multmin_solver(gsl_vector *x0) {
		struct my_f_params *params 	= (struct my_f_params *)orig_F.params;
		NumericVector s1			= (params->state_vec);
		Function Q_fn				= (params->Q_fn);
		NumericVector Inc_nodes 	= (params->Inc_nodes);	
		double lnzero				= (params->lnzero);
		NumericVector x1			= (params->x1);
		int n 						= orig_F.n; 
		double Inc 					= Inc_nodes[s1[1] - 1], Imax = x1(x1.length()-1); 
		
		double upper_y = Inc- exp(s1(3)); 
		double upper_Q = as<double>(Q_fn(upper_y));
		double upper_c = s1(0) + std::max(upper_Q, 0.0); 
		NumericVector out(n+1); 
		gsl_vector *x; 
		x = gsl_vector_alloc(n);
		
		if(upper_y <=0 && upper_c<=0){
			out.fill(0); 
			gsl_vector_set_all(x, 0); 
			out(n)	= orig_F.f(x, orig_F.params); 
		}else if (upper_y <= 0 && s1(0)>0 ){
			out		= opt1d(); 
			out(0)  = 0; 
			out(n) 	= -out(n); 
		}else{
			// Set up initial values at the first guess and boundary; 
			out(n) = 1e5; 
			gsl_vector *b; 
			gsl_matrix *init; 
			init	= gsl_matrix_alloc (n, 3);
			b		= gsl_vector_alloc(n); 
			
			// Set initial values; 
			gsl_vector_set(b, 0, log(upper_y)); 
			gsl_vector_set(b, 1, log(upper_c)); 
			gsl_matrix_set_col(init, 0, x0); 
			gsl_matrix_set_col(init, 1, b); 
			gsl_matrix_set(init, 0, 2, gsl_vector_get(b,0));
			gsl_matrix_set(init, 1, 2, gsl_vector_get(x0,1));
			
			// Solve for the minimizer for each initial values; 
			for(int i=0; i<3; i++){
				gsl_matrix_get_col(x0, init, i);
				NumericVector out1 = inner_multimin(x0);
				if(out1(n) < out(n)){
					out = clone(out1); 
				}
			}
				
			// Check boundary and convert back the variables; 		
			if(exp(out[0]) < bound_eps){
				out[0] = 0; 
			}else if(upper_y- exp(out[0]) < bound_eps){
				out[0] = upper_y; 	
			}else{
				out[0] = exp(out[0]); 
			}
	
			double Q = as<double>(Q_fn(out[0]));
			double lower_c = std::max(0.0, s1[0] + Q - Imax), upper_c1 = s1[0] + Q; 
			if(exp(out[1]) < bound_eps){
				out[1] = lower_c; 
			}else if(upper_c1 - lower_c - exp(out[1]) < bound_eps){
				out[1] = upper_c1; 
			}else{
				out[1] = lower_c + exp(out[1]); 
			}
	
			gsl_vector_set(x, 0, out[0]); 
			gsl_vector_set(x, 1, out[1]); 
			out(n) = orig_F.f(x, orig_F.params);	
			
			gsl_vector_free(b);
			gsl_matrix_free(init);
		}	
		gsl_vector_free(x);

		return(out);
	}
};
					
struct policy_fn{
	int n_act; NumericVector Inc_nodes; NumericMatrix Inc_weight; Function Q_fn;
	NumericVector x1; NumericVector x2; NumericVector x3; NumericVector x4; 
	struct spl *spl2d_y; struct spl *spl2d_c; 
	
	NumericVector policyEval(NumericVector s1, double y){
		int N1 = x1.length(), N2 = x2.length(), N3 = x3.length(), N4 = x4.length();
		NumericVector out(n_act + 1), y4(N4), c4(N4);
		double Q; 
		int x2_index = findInterval1(s1(1), x2);
		int x3_index = findInterval1(s1(2), x3);
		
		if(y == -999){
			// Interpolate I_next conditional on the 4th dimension; 			
			for(int j=0; j<N4; j++){
				int index = j*N2*N3 + (x3_index - 1)*N2 + x2_index - 1; 
				y4(j) = spl2d_y[index].splEval(s1(0));
				c4(j) = spl2d_c[index].splEval(s1(0));
			}
		
			// Create the spline interpolation along x4; 			
			struct spl y_spl_new = spl_init(x4, y4), c_spl_new = spl_init(x4, c4);
			double ystar	= y_spl_new.splEval(s1(3)); 
			double cstar	= c_spl_new.splEval(s1(3)); 
			y_spl_new.splfree();
			c_spl_new.splfree();
		
			// Adjust the interior solution to satisfy boundary conditions; 
			double Inc = Inc_nodes[s1(1) - 1];
			double upper_y = std::max(0.0, Inc - exp(s1(3)) );
			out(0) = std::max(0.0, std::min(ystar, upper_y));
			Q = as<double>(Q_fn(out(0)) );
			double lower_c = std::max(0.0, s1(0) + Q - x1[N1-1] );
			out(1) = std::max(lower_c, std::min(cstar, s1(0) + Q));
		}else{
			// Interpolate I_next conditional on the 4th dimension; 			
			for(int j=0; j<N4; j++){
				int index = j*N2*N3 + (x3_index - 1)*N2 + x2_index - 1; 
				c4(j) = spl2d_c[index].splEval(s1(0));
			}
		
			// Create the spline interpolation along x4; 			
			struct spl c_spl_new = spl_init(x4, c4);
			double cstar	= c_spl_new.splEval(s1(3)); 
			c_spl_new.splfree();
		
			// Adjust the interior solution to satisfy boundary conditions; 
			out(0) = y;
			Q = as<double>(Q_fn(out(0)) );
			double lower_c = std::max(0.0, s1(0) + Q - x1[N1-1] );
			out(1) = std::max(lower_c, std::min(cstar, s1(0) + Q));
		}
		out(2) = s1(0) + Q - out(1); 
		return(out);
	} 
}; 					

//-----------------------//;
// Define value function //;
//-----------------------//;
double my_valuefn(const gsl_vector *a, void *p){
	struct my_f_params *params 	= (struct my_f_params *)p;
	NumericVector param_vec		= (params->param_vec);	
	NumericVector state_vec		= (params->state_vec);	
	double beta					= (params->beta);	
	Function Q_fn				= (params->Q_fn);	
	Function omega_fn			= (params->omega_fn);
	NumericVector Inc_nodes 	= (params->Inc_nodes);	
	NumericMatrix Inc_weight	= (params->Inc_weight);	
	NumericVector lnZ_nodes		= (params->lnZ_nodes);	
	NumericVector lnZ_weight	= (params->lnZ_weight);
	NumericVector x1			= (params->x1);
	NumericVector x2			= (params->x2);
	NumericVector x3			= (params->x3);
	NumericVector x4			= (params->x4);	
	struct spl *spl2d			= (params->spl2d);
	
	// Assign parameter; 
	List param_out 				= param_assign(param_vec);
	double lambda 				= param_out["lambda"], tau = param_out["tau"], mu = param_out["mu"], sigma = param_out["sigma"];
	double I 					= state_vec[0], Z = exp(state_vec[3]); 
	int Inc_index 				= state_vec[1], Rc = state_vec[2];
	double Inc 					= Inc_nodes[Inc_index - 1], y = gsl_vector_get(a, 0), c = gsl_vector_get(a, 1);
	
	// Initiate working variables; 
	int n_Inc = Inc_nodes.length(), n_Z = lnZ_nodes.length(), 
		N1 = x1.length(), N2 = x2.length(), N3 = x3.length(), N4 = x4.length(); 
	double out = 0, ev = 0; 
	NumericVector s_next(4), y4(N4), v(n_Z), ev_Inc(n_Inc), lnZ_next = sqrt(2)*sigma*lnZ_nodes + mu; 
	gsl_set_error_handler_off();

	// Transition of the next state; 
	double Q = as<double>(Q_fn(y));
	s_next[0] = I + Q - c; 
	s_next[1] = 0; 
	s_next[2] = Rc;
	s_next[3] = lnZ_nodes[0];
	NumericVector Inc_weight1 = Inc_weight((n_Inc*Rc + Inc_index - 1),_);

	// Compute the expected option value; 
	for(int i=0; i<n_Inc; i++){
		s_next[1] = i + 1; 
		
		// Interpolate the value function conditional on s2, s3, and x4;
		int x2_index = findInterval1(s_next(1), x2);
		int x3_index = findInterval1(s_next(2), x3);
		for(int j=0; j<N4; j++){
			int index = j*N2*N3 + (x3_index - 1)*N2 + x2_index - 1; 
			y4(j) = spl2d[index].splEval(s_next(0));
		}
		
		// Create the spline interpolation along x4; 			
		struct spl myspl = spl_init(x4, y4);
		for(int j=0; j<n_Z; j++){
			v[j] = myspl.splEval(lnZ_next[j]);
		}
		myspl.splfree();
		ev_Inc(i) = sum(lnZ_weight*v)/sqrt(PI);				
	}
	ev = sum(ev_Inc * Inc_weight1); 

	// Compute the value function;
	double omega = as<double>(omega_fn(y));
	out = omega + lambda*log(c+0.01) + (Inc-y-Z) - tau * s_next[0] + beta*ev; 
	return(out);
}	

double my_valuefn_t(const gsl_vector *a_t, void *p){
	struct my_f_params *params 	= (struct my_f_params *)p;
	NumericVector state_vec		= (params->state_vec);
	Function Q_fn				= (params->Q_fn);
	double lnzero				= (params->lnzero);	
	NumericVector Inc_nodes		= (params->Inc_nodes);	
	NumericVector x1			= (params->x1);
	double I					= state_vec[0], Z = exp(state_vec[3]); 
	int Inc_index 				= state_vec[1], N1 = x1.length();
	double Inc 					= Inc_nodes[Inc_index - 1]; 

	// Variable transformation for the control variables and construct constraints; 
	double y = 0, c= 0;
	if( gsl_vector_get(a_t, 0) > lnzero){
		y = exp(gsl_vector_get(a_t, 0));
	}
	double Q = as<double>(Q_fn(y));
	double lowerc = std::max(0.0, I + Q - x1[N1-1]);
	if( gsl_vector_get(a_t, 1) > lnzero){
		c = lowerc + exp(gsl_vector_get(a_t, 1)); 
	}else{
		c = lowerc;
	}
	gsl_vector *a; 
	a = gsl_vector_alloc(2);
	gsl_vector_set(a, 0, y);
	gsl_vector_set(a, 1, c);
	
	// Set the constraints: g(x) <= 0; 
	NumericVector g(2), sign_g(2);
	g(0) = y - (Inc - Z)*1*(Inc>Z); 
	g(1) = c - (I + Q); 
	double my_penalty = 10000.0, v, out, tol_g = 1e-6; 
	sign_g.fill(0); 
	if(g(0) > tol_g) {sign_g(0) = 1; }
	if(g(1) > tol_g) {sign_g(1) = 1; }
	
	if(sum(sign_g)==0){
		v = my_valuefn(a, p);
		out = v; 
	}else{
// 		out = v + my_penalty * sum(g*sign_g);
// 		out = -my_penalty * sum(g*sign_g);
		out = -my_penalty;
	}
	gsl_vector_free(a);
	return(-out);
}

double my_value_1c(double c, void *p){	
	gsl_vector *a; 
	a = gsl_vector_alloc(2); 
	gsl_vector_set(a, 0, 0); 
	gsl_vector_set(a, 1, c); 
	double out = my_valuefn(a, p); 
	gsl_vector_free(a); 
	return(-out); 
}

// [[Rcpp::export]]
NumericVector my_fn(NumericMatrix policy, NumericVector V, List DP_list, List control_list){
	NumericVector param = DP_list["param"]; 
	NumericMatrix state = DP_list["state"]; 
	double beta = DP_list["beta"]; 
	Function omega_fn = DP_list["omega"], Q_fn = DP_list["Q_fn"]; 
	List Inc_kernel = DP_list["Inc_kernel"], lnZ_kernel = DP_list["lnZ_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"], 
				lnZ_nodes = lnZ_kernel["nodes"], lnZ_weight = lnZ_kernel["weight"]; 
	NumericMatrix Inc_weight = Inc_kernel["weight"]; 
	int ns = state.nrow();
	double lnzero = control_list["lnzero"]; 
	
	// Sort out the grids for each dimension in the states; 
	NumericVector x1 = unique(state(_,0)), x2 = unique(state(_,1)), x3 = unique(state(_,2)), x4 = unique(state(_,3));
	int N1 = x1.length(), N2 = x2.length(), N3 = x3.length(), N4 = x4.length();
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
	std::sort(x3.begin(), x3.end()); 
	std::sort(x4.begin(), x4.end()); 
	
	// Initiate the c-spline; 
	int nspl = N2*N3*N4;
	struct spl myspl[nspl];
	for(int i=0; i<nspl; i++){
		IntegerVector index	= i*N1 + seq_len(N1) - 1;
		NumericVector tmpy 	= V[index];
		myspl[i] = spl_init(x1, tmpy);
	}
	
	// Assign parameters for the value function; 
	struct my_f_params params = {param, state(0,_),beta, Q_fn, omega_fn, lnzero, 
								Inc_nodes, Inc_weight, lnZ_nodes, lnZ_weight, 
								 x1, x2, x3, x4, myspl};	
	gsl_vector *a; 
	a = gsl_vector_alloc(2); 
	NumericVector out(ns);

	double tmpy, tmpc, lower_c, Q; 
	for(int i=0; i<ns; i++){
		params.state_vec = state(i,_); 
		gsl_vector_set(a, 0, policy(i,0)); 
		gsl_vector_set(a, 1, policy(i,1));
// 		Rcout<<"state "<<i<<", y="<<tmpy<<", c="<<tmpc<<"\n";
		out(i) = my_valuefn(a, &params) ;
	}	  
	gsl_vector_free(a);
	return(out);
}

//-------------------------------------------------//;
// Solving Dynamic programming -- policy iteration //;
//-------------------------------------------------//;
// [[Rcpp::export]]
List policy_improveC(NumericVector V_old, List DP_list, List control_list, double NM_sizetol=0, int NM_iter = 0){
	// Unload parameters; 
	NumericVector param = DP_list["param"]; 
	NumericMatrix state = DP_list["state"]; 
	double beta 		= DP_list["beta"]; 
	Function omega_fn 	= DP_list["omega"], Q_fn = DP_list["Q_fn"]; 
	List Inc_kernel 	= DP_list["Inc_kernel"], lnZ_kernel = DP_list["lnZ_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"], 
			lnZ_nodes	= lnZ_kernel["nodes"], lnZ_weight = lnZ_kernel["weight"]; 
	NumericMatrix Inc_weight = Inc_kernel["weight"]; 
	int ns 				= state.nrow(), n_act = 2;
	double lnzero = control_list["lnzero"], Qk = DP_list["k"];
	if(NM_sizetol ==0) 	{NM_sizetol = control_list["NM_sizetol"]; }
	if(NM_iter ==0) 	{NM_iter = control_list["NM_iter"];}
	double NM_startsize = control_list["NM_startsize"], brent_tol = control_list["brent_tol"], bound_eps = control_list["bound_eps"]; 
	int brent_iter 		= control_list["brent_iter"]; 
	gsl_set_error_handler_off();
	
	// Sort out the grids for each dimension in the states; 
	NumericVector x1 = unique(state(_,0)), x2 = unique(state(_,1)), x3 = unique(state(_,2)), x4 = unique(state(_,3));
	int N1 = x1.length(), N2 = x2.length(), N3 = x3.length(), N4 = x4.length();
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
	std::sort(x3.begin(), x3.end()); 
	std::sort(x4.begin(), x4.end()); 
	
	// Initiate the c-spline; 
	int nspl = N2*N3*N4;
	struct spl myspl[nspl];
	for(int i=0; i<nspl; i++){
		IntegerVector index	= i*N1 + seq_len(N1) - 1;
		NumericVector tmpy 	= V_old[index];
		myspl[i] = spl_init(x1, tmpy);
	}
	
	// Assign parameters for the value function; 
	struct my_f_params params = {param, state(0,_),beta, Q_fn, omega_fn, lnzero, 
								Inc_nodes, Inc_weight, lnZ_nodes, lnZ_weight, 
								 x1, x2, x3, x4, myspl};	
	
	// Set up multimin; 
	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s;
	s = gsl_multimin_fminimizer_alloc (T, n_act);
	gsl_vector *m; 
	m = gsl_vector_alloc (n_act);
	
	// GSL objective function ;
	gsl_multimin_function F; 
	F.n			= n_act; 
	F.f 		= &my_valuefn_t;
	F.params	= &params;
	
	// Claim the original function; 
	gsl_multimin_function orig_func;
  	orig_func.n 	= n_act;
	orig_func.f 	= &my_valuefn;
	orig_func.params= &params;
	
	// Claim the 1d function when y=0; 
	const gsl_min_fminimizer_type *T1d;
	gsl_min_fminimizer *s1d;
	T1d = gsl_min_fminimizer_brent;
	s1d = gsl_min_fminimizer_alloc(T1d);
	gsl_function F1d; 
	F1d.function 	= &my_value_1c;
	F1d.params		= &params;
		
	// Initiate constrained optimization; 	
	struct mymin my_multmin = {s, F, s1d, F1d, orig_func, bound_eps, 
							   NM_sizetol, NM_startsize, NM_iter, brent_tol, brent_iter};  
	
	// Initiate working variables;
	NumericMatrix policy(ns, n_act); 
	NumericVector V(ns), min_f(n_act+1); 
	
	// Loop over each state to solve the maximization; 
	for(int i=0; i<ns; i++){
		params.state_vec = state(i,_); 
		
		// Starting values; 
		double Inc 		= Inc_nodes[state(i,1) - 1]; 
		double upper_y 	= Inc - exp(state(i,3)); 
		if(upper_y <= 0){
			gsl_vector_set(m, 0, lnzero); 
			gsl_vector_set(m, 1, std::max(lnzero, log(state(i,0))) ); 
		}else{
			gsl_vector_set(m, 0, log(.5*upper_y) );
			gsl_vector_set(m, 1, log(.5*(state(i,0) + .05*upper_y))) ;
		}	 
		
		min_f = my_multmin.multmin_solver(m);
		for(int j=0; j<n_act; j++){
			policy(i,j) = min_f(j); 
		}
		V(i)		= min_f(n_act);
	}
	gsl_vector_free(m);
	gsl_multimin_fminimizer_free (s);
	gsl_min_fminimizer_free(s1d);
	for(int i=0; i<nspl; i++){
		myspl[i].splfree();
	}
	return Rcpp::List::create( _["policy"]=policy, _["V"] = V); 	
}

// [[Rcpp::export]]
NumericVector policy_evalC(NumericMatrix policy, NumericVector V_old, NumericVector V0, List DP_list, List control_list, int M=0){
	// Unload parameters; 
	NumericVector param = DP_list["param"]; 
	NumericMatrix state = DP_list["state"]; 
	double beta = DP_list["beta"]; 
	Function omega_fn = DP_list["omega"], Q_fn = DP_list["Q_fn"]; 
	List Inc_kernel = DP_list["Inc_kernel"], lnZ_kernel = DP_list["lnZ_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"], 
				lnZ_nodes = lnZ_kernel["nodes"], lnZ_weight = lnZ_kernel["weight"]; 
	NumericMatrix Inc_weight = Inc_kernel["weight"]; 
	int ns = state.nrow(), n_act = 2;
	
	// Unload the iteration parameter;
	double tol = control_list["inner_tol"];
	if(M==0){
		M = control_list["inner_max"]; 
	}
	double lnzero = control_list["lnzero"];
	gsl_set_error_handler_off();
	
	// Sort out the grids for each dimension in the states; 
	NumericVector x1 = unique(state(_,0)), x2 = unique(state(_,1)), x3 = unique(state(_,2)), x4 = unique(state(_,3));
	int N1 = x1.length(), N2 = x2.length(), N3 = x3.length(), N4 = x4.length();
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
	std::sort(x3.begin(), x3.end()); 
	std::sort(x4.begin(), x4.end()); 
	int nspl = N2*N3*N4;
	
	// Initiate the iteration;
	NumericVector V_new(ns), V_oldl(ns); 
	V_new 	= clone(V0);
	V_oldl	= clone(V_old);
	double dif = max(abs(V_new - V_oldl)); 
	int k = 0; 
	gsl_vector *a; 
	a = gsl_vector_alloc(n_act);
	
	if(dif > tol){
		V_oldl = clone(V0); 
		do{
			k++; 
			V_new.fill(0); 

			// Initiate the c-spline; 
			struct spl myspl[nspl];
			for(int i=0; i<nspl; i++){
				IntegerVector index	= i*N1 + seq_len(N1) - 1;
				NumericVector tmpy 	= V_oldl[index];
				myspl[i] = spl_init(x1, tmpy);
			}	
			// Assign parameters for the value function; 
			struct my_f_params params = {param, state(0,_),beta, Q_fn, omega_fn, lnzero, 
							Inc_nodes, Inc_weight, lnZ_nodes, lnZ_weight, 
							 x1, x2, x3, x4, myspl};	
			
			// Compute the value function at each states; 
			for(int i=0; i<ns; i++){
				params.state_vec = state(i,_); 
				gsl_vector_set(a, 0, policy(i,0) ); 
				gsl_vector_set(a, 1, policy(i,1) );
				V_new(i) = my_valuefn(a, &params) ;
			}
			
			// Check distance; 
			dif = max(abs(V_new - V_oldl));
			V_oldl = clone(V_new); 
			// Free spline; 
			for(int i=0; i<nspl; i++){
				myspl[i].splfree();
			}			
		}
		while(k <= M && dif > tol);	
	}
	gsl_vector_free(a);
	return(V_new); 
}

// [[Rcpp::export]]
List policy_iterC(List DP_list, List control_list, int print_level=0){
	double tol 		= control_list["tol"];
	int max_iter	= control_list["max_iter"];
	NumericVector V0 = DP_list["value_fn"];
	int ns 			= V0.length(), n_act = 2;
	NumericVector V_new(ns); 
	NumericMatrix policy_new(ns, n_act);
	List DP_listout = clone(DP_list);
	
	// Compute the policy and value function at the initial guess of value function; 
	List tmp_sol = policy_improveC(V0, DP_list, control_list);
	NumericMatrix policy_old 	= tmp_sol["policy"]; 
	NumericVector V_old			= tmp_sol["V"]; 
	
	// Iterate policy function; 
	double dif = tol + 1, NM_sizetol; 
	int iter = 0, Ml, NM_iter;
	int M = control_list["inner_max"];  	
	do{
		iter++; 
		if(iter <=10 ) {
			Ml = .5*M;  
			NM_iter = 150; 
			NM_sizetol = 1e-4; 
		}else{ 
			Ml = M ; 
			NM_iter = 0; 
			NM_sizetol = 0; 
		}
		
		// Policy improvement step; 
		tmp_sol		= policy_improveC(V_old, DP_list, control_list, NM_sizetol, NM_iter);
		NumericMatrix policy1	= tmp_sol["policy"]; 
		policy_new	= clone(policy1); 
		V0 			= tmp_sol["V"]; 
		
		// Policy evaluation step;
		V_new 		= policy_evalC(policy_new, V_old, V0, DP_list, control_list, Ml); 
		dif			= max(abs(policy_new - policy_old));
		V_old		= clone(V_new); 
		policy_old	= clone(policy_new); 
		
		// Print the improvement; 
		if(print_level > 0){
			Rcout<<"Iteration "<<iter<<", norm_diff = "<<dif<<"\n";
		}
	}
	while( iter<=max_iter && dif > tol );
	
	// Return the policy and value function;
	DP_listout["value_fn"]	= policy_evalC(policy_new, V_new, V0, DP_list, control_list, 500);
	DP_listout["policy"]	= policy_new; 
	DP_listout["Iteration"] = iter;
	if(dif > tol){
		DP_listout["status"] = 2;
// 		if(print_level == 1){
			Rcout<<"Policy function iteration did not converge.\n";
// 		}
	}else{
		DP_listout["status"] = 0;	
	}
	return(DP_listout);
}

//------------------------------//;
// Simulate behavioral sequence //;
//------------------------------//;
List sim_seq1C(NumericVector init_state, int TT, NumericMatrix DataState, NumericMatrix choice_seq, struct policy_fn myplc, bool drawall=false){
	RNGScope scope;	
// 	struct policy_fn myplc = (struct policy_fn *)plcspl;
	NumericMatrix Inc_weight = myplc.Inc_weight; 
	int n_act = myplc.n_act;
	int nsd = init_state.length(), n_Inc = Inc_weight.ncol(); 
	NumericMatrix state_out(TT, nsd), choice_out(TT, n_act); 
	state_out	= clone(DataState); 
	choice_out 	= clone(choice_seq); 
	
// 	if(drawall){
// 		NumericVector lnZ_draw = rnorm(TT, mu, sigma);
// 		state_out(_,3) = lnZ_draw;
// 		state_out(_,2).fill(init_state(2));
// 	}
	for(int j=0; j<nsd; j++){
		state_out(0,j) = init_state(j);
	}
	
	// Simulate the policy; 
	for(int i=0; i<TT; i++){
		NumericVector out	= myplc.policyEval(state_out(i,_), choice_out(i,0)); 
		choice_out(i, 0) 	= out(0); 
		choice_out(i, 1)	= out(1); 
		if(i < TT-1){
			state_out(i+1,0)	= out(2); 
			if(drawall){
				state_out(i+1,1) = my_rmultinorm(Inc_weight(n_Inc*state_out(i,2) + state_out(i,1) -1 ,_));
			}
		}
	}
	return Rcpp::List::create(_["DataState"] = state_out, _["choice_seq"] = choice_out);
}

// [[Rcpp::export]]
List sim_hhseqC(NumericVector hh_index, IntegerVector TT_vec, NumericMatrix init_state, NumericMatrix DataState, NumericMatrix choice_seq, List DP_list, List control_list, bool scalez=true){
	// Assign parameters; 
	NumericVector param = DP_list["param"]; 
	List param_out 		= param_assign(param); 
	double mu 			= param_out["mu"], sigma = param_out["sigma"]; 
	NumericMatrix state = DP_list["state"]; 
	Function omega_fn 	= DP_list["omega"], Q_fn = DP_list["Q_fn"]; 
	List Inc_kernel 	= DP_list["Inc_kernel"], lnZ_kernel = DP_list["lnZ_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"], 
			lnZ_nodes 	= lnZ_kernel["nodes"], lnZ_weight = lnZ_kernel["weight"]; 
	NumericMatrix Inc_weight = Inc_kernel["weight"], policy	= DP_list["policy"];
	NumericVector policy_y = policy(_,0), policy_c = policy(_,1); 
	int n_act 			= 2, nsd = state.ncol();
	gsl_set_error_handler_off();
	
	// Initiate the working variables; 
	int nhh = hh_index.length(), N = choice_seq.nrow(); 
	NumericMatrix state_out(N, nsd), choice_out(N, n_act); 
	state_out 	= clone(DataState);
	if(scalez){
		state_out(_,3) = mu + sigma*DataState(_,3);
	}
	for(int i=0; i<n_act -1; i++){
		choice_out(_,i) = choice_seq(_,i); 
	}
		
	// Construct policy interpolation; 
	NumericVector x1 = unique(state(_,0)), x2 = unique(state(_,1)), x3 = unique(state(_,2)), x4 = unique(state(_,3));
	int N1 = x1.length(), N2 = x2.length(), N3 = x3.length(), N4 = x4.length();
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
	std::sort(x3.begin(), x3.end()); 
	std::sort(x4.begin(), x4.end()); 
	
	// Initiate the c-spline; 
	int nspl = N2*N3*N4;
	struct spl y_spl[nspl], c_spl[nspl];
	for(int i=0; i<nspl; i++){
		IntegerVector index	= i*N1 + seq_len(N1) - 1;
		NumericVector tmpy 	= policy_y[index];
		NumericVector tmpc 	= policy_c[index];
		y_spl[i] = spl_init(x1, tmpy);
		c_spl[i] = spl_init(x1, tmpc);
	}
	
	// Initiate policy function structure; 
	struct policy_fn myplc = {n_act, Inc_nodes, Inc_weight, Q_fn, x1, x2, x3, x4, y_spl, c_spl};
	
	// Loop over households; 
	int cnt = 0; 
	for(int h=0; h<nhh; h++){
		// Assign initial states; 
		for(int j=0; j<nsd; j++){
			state_out(cnt,j) = init_state(h,j);
		}		
		for(int i=cnt; i<cnt+TT_vec(h); i++){
			NumericVector out	= myplc.policyEval(state_out(i,_), choice_out(i,0)); 
			choice_out(i, 0) 	= out(0); 
			choice_out(i, 1)	= out(1); 
			if(i-cnt < TT_vec(h)-1){
				state_out(i+1,0)	= out(2); 
			}
		}
		cnt = cnt + TT_vec[h]; 
	}
	for(int i=0; i<nspl; i++){
		y_spl[i].splfree();
		c_spl[i].splfree();
	}
	return Rcpp::List::create(_["DataState"] = state_out, _["choice_seq"] = choice_out);
}

//------------//;
// Estimation //;
//------------//;
// [[Rcpp::export]]
NumericVector MM1C(NumericVector param, NumericMatrix choice_obs, NumericMatrix choice_seq, NumericMatrix DataState, IntegerVector TT_vec, List DP_list){
	// Assign parameter; 
	Function omega_j	= DP_list["omega_j"], Q_j = DP_list["Q_j"], Q_fn = DP_list["Q_fn"]; 
	double beta 		= DP_list["beta"];
	NumericMatrix state	= DP_list["state"]; 
	NumericVector x1	= state(_,0); 
	List Inc_kernel 	= DP_list["Inc_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"]; 
	List param_out 		= param_assign(param); 
	double lambda 		= param_out["lambda"], tau = param_out["tau"], mu = param_out["mu"], sigma = param_out["sigma"];
	NumericVector y 	= choice_seq(_,0), c = choice_seq(_,1); 
	int nhh 			= TT_vec.length(), N = DataState.nrow(), N1 = x1.length();
	std::sort(x1.begin(), x1.end());
	
	// Construct instrument matrix; 
// 	int n_inst = 2; 
// 	NumericMatrix inst(N-nhh, n_inst); 
// 	NumericVector Inc(N); 
// 	for(int i=0; i<N; i++) {Inc(i) = Inc_nodes(DataState(i,1) - 1); }
// 	double y_low 			= 0.0; 
// 	double domega_low		= as<double>(omega_j(0)); 
// 	double dQ_low			= as<double>(Q_j(0)); 
// 	NumericVector Q			= Q_fn(y), 
// 				  y_up		= pmax(0.0, Inc - exp(DataState(_,3))); 
// 	NumericVector c_low		= pmax(0.0, DataState(_,0) + Q - x1[N1-1]), 
// 				  c_up		= DataState(_,0) + Q, 
// 				  domega	= omega_j(y), 
// 				  domega_up	= omega_j(y_up), 
// 				  dQ 		= Q_j(y), 
// 				  dQ_up		= Q_j(y_up); 
// 	double dvy, dvy_low, dvy_up, dvc, dvc_low, dvc_up; 
// 	NumericVector my(N-nhh), mc(N-nhh); 
// 	int cnt = 0, cnt1 = 0; 
// 	for(int h=0; h<nhh; h++){
// 		for(int i=0; i<TT_vec(h) -1; i++){
// 			int j = cnt + i; 
// 			dvy				= domega(j) - 1- tau*dQ(j) + beta - beta*domega(j+1);
// 			dvy_low			= domega_low -1 - tau*dQ_low + beta - beta*domega(j+1);
// 			dvy_up			= domega_up(j) -1 - tau*dQ_up(j) + beta - beta*domega(j+1); 
// 			my(cnt1+i)		= dvy - std::min(dvy_low, std::max(0.0, dvy_up)); 
// 			
// 			dvc				= lambda/(c(j)+.01) + tau - beta*lambda/(c(j+1)+.01); 
// 			dvc_low			= lambda/(c_low(j)+.01) + tau - beta*lambda/(c(j+1)+.01); 
// 			dvc_up			= lambda/(c_up(j)+.01) + tau - beta*lambda/(c(j+1)+.01); 
// 			dvc				= lambda/(c(j)+.01) + tau - beta*(domega(j+1)-1)/dQ(j+1); 
// 			dvc_low			= lambda/(c_low(j)+.01) + tau - beta*(domega(j+1)-1)/dQ(j+1); 
// 			dvc_up			= lambda/(c_up(j)+.01) + tau - beta*(domega(j+1)-1)/dQ(j+1); 
// 			mc(cnt1+i) 		= dvc - std::min(dvc_low, std::max(0.0, dvc_up)); 			
// 			inst(cnt1+i, 0) = 1; 
// 			inst(cnt1+i, 1) = DataState(j,1); 
// 		}
// 		cnt 	= cnt + TT_vec(h); 
// 		cnt1 	= cnt1 + TT_vec(h) -1; 
// 	}
	NumericVector my(N), mc(N);
	my	= choice_obs(_,0) - choice_seq(_,0); 
	mc	= pow(my, 2); 	
	int n_inst = 2; 
	NumericMatrix inst(N, n_inst); 
	for(int i=0; i<N; i++){
		inst(i,0) = 1; 
	}
	inst(_,1) = DataState(_,1); 
	
	// Stacking the moments; 
	int n_m = 2; 
	int n_M = n_m*n_inst ; 
	NumericMatrix M(N-nhh, n_M); 
	NumericVector M_bar(n_M); 
	for(int i=0; i<n_inst; i++){
		M(_,i) 			= my * inst(_,i); 
		M(_,i+n_inst)	= mc * inst(_,i); 
	}
	for(int i=0; i<n_M -1; i++){
		M_bar(i) = sum(M(_,i))/N; 
	}
// 	M_bar(n_M-2) = sum(pow(my,2))/N; 
	
	return(M_bar);
}

// [[Rcpp::export]]
double GMM_objC(NumericVector param, NumericVector hh_index, IntegerVector TT_vec, NumericMatrix init_state, NumericMatrix DataState, NumericMatrix choice_obs, NumericMatrix W, List DP_init, NumericMatrix lnz_draw, List control_list){
	// Find the DP solution; 
	List DP_listout		= clone(DP_init);
	DP_listout["param"]	= param; 
	
	// Check the boundary of parameters; 
	List param_out = param_assign(param); 
	double mu = param_out["mu"], sigma = param_out["sigma"]; 
	NumericMatrix state = DP_listout["state"]; 
	double lnZ_min = min(state(_,3)), lnZ_max = max(state(_,3));
	if(mu - 3*sigma <lnZ_min || mu + 3*sigma > lnZ_max){
		return(NA_REAL); 
	}	
	DP_listout			= policy_iterC(DP_listout, control_list);
	
	// Simulate the whole sequence of states and choice variables at one set of z draw; 
	int nz 	= lnz_draw.ncol(), n_M = 4; 
	NumericMatrix DataState1(DataState.nrow(), DataState.ncol()); 
	NumericMatrix choice_seq(DataState.nrow(), 2); 
	choice_seq.fill(-999);
	NumericVector M_bar(n_M); 
	M_bar.fill(0); 
	DataState1 = clone(DataState);
	for(int i=0; i<nz; i++){
		DataState1(_,3) = lnz_draw(_,i); 
		List seqout 		= sim_hhseqC(hh_index, TT_vec, init_state, DataState1, choice_seq, DP_listout, control_list); 
		NumericMatrix state_out	= seqout["DataState"];
		NumericMatrix choice_out= seqout["choice_seq"];
		NumericVector m_new 	= MM1C(param, choice_obs, choice_out, state_out, TT_vec, DP_listout);
		M_bar					= M_bar + m_new; 
// 		Rcout<<i<<" z draw, m="<<m_new[0]<<","<<m_new[1]<<","<<m_new[2]<<","<<m_new[3]<<","<<m_new[4]<<".\n"; 
	}
	M_bar = M_bar/nz; 

	// Return the objective function; 
	double Q = 0;
	for(int i=0; i<n_M; i++){
		for(int j=0; j<n_M; j++){
			Q = Q + M_bar(i)*W(i,j)*M_bar(j); 
		}
	} 
	return(-Q);
}
