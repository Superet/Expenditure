#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
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

struct my_f_params {NumericVector param_vec; NumericVector state_vec; 
					double beta; Function Q_fn; Function omega_fn; 
					Function Q_j; Function omega_j; 
					double lnzero; 
					NumericVector Inc_nodes; NumericMatrix Inc_weight; 
					NumericVector lnZ_nodes; NumericVector lnZ_weight; 
					NumericVector x1; NumericVector x2; NumericVector x3; NumericVector x4; 
					struct spl *spl2d;};

struct constr_par {void *f_par; double *mu;}; 

struct pen_constrOptim{
	gsl_multimin_fdfminimizer *s; gsl_multimin_function_fdf F; gsl_multimin_function_fdf orig_F; 
	gsl_min_fminimizer *s1d; gsl_function F1d; 
	double outer_eps; int outer_iter; double mu_start; double mu_step; double mu_max; 
	double brent_tol; int brent_iter; double bound_eps; 
	
	double inner_multimin(gsl_vector *x, gsl_vector *par){
		int iter = 0, status, n=F.n; 
	
		gsl_multimin_fdfminimizer_set (s, &F, x, 0.001, 0.01);
		iter = 0; 
		status = 0; 
		status = gsl_multimin_fdfminimizer_iterate (s);

		do{
			iter++;
			status = gsl_multimin_fdfminimizer_iterate (s);
			if (status){	
				break; 
			}
			status = gsl_multimin_test_gradient (s->gradient, 1e-5);
		}while (status == GSL_CONTINUE && iter < 100);
		
		// Return the minimizer and optimal value; 
		gsl_vector_memcpy(par, s->x); 
		return(s->f );
	}
	
	NumericVector constr_optim(gsl_vector *x0){
		double obj, r, obj_old, r_old; 
		int n 				= F.n; 
		struct constr_par *par	= (struct constr_par *)F.params;
		double *mu_ptr 		= (par->mu); 
		double current_mu 	= mu_start; 
		*mu_ptr 			= mu_start; 
		
		// Initiate working variables; 
		NumericVector out(n+1); 
		gsl_vector *x; 
		x = gsl_vector_alloc(n); 
		gsl_vector_memcpy (x, x0);
		
		// Set the function values at the initial values; 
		obj			= orig_F.f(x0, orig_F.params); 
		r			= F.f(x0, F.params); 
		int check 	= 0; 
		int iter 	= 0; 
		do{
			iter++; 
			gsl_vector_memcpy (x0, x);
// 			obj_old 	= obj; 
			r_old		= r; 
			r = inner_multimin(x0, x); 
			if( gsl_finite(r) && gsl_finite(r_old) && abs(r - r_old)<(0.001 + abs(r))*outer_eps && current_mu >= mu_max){
// 			 
				check 	= 1; 
			}
// 			obj 		= orig_F.f(x, orig_F.params); 
			current_mu	= current_mu * mu_step; 
			*mu_ptr		= current_mu; 
// 			if(obj > obj_old) break; 
		}while(check==0 && iter < outer_iter);
		
		// Return the minimizer and optimal value; 
		for(int j=0; j<n; j++){
			out(j) = gsl_vector_get(x,j);
		}
		out(n) = orig_F.f(x, orig_F.params);
		gsl_vector_free(x); 
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
	
	NumericVector conv_constrOpitm(gsl_vector *x0){
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
			out(n) 	= -1e5; 
			gsl_vector *b; 
			gsl_matrix *init; 
			init	= gsl_matrix_alloc (n, 4);
			b		= gsl_vector_alloc(n); 
			
			// Set initial values; 
			gsl_vector_set(b, 0, log(upper_y)); 
			gsl_vector_set(b, 1, log(upper_c)); 
			gsl_matrix_set_col(init, 0, x0); 
			gsl_matrix_set_col(init, 1, b); 
			gsl_matrix_set(init, 0, 2, gsl_vector_get(b,0));
			gsl_matrix_set(init, 1, 2, gsl_vector_get(x0,1));
			gsl_matrix_set(init, 0, 3, gsl_vector_get(x0,0));
			gsl_matrix_set(init, 1, 3, gsl_vector_get(b,1));
			
			for(int i=0; i<4; i++){
				gsl_matrix_get_col(x0, init, i);
				NumericVector out1 = constr_optim(x0);
				// Check boundary and convert back the variables; 		
				if(exp(out1[0]) < bound_eps){
					out1[0] = 0; 
				}else if(upper_y- exp(out1[0]) < bound_eps){
					out1[0] = upper_y; 	
				}else{
					out1[0] = exp(out1[0]); 
				}
		
				double Q = as<double>(Q_fn(out1[0]));
				double lower_c = std::max(0.0, s1[0] + Q - Imax), upper_c1 = s1[0] + Q; 
				if(exp(out1[1]) < bound_eps){
					out1[1] = lower_c; 
				}else if(upper_c1 - lower_c - exp(out1[1]) < bound_eps){
					out1[1] = upper_c1; 
				}else{
					out1[1] = lower_c + exp(out1[1]); 
				}
		
				gsl_vector_set(x, 0, out1[0]); 
				gsl_vector_set(x, 1, out1[1]); 
				out1(n) = orig_F.f(x, orig_F.params);	
				if(out1(n) > out(n)){
					out = clone(out1); 
				}
			}
						
			// Check boundary and convert back the variables; 		
// 			if(exp(out[0]) < bound_eps){
// 				out[0] = 0; 
// 			}else if(upper_y- exp(out[0]) < bound_eps){
// 				out[0] = upper_y; 	
// 			}else{
// 				out[0] = exp(out[0]); 
// 			}
// 		
// 			double Q = as<double>(Q_fn(out[0]));
// 			double lower_c = std::max(0.0, s1[0] + Q - Imax), upper_c1 = s1[0] + Q; 
// 			if(exp(out[1]) < bound_eps){
// 				out[1] = lower_c; 
// 			}else if(upper_c1 - lower_c - exp(out[1]) < bound_eps){
// 				out[1] = upper_c1; 
// 			}else{
// 				out[1] = lower_c + exp(out[1]); 
// 			}
// 		
// 			gsl_vector_set(x, 0, out[0]); 
// 			gsl_vector_set(x, 1, out[1]); 
// 			out(n) = orig_F.f(x, orig_F.params);	
			gsl_vector_free(b); gsl_matrix_free(init); 	
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
	struct my_f_params *params = (struct my_f_params *)p;
	NumericVector param_vec		= (params->param_vec);	
	NumericVector state_vec		= (params->state_vec);	
	double beta					= (params->beta);	
	Function Q_fn				= (params->Q_fn);	
	Function omega_fn			= (params->omega_fn);
	double lnzero				= (params->lnzero);
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
	List param_out 	= param_assign(param_vec);
	double lambda	= param_out["lambda"], tau = param_out["tau"], mu = param_out["mu"], sigma = param_out["sigma"];
	double I 		= state_vec[0], Z = exp(state_vec[3]); 
	int Inc_index 	= state_vec[1], Rc = state_vec[2];
	double Inc 		= Inc_nodes[Inc_index - 1], y = gsl_vector_get(a, 0), c = gsl_vector_get(a, 1);
	
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

void my_valuefn_g(const gsl_vector *a, void *p, gsl_vector *df){
	struct my_f_params *params = (struct my_f_params *)p;
	NumericVector param_vec		= (params->param_vec);	
	NumericVector state_vec		= (params->state_vec);	
	double beta					= (params->beta);	
	Function Q_fn				= (params->Q_fn);	
	Function omega_fn			= (params->omega_fn);
	Function Q_j				= (params->Q_j);	
	Function omega_j			= (params->omega_j);	
	double lnzero				= (params->lnzero);
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
	List param_out 	= param_assign(param_vec);
	double lambda	= param_out["lambda"], tau = param_out["tau"], mu = param_out["mu"], sigma = param_out["sigma"];
	double I 		= state_vec[0], Z = exp(state_vec[3]); 
	int Inc_index 	= state_vec[1], Rc = state_vec[2];
	double Inc 		= Inc_nodes[Inc_index - 1], y = gsl_vector_get(a, 0), c = gsl_vector_get(a, 1);
	
	// Initiate working variables; 
	int n_Inc = Inc_nodes.length(), n_Z = lnZ_nodes.length(), 
		N1 = x1.length(), N2 = x2.length(), N3 = x3.length(), N4 = x4.length(); 
	double out = 0, e_dv_Inext = 0; 
	NumericVector s_next(4), y4(N4), v(n_Z), edv_Inc(n_Inc), lnZ_next = sqrt(2)*sigma*lnZ_nodes + mu; 
	gsl_set_error_handler_off();

	// Transition of the next state; 
	double Q = as<double>(Q_fn(y));
	s_next[0] = I + Q - c; 
	s_next[1] = 0; 
	s_next[2] = Rc;
	s_next[3] = lnZ_nodes[0];
	NumericVector Inc_weight1 = Inc_weight((n_Inc*Rc + Inc_index - 1),_);

	// Compute the expected value of dV(s')/dI'; 
	for(int i=0; i<n_Inc; i++){
		s_next[1] = i + 1; 
		int x2_index = findInterval1(s_next(1), x2);
		int x3_index = findInterval1(s_next(2), x3);
		
		// Take derivative w.r.t. I_next; 	
		for(int j=0; j<N4; j++){
			int index = j*N2*N3 + (x3_index - 1)*N2 + x2_index - 1; 
			y4(j) = spl2d[index].spldevEval(s_next(0));
		}
		
		// Create the spline interpolation along x4; 			
		struct spl myspl = spl_init(x4, y4);
		for(int j=0; j<n_Z; j++){
			v[j] = myspl.splEval(lnZ_next[j]);
		}
		myspl.splfree();
		edv_Inc(i) = sum(lnZ_weight*v)/sqrt(PI);				
	}
	e_dv_Inext = sum(edv_Inc * Inc_weight1); 

	// Compute the derivatives w.r.t. y and c; 
	double domega	= as<double>(omega_j(y)); 
	double dQ		= as<double>(Q_j(y)); 
	gsl_vector_set(df, 0, domega - 1 - tau*dQ + beta*dQ*e_dv_Inext); 
	gsl_vector_set(df, 1, lambda/(c+.01) + tau - beta*e_dv_Inext); 
}

double my_value_obj(const gsl_vector *a_t, void *p){
	struct constr_par *params 	= (struct constr_par *)p;
	void *vf_par				= (params->f_par); 
	struct my_f_params *f_par	= (struct my_f_params *)vf_par; 
	double *mu					= (params->mu); 
	double mu_n					= *mu; 	
	NumericVector state_vec		= (f_par->state_vec);
	Function Q_fn				= (f_par->Q_fn);
	double lnzero				= (f_par->lnzero);	
	NumericVector Inc_nodes		= (f_par->Inc_nodes);	
	NumericVector x1			= (f_par->x1);
	
	double I		= state_vec[0], Z = exp(state_vec[3]); 
	int Inc_index 	= state_vec[1], Rc = state_vec[2], N1 = x1.length();
	double Inc 		= Inc_nodes[Inc_index - 1]; 

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
	NumericVector g(2), pnty(2);
	g(0) = y - (Inc - Z)*1*(Inc>Z); 
	g(1) = c - (I + Q); 
	pnty(0)	= std::max(0.0, g(0)); 
	pnty(1) = std::max(0.0, g(1)); 
	double out; 
	
	out = my_valuefn(a, vf_par) - mu_n*sum(pow(pnty,2)); 
	gsl_vector_free(a);
	return(-out);		// Notice to return negative value function; 
}

void my_value_obj_g(const gsl_vector *a_t, void *p, gsl_vector *df){
	struct constr_par *params 	= (struct constr_par *)p;
	void *vf_par				= (params->f_par); 
	struct my_f_params *f_par	= (struct my_f_params *)vf_par; 
	double *mu					= (params->mu); 
	double mu_n					= *mu; 	
	NumericVector state_vec		= (f_par->state_vec);
	Function Q_fn				= (f_par->Q_fn);
	Function Q_j 				= (f_par->Q_j);
	double lnzero				= (f_par->lnzero);	
	NumericVector Inc_nodes		= (f_par->Inc_nodes);	
	NumericVector x1			= (f_par->x1);
	
	double I 		= state_vec[0], Z = exp(state_vec[3]); 
	int Inc_index 	= state_vec[1], Rc = state_vec[2], N1 = x1.length();
	double Inc 		= Inc_nodes[Inc_index - 1]; 	
	
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
	gsl_vector *a, *dp; 
	a = gsl_vector_alloc(2);
	dp= gsl_vector_alloc(2); 		
	gsl_vector_set(a, 0, y);
	gsl_vector_set(a, 1, c);
	
	// The derivative of penalty function; 
	double dQ = as<double>(Q_j(y));
	gsl_vector_set(dp,0, std::max(0.0, y - (Inc - Z)*1*(Inc>Z)) - dQ*std::max(0.0, c - (I + Q)) );
	gsl_vector_set(dp,1, std::max(0.0, c - (I + Q)) ) ;
	
	// The gradient of objective function;  
	my_valuefn_g(a, vf_par, df);
	gsl_blas_daxpy(-2*mu_n, dp, df); 
	// Multiply by dx/dtheta; 
	gsl_vector_set(a, 0, exp(gsl_vector_get(a_t, 0)));
	gsl_vector_set(a, 1, exp(gsl_vector_get(a_t, 1)));
// 	gsl_vector_mul(df, a); 
	gsl_vector_scale(df, -1.0); 
	
	// Free memory; 
	gsl_vector_free(a); gsl_vector_free(dp);
}

void my_value_obj_fdf(const gsl_vector *x, void *params, double *f, gsl_vector *df){
	*f = my_value_obj(x, params); 
	my_value_obj_g(x, params, df);
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

//-------------------------------------------------//;
// Solving Dynamic programming -- policy iteration //;
//-------------------------------------------------//;
// [[Rcpp::export]]
List policy_improveC(NumericVector V_old, List DP_list, List control_list, double multimin_tol=0, int multimin_max = 0){
	// Unload parameters; 
	NumericVector param = DP_list["param"]; 
	NumericMatrix state = DP_list["state"]; 
	double beta 		= DP_list["beta"]; 
	Function omega_fn 	= DP_list["omega"], Q_fn = DP_list["Q_fn"], 
			 omega_j 	= DP_list["omega_j"], Q_j = DP_list["Q_j"]; 
	List Inc_kernel 	= DP_list["Inc_kernel"], lnZ_kernel = DP_list["lnZ_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"], 
			lnZ_nodes 	= lnZ_kernel["nodes"], lnZ_weight = lnZ_kernel["weight"]; 
	NumericMatrix Inc_weight = Inc_kernel["weight"]; 
	int ns = state.nrow(), n_act = 2;
	size_t n_act1 		= 2; 
	double lnzero 		= control_list["lnzero"], Qk = DP_list["k"];
	if(multimin_tol ==0) {multimin_tol = control_list["multimin_tol"]; }
	if(multimin_max ==0) {multimin_max = control_list["multimin_max"];}
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
	struct my_f_params params = {param, state(0,_),beta, Q_fn, omega_fn, Q_j, omega_j, lnzero, 
								Inc_nodes, Inc_weight, lnZ_nodes, lnZ_weight, 
								 x1, x2, x3, x4, myspl};	
	double mu_n = 1; 
	struct constr_par constrp = {&params, &mu_n}; 
	
	// Initiate multimin minimizer; 
	const gsl_multimin_fdfminimizer_type *T;
  	gsl_multimin_fdfminimizer *s;	
	T = gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (T, n_act);
	
	// Claim minimizing function; 
  	gsl_multimin_function_fdf my_func;
  	my_func.n 		= n_act;
	my_func.f 		= my_value_obj;
	my_func.df 		= my_value_obj_g;
	my_func.fdf 	= my_value_obj_fdf;
	my_func.params 	= &constrp;
	
	// Claim the original function; 
	gsl_multimin_function_fdf orig_func;
  	orig_func.n 	= n_act;
	orig_func.f 	= my_valuefn;
	orig_func.df 	= my_valuefn_g;
	orig_func.params= &params;
	
	// Claim the 1d function when y=0; 
	const gsl_min_fminimizer_type *T1d;
	gsl_min_fminimizer *s1d;
	T1d = gsl_min_fminimizer_brent;
	s1d = gsl_min_fminimizer_alloc(T1d);

	gsl_function F1d; 
	F1d.function 	= my_value_1c;
	F1d.params		= &params;
		
	// Initiate constrained optimization; 	
	struct pen_constrOptim my_multmin = {s, my_func, orig_func, s1d, F1d, 1e-4, 100, .1, 10, 1000, 1e-4, 200, 1e-4};  
	
	// Initiate working variables;
	NumericMatrix policy(ns, n_act); 
	NumericVector V(ns), min_f(n_act+1); 
	gsl_vector *m; 
	m = gsl_vector_alloc(n_act); 
	
	// Loop over each state to solve the maximization; 
	for(int i=0; i<ns; i++){
		params.state_vec = state(i,_); 
		
		// Starting values; 
		double Inc = Inc_nodes[state(i,1) - 1]; 
		if(Inc - exp(state(i,3)) <= 0){
			gsl_vector_set(m, 0, lnzero); 
			gsl_vector_set(m, 1, std::max(lnzero, log(state(i,0))) ); 
		}else{
			gsl_vector_set(m, 0, log(.5*(Inc - exp(state(i,3)))) );
			gsl_vector_set(m, 1, log( (state(i,0) + .5*Qk*(Inc - exp(state(i,3)) )/2 )*.8) );
		}	 
				
		min_f = my_multmin.conv_constrOpitm(m);
// 		min_f = my_multmin.constr_optim(m);
		for(int j=0; j<n_act; j++){
			policy(i,j) = min_f(j); 
		}
		V(i)		= min_f(n_act);
	}
	
	// Free memory; 
	gsl_vector_free(m);
	gsl_multimin_fdfminimizer_free (s);
	gsl_min_fminimizer_free(s1d);
	for(int i=0; i<nspl; i++){
		myspl[i].splfree();
	}
	return Rcpp::List::create( _["policy"]=policy, _["V"] = V); 	
}

