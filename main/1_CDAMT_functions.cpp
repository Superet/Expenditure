#include <RcppGSL.h>
#include <Rcpp.h>
#include <time.h> 
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]

//------------------//;
// Define structure //;
//------------------//;
struct paramStr{double lambda; double lambda_o; double tau1; double tau2; double tau3; double tau4; }; 

struct paramStr param_assign(NumericVector param){
	struct paramStr out = {param(0), param(1), param(2), param(3), param(4), param(5)};
	return(out); 
}

struct spl{
	gsl_spline *spline; gsl_interp_accel *acc; 
	NumericVector x; NumericVector y; double slope_first; double slope_end; 
			
	inline double splEval(double x0){
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

inline struct spl spl_init(NumericVector x, NumericVector y){
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

struct DP_fun{
	NumericVector omega_y; NumericMatrix omega_val; struct spl *spl2d;
	NumericVector Q_coef; int nq; 
	
	inline double omega_fn(double y, int d){
		double out = 0; 
		if(d > 0){
			out = spl2d[d-1].splEval(y);
		} 
		return(out); 
	}
	
	inline double Q_fn(double y){
		double out		= log(y+1) * Q_coef(0); 
		return(out); 
	}
};

struct my_f_params {int *d; struct paramStr param_out; NumericVector state_vec; 
					double beta; struct DP_fun DPf; double lnzero; 
					NumericVector Inc_nodes; NumericMatrix Inc_weight; 
					NumericVector x1; NumericVector x2; NumericVector x3; struct spl *spl2d;};
					
struct myDP{NumericVector policy; NumericMatrix cond_policy;  NumericVector W; double V; double Q; };					
					
struct mymin{
// s is the NM minimizer, the function F is negative of value function. The parameters need to be converted when computing value function; 
// s1d is brent one-dimensional minimizer when y = 0, the function F1d is negative value function; 
	gsl_multimin_fminimizer *s; 	gsl_multimin_function F; gsl_min_fminimizer *s1d; gsl_function F1d; 
	gsl_multimin_function orig_F; 	double bound_eps; 
	double NM_sizetol; 				double NM_startsize;	int NM_iter; 
	double brent_tol; 				int brent_iter; 		int D; 
		
	inline NumericVector inner_multimin(gsl_vector *x){
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
	// Output: y*, c*, -v(y*,c*) conditional on d, Q; 
		// When y = 0; 
		struct my_f_params *params 	= (struct my_f_params *)orig_F.params;
		NumericVector s1			= params->state_vec; 
		NumericVector x1			= params->x1;			
		double lnzero				= params->lnzero;
		double 	lower_c 			= std::max(0.0, s1(0)+0.0-x1[x1.length()-1]), 
				upper_c 			= s1(0);
		double m = 0, v= 0, low = lower_c, up = upper_c;
		int sat, status;
		int iter	= 0;
		
		// Try initial values close to the boundaries;
		m		= std::max(upper_c - brent_tol,0.0);
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

		NumericVector out(F.n+2); 
		out(0) = lnzero; 
		out(1) = m; 
		out(2) = v; 
		out(3) = 0.0; 
		return(out); 	
	}	
			
	NumericVector multmin_solver(const gsl_vector *x0) {
	// Return y, c, V(d, y, c) and Q; 
		struct my_f_params *params 	= (struct my_f_params *)orig_F.params;
		NumericVector s1			= params->state_vec;
		struct DP_fun DPf			= params->DPf; 
		NumericVector Inc_nodes 	= params->Inc_nodes;	
		NumericVector x1			= params->x1;
		int n 						= orig_F.n; 
		double Imax 				= x1(x1.length()-1); 
		
		NumericVector out(n+2); 
		gsl_vector *x, *x0l, *b;
		x 				= gsl_vector_alloc(n);
		x0l 			= gsl_vector_alloc(n); 
		b				= gsl_vector_alloc(n); 	
		gsl_matrix *init; 
		init			= gsl_matrix_alloc (n, 3);
			
		// Set initial values; 
		out(n) = 1e6;
		double eps = NM_startsize ; 
		double upper_Q	= DPf.Q_fn(gsl_vector_get(x0,0));
		double upper_c	= s1(0) + std::max(upper_Q, 0.0); 
		gsl_vector_set(b, 0, log(gsl_vector_get(x0,0))); 
		gsl_vector_set(b, 1, log(upper_c)); 
		gsl_matrix_set_col(init, 0, x0); 
		gsl_matrix_set_col(init, 1, b); 
		gsl_vector_set_all(x0l, params->lnzero - NM_sizetol); 
		
		// Solve for the minimizer for each initial values; 
		for(int i=0; i<3; i++){
			if(gsl_vector_get(x0l,0)==gsl_matrix_get(init,0,i) && gsl_vector_get(x0l,1)==gsl_matrix_get(init,1,i)){
				continue; 
			}
			gsl_matrix_get_col(x0l, init, i);
			NumericVector out1 = inner_multimin(x0l);
			if(out1(n) < out(n)){
				out(0)	= out1(0); 
				out(1)	= out1(1); 
				out(2)	= out1(2); 
			}
			if(i<2){
				gsl_matrix_set(init, 0, i, out(0)); 
				gsl_matrix_set(init, 1, i, out(1)); 
			}
			if(i==1){
				gsl_matrix_set(init, 0, i+1, .5*(gsl_matrix_get(init, 0,0) + gsl_matrix_get(init, 0,1))); 
				gsl_matrix_set(init, 1, i+1, .5*(gsl_matrix_get(init, 1,0) + gsl_matrix_get(init, 1,1))); 
			}
		}
			
		// Check boundary and convert back the variables; 		
		if(exp(out(0)) < bound_eps){
			out(0) = 0; 
		}else{
			out(0) = exp(out(0)); 
		}

		double Q = DPf.Q_fn(out(0));
		double lower_c = std::max(0.0, s1[0] + Q - Imax), upper_c1 = s1(0) + std::max(Q,0.0); 
		if(exp(out(1)) < bound_eps){
			out(1) = lower_c; 
		}else if(upper_c1 - lower_c - exp(out(1)) < bound_eps){
			out(1) = upper_c1; 
		}else{
			out(1) = lower_c + exp(out(1)); 
		}

		gsl_vector_set(x, 0, out(0)); 
		gsl_vector_set(x, 1, out(1)); 
		out(n) 		= orig_F.f(x, orig_F.params);
		out(n+1)	= Q; 	
		
		gsl_vector_free(b);
		gsl_matrix_free(init);
		gsl_vector_free(x); gsl_vector_free(x0l); 

		return(out);
	}
	
	int action_solver(gsl_vector *x0, struct myDP *DP1_ptr){	
		int n = F.n; 
		struct my_f_params *params 	= (struct my_f_params *)orig_F.params;
		int *d_ptr 					= params->d; 
		NumericVector s1			= params->state_vec;
		struct DP_fun DPf			= params->DPf; 
		NumericVector Inc_nodes 	= params->Inc_nodes;	
		double lnzero				= params->lnzero;
		NumericVector x1			= params->x1;
		double Imax 				= x1(x1.length()-1); 
		NumericMatrix cond_sol(n, D+1); 
		cond_sol.fill(-999); 
		NumericVector multimin_out(n+2); 
		
		// When d=0, then y=0. We only solve for c; 
		// y, c, and value function conditional on d=0; 		
		*d_ptr	= 0; 
		(DP1_ptr->W).fill(0); 
		(DP1_ptr->policy).fill(0); 
		multimin_out	= opt1d();
		cond_sol(0,0) 	= 0.0; 
		cond_sol(1,0)	= multimin_out(1); 
		(DP1_ptr->W)[0]	= (-multimin_out(2)); 
		double v_max 	= (-multimin_out(2)), Q = 0; 
		int max_idx		= 0; 

		for(int i=1; i<=D; i++){
			*d_ptr	= i; 
			multimin_out	= multmin_solver(x0); 
			cond_sol(0,i)	= multimin_out(0); 
			cond_sol(1,i)	= multimin_out(1); 
			(DP1_ptr->W)[i]	= multimin_out(2); 
			if(multimin_out(n) > v_max && multimin_out(0) > 0 ){
				max_idx 	= i; 
				v_max		= multimin_out(n); 
				Q			= multimin_out(n+1); 
			}
		}			
		DP1_ptr->policy[0]	= max_idx; 
		DP1_ptr->policy[1]	= cond_sol(0, max_idx); 
		DP1_ptr->policy[2]	= cond_sol(1, max_idx); 
		DP1_ptr->Q			= Q; 
		
		DP1_ptr->cond_policy	= cond_sol; 
		DP1_ptr->V				= v_max; 
		return 0; 
	}
};

struct policy_fn{
	int n_act; NumericVector Inc_nodes; NumericMatrix Inc_weight; struct DP_fun DPf; 
	NumericVector x1; NumericVector x2; NumericVector x3; int D; 
	struct spl *spl2d_y; struct spl *spl2d_c; gsl_multimin_function orig_F; 
	
	NumericVector policySimall(NumericVector s1, NumericVector logit_draw){
	// This function simulates all the action variables; 
		struct my_f_params *params 	= (struct my_f_params *)orig_F.params;
		int N1 		= x1.length(), N2 = x2.length(), N3 = x3.length();
		int *d_ptr 	= (params->d); 	
		NumericVector out(n_act + 2);
		double Q; 
		int x2_index = s1(1);
		int x3_index = s1(2) + 1;
		gsl_vector *x; 
		x	= gsl_vector_alloc(n_act); 
		double v_max = -1e6, v; 
	
		// Interpolate I_next conditional on the other dimensions; 			
		int index 		= (x3_index - 1)*N2 + x2_index - 1; 		
		double ystar	= spl2d_y[index].splEval(s1(0));
		double cstar	= spl2d_c[index].splEval(s1(0));
		
		// Adjust the interior solution to satisfy boundary conditions; 
		double lower_c, tmpc, tmpQ; 
		out(1) = std::max(0.0, ystar);
		gsl_vector_set(x, 0, out(1)); 
		
		// Grid search for d; 
		if(out(1) ==0) {
			out(0) 	= 0; 
			Q 		= DPf.Q_fn(out[1]); 
			lower_c = std::max(0.0, s1(0) + Q - x1[N1-1] ); 
			out(2) 	= std::max(lower_c, std::min(cstar, s1(0) + Q));
		}else{
			for(int j=1; j<=D; j++){
				*d_ptr 	= j; 
				tmpQ 	= DPf.Q_fn(out(1));
				lower_c = std::max(0.0, s1(0) + tmpQ - x1[N1-1] );
				tmpc 	= std::max(lower_c, std::min(cstar, s1(0) + tmpQ));
				gsl_vector_set(x, 1, tmpc); 
				v = orig_F.f(x, orig_F.params) + logit_draw(j) ;
				if(v > v_max){
					out(0) 	= j; 
					out(2) 	= tmpc; 
					Q		= tmpQ; 
					v_max  	= v; 
				}
			}
		}
		gsl_vector_free(x); 
		out(3) = s1(0) + Q - out(2);
		return(out); 		
	}
	
	NumericVector policySim2(NumericVector s1, int d){
	// Output: (d, y, c, Inext); 
		struct my_f_params *params 	= (struct my_f_params *)orig_F.params;
		NumericVector Inc_nodes		= params->Inc_nodes; 
		int N1 = x1.length(), N2 = x2.length(), N3 = x3.length();
		NumericVector out(n_act + 2);
		double Q; 
		int x2_index = s1(1);
		int x3_index = s1(2) + 1;

		// Interpolate I_next; 			
		int index 		= (x3_index - 1)*N2 + x2_index - 1; 			
		double ystar	= spl2d_y[index].splEval(s1(0)); 
		double cstar	= spl2d_c[index].splEval(s1(0));

		// Adjust the interior solution to satisfy boundary conditions; 
		out(0) 			= d; 
		out(1) 			= std::max(0.0, ystar);
		Q 				= DPf.Q_fn(out(1));
		double lower_c 	= std::max(0.0, s1(0) + Q - x1[N1-1] );
		out(2) 			= std::max(lower_c, std::min(cstar, s1(0) + Q));
		out(3) 			= s1(0) + Q - out(2); 
		return(out);
	} 	
	
	NumericVector policySim1(NumericVector s1, int d, double y){
	// Output: (d, y, c, Inext); 
		int N1 = x1.length(), N2 = x2.length(), N3 = x3.length();
		NumericVector out(n_act + 2);
		double Q; 
		int x2_index = s1(1);
		int x3_index = s1(2) + 1;

		// Interpolate I_next conditional on the 4th dimension; 			
		int index		= (x3_index - 1)*N2 + x2_index - 1; 
		double cstar	= spl2d_c[index].splEval(s1(0));
	
		// Adjust the interior solution to satisfy boundary conditions; 
		out(0) 	= d; 
		out(1) 	= y;
		Q 		= DPf.Q_fn(out(1)); 
		double lower_c = std::max(0.0, s1(0) + Q - x1[N1-1] );
		out(2) = std::max(lower_c, std::min(cstar, s1(0) + Q));
		out(3) = s1(0) + Q - out(2); 
		return(out);
	}
	
	inline NumericMatrix policySim12(NumericVector s1, int d, double y, double Q_obs){
	// Output: matrix (2,4), col(d, y, c, Inext), row(sim1, sim2); 
		int N1 = x1.length(), N2 = x2.length(), N3 = x3.length();
		NumericMatrix out(2, n_act+2); 
		double Q; 
		int x2_index = s1(1);
		int x3_index = s1(2) + 1;

		// Interpolate I_next conditional on the 4th dimension; 			
		int index 		= (x3_index - 1)*N2 + x2_index - 1; 	
		double ystar	= spl2d_y[index].splEval(s1(0));
		double cstar	= spl2d_c[index].splEval(s1(0));
	
		//sim1: only simulate c conditional on y and d; 
		// Adjust the interior solution to satisfy boundary conditions; 
		double lower_c1	= std::max(0.0, s1(0) + Q_obs - x1[N1-1] );
		out(0,0)		= d; 
		out(0,1)		= y; 
		out(0,2) 		= std::max(lower_c1, std::min(cstar, s1(0) + Q_obs));
		out(0,3) 		= s1(0) + Q_obs - out(0,2); 
		
		//sim2: simulate y and c conditional on d; 
		out(1,0) 		= d; 		
		out(1,1) 		= std::max(0.0, ystar);
		Q 				= DPf.Q_fn(out(1,1));
		double lower_c2	= std::max(0.0, s1(0) + Q - x1[N1-1] );
		out(1,2)		= std::max(lower_c2, std::min(cstar, s1(0) + Q));
		out(1,3)		= s1(0) + Q - out(1,2); 
		return(out);
	}  
}; 		
					
//-----------------------//;
// Define value function //;
//-----------------------//;
double my_valuefn(const gsl_vector *a, void *p){
	struct my_f_params *params 	= (struct my_f_params *)p;
	int *d_ptr					= params->d;	
	struct paramStr param_out	= params->param_out;	
	NumericVector state_vec		= params->state_vec; 
	struct DP_fun DPf			= params->DPf; 
	NumericVector Inc_nodes 	= params->Inc_nodes;	
	NumericMatrix Inc_weight	= params->Inc_weight;	
	NumericVector x1			= params->x1;
	NumericVector x2			= params->x2;
	NumericVector x3			= params->x3;
	struct spl *spl2d			= params->spl2d;
	
	// Assign parameter; 
	double I 					= state_vec[0]; 
	int Inc_index 				= state_vec[1], Rc = state_vec[2];
	double Inc 					= Inc_nodes[Inc_index - 1], y = gsl_vector_get(a, 0), c = gsl_vector_get(a, 1);
	int n_Inc 					= Inc_nodes.length(), 
		N2 						= x2.length(), N3 = x3.length(); 
	
	// Transition of the next state; 
	int d_n		= 0; 
	if(y > 0)	{ d_n = *d_ptr; }
	NumericVector s_next(3); 
// 	double Q 	= DPf.Q_fn(y);
	double Q	= log(y + 1) * DPf.Q_coef(0); 
	s_next[0] 	= I + Q - c; 
	s_next[1] 	= Inc_index; 
	s_next[2] 	= Rc;
	
	// Markov income; 
	NumericVector Inc_weight0		= Inc_weight((n_Inc*Rc + Inc_index - 1),_);
	int n_Incnext					= sum(Inc_weight0>0), cnt = 0; 
	IntegerVector Inc_next_index(n_Incnext); 
	for(int i=0; i<n_Inc; i++){
		if(Inc_weight0(i) > 0){
			Inc_next_index(cnt) = i + 1; 
			cnt 				= cnt + 1; 
		}
	}
	NumericVector Inc_weight1		= Inc_weight0[Inc_next_index - 1]; 
	
	// Initiate working variables; 	
	double out = 0, ev = 0; 
	NumericVector ev_Inc(n_Incnext);
	int index =0; 
	
	// Compute the expected option value; 
	for(int i=0; i<n_Incnext; i++){
		s_next(1) 		= Inc_next_index(i); 
		index 			= s_next(2)*N2 + s_next(1) - 1; 
		ev_Inc(i) 		= spl2d[index].splEval(s_next(0));			
	}
	ev = sum(ev_Inc * Inc_weight1); 

	// Compute the value function;
// 	double omega = DPf.omega_fn(y, d_n); 
	double omega = 0; 
	if(d_n > 0){
		omega = (DPf.spl2d)[d_n-1].splEval(y);
	} 
	out = (param_out.lambda_o)*omega + (param_out.lambda) * log(c+1) + 0.1*(Inc-y) + (param_out.tau1) * s_next[0] 
		 + (param_out.tau2) * pow(s_next[0], 2) + (param_out.tau3) * s_next[0] * d_n + (param_out.tau4) * d_n + (params->beta)*ev; 
	if(!std::isfinite(out) || std::abs(out) > 999999){
		out = 999999;
	}
	return(out);
}	

double my_valuefn_t(const gsl_vector *a_t, void *p){
	struct my_f_params *params 	= (struct my_f_params *)p;
	int *d_ptr					= params->d; 
	struct DP_fun	DPf			= params->DPf; 
	double lnzero				= params->lnzero;	
	NumericVector Inc_nodes		= params->Inc_nodes;	
	NumericVector x1			= params->x1;
	double I					= (params->state_vec)(0); 
	int Inc_index 				= (params->state_vec)(1), N1 = x1.length();
	double Inc 					= Inc_nodes[Inc_index - 1]; 

	// Variable transformation for the control variables and construct constraints; 
	double y		= 0, c= 0;
	if( gsl_vector_get(a_t, 0) > lnzero){
		y 			= exp(gsl_vector_get(a_t, 0));
	}
	double d_n		= 0; 
	if(y > 0){ d_n = *d_ptr; }
	double Q 		= DPf.Q_fn(y);
	double lowerc	= std::max(0.0, I + Q - x1(N1-1));
	if( gsl_vector_get(a_t, 1) > lnzero){
		c 			= lowerc + exp(gsl_vector_get(a_t, 1)); 
	}else{
		c 			= lowerc;
	}
	gsl_vector *a; 
	a = gsl_vector_alloc(2);
	gsl_vector_set(a, 0, y);
	gsl_vector_set(a, 1, c);
	
	// Set the constraints: g(x) <= 0; 
	NumericVector g(2), sign_g(2), my_penalty(2); 
	sign_g.fill(0);
	g(0) 	= c - (I + Q); 
	g(1)	= y - Inc;
	double v, out, tol_g = 1e-6; 
	if(g(0) > tol_g) {sign_g(0) = 1; }
	if(g(1) > tol_g) {sign_g(1) = 1; }
	my_penalty(0) = 1000; 
	my_penalty(1) = 1000; 
		
	v = my_valuefn(a, p);
	if(sum(sign_g) == 0){
		out = v; 
	}else{
		out = v - sum(my_penalty * g*sign_g);
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

//-------------------------------------------------//;
// Solving Dynamic programming -- policy iteration //;
//-------------------------------------------------//;
List policy_improveC(NumericVector V_old, List DP_list, List control_list, NumericMatrix policy_inits, struct DP_fun DPf, bool inits_gen=false, double NM_sizetol=0, int NM_iter = 0, double NM_startsize = 0){
//NOTE: check the Q_fn when generating initial values; 
	// Unload parameters; 
	NumericVector param 	= DP_list["param"]; 
	NumericMatrix state 	= DP_list["state"]; 
	double beta 			= DP_list["beta"]; 
	List Inc_kernel 		= DP_list["Inc_kernel"];
	NumericVector Inc_nodes = Inc_kernel["nodes"];
	NumericMatrix Inc_weight= Inc_kernel["weight"]; 
	int ns 					= state.nrow(), n_act = 2, D = DP_list["D"];
	double lnzero 			= control_list["lnzero"];
	double brent_tol 		= control_list["brent_tol"], bound_eps = control_list["bound_eps"]; 
	int brent_iter 			= control_list["brent_iter"]; 
	if(NM_sizetol ==0) 		{NM_sizetol 	= control_list["NM_sizetol"]; }
	if(NM_iter ==0) 		{NM_iter 		= control_list["NM_iter"];}
	if(NM_startsize ==0)	{NM_startsize	= control_list["NM_startsize"]; }
	gsl_set_error_handler_off();
	
	// Sort out the grids for each dimension in the states; 
	NumericVector x1 = unique(state(_,0)), x2 = unique(state(_,1)), x3 = unique(state(_,2));
	int N1 = x1.length(), N2 = x2.length(), N3 = x3.length();
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
	std::sort(x3.begin(), x3.end()); 
	
	// Initiate the c-spline; 
	int nspl = N2*N3;
	struct spl *myspl	= new struct spl[nspl];
	for(int i=0; i<nspl; i++){
		IntegerVector index	= i*N1 + seq_len(N1) - 1;
		NumericVector tmpy 	= V_old[index];
		myspl[i] = spl_init(x1, tmpy);
	}

	// Assign parameters for the value function; 
	struct paramStr param_out = param_assign(param); 
	int d = 1; 
	int *d_ptr = &d; 
	struct my_f_params params = {d_ptr, param_out, state(0,_),beta, DPf, lnzero, 
								Inc_nodes, Inc_weight, x1, x2, x3, myspl};	
	
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
							   NM_sizetol, NM_startsize, NM_iter, brent_tol, brent_iter, D};  
	// Initiate working variables;
	NumericMatrix policy(ns, n_act + 1), cond_y(ns, D+1), cond_c(ns, D+1), cond_sol(n_act, D+1); 
	NumericVector V(ns), Q_vec(ns), W(D+1); 
	struct myDP DP1_sol = {policy(0,_), cond_sol, W, 0.0}; 

	// Loop over each state to solve the maximization; 
	double eps = NM_startsize; 
	double Inc, upper_y; 
	for(int i=0; i<ns; i++){
		Inc 			= Inc_nodes(state(i,1) - 1); 
		params.state_vec= state(i,_); 

		// Starting values; 
		if(inits_gen){			
			gsl_vector_set(m, 0, log(.1*Inc) );
			gsl_vector_set(m, 1, log(.5*(state(i,0) + .001))) ; 
		}else{
			double Q		= log(policy_inits(i,1) + 1) * DPf.Q_coef(0); 
			double lower_c 	= std::max(0.0, state(i,0) + Q - x1[N1-1]); 
			if(policy_inits(i,1)==0){
				gsl_vector_set(m, 0, lnzero); 
			}else{
				gsl_vector_set(m, 0, log(policy_inits(i,1))); 
			}
			
			if(policy_inits(i,2) == lower_c){
				gsl_vector_set(m, 1, lnzero );
			}else if(policy_inits(i,2) == state(i,0) + Q){
				gsl_vector_set(m, 1, log(policy_inits(i,2) - lower_c) );
			}else{
				gsl_vector_set(m, 1, log(policy_inits(i,2) - lower_c) );
			}
						 
		}

		my_multmin.action_solver(m, &DP1_sol);
		policy(i,_) 	= DP1_sol.policy; 
		V(i)			= DP1_sol.V;
		cond_y(i,_)		= DP1_sol.cond_policy(0,_); 
		cond_c(i,_)		= DP1_sol.cond_policy(1,_); 
		Q_vec(i)		= DP1_sol.Q; 
	}
	NumericVector c_ratio	= policy(_,2)/(state(_,0) + Q_vec); 
	c_ratio[is_na(c_ratio)] = 1; 
	
	gsl_vector_free(m);
	gsl_multimin_fminimizer_free (s);
	gsl_min_fminimizer_free(s1d);
	for(int i=0; i<nspl; i++){
		myspl[i].splfree();
	}
	delete[] myspl; 

	return Rcpp::List::create( _["policy"]=policy, _["V"] = V, _["c_ratio"]	= c_ratio, 
							   _["cond_y"]=cond_y, _["cond_c"]=cond_c); 
}

NumericVector policy_evalC(NumericMatrix policy, NumericMatrix policy_y, NumericMatrix policy_c, NumericVector V_old, NumericVector V0, List DP_list, List control_list, struct DP_fun DPf, int M=0){
// NOTE: policy_y and policy_c are policy function conditional on d; 	
	
	// Unload parameters; 
	NumericVector param 	= DP_list["param"]; 
	NumericMatrix state 	= DP_list["state"]; 
	double beta 			= DP_list["beta"]; 
	List Inc_kernel 		= DP_list["Inc_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"]; 
	NumericMatrix Inc_weight= Inc_kernel["weight"]; 
	int ns 					= state.nrow(), n_act = 2, D = DP_list["D"];
	struct paramStr param_out= param_assign(param);
	
	// Unload the iteration parameter;
	double tol = control_list["inner_tol"];
	if(M==0){
		M = control_list["inner_max"]; 
	}
	double lnzero = control_list["lnzero"];
	gsl_set_error_handler_off();
	
	// Sort out the grids for each dimension in the states; 
	NumericVector x1 = unique(state(_,0)), x2 = unique(state(_,1)), x3 = unique(state(_,2));
	int N1 = x1.length(), N2 = x2.length(), N3 = x3.length();
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
	std::sort(x3.begin(), x3.end()); 
	int nspl = N2*N3;
	
	// Initiate the iteration;
	NumericVector V_new(ns), V_oldl(ns); 
	V_new 		= clone(V0);
	V_oldl		= clone(V_old);
	double dif 	= max(abs(V_new - V_oldl)) ; 
	int k 		= 0; 
	gsl_vector *a; 
	a 			= gsl_vector_alloc(n_act);
	int d 		= 1;
	int *d_ptr 	= &d; 
	
	if(dif > tol){
		V_oldl = clone(V0); 
		do{
			k++; 
			V_new.fill(0); 
	
			// Initiate the c-spline; 
			struct spl *myspl	= new struct spl[nspl];
			for(int i=0; i<nspl; i++){
				IntegerVector index	= i*N1 + seq_len(N1) - 1;
				NumericVector tmpy 	= V_oldl[index];
				myspl[i] 			= spl_init(x1, tmpy);
			}	
			// Assign parameters for the value function; 
			struct my_f_params params = {d_ptr, param_out, state(0,_),beta, DPf, lnzero,
							Inc_nodes, Inc_weight, x1, x2, x3, myspl};	
			
			// Compute the value function at each states; 
			for(int i=0; i<ns; i++){
				params.state_vec = state(i,_); 
				*d_ptr = policy(i,0); 
		 		gsl_vector_set(a, 0, policy(i,1) ); 
				gsl_vector_set(a, 1, policy(i,2) );
				V_new(i) = my_valuefn(a, &params); 
// 				if(policy(i,0) > 0){
// 					V_new(i) = V_new(i) + 2.718282; 
// 				}
			}
			
			// Check distance; 
			dif = max(abs(V_new - V_oldl));
			V_oldl = clone(V_new); 
			
			// Free spline; 
			for(int i=0; i<nspl; i++){
				myspl[i].splfree();
			}
			delete[] myspl; 
// 			Rcout<<"Policy evaluation iter "<<k<<", dif="<<dif<<"\n"; 			
		}
		while(k <= M && dif > tol);	
	}
	gsl_vector_free(a);
		
	return(V_new); 
}

// [[Rcpp::export]]
NumericVector guessVC(List DP_list, List control_list){
// Unload parameters; 
	NumericVector param 	= DP_list["param"]; 
	NumericMatrix state 	= DP_list["state"]; 
	double beta 			= DP_list["beta"]; 
	List Inc_kernel 		= DP_list["Inc_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"]; 
	NumericMatrix Inc_weight= Inc_kernel["weight"]; 
	int ns 					= state.nrow(), n_act = 2, D = DP_list["D"];
	struct paramStr param_out= param_assign(param);
	
	List omega_fnv			= DP_list["omega"];	
	NumericVector omega_y	= omega_fnv["omega_y"], 
			Q_coef			= DP_list["Q_coef"]; 			
	NumericMatrix omega_val	= omega_fnv["omega"]; 		
	
	struct spl *omega_spl = new struct spl[D]; 
	for(int i=0; i<D; i++){
		omega_spl[i] 		= spl_init(omega_y, omega_val(i,_));
	}
	struct DP_fun DPf 		= {omega_y, omega_val, omega_spl, Q_coef, Q_coef.length()}; 
	
	NumericVector y(ns), c(ns), V(ns); 
	double Q, omega, Inext, Inc; 
	int d; 
	for(int i=0; i<ns; i++){
		d		= 3; 
		y(i)	= 100 - 10 * state(i,0); 
		Q 		= DPf.Q_fn(y(i)); 
		c(i)	= (.8 - state(i,2)*.3 )* (state(i,0) + Q); 
		Inext	= state(i,0) + Q - c(i); 
		omega 	= (DPf.spl2d)[d-1].splEval(y(i));
		Inc		= Inc_nodes(state(i,1)-1); 
		V(i) 	= (param_out.lambda_o)*omega + (param_out.lambda) * log(c(i)+1) + 0.1*(Inc-y(i)) + (param_out.tau1) * Inext 
		 		+ (param_out.tau2) * pow(Inext, 2) + (param_out.tau3) * Inext * d + (param_out.tau4) * d; 
		V(i) 	= V(i)/(1-beta); 		
	}
	for(int i=0; i<D; i++){
		omega_spl[i].splfree(); 
	}
	delete[] omega_spl; 
	return(V); 
}

// [[Rcpp::export]]
List policy_iterC(List DP_list, List control_list, int print_level=0, int init_V=1){
	double tol 				= control_list["tol"];
	int max_iter			= control_list["max_iter"], 
		D					= DP_list["D"], 
		NM_stop 			= control_list["NM_stop"];
	NumericVector V0		= DP_list["value_fn"], state = DP_list["state"];
	List omega_fnv			= DP_list["omega"];	
	NumericVector omega_y	= omega_fnv["omega_y"], 
			Q_coef			= DP_list["Q_coef"]; 			
	NumericMatrix omega_val	= omega_fnv["omega"]; 			
	int ns 					= V0.length(), n_act = 2;
	NumericVector V_new(ns); 
	NumericMatrix policy_new(ns, n_act+1);
	List DP_listout 		= clone(DP_list);

	struct spl *omega_spl	= new struct spl[D]; 
	for(int i=0; i<D; i++){
		omega_spl[i] 		= spl_init(omega_y, omega_val(i,_));
	}
	struct DP_fun DPf 		= {omega_y, omega_val, omega_spl, Q_coef, Q_coef.length()}; 
	
	// Set iteration parameters; 	
	double dif 	= tol + 1, NM_sizetol = 1e-3, NM_startsize = 0.5; 
	int iter 	= 0, Ml = 0, NM_iter = 150;
	int M 		= control_list["inner_max"], check_stop = 0;  	
	
	// Guess the initial value function; 
	if(init_V){
		V0 		= guessVC(DP_list, control_list); 
	}
	
	// Compute the policy and value function at the initial guess of value function; 
	List tmp_sol = policy_improveC(V0, DP_list, control_list, policy_new, DPf, true, NM_sizetol, NM_iter, NM_startsize);
	NumericMatrix policy_old 	= tmp_sol["policy"]; 
	NumericVector V_old			= tmp_sol["V"]; 
	NumericMatrix cond_y		= tmp_sol["cond_y"]; 
	NumericMatrix cond_c		= tmp_sol["cond_c"]; 
	NumericVector c_ratio		= tmp_sol["c_ratio"]; 
	
	bool inits_gen = false; 
	// Iterate policy function; 
	do{
		iter++; 
		double dif_old 	= dif; 
		if(iter <= 3) {
			Ml 			= 10;  
			NM_iter 	= 150; 
			NM_sizetol 	= 1e-2; 
			inits_gen	= true; 
		}else if(iter <= 5){
			Ml 			= .5*M;  
			NM_iter 	= 200; 
			NM_sizetol 	= 1e-3; 
			inits_gen	= false; 
			NM_startsize= .75; 
		}else if(iter <=8){
			NM_iter 	= 300; 
			NM_sizetol 	= 1e-5;
			NM_startsize = 1; 
		}else{ 
			Ml 			= M ; 
			NM_iter 	= 0; 
			NM_sizetol 	= 0; 
			NM_startsize= 0; 
		}
		
		// Policy improvement step; 
		tmp_sol					= policy_improveC(V_old, DP_list, control_list, policy_old, DPf, inits_gen, NM_sizetol, NM_iter, NM_startsize);
		NumericMatrix policy1	= tmp_sol["policy"]; 
		NumericMatrix cond_yl 	= tmp_sol["cond_y"], 
					  cond_cl 	= tmp_sol["cond_c"]; 
		policy_new				= clone(policy1); 
		V0 						= tmp_sol["V"]; 
		
		// Policy evaluation step;
		V_new 					= policy_evalC(policy_new, cond_yl, cond_cl, V_old, V0, DP_list, control_list, DPf, Ml); 
		dif						= max(abs(policy_new - policy_old));
		V_old					= clone(V_new); 
		policy_old				= clone(policy_new); 
		check_stop 				= 1*(iter > NM_stop && dif >= dif_old); 
		
		// Print the improvement; 
		if(print_level > 0){
			Rcout<<"Iteration "<<iter<<", norm_diff = "<<dif<<"\n";
		}
		if(iter > max_iter || dif <= tol || check_stop){
			cond_y = cond_yl; 
			cond_c = cond_cl; 
			c_ratio = tmp_sol["c_ratio"]; 
		}
	}
	while( iter<=max_iter && dif > tol && check_stop == 0);
	
	// Return the policy and value function;
	NumericVector value_fn 	= policy_evalC(policy_new, cond_y, cond_c, V_new, V0, DP_list, control_list, DPf, 500);
	DP_listout["value_fn"]	= value_fn; 
	DP_listout["policy"]	= policy_new; 
	DP_listout["c_ratio"]	= c_ratio; 
	DP_listout["Iteration"] = iter;
	double V_max			= max(value_fn); 
	
	// Check if consumption is positive; 
	NumericVector c_pos(policy_new.nrow()); 
	c_pos.fill(1); 
	for(int i=0; i<ns; i++){
		if(policy_new(i,2) <= 0 && state(i,0) > 0 && policy_new(i,1) > 0){
			c_pos(i) = 0; 
		}
	}
	
	if(dif > tol){
		DP_listout["status"] = 2;
		if(dif > .1){
			DP_listout["status"] = 3;
			Rcout<<"Policy function iteration did not converge.\n";
		}
// 		if(print_level == 1){
			
// 		}
	}else{
		DP_listout["status"] = 0;	
	}
	if(sum(c_ratio > .99) == ns){
		DP_listout["status"] = 4;
		Rcout<<"Policy function iteration generates bounded consumption.\n";
	} 
	if(V_max > 999999){
		DP_listout["status"] = 3;
		Rcout<<"Value function is out of bound.\n";
	}
	if(mean(c_pos) < 1){
		DP_listout["status"] = 4;
		Rcout<<"Policy function iteration generates zero consumption.\n";
	} 
	
	for(int i=0; i<D; i++){
		omega_spl[i].splfree(); 
	}
	delete[] omega_spl;
	
	return(DP_listout);
}

//------------------------------//;
// Simulate behavioral sequence //;
//------------------------------//;
// [[Rcpp::export]]
List sim_hhseqC(NumericVector hh_index, IntegerVector TT_vec, NumericMatrix init_state, NumericMatrix DataState, NumericMatrix choice_seq, List DP_list, List control_list, NumericVector Q, int simidx=12, bool logit_random = false, NumericMatrix logit_draw_mat = NumericMatrix(2,2)){
	// Assign parameters; 
	NumericVector param 	= DP_list["param"], 
				  V 		= DP_list["value_fn"]; 
	double	lnzero 			= control_list["lnzero"], 
			beta 			= DP_list["beta"]; 
	NumericMatrix state 	= DP_list["state"]; 
	List omega_fnv			= DP_list["omega"];
	NumericVector omega_y	= omega_fnv["omega_y"], 
			Q_coef			= DP_list["Q_coef"]; 
	NumericMatrix omega_val	= omega_fnv["omega"];
	List	Inc_kernel 		= DP_list["Inc_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"]; 
	NumericMatrix Inc_weight= Inc_kernel["weight"], 
					policy	= DP_list["policy"];
	NumericVector policy_y 	= policy(_,1), policy_c = policy(_,2); 
	int n_act 				= 2, nsd = state.ncol(), D = DP_list["D"];
	struct paramStr param_out=param_assign(param); 
	gsl_set_error_handler_off();

	// Initiate the working variables; 
	int nhh = hh_index.length(), N = choice_seq.nrow(); 
	NumericMatrix state_out(N, nsd), choice_out(N, n_act+1); 
	state_out 	= clone(DataState);
	for(int i=0; i<n_act; i++){
		choice_out(_,i) = choice_seq(_,i); 
	}

	// Construct policy interpolation; 
	NumericVector x1 = unique(state(_,0)), x2 = unique(state(_,1)), x3 = unique(state(_,2));
	int N1 = x1.length(), N2 = x2.length(), N3 = x3.length();
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
	std::sort(x3.begin(), x3.end()); 

	// Initiate the c-spline; 
	int nspl = N2*N3;
	struct spl *y_spl	= new struct spl[nspl]; 
	struct spl *c_spl	= new struct spl[nspl]; 
	struct spl *V_spl	= new struct spl[nspl]; 
	for(int i=0; i<nspl; i++){
		IntegerVector index	= i*N1 + seq_len(N1) - 1;
		NumericVector tmpy 	= policy_y[index];
		NumericVector tmpc 	= policy_c[index];
		NumericVector tmpv	= V[index]; 		
		y_spl[i] = spl_init(x1, tmpy);
		c_spl[i] = spl_init(x1, tmpc);
		V_spl[i] = spl_init(x1, tmpv);
	}
	
	struct spl *omega_spl	= new struct spl[D]; 
	for(int i=0; i<D; i++){
		omega_spl[i] = spl_init(omega_y, omega_val(i,_));
	}
	struct DP_fun DPf = {omega_y, omega_val, omega_spl, Q_coef, Q_coef.length()}; 
	
	// Assign parameters for the value function; 
	int d = 1; 
	struct my_f_params params = {&d, param_out, state(0,_),beta, DPf, lnzero, 
								Inc_nodes, Inc_weight, x1, x2, x3, V_spl};	
		
	// Claim the original function; 
	gsl_multimin_function orig_func;
  	orig_func.n 	= n_act;
	orig_func.f 	= &my_valuefn;
	orig_func.params= &params;
	
	// Initiate policy function structure; 
	struct policy_fn myplc = {n_act, Inc_nodes, Inc_weight, DPf, x1, x2, x3, D, y_spl, c_spl, orig_func};

	// Loop over households; 
	int cnt = 0; 	
	if(simidx == 3){
		// Set up Gumbel distribution draw; 
		const gsl_rng_type * Tran;
		gsl_rng * rdraw;
		gsl_rng_env_setup();
		Tran 	= gsl_rng_default;
		rdraw 	= gsl_rng_alloc (Tran);
		NumericVector logit_draw(D+1); 
		
		for(int h=0; h<nhh; h++){
			// Assign initial states; 
			for(int j=0; j<nsd; j++){
				state_out(cnt,j) = init_state(h,j);
			}
					
			for(int i=cnt; i<cnt+TT_vec(h); i++){
				// Take logit draw; 
				if(logit_random){
					for(int j=0; j<D+1; j++){
						logit_draw(j) =  gsl_ran_gumbel1(rdraw, 1, 1); 
					}
				}else{
					logit_draw		= logit_draw_mat(i,_); 
				}
				params.state_vec	= state_out(i,_); 
				NumericVector out	= myplc.policySimall(state_out(i,_), logit_draw); 
				choice_out(i, 0) 	= out(0); 
				choice_out(i, 1)	= out(1);
				choice_out(i, 2)	= out(2);
				if(i-cnt < TT_vec(h)-1){
					state_out(i+1,0)	= out(3); 
				}
			}
			cnt = cnt + TT_vec[h]; 
		}
		
		// Free memory; 
		gsl_rng_free (rdraw);	
		for(int i=0; i<nspl; i++){
			y_spl[i].splfree();
			c_spl[i].splfree();
			V_spl[i].splfree(); 
		}
		for(int i=0; i<D; i++){
			omega_spl[i].splfree(); 
		}
		return Rcpp::List::create(_["DataState"] = state_out, _["choice_seq"] = choice_out);
	}else{
		// Also output the simulated sequence conditional on d only; 
		NumericMatrix state_out2(N, nsd), choice_out2(N, n_act+1); 
		state_out2	= clone(state_out); 
		choice_out2	= clone(choice_out); 

		for(int h=0; h<nhh; h++){
			// Assign initial states; 
			for(int j=0; j<nsd-1; j++){
				state_out(cnt,j) 	= init_state(h,j);
				state_out2(cnt,j) 	= init_state(h,j);
			}		
			for(int i=cnt; i<cnt+TT_vec(h); i++){
				NumericMatrix out	= myplc.policySim12(state_out(i,_), choice_out(i,0), choice_out(i,1), Q(i)); 
				choice_out(i,2)		= out(0,2);
				
				// Simulate y as well; 
				choice_out2(i,1)	= out(1,1); 
				choice_out2(i,2)	= out(1,2); 
				
				if(i-cnt < TT_vec(h)-1){
					state_out(i+1,0)	= out(0,3); 
					state_out2(i+1,0)	= out(1,3); 
				}
			}
			cnt = cnt + TT_vec[h]; 
		}
	
		// Free memory
		for(int i=0; i<nspl; i++){
			y_spl[i].splfree();
			c_spl[i].splfree();
			V_spl[i].splfree(); 
		}
		for(int i=0; i<D; i++){
			omega_spl[i].splfree(); 
		}
		delete[] y_spl; 
		delete[] c_spl;
		delete[] V_spl;
		delete[] omega_spl; 
		return Rcpp::List::create(_["DataState"] = state_out, _["choice_seq"] = choice_out, 
								  _["DataState_sim2"] = state_out2, _["choice_seq_sim2"] = choice_out2);	
	}
}

//------------//;
// Estimation //;
//------------//;
// [[Rcpp::export]]
NumericVector MM1C(NumericVector param, NumericMatrix choice_obs, NumericMatrix choice_seq_sim2, NumericMatrix DataState, IntegerVector TT_vec, List DP_list, Function logit_pred, NumericVector Q, NumericMatrix omega_dif, NumericMatrix phat, bool nest_pred = false){
/*
param			...		parameter vector; 
choice_obs		...		N*3 matrix, cbind(d, y, simulated c conditional on d and y)		
choice_seq_sim2	...		N*3 matrix, cbind(d, simulated y, simulated c conditional on d)
DataState		... 	N*4 matrix, cbind(simulated I from choice_obs, Inc, Rc, scaled lnZ)
TT_vec			...		nhh vector of total periods for each household
DP_list			...		DP solution 
logit_pred		... 	A function that returns matrix of ln(P(d) - P(D))
Q				...		precomputed vector of Q = Q_fn(y)
omega_dif		... 	precomputed matrix omega = omega_fn(y, d) and omega_dif(,i) = omega(,i) - omega(,D); 
*/
	// Assign parameter; 
	NumericMatrix state			= DP_list["state"]; 
	NumericVector x1			= unique(state(_,0)); 
	List Inc_kernel 			= DP_list["Inc_kernel"]; 
	NumericVector Inc_nodes 	= Inc_kernel["nodes"]; 
	int D						= DP_list["D"]; 
	int nhh 					= TT_vec.length(), N = DataState.nrow(), N1 = x1.length();
	struct paramStr param_out	= param_assign(param); 
	std::sort(x1.begin(), x1.end());
	
	// Initiate moments; 
	NumericVector my(N), md(N); 
	NumericMatrix md_mat(N,D-1);
	md.fill(0);  
	
	// Moments for y: difference between the observed y and simulated y; 
	my		= choice_obs(_,1) - choice_seq_sim2(_,1); 
	
	// Moments for discrete choice d: U(d,y,c) - U(D,y,c) - ln(P(d)) + ln(P(D)); 
	// Compute phat = ln(P(d)) - ln(P(D)); 
	if(nest_pred){
		NumericMatrix phat1 = logit_pred(DataState, choice_obs(_,0)); 
		phat				= phat1; 
	}
	
	NumericVector Inext		= DataState(_,0) + Q - choice_obs(_,2); 
	for(int i=1; i<D; i++){
		NumericVector tmp_md 	= omega_dif(_,i-1)  + (param_out.tau3) * Inext*(i-D) + (param_out.tau4)*(i-D) - phat(_,i-1); 
		tmp_md[choice_obs(_,0)==0] = 0; 
		md_mat(_,i-1) 			= tmp_md; 
	}

	// Construct instrument matrix; 
	int n_inst = 2; 
	NumericMatrix inst(N,n_inst); 
	inst.fill(1); 
	inst(_,1) = DataState(_,1); 
	
	// Multiply moments with instruments; 
	int n_M = n_inst*(D+1); 
	NumericMatrix M(N, n_M); 
	for(int i=0; i<n_inst; i++){
		M(_,i) 			= my*inst(_,i); 
		M(_,i+n_inst)	= pow(my,2) * inst(_,i);
		for(int j=1; j<D; j++){
			M(_, i+n_inst*(j+1)) = md_mat(_,j-1) * inst(_,i); 
		}
	}
	
	// Averaging the moment equations over observations; 
	NumericVector M_bar(n_M); 
	for(int i=0; i<n_M; i++){
		M_bar(i) = sum(M(_,i))/N; 
	}
	
	return(M_bar);
}

// [[Rcpp::export]]
List GMM_objC(NumericVector param, NumericVector hh_index, IntegerVector TT_vec, NumericMatrix init_state, NumericMatrix DataState, NumericMatrix choice_obs, NumericMatrix W, List DP_init, List control_list, Function logit_pred, NumericVector Q, NumericMatrix omega_dif, NumericMatrix phat, bool nest_pred = false){	
	// Check the boundary of parameters; 
	List DP_listout				= clone(DP_init);
	struct paramStr param_out 	= param_assign(param); 
	int n_M 					= W.ncol(), status; 
	double Qobj 				= 0;
	NumericVector M_bar(n_M); 
	
	// Find the DP solution; 	
	DP_listout["param"]			= param; 	
	DP_listout					= policy_iterC(DP_listout, control_list);
	status 						= DP_listout["status"];
	if(status == 3){
		Qobj	= NA_REAL;
	}else{
		M_bar.fill(0); 
		List seqout 				= sim_hhseqC(hh_index, TT_vec, init_state, DataState, choice_obs, DP_listout, control_list, Q); 
		NumericMatrix state_out1	= seqout["DataState"];
		NumericMatrix choice_out1	= seqout["choice_seq"];
		NumericMatrix choice_out_sim2= seqout["choice_seq_sim2"];
		M_bar 		= MM1C(param, choice_out1, choice_out_sim2, state_out1, TT_vec, DP_listout, logit_pred, Q, omega_dif, phat, nest_pred); 

		// Return the objective function; 
		for(int i=0; i<n_M; i++){
			for(int j=0; j<n_M; j++){
				Qobj = Qobj + M_bar(i)*W(i,j)*M_bar(j); 
			}
		} 	
	}
	return Rcpp::List::create(_["obj"] = Qobj, _["status"] = DP_listout["status"], _["Mbar"] = M_bar, _["value"] = DP_listout["value_fn"]); 
}

// [[Rcpp::export]]
NumericMatrix comp_W(NumericVector param, NumericVector hh_index, IntegerVector TT_vec, NumericMatrix init_state, NumericMatrix DataState, NumericMatrix choice_obs, NumericMatrix W, List DP_init, List control_list, Function logit_pred, NumericVector Q, NumericMatrix omega_dif, NumericMatrix phat){
	// Assign parameter; 
	NumericMatrix state			= DP_init["state"]; 
	NumericVector x1			= unique(state(_,0)); 
	List Inc_kernel 			= DP_init["Inc_kernel"]; 
	NumericVector Inc_nodes 	= Inc_kernel["nodes"]; 
	int D						= DP_init["D"]; 
	int nhh 					= TT_vec.length(), N = DataState.nrow(), N1 = x1.length();
	std::sort(x1.begin(), x1.end());
	
	List DP_listout				= clone(DP_init);
	struct paramStr param_out 	= param_assign(param); 
	int n_M 					= W.ncol(), status; 
	NumericMatrix M_mat_bar(DataState.nrow(), n_M); 
	M_mat_bar.fill(0); 
	
	// Find the DP solution; 	
	DP_listout["param"]			= param; 	
	DP_listout					= policy_iterC(DP_listout, control_list);	
	List seqout 				= sim_hhseqC(hh_index, TT_vec, init_state, DataState, choice_obs, DP_listout, control_list, Q); 
	NumericMatrix state_out1	= seqout["DataState"];
	NumericMatrix choice_out1	= seqout["choice_seq"];
	NumericMatrix choice_out_sim2= seqout["choice_seq_sim2"];
	
	// Initiate moments; 
	NumericVector my(N); 
	NumericMatrix md_mat(N,D-1);

	// Moments for y: difference between the observed y and simulated y; 
	my		= choice_obs(_,1) - choice_out_sim2(_,1); 	
	NumericVector Inext		= state_out1(_,0) + Q - choice_out1(_,2); 
	for(int i=1; i<D; i++){
		NumericVector tmp_md 	= omega_dif(_,i-1)  + (param_out.tau3) * Inext*(i-D) + (param_out.tau4)*(i-D) - phat(_,i-1); 
		tmp_md[choice_obs(_,0)==0] = 0; 
		md_mat(_,i-1) 			= tmp_md; 
	}

	// Construct instrument matrix; 
	int n_inst = 2; 
	NumericMatrix inst(N,n_inst); 
	inst.fill(1); 
	inst(_,1) = DataState(_,1); 

	// Multiply moments with instruments; 
	NumericMatrix M(N, n_M); 
	for(int i=0; i<n_inst; i++){
		M(_,i) 			= my*inst(_,i); 
		M(_,i+n_inst)	= pow(my,2) * inst(_,i);
		for(int j=1; j<D; j++){
			M(_, i+n_inst*(j+1)) = md_mat(_,j-1) * inst(_,i); 
		}
	}
	
	return(M); 
}

// [[Rcpp::export]]
List comp_policy(NumericMatrix DataState, List DP_list, List control_list){
	// Assign parameters; 
	NumericVector param 	= DP_list["param"], 
				  V 		= DP_list["value_fn"]; 
	double	lnzero 			= control_list["lnzero"], 
			beta 			= DP_list["beta"]; 
	NumericMatrix state 	= DP_list["state"]; 
	List omega_fnv			= DP_list["omega"];
	NumericVector omega_y	= omega_fnv["omega_y"], 
			Q_coef			= DP_list["Q_coef"]; 
	NumericMatrix omega_val	= omega_fnv["omega"];
	List	Inc_kernel 		= DP_list["Inc_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"]; 
	NumericMatrix Inc_weight= Inc_kernel["weight"], 
					policy	= DP_list["policy"];
	NumericVector policy_y 	= policy(_,1), policy_c = policy(_,2); 
	int n_act 				= 2, nsd = state.ncol(), D = DP_list["D"];
	struct paramStr param_out=param_assign(param); 
	gsl_set_error_handler_off();
			
	// Construct policy interpolation; 
	NumericVector x1 = unique(state(_,0)), x2 = unique(state(_,1)), x3 = unique(state(_,2));
	int N1 = x1.length(), N2 = x2.length(), N3 = x3.length();
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
	std::sort(x3.begin(), x3.end()); 

	// Initiate the c-spline; 
	int nspl = N2*N3;
	struct spl *y_spl	= new struct spl[nspl]; 
	struct spl *c_spl	= new struct spl[nspl]; 
	struct spl *V_spl	= new struct spl[nspl]; 
	for(int i=0; i<nspl; i++){
		IntegerVector index	= i*N1 + seq_len(N1) - 1;
		NumericVector tmpy 	= policy_y[index];
		NumericVector tmpc 	= policy_c[index];
		NumericVector tmpv	= V[index]; 		
		y_spl[i] = spl_init(x1, tmpy);
		c_spl[i] = spl_init(x1, tmpc);
		V_spl[i] = spl_init(x1, tmpv);
	}
	
	struct spl *omega_spl	= new struct spl[D]; 
	for(int i=0; i<D; i++){
		omega_spl[i] = spl_init(omega_y, omega_val(i,_));
	}
	struct DP_fun DPf = {omega_y, omega_val, omega_spl, Q_coef, Q_coef.length()}; 
	
	// Assign parameters for the value function; 
	int d = 1; 
	struct my_f_params params = {&d, param_out, state(0,_),beta, DPf, lnzero, 
								Inc_nodes, Inc_weight, x1, x2, x3, V_spl};	
		
	// Claim the original function; 
	gsl_multimin_function orig_func;
  	orig_func.n 	= n_act;
	orig_func.f 	= &my_valuefn;
	orig_func.params= &params;
	
	// Initiate policy function structure; 
	struct policy_fn myplc = {n_act, Inc_nodes, Inc_weight, DPf, x1, x2, x3, D, y_spl, c_spl, orig_func};
	
	// Initiate output variables; 
	int N 				= DataState.nrow();
	NumericMatrix ccp(N, D), choice_out(N, n_act+1); 
	NumericVector draw_norm(D+1), v(D); 
	draw_norm.fill(0); 
	gsl_vector *a; 
	a 					= gsl_vector_alloc(n_act); 
	double vmax			= 0, vsum = 0; 
	
	for(int i=0; i<N; i++){
		params.state_vec	= DataState(i,_); 
		NumericVector out	= myplc.policySimall(DataState(i,_), draw_norm); 
		choice_out(i, 0) 	= out(0); 
		choice_out(i, 1)	= out(1);
		choice_out(i, 2)	= out(2);
		gsl_vector_set(a, 0, choice_out(i, 1)); 
		gsl_vector_set(a, 1, choice_out(i, 2));
		
		for(int j=1; j<=D; j++){
			d 				= j; 
			v(j-1) 			= my_valuefn(a, &params) ;
		}
		vmax 				= max(v); 
		v					= exp(v - vmax); 
		vsum				= sum(v); 
		
		for(int j=0; j<D; j++){
			ccp(i,j)		= v(j)/vsum; 
		}		
	}
	
	// Free memory
	for(int i=0; i<nspl; i++){
		y_spl[i].splfree();
		c_spl[i].splfree();
		V_spl[i].splfree(); 
	}
	for(int i=0; i<D; i++){
		omega_spl[i].splfree(); 
	}
	delete[] y_spl; 
	delete[] c_spl;
	delete[] V_spl;
	delete[] omega_spl; 
	return Rcpp::List::create(_["policy"] = choice_out, _["ccp"] = ccp); 
}


// [[Rcpp::export]]
NumericMatrix plotvalueC(IntegerVector plotstate, NumericVector c_grid,List DP_list, List control_list){
	// Assign parameters; 
	NumericVector param 	= DP_list["param"], 
				  V 		= DP_list["value_fn"]; 
	double	lnzero 			= control_list["lnzero"], 
			beta 			= DP_list["beta"]; 
	NumericMatrix state 	= DP_list["state"]; 
	List omega_fnv			= DP_list["omega"];
	NumericVector omega_y	= omega_fnv["omega_y"], 
			Q_coef			= DP_list["Q_coef"]; 
	NumericMatrix omega_val	= omega_fnv["omega"];
	List	Inc_kernel 		= DP_list["Inc_kernel"]; 
	NumericVector Inc_nodes = Inc_kernel["nodes"]; 
	NumericMatrix Inc_weight= Inc_kernel["weight"], 
					policy	= DP_list["policy"];
	int n_act 				= 2, nsd = state.ncol(), D = DP_list["D"];
	struct paramStr param_out=param_assign(param); 
	gsl_set_error_handler_off();
			
	// Construct policy interpolation; 
	NumericVector x1 = unique(state(_,0)), x2 = unique(state(_,1)), x3 = unique(state(_,2));
	int N1 = x1.length(), N2 = x2.length(), N3 = x3.length();
	std::sort(x1.begin(), x1.end()); 
	std::sort(x2.begin(), x2.end()); 
	std::sort(x3.begin(), x3.end()); 

	// Initiate the c-spline; 
	int nspl = N2*N3;
	struct spl *V_spl	= new struct spl[nspl];
	for(int i=0; i<nspl; i++){
		IntegerVector index	= i*N1 + seq_len(N1) - 1;
		NumericVector tmpv	= V[index]; 		
		V_spl[i] = spl_init(x1, tmpv);
	}
	
	struct spl *omega_spl	= new struct spl[D]; 
	for(int i=0; i<D; i++){
		omega_spl[i] = spl_init(omega_y, omega_val(i,_));
	}
	struct DP_fun DPf = {omega_y, omega_val, omega_spl, Q_coef, Q_coef.length()}; 
	
	// Assign parameters for the value function; 
	int d = 1; 
	struct my_f_params params = {&d, param_out, state(0,_),beta, DPf, lnzero, 
								Inc_nodes, Inc_weight, x1, x2, x3, V_spl};	
	gsl_vector *a; 
	a = gsl_vector_alloc(n_act);
		
// Compute the value function at each states; 
	int sel=0, nps = plotstate.length(), nc=c_grid.length(); 
	NumericMatrix out(nps, nc+1); 
	for(int i=0; i<nps; i++){
		sel = plotstate(i) - 1; 
		params.state_vec = state(sel,_); 
		d = policy(sel, 0); 
		gsl_vector_set(a, 0, policy(sel,1) ); 
		
		for(int j=0; j<nc; j++){
			double Q 	= DPf.Q_fn(policy(sel,1)); 
			double lower_c = std::max(state(sel,0) + Q - x1[N1-1], 0.0);
			double upper_c = state(sel,0) + Q;
			if(c_grid(j) >= lower_c && c_grid(j) <= upper_c){
				gsl_vector_set(a, 1, c_grid(j) );
				out(i,j) = my_valuefn(a, &params) ;
			}
		}
		gsl_vector_set(a, 1, policy(sel, 2)); 
		out(i,nc) = my_valuefn(a, &params) ;
	}
	gsl_vector_free(a);
	for(int i=0; i<nspl; i++){
		V_spl[i].splfree(); 
	}
	for(int i=0; i<D; i++){
		omega_spl[i].splfree();
	}
	delete[] V_spl; 
	delete[] omega_spl; 
	return(out);		
}

