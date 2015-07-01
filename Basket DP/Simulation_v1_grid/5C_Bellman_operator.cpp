#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppGSL)]]

struct my_f_params { double I_current; double lambda; double tau1; double tau3;
					double beta; 
					double q; NumericVector yweight_vec;
					const gsl_spline *myint1; gsl_interp_accel *acc1; NumericVector x1;NumericVector y1; 
					const gsl_spline *myint2; gsl_interp_accel *acc2; NumericVector x2;NumericVector y2;
					double slope1_first; double slope1_end; double slope2_first; double slope2_end; };

double my_valuefn( double cs, void *p){
    struct my_f_params *params = (struct my_f_params *)p;	
    double lambda				= (params->lambda);
    double tau1					= (params->tau1);
    double tau3					= (params->tau3);
	double beta					= (params->beta);
	double q					= (params->q);
	double I_current			= (params->I_current);
	NumericVector yweight_vec	= (params->yweight_vec);
	const gsl_spline *myint1    = (params->myint1);
    gsl_interp_accel *acc1      = (params->acc1);
    NumericVector x1           	= (params->x1);
    NumericVector y1            = (params->y1);
   	const gsl_spline *myint2    = (params->myint2);
    gsl_interp_accel *acc2      = (params->acc2);
    NumericVector x2           	= (params->x2);
    NumericVector y2            = (params->y2);
    double slope1_first			= (params->slope1_first);
    double slope1_end			= (params->slope1_end);
	double slope2_first			= (params->slope2_first);
    double slope2_end			= (params->slope2_end);
    int k = x1.length();
    	
	double I_next=0, out = 0, ev = 0;
	NumericVector v(2);
	I_next = I_current + q -cs;
	
	if(I_next>=x1(0) & I_next<=x1(k-1)){
		v(0) = gsl_spline_eval(myint1, I_next, acc1);
		v(1) = gsl_spline_eval(myint2, I_next, acc2);
	}else if(I_next < x1(0)){
		v(0) = y1(0) + slope1_first *(I_next - x1(0));
		v(1) = y2(0) + slope2_first *(I_next - x2(0));
	}else {
		v(0) = y1(k-1) + slope1_end *(I_next - x1(k-1));
		v(1) = y2(k-1) + slope2_end *(I_next - x2(k-1));
	}
	ev = sum(v*yweight_vec);
	out = lambda*log(cs+0.01) - tau1*I_next - tau3*q + beta*ev;
	return(-out);
}

// [[Rcpp::export]]
List Solve_c_fnC(NumericMatrix state, NumericVector x1, NumericVector V1, NumericVector x2, NumericVector V2,  
				int K, NumericVector param, double beta, List Q_kernal, List y_kernal, List control_list){
	const double lambda = param(0), tau1 = param(1), tau3 = param(2);
	NumericVector Q_grid	= Q_kernal[0], y_nodes = y_kernal[0];
	NumericMatrix yweight 	= y_kernal[1];
	List Q_weight 			= Q_kernal[1];
	double Imax 			= max(x1);
	int ns = state.nrow(), nQ = Q_grid.length(), ny = y_nodes.length(), nx = x1.length();

	// The object to store computation values; 
	NumericMatrix cstar(ns, nQ), Vc(ns, nQ), Vk(ns, K);
	NumericVector policy_k(ns), value_fn(ns);
	
	// Linear interpolation set up;
	gsl_interp_accel *accP1	= gsl_interp_accel_alloc() ;
	gsl_interp_accel *accP2	= gsl_interp_accel_alloc() ;
	gsl_spline *spline1 	= gsl_spline_alloc( gsl_interp_cspline , nx );
	gsl_spline *spline2 	= gsl_spline_alloc( gsl_interp_cspline , nx );
	gsl_spline_init( spline1, x1.begin(), V1.begin(), nx);
	gsl_spline_init( spline2, x2.begin(), V2.begin(), nx);
	
	// Compute the slope of linear extrapolating line;
	double eps = 0.0001;
	double slope1_first = (gsl_spline_eval (spline1, x1[0]+eps, accP1) - V1[0])/eps;
	double slope1_end 	= (V1[nx-1] - gsl_spline_eval (spline1, x1[nx-1]-eps, accP1))/eps;
	double slope2_first = (gsl_spline_eval (spline2, x2[0]+eps, accP1) - V2[0])/eps;
	double slope2_end 	= (V2[nx-1] - gsl_spline_eval (spline2, x2[nx-1]-eps, accP2))/eps;
	struct my_f_params params = {state(0,0), lambda, tau1, tau3, beta, Q_grid(0), yweight(0,_), 
								spline1, accP1, x1, V1, spline2, accP2, x2, V2, 
								slope1_first, slope1_end, slope2_first, slope2_end};
										
	// Optimizer setup;
	int status;
	int iter=0, max_iter = control_list["brent_max_iter"];
	double brent_tol	 = control_list["brent_tol"];
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;
	gsl_set_error_handler_off(); //TURN OFF GSL ERROR ERROR HANDLER (GSL JUST PRINT ERROR MESSAGE AND KILL THE PROGRAM IF FLAG IN ON)

	// GSL objective function ;
	gsl_function F; 
	F.function 	= &my_valuefn;
	F.params	= &params;
	
	//Create optimizer;
	T = gsl_min_fminimizer_quad_golden;
	s = gsl_min_fminimizer_alloc(T);
	double low = 0, up = Imax;
	double m = 0.2, v= 0;
	int sat;

	// Solve for c for every state in a loop;
	for(int i=0; i < ns; i++){
		int sel_y = state(i,1)-1;
		for(int j=0; j < nQ; j++){
			params.I_current	= state(i,0);
			params.q			= Q_grid(j);
			params.yweight_vec	= yweight(sel_y,_);
			if(state(i,0) + Q_grid(j)<=0){
				cstar(i,j) = 0;
				Vc(i,j) = -my_valuefn(0, &params); 
			}else{
				// Re-initiate minimizer and parameters;
				low		= 0;
				up		= state(i,0) + Q_grid(j);
				iter	= 0;
				
				// Try initial values close to the boundaries;
				m		= up - brent_tol;
				sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
				if(sat == GSL_EINVAL){
					m= brent_tol;
					sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
				}
// 				int m_iter	= 1;			
// 				do{
// 						m = m_iter*0.05*(state(i,0) + params.q);
// 						sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
// 						m_iter++;
// 					}
// 				while(sat == GSL_EINVAL && m_iter < 20);

				if(sat == GSL_EINVAL){
					cstar(i,j) =  up;
					Vc(i,j) =  -my_valuefn(cstar(i,j), &params);	 
				}else{
					//Optimization loop;
					do {
					iter++;
					status = gsl_min_fminimizer_iterate(s);
					m      = gsl_min_fminimizer_x_minimum(s);
					low    = gsl_min_fminimizer_x_lower(s);
					up     = gsl_min_fminimizer_x_upper(s);
					v	   = gsl_min_fminimizer_f_minimum(s);

					status = gsl_min_test_interval(low,up, brent_tol, 0.0);
					}	
					while (status == GSL_CONTINUE && iter < max_iter);
					cstar(i,j) = m;
					Vc(i,j) = -v;
				}			
			}
		}
		
		// Construct choice-specific value function Vk;
		NumericMatrix Qw = Q_weight[sel_y];
		for(int k = 0; k<K; k++){
			Vk(i,k) = sum(Vc(i,_)*Qw(k,_));
		}
		double value_max = max(Vk(i,_));
		for(int k = 0; k<K; k++){
			if(Vk(i,k)==value_max) { policy_k(i) = k+1;}
		}
		value_fn(i) = log(sum(exp(Vk(i,_) - value_max))) + value_max;
	}
	gsl_min_fminimizer_free(s);
    gsl_interp_accel_free (accP1);
    gsl_interp_accel_free (accP2);
	gsl_spline_free (spline1);
	gsl_spline_free (spline2);
	return Rcpp::List::create(_["policy"] = List::create(_["k"] = policy_k, _["c"] = cstar), 
							  _["value"]=value_fn, _["Vk"] = Vk, _["Vc"] = Vc);
}

// [[Rcpp::export]]
List Bellman_operatorC(List DP_list, List control_list){
	NumericMatrix state = DP_list["state"];
	NumericVector V = DP_list["value_fn"], state_y = state(_,1), state_I = state(_,0), param = DP_list["param"];
	NumericVector I_grid = state_I[state_y==1], V1 = V[state_y==1], V2 = V[state_y==2];
	int K = DP_list["K"];
	double	beta = DP_list["beta"];
	List Q_kernal = DP_list["Q_kernal"], y_kernal = DP_list["y_kernal"];
	List out;
	out  = Solve_c_fnC(state, I_grid, V1, I_grid, V2, K, param, beta, Q_kernal, y_kernal, control_list);
	return(out);
}

// [[Rcpp::export]]
List value_iteration_fnC(List DP_list, int print_level, Function Bellman_operator, List control_list){	
	double tol = control_list["tol"]; 
	int iter = 0, max_iter = control_list["value_max_iter"], display_freq = control_list["display_freq"];
	double norm_dif = 100;
	List DP_listout = clone(DP_list);

	do{
		iter++;
		NumericVector W_old 	= DP_listout["value_fn"];
		List sol				= Bellman_operator(DP_listout, control_list);
		NumericVector W_new		= sol["value"];
		List policy				= sol["policy"];
		norm_dif				= max(abs(W_new - W_old));
		DP_listout["value_fn"]	= W_new;
		DP_listout["policy"]	= policy;
		if(print_level == 1){
			if(iter%display_freq ==0 ){
				Rcout<<"Iteration "<<iter<<", norm_diff = "<<norm_dif<<"\n";
			}
		}
	}
	while(norm_dif > tol && iter < max_iter);
	DP_listout["Iteration"] = iter;
	if(norm_dif <= tol){
		DP_listout["status"] = 3;
		Rcout<<"Value function iteration has converged.\n";
	}else{
		DP_listout["status"] = 2;	
		Rcout<<"Value function iteration did not converge.\n";
	}
	return(DP_listout);
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
List CCP_fnC(NumericMatrix DataState, NumericVector x1, NumericVector V1, NumericVector x2, NumericVector V2,  
				int K, NumericVector param, double beta, List Q_kernal, List y_kernal, bool makedraw, List control_list){
	List Bellman_list = Solve_c_fnC(DataState, x1, V1, x2, V2, K, param, beta, Q_kernal, y_kernal, control_list);
	NumericMatrix Vk = Bellman_list["Vk"];
	List policy = Bellman_list["policy"];
	int ns = Vk.nrow();
	NumericVector u(K);
	NumericMatrix ccp(ns, K);
	IntegerVector draw(ns);
	
	for(int i = 0; i<ns; i++){
		NumericVector v = Vk(i,_);
		u = clone(v);
		u = exp(u - max(u));
		ccp(i,_) = u/sum(u);
		if(makedraw){
			draw(i) = my_rmultinorm(ccp(i,_));
		}
	}
	return Rcpp::List::create(_["ccp"] = ccp, _["policy"] = policy, _["draw"] = draw );
}

// [[Rcpp::export]]
NumericVector ll_fnC(NumericVector param, NumericMatrix DataState, IntegerVector choice_seq, 
					List DP_init, Function Bellman_operator, List control_list){
	// Construct DP_list and Bellman operator; 
	List DP_list = clone(DP_init);
	double deriv = sum(param);
// 	NumericVector param1 = clone(param);
// 	param1 = exp(param1);
	DP_list["param"] = param; 
	NumericMatrix state = DP_list["state"];
	NumericVector state_y = state(_,1), state_I = state(_,0);
	NumericVector I_grid = state_I[state_y==1];
	int K = DP_list["K"];
	double	beta = DP_list["beta"];
	List Q_kernal = DP_list["Q_kernal"], y_kernal = DP_list["y_kernal"];
	int ns = DataState.nrow(), choice_index = 0;
	NumericVector ll(ns);
	
	if(param[0] <=0 || param[1] <=0 ){
		for(int i=0; i<ns; i++){
			ll(i) = -100;
		}
		return(ll);
	}else{
		// Solve Dynamic programing at the given parameter; 
		List sol = value_iteration_fnC(DP_list, 0, Bellman_operator, control_list);
		int status = sol["status"];
		if(status == 2){
// 			return NumericVector::create(NA_REAL);
			return(ll);
		}else{
			NumericVector V = sol["value_fn"]; 
			NumericVector V1 = V[state_y==1], V2 = V[state_y==2];

			// Compute CCP at each state; 
			List ccp_ls = CCP_fnC(DataState, I_grid, V1, I_grid, V2, K, param, beta, Q_kernal, y_kernal, FALSE, control_list);
			NumericMatrix ccp = ccp_ls["ccp"];

			// Compute log-likelyhood of each state and choice data;
			for(int i=0; i<ns; i++){
				choice_index = choice_seq[i] - 1;
				ll(i) = log(ccp(i, choice_index)) ;
			}	
			return(ll);
		}
	}
}
