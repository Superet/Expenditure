#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]

//------------------//;
// Define structure //;
//------------------//;
int findInterval1(double x, NumericVector breaks);

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

// [[Rcpp::export]]
int findInterval1(double x, NumericVector breaks) {
	int out;
	out = std::distance(breaks.begin(), std::upper_bound(breaks.begin(), breaks.end(), x));
	return out;
}

struct my_f_params {NumericVector param_vec; double beta; 
					double I_current; double lny; double Inc; double Rc;
					int k; double q; double TR; double omega; 	
					NumericVector ynodes; NumericVector yweight; double price_min;
					NumericVector kappa; double sigma; double rho; NumericVector Inc_weight;
					NumericVector x1; NumericVector x2; NumericVector x3; NumericVector x4; struct spl *spl2d;};

//---------------------------//;
// Value function definition //;
//---------------------------//;
double my_valuefn(double cs, void *p){
	struct my_f_params *params = (struct my_f_params *)p;
	NumericVector param_vec		= (params->param_vec);	
	double beta					= (params->beta);
	double I_current			= (params->I_current);
	double lny					= (params->lny);
	double Inc					= (params->Inc);
	double Rc					= (params->Rc);
	int k						= (params->k);
	double q 					= (params->q);
	double TR 					= (params->TR);
	double omega				= (params->omega);
	NumericVector yweight		= (params->yweight);
	NumericVector ynodes		= (params->ynodes);
	double price_min			= (params->price_min);
	NumericVector kappa			= (params->kappa);
	double sigma				= (params->sigma);
	double rho					= (params->rho);
	NumericVector Inc_weight	= (params->Inc_weight);
	NumericVector x1			= (params->x1);
	NumericVector x2			= (params->x2);
	NumericVector x3			= (params->x3);
	NumericVector x4			= (params->x4);
	struct spl *spl2d			= (params->spl2d);
	
	// Assign parameters; 
	double lambda = param_vec(0), tau_b = param_vec(1), tau1 = param_vec(2), 
			tau2 = param_vec(3);
	
	//Initiate variables; 
	int  n 			= ynodes.length(), N2 = x2.length(), N3 = x3.length(), 
		n_Inc = Inc_weight.length();
	NumericVector lny_next(n), v(n), y2(N2), ev_Inc(n_Inc);
	double ev = 0;
	gsl_set_error_handler_off();
	
	// Assign values to the next state; 
	NumericVector  s_next(4);
	s_next(0) 	= I_current + q - cs; 
	s_next(1) 	= 0;			// Will change; 
	s_next(2) 	= Inc; 			// Will change; 
	s_next(3) 	= Rc; 
	
	// Compute the expected future value; 
	for(int i=0; i<n_Inc; i++){
		// Update the state values;
		s_next(2) 	= i + 1; 
	
		// Compute the location of next y; 
		double mu 	= kappa(i);
		lny_next	= sqrt(2)*sigma*ynodes + mu + rho*lny;
		
		// Interpolate the option value along the first dimension; 
		int x3_index = findInterval1(s_next(2), x3);		//Find which discrete value it is in the 3rd dimension; 
		int x4_index = findInterval1(s_next(3), x4);		//Find which discrete value it is in the 4th dimension; 	
		int index = (x4_index -1)*N2*N3 + (x3_index -1)*N2;
		for(int j=0; j<N2; j++){		
			y2(j) = spl2d[index+j].splEval(s_next(0));
		}
	
		// Create the spline in the second dimension;
		struct spl myspl = spl_init(x2, y2);
		for(int j=0; j<n; j++){
			v[j] = myspl.splEval(lny_next[j]);
		}
		myspl.splfree();
		ev_Inc(i) = sum(yweight*v)/sqrt(PI);
	}
	ev = sum(ev_Inc * Inc_weight);
	
	// Borrowing cost
	double borrow_cost = 0;
	if(exp(lny)/price_min <q){
		borrow_cost = tau_b * log(price_min*q - exp(lny));
	}
		
	// Compute the value function; 
	double out = lambda*log(cs+0.01) - tau1*s_next(0) - tau2*TR + omega - borrow_cost + beta*ev ;
	return(-out);
}


// [[Rcpp::export]]
NumericMatrix cal_valuefn(List DP_list, NumericVector c_grid, NumericMatrix DataState, int choice_k, bool optimalc, NumericVector cstar){
	NumericMatrix state = DP_list["state"], Inc_weight = DP_list["Inc_weight"];
	NumericVector	V = DP_list["value_fn"], param = DP_list["param"], Q_grid = DP_list["Q_grid"], 
				 	TR_grid = DP_list["TR_grid"], price=DP_list["price"];
	double	beta = DP_list["beta"], price_min = min(price);
 	List y_kernal = DP_list["y_kernal"];
	Function omega_fn	= DP_list["omega"];
	
	NumericVector 	ynodes = y_kernal["nodes"], yweight = y_kernal["weight"], 
					ysigma = y_kernal["sigma"], yrho = y_kernal["rho"];
	NumericMatrix ykappa = y_kernal["kappa"];
	int n_Inc = Inc_weight.ncol(), ns1 = DataState.nrow(), nc = c_grid.length();
	gsl_set_error_handler_off();
	NumericMatrix out(ns1, nc+1);
	
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
	struct my_f_params params = {param, beta, 
								 state(1,0), state(1,1), state(1,2), state(1,3),
								 choice_k, Q_grid(choice_k-1), TR_grid(choice_k-1), 0, 
								 ynodes, yweight, price_min,
								 ykappa(_,0), ysigma(0), yrho(0), Inc_weight(0,_),
								 x1, x2, x3, x4, myspl};
	for(int i=0; i<ns1; i++){
		for(int j=0; j<nc; j++){
			double tmp_omega 	= as<double>(omega_fn(choice_k, DataState(i,1)));
			int Rc_int 			= DataState(i,3); 
			params.I_current	= DataState(i,0);
			params.lny			= DataState(i,1);
			params.Inc			= DataState(i,2);
			params.Rc			= DataState(i,3);
			params.omega		= tmp_omega;
			params.kappa		= ykappa(_, Rc_int);
			params.sigma		= ysigma(Rc_int);
			params.rho			= yrho(Rc_int);
			params.Inc_weight	= Inc_weight(n_Inc*Rc_int + DataState(i,2) -1 ,_);
			out(i,j) 			= -my_valuefn(c_grid(j), &params);
		}
		if(optimalc){
			if(cstar(i) == NA_REAL){
				out(i, nc) = NA_REAL;
			}else{
				out(i, nc) = -my_valuefn(cstar(i), &params);
			}
		}
	}	
	for(int i=0; i<nspl; i++){
		myspl[i].splfree();
	}

	return(out);
}

// [[Rcpp::export]]
NumericVector V_initC(NumericVector param, List DP_list, IntegerVector policy_k, NumericVector policy_c){
	NumericMatrix state = DP_list["state"]; 
	NumericVector price = DP_list["state"], Q_grid = DP_list["Q_grid"], TR_grid = DP_list["TR_grid"];
	double lambda = param(0), tau_b = param(1), tau1 = param(2), price_min = min(price), 
			tau2 = param(3),
			beta = DP_list["beta"], borrow_cost = 0; 
	int ns = state.nrow();
	Function omega_fn = DP_list["omega"];
	NumericVector out(ns);
	
	for(int i=0; i<ns; i++){
// 		double tripcost = param(2 + policy_k(i)); 
		double q = Q_grid(policy_k(i)-1), TR = TR_grid(policy_k(i)-1), lny = state(i,1);
		double I_next = state(i,0) + q - policy_c(i);
		double omega = as<double>(omega_fn(policy_k(i), lny));
		if(exp(lny)/price_min <q){
			borrow_cost = tau_b * log(price_min*q - exp(lny));
		}	
		out(i) = lambda*log(policy_c(i)+0.01) - tau1*I_next - tau2*TR - borrow_cost + omega; 
	} 
	return(out/(1-beta));
}

//------------------//;
// Bellman operator //;
//------------------//;
// [[Rcpp::export]]
List Solve_DP_fnC(NumericMatrix state, NumericVector V_old, NumericMatrix omega_mat, int K,
				  NumericVector param, double beta, NumericVector Q_grid, NumericVector TR_grid, 
				  NumericMatrix Inc_weight, List y_kernal, NumericVector price, List control_list){
	NumericVector 	ynodes = y_kernal["nodes"], yweight = y_kernal["weight"], 
					ysigma = y_kernal["sigma"], yrho = y_kernal["rho"];
	NumericMatrix ykappa = y_kernal["kappa"];
	int ns = state.nrow(), n_Inc = Inc_weight.ncol();
	double price_min	= min(price), Imax = max(state(_,0));
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
	struct my_f_params params = {param, beta, 
								 state(1,0), state(1,1), state(1,2), state(1,3),
								 1, Q_grid(1), TR_grid(1), omega_mat(0,0), 
								 ynodes, yweight, price_min,
								 ykappa(_,0), ysigma(0), yrho(0), Inc_weight(0,_),
								 x1, x2, x3, x4, myspl};
										
	// Optimizer setup;
	int status;
	int iter=0, max_iter = control_list["brent_max_iter"];
	double brent_tol	 = control_list["brent_tol"];
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;

	// GSL objective function ;
	gsl_function F; 
	F.function 	= &my_valuefn;
	F.params	= &params;
	
	//Create optimizer;
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	double low = 0, up = state(0,0);
	double m = 0.2, v= 0;
	int sat;
	
	// The object to store computation values; 
	NumericMatrix cstar(ns, K), vc(ns, K), ccp(ns, K);
	NumericVector V_new(ns);
					
	// Solve for c for every state in a loop;
	for(int i=0; i<ns; i++){
		for(int k=0; k<K; k++){
			int Rc_int 			= state(i,3); 
			params.I_current	= state(i,0);
			params.lny			= state(i,1);
			params.Inc			= state(i,2);
			params.Rc			= state(i,3);
			params.k			= k + 1; 
			params.q			= Q_grid(k);
			params.TR			= TR_grid(k);
			params.omega		= omega_mat(i,k);
			params.kappa		= ykappa(_, Rc_int);
			params.sigma		= ysigma(Rc_int);
			params.rho			= yrho(Rc_int);
			params.Inc_weight	= Inc_weight(n_Inc*Rc_int + state(i,2) -1 ,_);
			
			if(state(i,0)+Q_grid[k] <= 0){
				cstar(i,k)	= 0;
				vc(i,k)		= -my_valuefn(0, &params);			
			}else{
				// Re-initiate minimizer and parameters;
				if(state(i,0) + params.q - Imax > 0){
					low 	= state(i,0) + params.q - Imax; 
				}else{
					low		= 0; 
				}
				up			= state(i,0) + params.q;
				iter		= 0;
				// Try initial values close to the boundaries;
				m		= up - brent_tol;
				sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
				if(sat == GSL_EINVAL){
					m= low + brent_tol;
					sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
				}

				if(sat == GSL_EINVAL){
					cstar(i,k) 	=  up;
					vc(i,k) 	=  -my_valuefn(up, &params);	 
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
					cstar(i,k) = m;
					vc(i,k) = -v;
				}
			}
		}
		double vmax = max(vc(i,_));
		NumericVector tmpu = vc(i,_);
		NumericVector tmpv = clone(tmpu);
		tmpv = exp(tmpv- vmax);
		V_new(i) = log(sum(tmpv)) + vmax;
		ccp(i,_) = tmpv/sum(tmpv);
	}
	
	// Release memory of GSL pointers;
	gsl_min_fminimizer_free(s);
	for(int i=0; i<nspl; i++){
		myspl[i].splfree();
	}

	return Rcpp::List::create(_["policy"] = List::create(_["c"] = cstar), 
							  _["value"]=V_new, _["ccp"] = ccp, _["vc"] = vc);
}

// [[Rcpp::export]]
List Bellman_operatorC(List DP_list, List control_list, bool is_data=false, NumericMatrix DataState = NumericMatrix(2,2) ){
	NumericMatrix state = DP_list["state"], Inc_weight = DP_list["Inc_weight"];
	NumericVector V_old = DP_list["value_fn"], param = DP_list["param"], Q_grid = DP_list["Q_grid"], 
				  TR_grid = DP_list["TR_grid"], price=DP_list["price"];
	int K = DP_list["K"];
	double	beta = DP_list["beta"];
 	List y_kernal = DP_list["y_kernal"];
	Function omega_fn = DP_list["omega"];
	
	List out;
	if(is_data){
		int ns = DataState.nrow();
		NumericMatrix omega(ns, K);
		for(int k=0; k<K; k++){
			NumericVector tmp = omega_fn(k+1, DataState(_,1));
			omega(_,k) = tmp;
		}
		out  = Solve_DP_fnC(DataState, V_old, omega, K, param, beta, Q_grid, TR_grid, Inc_weight, y_kernal, price, control_list);
	}else{
		int ns = state.nrow();
		NumericMatrix omega(ns, K);
		for(int k=0; k<K; k++){
			NumericVector tmp = omega_fn(k+1, state(_,1));
			omega(_,k) = tmp;
		}
		out  = Solve_DP_fnC(state, V_old, omega, K, param, beta, Q_grid, TR_grid, Inc_weight, y_kernal, price, control_list); 
	}
	return(out);
}

//--------------------------//;
// Value function iteration //;
//--------------------------//;
// [[Rcpp::export]]
List value_iteration_fnC(List DP_list, Function Bellman_operator, int print_level, List control_list){	
	double tol = control_list["tol"]; 
	int iter = 0, max_iter = control_list["value_max_iter"], display_freq = control_list["display_freq"];
	double norm_dif = 10;
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
		if(print_level == 1){
			Rcout<<"Value function iteration has converged.\n";
		}
	}else{
		DP_listout["status"] = 2;	
		//if(print_level == 1){
			Rcout<<"Value function iteration did not converge.\n";
		//}
	}
	return(DP_listout);
}

//-------------------------------------//;
// Simulate a sequence of dynamic path //;
//-------------------------------------//;
int my_rmultinorm(NumericVector probs) {
    int k = probs.size();
    IntegerVector ans(k);
    rmultinom(1, probs.begin(), k, ans.begin());
    IntegerVector index = seq_len(k);
    IntegerVector out = index[ans==1];
    return(out[0]);
}

// [[Rcpp::export]]
List simulate_seq1_fnC( NumericVector init_state, int TT, bool draw_all, 
						NumericVector lny_seq, NumericVector Inc_seq, NumericVector Rc_seq,
						IntegerVector k_seq, 
						List DP_list, List control_list){
	RNGScope scope;		// ensure RNG gets set/reset;
	NumericMatrix state = DP_list["state"], Inc_weight = DP_list["Inc_weight"];
	NumericVector	V = DP_list["value_fn"], param = DP_list["param"], Q_grid = DP_list["Q_grid"], 
				 	TR_grid = DP_list["TR_grid"], price=DP_list["price"];
	int K = DP_list["K"];
	double	beta = DP_list["beta"], price_min = min(price), Imax = max(state(_,0));
 	List y_kernal = DP_list["y_kernal"];
	Function omega_fn = DP_list["omega"];
	int ms = state.ncol(), n_Inc = Inc_weight.ncol();
	
	NumericVector 	ynodes = y_kernal["nodes"], yweight = y_kernal["weight"], 
					ysigma = y_kernal["sigma"], yrho = y_kernal["rho"];
	NumericMatrix ykappa = y_kernal["kappa"];
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
		NumericVector tmpy 	= V[index];
		myspl[i] = spl_init(x1, tmpy);
	}
	
	// Assign parameters for the value function; 
	struct my_f_params params = {param, beta, 
								 state(1,0), state(1,1), state(1,2), state(1,3),
								 1, Q_grid(1), TR_grid(1), 0, 
								 ynodes, yweight, price_min,
								 ykappa(_,0), ysigma(0), yrho(0), Inc_weight(0,_),
								 x1, x2, x3, x4, myspl};
																				
	// Optimizer setup;
	int status;
	int iter=0, max_iter = control_list["brent_max_iter"];
	double brent_tol	 = control_list["brent_tol"];
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;

	// GSL objective function ;
	gsl_function F; 
	F.function 	= &my_valuefn;
	F.params	= &params;
	
	//Create optimizer;
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	double low = 0, up = state(0,0);
	double m = 0.2, v= 0;
	int sat;
	
	// Initiate result storage;
	NumericMatrix DataState(TT, state.ncol()), ccp(TT, K);
	NumericVector c_seq(TT), ccp_vec(TT), Q_seq(TT);
	
	// Initiate the states; 
	for(int i=0; i<ms; i++){
		DataState(0,i) = init_state(i);
	}
	if(!draw_all){
		DataState(_,1) = lny_seq; 
		DataState(_,2) = Inc_seq; 
		DataState(_,3) = Rc_seq; 
	}
	
	// Solve for q, c for every period in a loop;
	for(int i=0; i<TT; i++){
		// Working objects;
		NumericVector vc(K), cstar(K);
		
		// Obtain ccp function;
		for(int k=0; k<K; k++){
			double tmp_omega 	= as<double>(omega_fn(k+1, DataState(i,1)));
			int Rc_int 			= DataState(i,3); 
			params.I_current	= DataState(i,0);
			params.lny			= DataState(i,1);
			params.Inc			= DataState(i,2);
			params.Rc			= DataState(i,3);
			params.k			= k + 1; 
			params.q			= Q_grid(k);
			params.TR			= TR_grid(k);
			params.omega		= tmp_omega;
			params.kappa		= ykappa(_, Rc_int);
			params.sigma		= ysigma(Rc_int);
			params.rho			= yrho(Rc_int);
			params.Inc_weight	= Inc_weight(n_Inc*Rc_int + DataState(i,2) -1 ,_);
	
			if(DataState(i,0)+Q_grid(k) <= 0){
				cstar(k)	= 0;
				vc(k)		= -my_valuefn(0, &params);
			}else{
				// Re-initiate minimizer and parameters;
				if(DataState(i,0) + params.q - Imax >0){
					low 	= DataState(i,0) + params.q - Imax;
				}else{
					low 	= 0; 
				}
				up			= DataState(i,0) + params.q;
				iter		= 0;
				// Try initial values close to the boundaries;
				m		= up - brent_tol;
				sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
				if(sat == GSL_EINVAL){
					m= low + brent_tol;
					sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
				}

				if(sat == GSL_EINVAL){
					cstar(k) =  up;
					vc(k) =  -my_valuefn(up, &params);	 
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
					cstar(k) = m;
					vc(k) = -v;
				}
			}
		}
		double vmax = max(vc);
		NumericVector tmpv = exp(vc - vmax); 
		ccp(i,_) = tmpv/sum(tmpv);		
		
		if(draw_all){
			k_seq(i) 	= my_rmultinorm(ccp(i,_));
			Q_seq(i)	= Q_grid(k_seq(i)-1);
		}
		ccp_vec(i) 	= ccp(i,k_seq[i] - 1);
		c_seq(i)	= cstar(k_seq[i] - 1);
		
		// Evolve to the next state;
		if(i<=TT-2){
			DataState(i+1, 0) = DataState(i,0) + Q_grid(k_seq[i] - 1) - cstar(k_seq[i] - 1);			
			if(draw_all){
				DataState(i+1,3) = DataState(i, 3); 	// Assuming recession is also constant; 
				int Rc_int 	= DataState(i+1, 3);
				DataState(i+1,2) = my_rmultinorm(Inc_weight(n_Inc*Rc_int + DataState(i,2) -1 ,_)); 	
				int s_int 	= DataState(i+1, 2);
				double mu	= ykappa(s_int - 1, Rc_int);	
				DataState(i+1, 1) =  rnorm(1, mu + yrho(Rc_int)*DataState(i,1), ysigma(Rc_int))[0] ;					
			}
		}
	}
	// Release memory of GSL pointers;
	gsl_min_fminimizer_free(s);
	for(int i=0; i<nspl; i++){
		myspl[i].splfree();
	}
	
	return Rcpp::List::create(Rcpp::DataFrame::create(_["t"] = seq_len(TT), _["I"]=DataState(_,0),
								_["lny"] = DataState(_,1), _["Inc"] = DataState(_,2), _["Rc"] = DataState(_,3), 
								_["k"] = k_seq, _["Q"]=Q_seq, _["c"]=c_seq), 
		    				  	_["ccp"]=ccp, _["ccp_vec"] = ccp_vec);
}

// [[Rcpp::export]]
List simulate_hhseq1_fnC(NumericVector hh_index, IntegerVector TT_vec, NumericMatrix init_state,
						NumericVector lny_seq, NumericVector Inc_seq, NumericVector Rc_seq,
					    IntegerVector k_seq, 
						List DP_list, List control_list){
	int nhh = hh_index.length(), N = k_seq.length(), K = DP_list["K"];
	NumericMatrix ccp(N,K);
	NumericVector ccp_vec(N);
	
	// Loop over all the household; 
	int cnt = 0;				// Record the position of the last record for an household; 
	for(int i=0; i<nhh; i++){
		IntegerVector index = seq_len(TT_vec[i]) + cnt - 1;
		List sim				= simulate_seq1_fnC(init_state(i,_), TT_vec(i), false, lny_seq[index], 
									Inc_seq[index], Rc_seq[index], k_seq[index], DP_list, control_list);
		NumericMatrix tmp_ccp	= sim["ccp"];
		NumericVector tmp_ccpv 	= sim["ccp_vec"]; 
// 		ccp(index,_) 			= tmp_ccp;
		ccp_vec[index]			= tmp_ccpv;
		cnt			= cnt + TT_vec[i];
	}
	return Rcpp::List::create( _["ccp_vec"] = ccp_vec); 							
}


// [[Rcpp::export]]
List simulate_hhseq_fnC(NumericVector hh_index, IntegerVector TT_vec, NumericMatrix init_state,
						NumericMatrix DataState, IntegerVector choice_seq, 
						List DP_list, List control_list, bool draw_all){
	RNGScope scope;		// ensure RNG gets set/reset;
	// Initiate result storage;
	int nhh = hh_index.length(), N = choice_seq.length(), K = DP_list["K"];
	NumericMatrix ccp(N,K);
	NumericVector ccp_vec(N), c_seq(N);
	
	// Extract variables and parameters; 
	NumericMatrix state = DP_list["state"], Inc_weight = DP_list["Inc_weight"];
	NumericVector	V = DP_list["value_fn"], param = DP_list["param"], Q_grid = DP_list["Q_grid"], 
				 	TR_grid = DP_list["TR_grid"], price=DP_list["price"];
	double	beta = DP_list["beta"], price_min = min(price), Imax = max(state(_,0));
 	List y_kernal = DP_list["y_kernal"];
	Function omega_fn = DP_list["omega"];
	int n_Inc = Inc_weight.ncol(), ms = state.ncol();
	
	NumericVector 	ynodes = y_kernal["nodes"], yweight = y_kernal["weight"], 
					ysigma = y_kernal["sigma"], yrho = y_kernal["rho"];
	NumericMatrix ykappa = y_kernal["kappa"];
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
		NumericVector tmpy 	= V[index];
		myspl[i] = spl_init(x1, tmpy);
	}
	
	// Assign parameters for the value function; 
	struct my_f_params params = {param, beta, 
								 state(1,0), state(1,1), state(1,2), state(1,3),
								 1, Q_grid(1), TR_grid(1), 0, 
								 ynodes, yweight, price_min,
								 ykappa(_,0), ysigma(0), yrho(0), Inc_weight(0,_),
								 x1, x2, x3, x4, myspl};
																				
	// Optimizer setup;
	int status;
	int iter=0, max_iter = control_list["brent_max_iter"];
	double brent_tol	 = control_list["brent_tol"];
	const gsl_min_fminimizer_type *T;
	gsl_min_fminimizer *s;

	// GSL objective function ;
	gsl_function F; 
	F.function 	= &my_valuefn;
	F.params	= &params;
	
	//Create optimizer;
	T = gsl_min_fminimizer_brent;
	s = gsl_min_fminimizer_alloc(T);
	double low = 0, up = state(0,0);
	double m = 0.2, v= 0;
	int sat;
	
	// Loop over all the household; 
	int cnt = 0;				// Record the first position for an household; 
	for(int h=0; h<nhh; h++){
		// Assign initial values; 
		for(int j = 0; j < ms; j++){
			DataState(cnt,j) 	= init_state(h,j);
		}
		
		// Solve for q, c for every period in a loop;
		for(int i=cnt; i< cnt + TT_vec(h); i++){
			// Working objects;
			NumericVector vc(K), cstar(K);
		
			// Obtain ccp function;
			for(int k=0; k<K; k++){
				double tmp_omega 	= as<double>(omega_fn(k+1, DataState(i,1)));
				int Rc_int 			= DataState(i,3); 
				params.I_current	= DataState(i,0);
				params.lny			= DataState(i,1);
				params.Inc			= DataState(i,2);
				params.Rc			= DataState(i,3);
				params.k			= k + 1; 
				params.q			= Q_grid(k);
				params.TR			= TR_grid(k);
				params.omega		= tmp_omega;
				params.kappa		= ykappa(_, Rc_int);
				params.sigma		= ysigma(Rc_int);
				params.rho			= yrho(Rc_int);
				params.Inc_weight	= Inc_weight(n_Inc*Rc_int + DataState(i,2) -1 ,_);
	
				if(DataState(i,0)+Q_grid(k) <= 0){
					cstar(k)	= 0;
					vc(k)		= -my_valuefn(0, &params);
				}else{
				// Re-initiate minimizer and parameters;
					if(DataState(i,0) + params.q - Imax >0){
						low 	= DataState(i,0) + params.q - Imax;
					}else{
						low 	= 0; 
					}
					up			= DataState(i,0) + params.q;
					iter		= 0;
					
					// Try initial values close to the boundaries;
					m		= up - brent_tol;
					sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
					if(sat == GSL_EINVAL){
						m= low + brent_tol;
						sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
					}

					if(sat == GSL_EINVAL){
						cstar(k) =  up;
						vc(k) =  -my_valuefn(up, &params);	 
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
						cstar(k) = m;
						vc(k) = -v;
					}
				}
			}
			double vmax = max(vc);
			NumericVector tmpv = exp(vc - vmax); 
			ccp(i,_) = tmpv/sum(tmpv);		
			if(draw_all){
				choice_seq(i) 	= my_rmultinorm(ccp(i,_));
			}
			ccp_vec(i) 	= ccp(i,choice_seq[i] - 1);
			c_seq(i)	= cstar(choice_seq[i] - 1);
		
			// Evolve to the next state;
			if(i<= cnt + TT_vec(h)-2){
				DataState(i+1, 0) = DataState(i,0) + Q_grid(choice_seq[i] - 1) - cstar(choice_seq[i] - 1);
				if(draw_all){
					DataState(i+1,3) = DataState(i, 3); 	// Assuming recession is also constant; 
					int Rc_int 	= DataState(i+1, 3);
					DataState(i+1,2) = my_rmultinorm(Inc_weight(n_Inc*Rc_int + DataState(i,2) -1 ,_)); 	
					int s_int 	= DataState(i+1, 2);
					double mu	= ykappa(s_int - 1, Rc_int);	
					DataState(i+1, 1) =  rnorm(1, mu + yrho(Rc_int)*DataState(i,1), ysigma(Rc_int))[0] ;					
				}
			}
		}
		cnt			= cnt + TT_vec[h];
	}
	
	// Release memory of GSL pointers;
	gsl_min_fminimizer_free(s);
	for(int i=0; i<nspl; i++){
		myspl[i].splfree();
	}
	
	return Rcpp::List::create( _["ccp_vec"] = ccp_vec, _["ccp"] = ccp, _["c_seq"] = c_seq, 
							   _["DataState"] = DataState, _["choice_seq"] = choice_seq); 							
}


//-------------------------//;
// Log likelihood function //;
//-------------------------//;
// [[Rcpp::export]]
NumericVector ll_fnC(NumericVector param, NumericVector hh_index, IntegerVector TT_vec, NumericMatrix init_state,
					NumericMatrix DataState, IntegerVector choice_seq, 
					List DP_init, Function Bellman_operator, List control_list, IntegerVector policy_k_guess){
	// Construct DP_list and Bellman operator; 
	List DP_list = clone(DP_init);
	DP_list["param"] 	= param; 
// 	NumericMatrix state = DP_init["state"];
// 	int ns				= state.nrow();
// 	DP_list["value_fn"]	= V_initC(param, DP_init, policy_k_guess, state(_,0));
	int n = choice_seq.length();
	NumericVector ll(n, -100.0);
	
	if((param[0] <=0) | (param[2] <0) ){
		return(ll);
	}else{
		List sol = value_iteration_fnC(DP_list, Bellman_operator, 0, control_list);
		int status = sol["status"];
		if(status == 2){
			return(ll);
		}else{
			List sim_seq = simulate_hhseq_fnC(hh_index, TT_vec, init_state, DataState, 
							choice_seq, DP_list, control_list, FALSE);
			NumericVector ccp_vec = sim_seq["ccp_vec"];	
			ll = log(ccp_vec);
			return(ll);
		}
	}
}

