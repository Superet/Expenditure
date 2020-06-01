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
int findInterval1(double x, NumericVector breaks) {
	int out;
	out = std::distance(breaks.begin(), std::upper_bound(breaks.begin(), breaks.end(), x));
	return out;
}

struct spl4d{
	NumericMatrix x; NumericVector y;
	
	double spl4dEval(NumericVector xnew){
		NumericVector x1 = unique(x(_,0)), x2 = unique(x(_,1)), x3 = unique(x(_,2)), x4 = unique(x(_,3));
		std::sort(x1.begin(), x1.end()); 
		std::sort(x2.begin(), x2.end()); 
		std::sort(x3.begin(), x3.end()); 
		std::sort(x4.begin(), x4.end()); 
		int N1 = x1.length(), N2 = x2.length(), N3 = x3.length(), N4 = x4.length();
		int x3_index = findInterval1(xnew(2), x3);		//Find which discrete value it is in the 3rd dimension; 
		int x4_index = findInterval1(xnew(3), x4);		//Find which discrete value it is in the 3rd dimension; 
		
		// Interpolating the 1st value conditional on the rest; 
		NumericVector y2(N2);
		for(int i=0; i<N2; i++){
			IntegerVector index = (x4_index -1)*N1*N2*N3 + (x3_index -1)*N1*N2 + i*N1 + seq_len(N1) - 1;
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

struct my_f_params { double lambda; double tau1; double tau2;double beta; 
					double I_current; double lny; double Inc; double Rc;
					double q; double omega;		
					NumericVector ynodes; NumericVector yweight;
					NumericVector kappa; double sigma; double rho; 
					NumericMatrix state; NumericVector V;};

//---------------------------//;
// Value function definition //;
//---------------------------//;
double my_valuefn(double cs, void *p){
	struct my_f_params *params = (struct my_f_params *)p;	
    double lambda				= (params->lambda);
    double tau1					= (params->tau1);
    double tau2					= (params->tau2);
	double beta					= (params->beta);
	double I_current			= (params->I_current);
	double lny					= (params->lny);
	double Inc					= (params->Inc);
	double Rc					= (params->Rc);
	double q 					= (params->q);
	double omega				= (params->omega);
	NumericVector yweight		= (params->yweight);
	NumericVector ynodes		= (params->ynodes);
	NumericVector kappa			= (params->kappa);
	double sigma				= (params->sigma);
	double rho					= (params->rho);
	NumericMatrix state			= (params->state);
	NumericVector V				= (params->V);
	
	int  n 			= ynodes.length();
	NumericVector lny_next(n), v(n);
	double I_next 	= I_current + q - cs, ev = 0;
	
	// Assign values to the next state and ; 
	NumericVector  s_next(4);
	s_next(0) 	= I_next; 
	s_next(1) 	= lny_next(0);
	s_next(2) 	= Inc; 
	s_next(3) 	= Rc; 
	int s_int 	= s_next(2);
	double mu	= kappa(0);
	if(s_int>1){
		mu 	= mu + kappa(s_int-1);
	}
	lny_next 		= sqrt(2)*sigma*ynodes + mu + rho*lny;
	
	// Interpolate the value function for the next state; 
	gsl_set_error_handler_off();
	struct spl4d myspl = {state, V};
	for(int i=0; i<n; i++){
		s_next(1) 	= lny_next(i); 
		v(i) 		= myspl.spl4dEval(s_next);
	}
	
	// Compute the value function; 
	ev = sum(yweight*v)/sqrt(PI);
	double out = lambda*log(cs+1) - tau1*I_next - tau2*q + omega + beta*ev;
	return(-out);
}

//------------------//;
// Bellman operator //;
//------------------//;
// [[Rcpp::export]]
List Solve_DP_fnC(NumericMatrix state, NumericVector V_old, NumericMatrix omega_mat, int K,
				  NumericVector param, double beta, NumericVector Q_grid, List y_kernal, NumericVector price, List control_list){
	const double lambda = param(0), tau1 = param(1), tau2 = param(2);
	NumericVector 	ynodes = y_kernal["nodes"], yweight = y_kernal["weight"], 
					ysigma = y_kernal["sigma"], yrho = y_kernal["rho"];
	NumericMatrix ykappa = y_kernal["kappa"];
	int ns = state.nrow();
	gsl_set_error_handler_off();

	// Assign parameters for the value function; 
	struct my_f_params params = {lambda, tau1, tau2, beta, 
								 state(1,0), state(1,1), state(1,2), state(1,3),
								 Q_grid(1), omega_mat(0,1),
								 ynodes, yweight, ykappa(_,0), ysigma(0), yrho(0),state, V_old};
										
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
	
	double checkv = my_valuefn(0, &params);
					
	// Solve for c for every state in a loop;
	for(int i=0; i<ns; i++){
		for(int k=0; k<K; k++){
			int Rc_int 			= state(i,3); 
			params.I_current	= state(i,0);
			params.lny			= state(i,1);
			params.Inc			= state(i,2);
			params.Rc			= state(i,3);
			params.q			= Q_grid(k);
			params.omega		= omega_mat(i,k);
			params.kappa		= ykappa(_, Rc_int);
			params.sigma		= ysigma(Rc_int);
			params.rho			= yrho(Rc_int);
			
			vc(i,k) = my_valuefn(0, &params);
			if(state(i,0)+Q_grid[k] <= 0){
				cstar(i,k)	= 0;
				vc(i,k)		= -my_valuefn(0, &params);
			}else if(Q_grid(k) > exp(state(i,1))/min(price)){
				cstar(i,k) 	= NA_REAL;
				vc(i,k)		= -500;				
			}else{
				// Re-initiate minimizer and parameters;
				low			= 0;
				up			= state(i,0) + params.q;
				iter		= 0;
				// Try initial values close to the boundaries;
				m		= up - brent_tol;
				sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
				if(sat == GSL_EINVAL){
					m= brent_tol;
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

	return Rcpp::List::create(_["policy"] = List::create(_["c"] = cstar), 
							  _["value"]=V_new, _["ccp"] = ccp, _["vc"] = vc);
}


// [[Rcpp::export]]
List Bellman_operatorC(List DP_list, List control_list, bool is_data=false, NumericMatrix DataState = NumericMatrix(2,2) ){
	NumericMatrix state = DP_list["state"];
	NumericVector V_old = DP_list["value_fn"], param = DP_list["param"], Q_grid = DP_list["Q_grid"], price=DP_list["price"];
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
		out  = Solve_DP_fnC(DataState, V_old, omega, K, param, beta, Q_grid, y_kernal, price, control_list); 
	}else{
		int ns = state.nrow();
		NumericMatrix omega(ns, K);
		for(int k=0; k<K; k++){
			NumericVector tmp = omega_fn(k+1, state(_,1));
			omega(_,k) = tmp;
		}
		out  = Solve_DP_fnC(state, V_old, omega, K, param, beta, Q_grid, y_kernal, price, control_list); 
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
		Rcout<<"Value function iteration has converged.\n";
	}else{
		DP_listout["status"] = 2;	
		Rcout<<"Value function iteration did not converge.\n";
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
						NumericVector Q_seq, IntegerVector k_seq, 
						List DP_list, List control_list){
	RNGScope scope;		// ensure RNG gets set/reset
	NumericMatrix state = DP_list["state"];
	NumericVector V = DP_list["value_fn"], param = DP_list["param"], Q_grid = DP_list["Q_grid"], price=DP_list["price"];
	int K = DP_list["K"];
	double	beta = DP_list["beta"];
 	List y_kernal = DP_list["y_kernal"];
	Function omega_fn = DP_list["omega"];
	int ns = state.nrow(), ms = state.ncol();
	
	const double lambda = param(0), tau1 = param(1), tau2 = param(2);
	NumericVector 	ynodes = y_kernal["nodes"], yweight = y_kernal["weight"], 
					ysigma = y_kernal["sigma"], yrho = y_kernal["rho"];
	NumericMatrix ykappa = y_kernal["kappa"];
	gsl_set_error_handler_off();

	// Assign parameters for the value function; 
	struct my_f_params params = {lambda, tau1, tau2, beta, 
								 state(0,0), state(0,1), state(0,2), state(0,3), 
								 Q_grid(0), 0,
								 ynodes, yweight, ykappa(_,0), ysigma(0), yrho(0), state, V};
																				
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
	NumericVector c_seq(TT), ccp_vec(TT);
	
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
			params.q			= Q_grid(k);
			params.omega		= tmp_omega;
			params.kappa		= ykappa(_, Rc_int);
			params.sigma		= ysigma(Rc_int);
			params.rho			= yrho(Rc_int);
	
			if(DataState(i,0)+Q_grid(k) < 0){
				cstar(k)	= 0;
				vc(k)		= -my_valuefn(0, &params);
			}else if(Q_grid(k) > DataState(i,1)/min(price)){
				cstar(k) 	= NA_REAL;
				vc(k)		= -500;				
			}else{
				// Re-initiate minimizer and parameters;
				low			= 0;
				up			= DataState(i,0) + params.q;
				iter		= 0;
				// Try initial values close to the boundaries;
				m		= up - brent_tol;
				sat = gsl_min_fminimizer_set ( s, &F, m, low, up );
				if(sat == GSL_EINVAL){
					m= brent_tol;
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
			DataState(i+1, 0) = DataState(i,0) + Q_seq(i) - c_seq(i);
			if(draw_all){
				DataState(i+1,2) = DataState(i, 2); 	// Assuming income level does not change; 
				DataState(i+1,3) = DataState(i, 3); 	// Assuming recession is also constant; 
				int Rc_int 	= DataState(i+1, 3);
				int s_int 	= DataState(i+1, 2);
				double mu	= ykappa(0, Rc_int);
				if(s_int>1){
					mu 	= mu + ykappa(s_int-2, Rc_int);
				}
							
				DataState(i+1, 1) =  rnorm(1, mu + yrho(Rc_int)*DataState(i,1), ysigma(Rc_int))[0] ;					
			}
		}
	}
	// Release memory of GSL pointers;
	gsl_min_fminimizer_free(s);
	
	return Rcpp::List::create(Rcpp::DataFrame::create(_["t"] = seq_len(TT), _["I"]=DataState(_,0),
								_["y"] = DataState(_,1),_["k"] = k_seq, _["Q"]=Q_seq, _["c"]=c_seq), 
		    				  	_["ccp"]=ccp, _["ccp_vec"] = ccp_vec);
}

//-------------------------//;
// Log likelihood function //;
//-------------------------//;