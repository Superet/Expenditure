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

struct my_f_params { double lambda; double tau1; double tau2;double beta; 
					double I_current; double y; double q; double omega;		
					NumericVector y_grid; NumericVector ynodes; NumericVector yweight;
					double mu; double sigma; double rho; struct spl *spl2d;};
					

// [[Rcpp::export]]					
double uP_fnC(NumericVector e, NumericVector psi, NumericVector gamma, NumericVector price, int S, double qz_cons){
	NumericVector u(S+2);
	for(int i=0; i<S; i++){
		u[i] = log(e[i]/(gamma[i]*price[i]) + 1);
	}
	u[S] 	= log(e[S]/gamma[S]);
	u[S+1]	= log(e[S+1]/gamma[S+1] + qz_cons);
	double out = sum(psi*gamma*u);
	return(out);
}					

double my_valuefn(double cs, void *p){
	struct my_f_params *params = (struct my_f_params *)p;	
    double lambda				= (params->lambda);
    double tau1					= (params->tau1);
    double tau2					= (params->tau2);
	double beta					= (params->beta);
	double q					= (params->q);
	double I_current			= (params->I_current);
	double y					= (params->y);
	double omega				= (params->omega);
	NumericVector y_grid		= (params->y_grid);
	NumericVector yweight		= (params->yweight);
	NumericVector ynodes		= (params->ynodes);
	double mu					= (params->mu);
	double sigma				= (params->sigma);
	double rho					= (params->rho);
	struct spl *spl2d			= (params->spl2d);
	
	int  n = ynodes.length(), ny = y_grid.length();
	double I_next = I_current + q - cs, ev = 0, Iq = 0;
	NumericVector y_next(n), v1(ny), v2(n);
	y_next = exp(sqrt(2)*sigma*ynodes + mu + rho*log(y));
	if(q>0){Iq = 1; }
	
	// Interpolate I_next along conditional on y;
	for(int i=0; i<ny; i++){
		v1[i] = spl2d[i].splEval(I_next);
	}
	
	// Create the spline in the second dimension;
	struct spl myspl;
	double eps 		= 0.0001;
	myspl.acc 		= gsl_interp_accel_alloc();    	
	myspl.spline	= gsl_spline_alloc( gsl_interp_cspline , ny );
	myspl.x			= y_grid;
	myspl.y			= v1;
	gsl_spline_init( myspl.spline, y_grid.begin(), v1.begin(), ny);	
	myspl.slope_first= (gsl_spline_eval (myspl.spline, y_grid[0]+eps, myspl.acc) - v1[0])/eps;
	myspl.slope_end	 = (v1[ny-1] - gsl_spline_eval (myspl.spline, y_grid[ny-1]-eps, myspl.acc))/eps;
	
	// Interopolate the value function at quadrature nodeds;
	for(int i=0; i<n; i++){
		v2[i] = myspl.splEval(y_next[i]);
	}
	ev = sum(yweight*v2)/sqrt(PI);
	double out = lambda*log(cs+0.01) - tau1*I_next - tau2*Iq + omega + beta*ev;
	myspl.splfree();
	return(-out);
}

// [[Rcpp::export]]
List Solve_DP_fnC(NumericMatrix state, NumericVector I_grid, NumericVector y_grid, NumericMatrix V_old, NumericMatrix omega_mat, int K,
				  NumericVector param, double beta, NumericVector Q_grid, List y_kernal, NumericVector price, List control_list){
	const double lambda = param(0), tau1 = param(1), tau2 = param(2);
	NumericVector ynodes  = y_kernal["nodes"], yweight = y_kernal["weight"];
	double ymu = y_kernal["mu"], ysigma = y_kernal["sigma"], yrho = y_kernal["rho"];
	int ns = state.nrow(), nx = I_grid.length(), ny = y_grid.length();
	
	// Spline interpolation set up;
	gsl_set_error_handler_off();
	struct spl myspl[ny];
	double eps = 0.0001;
	for(int i=0; i<ny; i++){
		NumericVector tmpz = V_old(_,i);
		myspl[i].acc 	= gsl_interp_accel_alloc();    	
		myspl[i].spline = gsl_spline_alloc( gsl_interp_cspline , nx );
		myspl[i].x		= I_grid;
		myspl[i].y		= tmpz;
		gsl_spline_init( myspl[i].spline, I_grid.begin(), tmpz.begin(), nx);	
		myspl[i].slope_first = (gsl_spline_eval (myspl[i].spline, I_grid[0]+eps, myspl[i].acc) - tmpz[0])/eps;
		myspl[i].slope_end	 = (tmpz[nx-1] - gsl_spline_eval (myspl[i].spline, I_grid[nx-1]-eps, myspl[i].acc))/eps;
	}
	
	struct my_f_params params = {lambda, tau1, tau2, beta, state(0,0), state(0,1), Q_grid(0), omega_mat(0,0),
								 y_grid, ynodes, yweight, ymu, ysigma, yrho, myspl};
										
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
	double low = 0, up = I_grid[nx-1];
	double m = 0.2, v= 0;
	int sat;
	
	// The object to store computation values; 
	NumericMatrix cstar(ns, K), vc(ns, K), ccp(ns, K);
	NumericVector V_new(ns), policy_k(ns);

	// Solve for c for every state in a loop;
	for(int i=0; i<ns; i++){
		for(int k=0; k<K; k++){
			params.I_current	= state(i,0);
			params.y			= state(i,1);
			params.q			= Q_grid(k);
			params.omega		= omega_mat(i,k);
			if(state(i,0)+Q_grid[k] <= 0){
				cstar(i,k) = 0;
				vc(i,k)		= -my_valuefn(0, &params);
			}else if(Q_grid(k) > exp(state(i,1))/min(price)){
				cstar(i,k) = NA_REAL;
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
					cstar(i,k) =  up;
					vc(i,k) =  -my_valuefn(up, &params);	 
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
		for(int k = 0; k<K; k++){
			if(vc(i,k)==vmax) { policy_k(i) = k+1;}
		}
		ccp(i,_) = tmpv/sum(tmpv);
	}
	
	// Release memory of GSL pointers;
	gsl_min_fminimizer_free(s);
	for(int i=0; i<ny; i++){
		myspl[i].splfree();
	}

	return Rcpp::List::create(_["policy"] = List::create(_["k"] = policy_k, _["c"] = cstar), 
							  _["value"]=V_new, _["ccp"] = ccp, _["vc"] = vc);
}

// [[Rcpp::export]]
List Bellman_operatorC(List DP_list, List control_list, bool is_data=false, NumericMatrix DataState = NumericMatrix(2,2) ){
	NumericMatrix state = DP_list["state"];
	NumericVector V = DP_list["value_fn"], param = DP_list["param"], Q_grid = DP_list["Q_grid"], price=DP_list["price"], 
				  I_grid = unique(state(_,0)), y_grid = unique(state(_,1));
	int K = DP_list["K"];
	double	beta = DP_list["beta"];
 	List y_kernal = DP_list["y_kernal"];
	Function omega_fn = DP_list["omega"];
	std::sort(I_grid.begin(), I_grid.end()); 
	std::sort(y_grid.begin(), y_grid.end()); 
	int nx = I_grid.length(), ny = y_grid.length(), ns;
	NumericMatrix V_old(nx, ny, V.begin());
	
	List out;
	if(is_data){
		ns = DataState.nrow();
		NumericMatrix omega(ns, K);
		for(int k=0; k<K; k++){
			NumericVector tmp = omega_fn(k+1, DataState(_,1));
			omega(_,k) = tmp;
		}
		out  = Solve_DP_fnC(DataState, I_grid, y_grid, V_old, omega, K, param, beta, Q_grid, y_kernal, price, control_list);			
	}else{
		ns = state.nrow();
		NumericMatrix omega(ns, K);
		for(int k=0; k<K; k++){
			NumericVector tmp = omega_fn(k+1, state(_,1));
			omega(_,k) = tmp;
		}
		out  = Solve_DP_fnC(state, I_grid, y_grid, V_old, omega, K, param, beta, Q_grid, y_kernal, price, control_list);
	}
	return(out);
}

// [[Rcpp::export]]
List value_iteration_fnC(List DP_list, Function Bellman_operator,  int print_level, List control_list){	
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

int my_rmultinorm(NumericVector probs) {
    int k = probs.size();
    IntegerVector ans(k);
    rmultinom(1, probs.begin(), k, ans.begin());
    IntegerVector index = seq_len(k);
    IntegerVector out = index[ans==1];
    return(out[0]);
}

// [[Rcpp::export]]
List simulate_seq_fnC( NumericVector init_state, int TT, bool draw_all, 
							NumericVector y_seq, NumericVector Q_seq, IntegerVector k_seq, 
							List DP_list, List control_list){
	RNGScope scope;		// ensure RNG gets set/reset
	NumericMatrix state = DP_list["state"];
	NumericVector V = DP_list["value_fn"], param = DP_list["param"], Q_grid = DP_list["Q_grid"], price=DP_list["price"], 
				  I_grid = unique(state(_,0)), y_grid = unique(state(_,1));
	int K = DP_list["K"];
	double	beta = DP_list["beta"];
 	List y_kernal = DP_list["y_kernal"];
	Function omega_fn = DP_list["omega"];
	std::sort(I_grid.begin(), I_grid.end()); 
	std::sort(y_grid.begin(), y_grid.end()); 
	int nx = I_grid.length(), ny = y_grid.length(), ns = state.nrow();
	NumericMatrix V_old(nx, ny, V.begin());
	
	const double lambda = param(0), tau1 = param(1), tau2 = param(2);
	NumericVector ynodes  = y_kernal["nodes"], yweight = y_kernal["weight"];
	double ymu = y_kernal["mu"], ysigma = y_kernal["sigma"], yrho = y_kernal["rho"];
	
	// Spline interpolation set up;
	gsl_set_error_handler_off();
	struct spl myspl[ny];
	double eps = 0.0001;
	for(int i=0; i<ny; i++){
		NumericVector tmpz = V_old(_,i);
		myspl[i].acc 	= gsl_interp_accel_alloc();    	
		myspl[i].spline = gsl_spline_alloc( gsl_interp_cspline , nx );
		myspl[i].x		= I_grid;
		myspl[i].y		= tmpz;
		gsl_spline_init( myspl[i].spline, I_grid.begin(), tmpz.begin(), nx);	
		myspl[i].slope_first = (gsl_spline_eval (myspl[i].spline, I_grid[0]+eps, myspl[i].acc) - tmpz[0])/eps;
		myspl[i].slope_end	 = (tmpz[nx-1] - gsl_spline_eval (myspl[i].spline, I_grid[nx-1]-eps, myspl[i].acc))/eps;
	}
	
	struct my_f_params params = {lambda, tau1, tau2, beta, state(0,0), state(0,1), Q_grid(0), 0,
								 y_grid, ynodes, yweight, ymu, ysigma, yrho, myspl};
										
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
	double low = 0, up = I_grid[nx-1];
	double m = 0.2, v= 0;
	int sat;
	
	// Initiate result storage;
	NumericMatrix DataState(TT, state.ncol()), ccp(TT, K);
	NumericVector c_seq(TT), ccp_vec(TT);
	
	// Solve for q, c for every period in a loop;
	DataState(0,0) = init_state(0);
 	DataState(0,1) = init_state(1);
	for(int i=0; i<TT; i++){
		// Working objects;
		NumericVector vc(K), cstar(K);
		
		// Obtain ccp function;
		for(int k=0; k<K; k++){
			double tmp_omega 	= as<double>(omega_fn(k+1, DataState(i,1)));
			params.I_current	= DataState(i,0);
			params.y			= DataState(i,1);
			params.q			= Q_grid(k);
			params.omega		= tmp_omega;
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
				DataState(i+1, 1) =  exp(rnorm(1, ymu + yrho*log(DataState(i,1)), ysigma)[0]);	
			}else{
				DataState(i+1, 1) = y_seq(i+1);	
			}
		}
	}
	// Release memory of GSL pointers;
	gsl_min_fminimizer_free(s);
	for(int i=0; i<ny; i++){
		myspl[i].splfree();
	}
	
	return Rcpp::List::create(Rcpp::DataFrame::create(_["t"] = seq_len(TT), _["I"]=DataState(_,0),
								 _["y"] = DataState(_,1),_["k"] = k_seq, _["Q"]=Q_seq, _["c"]=c_seq), 
		    				  _["ccp"]=ccp, _["ccp_vec"] = ccp_vec);
}

// [[Rcpp::export]]
NumericVector ll_fnC(NumericVector param, NumericMatrix DataState, IntegerVector choice_seq, 
					NumericVector Q_seq, List DP_init, Function Bellman_operator, List control_list){
	// Construct DP_list and Bellman operator; 
	List DP_list = clone(DP_init);
	DP_list["param"] = param; 
	NumericVector init_state = DataState(0,_), y_seq = DataState(_,1);
	int n = choice_seq.length();
	NumericVector ll(n, -100.0);
	
	if(param[0] <=0 || param[1] <=0 ){
		return(ll);
	}else{
		List sol = value_iteration_fnC(DP_list, Bellman_operator, 0, control_list);
		int status = sol["status"];
		if(status == 2){
			return(ll);
		}else{
			List sim_seq = simulate_seq_fnC(init_state, n, false, y_seq, Q_seq, choice_seq, sol, control_list);
			NumericVector ccp_vec = sim_seq["ccp_vec"];	
			ll = log(ccp_vec);
			return(ll);
		}
	}
}
