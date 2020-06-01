#include <RcppArmadillo.h>
#include <algorithm>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

// int I_idx=0, y_idx = 1;

// [[Rcpp::export]]
double singleu_fn(double c, double I, double Q, NumericVector param){
	double out=0, I_next = 0;
	double lambda 	= param[0];
	double tau1 	= param[1];
	double tau2	 	= param[2];
	double tau3 	= param[3];
	I_next = I + Q - c;
	out = lambda*log(c+.01) - tau1*(I+Q-c) -tau2*pow(I+Q-c,2) - tau3*Q;
	return out;
}

// [[Rcpp::export]]
colvec flow_Eu_fn(NumericMatrix state, List policy, NumericVector param, List Q_kernal, List y_kernal){
	NumericVector Qnodes 	= Q_kernal[0], ynodes = y_kernal[0];
	IntegerVector policy_k	= policy[0];
	NumericMatrix policy_c	= policy[1];
	List Qweight = Q_kernal[1];
	int ns = state.nrow(), nQ = Qnodes.size(), ;
	mat u(ns,nQ);
	int k = policy_k[0];
	mat w = as<mat>(Qweight[0]);
	int yidx = state(0,1);
	
	for(int i=0; i<ns; ++i){
		k = policy_k[i]-1;
		w = as<mat>(Qweight[k]);
		yidx = state(i,1)-1;
		for(int j=0; j<nQ; ++j){
			u(i,j) = singleu_fn(policy_c(i,j), state(i,0), Qnodes[j], param)*w(yidx,j);
		}
	}
	colvec Eu = sum(u,1);
	return (Eu);
}

// [[Rcpp::export]]
int findInterval1(double x, NumericVector breaks) {
	int out;
	out = std::distance(breaks.begin(), std::upper_bound(breaks.begin(), breaks.end(), x));
	return out;
}

// [[Rcpp::export]]
IntegerVector findState_Interval(double I, NumericVector I_grid, IntegerMatrix s_idx){
	int sel = std::max(1,findInterval1(I, I_grid));
	IntegerVector Ivec = s_idx(_,0);
	IntegerVector index = seq_len(Ivec.size());
	IntegerVector out = index[Ivec==sel];
	return(out);
}

// [[Rcpp::export]]
IntegerVector findState_Discrete(int y, IntegerMatrix s_idx){
	IntegerVector yvec = s_idx(_,1);
	IntegerVector index = seq_len(yvec.size());
	IntegerVector out = index[yvec==y];
	return(out);
}

// [[Rcpp::export]]
IntegerVector findState(NumericVector state_single, IntegerMatrix s_idx, NumericVector I_grid){
	IntegerVector sel_I = findState_Interval(state_single[0], I_grid, s_idx);
	IntegerVector sel_y = findState_Discrete(state_single[1], s_idx);
	IntegerVector out = intersect(sel_I, sel_y);
	return(out);
}

// [[Rcpp::export]]
NumericMatrix transition_fn(NumericMatrix state, IntegerMatrix s_idx, List policy, List Q_kernal, List y_kernal, NumericVector I_grid){
	NumericVector Qnodes 	= Q_kernal[0], ynodes = y_kernal[0];
	IntegerVector policy_k	= policy[0];
	NumericMatrix policy_c	= policy[1];
	List Qweight 			= Q_kernal[1];
	NumericMatrix yweight	= y_kernal[1];
	IntegerVector I_idxvec	= s_idx(_,0);
	int ns = state.nrow(), nQ = Qnodes.size(), ny = ynodes.size();
	
	double I_next = 0;
	int sel_I 	= 1;
	int sel_y	= 1;
	NumericMatrix out(ns,ns) ;	
	
	for(int i=0; i<ns; i++){
		int k = policy_k[i]-1;
		NumericMatrix w = Qweight[k];
		sel_y =  state(i,1)-1;
		
		// Transition probability of inventory;
		for(int j=0; j<nQ; j++){
			I_next 	= state(i,0) + Qnodes[j] - policy_c(i,j);
			sel_I	= std::max(1,findInterval1(I_next, I_grid));
			for(int l=0; l<ns; l++){
				if(I_idxvec[l]==sel_I) {out(i,l) += w(sel_y,j);}
			}			
		}
		
		// Transition probability of budget;
		for(int j=0; j<ny; j++){
			for(int l=0; l<ns; l++){
				if(state(l,1)==ynodes[j]){ out(i,l) *= yweight(sel_y,j); }
			}
		}
	}
	return(out);
}

// [[Rcpp::export]]
List solve_kc_fn(List DP_list, List Q_kernal, List y_kernal, NumericVector I_grid, NumericVector c_grid, arma::mat s_idx, IntegerVector plot_state){
	NumericMatrix state 	= DP_list["state"];
	NumericVector param		= DP_list["param"];
	NumericVector Qnodes	= Q_kernal[0], ynodes = y_kernal[0];
	mat yweight				= as<mat>(y_kernal[1]);
	vec I_idxvec			= s_idx.col(0);
	List Qweight 			= Q_kernal[1];
	vec W 					= as<vec>(DP_list["value_fn"]);	
	double beta 			= DP_list["beta"];
	int ns = state.nrow(), nQ = Qnodes.size(), ny = ynodes.size(), nc=c_grid.size(), K = Qweight.size();
	
	// Initiate the return policy;
	vec policy_k	= zeros(ns);
	mat policy_c	= zeros(ns,nQ);
	vec value_max	= zeros(ns);
	vec vc 			= zeros(nc);
	mat V			= zeros(ns,nQ);
	mat V_choice	= zeros(ns,K);
	
	int ploti = 0;
	List plot_data(plot_state.size());
	for(int i=0; i<ns; i++){
		int sel_y =  state(i,1)-1;
		mat tmp_plot = zeros(nc, nQ);
		for(int j=0; j<nQ; j++){			
			for(int l=0; l<nc; l++){
				double I_next = state(i,0) + Qnodes[j] - c_grid[l];
				int sel_I = findInterval1(I_next, I_grid);
				if(sel_I<1){
					vc.row(l) = NA_REAL;
				}else{
					vc.row(l) = singleu_fn(c_grid[l], state(i,0), Qnodes[j], param) +
								beta* (yweight.row(sel_y) * W.elem(arma::find(I_idxvec==sel_I)));
				}
			}
			
			uword tmp;
			double tmp1 = vc.max(tmp);
			if(tmp<=nc){
				V(i,j) 			= tmp1;
				policy_c(i,j)	= c_grid[tmp];
			}else{
				V(i,j) 			= NA_REAL;
				policy_c(i,j)	= c_grid[0];
			}
			
			if(i==plot_state[ploti]){
				tmp_plot.col(j) = vc;
			}
		}		
		for(int k=0; k<K; k++){
			mat Qw = as<mat>(Qweight[k]);
			V_choice(i,k) = sum(V.row(i)*trans(Qw.row(sel_y)));
		}
		uword tmp;
		double tmp1 = V_choice.row(i).max(tmp);
		if(tmp<=K){
			policy_k.row(i)	= tmp+1;
			value_max.row(i) = tmp1;
		}else{
			policy_k.row(i)	= K;
			value_max.row(i) = NA_REAL;
		}
		
		if(i==plot_state[ploti]){
			plot_data[ploti] = tmp_plot;
			ploti++;
		}
	}

	return Rcpp::List::create(
		Rcpp::Named("policy") = Rcpp::List::create(Rcpp::Named("k")=policy_k, Rcpp::Named("c")=policy_c),
		Rcpp::Named("Value") = value_max,
		Rcpp::Named("plot_data") = plot_data
	);
}

// [[Rcpp::export]]
NumericVector CCP_single_fn(NumericVector state_single, List DP_list, NumericVector I_grid, List Q_kernal, List y_kernal){
	NumericMatrix state 	= DP_list["state"];
	NumericVector param		= DP_list["param"];
	IntegerMatrix s_idx		= DP_list["state_idx"];
	NumericVector Qnodes	= Q_kernal[0], ynodes = y_kernal[0];
	mat yweight				= as<mat>(y_kernal[1]);
	List Qweight 			= Q_kernal[1];
	vec W 					= as<vec>(DP_list["value_fn"]);	
	List policy				= DP_list["policy"];
	NumericMatrix policy_c	= policy["c"];
	double beta 			= DP_list["beta"];
	int nQ = Qnodes.size(), K = Qweight.size();
	
	vec V		= zeros(nQ);
	vec EV		= zeros(K);
	
	int sel 	= findState(state_single, s_idx, I_grid)[0]-1;
	int sel_y 	= state_single[1]-1;
	for(int i=0; i<nQ; i++){
		double I_next = state_single[0] + Qnodes[i] - policy_c(sel,i);
		uvec idx = as<uvec>(findState_Interval(I_next, I_grid, s_idx));
		idx = idx - 1;
		V.row(i) = singleu_fn(policy_c(sel,i), state_single[0], Qnodes[i], param) +
					beta* (yweight.row(sel_y) * W.elem(idx));
	}
	for(int i=0; i<K; i++){
		mat Qw = Qweight[i];
		EV.row(i) = Qw.row(sel_y) * V;
	}
	vec EV_exp 	= exp(EV - EV.max());
	vec ccp 	= EV_exp/sum(EV_exp);
	NumericVector out = NumericVector(ccp.begin(),ccp.end());
	return(out);
}

// [[Rcpp::export]]
mat CCP_fullstate_fn(List DP_list, NumericVector I_grid, List Q_kernal, List y_kernal){
	NumericMatrix state 	= DP_list["state"];
	List Qweight			= Q_kernal[1];
	int ns = state.nrow(), K = Qweight.size();
	mat out(ns,K);
	for(int i=0; i<ns; i++){ 	
		out.row(i) 	= as<rowvec>(CCP_single_fn(state(i,_), DP_list, I_grid, Q_kernal, y_kernal));
	}
	return(out);
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
NumericVector ll_fullstate_fn(NumericVector theta, NumericMatrix state, IntegerVector choice_seq, List DP_init, 
							Function policy_iteration_fn, Function u_fn, Function transition_fn, Function policy_improve_fn, 
							List set_list,List Q_kernal, List y_kernal, NumericVector I_grid, NumericVector c_grid){
	NumericVector param(theta.size());
	List est_list 		= DP_init;
	param[0]			= exp(theta[0]);
	param[1]			= exp(theta[1]);
	param[2]			= -exp(theta[2]);
	param[3]			= exp(theta[3]);
	est_list["param"] 	= param;
	int ns 				= state.nrow();
	IntegerMatrix s_idx	= DP_init["state_idx"];
	List sol 			= policy_iteration_fn(est_list,set_list,u_fn,transition_fn,policy_improve_fn,
							Q_kernal,y_kernal,I_grid,c_grid);
	List sol_DP			= sol[0];
	mat ccp_mat 		= CCP_fullstate_fn(sol_DP, I_grid, Q_kernal, y_kernal);
	NumericVector out(ns);

	for(int i=0; i<ns; i++){
		int sel = findState(state(i,_), s_idx, I_grid)[0] - 1;
		out(i)	= log(ccp_mat.at(sel,choice_seq[i]-1));		
	}
	return(out);	
}





