#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]

struct constr_par{const gsl_matrix *ui; const gsl_vector *ci; double *mu; void *p ; };

struct constr_multimin{
	gsl_multimin_fdfminimizer *s; gsl_multimin_function_fdf F; gsl_multimin_function_fdf orig_F; 
	double outer_eps; int outer_iter; double mu_max;  
	
	double inner_multimin(gsl_vector *x, gsl_vector *par){
		int iter = 0, status, n=F.n; 
		NumericVector out(n+1);
	
		gsl_multimin_fdfminimizer_set (s, &F, x, 0.001, 0.01);
		iter = 0; 
		status = 0; 
		status = gsl_multimin_fdfminimizer_iterate (s);

		do{
			iter++;
			status = gsl_multimin_fdfminimizer_iterate (s);
			if (status){	
				Rcout<<"GSL_error, f="<< gsl_multimin_fdfminimizer_minimum(s)<<"\n"; 
				break; 
			}
			status = gsl_multimin_test_gradient (s->gradient, 1e-6);
		}while (status == GSL_CONTINUE && iter < 500);
		
		// Return the minimizer and optimal value; 
		Rcout<<"gradient="<<gsl_vector_get(s->gradient, 0)<<","<<gsl_vector_get(s->gradient, 1)<<
			", "<<gsl_vector_get(s->gradient, 2)<<"\n";
		gsl_vector_memcpy(par, s->x); 
		return(s->f );
	}
	
	NumericVector constr_optim(gsl_vector *x0){
		double obj, r, obj_old, r_old; 
		int n = F.n; 
		NumericVector out(n+1); 
		gsl_vector *x; 
		x = gsl_vector_alloc(F.n); 
		gsl_vector_memcpy (x, x0);
		struct constr_par *par	= (struct constr_par *)F.params;
		double *mu_ptr = (par->mu); 
		
		obj			= orig_F.f(x0, orig_F.params); 
		r			= F.f(x0, F.params); 
		int check 	= 0; 
		int iter = 0; 
		Rcout<<"Initial "<<iter<<": r="<<r<<", obj="<<obj<<".\n";
		do{
			iter++; 
			gsl_vector_memcpy (x0, x);
			obj_old 	= obj; 
			r_old		= r; 
// 			constr_parl.x_old = x0; 
			r = inner_multimin(x0, x); 
			if( gsl_finite(r) && gsl_finite(r_old) && abs(r - r_old)<(0.001 + abs(r))*outer_eps && *mu_ptr <= mu_max){
				check = 1; 
			}
			obj = orig_F.f(x, orig_F.params); 
			*(mu_ptr) = *(mu_ptr)/5; 
// 			Rcout<<"Outer iteration"<<iter<<": r="<<r<<", r_old="<<r_old<<", obj="<<obj<<", obj_old"<<obj_old<<".\n"; 
// 			if(obj > obj_old) break; 
		}while(check==0 && iter < outer_iter);
		
		// Return the minimizer and optimal value; 
		for(int j=0; j<n; j++){
			out(j) = gsl_vector_get(x,j);
		}
		out(n) = obj;
		gsl_vector_free(x); 
		return(out); 
	}
};

// Define objective function and gradients; 
double FC(const gsl_vector *x, void *params){	
	NumericVector b(3); 
	double *p = (double *)params;
	b(0) = gsl_vector_get(x,0); 
	b(1) = gsl_vector_get(x,1); 
	b(2) = gsl_vector_get(x,2); 
	return( - p[0]*b[0] - p[1]*b[1] - p[2]*b[2] + p[3]*sum(pow(b,2))); 
}

void gradC(const gsl_vector *x, void *params, gsl_vector *df){
	NumericVector b(3); 
	double *p = (double *)params;
	b(0) = gsl_vector_get(x,0); 
	b(1) = gsl_vector_get(x,1); 
	b(2) = gsl_vector_get(x,2); 
	gsl_vector_set(df, 0, -p[0] + 2*p[3]*b(0) ); 
	gsl_vector_set(df, 1, -p[1] + 2*p[3]*b(1) ); 
	gsl_vector_set(df, 2, -p[2] + 2*p[3]*b(2) ); 
}

// Penalty function; 
double pnty(const gsl_vector *pgi){
	double out = 0; 
	for(size_t i=0; i<pgi->size; i++){
		// Log penalty; 
// 		out = out - log(gsl_vector_get(pgi, i)); 
		// quadratic function; 
		out = out + 1/gsl_vector_get(pgi, i); 
	}
	return(out);
}

void pnty_g(const gsl_vector *pgi, gsl_vector *dpnty){
	for(size_t i=0; i<pgi->size; i++){
// 		gsl_vector_set(dpnty, i, -1/gsl_vector_get(pgi,i)); 
		gsl_vector_set(dpnty, i, -1/pow(gsl_vector_get(pgi,i), 2)); 
	}
}

// Define Lagrangian function and gradients; 
double RC(const gsl_vector *x, void *params){
	struct constr_par *par	= (struct constr_par *)params;
	const gsl_matrix *ui	= (par->ui);
	const gsl_vector *ci	= (par->ci);
	void *p1				= (par->p);
	double *mu				= (par->mu);
	double mu_n				= *mu; 
	double out=0; 
	
	size_t n 	= ui->size2, m = ui->size1; 
	gsl_vector  *gi;
	gi			= gsl_vector_alloc(m); 
	double bar = 0 ; 
	
	gsl_vector_memcpy(gi, ci); 
	gsl_blas_dgemv(CblasNoTrans, 1.0, ui, x, -1.0, gi); 
	
	bar = pnty(gi); 
	out = FC(x, p1) + mu_n*bar; 
	gsl_vector_free(gi); 
	return(out);
}

void dRC(const gsl_vector *x, void *params, gsl_vector *df){
	struct constr_par *par	= (struct constr_par *)params;
	const gsl_matrix *ui	= (par->ui);
	const gsl_vector *ci	= (par->ci);
	void *p1				= (par->p);	
	double *mu				= (par->mu);
	double mu_n 			= *mu; 
	
	size_t n 	= ui->size2, m = ui->size1; 
	gsl_vector  *gi, *dbar;
	gi			= gsl_vector_alloc(m); 
	dbar		= gsl_vector_alloc(m); 
	gsl_vector_memcpy(gi, ci); 
	gsl_blas_dgemv(CblasNoTrans, 1.0, ui, x, -1.0, gi); 
	
	gradC(x, p1, df);
	pnty_g(gi, dbar); 
	gsl_blas_dgemv(CblasTrans, mu_n, ui, dbar, 1.0, df); 
		
	gsl_vector_free(dbar); 
	gsl_vector_free(gi); 
}

void RdRC(const gsl_vector *x, void *params, double *f, gsl_vector *df){
	*f = RC(x, params); 
	dRC(x, params, df);
}

// A wrapper of constrained optimization; 
// [[Rcpp::export]]
double my_fn(NumericVector init, NumericMatrix ui, NumericVector ci){
	double p[4] = {0, 5, 0, .5};
	int n = init.length(), m=ci.length(); 
	size_t nt = init.length(), mt = ci.length();  
	
	RcppGSL::vector<double> x(init);
	RcppGSL::matrix<double> ui1(ui); 
	RcppGSL::vector<double> ci1(ci); 
	double mu_n= 1; 
	double *mu = &mu_n; 
	struct constr_par mypar={ui1, ci1, mu, p}; 
	gsl_vector *df; 
	df = gsl_vector_alloc(nt); 
		
	double out = RC(x, &mypar); 
	dRC(x, &mypar, df); 
	for(size_t i=0; i<n; i++){
		Rcout<<gsl_vector_get(df, i)<<","; 
	}
	
	gsl_vector_free(ci1); gsl_vector_free(x); gsl_vector_free(df); 
	gsl_matrix_free(ui1); 
	return(out);
}

// [[Rcpp::export]]
NumericVector constrOptimC(NumericVector init, NumericMatrix ui, NumericVector ci){
	double p[4] = {0, 5, 0, .5};
	int n = init.length(), m=ci.length(); 
	size_t nt = init.length(), mt = ci.length();  
	
	RcppGSL::vector<double> x(init);
	RcppGSL::matrix<double> ui1(ui); 
	RcppGSL::vector<double> ci1(ci); 
	double mu_n= 1; 
	double *mu = &mu_n; 
	struct constr_par mypar={ui1, ci1, mu, p}; 
	gsl_set_error_handler_off();
	
	// Initiate multimin minimizer; 
	const gsl_multimin_fdfminimizer_type *T;
  	gsl_multimin_fdfminimizer *s;	
	T = gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (T, n);
	
	// Claim minimizing function; 
  	gsl_multimin_function_fdf my_func;
  	my_func.n = n;
	my_func.f = RC;
	my_func.df = dRC;
	my_func.fdf = RdRC;
	my_func.params = &mypar;
	
	// Claim the original function; 
	gsl_multimin_function_fdf orig_func;
  	orig_func.n = n;
	orig_func.f = FC;
	orig_func.df = gradC;
	orig_func.params = p;
	
	// Initiate the constrained optimization structure; 
	struct constr_multimin my_constrOpitm = {s, my_func, orig_func, 1e-4, 100, 1e-05}; 
	NumericVector out(n+1); 
	gsl_vector *opt;
	opt = gsl_vector_alloc(nt); 
	double f; 
	f = my_constrOpitm.inner_multimin(x, opt);
	out = my_constrOpitm.constr_optim(x); 
	
	gsl_vector_free(ci1); gsl_vector_free(x); gsl_vector_free(opt); 
	gsl_matrix_free(ui1); 
	return(out);
}