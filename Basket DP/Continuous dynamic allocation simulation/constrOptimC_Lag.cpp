#include <RcppGSL.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppGSL)]]

struct constr_par{const gsl_matrix *ui; const gsl_vector *ci; void *p ; const double mu;};

struct constr_multimin{
	gsl_multimin_fdfminimizer *s; gsl_multimin_function_fdf F; gsl_multimin_function_fdf orig_F; 
	double tol; int max_iter;  struct constr_par constr_parl; 
	
	double inner_multimin(gsl_vector *x, gsl_vector *par){
		int iter = 0, status; 
	
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
			Rcout<<"Iteration"<<iter<<": gradient="<<gsl_vector_get(s->gradient, 0)<<","<<gsl_vector_get(s->gradient, 1)<<
			", "<<gsl_vector_get(s->gradient, 2)<<","<<gsl_vector_get(s->gradient, 3)<<","<<gsl_vector_get(s->gradient, 4)<<"\n";
		}while (status == GSL_CONTINUE && iter < 200);
		
		// Return the minimizer and optimal value; 
		Rcout<<"gradient="<<gsl_vector_get(s->gradient, 0)<<","<<gsl_vector_get(s->gradient, 1)<<
			", "<<gsl_vector_get(s->gradient, 2)<<","<<gsl_vector_get(s->gradient, 3)<<","<<gsl_vector_get(s->gradient, 4)<<"\n";
		gsl_vector_memcpy(par, s->x); 
		return(s->f );
	}
	
	NumericVector constr_optim(gsl_vector *x0){
		size_t n = orig_F.n, n1 = F.n; 
		size_t m = n1 - n; 
		gsl_vector *x1, *x1star; 
		double r=0; 
		NumericVector out(n+1); 
		
		// Assign initial values; 
		x1 = gsl_vector_alloc(n1);
		x1star = gsl_vector_alloc(n1); 
		for(size_t i=0; i<n; i++){
			gsl_vector_set(x1, i, gsl_vector_get(x0, i)); 
		}
		for(size_t i=n; i<n1; i++){
			gsl_vector_set(x1, i, -10); 
		}
		
		// Lagrangian optimization; 
		r = inner_multimin(x1, x1star); 
		Rcout<<"log(lambda)="<<gsl_vector_get(x1star,n+1)<<","<<gsl_vector_get(x1star,n+2)<<"\n"; 
		
		// Return optimal values; 
		for(int i=0; i<n; i++){
			gsl_vector_set(x0, i, gsl_vector_get(x1star,i)); 
			out(i) = gsl_vector_get(x1star, i);
		}
		out(n) = orig_F.f(x0, orig_F.params);
		gsl_vector_free(x1); gsl_vector_free(x1star); 
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

// Define Lagrangian function and gradients; 
double RC(const gsl_vector *x, void *params){
	struct constr_par *par	= (struct constr_par *)params;
	const gsl_matrix *ui	= (par->ui);
	const gsl_vector *ci	= (par->ci);
	void *p1				= (par->p);
	double mu 				= (par->mu);
	double out=0; 
	
	size_t n 	= ui->size2, m = ui->size1; 
	size_t n1 	= n + m; 
	gsl_vector *x0, *lambda, *gi;
	x0			= gsl_vector_alloc(n); 
	lambda		= gsl_vector_alloc(m); 
	gi			= gsl_vector_alloc(m); 
	
	for(size_t i=0; i<n; i++){
		gsl_vector_set(x0, i, gsl_vector_get(x, i)); 
	}
	for(size_t i=0; i<m; i++){
		gsl_vector_set(lambda, i, exp(gsl_vector_get(x,i+n))); 
	}
	
	gsl_vector_memcpy(gi, ci); 
	gsl_blas_dgemv(CblasNoTrans, 1.0, ui, x0, -1.0, gi); 
	for(size_t i=0; i<m; i++){
		double d = gsl_vector_get(gi,i)*log(gsl_vector_get(gi,i)); 
		gsl_vector_set(gi,i, d); 
	}
	
	double bar; 
	gsl_blas_ddot(lambda, gi, &bar);
	out = FC(x0, p1) - bar; 
	
	gsl_vector_free(x0); gsl_vector_free(lambda); 
	gsl_vector_free(gi); 
	return(out);
}

void dRC(const gsl_vector *x, void *params, gsl_vector *df){
	struct constr_par *par	= (struct constr_par *)params;
	const gsl_matrix *ui	= (par->ui);
	const gsl_vector *ci	= (par->ci);
	void *p1				= (par->p);	
	double mu 				= (par->mu);
	
	size_t n 	= ui->size2, m = ui->size1; 
	size_t n1 	= m + n ; 
	gsl_vector *x0, *lambda, *gi, *dfx0;
	x0			= gsl_vector_alloc(n); 
	lambda		= gsl_vector_alloc(m); 
	gi			= gsl_vector_alloc(m); 
	dfx0		= gsl_vector_alloc(n); 
	
	for(size_t i=0; i<n; i++){
		gsl_vector_set(x0, i, gsl_vector_get(x, i)); 
	}
	for(size_t i=0; i<m; i++){
		gsl_vector_set(lambda, i, exp(gsl_vector_get(x,i+n))); 
	}
	
	gsl_vector_memcpy(gi, ci); 
	gsl_blas_dgemv(CblasNoTrans, 1.0, ui, x0, -1.0, gi); 
	
	// Compute gradients; 
	gradC(x, p1, dfx0); 
	gsl_blas_dgemv(CblasTrans, -1.0, ui, lambda, 1, dfx0);
	
	for(size_t i=0; i<n; i++){
		gsl_vector_set(df, i, gsl_vector_get(dfx0, i));  
	}
	for(size_t i=0; i<m; i++){
// 		gsl_vector_set(df, i+n, (gsl_vector_get(gi,i)*gsl_vector_get(lambda,i) - mu)* gsl_vector_get(lambda,i)); 
		gsl_vector_set(df, i+n, (gsl_vector_get(gi,i)*gsl_vector_get(lambda,i) - mu)); 
	}
	
	gsl_vector_free(x0); gsl_vector_free(lambda); gsl_vector_free(gi); gsl_vector_free(dfx0); 
}

void RdRC(const gsl_vector *x, void *params, double *f, gsl_vector *df){
	*f = RC(x, params); 
	dRC(x, params, df);
}

// A wrapper of constrained optimization; 
// [[Rcpp::export]]
double my_fn(NumericVector init, NumericMatrix ui, NumericVector ci){
	gsl_set_error_handler_off();
	double p[4] = {0, 5, 0, .5};
	size_t n = init.length(), m = ci.length();  
	size_t n1 = n+m; 
	
	RcppGSL::vector<double> x(init);
	RcppGSL::matrix<double> ui1(ui); 
	RcppGSL::vector<double> ci1(ci); 
	struct constr_par mypar={ui1, ci1, p, 1e-4}; 
	gsl_vector *x1, *df; 
	x1 = gsl_vector_alloc(n1); 
	df = gsl_vector_alloc(n1); 
	
	for(size_t i=0; i<n; i++){
		gsl_vector_set(x1, i, gsl_vector_get(x, i)); 
	}
	for(size_t i=n; i<n1; i++){
		gsl_vector_set(x1, i, -50); 
	}
	
	double out = RC(x1, &mypar); 
	dRC(x1, &mypar, df); 
	for(size_t i=0; i<n1; i++){
		Rcout<<gsl_vector_get(df, i)<<","; 
	}
	
	gsl_vector_free(ci1); gsl_vector_free(x); gsl_vector_free(df); gsl_vector_free(x1);
	gsl_matrix_free(ui1); 
	return(out);
}

// [[Rcpp::export]]
NumericVector constrOptimC(NumericVector init, NumericMatrix ui, NumericVector ci){
	gsl_set_error_handler_off();
	double p[4] = {0, 5, 0, .5};
	size_t n = init.length(), m = ci.length();  
	size_t n1 = n+m; 
	
	RcppGSL::vector<double> x(init);
	RcppGSL::matrix<double> ui1(ui); 
	RcppGSL::vector<double> ci1(ci); 
	struct constr_par mypar={ui1, ci1, p, 1e-4}; 

	// Initiate multimin minimizer; 
	const gsl_multimin_fdfminimizer_type *T;
  	gsl_multimin_fdfminimizer *s;	
	T = gsl_multimin_fdfminimizer_vector_bfgs;
	s = gsl_multimin_fdfminimizer_alloc (T, n1);
	
	// Claim minimizing function; 
  	gsl_multimin_function_fdf my_func;
  	my_func.n = n1;
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
	struct constr_multimin my_constrOpitm = {s, my_func, orig_func, 1e-4, 100, mypar}; 
	NumericVector out(n+1); 
// 	gsl_vector *opt;
// 	opt = gsl_vector_alloc(n1); 
	out = my_constrOpitm.constr_optim(x); 
	
	gsl_vector_free(ci1); gsl_vector_free(x); 
// 	gsl_vector_free(opt); 
	gsl_matrix_free(ui1); 
	gsl_multimin_fdfminimizer_free(s);
	
	return(out);
}