/*
Estimation of linear quantile models and linear quantile mixed models (ver 1.0)

Marco Geraci, University College London with contributions by Matteo Bottai, Karolinska Institutet

This suite of routines is part of the R package lqmm. The optimization algorithm is described in "A gradient search maximization algorithm for Laplace likelihood" (Bottai, Orsini and Geraci, 2011). The models stem from the work by Geraci (Ph.D. Dissertation, University of Florence, 2005) and by Geraci and Bottai ("Quantile regression for longitudinal data using the asymmetric Laplace distribution", Biostatistics 8, 2007) in which the asymmetric Laplace likelihood is linked to the estimation of conditional quantiles of a response variable given covariates and subject-specific random intercepts. The estimation of its parameters entails unconstrained maximization of a concave and nondifferentiable function over the real space. The algorithm is based on the gradient of the log-likelihood that generates a finite sequence of parameter values along which the likelihood increases. Its application to quantile mixed models is described in "Linear quantile mixed models" (Geraci and Bottai, 2011, Unpublished manuscript).
*/

#include <Rinternals.h>
#include <Rmath.h>
#include <R.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/RS.h>

double MIN(double x[], int length)
{
	int i;
	double minv = x[0];
	
	for (i = 1; i < length; i++)
	{
		if (minv > x[i])
			minv = x[i];
	}
	return minv;
}

double MAX(double x[], int length)
{
	int i;
	double maxv = x[0];
	
	for (i = 1; i < length; i++)
	{
		if (maxv < x[i])
			maxv = x[i];
	}
	return maxv;
}

/*

GRADIENT SEARCH FOR LAPLACE LIKELIHOOD

INDEPENDENT DATA

*/

double ll_s_d(double *theta, double *x, double *y, float *quantile, int *N, int *p, double *deriv)
{

/*  Likelihood (objective) function and gradient for independent data
	
	theta = parameter vector
	x = design matrix fixeff (vectorized by column)	
	y = response
	quantile = tau
	p = length theta_x
	N = total number of observations
*/
	
	int i, loc_N;
	loc_N = *N;
	double I[loc_N], resid[loc_N];
	double ans = 0;
	double alpha = -1.0; double beta = 1.0; int inc = 1;
	char *transn = "N";
	char *transt = "T";

	// resid <-  y
	F77_CALL(dcopy)(N, y, &inc, resid, &inc);
	// resid <- resid - X%*%theta
	F77_CALL(dgemv)(transn, N, p, &alpha, x, N, theta, &inc, &beta, resid, &inc);
   
	for (i = 0; i < (*N); i++){
		if (resid[i] < 0)
			{I[i] = *quantile - 1;}
		else
			{I[i] = *quantile;}
		ans += resid[i]*I[i];
	}
	// deriv <- -t(X)%*%I (p times 1)
	beta = 0.0;
	F77_CALL(dgemv)(transt, N, p, &alpha, x, N, I, &inc, &beta, deriv, &inc);

	return ans;
}

void gradientSd_s(double theta[], double *x, double *y, float *quantile, int *N, int *p, double *STEP, double *beta, double *gamma, int *RESET_STEP, double *TOL_LL, double *TOL_THETA, int *CHECK_THETA, int *MAXIT, int *verbose, int CONVERGE[], double deriv[], double optimum[])
{

/*
Optimization algorithm for quantile models with independent data
*/
	
	double f_0, f_1, tol_ll, tol_theta, step;
	int j, loc_p, iter, maxit, check_theta;
	int inc = 1;
	loc_p = *p;
	double theta_tmp[loc_p], grad[loc_p], test[loc_p];
  
	CONVERGE[0] = -2;
	maxit = *MAXIT;
	tol_ll = *TOL_LL;
	tol_theta = *TOL_THETA;
	step = *STEP;
	iter = 0;
	check_theta = 0;
	
	//initialize loglikelihood and gradient with theta_0
	f_0 = ll_s_d(theta, x, y, quantile, N, p, deriv);

	for(j = 0; j < loc_p; j++){
    grad[j] = deriv[j];
	}

	//NOTE: this loop MINIMIZES the negative log-likelihood
	for(iter = 0; iter < maxit; iter++){//loop
	  
		//update parameter
		for(j = 0; j < loc_p; j++){
		  theta_tmp[j] = theta[j] - grad[j]*step;
		}

		//update likelihood and gradient
		f_1 = ll_s_d(theta_tmp, x, y, quantile, N, p, deriv);
		if(*verbose == 1) {Rprintf("  (%i",iter+1); Rprintf(")"); Rprintf(" logLik = %5.12f\n", f_1);}
		
		//line search
		if (f_1 > f_0) {
			step = step * (*beta);
			if(*verbose == 1) Rprintf("  Decreasing step...\n");
		} else		
		{
		
			if(*CHECK_THETA == 1){
				for(j = 0; j < loc_p; j++){
				test[j] = fabs(theta_tmp[j]-theta[j]);
				}
				if (MAX(test, loc_p) < tol_theta) {check_theta = 0;} else {check_theta = 1;}
			}

			if (fabs(f_1 - f_0) < tol_ll && check_theta == 0) {goto OUT_FOR;}//exit loop
				
			F77_CALL(dcopy)(&loc_p, theta_tmp, &inc, theta, &inc);
			F77_CALL(dcopy)(&loc_p, deriv, &inc, grad, &inc);
			f_0 = f_1;

			if(*RESET_STEP == 1) {step = *STEP;} else {step = step * (*gamma);}
		}
	} //end loop1
  


	OUT_FOR:;
	if(iter < maxit) {CONVERGE[0] = iter + 1;}
	if((iter == maxit) & (maxit > 0)){CONVERGE[0] = -1;}
   
	 optimum[0] = f_0;
}


/*

HIERARCHICAL DATA

*/

double psi_mat(double *theta, int *Tq, int *p, int *q, int *m, int *cov_type, double *psi){

/*
cov_type
0 pdIdent
1 pdDiag
2 pdCompSymm
3 pdSymm
*/

	if(*cov_type == 0){
	int j;
		for (j = 0; j < (*m); j++){
			psi[j] = theta[j + (*p)];
		}
	}

	if(*cov_type == 1){
	int j;
		for (j = 0; j < (*m); j++){
			psi[j] = theta[j + (*p)];
		}
	}

	if(*cov_type == 2){
	int i, j, qsq;
	qsq = (*q) * (*q);
		for (i = 0; i < qsq; i++){
			psi[i] = 0;
			for (j = 0; j < (*m); j++){
				psi[i] += Tq[i + qsq * j] * theta[j + (*p)];
				}
		}
	}
	
	if(*cov_type == 3){
	int i, j, qsq;
	qsq = (*q) * (*q);
		for (i = 0; i < qsq; i++){
			psi[i] = 0;
			for (j = 0; j < (*m); j++){
				psi[i] += Tq[i + qsq * j] * theta[j + (*p)];
				}
		}
	}
	
return 0;
}



double pdIdent(double *z, double *V, double *psi, int *q, int *N, int *Kq, int k, int j){

int l;
double val = 0;

	for (l = 0; l < (*q); l++){
		val += V[(k) + (*Kq) * l] * z[(j) + (*N) * l] * psi[0];	
	}

	return val;

}

double pdDiag(double *z, double *V, double *psi, int *q, int *N, int *Kq, int k, int j){

int l;
double val = 0;

	for (l = 0; l < (*q); l++){
		val += V[(k) + (*Kq) * l] * z[(j) + (*N) * l] * psi[l];	
	}

return val;

}

double pdCompSymm(double *z, double *V, double *psi, int *q, int *N, int *Kq, int k, int j){

int STEPQ, l;
double val = 0;

	for (STEPQ = 0; STEPQ < (*q); STEPQ++){
		for (l = 0; l < (*q); l++){
			val += V[(k) + (*Kq) * STEPQ] * z[(j) + (*N) * l] * psi[l + (*q) * STEPQ];
		}
	}

return val;

}

double pdSymm(double *z, double *V, double *psi, int *q, int *N, int *Kq, int k, int j){

int STEPQ, l;
double val = 0;

	for (STEPQ = 0; STEPQ < (*q); STEPQ++){
		for (l = 0; l < (*q); l++){
			val += V[(k) + (*Kq) * STEPQ] * z[(j) + (*N) * l] * psi[l + (*q) * STEPQ];
		}
	}


return val;

}


double lin_pred_ll(int *cov_type, double *z, double *V, double *psi, int *q, int *N, int *Kq, int k, int j){

double val = 0;

	if(*cov_type == 0){
	val = pdIdent(z, V, psi, q, N, Kq, k, j);
	}
	if(*cov_type == 1){
	val = pdDiag(z, V, psi, q, N, Kq, k, j);
	}
	if(*cov_type == 2){
	val = pdCompSymm(z, V, psi, q, N, Kq, k, j);
	}
	if(*cov_type == 3){
	val = pdSymm(z, V, psi, q, N, Kq, k, j);
	}

return val;

}



/*

GRADIENT SEARCH FOR LAPLACE LIKELIHOOD

HIERARCHICAL DATA

*/

void ll_h_R(double theta[], double *x, double *y, double *z, double *weights, int *Tq, int *cov_type, double *V, double *W, double *sigma, float *quantile, int *p, int *q, int *m, int *M, int *N, int *Kq, int *start, int *end, double *loglik)
{

/*  Likelihood (objective) function for hierarchical data - type void to be called from within R
	
	theta = parameter vector
	x = design matrix fixeff (vectorized by column)	
	y = response
	z = Z design matrix raneff (vectorized by column)
	Tq = duplication matrix (vectorized by column)
	V = matrix Kq x q nodes (vectorized by column)
	W = product weights in a_ik
	sigma = scale residuals
	quantile = tau
	p = length theta_x
	q = number raneff
	m = length theta_z
	M = number of clusters
	N = total number of observations
	Kq = K^q (no. of knots ^ number of random effects)
	start = vector of indices (e.g., start[i] = beginning observations cluster i)
	end = vector of indices (e.g., end[i] = end observations cluster i)

*/
	
	int i, k, j, ls, le, NN, KKq;
	double _cons, resid, val, mvalk, MAXVALK, TMP_VAR, tmp_psi;
	
	_cons = log((*quantile) * (1-(*quantile)) / *sigma);
	NN = *N;
	KKq = *Kq;
	
	int dim_psi = 1;
		if((*cov_type == 0) || (*cov_type == 1)){
			dim_psi = *m;
		}
		if((*cov_type == 1) || (*cov_type == 2)){
			dim_psi = (*q)*(*q);
		}
	
	double psi[dim_psi], w[NN], valk[KKq];

		
// y - X*theta_x
	
	for (i = 0; i < NN; i++){
		w[i] = 0;
		for (j = 0; j < *p; j++){
			w[i] += x[i + NN*j]*theta[j];
			}
		w[i] = y[i] - w[i];
	}

// Psi.sqroot = Tq*theta_z
	
	tmp_psi = psi_mat(theta, Tq, p, q, m, cov_type, psi);

// Asym. Laplace Likelihood	for M clusters
	
	for (i = 0; i < *M; i++)
	{
	mvalk = 0;
	MAXVALK = -1000000000;
	ls = start[i];
	le = end[i];
		for (k = 0; k < KKq; k++)
		{
		val = 0;
		//valk[k] = 0;
			for (j = ls; j < le; j++)
			{
				TMP_VAR = 0;
				TMP_VAR = lin_pred_ll(cov_type, z, V, psi, q, N, Kq, k, j);
				//Rprintf("TMP_VAR = %5.5f\n", TMP_VAR);
				resid = w[j] - TMP_VAR;
				val += _cons - (1 / (2 * (*sigma))) * (fabs(resid) + (2 * (*quantile) - 1) * resid);
			}
		//valk += exp(val + log(W[k]));
		valk[k] = val + log(W[k]);
		if(valk[k] > MAXVALK)	MAXVALK = valk[k];
		}

		for (k = 0; k < KKq; k++){
			mvalk += exp(valk[k] - MAXVALK);
		}
	
	loglik[0] += - weights[i] * (MAXVALK + log(mvalk));
	}
}


double ll_h_d(double *theta, double *x, double *y, double *z, double *weights, int *Tq, int *cov_type, double *V, double *W, double *sigma, float *quantile, int *p, int *q, int *m, int *M, int *N, int *Kq, int *start, int *end, double *deriv, double *loglik)
{

/*  Likelihood (objective) function and gradient for hierarchical data
	
	theta = parameter vector
	x = design matrix fixeff (vectorized by column)	
	y = response
	z = Z design matrix raneff (vectorized by column)
	Tq = duplication matrix (vectorized by column)
	V = matrix Kq x q nodes (vectorized by column)
	W = product weights in a_ik
	sigma = scale residuals
	quantile = tau
	p = length theta_x
	q = number raneff
	m = length theta_z
	M = number of clusters
	N = total number of observations
	Kq = K^q (no. of knots ^ number of random effects)
	start = vector of indices (e.g., start[i] = beginning observations cluster i)
	end = vector of indices (e.g., end[i] = end observations cluster i)
*/
	
	int i, k, j, l, ls, le, STEPQ, STEPM, NN, KKq, n, qsq;
	double _cons, resid, val, mvalk, phi, TMP_VAR, tmp_psi, MAXVALK;
	
	_cons = log((*quantile) * (1-(*quantile)) / *sigma);
	NN = *N;
	KKq = *Kq;
	n = *p + *m;
	qsq = (*q)*(*q);

	int dim_psi = 1;
	if((*cov_type == 0) || (*cov_type == 1)){
		dim_psi = *m;
	}
	if((*cov_type == 1) || (*cov_type == 2)){
		dim_psi = (*q)*(*q);
	}
	
	double psi[dim_psi], w[NN], GRAD[n], DERIV[n], valk[KKq];

	for (l = 0; l < n; l ++){
	 GRAD[l] = 0;
	 DERIV[l] = 0;
	 deriv[l] = 0;
	}
        	
// y - X*theta_x
	
	for (i = 0; i < (*N); i++){
		w[i] = 0;
		for (j = 0; j < *p; j++){
			w[i] += x[i + (*N)*j]*theta[j];
			}
		w[i] = y[i] - w[i];
	}

// Psi.sqroot = Tq*theta_z
	
	tmp_psi = psi_mat(theta, Tq, p, q, m, cov_type, psi);
	
// Asym. Laplace Likelihood	for M clusters
	
	loglik[0] = 0;
	//mvalk = 0;
	for (i = 0; i < *M; i++)
	{
		mvalk = 0;
		MAXVALK = -1000000000;
		ls = start[i];
		le = end[i];
		for (k = 0; k < KKq; k++)
		{
		 val = 0;
			for (j = ls; j < le; j++)
			{
				TMP_VAR = 0;
				TMP_VAR = lin_pred_ll(cov_type, z, V, psi, q, N, Kq, k, j);
				//Rprintf("TMP_VAR = %5.5f\n", TMP_VAR);
				resid = w[j] - TMP_VAR;

				phi = *quantile;
				if (resid <= 0) phi = *quantile - 1;
				val += _cons - (1 / (*sigma)) * resid * phi;

				//calculate numerator derivative    
				for (l = 0; l < *p; l++){
				 GRAD[l] += x[j + (*N) * l]*phi;
				}
 				
				for (STEPM = 0; STEPM < *m; STEPM++){
					TMP_VAR = 0;
					for (STEPQ = 0; STEPQ < *q; STEPQ++){
						for (l = 0; l < *q; l++){
							TMP_VAR += V[k + (*Kq)*STEPQ]*z[j + (*N) * l]*Tq[l + (*q) * STEPQ + STEPM * (qsq)];
						}
					}
					GRAD[STEPM + (*p)] += TMP_VAR*phi;
				}

			} //j
			valk[k] = val + log(W[k]);

			for (l = 0; l < n; l++) {
				DERIV[l] += GRAD[l] * exp(valk[k])/(*sigma);
				GRAD[l] = 0;
			}

			if(valk[k] > MAXVALK)	MAXVALK = valk[k];
		} //k
  
		for (k = 0; k < KKq; k++){
			mvalk += exp(valk[k] - MAXVALK);
		}


		loglik[0] += - weights[i] * (MAXVALK + log(mvalk));

		for (l = 0; l < n; l++) {
			deriv[l] += - weights[i] * DERIV[l]/(mvalk*exp(MAXVALK));
			DERIV[l] = 0;
			//Rprintf(" logLik = %5.5f\n", deriv[l]);
		}
	
	}//i
return loglik[0];
}


void gradientSd_h(double theta[], double *x, double *y, double *z, double *weights, int *Tq, int *cov_type, double *V, double *W, double *sigma, float *quantile, int *p, int *q, int *m, int *M, int *N, int *Kq, int *start, int *end, double *STEP, double *beta, double *gamma, int *RESET_STEP, double *TOL_LL, double *TOL_THETA, int *CHECK_THETA, int *MAXIT, int *verbose, int CONVERGE[], double loglik[], double deriv[], double optimum[])
{

/*
Optimization algorithm for quantile models with hierarchical data
*/
	
	double f_0, f_1, tol_ll, tol_theta, step;
	int j, n, iter, maxit, check_theta;
	n = *p + *m;
	double theta_tmp[n], grad[n], test[n];

	CONVERGE[0] = -2;
	maxit = *MAXIT;
	tol_ll = *TOL_LL;
	tol_theta = *TOL_THETA;
	step = *STEP;
	iter = 0;
	check_theta = 0;

	//initialize loglikelihood and gradient with theta_0
	f_0 = ll_h_d(theta, x, y, z, weights, Tq, cov_type, V, W, sigma, quantile, p, q, m, M, N, Kq, start, end, deriv, loglik);

	for(j = 0; j < n; j++){
		grad[j] = deriv[j];
	}
  
	if(*verbose == 1) {
		Rprintf(" Lower loop\n");
	}

	//NOTE: this loop MINIMIZES the negative log-likelihood
	for(iter = 0; iter < maxit; iter++){//loop
  
		for(j = 0; j < n; j++){
		theta_tmp[j] = theta[j] - grad[j]*step;
		}

		f_1 = ll_h_d(theta_tmp, x, y, z, weights, Tq, cov_type, V, W, sigma, quantile, p, q, m, M, N, Kq, start, end, deriv, loglik);
    
		if(*verbose == 1) {
			Rprintf("  (%i",iter+1); Rprintf(")"); Rprintf(" logLik = %5.5f\n", f_1);
		}

		if (f_1 > f_0) {
			step = step * (*beta);
			if(*verbose == 1) Rprintf("  Decreasing step...\n");
		} else
		{

			if(*CHECK_THETA == 1){
				for(j = 0; j < n; j++){
				test[j] = fabs(theta_tmp[j]-theta[j]);
				}
				if (MAX(test,n) < tol_theta) {check_theta = 0;} else {check_theta = 1;}
			}

		
		if (fabs(f_1 - f_0) < tol_ll && check_theta == 0) {goto OUT_FOR;}//exit loop

		for(j = 0; j < n; j++){
			theta[j] = theta_tmp[j];
		}
		f_0 = f_1;
		for(j = 0; j < n; j++){
			grad[j] = deriv[j];
		}
		if(*RESET_STEP == 1) {step = *STEP;} else {step = step * (*gamma);}
		}
	}//end loop


	OUT_FOR:;

	if(iter < maxit) {CONVERGE[0] = iter + 1;}
	if((iter == maxit) & (maxit > 0)){CONVERGE[0] = -1;}
	
	optimum[0] = f_0;
}
