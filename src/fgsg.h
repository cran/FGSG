#ifndef FGSG_H
#define FGSG_H

#include <R.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include "lafunc.h"


typedef struct{
	integer n,p,g,maxIter, aMaxIter;
	double rho,tol,tau, aTol;
	double* wt;
}Parm;

void computeTTx(integer *edge, double *x, double *TTx, integer g, double scale, double *wt);

void computeDegree(integer *edge, double *d, integer g, double *wt);

/* regularization version */
void goscarSub(double *x, double* Ch_R, double* bv, integer* edge, double scale,double s1, double s2,Parm opts);

void ncTFGS(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts);

void ncFGS(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts);

void goscar(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts);

void gflassoSub(double *x, double* Ch_R, double* bv, integer* edge, double s1, double s2,Parm opts,  ptrdiff_t *sign);

void gflasso(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts);

void ncTL(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts);

void ncTF(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts);

void ncTLF(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts);


/* functions */
/* regularization functions */
void computeTTx(integer *edge, double *x, double *TTx, integer g, double scale, double *wt)
{
	integer i;
	for (i = 0; i< g; i++)
	{
		TTx[edge[2*i]] += (x[2*i] + x[2*i + 1])*scale*wt[i];
		TTx[edge[2*i+1]] += (x[2*i + 1] - x[2*i])*scale*wt[i];  
	} 
}

void computeDegree(integer *edge, double *d, integer g, double *wt)
{
	integer i;
	for (i = 0; i< 2*g; i++)
		d[edge[i]] +=wt[(integer)i/2]*wt[(integer)i/2];
}

void goscarSub(double *x, double* Ch_R, double* bv, integer* edge, double scale, double s1, double s2, Parm opts)
{
	integer iter,i;
	ptrdiff_t info;
	integer nrhs;
	double tv,td,tz;
	double *x_p,*q, *pk, *v, *u, *Tx,*funv,*TTx;
	char UPLO = 'U';
	nrhs = 1;

	x_p = (double*)R_alloc(opts.p,sizeof(double));
	q   = (double*)R_alloc(opts.p,sizeof(double));
	pk  = (double*)R_alloc(2*opts.g,sizeof(double));
	v   = (double*)R_alloc(2*opts.g,sizeof(double));
	u   = (double*)R_alloc(opts.p,sizeof(double));
	Tx  = (double*)R_alloc(2*opts.g,sizeof(double));
	funv = (double*)R_alloc(opts.aMaxIter,sizeof(double));
	TTx = (double*)R_alloc(opts.p,sizeof(double));

	memset(x_p, 0, opts.p*sizeof(double));
	memset(q,   0, opts.p*sizeof(double));
	memset(pk,  0, 2*opts.g*sizeof(double));
	memset(v,   0, 2*opts.g*sizeof(double));
	memset(u,   0, opts.p*sizeof(double));
	memset(funv,0, opts.aMaxIter*sizeof(double));

	for(iter = 0; iter < opts.aMaxIter; iter ++)
	{
		memset(TTx,0,opts.p*sizeof(double));
		computeTTx(edge,pk,TTx,opts.g,opts.rho*scale, opts.wt);
		computeTTx(edge,v,TTx,opts.g,-scale, opts.wt);
		for(i = 0; i< opts.p; i++)
			x[i] = bv[i] - u[i] + opts.rho*q[i] + TTx[i];

		dpotrs_(&UPLO, &opts.p, &nrhs, Ch_R, &opts.p, x, &opts.p, &info);

		for (i = 0; i < opts.p ; i ++)
		{
			td = x[i] + u[i]/opts.rho; 
			if (td < -s1/opts.rho)
				q[i] = td + s1/opts.rho;
			else if(td > s1/opts.rho)
				q[i] = td - s1/opts.rho;
			else
				q[i]= 0;

			u[i] += opts.rho*(x[i] - q[i]);
		}

		for (i = 0; i< opts.g; i++)
		{
			td = x[edge[2*i]];
			tv = x[edge[2*i+1]];
			Tx[2*i] = (td - tv)*scale*opts.wt[i];
			Tx[2*i+1] = (td + tv)*scale*opts.wt[i];
		}

		for (i = 0; i<2*opts.g; i++)
		{
			td = Tx[i] + v[i]/opts.rho;
			if (td < -s2/opts.rho)
				pk[i] = td + s2/opts.rho;
			else if(td > s2/opts.rho)
				pk[i] = td - s2/opts.rho;
			else
				pk[i]= 0;

			v[i] += opts.rho*(Tx[i] - pk[i]);
		}

		tv = 0; td = 0;
		for (i = 0; i < opts.p; i++)
		{
			tv += pow(x[i] - q[i],2);
			td += pow(x[i] - x_p[i], 2);
		}
		funv[iter] = sqrt(tv);
		/*printf("The Objective value in the %td-th iteration is %f\opts.n", iter+1, (float)funv[iter]);*/
		td = sqrt(td);
		/*printf("The gap between x and x_p in the %td-th iteration is %f\opts.n", iter+1, (float)td); */

		if (iter >0)
		{
			tv = funv[iter-1] < 1 ? 1 : funv[iter-1];
			tz = td<1 ? 1:td;
			if (funv[iter] < opts.aTol *tv && td < opts.aTol *tz)
				break;
		}
		memcpy(x_p,x,opts.p*sizeof(double));
	}
    
	memcpy(x,q,opts.p*sizeof(double));
/*	free(x_p);
	free(q);
	free(pk);
	free(v);
	free(u);
	free(Tx);
	free(funv);
	free(TTx);*/
}

void ncTFGS(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts)
{
	integer i,j,k,iter; 
	double *bv, *Ax, *cp, *Ch_R,*funVal,*degree;
	double f1,td,tv,tz;
	double scale = 1;
	char trans = 'T';
	char Tra = 'T';
	char Trb= 'N';
	char UPLO = 'U';

	ptrdiff_t info = 0;
	integer incx = 1;
	integer incy = 1;
	double alpha = 1.0;
	double beta = 0;

	degree = (double*)R_alloc(opts.p,sizeof(double));
	memset(degree,0,opts.p*sizeof(double));
	computeDegree(edge,degree,opts.g, opts.wt);

	Ax = (double*)R_alloc(opts.n,sizeof(double));
	bv = (double*)R_alloc(opts.p,sizeof(double));
	memset(bv, 0, opts.p*sizeof(double));
	dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, y, &incx, &beta, bv, &incy);

	/* compute the chol factoration */
	Ch_R = (double*) R_alloc(opts.p*opts.p,sizeof(double));
	memset(Ch_R, 0, opts.p*opts.p*sizeof(double));

	/* set Ch_R as rho * eye(p) + rho*T'*T */
	for(i = 0; i< opts.p; i++)
		Ch_R[i*opts.p+i] = opts.rho*(1 + degree[i]*2);

	/* compute Ch_R = (A'*A) + rho * eye(p) + rho*T'*T;*/
	beta = 1.0;
	alpha = 1.0;
	dgemm_(&Tra, &Trb,&opts.p, &opts.p, &opts.n, &alpha, A, &opts.n, A, &opts.n, &beta, Ch_R, &opts.p);
	dpotrf_(&UPLO, &opts.p, Ch_R, &opts.p, &info);

	if (info !=0)
	{
		
		/*free(bv);
		free(Ax);
		free(degree);
		free(Ch_R);*/
		error("Chol factorization fails");
		return;
	}  

	cp = (double*)R_alloc(opts.p,sizeof(double));

	/* compute Ax - y */
	trans = 'N';
	alpha = 1.0;
	beta = 0.0;
	/* assign memeory for funVal */
	funVal = (double*) R_alloc(opts.maxIter,sizeof(double));
	memset(funVal,0,opts.maxIter*sizeof(double));

	s1 /= opts.tau;
	s2 /= opts.tau;

	for(iter = 0; iter<opts.maxIter; iter++)
	{
		memset(cp,0,opts.p*sizeof(double));
		for (i = 0; i<opts.g; i++)
		{
			j = edge[2*i];
			k = edge[2*i+1];

			if (fabs(x[j]) - fabs(x[k]) > opts.tau)
				cp[j] += 2.0*s2*opts.wt[i];
			else if (fabs(x[k])-fabs(x[j]) > opts.tau)
				cp[k] += 2.0*s2*opts.wt[i];

			if (fabs(fabs(x[j]) - fabs(x[k])) < opts.tau)
			{
				cp[j] += s2*opts.wt[i];
				cp[k] += s2*opts.wt[i];
			}
		}
		for (i = 0; i<opts.p; i++)
		{
			if (x[i] == 0)
			{
				cp[i] = bv[i];
				continue;
			}

			if (x[i]<0)
				cp[i] *= -1.0;

			if (x[i]>opts.tau)
				cp[i] += s1;
			else if (x[i]<-opts.tau)
				cp[i] -= s1;

			cp[i] += bv[i];
		}

		goscarSub(x,Ch_R,cp,edge,scale,s1,s2,opts);

		dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, x, &incx, &beta, Ax, &incy);

		f1 =0; td = 0;
		/* compute f1 */
		for (i = 0; i<opts.n; i++)
			f1 += (Ax[i]-y[i])*(Ax[i]-y[i]);
		f1 /= 2.0;
		
		for (i = 0; i<opts.p; i++)
		{
		   tv = fabs(x[i]) < opts.tau ? fabs(x[i]):opts.tau;
		   td +=tv;
		}
		
		f1 += s1*td;
		td = 0;
		for (i = 0; i<opts.g; i++)
		{
		   tz = fabs(fabs(x[edge[2*i]]) - fabs(x[edge[2*i+1]]));
		   if (tz > opts.tau)
		       tz = opts.tau;
		   td +=tz*opts.wt[i];
		}
		f1 += s2*td;
		
        funVal[iter]  = f1;
		if (iter<1)
			continue;
	    tv = funVal[iter-1] < 1 ? 1 : funVal[iter-1];
		/*if (fabs(funVal[iter] - funVal[iter-1])/funVal[iter] < opts.tol)*/
		if (fabs(funVal[iter] - funVal[iter-1]) < opts.tol * tv)
			break;
	}
	
/*	free(bv);
	free(Ax);
	free(cp);
	free(Ch_R);
	free(funVal);
	free(degree);*/
}

void ncFGS(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts)
{
	integer i,iter,j; 
	double *bv, *Ax, *cp, *Ch_R, *funVal,*degree;
	double f1,f2,td,tv;
	double scale = 1;
	char trans = 'T';
	char Tra = 'T';
	char Trb = 'N';
	char UPLO = 'U';

	ptrdiff_t info = 0;
	integer incx = 1;
	integer incy = 1;
	double alpha = 1.0;
	double beta = 0;

	degree = (double*)R_alloc(opts.p,sizeof(double));
	memset(degree,0,opts.p*sizeof(double));
	computeDegree(edge,degree,opts.g,opts.wt);

	Ax = (double*) R_alloc(opts.n,sizeof(double));
	bv = (double*) R_alloc(opts.p,sizeof(double));
	memset(bv, 0, opts.p*sizeof(double));
	dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, y, &incx, &beta, bv, &incy);
	/* compute the chol factoration */
	Ch_R = (double*) R_alloc(opts.p*opts.p,sizeof(double));
	memset(Ch_R, 0, opts.p*opts.p*sizeof(double));

	/* set Ch_R as  opts.rho * eye(opts.p) + opts.rho*T'*T */
	for(i = 0; i< opts.p; i++)
		Ch_R[i*opts.p+i] = opts.rho*(1 + degree[i]*2);
	/* compute Ch_R = (A'*A) + opts.rho * eye(opts.p) + opts.rho*T'*T; */
	beta = 1.0;
	alpha = 1.0;
	dgemm_(&Tra, &Trb,&opts.p, &opts.p, &opts.n, &alpha, A, &opts.n, A, &opts.n, &beta, Ch_R, &opts.p);

	dpotrf_(&UPLO, &opts.p, Ch_R, &opts.p, &info);

	if (info != 0)
	{

	/*	free(bv);
		free(Ax);
		free(degree);
		free(Ch_R);*/
		error("Chol factorization fails");
		return;
	}  

	cp = (double*)R_alloc(opts.p,sizeof(double));
	memset(cp,0,opts.p*sizeof(double));

	/* compute Ax - y */
	trans = 'N';
	alpha = 1.0;
	beta = 0.0;
	/* assign memeory for funVal */
	funVal = (double*) R_alloc(opts.maxIter,sizeof(double));
	memset(funVal,0,opts.maxIter*sizeof(double));

	for(iter = 0; iter<opts.maxIter; iter++)
	{
        memcpy(cp,bv,opts.p*sizeof(double));
		
		for (i = 0; i<2*opts.g; i++)
		{
			j = edge[i];

			if (x[j] > 0)
				cp[j] += s2*opts.wt[i/2];
			else if (x[j] < 0)
				cp[j] -= s2*opts.wt[i/2];
		}

		goscarSub(x,Ch_R,cp,edge,scale,s1,s2,opts);
		dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, x, &incx, &beta, Ax, &incy);

		f1 =0; f2 = 0; td = 0;
		/* compute f1 */
		for (i = 0; i<opts.n; i++)
			f1 += (Ax[i]-y[i])*(Ax[i]-y[i]);
		f1 /= 2.0;
		for (i = 0; i<opts.p; i++)
			td += fabs(x[i]);
		f1 += s1*td;
		td = 0; tv = 0;
		for (i=0; i < opts.g; i++)
		{
			td += fabs(x[edge[2*i]] - x[edge[2*i+1]])*opts.wt[i];
			td += fabs(x[edge[2*i]] + x[edge[2*i+1]])*opts.wt[i];
			tv += (fabs(x[edge[2*i]]) + fabs(x[edge[2*i+1]]) ) *opts.wt[i];
		}
		f1 += s2*td; 
		f2 = s2*tv;
		funVal[iter] = f1 - f2;

		if (iter<1)
		  continue;
		tv = funVal[iter-1] < 1 ? 1 : funVal[iter-1];  
		if (fabs(funVal[iter] - funVal[iter-1]) < opts.tol*tv)
				break;
			/*if (fabs(funVal[iter] - funVal[iter-1])/funVal[iter] < opts.tol)*/	
	}
	
/*	free(degree);
	free(bv);
	free(Ax);
	free(cp);
	free(Ch_R);
	free(funVal);*/
}

void goscar(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts)
{
	integer i,k;
	double *bv, *Ch_R, *degree;
	char trans = 'T';
	char Tra = 'T';
	char Trb = 'N';
	char UPLO = 'U';

	ptrdiff_t info = 0;   
	integer incx = 1;
	integer incy = 1;
	double alpha = 1.0;
	double beta = 0;

	degree = (double*)R_alloc(opts.p,sizeof(double));
	memset(degree,0,opts.p*sizeof(double));
	computeDegree(edge,degree,opts.g,opts.wt);

	bv = (double*) R_alloc(opts.p,sizeof(double));
	memset(bv, 0, opts.p*sizeof(double));
	dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, y, &incx, &beta, bv, &incy);
	/* compute the chol factoration */
	Ch_R = (double*) R_alloc(opts.p*opts.p,sizeof(double));
	memset(Ch_R, 0, opts.p*opts.p*sizeof(double));

	/* set Ch_R as rho * eye(p) + rho*T'*T*/
	for(i = 0; i< opts.p; i++)
		Ch_R[i*opts.p+i] = opts.rho*(1 + degree[i]*0.5);

	/* compute Ch_R = (A'*A) + rho * eye(p) + rho*T'*T; */
	beta = 1.0;
	alpha = 1.0;
	dgemm_(&Tra, &Trb,&opts.p, &opts.p, &opts.n, &alpha, A, &opts.n, A, &opts.n, &beta, Ch_R, &opts.p);
	alpha = opts.rho;
	k = 2*opts.g;

	dpotrf_(&UPLO, &opts.p, Ch_R, &opts.p, &info);

	if (info !=0)
	{
/*		free(bv);
		free(degree);
		free(Ch_R);*/
		error("Chol factorization fails");
		return;
	}  

	goscarSub(x,Ch_R,bv,edge,0.5,s1,s2,opts);

/*	free(Ch_R);
	free(bv);
	free(degree);*/
}

void gflassoSub(double *x, double* Ch_R, double* bv, integer* edge, double s1, double s2, Parm opts, ptrdiff_t *sign)
{
	double *x_p,*q, *pk, *v, *u, *Tx,*TTx;
	integer i,iter,nrhs;
	double tv,td,tz;
	double funValp = -1E37;
	double funVal = 0;
	char UPLO = 'U';
	ptrdiff_t info = 0;  
	nrhs = 1;

	x_p = (double*)R_alloc(opts.p,sizeof(double));
	q   = (double*)R_alloc(opts.p,sizeof(double));
	pk  = (double*)R_alloc(opts.g,sizeof(double));
	v   = (double*)R_alloc(opts.g,sizeof(double));
	u   = (double*)R_alloc(opts.p,sizeof(double));
	Tx  = (double*)R_alloc(opts.g,sizeof(double));
	TTx = (double*)R_alloc(opts.p,sizeof(double));

	memset(x_p, 0, opts.p*sizeof(double));
	memset(q,   0, opts.p*sizeof(double));
	memset(pk,  0, opts.g*sizeof(double));
	memset(v,   0, opts.g*sizeof(double));
	memset(u,   0, opts.p*sizeof(double));

	for(iter = 0; iter < opts.aMaxIter; iter ++)
	{
		memset(TTx,0,opts.p*sizeof(double));

		for (i = 0; i< opts.g; i++)
		{
			TTx[edge[2*i]] += opts.wt[i]*(opts.rho*pk[i] - v[i]);
			TTx[edge[2*i+1]] += opts.wt[i]*sign[i]*(opts.rho*pk[i] - v[i]);
		}

		for(i = 0; i< opts.p; i++)
			x[i] = bv[i] - u[i] + opts.rho*q[i] + TTx[i];

		dpotrs_(&UPLO, &opts.p, &nrhs, Ch_R, &opts.p, x, &opts.p, &info);

		for (i = 0; i < opts.p ; i ++)
		{
			td = x[i] + u[i]/opts.rho; 
			if (td < -s1/opts.rho)
				q[i] = td + s1/opts.rho;
			else if(td > s1/opts.rho)
				q[i] = td - s1/opts.rho;
			else
				q[i]= 0;
			u[i] += opts.rho*(x[i] - q[i]);
		}

		for (i = 0; i< opts.g; i++)
			Tx[i] = opts.wt[i] * (x[edge[2*i]] + sign[i]*x[edge[2*i+1]]);

		for (i = 0; i<opts.g; i++)
		{
			td = Tx[i] + v[i]/opts.rho;
			if (td < -s2/opts.rho)
				pk[i] = td + s2/opts.rho;
			else if(td > s2/opts.rho)
				pk[i] = td - s2/opts.rho;
			else
				pk[i]= 0;

			v[i] += opts.rho*(Tx[i] - pk[i]);
		}

		tv = 0; td = 0;
		for (i = 0; i < opts.p; i++)
		{
			tv += pow(x[i] - q[i],2);
			td += pow(x[i] - x_p[i], 2);
		}
		funVal = sqrt(tv);
		/*printf("The Objective value in the %td-th iteration is %f\n", iter+1, (float)funv[iter]);*/
		td = sqrt(td);
		/*printf("The gap between x and x_p in the %td-th iteration is %f\n", iter+1, (float)td); */

		tv = funValp < 1 ? 1 : funValp;
		tz = td<1 ? 1:td;
		if (funVal < opts.aTol *tv && td < opts.aTol *tz)
			break;
		memcpy(x_p,x,opts.p*sizeof(double));
		funValp = funVal;
	}

	memcpy(x,q,opts.p*sizeof(double));
/*	free(x_p);
	free(q);
	free(pk);
	free(v);
	free(u);
	free(Tx);
	free(TTx);*/
}


void gflasso(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts)
{
	integer i,k, nrhs,incx,incy;
	double *bv, *Ch_R, tv, td;
	double alpha = 1.0;
	double beta = 0;

	char trans = 'T';
	char Tra = 'T';
	char Trb = 'N';
	char UPLO = 'U';

	ptrdiff_t *sign;
	ptrdiff_t info = 0;  

	incx = 1;
	incy = 1;
	nrhs = 1;

	sign = (ptrdiff_t*)R_alloc(opts.g,sizeof(ptrdiff_t));
	memset(sign,0,opts.g*sizeof(ptrdiff_t));
	for (i = 0; i<opts.g; i++)
		if (opts.wt[i]<0)
		{ 
			sign[i] = 1;
			opts.wt[i] = -opts.wt[i];
		}
		else if (opts.wt[i]>0)
			sign[i] = -1;
		else
			sign[i] = 0;

	bv = (double*) R_alloc(opts.p,sizeof(double));
	memset(bv, 0, opts.p*sizeof(double));
	dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, y, &incx, &beta, bv, &incy);
		/* compute the chol factoration */
	Ch_R = (double*) R_alloc(opts.p*opts.p,sizeof(double));
	memset(Ch_R, 0, opts.p*opts.p*sizeof(double));

		/* set Ch_R as rho * eye(p) + rho*T'*T*/
	for(i = 0; i< opts.p; i++)
		Ch_R[i*opts.p+i] = opts.rho;

	for (i = 0; i<opts.g;i++)
	{
		td = opts.rho*opts.wt[i]*opts.wt[i];
		tv = td*sign[i]; 
		Ch_R[edge[2*i]*opts.p + edge[2*i+1]] += tv;
		Ch_R[edge[2*i+1]*opts.p + edge[2*i]] += tv;
		Ch_R[edge[2*i]*opts.p + edge[2*i]] += td;
		Ch_R[edge[2*i+1]*opts.p + edge[2*i+1]] +=td;
	}

		/* compute Ch_R = (A'*A) + rho * eye(p) + rho*T'*T; */
	beta = 1.0;
	alpha = 1.0;
	dgemm_(&Tra, &Trb,&opts.p, &opts.p, &opts.n, &alpha, A, &opts.n, A, &opts.n, &beta, Ch_R, &opts.p);
	alpha = opts.rho;
	k = 2*opts.g;

	dpotrf_(&UPLO, &opts.p, Ch_R, &opts.p, &info);
	if (info !=0)
	{
/*		free(bv);
		free(sign);
		free(Ch_R);*/
		error("Chol factorization fails");
		return;
	}  

	gflassoSub(x,Ch_R,bv,edge,s1,s2,opts,sign);

/*	free(Ch_R);
	free(bv);
	free(sign);*/
}


void ncTL(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts)
{
	integer i,k, iter, nrhs,incx,incy;
	double *bv, *Ch_R, *cp, *Ax;
	double ofunVal, funVal, tv, td, f1, tz;
	double alpha = 1.0;
	double beta = 0;

	char trans = 'T';
	char Tra = 'T';
	char Trb = 'N';
	char UPLO = 'U';

	ptrdiff_t *sign;
	ptrdiff_t info = 0;  

	incx = 1;
	incy = 1;
	nrhs = 1;

	ofunVal = -1E20;

	sign = (ptrdiff_t*)R_alloc(opts.g,sizeof(ptrdiff_t));
	memset(sign,0,opts.g*sizeof(ptrdiff_t));
	for (i = 0; i<opts.g; i++)
		if (opts.wt[i]<0)
		{ 
			sign[i] = 1;
			opts.wt[i] = -opts.wt[i];
		}
		else if (opts.wt[i]>0)
			sign[i] = -1;
		else
			sign[i] = 0;

	bv = (double*) R_alloc(opts.p,sizeof(double));
	memset(bv, 0, opts.p*sizeof(double));
	dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, y, &incx, &beta, bv, &incy);
		/* compute the chol factoration */
	Ch_R = (double*) R_alloc(opts.p*opts.p,sizeof(double));
	memset(Ch_R, 0, opts.p*opts.p*sizeof(double));

		/* set Ch_R as rho * eye(p) + rho*T'*T*/
	for(i = 0; i< opts.p; i++)
		Ch_R[i*opts.p+i] = opts.rho;

	for (i = 0; i<opts.g;i++)
	{
		td = opts.rho*opts.wt[i]*opts.wt[i];
		tv = td*sign[i]; 
		Ch_R[edge[2*i]*opts.p + edge[2*i+1]] += tv;
		Ch_R[edge[2*i+1]*opts.p + edge[2*i]] += tv;
		Ch_R[edge[2*i]*opts.p + edge[2*i]] += td;
		Ch_R[edge[2*i+1]*opts.p + edge[2*i+1]] +=td;
	}

		/* compute Ch_R = (A'*A) + rho * eye(p) + rho*T'*T; */
	beta = 1.0;
	alpha = 1.0;
	dgemm_(&Tra, &Trb,&opts.p, &opts.p, &opts.n, &alpha, A, &opts.n, A, &opts.n, &beta, Ch_R, &opts.p);
	alpha = opts.rho;
	k = 2*opts.g;

	dpotrf_(&UPLO, &opts.p, Ch_R, &opts.p, &info);
	if (info !=0)
	{
/*		free(bv);
		free(sign);
		free(Ch_R);*/
		error("Chol factorization fails");
		return;
	}  

	cp = (double*) R_alloc(opts.p,sizeof(double));
	Ax = (double*) R_alloc(opts.n,sizeof(double));

		/* Ax - y*/
	trans = 'N';
	alpha = 1.0;
	beta = 0.0;

	s1 /= opts.tau;

	for (iter = 0; iter < opts.maxIter; iter++)
	{
		for (i = 0; i<opts.p; i++)
			if (x[i]>opts.tau)
				cp[i] = s1 + bv[i];
			else if (x[i]<-opts.tau)
				cp[i] = bv[i] - s1;
            else
                cp[i] = bv[i]; 

		gflassoSub(x,Ch_R,cp,edge,s1,s2,opts,sign);

		dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, x, &incx, &beta, Ax, &incy);

		f1 =0; td = 0;
			/* compute f1 */
		for (i = 0; i<opts.n; i++)
			f1 += (Ax[i]-y[i])*(Ax[i]-y[i]);
	    f1 /= 2.0;

		for (i = 0; i<opts.p; i++)
		{
			tv = fabs(x[i]) < opts.tau ? fabs(x[i]):opts.tau;
			td +=tv;
		}

		f1 += s1*td;

		td = 0;
		for (i = 0; i<opts.g; i++)
		{
			tz = fabs(x[edge[2*i]] + sign[i]*x[edge[2*i+1]]);
			td +=tz*opts.wt[i];
		}
		f1 += s2*td;
		funVal = f1;

		/*printf("%d iteration is %f\n", iter, funVal);*/
		tv = ofunVal < 1 ? 1 : ofunVal;
		if (fabs(funVal - ofunVal) < opts.tol*tv)
			break;
		ofunVal = funVal;
	}

/*	free(Ax);
	free(cp);	
	free(Ch_R);
	free(bv);
	free(sign);*/

}

void ncTF(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts)
{
	integer i,j,k, iter, nrhs,incx,incy;
	double *bv, *Ch_R, *cp, *Ax;
	double ofunVal, funVal,tv,td,f1,tz;
	double alpha = 1.0;
	double beta = 0;

	char trans = 'T';
	char Tra = 'T';
	char Trb = 'N';
	char UPLO = 'U';

	ptrdiff_t *sign;
	ptrdiff_t info = 0;  

	incx = 1;
	incy = 1;
	nrhs = 1;

	ofunVal = -1E20;

	sign = (ptrdiff_t*)R_alloc(opts.g,sizeof(ptrdiff_t));
	memset(sign,0,opts.g*sizeof(ptrdiff_t));
	for (i = 0; i<opts.g; i++)
		if (opts.wt[i]<0)
		{ 
			sign[i] = 1;
			opts.wt[i] = -opts.wt[i];
		}
		else if (opts.wt[i]>0)
			sign[i] = -1;
		else
			sign[i] = 0;

	bv = (double*) R_alloc(opts.p,sizeof(double));
	memset(bv, 0, opts.p*sizeof(double));
	dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, y, &incx, &beta, bv, &incy);
	/* compute the chol factoration */
	Ch_R = (double*) R_alloc(opts.p*opts.p,sizeof(double));
	memset(Ch_R, 0, opts.p*opts.p*sizeof(double));

	/* set Ch_R as rho * eye(p) + rho*T'*T*/
	for(i = 0; i< opts.p; i++)
		Ch_R[i*opts.p+i] = opts.rho;

	for (i = 0; i<opts.g;i++)
	{
		td = opts.rho*opts.wt[i]*opts.wt[i];
		tv = td*sign[i]; 
		Ch_R[edge[2*i]*opts.p + edge[2*i+1]] += tv;
		Ch_R[edge[2*i+1]*opts.p + edge[2*i]] += tv;
		Ch_R[edge[2*i]*opts.p + edge[2*i]] += td;
		Ch_R[edge[2*i+1]*opts.p + edge[2*i+1]] +=td;
	}

		/* compute Ch_R = (A'*A) + rho * eye(p) + rho*T'*T; */
	beta = 1.0;
	alpha = 1.0;
	dgemm_(&Tra, &Trb,&opts.p, &opts.p, &opts.n, &alpha, A, &opts.n, A, &opts.n, &beta, Ch_R, &opts.p);
	alpha = opts.rho;
	k = 2*opts.g;

	dpotrf_(&UPLO, &opts.p, Ch_R, &opts.p, &info);
	if (info !=0)
	{
/*		free(bv);
		free(sign);
		free(Ch_R);*/
		error("Chol factorization fails");
		return;
	}  

	cp = (double*)R_alloc(opts.p,sizeof(double));
	Ax = (double*) R_alloc(opts.n,sizeof(double));

		/* Ax - y*/
	trans = 'N';
	alpha = 1.0;
	beta = 0.0;

	s2 /= opts.tau;

	for (iter = 0; iter < opts.maxIter; iter++)
	{
		memcpy(cp,bv,opts.p*sizeof(double));
		for (i = 0; i<opts.g; i++)
		{
			j = edge[2*i];
			k = edge[2*i+1];

			if (x[j] + sign[i]*x[k]> opts.tau)
			{
				cp[j] += s2*opts.wt[i];
				cp[k] += s2*sign[i]*opts.wt[i];
			}
			else if (x[j] + sign[i]*x[k] <- opts.tau)
			{
				cp[j] -= s2*opts.wt[i];
				cp[k] -= s2*sign[i]*opts.wt[i];
			}
		}

		gflassoSub(x,Ch_R,cp,edge,s1,s2,opts,sign);

		dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, x, &incx, &beta, Ax, &incy);

		f1 =0; td = 0;
		/* compute f1 */
		for (i = 0; i<opts.n; i++)
			f1 += (Ax[i]-y[i])*(Ax[i]-y[i]);
		f1 /= 2.0;

		for (i = 0; i<opts.p; i++)
			td +=fabs(x[i]);

		f1 += s1*td;

		td = 0;
		for (i = 0; i<opts.g; i++)
		{
			tz = fabs(x[edge[2*i]] + sign[i]*x[edge[2*i+1]]);
			if (tz > opts.tau)
			   tz = opts.tau;

			td +=tz*opts.wt[i];
		}
		f1 += s2*td;
		funVal = f1;

		/*printf("%d iteration is %f\n", iter, funVal);*/
		tv = ofunVal < 1 ? 1 : ofunVal;
		if (fabs(funVal - ofunVal) < opts.tol*tv)
			break;
		ofunVal = funVal;
	}

/*	free(Ax);
	free(cp);	
	free(Ch_R);
	free(bv);
	free(sign);*/
}

void ncTLF(double *A, double *y, double *x, integer *edge, double s1, double s2, Parm opts)
{
	integer i,k,j, iter, nrhs,incx,incy;
	double *bv, *Ch_R, *cp, *Ax;
	double ofunVal, funVal, tv, td, f1, tz;
	double alpha = 1.0;
	double beta = 0;

	char trans = 'T';
	char Tra = 'T';
	char Trb = 'N';
	char UPLO = 'U';

	ptrdiff_t *sign;
	ptrdiff_t info = 0;  

	incx = 1;
	incy = 1;
	nrhs = 1;

	ofunVal = -1E20;

	sign = (ptrdiff_t*)R_alloc(opts.g,sizeof(ptrdiff_t));
	memset(sign,0,opts.g*sizeof(ptrdiff_t));
	for (i = 0; i<opts.g; i++)
		if (opts.wt[i]<0)
		{ 
			sign[i] = 1;
			opts.wt[i] = -opts.wt[i];
		}
		else if (opts.wt[i]>0)
			sign[i] = -1;
		else
			sign[i] = 0;

	bv = (double*) R_alloc(opts.p,sizeof(double));
	memset(bv, 0, opts.p*sizeof(double));
	dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, y, &incx, &beta, bv, &incy);
		/* compute the chol factoration */
	Ch_R = (double*) R_alloc(opts.p*opts.p,sizeof(double));
	memset(Ch_R, 0, opts.p*opts.p*sizeof(double));

		/* set Ch_R as rho * eye(p) + rho*T'*T*/
	for(i = 0; i< opts.p; i++)
		Ch_R[i*opts.p+i] = opts.rho;

	for (i = 0; i<opts.g;i++)
	{
		td = opts.rho*opts.wt[i]*opts.wt[i];
		tv = td*sign[i]; 
		Ch_R[edge[2*i]*opts.p + edge[2*i+1]] += tv;
		Ch_R[edge[2*i+1]*opts.p + edge[2*i]] += tv;
		Ch_R[edge[2*i]*opts.p + edge[2*i]] += td;
		Ch_R[edge[2*i+1]*opts.p + edge[2*i+1]] +=td;
	}

		/* compute Ch_R = (A'*A) + rho * eye(p) + rho*T'*T; */
	beta = 1.0;
	alpha = 1.0;
	dgemm_(&Tra, &Trb,&opts.p, &opts.p, &opts.n, &alpha, A, &opts.n, A, &opts.n, &beta, Ch_R, &opts.p);
	alpha = opts.rho;
	k = 2*opts.g;

	dpotrf_(&UPLO, &opts.p, Ch_R, &opts.p, &info);
	if (info !=0)
	{
/*		free(bv);
		free(sign);
		free(Ch_R);*/
		error("Chol factorization fails");
		return;
	}  

	cp = (double*)R_alloc(opts.p,sizeof(double));
	Ax = (double*) R_alloc(opts.n,sizeof(double));

		/* Ax - y*/
	trans = 'N';
	alpha = 1.0;
	beta = 0.0;

	s2 /= opts.tau;
	s1 /= opts.tau;

	for (iter = 0; iter < opts.maxIter; iter++)
	{
		for (i = 0; i<opts.p; i++)
			if (x[i]>opts.tau)
				cp[i] = s1 + bv[i];
			else if (x[i]<-opts.tau)
				cp[i] = bv[i] - s1;
            else
                cp[i] = bv[i]; 

		for (i = 0; i<opts.g; i++)
		{
			j = edge[2*i];
			k = edge[2*i+1];

			if (x[j] + sign[i]*x[k]> opts.tau)
			{
				cp[j] += s2*opts.wt[i];
				cp[k] += s2*sign[i]*opts.wt[i];
			}
			else if (x[j] + sign[i]*x[k] <- opts.tau)
			{
				cp[j] -= s2*opts.wt[i];
				cp[k] -= s2*sign[i]*opts.wt[i];
			}
		}

		gflassoSub(x,Ch_R,cp,edge,s1,s2,opts,sign);

		dgemv_(&trans, &opts.n, &opts.p, &alpha, A, &opts.n, x, &incx, &beta, Ax, &incy);

		f1 =0; td = 0;
		/* compute f1 */
		for (i = 0; i<opts.n; i++)
			f1 += (Ax[i]-y[i])*(Ax[i]-y[i]);
		f1 /= 2.0;

		for (i = 0; i<opts.p; i++)
		{
			tv = fabs(x[i]) < opts.tau ? fabs(x[i]):opts.tau;
			td +=tv;
		}

		f1 += s1*td;

		td = 0;
		for (i = 0; i<opts.g; i++)
		{
			tz = fabs(x[edge[2*i]] + sign[i]*x[edge[2*i+1]]);
			if (tz > opts.tau)
				tz = opts.tau;

			td +=tz*opts.wt[i];
		}
		f1 += s2*td;
		funVal = f1;
		
		/*printf("%d iteration is %f\n", iter, funVal);*/
		tv = ofunVal < 1 ? 1 : ofunVal;
		if (fabs(funVal - ofunVal) < opts.tol*tv)
			break;
		ofunVal = funVal;
	}

/*	free(Ax);
	free(cp);	
	free(Ch_R);
	free(bv);
	free(sign);*/
}

#endif
