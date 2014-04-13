#include "fgsg.h"


void do_gflasso(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol);
void do_goscar(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol);
void do_ncFGS(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol);
void do_ncTF(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol);
void do_ncTFGS(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol);
void do_ncTL(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol);
void do_ncTLF(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol);

void do_gflasso(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol){
	int i;
	Parm opts;
	opts.n=opts_n[0];
	opts.p=opts_p[0];
	opts.g=opts_g[0];
	opts.maxIter=(integer)RmaxIter[0];
	opts.aMaxIter=(integer)RaMaxIter[0];
	opts.rho=Rrho[0];
	opts.tau=Rtau[0];
	opts.tol=Rtol[0];
	opts.aTol=Ratol[0];
	opts.wt=Rwt;

	integer* edge = (integer*)R_alloc(opts.g,sizeof(integer));
	for (i = 0; i<opts.g; i++){
		edge[i] = (integer)tp[i] - 1;
	}
	opts.g=opts.g/2;
	gflasso(A,y,x,edge,s1[0],s2[0],opts);

    //free(edge);
}

void do_goscar(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol){
	int i;
	Parm opts;
	opts.n=opts_n[0];
	opts.p=opts_p[0];
	opts.g=opts_g[0];
	opts.maxIter=(integer)RmaxIter[0];
	opts.aMaxIter=(integer)RaMaxIter[0];
	opts.rho=Rrho[0];
	opts.tau=Rtau[0];
	opts.tol=Rtol[0];
	opts.aTol=Ratol[0];
	opts.wt=Rwt;
	for (i = 0; i<opts.g; i++){
		opts.wt[i] = fabs(opts.wt[i]);
	}

	integer* edge = (integer*)R_alloc(opts.g,sizeof(integer));
	for (i = 0; i<opts.g; i++){
		edge[i] = (integer)tp[i] - 1;
	}
	opts.g=opts.g/2;
	goscar(A,y,x,edge,s1[0],s2[0],opts);

    //free(edge);
}

void do_ncFGS(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol){
	int i;
	Parm opts;
	opts.n=opts_n[0];
	opts.p=opts_p[0];
	opts.g=opts_g[0];
	opts.maxIter=(integer)RmaxIter[0];
	opts.aMaxIter=(integer)RaMaxIter[0];
	opts.rho=Rrho[0];
	opts.tau=Rtau[0];
	opts.tol=Rtol[0];
	opts.aTol=Ratol[0];
	opts.wt=Rwt;
	for (i = 0; i<opts.g; i++){
		opts.wt[i] = fabs(opts.wt[i]);
	}
	integer* edge = (integer*)R_alloc(opts.g,sizeof(integer));

	for (i = 0; i<opts.g; i++){
		edge[i] = (integer)tp[i] - 1;
	}
	opts.g=opts.g/2;

	ncFGS(A,y,x,edge,s1[0],s2[0],opts);

    //free(edge);
}

void do_ncTF(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol){
	int i;
	Parm opts;
	opts.n=opts_n[0];
	opts.p=opts_p[0];
	opts.g=opts_g[0];
	opts.maxIter=(integer)RmaxIter[0];
	opts.aMaxIter=(integer)RaMaxIter[0];
	opts.rho=Rrho[0];
	opts.tau=Rtau[0];
	opts.tol=Rtol[0];
	opts.aTol=Ratol[0];
	opts.wt=Rwt;
	for (i = 0; i<opts.g; i++){
		opts.wt[i] = fabs(opts.wt[i]);
	}

	integer* edge = (integer*)R_alloc(opts.g,sizeof(integer));
	for (i = 0; i<opts.g; i++){
		edge[i] = (integer)tp[i] - 1;
	}
	opts.g=opts.g/2;

	ncTF(A,y,x,edge,s1[0],s2[0],opts);

   // free(edge);
}

void do_ncTFGS(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol){
	int i;
	Parm opts;
	opts.n=opts_n[0];
	opts.p=opts_p[0];
	opts.g=opts_g[0];
	opts.maxIter=(integer)RmaxIter[0];
	opts.aMaxIter=(integer)RaMaxIter[0];
	opts.rho=Rrho[0];
	opts.tau=Rtau[0];
	opts.tol=Rtol[0];
	opts.aTol=Ratol[0];
	opts.wt=Rwt;
	for (i = 0; i<opts.g; i++){
		opts.wt[i] = fabs(opts.wt[i]);
	}

	integer* edge = (integer*)R_alloc(opts.g,sizeof(integer));
	for (i = 0; i<opts.g; i++){
		edge[i] = (integer)tp[i] - 1;
	}
	opts.g=opts.g/2;

	ncTFGS(A,y,x,edge,s1[0],s2[0],opts);

    //free(edge);
}

void do_ncTL(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol){
	int i;
	Parm opts;
	opts.n=opts_n[0];
	opts.p=opts_p[0];
	opts.g=opts_g[0];
	opts.maxIter=(integer)RmaxIter[0];
	opts.aMaxIter=(integer)RaMaxIter[0];
	opts.rho=Rrho[0];
	opts.tau=Rtau[0];
	opts.tol=Rtol[0];
	opts.aTol=Ratol[0];
	opts.wt=Rwt;

	integer* edge = (integer*)R_alloc(opts.g,sizeof(integer));
	for (i = 0; i<opts.g; i++){
		edge[i] = (integer)tp[i] - 1;
	}
	opts.g=opts.g/2;

	ncTL(A,y,x,edge,s1[0],s2[0],opts);

   // free(edge);
}

void do_ncTLF(double* x, double* A, double* y, int* tp, double* Rwt, double* s1, double* s2, int* opts_n,
				int* opts_p, int* opts_g, int* RmaxIter, int* RaMaxIter, int* Rrho, int* Rtau, int* Rtol, int* Ratol){
	int i;
	Parm opts;
	opts.n=opts_n[0];
	opts.p=opts_p[0];
	opts.g=opts_g[0];
	opts.maxIter=(integer)RmaxIter[0];
	opts.aMaxIter=(integer)RaMaxIter[0];
	opts.rho=Rrho[0];
	opts.tau=Rtau[0];
	opts.tol=Rtol[0];
	opts.aTol=Ratol[0];
	opts.wt=Rwt;

	integer* edge = (integer*)R_alloc(opts.g,sizeof(integer));
	for (i = 0; i<opts.g; i++){
		edge[i] = (integer)tp[i] - 1;
	}
	opts.g=opts.g/2;

	ncTLF(A,y,x,edge,s1[0],s2[0],opts);

    //free(edge);
}
