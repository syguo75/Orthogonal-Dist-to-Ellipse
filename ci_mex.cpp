#include "mex.h"
#include <math.h>

#ifdef __cplusplus
#define EXTERN_C extern "C"
#else
#define EXTERN_C
#endif


inline double sign(double v)
{
	return v > 0 ? 1 : (v < 0 ? -1:0);
}

#define xIn prhs[0]
#define yIn prhs[1]
#define aIn prhs[2]
#define bIn prhs[3]
#define tolIn prhs[4]
#define sOut plhs[0]

EXTERN_C void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double a = mxGetScalar(aIn);
	double b = mxGetScalar(bIn);
	double d = a*a - b*b;

	double tol = mxGetScalar(tolIn);
	double stol = tol / a;

	int n = mxGetNumberOfElements(xIn);

	sOut = mxCreateDoubleMatrix(n,1,mxREAL);
	
	const double *x = mxGetPr(xIn);
	const double *y = mxGetPr(yIn);
	double *s = mxGetPr(sOut);
	
	for (int i = 0; i < n; i++) {
		if (x[i]==0) {
			s[i] = 1;
			continue;
		}

		if (y[i]==0) {
			if (x[i] >= d/a)
				s[i] = 0;
			else
				s[i] = sqrt(1 - (a*a*x[i]*x[i])/(d*d));
			continue;
		}

		s[i] = 1;
		
		while (true) {
			double a1 = d*d;
			double a2 = b*d*y[i];
			double a3 = a*a*x[i]*x[i] + b*b*y[i]*y[i] - d*d;
			double a4 = b*b*y[i]*y[i];
			
			double num = ((3*a1*s[i] + 4*a2)*s[i] + a3)*s[i]*s[i] + a4;
			double den = ((4*a1*s[i] + 6*a2)*s[i] + 2*a3)*s[i] - 2*a2;

			if (den==0)
				break;

			double ns = num / den;
			double delta = abs(ns-s[i]);
			s[i] = ns;

			if (delta<stol)
				break;
		}
	}
}
