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
#define xeOut plhs[0]
#define yeOut plhs[1]

#define MAX_ITER 50
double NaN = mxGetNaN();

EXTERN_C void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double a = mxGetScalar(aIn);
	double b = mxGetScalar(bIn);
	double a_b = a*a - b*b;

	double tol = mxGetScalar(tolIn);

	int n = mxGetNumberOfElements(xIn);

	xeOut = mxCreateDoubleMatrix(n,1,mxREAL);
	yeOut = mxCreateDoubleMatrix(n,1,mxREAL);

	const double *x = mxGetPr(xIn);
	const double *y = mxGetPr(yIn);
	double *xe = mxGetPr(xeOut);
	double *ye = mxGetPr(yeOut);

	for (int i = 0; i < n; i++) {
		double k1 = (a*b) / sqrt(b*b*x[i]*x[i] + a*a*y[i]*y[i]);
		double xk1 = x[i]*k1;
		double yk1 = y[i]*k1;

		double xk2, yk2;
		if (abs(x[i] < a)) {
			xk2 = x[i];
			yk2 = sign(y[i]) * (b*sqrt(a*a - x[i]*x[i])/a);
		} else {
			xk2 = sign(x[i]) * a;
			yk2 = 0;
		}

		xe[i] = (xk1+xk2) / 2;
		ye[i] = (yk1+yk2) / 2;

		int nIter = 0;

		while (true) {
			double q1 = b*b*xe[i];
			double q2 = a*a*ye[i];
			double q3 = a_b*ye[i] + b*b*y[i];
			double q4 = a_b*xe[i] - a*a*x[i];
			double f1 = (a*a*ye[i]*ye[i] + b*b*xe[i]*xe[i] - a*a*b*b) / 2;
			double f2 = b*b*xe[i]*(y[i]-ye[i]) - a*a*ye[i]*(x[i]-xe[i]);

			double dx = (f2*q2 - f1*q4) / (q1*q4 - q2*q3);
			double dy = (f1*q3 - f2*q1) / (q1*q4 - q2*q3);

			xe[i] += dx;
			ye[i] += dy;

			if (abs(dx)<tol && abs(dy)<tol)
				break;

			nIter++;
			if (nIter == MAX_ITER) {
				xe[i] = NaN;
				ye[i] = NaN;
				break;
			}

		}
	}
}
