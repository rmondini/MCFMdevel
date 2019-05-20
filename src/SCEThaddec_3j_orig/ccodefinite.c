#include <stdio.h>
#include "wrpfinite.h"

double hplevalfinite_(double*yy, double* zz, double* lm, int* nn){

	double y = *yy;
	double z = *zz;
	double logm = *lm;
	int n = *nn;

//	printf("y = %f \n",y);
//	printf("z = %f \n",z);
//	printf("n = %i \n",n);

	return wrpfinite(y,z,logm,n);
}
