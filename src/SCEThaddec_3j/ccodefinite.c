#include <stdio.h>
#include "wrpfinite.h"

double hplevalfinite_(double*yy, double* zz, double* lm){

	double y = *yy;
	double z = *zz;
	double logm = *lm;

//	printf("y = %f \n",y);
//	printf("z = %f \n",z);

	return wrpfinite(y,z,logm);
}
