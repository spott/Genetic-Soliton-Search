#include <stdio.h>
#include "../config.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>

/*  This file contains all the curve generation rutines.
 */ 


int curve_print(unsigned int length, double dx, double *ap)
{
    int i;
    for(i=0;i<length;i++)
    {
        printf("%g %g\n",(double)i*dx,*(ap+i));
    }
    return 0;
}
