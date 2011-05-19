/* This contains the main line of execution for the program.
 */

#include <stdio.h>
#include "../config.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>

#include "curve.h"
#include "propagate.h"

#define MAX_GENE_SIZE 100
#define CURVE_SIZE 2048


int main (void)
{
    double width = 2*M_PI;
    double dx    = width/CURVE_SIZE;
    double dt    = pow(10,-3);
    int i;
    unsigned int chromosome_size = 0;
    gsl_vector         * chromosome_x  = gsl_vector_alloc (MAX_GENE_SIZE);
    gsl_vector_complex * chromosome_y  = gsl_vector_complex_alloc (MAX_GENE_SIZE);
    gsl_vector         * curve_x = gsl_vector_alloc (CURVE_SIZE);
    gsl_vector_complex * curve_y = gsl_vector_complex_alloc (CURVE_SIZE);
    gsl_vector_set_zero (chromosome_x);
    gsl_vector_complex_set_zero (chromosome_y);
    gsl_vector_set_all (curve_x,1.0);
    gsl_vector_complex_set_zero (curve_y);

    //Fill curve_x with descrete steps
    for(i=0;i<curve_x->size;i++)
        gsl_vector_set (curve_x,i,(-width/2)+(double)i*dx);

    // Set Pinning on LHS.
    // We only plot points x \in [-pi/4,pi/4].
    // This way we can propagate the cure.
    int chromosome_index = -1;
    chromosome_index++;
    gsl_vector_set (chromosome_x, chromosome_index,-M_PI_4);
    gsl_vector_complex_set (chromosome_y, chromosome_index,gsl_complex_rect(0.0,0.0));
    chromosome_index++;
    gsl_vector_set (chromosome_x, chromosome_index,-M_PI_4+dx);
    gsl_vector_complex_set (chromosome_y, chromosome_index,gsl_complex_rect(0.0,0.0));
    // Set points
    chromosome_index++;
    gsl_vector_set (chromosome_x, chromosome_index,-M_PI_4/2);
    gsl_vector_complex_set (chromosome_y, chromosome_index,gsl_complex_rect(-1.0,0));
    chromosome_index++;
    gsl_vector_set (chromosome_x, chromosome_index,0.0);
    gsl_vector_complex_set (chromosome_y, chromosome_index,gsl_complex_rect(1.0,0.7));
    chromosome_index++;
    gsl_vector_set (chromosome_x, chromosome_index,M_PI_4/2);
    gsl_vector_complex_set (chromosome_y, chromosome_index,gsl_complex_rect(-1.0,0));
    // Set Pinning on LHS
    chromosome_index++;
    gsl_vector_set (chromosome_x, chromosome_index,M_PI_4-dx);
    gsl_vector_complex_set (chromosome_y, chromosome_index,gsl_complex_rect(0.0,0.0));
    chromosome_index++;
    gsl_vector_set (chromosome_x, chromosome_index,M_PI_4);
    gsl_vector_complex_set (chromosome_y, chromosome_index,gsl_complex_rect(0.0,0.0));
    
    //printf("chromosome_index=%i\n",chromosome_index);
    //for(i=0;i<chromosome_index;i++)
    //    printf("%g %g+%gi\n",gsl_vector_get (chromosome_x,i), GSL_REAL(gsl_vector_complex_get(chromosome_y,i)), GSL_IMAG(gsl_vector_complex_get(chromosome_y,i)));
    

    // Convert to curve
    chromosome2curve(dx,chromosome_index+1,chromosome_x,chromosome_y,curve_x,curve_y);

    print_curve(curve_x,curve_y);

    //propagate_curve(dx,dt,curve_x,curve_y);


    //Clear out used memory.
    gsl_vector_free(chromosome_x);
    gsl_vector_complex_free(chromosome_y);
    gsl_vector_free(curve_x);
    gsl_vector_complex_free(curve_y);
    

    return 0;
}
