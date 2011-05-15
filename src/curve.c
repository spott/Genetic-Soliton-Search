#include <stdio.h>
#include <stdlib.h>
#include "../config.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*  curve.c: Contains all the curve generation rutines.
 */



/* We generate the curve from a curve-gene. 
 * The gene is a two dimensional array that has the location in gene_x
 * and the height in gene_y. As we add more points we refine our curve.
 */

int gene2curve ( double width, gsl_vector * gene_x, gsl_vector * gene_y, gsl_vector * curve_x, gsl_vector * curve_y)
{
    unsigned int size       = gene_x->size;
    unsigned int curve_size = curve_x->size;
    double dx         = width/curve_size;
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, size);
    gsl_spline_init (spline, gene_x->data, gene_y->data , size);
    int i;
    for(i=0; i < curve_x->size; i++)
    {
        gsl_vector_set (curve_x,i,(double)i*dx);
        gsl_vector_set (curve_y,i,gsl_spline_eval (spline,(double)i*dx,acc));
    }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    return 0;
}

int print_curve (gsl_vector * curve_x, gsl_vector * curve_y)
{
    int i;
    for(i=0;i < curve_x->size ;i++)
    {
        printf("%g %g\n",gsl_vector_get (curve_x,i), gsl_vector_get (curve_y,i) );
    }
    return 0;
}
