#include <stdio.h>
#include <stdlib.h>
#include "../config.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>


/*  curve.c: Contains all the curve chromosomeration rutines.
 */




int chromosome2curve ( double dx,
                       unsigned int chromosome_size,
                       gsl_vector * chromosome_x,
                       gsl_vector_complex * chromosome_y,
                       gsl_vector * curve_x,
                       gsl_vector_complex * curve_y )
{
    printf("#Entering chromosome2curve\n");
    /* We chromosomerate the curve from a curve-chromosome. 
     * The chromosome is a two dimensional array that has the location in chromosome_x
     * and the height in chromosome_y. As we add more points we refine our curve.
     */
    unsigned int curve_size = curve_x->size;
    int i;

    gsl_vector_complex_set_zero (curve_y);

    gsl_interp_accel *acc      = gsl_interp_accel_alloc ();
    gsl_spline *spline_real    = gsl_spline_alloc (gsl_interp_cspline, chromosome_size);
    gsl_spline *spline_imag    = gsl_spline_alloc (gsl_interp_cspline, chromosome_size);
    //Let us do the real parts first
    printf("#chromosome2curve: Finding Spline on real part\n");
    gsl_vector * chromosome_y_temp = gsl_vector_alloc (chromosome_size);
    for(i=0; i < chromosome_size;i++)
        gsl_vector_set (chromosome_y_temp,i,GSL_REAL(gsl_vector_complex_get(chromosome_y,i)));

    gsl_spline_init (spline_real, chromosome_x->data, chromosome_y_temp->data , chromosome_size);

    //Now for the imag parts
    printf("#chromosome2curve: Finding Spline on imag part\n");
    for(i=0; i < chromosome_size;i++)
        gsl_vector_set (chromosome_y_temp,i,GSL_IMAG(gsl_vector_complex_get(chromosome_y,i)));

    gsl_spline_init (spline_imag, chromosome_x->data, chromosome_y_temp->data , chromosome_size);
    
    printf("#chromosome2curve: Writing Spline data to curve\n");
    for(i=0; i < curve_x->size; i++)
    {
        gsl_vector_complex_set (curve_y,i,
            gsl_complex_rect (
                gsl_spline_eval (spline_real,(double)i*dx,acc),
                gsl_spline_eval (spline_imag,(double)i*dx,acc)
            )
        );
    }
    printf("#chromosome2curve: Freeing allocated workspace\n");
    gsl_spline_free (spline_real);
    gsl_spline_free (spline_imag);
    gsl_interp_accel_free (acc);
    gsl_vector_free (chromosome_y_temp);

    return 0;
}

void print_curve (gsl_vector * curve_x, gsl_vector_complex * curve_y)
{
    printf("#Entering print_curve\n");
    int i;
    for(i=0;i < curve_x->size ;i++)
    {
        printf("%g %g\n",gsl_vector_get (curve_x,i), GSL_REAL(gsl_vector_complex_get (curve_y,i)) );
    }
    printf("\n");
    for(i=0;i < curve_x->size ;i++)
    {
        printf("%g %g\n",gsl_vector_get (curve_x,i), GSL_IMAG(gsl_vector_complex_get (curve_y,i)) );
    }
}

void curve2arena(gsl_vector_complex * arena , gsl_vector_complex * curve)
{
    int i;
    gsl_vector_complex_set_zero (arena);
    int center = floor ((arena->size)/2);
    for (i=0; i< curve->size;i++)
        gsl_vector_complex_set(arena,center-floor((curve->size)/2)+i,gsl_vector_complex_get(curve,i));
}
