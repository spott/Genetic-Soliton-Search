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


/*  curve.c: Contains all the curve generation rutines.
 */




int gene2curve ( double dx, unsigned int gene_size, gsl_vector * gene_x, gsl_vector_complex * gene_y, gsl_vector * curve_x, gsl_vector_complex * curve_y)
{
    printf("#Entering gene2curve\n");
    /* We generate the curve from a curve-gene. 
     * The gene is a two dimensional array that has the location in gene_x
     * and the height in gene_y. As we add more points we refine our curve.
     */
    unsigned int curve_size = curve_x->size;
    int i;
    gsl_interp_accel *acc      = gsl_interp_accel_alloc ();
    gsl_spline *spline_real    = gsl_spline_alloc (gsl_interp_cspline, gene_size);
    gsl_spline *spline_imag    = gsl_spline_alloc (gsl_interp_cspline, gene_size);
    //Let us do the real parts first
    printf("#gene2curve: Finding Spline on real part\n");
    gsl_vector * gene_y_temp = gsl_vector_alloc (gene_size);
    for(i=0; i < gene_size;i++)
        gsl_vector_set (gene_y_temp,i,GSL_REAL(gsl_vector_complex_get(gene_y,i)));

    gsl_spline_init (spline_real, gene_x->data, gene_y_temp->data , gene_size);

    //Now for the imag parts
    printf("#gene2curve: Finding Spline on imag part\n");
    for(i=0; i < gene_size;i++)
        gsl_vector_set (gene_y_temp,i,GSL_IMAG(gsl_vector_complex_get(gene_y,i)));

    gsl_spline_init (spline_imag, gene_x->data, gene_y_temp->data , gene_size);
    
    printf("#gene2curve: Writing Spline data to curve\n");
    for(i=0; i < curve_x->size; i++)
    {
        gsl_vector_complex_set (curve_y,i,
            gsl_complex_rect (
                gsl_spline_eval (spline_real,(double)i*dx,acc),
                gsl_spline_eval (spline_imag,(double)i*dx,acc)
            )
        );
    }
    printf("#gene2curve: Freeing allocated workspace\n");
    gsl_spline_free (spline_real);
    gsl_spline_free (spline_imag);
    gsl_interp_accel_free (acc);

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
