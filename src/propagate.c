/* This file contains the code to propagate the wave. This is also 
 * where we set the differential equation coefficients.
 */

#include <stdio.h>
#include "../config.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_fft_complex.h>

void propagate_curve(double dx, double dt, gsl_vector * curve_x, gsl_vector_complex * curve_y)
{
    int i;
    int debug = 0;
    gsl_vector_complex * fft_able = gsl_vector_complex_alloc(curve_y->size);
    gsl_vector_complex_memcpy (fft_able,curve_y);

    gsl_fft_complex_wavetable * wavetable = gsl_fft_complex_wavetable_alloc (curve_y->size);
    gsl_fft_complex_workspace * workspace = gsl_fft_complex_workspace_alloc (curve_y->size);

    gsl_fft_complex_forward (fft_able->data, 1, curve_y->size, wavetable, workspace);
    
    // Print out fft coeffs
    //for (i = 1015; i < curve_y->size; i++)
    //{
    //       printf ("%d %e\n", i, GSL_REAL(gsl_vector_complex_get(fft_able,i)));
    //}
    //printf("\n");
    //for (i = 1015; i < curve_y->size; i++)
    //{
    //       printf ("%d %e\n", i, GSL_IMAG(gsl_vector_complex_get(fft_able,i)));
    //}

    gsl_fft_complex_wavetable_free (wavetable);
    gsl_fft_complex_workspace_free (workspace);
}

