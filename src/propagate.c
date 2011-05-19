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
#include <gsl/gsl_errno.h>

void propagate_curve(double dx, double tf, gsl_vector * curve_x, gsl_vector_complex * curve_y)
{
    int i;
    int debug = 0;
    gsl_vector_complex * fft_able = gsl_vector_complex_alloc(curve_y->size);
    gsl_vector_complex_memcpy (fft_able,curve_y);

    gsl_fft_complex_wavetable * wavetable = gsl_fft_complex_wavetable_alloc (curve_y->size);
    gsl_fft_complex_workspace * workspace = gsl_fft_complex_workspace_alloc (curve_y->size);

    gsl_fft_complex_forward (fft_able->data, 1, curve_y->size, wavetable, workspace);
    
    // Print out fft coeffs
    for (i = 1900; i < curve_y->size; i++)
    {
        printf ("%e %e\n", (double)i/curve_y->size, GSL_REAL(gsl_vector_complex_get(fft_able,i)));
    }
    printf("\n");
    for (i = 1900; i < curve_y->size; i++)
    {
        printf ("%e %e\n", (double)i/curve_y->size, GSL_IMAG(gsl_vector_complex_get(fft_able,i)));
    }
    
    // Initialize vectors for RK4
    gsl_vector_complex * rk4_1 = gsl_vector_complex_alloc(curve_y->size);
    gsl_vector_complex * rk4_2 = gsl_vector_complex_alloc(curve_y->size);
    gsl_vector_complex * rk4_3 = gsl_vector_complex_alloc(curve_y->size);
    gsl_vector_complex * rk4_4 = gsl_vector_complex_alloc(curve_y->size);

    // Initialize vectors for stepping
    gsl_vector_complex * old_step = gsl_vector_complex_alloc(curve_y->size);
    gsl_vector_complex * new_step = gsl_vector_complex_alloc(curve_y->size);

    //Set starting point for stepping variable.
    double h = pow(10,-5);
    double t = 0.0;
    gsl_complex * II = gsl_complex_rect(0.0,1.0);

    //Set coefficients for GNLSE
    gsl_complex * q1 = gsl_complex_rect(1,0);
    gsl_complex * q2 = gsl_complex_rect(0,0);
    gsl_complex * q3 = gsl_complex_rect(0,0);
    gsl_complex * q4 = gsl_complex_rect(0,0);

    gsl_complex * term1;
    gsl_complex * term2;
    gsl_complex * term3;
    gsl_complex * term4;
    
    gsl_vector_complex_memcp(old_step,curve_y);

    while(t<tf)
    {
        for(i = 0;i<curve_y->size;i++)
        {
            term1 = gsl_complex_mul(II,gsl_complex_mul(q1, gsl_complex_mul(gsl_complex_abs2(gsl_vector_complex_get(old_step,i)) gsl_vector_complex_get(old_step,i))));
            term2 = gsl_complex_mul(II,gsl_complex_mul(q2, gsl_complex_mul(pow(gsl_complex_abs2(gsl_vector_complex_get(old_step,i)),2) gsl_vector_complex_get(old_step,i))));
            gsl_vector_complex_set(rk4_1,i) = gsl_complex_add(term1,gsl_complex_add(term2,gsl_complex_add(term3,term4));
        }
    }



    gsl_vector_complex_free (rk4_1);
    gsl_vector_complex_free (rk4_2);
    gsl_vector_complex_free (rk4_3);
    gsl_vector_complex_free (rk4_4);
    gsl_vector_complex_free (old_step);
    gsl_vector_complex_free (new_step);
    gsl_vector_complex_free (fft_able);
    gsl_fft_complex_wavetable_free (wavetable);
    gsl_fft_complex_workspace_free (workspace);
}

