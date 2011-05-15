/* This contains the main line of execution for the program.
 */

#include <stdio.h>
#include "../config.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>

#define MAX_GENE_SIZE 100
#define CURVE_SIZE 1000


int main (void)
{
    double width = 1.0;
    double dx    = width/CURVE_SIZE;
    int i;
    unsigned int gene_size = 0;
    gsl_vector         * gene_x  = gsl_vector_alloc (MAX_GENE_SIZE);
    gsl_vector_complex * gene_y  = gsl_vector_complex_alloc (MAX_GENE_SIZE);
    gsl_vector         * curve_x = gsl_vector_alloc (CURVE_SIZE);
    gsl_vector_complex * curve_y = gsl_vector_complex_alloc (CURVE_SIZE);
    gsl_vector_set_zero (gene_x);
    gsl_vector_complex_set_zero (gene_y);
    gsl_vector_set_all (curve_x,1.0);
    gsl_vector_complex_set_zero (curve_y);

    //Fill curve_x with descrete steps
    for(i=0;i<curve_x->size;i++)
        gsl_vector_set (curve_x,i,(double)i*dx);

    // Set Pinning on LHS
    int gene_index = -1;
    gene_index++;
    gsl_vector_set (gene_x, gene_index,0);
    gsl_vector_complex_set (gene_y, gene_index,gsl_complex_rect(0.0,0.0));
    gene_index++;
    gsl_vector_set (gene_x, gene_index,dx);
    gsl_vector_complex_set (gene_y, gene_index,gsl_complex_rect(0.0,0.0));
    // Set points
    gene_index++;
    gsl_vector_set (gene_x, gene_index,0.5);
    gsl_vector_complex_set (gene_y, gene_index,gsl_complex_rect(1.0,0.7));
    // Set Pinning on LHS
    gene_index++;
    gsl_vector_set (gene_x, gene_index,width-dx);
    gsl_vector_complex_set (gene_y, gene_index,gsl_complex_rect(0.0,0.0));
    gene_index++;
    gsl_vector_set (gene_x, gene_index,width);
    gsl_vector_complex_set (gene_y, gene_index,gsl_complex_rect(0.0,0.0));
    
    //printf("gene_index=%i\n",gene_index);
    //for(i=0;i<gene_index;i++)
    //    printf("%g %g+%gi\n",gsl_vector_get (gene_x,i), GSL_REAL(gsl_vector_complex_get(gene_y,i)), GSL_IMAG(gsl_vector_complex_get(gene_y,i)));
    

    // Convert to curve
    gene2curve(dx,gene_index+1,gene_x,gene_y,curve_x,curve_y);

    print_curve(curve_x,curve_y);

    return 0;
}
