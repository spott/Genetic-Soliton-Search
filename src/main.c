/* This contains the main line of execution for the program.
 */

#include <stdio.h>
#include "../config.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define MAX_GENE_SIZE 5
#define CURVE_SIZE 1000


int main (void)
{
    double width = 1.0;
    double dx    = width/CURVE_SIZE;
    int i;
    gsl_vector * gene_x  = gsl_vector_alloc (MAX_GENE_SIZE);
    gsl_vector * gene_y  = gsl_vector_alloc (MAX_GENE_SIZE);
    gsl_vector * curve_x = gsl_vector_alloc (CURVE_SIZE);
    gsl_vector * curve_y = gsl_vector_alloc (CURVE_SIZE);
    for(i=0; i < MAX_GENE_SIZE; i++)
    {
        gsl_vector_set (gene_x, i,1);
        gsl_vector_set (gene_y, i,0);
    }
    for(i=0; i < CURVE_SIZE; i++)
    {
        gsl_vector_set (curve_x, i,0);
        gsl_vector_set (curve_y, i,0);
    }
    // Set Pinning on LHS
    gsl_vector_set (gene_x, 0,0);
    gsl_vector_set (gene_y, 0,0);
    gsl_vector_set (gene_x, 1,dx);
    gsl_vector_set (gene_y, 1,0);
    // Set points
    gsl_vector_set (gene_x, 2,0.5);
    gsl_vector_set (gene_y, 2,1.0);
    // Set Pinning on LHS
    gsl_vector_set (gene_x, 3,width-dx);
    gsl_vector_set (gene_y, 3,0);
    gsl_vector_set (gene_x, 4,width);
    gsl_vector_set (gene_y, 4,0);

    /*for(i=0;i<6;i++)
        printf("%g %g\n",gsl_vector_get (gene_x,i), gsl_vector_get(gene_y,i));
    */

    // Convert to curve
    gene2curve(1.0,gene_x,gene_y,curve_x,curve_y);

    print_curve(curve_x,curve_y);
    return 0;
}
