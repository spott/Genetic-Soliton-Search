int chromosome2curve ( double dx,
                       unsigned int chromosome_size,
                       gsl_vector * chromosome_x,
                       gsl_vector_complex * chromosome_y,
                       gsl_vector * curve_x,
                       gsl_vector_complex * curve_y );
void print_curve (gsl_vector * curve_x, gsl_vector_complex * curve_y);
void curve2arena(gsl_vector_complex * arena , gsl_vector_complex * curve);
