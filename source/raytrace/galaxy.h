class Galaxy {

 public:

    double **F_all;
    double **C_all;
    double *D_all;
    double **F_all_2d;
    double **C_all_2d;
    double *D_all_2d;

    int sersic(double a,double b,double c,double alpha,double beta,
               double n,double *x_out,double *y_out);

    double sample_sersic(char*);

    int sersic2d(double a,double b,double beta,double n,
                 double *x_out,double *y_out);

    double sample_sersic_2d(char*);

};
