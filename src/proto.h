/*
 * Decomp
 * 
 * Decompose the optical emissions from disk and jet in 3C 273 with reveberation mapping data.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 * Apr 29, 2019
 */

#include <stdio.h>
#include <stdlib.h>

/* read */
void read_data();
void read_parset();

/* run */
void begin_run();
void end_run();

/* data */
void allocate_memory_data();
void free_memory_data();
void scale_data();
void allocate_memory();
void free_memory();
void init();
double Slc(double tcon, double tline, const void *model);
double Sll(double tcon, double tline, const void *model);
double Src(double tcon, double tline, const void *model);
double Srr(double tcon, double tline, const void *model);
void set_covar_Pmat_data(const void *model);
void set_covar_Umat_rec(const void *model);
void set_covar_Amat_rec(const void *model);
void set_covar_PSmat_data(const void *model);
void set_covar_Umat_conjd(const void *model);
void set_covar_Amat_conjd(const void *model);
void reconstruct_all(const void *model);
void reconstruct_all2(const void *model);
void reconstruct_conjd(const void *model);

/* mcmc */
void set_par_range();
double prob(const void *model);
double prob2(const void *model);
double prob3(const void *model);
void recon();
void recon_init();
void recon_end();
double mc_dnest(int argc, char **argv);
void recon_postprocess();

/* time */
double second();
double timediff(double t0, double t1);
void get_hms(double dt, int *h, int *m, double *s);

/* help */
void print_help();

/* matrix operations */
void inverse_mat(double *a, int n, int *info);
void inverse_symat_lndet(double * a, int n, double *lndet, int *info, int *ipiv);
void inverse_symat_lndet_sign(double * a, int n, double *lndet, int *info, int *sign, int *ipiv);
double det_mat(double *a, int n, int *info);
double lndet_mat(double *a, int n, int *info);
double lndet_mat2(double *a, int n, int *info, int *sign);
double lndet_mat3(double *a, int n, int *info, int *sign);
void display_mat(double *a, int m, int n);
void multiply_mat(double * a, double *b, double *c, int n);
void multiply_mat_transposeA(double * a, double *b, double *c, int n);
void multiply_mat_transposeB(double * a, double *b, double *c, int n);
void multiply_matvec(double *a, double *x, int n, double *y);
void multiply_matvec_transposeA(double *a, double *x, int n, double *y);
void multiply_matvec_MN(double * a, int m, int n, double *x, double *y);
void multiply_mat_MN(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_transposeA(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_transposeB(double * a, double *b, double *c, int m, int n, int k);
int multiply_mat_MN_inverseA(double * a, double *b, int m, int n, int *ipiv);
int multiply_symat_MN_inverseA(double * a, double *b, int m, int n, int *ipiv);
int multiply_pomat_MN_inverseA(double * a, double *b, int m, int n, int *ipiv);
void multiply_vec2mat(double * x, double * a, int n);
void eigen_sym_mat(double *a, int n, double *val, int *info);
void Chol_decomp_U(double *a, int n, int *info);
void Chol_decomp_L(double *a, int n, int *info);
double ** matrix_malloc(int n1, int n2);
double * array_malloc(int n);
void test_mathfun();
void compute_semiseparable_drw(double *t, int n, double a1, double c1, double *sigma, double syserr,  double *W, double *D, double *phi);
void multiply_matvec_semiseparable_drw(double *y, double  *W, double *D, double *phi, int n, double a1, double *z);
void multiply_mat_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z);
void multiply_mat_transposeB_semiseparable_drw(double *Y, double  *W, double *D, double *phi, int n, int m, double a1, double *Z);