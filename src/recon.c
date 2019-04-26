#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_sf_erf.h>

#include "dnestvars.h"
#include "allvars.h"
#include "proto.h"

#include "mc_dnest.h"

void recon()
{
  int i, argc=0;
  char **argv;
  double logz;

  // configure restart of dnest 
  argv = malloc(9*sizeof(char *));
  for(i=0; i<9; i++)
  {
    argv[i] = malloc(BRAINS_MAX_STR_LENGTH*sizeof(char));
  }
  //setup argc and argv
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc], parset.file_dir);
  strcat(argv[argc++], "/data/restart_dnest.txt");

  if(parset.flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restart_dnest.txt");
  }
  if(parset.flag_postprc == 1)
  {
    strcpy(argv[argc++], "-p");
  }
  if(parset.flag_temp == 1)
  {
    sprintf(argv[argc++], "-t%f", parset.temperature);
  }

  strcpy(argv[argc++], "-l");

  recon_init();

  logz = mc_dnest(argc, argv);

  recon_end();
  
  //clear up argv 
  for(i=0; i<9; i++)
  {
    free(argv[i]);
  }
  free(argv);
}

void recon_init()
{
  /* option file */
  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONS");
}

void recon_end()
{
  int i;
  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
  }
  free(par_range_model);
}

/* 
 *  likelehood probability 
 *  note that matrix operation A^-1 x B is implemented by calling "multiply_mat_MN_inverseA()".
 */
double prob(const void *model)
{
  int i, nq=3, info, sign;
  double prob=0.0, lndet, lndet_ICq;
  double *ICq, *yq, *ybuf, *y, *yave;

  ICq = workspace;
  yq = ICq + nq*nq;
  yave = yq + nq;
  y = yave + n_all_data;
  ybuf = y + n_all_data;

  set_covar_Pmat_data(model);

  memcpy(Tmat1, PCmat_data, n_all_data*n_all_data*sizeof(double));
  memcpy(Tmat2, Larr_data, n_all_data*nq*sizeof(double));

  multiply_mat_MN_inverseA(Tmat1, Tmat2, n_all_data, nq); // Tmat2 = C^-1 * L;  NxNq
  
  multiply_mat_MN_transposeA(Larr_data, Tmat2, ICq, nq, nq, n_all_data); // ICq = L^T*C^-1*L; NqxNq
  multiply_mat_MN_transposeA(Tmat2, Fall_data, yq, nq, 1, nq); // yq = L^T*C^-1*y;  Nqx1
  memcpy(Tmat1, ICq, nq*nq*sizeof(double));
  multiply_mat_MN_inverseA(Tmat1, yq, nq, 1); // yq = (L^T*C^-1*L)^-1 * L^T*C^-1*y; Nqx1

  multiply_mat_MN(Larr_data, yq, yave, n_all_data, 1, nq); // yave = L * q; Nx1

  for(i=0; i<n_all_data; i++)y[i] = Fall_data[i] - yave[i];
  memcpy(Tmat1, PCmat_data, n_all_data*n_all_data*sizeof(double));
  memcpy(ybuf, y, n_all_data*sizeof(double));
  multiply_mat_MN_inverseA(Tmat1, ybuf, n_all_data, 1); // ybuf = C^-1 * y; Nx1

  prob = -0.5*cblas_ddot(n_all_data, y, 1, ybuf, 1); // y^T * C^-1 * y
  if(prob > 0.0 )  // check if prob is positive
  {
    prob = -DBL_MAX;
    printf("prob >0!\n");
    return prob;
  }

  lndet = lndet_mat3(PCmat_data, n_all_data, &info, &sign);
  if(info!=0|| sign==-1)
  {
    prob = -DBL_MAX;
    printf("lndet_C %f %d!\n", lndet, sign);
    return prob;
  }

  lndet_ICq = lndet_mat3(ICq, nq, &info, &sign);
  if(info!=0 || sign==-1 )
  {
    prob = -DBL_MAX;
    printf("lndet_ICq %f %d!\n", lndet_ICq, sign);
    return prob;
  }

  prob += -0.5*lndet - 0.5*lndet_ICq;
  
  return prob;
}


void set_covar_Pmat_data(const void *model)
{
  int i, j, np;
  double *pm = (double *)model;
  double syserr_con, syserr_line, syserr_radio, sig_d, tau_d, sig_j, tau_j;
  double t1, t2, error;

  syserr_con = (exp(pm[0]) - 1.0) * con_error_mean;
  syserr_line = (exp(pm[1]) - 1.0) * line_error_mean;
  syserr_radio = (exp(pm[2]) - 1.0) * radio_error_mean;
  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);

  /* con - con */
  for(i=0; i<n_con_data; i++)
  {
    t1 = Tcon_data[i];
    for(j=0; j<i; j++)
    {
      t2 = Tcon_data[j];
      PCmat_data[i*n_all_data + j] = PCmat_data[j*n_all_data + i] = 
            sig_d*sig_d * exp(-fabs(t2-t1)/tau_d) + sig_j*sig_j * exp(-fabs(t2-t1)/tau_j);
    }
    error = Fcerrs_data[i]*Fcerrs_data[i] + syserr_con*syserr_con;
    PCmat_data[i*n_all_data + i] = sig_d*sig_d + sig_j*sig_j + error;
  }

  /* con - line */
  np = n_con_data;
  for(i=0; i<n_con_data; i++)
  {
    t1 = Tcon_data[i];
    for(j=0; j<n_line_data; j++)
    {
      t2 = Tline_data[j];
      PCmat_data[i*n_all_data + j+np] = PCmat_data[(j+np)*n_all_data + i] = Slc(t1, t2, model);
    }
  }

  /* con - radio */
  np = n_con_data + n_line_data;
  for(i=0; i<n_con_data; i++)
  {
    t1 = Tcon_data[i];
    for(j=0; j<n_radio_data; j++)
    {
      t2 = Tradio_data[j];
      PCmat_data[i*n_all_data + j+np] = PCmat_data[(j+np)*n_all_data + i] = Src(t1, t2, model);
    }
  }

  /* line - line */
  np = n_con_data;
  for(i=0; i<n_line_data; i++)
  {
    t1 = Tline_data[i];
    for(j=0; j<i; j++)
    {
      t2 = Tline_data[j];
      PCmat_data[(np+i)*n_all_data + (np+j)] = PCmat_data[(np+j)*n_all_data + (np+i)] = Sll(t1, t2, model);
    }
    error = Flerrs_data[i]*Flerrs_data[i] + syserr_line*syserr_line;
    PCmat_data[(np+i)*n_all_data + (np+i)] = Sll(t1, t1, model) + error;
  }

  /* line - radio */
  np = n_con_data + n_line_data; 
  for(i=0; i<n_line_data; i++)
  {
    for(j=0; j<n_radio_data; j++)
    {
      PCmat_data[(n_con_data+i)*n_all_data + (np+j)] = PCmat_data[(np+j)*n_all_data + (n_con_data+i)] = 0.0;
    }
  }

  /* radio - radio */
  np = n_con_data + n_line_data;
  for(i=0; i<n_radio_data; i++)
  {
    t1 = Tradio_data[i];
    for(j=0; j<i; j++)
    {
      t2 = Tradio_data[j];
      PCmat_data[(np+i)*n_all_data + (np+j)] = PCmat_data[(np+j)*n_all_data + (np+i)] = Srr(t1, t2, model);
    }
    error = Frerrs_data[i] * Frerrs_data[i] + syserr_radio*syserr_radio;
    PCmat_data[(np+i)*n_all_data + (np+i)] = Srr(t1, t1, model) + error;
  }
  
  return;
}

double Slc(double tcon, double tline, const void *model)
{
  double *pm = (double *)model;
  double St;
  double sig_d, tau_d, fg, taug, wg, Dt, DT;
  
  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  
  Dt = tline - tcon;
  fg = exp(pm[7]);
  taug = pm[8];
  wg = exp(pm[9]);
  
  DT = Dt - taug;
  St = exp(-DT/tau_d + gsl_sf_log_erfc( -(DT/wg - wg/tau_d)/sqrt(2.0) ) +  wg*wg/2.0/tau_d/tau_d )
      +exp( DT/tau_d + gsl_sf_log_erfc(  (DT/wg + wg/tau_d)/sqrt(2.0) ) +  wg*wg/2.0/tau_d/tau_d );
  
  St *= 1.0/2.0 * fg * sig_d*sig_d;
  return St;
}

double Sll(double t1, double t2, const void *model)
{
  double *pm = (double *)model;
  double St;
  double sig_d, tau_d, fg, taug, wg, Dt, DT;

  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  
  Dt = t2 - t1;
  fg = exp(pm[7]);
  taug = pm[8];
  wg = exp(pm[9]);
  
  DT = Dt;

  St = exp( -DT/tau_d + gsl_sf_log_erfc( -(DT/wg/2.0 - wg/tau_d) ) + wg*wg/tau_d/tau_d )
      +exp(  DT/tau_d + gsl_sf_log_erfc(  (DT/wg/2.0 + wg/tau_d) ) + wg*wg/tau_d/tau_d ) ;

  St *= 1.0/2.0 * fg*fg * sig_d*sig_d;
  return St;
}

double Src(double tcon, double tradio, const void *model)
{
  double *pm = (double *)model;
  double St;
  double sig_j, tau_j, fg, taug, wg, Dt, DT;
  
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  
  Dt = tradio - tcon;
  fg = exp(pm[10]);
  taug = pm[11];
  wg = exp(pm[12]);
  
  DT = Dt - taug;
  St = exp(-DT/tau_j + gsl_sf_log_erfc( -(DT/wg - wg/tau_j)/sqrt(2.0) ) +  wg*wg/2.0/tau_j/tau_j )
      +exp( DT/tau_j + gsl_sf_log_erfc(  (DT/wg + wg/tau_j)/sqrt(2.0) ) +  wg*wg/2.0/tau_j/tau_j );
  
  St *= 1.0/2.0 * fg * sig_j*sig_j;

  return St;
}

double Srr(double t1, double t2, const void *model)
{
  double *pm = (double *)model;
  double St;
  double sig_j, tau_j, fg, taug, wg, Dt, DT;

  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  
  Dt = t2 - t1;
  fg = exp(pm[10]);
  taug = pm[11];
  wg = exp(pm[12]);
  
  DT = Dt;

  St = exp( -DT/tau_j + gsl_sf_log_erfc( -(DT/wg/2.0 - wg/tau_j) ) + wg*wg/tau_j/tau_j )
      +exp(  DT/tau_j + gsl_sf_log_erfc(  (DT/wg/2.0 + wg/tau_j) ) + wg*wg/tau_j/tau_j ) ;

  St *= 1.0/2.0 * fg*fg * sig_j*sig_j;
  return St;
}