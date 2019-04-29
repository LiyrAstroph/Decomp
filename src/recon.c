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
#include <math.h>
#include <float.h>
#include <string.h>
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

  //strcpy(argv[argc++], "-l");

  recon_init();

  mc_dnest(argc, argv);

  recon_postprocess();

  recon_end();
  
  //clear up argv 
  for(i=0; i<9; i++)
  {
    free(argv[i]);
  }
  free(argv);
}

void recon_postprocess()
{
  char posterior_sample_file[BRAINS_MAX_STR_LENGTH];
  double *pm, *pmstd;
  int num_ps, i, j, np;
  void *posterior_sample, *post_model;
  int size_of_modeltype = num_params * sizeof(double);

  if(thistask == roottask)
  {
    char fname[200];
    FILE *fp, *fcon, *fline, *fradio, *fconj, *fcond;

    /* get file name of posterior sample file */
    get_posterior_sample_file(dnest_options_file, posterior_sample_file);

    /* open file for posterior sample */
    fp = fopen(posterior_sample_file, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }
    
    /* open file for continuum reconstruction */
    sprintf(fname, "%s/%s", parset.file_dir, "data/con_rec.txt");
    fcon = fopen(fname, "w");
    if(fcon == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }
    sprintf(fname, "%s/%s", parset.file_dir, "data/line_rec.txt");
    fline = fopen(fname, "w");
    if(fline == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }
    sprintf(fname, "%s/%s", parset.file_dir, "data/radio_rec.txt");
    fradio = fopen(fname, "w");
    if(fradio == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }
    sprintf(fname, "%s/%s", parset.file_dir, "data/conj_rec.txt");
    fconj = fopen(fname, "w");
    if(fconj == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }
    sprintf(fname, "%s/%s", parset.file_dir, "data/cond_rec.txt");
    fcond = fopen(fname, "w");
    if(fcond == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }

    /* read number of points in posterior sample */
    if(fscanf(fp, "# %d", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(0);
    }
    printf("# Number of points in posterior sample: %d\n", num_ps);

    post_model = malloc(size_of_modeltype);
    posterior_sample = malloc(num_ps * size_of_modeltype);

    for(i=0; i<num_ps; i++)
    {
      for(j=0; j<num_params; j++)
      {
        if(fscanf(fp, "%lf", (double *)post_model + j) < 1)
        {
          fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
          exit(0);
        }
      }
      fscanf(fp, "\n");

      memcpy(posterior_sample+i*size_of_modeltype, post_model, size_of_modeltype);
      
      reconstruct_all(post_model);
      reconstruct_conjd(post_model);

      for(j=0; j<n_con_rec; j++)
      {
        fprintf(fcon, "%f %f %f\n", Tall_rec[j], Fall_rec[j], Feall_rec[j]/con_scale);
      }
      fprintf(fcon, "\n");
      
      np = n_con_rec;
      for(j=0; j<n_line_rec; j++)
      {
        fprintf(fline, "%f %f %f\n", Tall_rec[np+j], Fall_rec[np+j]/line_scale, Feall_rec[np+j]/line_scale);
      }
      fprintf(fline, "\n");
      
      np = n_con_rec + n_line_rec;
      for(j=0; j<n_radio_rec; j++)
      {
        fprintf(fradio, "%f %f %f\n", Tall_rec[np+j], Fall_rec[np+j]/radio_scale, Feall_rec[np+j]/radio_scale);
      }
      fprintf(fradio, "\n");

      np = 0;
      for(j=0; j<n_con_rec; j++)
      {
        fprintf(fcond, "%f %f %f\n", Tconjd_rec[np+j], Fconjd_rec[np+j], Feconjd_rec[np+j]/con_scale);
      }
      fprintf(fcond, "\n");

      np = n_con_rec;
      for(j=0; j<n_radio_rec; j++)
      {
        fprintf(fconj, "%f %f %f\n", Tconjd_rec[np+j], Fconjd_rec[np+j], Feconjd_rec[np+j]/con_scale);
      }
      fprintf(fconj, "\n");
    }
    
    free(post_model);
    free(posterior_sample);

    fclose(fp);
    fclose(fcon);
    fclose(fline);
    fclose(fradio);
    fclose(fconj);
    fclose(fcond);
  }

  return;
}

void recon_init()
{
  int i, np;
  double fac, Tspan;
  /* option file */
  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONS");

  fac = 2.0;
  n_con_rec = (int)(n_con_data * fac);
  n_line_rec = (int)(n_line_data * fac);
  n_radio_rec = (int)(n_radio_data * fac);
  n_all_rec = n_con_rec + n_line_rec + n_radio_rec;

  n_conjd_rec = n_con_rec + n_radio_rec;

  Larr_rec = malloc(n_all_rec * 3 * sizeof(double));
  Larr_conjd = malloc(n_conjd_rec * 3 * sizeof(double));
  USmat_rec = malloc(n_all_rec * n_all_data * sizeof(double));
  USmatT_rec = malloc(n_all_rec * n_all_data * sizeof(double));
  ASmat_rec = malloc(n_all_rec * n_all_rec * sizeof(double));
  Tall_rec = malloc(n_all_rec * sizeof(double));
  Fall_rec = malloc(n_all_rec * sizeof(double));
  Feall_rec = malloc(n_all_rec * sizeof(double));
  Tconjd_rec = malloc(n_conjd_rec * sizeof(double));
  Fconjd_rec = malloc(n_conjd_rec * sizeof(double));
  Feconjd_rec = malloc(n_conjd_rec * sizeof(double));
  PEmat1 = malloc(n_all_rec * n_all_rec * sizeof(double));
  PEmat2 = malloc(n_all_rec * n_all_rec * sizeof(double));

  /* con grid */
  np = 0;
  Tspan = Tcon_data[n_con_data -1] - Tcon_data[0];
  for(i=0; i<n_con_rec; i++)
  {
    Tall_rec[np+i] = (Tcon_data[0] - 0.2*Tspan) + (Tspan + 0.4*Tspan)/(n_con_rec - 1.0) * i;
  }

  /* line grid */
  np = n_con_rec;
  Tspan = Tline_data[n_line_data -1] - Tline_data[0];
  for(i=0; i<n_line_rec; i++)
  {
    Tall_rec[np+i] = (Tline_data[0] - 0.2*Tspan) + (Tspan + 0.4*Tspan)/(n_line_rec - 1.0) * i;
  }

  /* radio grid */
  np += n_line_rec;
  Tspan = Tradio_data[n_radio_data -1] - Tradio_data[0];
  for(i=0; i<n_radio_rec; i++)
  {
    Tall_rec[np+i] = (Tradio_data[0] - 0.2*Tspan) + (Tspan + 0.4*Tspan)/(n_radio_rec - 1.0) * i;
  }

  /* setup Larr_data */
  for(i=0; i<n_con_rec; i++)
  {
    Larr_rec[i*3 + 0] = 1.0;
    Larr_rec[i*3 + 1] = 0.0;
    Larr_rec[i*3 + 2] = 0.0;
  }
  np = n_con_rec;
  for(i=0; i<n_line_rec; i++)
  {
    Larr_rec[(np+i)*3 + 0] = 0.0;
    Larr_rec[(np+i)*3 + 1] = 1.0;
    Larr_rec[(np+i)*3 + 2] = 0.0;
  }
  np = n_con_rec + n_line_rec;
  for(i=0; i<n_radio_rec; i++)
  {
    Larr_rec[(np+i)*3 + 0] = 0.0;
    Larr_rec[(np+i)*3 + 1] = 0.0;
    Larr_rec[(np+i)*3 + 2] = 1.0;
  }

  /* disk and jet emissions at con  */
  memcpy(Tconjd_rec, Tall_rec, n_con_rec*sizeof(double));
  memcpy(Tconjd_rec + n_con_rec, Tall_rec + n_con_rec+n_line_rec, n_radio_rec*sizeof(double));

  /* setup Larr_data */
  for(i=0; i<n_con_rec; i++)
  {
    Larr_conjd[i*3 + 0] = 0.5;
    Larr_conjd[i*3 + 1] = 0.0;
    Larr_conjd[i*3 + 2] = 0.0;
  }
  np = n_con_rec;
  for(i=0; i<n_radio_rec; i++)
  {
    Larr_conjd[(np+i)*3 + 0] = 0.5;
    Larr_conjd[(np+i)*3 + 1] = 0.0;
    Larr_conjd[(np+i)*3 + 2] = 0.0;
  }

  return;
}

void recon_end()
{
  int i;
  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
  }
  free(par_range_model);

  free(Larr_rec);
  free(Larr_conjd);
  free(USmat_rec);
  free(USmatT_rec);
  free(ASmat_rec);
  free(Tall_rec);
  free(Fall_rec);
  free(Feall_rec);
  free(Tconjd_rec);
  free(Fconjd_rec);
  free(Feconjd_rec);
  free(PEmat1);
  free(PEmat2);
}

void reconstruct_all(const void *model)
{
  int i, j, nq=3;
  double *ICq, *yq, *ybuf, *y, *yave, *yave_rec;

  ICq = workspace;
  yq = ICq + nq*nq;
  yave = yq + nq;
  y = yave + n_all_data;
  ybuf = y + n_all_data;
  yave_rec = ybuf + n_all_data;

  set_covar_Pmat_data(model);
  set_covar_Umat_rec(model);

  memcpy(Tmat1, PCmat_data, n_all_data*n_all_data*sizeof(double));
  memcpy(Tmat2, Larr_data, n_all_data*nq*sizeof(double));

  multiply_mat_MN_inverseA(Tmat1, Tmat2, n_all_data, nq); // Tmat2 = C^-1 * L;  NxNq
  
  multiply_mat_MN_transposeA(Larr_data, Tmat2, ICq, nq, nq, n_all_data); // ICq = L^T*C^-1*L; NqxNq
  multiply_mat_MN_transposeA(Tmat2, Fall_data, yq, nq, 1, n_all_data); // yq = L^T*C^-1*y;  Nqx1
  memcpy(Tmat1, ICq, nq*nq*sizeof(double));
  multiply_mat_MN_inverseA(Tmat1, yq, nq, 1); // yq = (L^T*C^-1*L)^-1 * L^T*C^-1*y; Nqx1

  multiply_mat_MN(Larr_data, yq, yave, n_all_data, 1, nq); // yave = L * q; Nx1

  for(i=0; i<n_all_data; i++)y[i] = Fall_data[i] - yave[i];
  memcpy(Tmat1, PCmat_data, n_all_data*n_all_data*sizeof(double));
  memcpy(ybuf, y, n_all_data*sizeof(double));
  multiply_mat_MN_inverseA(Tmat1, ybuf, n_all_data, 1); // ybuf = C^-1 * y; Nx1

  multiply_matvec_MN(USmat_rec, n_all_rec, n_all_data, ybuf, Fall_rec); // Fall_rec = S*C^-1*y
  multiply_matvec_MN(Larr_rec, n_all_rec, nq, yq, yave_rec);   // yave_rec = L*yq
  
  for(i=0; i<n_all_rec; i++)
  {
    Fall_rec[i] += yave_rec[i];
  }

  // get errors
  /* S x C^-1 x S */
  /* Transpose of USmat */
  for(i=0; i<n_all_rec; i++)
    for(j=0; j<n_all_data; j++)
      USmatT_rec[j*n_all_rec + i] = USmat_rec[i*n_all_data + j];
  
  memcpy(PEmat1, PCmat_data, n_all_data*n_all_data*sizeof(double));
  memcpy(PEmat2, USmatT_rec, n_all_rec*n_all_data*sizeof(double));

  multiply_mat_MN_inverseA(PEmat1, PEmat2, n_all_data, n_all_rec); // C^-1 x S; NdxN
  multiply_mat_MN(USmat_rec, PEmat2, PEmat1, n_all_rec, n_all_rec, n_all_data); // S x C^-1 x S; NxN

  set_covar_Amat_rec(model);

  for(i=0; i<n_all_rec; i++)
  {
    Feall_rec[i] = sqrt(ASmat_rec[i*n_all_rec + i] - PEmat1[i*n_all_rec + i]);
  }
}

void reconstruct_all2(const void *model)
{
  int i, nq=3, info;
  double *Cq, *yq, *ybuf, *y, *yave, *yave_rec;

  Cq = workspace;
  yq = Cq + nq*nq;
  yave = yq + nq;
  y = yave + n_all_data;
  ybuf = y + n_all_data;
  yave_rec = ybuf + n_all_data * nq;

  set_covar_Pmat_data(model);
  set_covar_Umat_rec(model);

  /* C^-1 */
  memcpy(IPCmat_data, PCmat_data, n_all_data*n_all_data*sizeof(double));
  inverse_mat(IPCmat_data, n_all_data, &info); 

  /* L^T*C^-1*L */
  multiply_mat_MN(IPCmat_data, Larr_data, ybuf, n_all_data, nq, n_all_data); // ybuf = C^-1*L; nd*nq
  multiply_mat_MN_transposeA(Larr_data, ybuf, Cq, nq, nq, n_all_data); // Cq = L^T*C^-1*L; nq*nq
  
  /* L^T*C^-1*y */
  multiply_matvec(IPCmat_data, Fall_data, n_all_data, ybuf);   // ybuf = C^-1*Fall_data; nd*1
  multiply_mat_MN_transposeA(Larr_data, ybuf, yq, nq, 1, n_all_data); // yq = L^T*C^-1*Fall_data; nq*1

  /* (L^T*C^-1*L)^-1 * L^T*C^-1*y */
  inverse_mat(Cq, nq, &info);
  multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq); // ybuf = (L^T*C^-1*L)^-1 * L^T*C^-1*y; nq*1
  
  multiply_matvec_MN(Larr_data, n_all_data, nq, ybuf, yave); // yave = L*ybuf; nd*1
  for(i=0; i<n_all_data; i++)
  {
    y[i] = Fall_data[i] - yave[i];
  }
  multiply_matvec(IPCmat_data, y, n_all_data, yave);  // yave = C^-1*y; nd*1
  
  /* S*C^-1*y */
  multiply_matvec_MN(USmat_rec, n_all_rec, n_all_data, yave, Fall_rec); // Fall_rec =  S*C^-1*y; nr*1
  multiply_matvec_MN(Larr_rec, n_all_rec, nq, ybuf, yave_rec);   

  for(i=0; i<n_all_rec; i++)
  {
    Fall_rec[i] += yave_rec[i];
  }

  // get errors
  /* S x C^-1 x S */
  set_covar_Amat_rec(model);
  multiply_mat_MN(USmat_rec, IPCmat_data, PEmat1, n_all_rec, n_all_data, n_all_data);
  multiply_mat_MN_transposeB(PEmat1, USmat_rec, PEmat2, n_all_rec, n_all_rec, n_all_data);
  
  for(i=0; i<n_all_rec; i++)
  {
    Feall_rec[i] = sqrt(ASmat_rec[i*n_all_rec + i] - PEmat2[i*n_all_rec + i]);
  }

  return;
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
  double *pm = (double *)model, sig_j, sig_d, sig2_all;
  
  sig_d = exp(pm[3] + 0.5*pm[4]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  sig2_all = sig_d*sig_d + sig_j*sig_j;

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
  multiply_mat_MN_transposeA(Tmat2, Fall_data, yq, nq, 1, n_all_data); // yq = L^T*C^-1*y;  Nqx1
  memcpy(Tmat1, ICq, nq*nq*sizeof(double));
  multiply_mat_MN_inverseA(Tmat1, yq, nq, 1); // yq = (L^T*C^-1*L)^-1 * L^T*C^-1*y; Nqx1

  multiply_mat_MN(Larr_data, yq, yave, n_all_data, 1, nq); // yave = L * q; Nx1

  for(i=0; i<n_all_data; i++)y[i] = Fall_data[i] - yave[i];
  memcpy(Tmat1, PCmat_data, n_all_data*n_all_data*sizeof(double));
  memcpy(ybuf, y, n_all_data*sizeof(double));
  multiply_mat_MN_inverseA(Tmat1, ybuf, n_all_data, 1); // ybuf = C^-1 * y; Nx1

  prob = -0.5*cblas_ddot(n_all_data, y, 1, ybuf, 1)/sig2_all; // y^T * C^-1 * y
  if(prob > 0.0 )  // check if prob is positive
  {
    prob = -DBL_MAX;
    printf("prob >0!\n");
    return prob;
  }

  lndet = lndet_mat3(PCmat_data, n_all_data, &info, &sign) + n_all_data * log(sig2_all);
  if(info!=0|| sign==-1)
  {
    prob = -DBL_MAX;
    printf("lndet_C %f %d!\n", lndet, sign);
    return prob;
  }

  lndet_ICq = lndet_mat3(ICq, nq, &info, &sign) + nq*log(sig2_all);
  if(info!=0 || sign==-1 )
  {
    prob = -DBL_MAX;
    printf("lndet_ICq %f %d!\n", lndet_ICq, sign);
    return prob;
  }

  prob += -0.5*lndet - 0.5*lndet_ICq;
  
  return prob;
}

/*!
 * matrix operation A^-1 x B is implemented by calling functions
 *    inverse_mat()
 *    multiply_mat_MN()
 */
double prob2(const void *model)
{
  int i, nq=3, info, sign;
  double prob=0.0, lndet, lndet_ICq;
  double *Cq, *ICq, *yq, *ybuf, *y, *yave;
  double *pm = (double *)model, sig_j, sig_d, sig2_all;

  sig_d = exp(pm[3] + 0.5*pm[4]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  sig2_all = sig_d*sig_d + sig_j*sig_j;
  
  Cq = workspace;
  ICq = Cq + nq*nq;
  yq = ICq + nq*nq;
  yave = yq + nq;
  y = yave + n_all_data;
  ybuf = y + n_all_data;

  set_covar_Pmat_data(model);

  /* C^-1 */
  memcpy(IPCmat_data, PCmat_data, n_all_data*n_all_data*sizeof(double));
  inverse_mat(IPCmat_data, n_all_data, &info); 

  /* L^T*C^-1*L */
  multiply_mat_MN(IPCmat_data, Larr_data, ybuf, n_all_data, nq, n_all_data);
  multiply_mat_MN_transposeA(Larr_data, ybuf, Cq, nq, nq, n_all_data);
  memcpy(ICq, Cq, nq*nq*sizeof(double));

  /* L^T*C^-1*y */
  multiply_matvec(IPCmat_data, Fall_data, n_all_data, ybuf);
  multiply_mat_MN_transposeA(Larr_data, ybuf, yq, nq, 1, n_all_data);

  /* (L^T*C^-1*L)^-1 * L^T*C^-1*y */
  inverse_mat(Cq, nq, &info);
  multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);

  multiply_matvec_MN(Larr_data, n_all_data, nq, ybuf, yave);
  for(i=0; i<n_all_data; i++)
  {
    y[i] = Fall_data[i] - yave[i];
  }

  /* y^T x C^-1 x y */
  multiply_matvec(IPCmat_data, y, n_all_data, ybuf);
  prob = -0.5 * cblas_ddot(n_all_data, y, 1, ybuf, 1) / sig2_all;
  
  if(prob > 0.0 )  // check if prob is positive
  { 
    prob = -DBL_MAX;
    printf("prob >0!\n");
    return prob;
  }
  lndet = lndet_mat3(PCmat_data, n_all_data, &info, &sign) + n_all_data*log(sig2_all);;
  if(info!=0|| sign==-1)
  {
    prob = -DBL_MAX;
    printf("lndet_C %f %d!\n", lndet, sign);
    return prob;
  }
  lndet_ICq = lndet_mat3(ICq, nq, &info, &sign) - nq*log(sig2_all);
  if(info!=0 || sign==-1 )
  {
    prob = -DBL_MAX;
    printf("lndet_ICq %f %d!\n", lndet_ICq, sign);
    return prob;
  }
  prob += prob - 0.5*lndet - 0.5*lndet_ICq;
  
  return prob;
}

void set_covar_Pmat_data(const void *model)
{
  int i, j, np;
  double *pm = (double *)model;
  double syserr_con, syserr_line, syserr_radio, sig_d, tau_d, sig_j, tau_j;
  double t1, t2, error, sig2_all;

  syserr_con = (exp(pm[0]) - 1.0) * con_error_mean;
  syserr_line = (exp(pm[1]) - 1.0) * line_error_mean;
  syserr_radio = (exp(pm[2]) - 1.0) * radio_error_mean;
  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);

  sig2_all = sig_d*sig_d + sig_j*sig_j;

  for(i=0; i<n_con_data; i++)
  {
    t1 = Tcon_data[i];

    /* con - con */
    for(j=0; j<i; j++)
    {
      t2 = Tcon_data[j];
      PCmat_data[i*n_all_data + j] = PCmat_data[j*n_all_data + i] = 
            ( sig_d*sig_d * exp(-fabs(t2-t1)/tau_d) + sig_j*sig_j * exp(-fabs(t2-t1)/tau_j) )/sig2_all ;
    }
    error = Fcerrs_data[i]*Fcerrs_data[i] + syserr_con*syserr_con;
    PCmat_data[i*n_all_data + i] = 1.0 + error/sig2_all;

    /* con - line */
    np = n_con_data;
    for(j=0; j<n_line_data; j++)
    {
      t2 = Tline_data[j];
      PCmat_data[i*n_all_data + j+np] = PCmat_data[(j+np)*n_all_data + i] = Slc(t1, t2, model);
    }

    /* con - radio */
    np += n_line_data;
    for(j=0; j<n_radio_data; j++)
    {
      t2 = Tradio_data[j];
      PCmat_data[i*n_all_data + j+np] = PCmat_data[(j+np)*n_all_data + i] = Src(t1, t2, model);
    }
  }

  for(i=0; i<n_line_data; i++)
  {
    t1 = Tline_data[i];

    /* line - line */
    np = n_con_data;
    for(j=0; j<i; j++)
    {
      t2 = Tline_data[j];
      PCmat_data[(np+i)*n_all_data + (np+j)] = PCmat_data[(np+j)*n_all_data + (np+i)] = Sll(t1, t2, model);
    }
    error = Flerrs_data[i]*Flerrs_data[i] + syserr_line*syserr_line;
    PCmat_data[(np+i)*n_all_data + (np+i)] = Sll(t1, t1, model) + error/sig2_all;

    /* line - radio */
    np += n_line_data; 
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
    PCmat_data[(np+i)*n_all_data + (np+i)] = Srr(t1, t1, model) + error/sig2_all;
  }
  
  return;
}

/* covariance between con and line data */
double Slc(double tcon, double tline, const void *model)
{
  double *pm = (double *)model;
  double St;
  double sig_d, tau_d, sig_j, tau_j, fg, taug, wg, Dt, DT;
  
  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  
  Dt = tline - tcon;
  fg = exp(pm[7]);
  taug = pm[8];
  wg = exp(pm[9]);
  
  DT = Dt - taug;
  St = exp(-DT/tau_d + gsl_sf_log_erfc( -(DT/wg - wg/tau_d)/sqrt(2.0) ) +  wg*wg/2.0/tau_d/tau_d )
      +exp( DT/tau_d + gsl_sf_log_erfc(  (DT/wg + wg/tau_d)/sqrt(2.0) ) +  wg*wg/2.0/tau_d/tau_d );
  
  St *= 1.0/2.0 * fg * sig_d*sig_d/(sig_d*sig_d+sig_j*sig_j);
  return St;
}

/* covariance of line data */
double Sll(double t1, double t2, const void *model)
{
  double *pm = (double *)model;
  double St;
  double sig_d, tau_d, sig_j, tau_j, fg, taug, wg, Dt, DT;

  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  
  Dt = t2 - t1;
  fg = exp(pm[7]);
  taug = pm[8];
  wg = exp(pm[9]);
  
  DT = Dt;

  St = exp( -DT/tau_d + gsl_sf_log_erfc( -(DT/wg/2.0 - wg/tau_d) ) + wg*wg/tau_d/tau_d )
      +exp(  DT/tau_d + gsl_sf_log_erfc(  (DT/wg/2.0 + wg/tau_d) ) + wg*wg/tau_d/tau_d ) ;

  St *= 1.0/2.0 * fg*fg * sig_d*sig_d/(sig_d*sig_d+sig_j*sig_j);
  return St;
}

/* covariance between con and radio data */
double Src(double tcon, double tradio, const void *model)
{
  double *pm = (double *)model;
  double St;
  double sig_d, tau_d, sig_j, tau_j, fg, taug, wg, Dt, DT;
  
  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  
  Dt = tradio - tcon;
  fg = exp(pm[10]);
  taug = pm[11];
  wg = exp(pm[12]);
  
  DT = Dt - taug;
  St = exp(-DT/tau_j + gsl_sf_log_erfc( -(DT/wg - wg/tau_j)/sqrt(2.0) ) +  wg*wg/2.0/tau_j/tau_j )
      +exp( DT/tau_j + gsl_sf_log_erfc(  (DT/wg + wg/tau_j)/sqrt(2.0) ) +  wg*wg/2.0/tau_j/tau_j );
  
  St *= 1.0/2.0 * fg * sig_j*sig_j/(sig_d*sig_d+sig_j*sig_j);

  return St;
}

/* covariance of radio data */
double Srr(double t1, double t2, const void *model)
{
  double *pm = (double *)model;
  double St;
  double sig_d, tau_d, sig_j, tau_j, fg, taug, wg, Dt, DT;

  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  
  Dt = t2 - t1;
  fg = exp(pm[10]);
  taug = pm[11];
  wg = exp(pm[12]);
  
  DT = Dt;

  St = exp( -DT/tau_j + gsl_sf_log_erfc( -(DT/wg/2.0 - wg/tau_j) ) + wg*wg/tau_j/tau_j )
      +exp(  DT/tau_j + gsl_sf_log_erfc(  (DT/wg/2.0 + wg/tau_j) ) + wg*wg/tau_j/tau_j ) ;

  St *= 1.0/2.0 * fg*fg * sig_j*sig_j/(sig_d*sig_d+sig_j*sig_j);
  return St;
}


/*
 * covariances between reconstructed and data time points.
 * 
 */
void set_covar_Umat_rec(const void *model)
{
  int i, j, np_data, np_rec;
  double *pm = (double *)model;
  double sig_d, tau_d, sig_j, tau_j, sig2_all;
  double t1, t2;

  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  sig2_all = sig_d*sig_d + sig_j*sig_j;

  np_rec = 0;
  for(i=0; i<n_con_rec; i++)
  {
    t1 = Tall_rec[np_rec + i];

    /* con - con */
    for(j=0; j<n_con_data; j++)
    {
      t2 = Tcon_data[j];
      USmat_rec[(np_rec+i)*n_all_data + j] = 
         (sig_d*sig_d * exp(-fabs(t2-t1)/tau_d) + sig_j*sig_j * exp(-fabs(t2-t1)/tau_j))/sig2_all;
    }
    /* con - line */
    np_data = n_con_data;
    for(j=0; j<n_line_data; j++)
    {
      t2 = Tline_data[j];
      USmat_rec[(np_rec+i)*n_all_data + j+np_data] = Slc(t1, t2, model);
    }
    
    /* con - radio */
    np_data += n_line_data;
    for(j=0; j<n_radio_data; j++)
    {
      t2 = Tradio_data[j];
      USmat_rec[(np_rec +i)*n_all_data + j+np_data] = Src(t1, t2, model);
    }
  }


  /* line  */
  np_rec = n_con_rec;
  for(i=0; i<n_line_rec; i++)
  {
    t1 = Tall_rec[np_rec + i];

    /* line - con */
    np_data = 0;
    for(j=0; j<n_con_data; j++)
    {
      t2 = Tcon_data[j];
      USmat_rec[(np_rec+i)*n_all_data + (np_data+j)] = Slc(t2, t1, model);
    }

    /* line - line */
    np_data += n_con_data; 
    for(j=0; j<n_line_data; j++)
    {
      t2 = Tline_data[j];
      USmat_rec[(np_rec+i)*n_all_data + (np_data+j)] = Sll(t1, t2, model);
    }
    
    /* line - radio */
    np_data += n_line_data; 
    for(j=0; j<n_radio_data; j++)
    {
      USmat_rec[(np_rec+i)*n_all_data + (np_data+j)] = 0.0;
    }
  }


  /* radio */
  np_rec += n_line_rec;
  for(i=0; i<n_radio_rec; i++)
  {
    t1 = Tall_rec[np_rec + i];

    /* radio - con */
    np_data = 0;
    for(j=0; j<n_con_data; j++)
    {
      t2 = Tcon_data[j];
      USmat_rec[(np_rec+i)*n_all_data + (np_data+j)] = Src(t2, t1, model);
    }
    
    /* radio - line */
    np_data += n_con_data;
    for(j=0; j<n_line_data; j++)
    {
      USmat_rec[(np_rec+i)*n_all_data + (np_data+j)] = 0.0;
    }
    
    /* radio - radio */
    np_data += n_line_data;
    for(j=0; j<n_radio_data; j++)
    {
      t2 = Tradio_data[j];
      USmat_rec[(np_rec+i)*n_all_data + (np_data+j)] = Srr(t1, t2, model);
    }
  }

  return;
}

/*
 * covariances between reconstructed time points.
 *  
 */
void set_covar_Amat_rec(const void *model)
{
  int i, j, np;
  double *pm = (double *)model;
  double syserr_con, syserr_line, syserr_radio, sig_d, tau_d, sig_j, tau_j, sig2_all;
  double t1, t2, error;
  
  syserr_con = (exp(pm[0]) - 1.0) * con_error_mean;
  syserr_line = (exp(pm[1]) - 1.0) * line_error_mean;
  syserr_radio = (exp(pm[2]) - 1.0) * radio_error_mean;
  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  sig2_all = sig_d*sig_d + sig_j*sig_j;

  for(i=0; i<n_con_rec; i++)
  {
    t1 = Tall_rec[i];

    /* con - con */
    for(j=0; j<i; j++)
    {
      t2 = Tall_rec[j];
      ASmat_rec[i*n_all_rec + j] = ASmat_rec[j*n_all_rec + i] = 
            ( sig_d*sig_d * exp(-fabs(t2-t1)/tau_d) + sig_j*sig_j * exp(-fabs(t2-t1)/tau_j) ) / sig2_all;
    }
    error = syserr_con*syserr_con;
    ASmat_rec[i*n_all_rec + i] = 1.0 + error/sig2_all;

    /* con - line */
    np = n_con_rec;
    for(j=0; j<n_line_rec; j++)
    {
      t2 = Tall_rec[np+j];
      ASmat_rec[i*n_all_rec + j+np] = ASmat_rec[(j+np)*n_all_rec + i] = Slc(t1, t2, model);
    }

    /* con - radio */
    np += n_line_rec;
    for(j=0; j<n_radio_rec; j++)
    {
      t2 = Tall_rec[np+j];
      ASmat_rec[i*n_all_rec + j+np] = ASmat_rec[(j+np)*n_all_rec + i] = Src(t1, t2, model);
    }
  }

  for(i=0; i<n_line_rec; i++)
  {
    t1 = Tall_rec[n_con_rec+i];

    /* line - line */
    np = n_con_rec;
    for(j=0; j<i; j++)
    {
      t2 = Tall_rec[np+j];
      ASmat_rec[(np+i)*n_all_rec + (np+j)] = ASmat_rec[(np+j)*n_all_rec + (np+i)] = Sll(t1, t2, model);
    }
    error = syserr_line*syserr_line;
    ASmat_rec[(np+i)*n_all_rec + (np+i)] = Sll(t1, t1, model) + error/sig2_all;

    /* line - radio */
    np += n_line_rec; 
    for(j=0; j<n_radio_rec; j++)
    {
      ASmat_rec[(n_con_rec+i)*n_all_rec + (np+j)] = ASmat_rec[(np+j)*n_all_rec + (n_con_rec+i)] = 0.0;
    }
  }

  /* radio - radio */
  np = n_con_rec + n_line_rec;
  for(i=0; i<n_radio_rec; i++)
  {
    t1 = Tall_rec[np+i];
    for(j=0; j<i; j++)
    {
      t2 = Tall_rec[np+j];
      ASmat_rec[(np+i)*n_all_rec + (np+j)] = ASmat_rec[(np+j)*n_all_rec + (np+i)] = Srr(t1, t2, model);
    }
    error = syserr_radio*syserr_radio;
    ASmat_rec[(np+i)*n_all_rec + (np+i)] = Srr(t1, t1, model) + error/sig2_all;
  }
  
  return;
}

/*
 * covariances of data time points.
 * not include measurement noises.
 * 
 */ 
void set_covar_PSmat_data(const void *model)
{
  int i, j, np;
  double *pm = (double *)model;
  double syserr_con, syserr_line, syserr_radio, sig_d, tau_d, sig_j, tau_j;
  double t1, t2, error, sig2_all;

  syserr_con = (exp(pm[0]) - 1.0) * con_error_mean;
  syserr_line = (exp(pm[1]) - 1.0) * line_error_mean;
  syserr_radio = (exp(pm[2]) - 1.0) * radio_error_mean;
  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);

  sig2_all = sig_d*sig_d + sig_j*sig_j;

  for(i=0; i<n_con_data; i++)
  {
    t1 = Tcon_data[i];

    /* con - con */
    for(j=0; j<i; j++)
    {
      t2 = Tcon_data[j];
      PCmat_data[i*n_all_data + j] = PCmat_data[j*n_all_data + i] = 
            ( sig_d*sig_d * exp(-fabs(t2-t1)/tau_d) + sig_j*sig_j * exp(-fabs(t2-t1)/tau_j) )/sig2_all ;
    }
    PCmat_data[i*n_all_data + i] = 1.0;

    /* con - line */
    np = n_con_data;
    for(j=0; j<n_line_data; j++)
    {
      t2 = Tline_data[j];
      PCmat_data[i*n_all_data + j+np] = PCmat_data[(j+np)*n_all_data + i] = Slc(t1, t2, model);
    }

    /* con - radio */
    np += n_line_data;
    for(j=0; j<n_radio_data; j++)
    {
      t2 = Tradio_data[j];
      PCmat_data[i*n_all_data + j+np] = PCmat_data[(j+np)*n_all_data + i] = Src(t1, t2, model);
    }
  }

  for(i=0; i<n_line_data; i++)
  {
    t1 = Tline_data[i];

    /* line - line */
    np = n_con_data;
    for(j=0; j<i; j++)
    {
      t2 = Tline_data[j];
      PCmat_data[(np+i)*n_all_data + (np+j)] = PCmat_data[(np+j)*n_all_data + (np+i)] = Sll(t1, t2, model);
    }
    PCmat_data[(np+i)*n_all_data + (np+i)] = Sll(t1, t1, model);

    /* line - radio */
    np += n_line_data; 
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
    PCmat_data[(np+i)*n_all_data + (np+i)] = Srr(t1, t1, model) + error/sig2_all;
  }
  
  return;
}

/*
 * covariances between reconstructed and data time points, for disk and jet emissions.
 * 
 */
void set_covar_Umat_conjd(const void *model)
{
  int i, j, np_data, np_rec;
  double *pm = (double *)model;
  double sig_d, tau_d, sig_j, tau_j, sig2_all;
  double t1, t2;

  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  sig2_all = sig_d*sig_d + sig_j*sig_j;
  
  /* disk emission */
  np_rec = 0;
  for(i=0; i<n_con_rec; i++)
  {
    t1 = Tconjd_rec[np_rec + i];
    /* con - con */
    for(j=0; j<n_con_data; j++)
    {
      t2 = Tcon_data[j];
      USmat_rec[(np_rec+i)*n_all_data + j] = (sig_d*sig_d * exp(-fabs(t2-t1)/tau_d))/sig2_all;
    }
    /* con - line */
    np_data = n_con_data;
    for(j=0; j<n_line_data; j++)
    {
      t2 = Tline_data[j];
      USmat_rec[(np_rec+i)*n_all_data + j+np_data] = Slc(t1, t2, model);
    }
    
    /* con - radio */
    np_data += n_line_data;
    for(j=0; j<n_radio_data; j++)
    {
      USmat_rec[(np_rec +i)*n_all_data + j+np_data] = 0.0;
    }
  }


  /* jet emission */
  np_rec += n_con_rec;
  for(i=0; i<n_radio_rec; i++)
  {
    t1 = Tconjd_rec[np_rec + i];

    /* radio - con */
    np_data = 0;
    for(j=0; j<n_con_data; j++)
    {
      t2 = Tcon_data[j];
      USmat_rec[(np_rec+i)*n_all_data + (np_data+j)] = (sig_j*sig_j * exp(-fabs(t2-t1)/tau_j))/sig2_all;;
    }
    
    /* radio - line */
    np_data += n_con_data;
    for(j=0; j<n_line_data; j++)
    {
      USmat_rec[(np_rec+i)*n_all_data + (np_data+j)] = 0.0;
    }
    
    /* radio - radio */
    np_data += n_line_data;
    for(j=0; j<n_radio_data; j++)
    {
      t2 = Tradio_data[j];
      USmat_rec[(np_rec+i)*n_all_data + (np_data+j)] = Src(t1, t2, model);
    }
  }

  return;
}

/*
 * covariances between reconstructed time points for disk and jet.
 *  
 */
void set_covar_Amat_conjd(const void *model)
{
  int i, j, np;
  double *pm = (double *)model;
  double syserr_con, syserr_line, syserr_radio, sig_d, tau_d, sig_j, tau_j, sig2_all;
  double t1, t2, error;
  
  syserr_con = (exp(pm[0]) - 1.0) * con_error_mean;
  syserr_line = (exp(pm[1]) - 1.0) * line_error_mean;
  syserr_radio = (exp(pm[2]) - 1.0) * radio_error_mean;
  tau_d = exp(pm[4]);
  sig_d = exp(pm[3] + 0.5*pm[4]);
  tau_j = exp(pm[6]);
  sig_j = exp(pm[5] + 0.5*pm[6]);
  sig2_all = sig_d*sig_d + sig_j*sig_j;

  for(i=0; i<n_con_rec; i++)
  {
    t1 = Tconjd_rec[i];

    /* disk - disk */
    for(j=0; j<i; j++)
    {
      t2 = Tconjd_rec[j];
      ASmat_rec[i*n_conjd_rec + j] = ASmat_rec[j*n_conjd_rec + i] = 
            ( sig_d*sig_d * exp(-fabs(t2-t1)/tau_d)) / sig2_all;
    }
    error = syserr_con*syserr_con;
    ASmat_rec[i*n_conjd_rec + i] = (sig_d*sig_d + error)/sig2_all;

    /* disk - jet */
    np = n_con_rec;
    for(j=0; j<n_radio_rec; j++)
    {
      t2 = Tconjd_rec[np+j];
      ASmat_rec[i*n_conjd_rec + j+np] = ASmat_rec[(j+np)*n_conjd_rec + i] = 0.0;
    }
  }

  /* radio - radio */
  np = n_con_rec;
  for(i=0; i<n_radio_rec; i++)
  {
    t1 = Tconjd_rec[np+i];
    for(j=0; j<i; j++)
    {
      t2 = Tconjd_rec[np+j];
      ASmat_rec[(np+i)*n_conjd_rec + (np+j)] = ASmat_rec[(np+j)*n_conjd_rec + (np+i)] = 
          ( sig_j*sig_j * exp(-fabs(t2-t1)/tau_j)) / sig2_all;
    }
    error = syserr_radio*syserr_radio;
    ASmat_rec[(np+i)*n_conjd_rec + (np+i)] = (sig_j*sig_j + error)/sig2_all;
  }
  
  return;
}

/* reconstruct disk and jet emission */
void reconstruct_conjd(const void *model)
{
  int i, nq=3, info;
  double *Cq, *yq, *ybuf, *y, *yave, *yave_rec;

  Cq = workspace;
  yq = Cq + nq*nq;
  yave = yq + nq;
  y = yave + n_all_data;
  ybuf = y + n_all_data;
  yave_rec = ybuf + n_all_data * nq;

  set_covar_Pmat_data(model);
  set_covar_Umat_conjd(model);

  /* C^-1 */
  memcpy(IPCmat_data, PCmat_data, n_all_data*n_all_data*sizeof(double));
  inverse_mat(IPCmat_data, n_all_data, &info); 

  /* L^T*C^-1*L */
  multiply_mat_MN(IPCmat_data, Larr_data, ybuf, n_all_data, nq, n_all_data); // ybuf = C^-1*L; nd*nq
  multiply_mat_MN_transposeA(Larr_data, ybuf, Cq, nq, nq, n_all_data); // Cq = L^T*C^-1*L; nq*nq
  
  /* L^T*C^-1*y */
  multiply_matvec(IPCmat_data, Fall_data, n_all_data, ybuf);   // ybuf = C^-1*Fall_data; nd*1
  multiply_mat_MN_transposeA(Larr_data, ybuf, yq, nq, 1, n_all_data); // yq = L^T*C^-1*Fall_data; nq*1

  /* (L^T*C^-1*L)^-1 * L^T*C^-1*y */
  inverse_mat(Cq, nq, &info);
  multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq); // ybuf = (L^T*C^-1*L)^-1 * L^T*C^-1*y; nq*1
  
  multiply_matvec_MN(Larr_data, n_all_data, nq, ybuf, yave); // yave = L*ybuf; nd*1
  for(i=0; i<n_all_data; i++)
  {
    y[i] = Fall_data[i] - yave[i];
  }
  multiply_matvec(IPCmat_data, y, n_all_data, yave);  // yave = C^-1*y; nd*1
  
  /* S*C^-1*y */
  multiply_matvec_MN(USmat_rec, n_conjd_rec, n_all_data, yave, Fconjd_rec); // Fall_rec =  S*C^-1*y; nr*1
  multiply_matvec_MN(Larr_conjd, n_conjd_rec, nq, ybuf, yave_rec);   

  for(i=0; i<n_conjd_rec; i++)
  {
    Fconjd_rec[i] += yave_rec[i];
  }

  // get errors
  /* S x C^-1 x S */
  set_covar_Amat_conjd(model);
  multiply_mat_MN(USmat_rec, IPCmat_data, PEmat1, n_conjd_rec, n_all_data, n_all_data);
  multiply_mat_MN_transposeB(PEmat1, USmat_rec, PEmat2, n_conjd_rec, n_conjd_rec, n_all_data);
  
  for(i=0; i<n_conjd_rec; i++)
  {
    Feconjd_rec[i] = sqrt(ASmat_rec[i*n_conjd_rec + i] - PEmat2[i*n_conjd_rec + i]);
  }

  return;
}

