/*
 * Decomp
 * 
 * Decompose the optical emissions from disk and jet in 3C 273 with reveberation mapping data.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 * Apr 29, 2019
 */

#ifndef _ALLVARS_H
#define _ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define GRAVITY       6.672e-8
#define SOLAR_MASS    1.989e33
#define C             2.9979e10
#define SEC_PER_YEAR  3.155e7
#define CM_PER_LD    (C*8.64e4)

#define PI            M_PI
#define BRAINS_MAX_STR_LENGTH  (100)

/* variables for MPICH */
extern int thistask, totaltask, namelen;
extern int roottask;
extern char proc_name[MPI_MAX_PROCESSOR_NAME];

/* data */ 
extern int n_con_data, n_line_data, n_radio_data, n_all_data;
extern double *Tcon_data, *Fcon_data, *Fcerrs_data;
extern double *Tline_data, *Fline_data, *Flerrs_data;
extern double *Tradio_data, *Fradio_data, *Frerrs_data;
extern double *Fall_data, *Fall_rec;
extern double *workspace;
extern int *workspace_ipiv;

extern int n_con_rec, n_line_rec, n_radio_rec, n_all_rec, n_conjd_rec;
extern double *Tall_rec, *Fall_rec, *Feall_rec;
extern double *Tconjd_rec, *Fconjd_rec, *Feconjd_rec;

extern double con_error_mean, line_error_mean, radio_error_mean;
extern double con_scale, line_scale, radio_scale;

extern double *Larr_data, *Larr_rec, *Larr_conjd;
extern double *PCmat_data, *IPCmat_data, *USmat_rec, *USmatT_rec, *ASmat_rec;
extern double *Tmat1, *Tmat2, *PEmat1, *PEmat2;

extern double *mean;

/* MCMC */
extern int num_params;
extern double **par_range_model;
extern int *par_fix, npar_fix;
extern double *par_fix_val;

extern char dnest_options_file[BRAINS_MAX_STR_LENGTH];

/*! \struct PARSET
 *  \brief the configuration parameters.
 */
typedef struct 
{
  char param_file[BRAINS_MAX_STR_LENGTH];
  char continuum_file[BRAINS_MAX_STR_LENGTH], 
       line_file[BRAINS_MAX_STR_LENGTH], 
       radio_file[BRAINS_MAX_STR_LENGTH];
       
  char file_dir[BRAINS_MAX_STR_LENGTH];
  
  int flag_temp;
  double temperature;

  int flag_restart, flag_postprc;
  int flag_rng_seed, rng_seed;
  int flag_con_sys_err, flag_line_sys_err, flag_radio_sys_err;

  double lag_line_low, lag_line_upp, lag_radio_low, lag_radio_upp;

}PARSET;
extern PARSET parset;

#endif