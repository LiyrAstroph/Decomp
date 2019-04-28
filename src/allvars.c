#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"

/* variables for MPICH */
int thistask, totaltask, namelen;
int roottask = 0;
char proc_name[MPI_MAX_PROCESSOR_NAME];

/* data */
int n_con_data, n_line_data, n_radio_data, n_all_data;
double *Tcon_data, *Fcon_data, *Fcerrs_data;
double *Tline_data, *Fline_data, *Flerrs_data;
double *Tradio_data, *Fradio_data, *Frerrs_data;
double *Fall_data;
double *workspace;
int *workspace_ipiv;

int n_con_rec, n_line_rec, n_radio_rec, n_all_rec;
double *Tall_rec, *Fall_rec, *Feall_rec;

double con_error_mean, line_error_mean, radio_error_mean;
double con_scale, line_scale, radio_scale;

double *Larr_data, *Larr_rec;
double *PCmat_data, *IPCmat_data, *USmat_rec, *USmatT_rec, *ASmat_rec;
double *Tmat1, *Tmat2, *PEmat1, *PEmat2;

// MCMC
int num_params;
double **par_range_model;
int *par_fix, npar_fix;
double *par_fix_val;

char dnest_options_file[BRAINS_MAX_STR_LENGTH];

PARSET parset;