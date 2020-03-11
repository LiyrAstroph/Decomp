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
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>
// header file for DNEST
#include "dnestvars.h"

#include "allvars.h"
#include "mc_dnest.h"
#include "proto.h"

DNestFptrSet *fptrset;

double mc_dnest(int argc, char **argv)
{
  int i;
  double logz;

  num_params = 13;

  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  
  fptrset = dnest_malloc_fptrset();

  /* setup functions used for dnest*/
  fptrset->from_prior = from_prior_mc;
  fptrset->perturb = perturb_mc;
  fptrset->print_particle = print_particle_mc;
  fptrset->restart_action = restart_action_mc;
  fptrset->accept_action = accept_action_mc;
  fptrset->kill_action = kill_action_mc;
  fptrset->read_particle = read_particle_mc;
  fptrset->log_likelihoods_cal = log_likelihoods_cal_mc;
  fptrset->log_likelihoods_cal_initial = log_likelihoods_cal_mc;
  fptrset->log_likelihoods_cal_restart = log_likelihoods_cal_mc;
  
  set_par_range();
  
  /* setup fixed parameters */
  for(i=0; i<num_params; i++)
    par_fix[i] = 0;

  /* fix systematic error of continuum */
  if(parset.flag_con_sys_err != 1)
  {
    par_fix[0] = 1;
    par_fix_val[0] = log(1.0);
  }
  /* fix systematic error of line */
  if(parset.flag_line_sys_err != 1)
  {
    par_fix[1] = 1;
    par_fix_val[1] = log(1.0);
  }
  /* fix systematic error of radio */
  if(parset.flag_radio_sys_err != 1)
  {
    par_fix[2] = 1;
    par_fix_val[2] = log(1.0);
  }

  logz = dnest(argc, argv, fptrset, num_params, dnest_options_file);
  
  dnest_free_fptrset(fptrset);
  return logz;
}

void set_par_range()
{
  int i = 0;

  /* systemactic errors, con, line, radio */
  par_range_model[i][0] = log(1.0);
  par_range_model[i++][1] = log(1.0+10.0);
  
  par_range_model[i][0] = log(1.0);
  par_range_model[i++][1] = log(1.0+10.0);

  par_range_model[i][0] = log(1.0);
  par_range_model[i++][1] = log(1.0+10.0);

  /* con variability */
  par_range_model[i][0] = log(1.0e-5);
  par_range_model[i++][1] = log(0.1);

  par_range_model[i][0] = log(1.0);
  par_range_model[i++][1] = log(Tcon_data[n_con_data-1] - Tcon_data[0]);

  /* radio variability */  
  par_range_model[i][0] = log(1.0e-5);
  par_range_model[i++][1] = log(0.1);

  par_range_model[i][0] = log(1.0);
  par_range_model[i++][1] = log(Tcon_data[n_con_data-1] - Tcon_data[0]);

  /* Hbeta transfer function */
  par_range_model[i][0] = log(1.0e-3);
  par_range_model[i++][1] = log(1.0e3);
  
  par_range_model[i][0] = 100.0;
  par_range_model[i++][1] = 300.0;

  par_range_model[i][0] = log(1.0);
  par_range_model[i++][1] = log(1.0e3);

  /* radio transfer function */
  par_range_model[i][0] = log(1.0e-3);
  par_range_model[i++][1] = log(1.0e3);
  
  par_range_model[i][0] = 200.0;
  par_range_model[i++][1] = 800.0;

  par_range_model[i][0] = log(1.0);
  par_range_model[i++][1] = log(1.0e3);

  return;
}

void from_prior_mc(void *model)
{
  int i;
  double *pm = (double *)model;
  
  for(i=0; i<num_params; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand()*(par_range_model[i][1] - par_range_model[i][0]);
  }

  for(i=0; i<num_params; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }
}

double perturb_mc(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, width, limit1, limit2;
  int which, which_level, which_level_update;
  
  /* sample variability parameters more frequently */
  do
  {
    which = dnest_rand_int(num_params);
  
  }while(par_fix[which] == 1);
  
  /*which_level_update = dnest_get_which_level_update();
  which_level = which_level_update > (size_levels - 10)?(size_levels -10):which_level_update;

  if( which_level > 0)
  {
    limit1 = limits[(which_level-1) * num_params *2 + which *2];
    limit2 = limits[(which_level-1) * num_params *2 + which *2 + 1];
    width = limit2 - limit1;
  }
  else
  {
    width = ( par_range_model[which][1] - par_range_model[which][0] );
  }*/

  width = ( par_range_model[which][1] - par_range_model[which][0] );
  pm[which] += dnest_randh() * width;
  dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  
  return logH;
}

/*!
 * This function calculate log likelihood probability.
 */
double log_likelihoods_cal_mc(const void *model)
{
  double logL;
  logL = prob(model);
  return logL;
}

/*!
 * This function print the particle into the file.
 */
void print_particle_mc(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%e ", pm[i] );
  }
  fprintf(fp, "\n");
}

/*!
 * This function read the particle from the file.
 */
void read_particle_mc(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < dnest_num_params; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "#Error: Cannot read file %s.\n", options.sample_file);
      exit(0);
    }
  }
  return;
}

void restart_action_mc(int iflag)
{

  return;
}

void accept_action_mc()
{
  return;
}

void kill_action_mc(int i, int i_copy)
{
  return;
}


