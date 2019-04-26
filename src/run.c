#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

void begin_run()
{
  read_parset();

  read_data();
  scale_data();
  
  init();

  recon();
}

void end_run()
{
  free_memory_data();
  free_memory();
}

void init()
{
  int i, np;
  allocate_memory();

  /* setup Fall_data */
  memcpy(Fall_data, Fcon_data, n_con_data*sizeof(double));
  memcpy(Fall_data+n_con_data, Fline_data, n_line_data*sizeof(double));
  memcpy(Fall_data+n_con_data+n_line_data, Fradio_data, n_radio_data*sizeof(double));
  
  /* setup Larr_data */
  for(i=0; i<n_con_data; i++)
  {
    Larr_data[i*3 + 0] = 1.0;
    Larr_data[i*3 + 1] = 0.0;
    Larr_data[i*3 + 2] = 0.0;
  }
  np = n_con_data;
  for(i=0; i<n_line_data; i++)
  {
    Larr_data[(np+i)*3 + 0] = 0.0;
    Larr_data[(np+i)*3 + 1] = 1.0;
    Larr_data[(np+i)*3 + 2] = 0.0;
  }
  np = n_con_data + n_line_data;
  for(i=0; i<n_radio_data; i++)
  {
    Larr_data[(np+i)*3 + 0] = 0.0;
    Larr_data[(np+i)*3 + 1] = 0.0;
    Larr_data[(np+i)*3 + 2] = 1.0;
  }
}

void allocate_memory()
{
  workspace = malloc(n_all_data * 10 * sizeof(double));
  workspace_ipiv = malloc(n_all_data * 5 *sizeof(double));
  Fall_data = malloc(n_all_data * sizeof(double));
  Larr_data = malloc(n_all_data * 3 * sizeof(double));
  PCmat_data = malloc(n_all_data * n_all_data * sizeof(double));

  Tmat1 = malloc(n_all_data*n_all_data*sizeof(double));
  Tmat2 = malloc(n_all_data*n_all_data*sizeof(double));
}

void free_memory()
{
  free(workspace);
  free(workspace_ipiv);
  free(Fall_data);
  free(Larr_data);
  free(PCmat_data);
  free(Tmat1);
  free(Tmat2);
}