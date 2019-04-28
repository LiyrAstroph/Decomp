#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

/*!
 * read parameter set from parameter file.
 */
void read_parset()
{
  if(thistask == roottask)
  {
    #define MAXTAGS 300
    #define DOUBLE 1
    #define STRING 2
    #define INT 3

    FILE *fparam;
    int i, j, nt;
    char str[200], buf1[200], buf2[200], buf3[200];
    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];

    nt = 0;
    strcpy(tag[nt], "FileDir");
    addr[nt] = &parset.file_dir;
    id[nt++] = STRING;

    strcpy(tag[nt], "ContinuumFile");
    addr[nt] = &parset.continuum_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "LineFile");
    addr[nt] = &parset.line_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "RadioFile");
    addr[nt] = &parset.radio_file;
    id[nt++] = STRING;

    strcpy(tag[nt], "FlagConSysErr");
    addr[nt] = &parset.flag_con_sys_err;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagLineSysErr");
    addr[nt] = &parset.flag_line_sys_err;
    id[nt++] = INT;

    strcpy(tag[nt], "FlagRadioSysErr");
    addr[nt] = &parset.flag_radio_sys_err;
    id[nt++] = INT;

    char fname[200];
    sprintf(fname, "%s", parset.param_file);
    
    fparam = fopen(fname, "r");
    if(fparam == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    while(!feof(fparam))
    {
      sprintf(str,"empty");

      fgets(str, 200, fparam);
      if(sscanf(str, "%s%s%s", buf1, buf2, buf3)<2)
        continue;
      if(buf1[0]=='%')
        continue;
      for(i=0, j=-1; i<nt; i++)
        if(strcmp(buf1, tag[i]) == 0)
        {
          j = i;
          tag[i][0] = 0;
          //printf("%s %s\n", buf1, buf2);
          break;
        }
      if(j >=0)
      {
        switch(id[j])
        {
          case DOUBLE:
            *((double *) addr[j]) = atof(buf2);
            break;
          case STRING:
            strcpy(addr[j], buf2);
            break;
          case INT:
            *((int *)addr[j]) = (int) atof(buf2);
            break;
        }
      }
      else
      {
        fprintf(stderr, "# Error in file %s: Tag '%s' is not allowed or multiple defined.\n", 
                      parset.param_file, buf1);
        exit(0);
      }
    }
    fclose(fparam);
  }
  
  MPI_Bcast(&parset, sizeof(parset), MPI_BYTE, roottask, MPI_COMM_WORLD);
  return;
}

void read_data()
{
  FILE *fp;
  int i;
  char buf[200], fname[200];

  // first need to determine the number of data points 
  if(thistask == roottask)
  {
    int count;

    { 
      // continuum file
      sprintf(fname, "%s/%s", parset.file_dir, parset.continuum_file);
      fp = fopen(fname, "r");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      // count the number of lines
      count = 0;
      while(1)
      {
        fgets(buf, 200, fp);
        if(feof(fp)!=0)
          break;
        count++;
      }
      fclose(fp);
      n_con_data = count;
    
      printf("continuum data points: %d\n", n_con_data);
    }

    {
      sprintf(fname, "%s/%s", parset.file_dir, parset.line_file);
    // emission flux line
      fp = fopen(fname, "r");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }

      // count the number of lines
      count = 0;
      while(1)
      {
        fgets(buf, 200, fp);
        if(feof(fp)!=0)
          break;
        count++;
      }
      fclose(fp);
      n_line_data = count;
      printf("line data points: %d\n", n_line_data);
    }

    {
      sprintf(fname, "%s/%s", parset.file_dir, parset.radio_file);
      fp = fopen(fname, "r");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      // count the number of lines
      count = 0;
      while(1)
      {
        fgets(buf, 200, fp);
        if(feof(fp)!=0)
          break;
        count++;
      }
      fclose(fp);
      n_radio_data = count;
      printf("radio data points: %d\n", n_radio_data);
    }
  }

  MPI_Bcast(&n_con_data, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&n_line_data, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&n_radio_data, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  // now allocate memory for data
  allocate_memory_data();

  // now read data
  {
    if(thistask == roottask)
    {
      // continuum data
      sprintf(fname, "%s/%s", parset.file_dir, parset.continuum_file);
      fp = fopen(fname, "r");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      for(i=0; i<n_con_data; i++)
      {
        fscanf(fp, "%lf %lf %lf \n", &Tcon_data[i], &Fcon_data[i], &Fcerrs_data[i]);
      }
      fclose(fp);

      /* cal mean continuum error */
      con_error_mean = 0.0;
      for(i=0; i<n_con_data; i++)
      {
        con_error_mean += Fcerrs_data[i];
      }
      con_error_mean /= n_con_data;
    }

    MPI_Bcast(Tcon_data, n_con_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Fcon_data, n_con_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Fcerrs_data, n_con_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(&con_error_mean, 1, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  }

  {
    if(thistask == roottask)
    {
      // line data
      sprintf(fname, "%s/%s", parset.file_dir, parset.line_file);
      fp = fopen(fname, "r");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      for(i=0; i<n_line_data; i++)
      {
        fscanf(fp, "%lf %lf %lf \n", &Tline_data[i], &Fline_data[i], &Flerrs_data[i]);
      }
      fclose(fp);

      /* cal mean line error */
      line_error_mean = 0.0;
      for(i=0; i<n_line_data; i++)
      {
        line_error_mean += Flerrs_data[i];
      }
      line_error_mean /= n_line_data;
    }

    MPI_Bcast(Tline_data, n_line_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Fline_data, n_line_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Flerrs_data, n_line_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(&line_error_mean, 1, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  }

  {
    if(thistask == roottask)
    {
      // radio data
      sprintf(fname, "%s/%s", parset.file_dir, parset.radio_file);
      fp = fopen(fname, "r");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      for(i=0; i<n_radio_data; i++)
      {
        fscanf(fp, "%lf %lf %lf \n", &Tradio_data[i], &Fradio_data[i], &Frerrs_data[i]);
      }
      fclose(fp);

      /* cal mean radio error */
      radio_error_mean = 0.0;
      for(i=0; i<n_radio_data; i++)
      {
        radio_error_mean += Frerrs_data[i];
      }
      radio_error_mean /= n_radio_data;
    }

    MPI_Bcast(Tradio_data, n_radio_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Fradio_data, n_radio_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(Frerrs_data, n_radio_data, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    MPI_Bcast(&radio_error_mean, 1, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  }
  
  n_all_data = n_con_data + n_line_data + n_radio_data;
  
  if(thistask == roottask)
  {
    printf("con data: N=%d, E=%f\n", n_con_data, con_error_mean);
    printf("line data : N=%d, E=%f\n", n_line_data, line_error_mean);
    printf("radio data: N=%d, E=%f\n", n_radio_data, radio_error_mean);
  }
}

void scale_data()
{
  int i;
  double ave_con, ave_line, ave_radio;
  
  con_scale = 1.0;
  line_scale = 1.0;
  radio_scale = 1.0;

  ave_con = 0.0;
  ave_line = 0.0;
  ave_radio = 0.0;

  /* con */
  for(i=0; i<n_con_data; i++)
  {
    ave_con += Fcon_data[i];
  }

  ave_con /= n_con_data;
  con_scale = 1.0/ave_con;

  for(i=0; i<n_con_data; i++)
  {
    Fcon_data[i] *=con_scale;
    Fcerrs_data[i] *=con_scale;
  }
  con_error_mean *= con_scale;

  if(thistask == roottask)
    printf("task %d con scale: %e\t%e\n", thistask, con_scale, ave_con);
  
  /* line */
  for(i=0; i<n_line_data; i++)
  {
    ave_line += Fline_data[i];
  }
  ave_line /=n_line_data;
  line_scale = 1.0/ave_line;
  
  for(i=0; i<n_line_data; i++)
  {
    Fline_data[i] *= line_scale;
    Flerrs_data[i] *= line_scale;
  }
  line_error_mean *= line_scale;

  if(thistask == roottask)
    printf("task %d line scale: %e\t%e\n", thistask, line_scale, ave_line);

  /* radio */
  for(i=0; i<n_radio_data; i++)
  {
    ave_radio += Fradio_data[i];
  }
  ave_radio /=n_radio_data;
  radio_scale = 1.0/ave_radio;
  
  for(i=0; i<n_radio_data; i++)
  {
    Fradio_data[i] *= radio_scale;
    Frerrs_data[i] *= radio_scale;
  }
  radio_error_mean *= radio_scale;

  if(thistask == roottask)
    printf("task %d radio scale: %e\t%e\n", thistask, radio_scale, ave_radio);
    
  return;
}

void allocate_memory_data()
{

  Tcon_data = malloc(n_con_data * sizeof(double));
  Fcon_data = malloc(n_con_data * sizeof(double));
  Fcerrs_data = malloc(n_con_data * sizeof(double));

  Tline_data = malloc(n_line_data * sizeof(double));
  Fline_data = malloc(n_line_data * sizeof(double));
  Flerrs_data = malloc(n_line_data * sizeof(double));

  Tradio_data = malloc(n_radio_data * sizeof(double));
  Fradio_data = malloc(n_radio_data * sizeof(double));
  Frerrs_data = malloc(n_radio_data * sizeof(double));
  
  return;
}

void free_memory_data()
{
  free(Tcon_data);
  free(Fcon_data);
  free(Fcerrs_data);

  free(Tline_data);
  free(Fline_data);
  free(Flerrs_data);

  free(Tradio_data);
  free(Fradio_data);
  free(Frerrs_data);

  return;
}

/*!
 * get file name of posterior sample. 
 */
void get_posterior_sample_file(char *fname, char *samplefile)
{
  FILE *fp;
  char buf[BRAINS_MAX_STR_LENGTH], buf1[BRAINS_MAX_STR_LENGTH];
  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  buf[0]='#';
  while(buf[0]=='#')
  {
    fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
    if(sscanf(buf, "%s", buf1) < 1)  // a blank line
    {
      buf[0] = '#';
    }
  }
  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.new_level_interval);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.save_interval);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.thread_steps);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.max_num_levels);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%lf", &options.lambda);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%lf", &options.beta);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%d", &options.max_num_saves);

  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.sample_file);
  
  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.sample_info_file);
  
  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.levels_file);
  
  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
//  sscanf(buf, "%s", options.sampler_state_file);
  
  fgets(buf, BRAINS_MAX_STR_LENGTH, fp);
  sscanf(buf, "%s", samplefile);
  fclose(fp);
}
