#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

int main(int argc, char **argv)
{
  double t0=0.0, t1=0.0, dt;
  int opt, flag_help=0, flag_end=0;
  /* initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);
  
  if(thistask == roottask)
  {
    t0 = second();
    printf("===============BRAINS==================\n");
    printf("Starts to run...\n");
    printf("%d cores used.\n", totaltask);
  }
  /* cope with command options. */
  if(thistask == roottask)
  {
    opterr = 0; /* reset getopt. */
    optind = 0; /* reset getopt. */
    
    parset.flag_con_sys_err = 0;
    parset.flag_line_sys_err = 0;
    parset.flag_radio_sys_err = 0;
    parset.flag_postprc = 0;
    parset.flag_temp = 0;
    parset.flag_restart = 0;
    parset.flag_rng_seed = 0;

    
    while( (opt = getopt(argc, argv, "pt:rs:h")) != -1)
    {
      switch(opt)
      {
        case 'p':  /* only do postprocessing */
          parset.flag_postprc = 1;
          parset.temperature = 1.0;
          printf("# MCMC samples available, only do post-processing.\n");
          break;
        case 't': /* temperature for postprocessing */
          parset.flag_temp = 1;
          parset.temperature = atof(optarg);
          printf("# Set a temperature %f.\n", parset.temperature);
          if(parset.temperature == 0.0)
          {
            printf("# Incorrect option -t %s.\n", optarg);
            exit(0);
          }
          if(parset.temperature < 1.0)
          {
            printf("# Temperature should >= 1.0\n");
            exit(0);
          }
          break;

        case 'r':   /* restart from restored file */
          parset.flag_restart = 1;
          printf("# Restart run.\n");
          break;

        case 's':  /* set random number generator seed */
          parset.flag_rng_seed = 1;
          parset.rng_seed = atoi(optarg);
          printf("# Set random seed %d.\n", parset.rng_seed);
          break;

        case 'h':  /* print help */
          flag_help = 1;
          print_help();
          break;

        case '?':
          printf("# Incorrect option -%c %s.\n", optopt, optarg);
          exit(0);
          break;

        default:
          break;
      }
    }
    
    if(parset.flag_postprc == 1)
      parset.flag_restart = 0;
    
    if(flag_help == 0) // not only print help.
    {
      if(argv[optind] != NULL) // parameter file is specified 
        strcpy(parset.param_file, argv[optind]); /* copy input parameter file */
      else
      {
        flag_end = 1;
        fprintf(stderr, "# Error: No parameter file specified!\n");
      }
    }
    
  }
  
  MPI_Bcast(&flag_help, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(&flag_end, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  
  if(flag_end == 1 && flag_help ==0 )
  {
    if(thistask == roottask)
    {
      fprintf(stdout, "Ends incorrectly.\n");
    }

    MPI_Finalize();
    return 0;
  }

  if(flag_help == 0)
  {
    begin_run();    /* implementation run */

    end_run();      /* end run */
  }

  MPI_Finalize();   /* clean up and finalize MPI */
  if(thistask == roottask)
  {
    int ht, mt;
    double st;
    t1 = second();
    dt = timediff(t0, t1);
    get_hms(dt, &ht, &mt, &st);
    printf("Time used: %dh %dm %fs.\n", ht, mt, st);
    printf("Ends successfully.\n");
    printf("===============BRAINS==================\n");
  }
  return 0;
}