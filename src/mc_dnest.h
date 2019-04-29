/*
 * Decomp
 * 
 * Decompose the optical emissions from disk and jet in 3C 273 with reveberation mapping data.
 * 
 * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
 * 
 * Apr 29, 2019
 */

#ifndef _MC_DNEST_CON_H

#define _MC_DNEST_CON_H

#include <stdbool.h>

/* functions */
void from_prior_mc(void *model);
void print_particle_mc(FILE *fp, const void *model);
void read_particle_mc(FILE *fp, void *model);
void restart_action_mc(int iflag);
void accept_action_mc();
void kill_action_mc(int i, int i_copy);
double log_likelihoods_cal_mc(const void *model);
double log_likelihoods_cal_initial_mc(const void *model);
double log_likelihoods_cal_restart_mc(const void *model);
double perturb_mc(void *model);

#endif
