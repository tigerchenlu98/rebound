/**
 * Violent scattering system statistics with TRACE. Figures 4 and 5 in Lu, Hernandez & Rein (2024)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "rebound.h"

double e_start; // initial energy
double tmax = 1e7*2*M_PI;
int nbodies = 3;
int ind;

char title_stats[100] = "error_and_time_stats.txt";
char element_stats[100] = "orbital_elements.txt";

int main(int argc, char* argv[]){

    // Initialize masses

    struct reb_simulation* r = reb_simulation_create();
    struct reb_particle star = {0};
    star.m = 1;
    star.r = 0.005;
    reb_simulation_add(r, star);
    double planet_m = 9.55e-4;
    double planet_r = 0.000477895;

    double sma = 5.;
    double add = 0.;
    double delta= 3.;

    ind = 0;

    // The random seed is passed as a command line argument
    if (argc == 2){
      r->rand_seed = atoi(argv[1]);
    }

    double smas[3] = {};

    for (int i = 0; i < nbodies; i++){
      smas[i] = sma;
      reb_simulation_add_fmt(r, "m r a e inc hash", planet_m, planet_r, sma, 0.05, (double)i * M_PI / 180., i+1);
      double num = -pow(2., 1./3.) * pow(3., 1./3.) * sma - pow((planet_m / star.m), 1./3.) * delta * sma;
      double denom = -pow(2., 1./3.) * pow(3., 1./3.) + pow((planet_m / star.m), 1./3.) * delta;
      sma = num/denom;
    }

    // Shift outermost planet
    add = reb_random_uniform(r, -1e-12, 1e-12);
    struct reb_particle* pouter = &r->particles[nbodies];
    pouter->x += add;

    // Select Timestep 
    struct reb_particle* sun = &r->particles[0];
    double final_a = smas[0] * smas[1] * smas[2] / (smas[0]*smas[1] + smas[1]*smas[2] + smas[0]*smas[2]);
    double final_ts = sqrt(4 * M_PI * M_PI / (r->G * star.m) * final_a * final_a * final_a) / 15.12345;

    reb_simulation_move_to_com(r);

    r->integrator = REB_INTEGRATOR_TRACE;
    r->dt = final_ts;
    r->exact_finish_time = 0;
    r->exit_max_distance = 1e4;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;

    e_start = reb_simulation_energy(r);
    clock_t begin = clock();

    // Integrate
    while (r->t < tmax){
       int retval = reb_simulation_integrate(r, tmax);
       if (retval == REB_STATUS_ESCAPE){
          // Track energy offset
          double Ei = reb_simulation_energy(r);

          // Find and remove escaping  particles
          int remove_ind;
          for (int i = 1; i < nbodies+1; i++){

              struct reb_particle* p = reb_simulation_particle_by_hash(r, i);
              if (p != NULL){
                double dx = p->x;
                double dy = p->y;
                double dz = p->z;
                double d2 = dx*dx+dy*dy+dz*dz;

                if (d2>r->exit_max_distance * r->exit_max_distance){
                    remove_ind = i;
                }
              }
          }
          reb_simulation_remove_particle_by_hash(r, remove_ind, 1);
          reb_simulation_move_to_com(r);

          r->energy_offset += Ei - reb_simulation_energy(r);
       }
    }
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    FILE* tf = fopen(title_stats, "a");
    fprintf(tf, "%d,%d,%e,%e\n", ind, r->N-1, fabs((reb_simulation_energy(r) - e_start)/e_start), time_spent);
    fclose(tf);

    struct reb_particle* s = &r->particles[0];
    FILE* ef = fopen(element_stats, "a");
    fprintf(ef, "%d", ind);
    double tot_m = r->particles[0].m;
    for (unsigned int i = 1; i < nbodies+1; i++){
        struct reb_particle* p = reb_simulation_particle_by_hash(r, i);
        if (p != NULL){
          struct reb_orbit o = reb_orbit_from_particle(r->G, *p, *s);
          fprintf(ef, ",%f,%f,%f", o.a, o.e, o.inc);
          tot_m += p->m;
        }
        else{
          fprintf(ef, ", , , ");
        }
    }
    int num_collisions = (int) ((tot_m -1.)/ planet_m) - (r->N - 1);
    fprintf(ef, ",%d\n", num_collisions);
    fclose(ef);

    reb_simulation_free(r);
}
