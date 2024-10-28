#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "rebound.h"
#include "particle.h"

// Code for 3-body scattering. Used to generate Figures 2 & 3 in Lu+ (2024)

double MASSC = 2.68 * 9.55e-4;
double tmax = 1e7*2*M_PI;
int nbodies = 3;
int ind;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    struct reb_particle star = {0};
    star.m = 0.81;
    reb_simulation_add(r, star);

    // Random seed passed as command line argument
    ind = 0;
    if (argc == 2){
      ind = atoi(argv[1]);
    }

    r->rand_seed = ind;

    double smas[3] = {8.759, 12.930, 19.086};
    double masses[3] = {MASSC, MASSC, MASSC};

    double rad_b = 4.36 * 4.26e-5;
    double rho = 1.33 * pow(1.496e13, 3.) / (1.989e33);
    double rad_c = pow(((3. * MASSC) / (4. * M_PI * rho)), 1./3.);

    double rads[3] = {rad_c, rad_c, rad_c};

    for (int i = 0; i < nbodies; i++){
      reb_simulation_add_fmt(r, "m a r e inc omega Omega f hash", masses[i], smas[i], rads[i], reb_random_uniform(r, 0.01, 0.05), reb_random_uniform(r, 0.0, 2.*M_PI/180.), reb_random_uniform(r, 0.0, 2.*M_PI), reb_random_uniform(r, 0.0, 2.*M_PI), reb_random_uniform(r, 0.0, 2.*M_PI), i+1);
    }

    reb_simulation_move_to_com(r);
    r->integrator = REB_INTEGRATOR_IAS15;
    r->exit_max_distance = 1e4;
    r->track_energy_offset = 1;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;

    double e_init = reb_simulation_energy(r);

    while (r->t < tmax){
       int retval = reb_simulation_integrate(r, tmax);
       if (retval == REB_STATUS_ESCAPE){
          double Ei = reb_simulation_energy(r);

          // Find and remove the particle
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
    reb_simulation_free(r);
}
