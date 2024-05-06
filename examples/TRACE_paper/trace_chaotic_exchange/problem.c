/**
 * Chaotic Exchange Orbit, Figure 2 in Lu, Hernandez & Rein (2024)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double j_init; // initial Jacobi constant
double tmax = 5000. * 11.86 * 2. * M_PI;

double jacobi_dh(struct reb_simulation* r){
  // Jacobi constant for the restricted 3-body problem
  // r needs to be interpreted in the rotating frame
  struct reb_particle* sun = &r->particles[0];
  struct reb_particle* jup = &r->particles[1];
  struct reb_particle* test = &r->particles[2];

  struct reb_orbit o = reb_orbit_from_particle(r->G, *jup, *sun);

  struct reb_vec3d x = {.x = test->x, .y=test->y, .z=test->z};
  struct reb_vec3d v = {.x = test->vx, .y=test->vy, .z=test->vz};

  double kinetic = 0.5 * reb_vec3d_length_squared(v);

  double d0 = sqrt(reb_vec3d_length_squared((struct reb_vec3d){.x = test->x - sun->x, .y = test->y - sun->y, .z =test->z - sun->z}));
  double d1 = sqrt(reb_vec3d_length_squared((struct reb_vec3d){.x = test->x - jup->x, .y = test->y - jup->y, .z =test->z - jup->z}));

  double U = -((r->G * sun->m) / d0) - ((r->G * jup->m) / d1);

  double omega = sqrt((r->G * (sun->m + jup->m) / (o.a * o.a * o.a)));
  double L =  x.x * v.y - x.y * v.x;

  double jac = kinetic + U - omega * L;
  return jac;
}


int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    r->dt = (8./365.) * 2 * M_PI;                // 8 days
    r->integrator = REB_INTEGRATOR_TRACE;
    r->ri_trace.r_crit_hill = 4. * 1.21;         // Value used in Lu Hernandez & Rein (2024)
    
    r->N_active = 2;

    // Initialize masses
    struct reb_particle star = {0};
    star.m = 1;
    struct reb_particle jup = {0};
    jup.m = 0.01 / (star.m - 0.01);

    // velocities
    double a = 5.2;
    double e = 0.0;
    star.x = 0.0;
    star.vy = 0.0;
    reb_simulation_add(r, star);

    jup.x = a;
    jup.vy = sqrt((r->G * (star.m + jup.m) / a) * ((1 - e) / (1 + e)));
    reb_simulation_add(r, jup);

    // Test particle
    struct reb_particle test = {0};
    double xhel = 4.42;
    double vhel = 0.0072 * (365.25) * (1 / (2 * M_PI)); // days to REBOUND years

    test.x = xhel + star.x;
    test.vy = vhel + star.vy;
    reb_simulation_add(r, test);

    reb_simulation_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    j_init = jacobi_dh(r);

    r->heartbeat  = heartbeat;
    system("rm -rf energy.txt");

    clock_t begin = clock();
    reb_simulation_integrate(r, tmax);
    clock_t end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    
    printf("\n%e %e\n", fabs((jacobi_dh(r) - j_init)/j_init), time_spent);
    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 10.*2.*M_PI)){
        reb_simulation_output_timing(r, tmax);
    }
    if (reb_simulation_output_check(r, 100. * 2.*M_PI)){
        // Once per 100 years, output the relative energy error to a text file
        FILE* f = fopen("energy.txt","a");

        double j = jacobi_dh(r);
        fprintf(f,"%e,%e\n",r->t, (j-j_init)/j_init);
        fclose(f);
    }
}
