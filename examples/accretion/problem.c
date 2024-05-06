/**
 * Accretion of the Moon. Figure 7 in Lu, Hernandez & Rein (2024)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
int Nparticles = 1000;
double tmax = 6 * M_PI;

// Save Files
char title[100] = "energy.txt";
char title_remove[100] = "rm -rf energy.txt";

char remove_snapshots[100] = "rm -rf *snapshots_*";
char snapshot_1[100] = "snapshot_1.txt";
char snapshot_2[100] = "snapshot_2.txt";
char snapshot_3[100] = "snapshot_3.txt";
char snapshot_4[100] = "snapshot_4.txt";

double snap1_time = 0.0;
double snap2_time = 10.0 * 2 * M_PI;
double snap3_time = 500. * 2 * M_PI;
int snap1=1;
int snap2=1;
int snap3=1;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();
    
    r->integrator = REB_INTEGRATOR_TRACE;
    r->dt = 0.1;
    r->softening = 3e-8;
    r->rand_seed = 1;
    r->heartbeat = heartbeat;
   
    // Collisions
    r->collision            = REB_COLLISION_DIRECT;
    r->collision_resolve    = reb_collision_resolve_merge;

    // Initialize masses
    struct reb_particle earth = {0};
    earth.m = 1;
    double earth_r = 1./2.9; // in units of Roche radius
    earth.r = earth_r;
    reb_simulation_add(r, earth);

    double lunar_mass = 0.0123;

    // Add Disk Particles
    for (unsigned int i = 0; i < Nparticles; i++){
      double m = reb_random_powerlaw(r, 3.2e-7, 3.2e-4, -1.);
      double rad = pow(m, 1./3.) * 1.185 * earth_r;
      double a = reb_random_powerlaw(r, earth_r, 1.5, -1.);
      double e = reb_random_uniform(r, 0., 0.95);
      double inc = reb_random_uniform(r, 0, 50. * M_PI / 180.);
      double Omega = reb_random_uniform(r, 0, 2 * M_PI);
      double omega = reb_random_uniform(r, 0, 2 * M_PI);
      double f = reb_random_uniform(r, 0, 2 * M_PI);
      reb_simulation_add_fmt(r, "primary m r a e inc Omega omega f", earth, m, rad, a, e, inc, Omega, omega, f);
    }
    
    reb_simulation_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    
    system(title_remove);
    system(remove_snapshots);



    FILE* f = fopen(snapshot_1,"a");
    struct reb_particle* e = &r->particles[0];
    for (unsigned int i = 1; i < r->N; i++){
      struct reb_particle* p = &r->particles[i];

      double dx = p->x - e->x;
      double dy = p->y - e->y;
      double dz = p->z - e->z;

      double r = sqrt(dx*dx+dy*dy);
      fprintf(f, "%f,%f,%f\n",p->m,r,dz);
    }
    fclose(f);


    e_init = reb_simulation_energy(r);
    clock_t begin = clock();
    reb_simulation_integrate(r, tmax);
    clock_t end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

    FILE* f4 = fopen(snapshot_4,"a");
    e = &r->particles[0];
    for (unsigned int i = 1; i < r->N; i++){
      struct reb_particle* p = &r->particles[i];
      struct reb_orbit o = reb_orbit_from_particle(r->G, *p, *e);

      double dx = p->x - e->x;
      double dy = p->y - e->y;
      double dz = p->z - e->z;

      double r = sqrt(dx*dx+dy*dy);
      fprintf(f4, "%f,%f,%f,%f,%f,%f\n",p->m,r,dz,o.a,o.e,o.inc);
    }
    fclose(f4);

    reb_simulation_free(r);
}

void heartbeat(struct reb_simulation* r){
  if (reb_simulation_output_check(r, 0.01*2.*M_PI)){
      reb_simulation_output_timing(r, tmax);
  }

  // Take Snapshots
  if (snap2 && r->t > snap2_time){
    struct reb_particle* e = &r->particles[0];
    snap2 = 0;
    FILE* f = fopen(snapshot_2,"a");
    for (unsigned int i = 1; i < r->N; i++){
      struct reb_particle* p = &r->particles[i];

      double dx = p->x - e->x;
      double dy = p->y - e->y;
      double dz = p->z - e->z;

      double r = sqrt(dx*dx+dy*dy);
      fprintf(f, "%f,%f,%f\n",p->m,r,dz);
    }
    fclose(f);
  }

  if (snap3 && r->t > snap3_time){
    struct reb_particle* e = &r->particles[0];
    snap3 = 0;
    FILE* f = fopen(snapshot_3,"a");
    for (unsigned int i = 1; i < r->N; i++){
      struct reb_particle* p = &r->particles[i];

      double dx = p->x - e->x;
      double dy = p->y - e->y;
      double dz = p->z - e->z;

      double r = sqrt(dx*dx+dy*dy);
      fprintf(f, "%f,%f,%f\n",p->m,r,dz);
    }
    fclose(f);
  }

  // Track number of particles
  if (reb_simulation_output_check(r, 0.1)){
    clock_t end = clock();
    FILE* f = fopen(title,"a");
    fprintf(f, "%f,%e,%d\n", r->t, fabs((reb_simulation_energy(r) - e_init) / e_init), r->N-1);
    fclose(f);
  }
}
