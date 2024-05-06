/**
 * Highly Eccentric Orbits with TRACE. Figure 3 in Lu Hernandez & Rein (2024)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "rebound.h"

void heartbeat(struct reb_simulation* r);

double e_init; // initial energy
double tmax = 300. * 2. * M_PI * 29.4;

int main(int argc, char* argv[]){

    // Initialize masses
    struct reb_simulation* r = reb_simulation_create();
    
    r->integrator = REB_INTEGRATOR_TRACE;
    r->dt = 0.15 * 2 * M_PI;
    r->heartbeat = heartbeat;
    r->exact_finish_time = 0;

    struct reb_particle star = {0};
    star.m = 1.0;

    // Jupiter
    double jm = 9.55e-4;
    double ja = 5.2;
    double je = 0.05;
    double jmse = 1.-je;

    // Saturn
    double sm = 2.857e-4;
    double sa = 9.58;
    double se = 0.99;
    double omse = 1. - se;
    double si = M_PI / 2.;
    
    reb_simulation_add(r, star);
    reb_simulation_add_fmt(r, "m a e", jm, ja, je);
    reb_simulation_add_fmt(r, "primary m a e inc omega f", star, sm, sa, se, si, M_PI/2., 0.0);

    reb_simulation_move_to_com(r);                // This makes sure the planetary systems stays within the computational domain and doesn't drift.
    e_init = reb_simulation_energy(r);
    system("rm -rf energy.txt");
    
    clock_t begin = clock();
    reb_simulation_integrate(r, tmax);
    clock_t end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    
    printf("%e\t%e\n", fabs((reb_simulation_energy(r) - e_init) / e_init), time_spent);
    reb_simulation_free(r);
}


void heartbeat(struct reb_simulation* r){
    if (reb_simulation_output_check(r, 10.*2.*M_PI)){
        reb_simulation_output_timing(r, tmax);
    }
    
    if (reb_simulation_output_check(r, 1. * 2.*M_PI)){
        // Once per year, output the relative energy error to a text file
        FILE* f = fopen("energy.txt","a");

        double e = fabs((reb_simulation_energy(r) - e_init) / e_init);
        fprintf(f,"%e,%e\n", r->t, e);
        fclose(f);
    }
}
