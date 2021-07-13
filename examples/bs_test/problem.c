/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
 */
#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m", 1.);                // Central object
    reb_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
    reb_move_to_com(r);
    double E0 = reb_tools_energy(r);

    r->integrator = REB_INTEGRATOR_BS;
    r->ri_bs.scalAbsoluteTolerance = 1e-14;
    r->ri_bs.scalRelativeTolerance = 1e-14;
    r->dt = 1e-4;

    for(int i=0; i<100;i ++){
        printf("==============================\n");
        reb_step(r);
        double E1 = reb_tools_energy(r);
        printf("t  = %20.16f     DE/E = %.3e\n",r->t, fabs((E1-E0)/E0));
        printf("dt = %20.16f  \n",r->dt);
    }

    reb_free_simulation(r);
}

