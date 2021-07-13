/**
 * A very simple test problem
 * 
 * We first create a REBOUND simulation, then we add 
 * two particles and integrate the system for 100 time 
 * units.
 */
#include "rebound.h"
#include "integrator_bs.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

void computeDerivatives_ho(double* const yDot, const double* const y, double const t, void * ref){
    yDot[0] = y[1];
    yDot[1] = -y[0];
}

int main(int argc, char* argv[]) {
    struct reb_simulation_integrator_bs ri_bs = {0};
    reb_integrator_bs_reset_struct(&ri_bs);


    ri_bs.scalAbsoluteTolerance= 1e-8;
    ri_bs.scalRelativeTolerance= 1e-8;

    ri_bs.state.y = malloc(sizeof(double)*2);
    ri_bs.state.y[0] = 0;
    ri_bs.state.y[1] = 1;
    ri_bs.state.length = 2;
    ri_bs.hNew = 1e-3;

    ri_bs.state.derivatives = computeDerivatives_ho;
    FILE* f = fopen("output.txt", "w");
    for(int i=0;i<100;i++){
        double t = ri_bs.state.t;
        double sol = sin(t);
        fprintf(f,"%.10e  %.10e %.20f %.20f\n", t, ri_bs.hNew, ri_bs.state.y[0], sol);
        reb_integrator_bs_step(&ri_bs);

    }
    fclose(f);

    reb_integrator_bs_reset_struct(&ri_bs);
}

