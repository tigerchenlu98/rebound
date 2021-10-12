#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

// Right hand side for Harmonic Oscillator written as a coupled set of first order equations
// (f_0)' = f_1
// (f_1)' = - omega**2 * f_0
void rhs(double t, double* f, double* fp){
    const double omega = 0.2;
    fp[0] = f[1];
    fp[1] = -f[0] * omega*omega;
}

void dummy_forces(struct reb_simulation* r){
    struct reb_particle* dummy1 = &(r->particles[2]);
    struct reb_particle* dummy2 = &(r->particles[3]);
    // Set the overall scale of the problem. 
    // Required for convergence criteria only.
    dummy1->x =  1.; 
    dummy1->y =  1.; 
    dummy1->z =  1.; 
    dummy2->x =  1.; 
    dummy2->y =  1.; 
    dummy2->z =  1.; 
    
    // Overwrite gravity acceleration
    dummy1->ax = 0.; 
    dummy1->ay = 0.;
    dummy1->az = 0.;
    dummy2->ax = 0.; 
    dummy2->ay = 0.;
    dummy2->az = 0.;

    // Add HO force
    double f[2];
    f[0] = dummy1->vx;
    f[1] = dummy2->vx;
    double fp[2];
    rhs(r->t, f, fp);

    dummy1->ax = fp[0]; 
    dummy2->ax = fp[1]; 
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();

    reb_add_fmt(r, "m", 1.);                // Central object
    reb_add_fmt(r, "m a e", 1e-3, 1., 0.1); // Jupiter mass planet
    reb_move_to_com(r);

    r->N_active = 2;
    reb_add_fmt(r, "m vx", 0., 1.);         // Dummy particle
                                            // Mass is not relevant
                                            // vx is set to 1. (this is the initial value of the ODE)
    reb_add_fmt(r, "m vx", 0., 0.);         // Dummy particle
                                            // Mass is not relevant
                                            // vx is set to 0. (this is the initial value of the ODE)

    r->additional_forces = dummy_forces;
    r->integrator = REB_INTEGRATOR_IAS15;
    r->force_is_velocity_dependent = 1;


   FILE* f = fopen("out.txt","w"); 

    for (int i=0; i<300; i++){
        reb_integrate(r,r->t+0.3);

        struct reb_particle p = r->particles[1];
        struct reb_particle dummy1 = r->particles[2];
        struct reb_particle dummy2 = r->particles[3];
        printf("t=%.4f\t planet = %6.3f %6.3f   \t HO = %6.3f %6.3f\n", r->t, p.x, p.y, dummy1.vx, dummy2.vx);
        fprintf(f, "%e %e %e\n", r->t, p.x, dummy1.vx);
    }
    fclose(f);

    reb_free_simulation(r);
}

