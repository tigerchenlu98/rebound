#include "rebound.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double E0;

void heartbeat(struct reb_simulation* r){
    //printf("%e    %e %e    %e %e  \n",r->t, r->particles[0].x, r->particles[0].y, r->particles[0].vx, r->particles[0].vy);
    //printf("%e    %e %e    %e %e  \n",r->t, r->particles[1].x, r->particles[1].y, r->particles[1].vx, r->particles[1].vy);
    for(int i=0;i<r->ri_mercurana.Nmaxshells;i++){
        if (r->ri_mercurana.shellN_encounter){
        printf("%2d dom=%4d sub=%4d enc=%4d\n",i, r->ri_mercurana.shellN_dominant[i], r->ri_mercurana.shellN_subdominant[i], r->ri_mercurana.shellN_encounter[i]);
        }
    }
    printf("-------------\n");
    printf("maxshells %d\n", r->ri_mercurana.Nmaxshellsused);
    double E1 = reb_tools_energy(r);
    printf("dE/E = %e\n",fabs((E0-E1)/E0));
    printf("-------------\n");
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();
    r->exact_finish_time = 0;
    r->dt = 0.012;
    r->heartbeat = heartbeat;
    r->integrator = REB_INTEGRATOR_MERCURANA;
    r->ri_mercurana.kappa = 1e-4;
    r->ri_mercurana.N_dominant = 1;
    r->ri_mercurana.Nmaxshells = 20;
   
    if (0){ 
        struct reb_particle p1 = {0}; 
        p1.m = 1.;
        reb_add(r, p1);  
        
        struct reb_particle p2 = {0};
        p2.x = 1;
        p2.vy = 1;
        p2.m = 1e-3;
        reb_add(r, p2); 
        
        struct reb_particle p3 = {0};
        p3.x = 2;
        p3.vy = 0.74;
        p3.m = 1e-3;
        reb_add(r, p3); 
    }

    if (1){
        struct reb_particle p1 = {0}; 
        p1.m = 1.;
        reb_add(r, p1);  
        
        struct reb_particle p2 = {0};
        p2.x = 1;
        p2.vy = 1;
        p2.m = 1e-3;
        reb_add(r, p2); 
        
        struct reb_particle p3 = {0};
        p3.x = 2;
        p3.vy = 0.74;
        p3.m = 1e-3;
        reb_add(r, p3); 

        for (int i=0; i<100;i++){
            double a = reb_random_uniform(0.8,1.2);
            double omega = reb_random_uniform(0.,M_PI*2);
            double f = reb_random_uniform(0.,M_PI*2.);
            struct reb_particle p = reb_tools_orbit2d_to_particle(1.,p1,1e-6,a,0.1,omega,f);
            reb_add(r,p);
        }

        reb_move_to_com(r);
    }

    E0 = reb_tools_energy(r);
    reb_integrate(r,10000.);
    double E1 = reb_tools_energy(r);
    printf("dE/E = %e\n",fabs((E0-E1)/E0));
}

