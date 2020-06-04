#include "rebound.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXSHELLS 100

extern unsigned long rebd_drift[MAXSHELLS];
extern unsigned long rebd_kick[MAXSHELLS];
extern unsigned long rebd_md1a[MAXSHELLS];
extern unsigned long rebd_md1b[MAXSHELLS];
extern unsigned long rebd_md1c[MAXSHELLS];
extern unsigned long rebd_md2a[MAXSHELLS];
extern unsigned long rebd_md2b[MAXSHELLS];
extern unsigned long rebd_md2c[MAXSHELLS];
extern unsigned long rebd_md3a[MAXSHELLS];
extern unsigned long rebd_md3b[MAXSHELLS];
extern unsigned long rebd_md3c[MAXSHELLS];
extern unsigned long rebd_md4a[MAXSHELLS];
extern unsigned long rebd_md4b[MAXSHELLS];
extern unsigned long rebd_md4c[MAXSHELLS];
double E0;

void heartbeat(struct reb_simulation* r){
    return;
    if (r->steps_done%10!=0) return;
    //printf("%e    %e %e    %e %e  \n",r->t, r->particles[0].x, r->particles[0].y, r->particles[0].vx, r->particles[0].vy);
    //printf("%e    %e %e    %e %e  \n",r->t, r->particles[1].x, r->particles[1].y, r->particles[1].vx, r->particles[1].vy);
    for(int i=0;i<r->ri_mercurana.Nmaxshells;i++){
        if (r->ri_mercurana.shellN_encounter){
        printf("%2d dom=%4d sub=%4d enc=%4d ",i, r->ri_mercurana.shellN_dominant[i], r->ri_mercurana.shellN_subdominant[i], r->ri_mercurana.shellN_encounter[i]);
        printf("%12lu ", rebd_drift[i]);
        printf("%12lu ", rebd_kick[i]);
        printf("%7lu ", rebd_md1a[i]);
        printf("%7lu ", rebd_md1b[i]);
        printf("%7lu ", rebd_md1c[i]);
        printf("%7lu ", rebd_md2a[i]);
        printf("%7lu ", rebd_md2b[i]);
        printf("%7lu ", rebd_md2c[i]);
        printf("%7lu ", rebd_md3a[i]);
        printf("%7lu ", rebd_md3b[i]);
        printf("%7lu ", rebd_md3c[i]);
        printf("%7lu ", rebd_md4a[i]);
        printf("%7lu ", rebd_md4b[i]);
        printf("%7lu ", rebd_md4c[i]);
        printf("\n");
        }
    }
    printf("-------------\n");
    printf("maxshells %d\n", r->ri_mercurana.Nmaxshellsused);
    printf("moved %d\n", r->ri_mercurana.moved_particles);
    double E1 = reb_tools_energy(r);
    printf("dE/E = %e\n",fabs((E0-E1)/E0));
    printf("N    = %d\n",r->N);
    printf("-------------\n");
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();
    srand(5);
    r->exact_finish_time = 0;
    r->dt = 0.02;
    r->heartbeat = heartbeat;
    r->integrator = REB_INTEGRATOR_MERCURANA;
    r->ri_mercurana.kappa = 1e-5;
    r->ri_mercurana.N_dominant = 1;
    r->ri_mercurana.Nmaxshells = 30;
    int rad = 1; 
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;

    if (0){ 
        struct reb_particle p1 = {0}; 
        p1.m = 1.;
        if (rad) p1.r = 0.0046524726;
        reb_add(r, p1);  
        
        struct reb_particle p2 = {0};
        p2.x = 1;
        p2.vy = 1;
        p2.m = 1e-3;
        if (rad) p2.r = 0.00046732617;
        reb_add(r, p2); 
        
        struct reb_particle p3 = {0};
        p3.x = 2;
        p3.vy = 0.74;
        p3.m = 1e-3;
        if (rad) p3.r = 0.00046732617;
        reb_add(r, p3); 
    }

    if (1){
        struct reb_particle p1 = {0}; 
        p1.m = 1.;
        if (rad) p1.r = 0.0046524726;
        reb_add(r, p1);  
        
        struct reb_particle p2 = {0};
        p2.x = 1;
        p2.vy = 1;
        p2.m = 1e-3;
        if (rad) p2.r = 0.00046732617;
        reb_add(r, p2); 
        
        struct reb_particle p3 = {0};
        p3.x = 2;
        p3.vy = 0.74;
        p3.m = 1e-3;
        if (rad) p3.r = 0.00046732617;
        reb_add(r, p3); 

        for (int i=0; i<100;i++){
            double a = reb_random_uniform(0.8,1.2);
            double omega = reb_random_uniform(0.,M_PI*2);
            double f = reb_random_uniform(0.,M_PI*2.);
            struct reb_particle p = reb_tools_orbit2d_to_particle(1.,p1,1e-6,a,0.1,omega,f);
            if (rad) p.r = 0.000046;
            reb_add(r,p);
        }

        reb_move_to_com(r);
    }


    E0 = reb_tools_energy(r);
    reb_integrate(r,10000.);
    double E1 = reb_tools_energy(r);
    printf("dE/E = %e\n",fabs((E0-E1)/E0));
}

