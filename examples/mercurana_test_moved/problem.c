#include "rebound.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double E0;

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();
    srand(5);
    r->exact_finish_time = 0;
    r->dt = 0.01;
    r->integrator = REB_INTEGRATOR_MERCURANA;
    r->ri_mercurana.kappa = 1e-5;
    r->ri_mercurana.N_dominant = 0;
    r->ri_mercurana.Nmaxshells = 30;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    //r->usleep = 100;

    struct reb_particle p1 = {0}; 
    p1.m = 1.;
    p1.r = 1e-2;
    reb_add(r, p1);  

    struct reb_particle p3 =  reb_tools_orbit2d_to_particle(1.,p1,1,4.,0.,0,0);      
    p3.r = 1e-2;
    reb_add(r, p3); 

    r->N_active = r->N;

    int Np = 1000;
    for(int i=0;i<1000;i++){
    struct reb_particle p2 =  reb_tools_orbit2d_to_particle(1.,p1,0,1.,0.49,((double)i)/(double)Np*2.*M_PI,-0.21);      
    reb_add(r, p2); 
    }
    
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    reb_integrate(r, 1000);
    


    double E1 = reb_tools_energy(r);
    printf("dE/E = %e\n",fabs((E0-E1)/E0));
    printf("dcrit = %.20f\n",r->ri_mercurana.dcrit[0][0]);
    printf("dcrit = %.20f\n",r->ri_mercurana.dcrit[0][0]);
    printf("xyz %.20f %.20f %.20f\n",r->particles[0].x,r->particles[0].y,r->particles[0].z);
    printf("vxyz %.20f %.20f %.20f\n",r->particles[0].vx,r->particles[0].vy,r->particles[0].vz);
    printf("xyz %.20f %.20f %.20f\n",r->particles[1].x,r->particles[1].y,r->particles[1].z);
    printf("vxyz %.20f %.20f %.20f\n",r->particles[1].vx,r->particles[1].vy,r->particles[1].vz);
    printf("xyz %.20f %.20f %.20f\n",r->particles[2].x,r->particles[2].y,r->particles[2].z);
    printf("vxyz %.20f %.20f %.20f\n",r->particles[2].vx,r->particles[2].vy,r->particles[2].vz);
    printf("shellsused  %d\n",r->ri_mercurana.Nmaxshellsused);
    printf("moved %d\n",r->ri_mercurana.moved_particles);
}

