#include "rebound.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MAXSHELLS 100

extern unsigned long rebd_drift[MAXSHELLS];
extern unsigned long rebd_viol1[MAXSHELLS];
extern unsigned long rebd_viol2[MAXSHELLS];
double E0;

#define Np 9
double J0[Np];

double jacobi(struct reb_simulation* r, unsigned int i){
    struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[1],r->particles[0]);
    double n = o.n;
    struct reb_particle p = r->particles[i];
    double dx0 = p.x - r->particles[0].x;
    double dy0 = p.y - r->particles[0].y;
    double dx1 = p.x - r->particles[1].x;
    double dy1 = p.y - r->particles[1].y;
    double p0d = sqrt(dx0*dx0+dy0*dy0);
    double p1d = sqrt(dx1*dx1+dy1*dy1);
    return (p.vx*p.vx + p.vy*p.vy)/2. - r->particles[0].m/p0d - r->particles[1].m/p1d - n*(p.x*p.vy - p.y*p.vx);
}

void heartbeat(struct reb_simulation* r){
    //if (r->steps_done==4) exit(0);
    for(int s=0;s<r->ri_mercurana.Nmaxshells;s++){
        printf("%2d", s);
        if (s+1==r->ri_mercurana.Nmaxshellsused){
            printf("* "); 
        }else{
            printf("  ");
        }
        for (int ptype=0; ptype<5; ptype++){
            if (r->ri_mercurana.pisd[ptype].shellN){
                printf("pt(%d)=%04d ", ptype, r->ri_mercurana.pisd[ptype].shellN[s]);
            }
        }
        //printf("%12lu ", rebd_drift[i]);
        //printf("%12lu ", rebd_viol1[i]);
        //printf("%12lu ", rebd_viol2[i]);
        printf("\n");
    }
    printf("-------------\n");
    double E1 = reb_tools_energy(r);
    printf("dE/E= %e (offset=%e)  N= %d  N_active= %d  moved= %d\n",fabs((E0-E1)/E0), r->energy_offset, r->N, r->N_active, r->ri_mercurana.Nmoved);
    double dJ_max = 0;
    int maxid = -1;
    for (int i=2;i<r->N;i++){
        double J = jacobi(r,i);
        double dJ = fabs((J-J0[i-2])/J0[i-2]);
        if (dJ>dJ_max) { 
            dJ_max = dJ;
            maxid = i;
        }
    }
    printf("dJ/J = %e  i=%d   steps=%lld\n",dJ_max,maxid,r->steps_done);
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();
    srand(1);
    r->heartbeat = heartbeat;
    r->exact_finish_time = 0;
    r->dt = 0.129;
    r->integrator = REB_INTEGRATOR_MERCURANA;
    r->ri_mercurana.kappa = 1e-4;
    //r->ri_mercurana.Gm0r0 = 1;
    r->ri_mercurana.N_dominant = 2;
    r->ri_mercurana.Nmaxshells = 30;
    //r->collision = REB_COLLISION_DIRECT;
    //r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    //r->usleep = 10000;

    struct reb_particle p1 = {0}; 
    p1.m = 1.;
    p1.r = 1e-2;
    reb_add(r, p1);  

    struct reb_particle p3 =  reb_tools_orbit2d_to_particle(1.,p1,1,4.,0.,0,0);      
    p3.r = 1e-2;
    reb_add(r, p3); 
    
    //struct reb_particle p4 =  reb_tools_orbit2d_to_particle(1.,p1,1e-3,0.3,0.,0,0);      
    //p4.r = 1e-3;
    //reb_add(r, p4); 
    
    //struct reb_particle p5 =  reb_tools_orbit2d_to_particle(1.,p3,1e-3,0.3,0.,0,0);      
    //p5.r = 1e-3;
    //reb_add(r, p5); 

    r->N_active = r->N;

    for(int i=0;i<Np;i++){
        if (i+2!=5) continue;
        struct reb_particle p2 =  reb_tools_orbit2d_to_particle(1.,p1,0,1.,0.49,((double)i)/(double)Np*2.*M_PI,-0.21);           
        reb_add(r, p2); 
        
    }
   
    reb_move_to_com(r);
    for(int i=0;i<Np;i++){
        J0[i] = jacobi(r,i+2); 
    }
    E0 = reb_tools_energy(r);
    reb_integrate(r, 1000000);
    


    double E1 = reb_tools_energy(r);
    printf("dE/E = %e\n",fabs((E0-E1)/E0));
}

