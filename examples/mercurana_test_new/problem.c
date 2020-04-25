#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

void heartbeat(struct reb_simulation* r){
    printf("%e    %e %e    %e %e  \n",r->t, r->particles[0].x, r->particles[0].y, r->particles[0].vx, r->particles[0].vy);
    printf("%e    %e %e    %e %e  \n",r->t, r->particles[1].x, r->particles[1].y, r->particles[1].vx, r->particles[1].vy);
    printf("-------------\n");
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();
    r->dt = 0.1;
    r->heartbeat = heartbeat;
    r->integrator = REB_INTEGRATOR_MERCURANA;

    struct reb_particle p1 = {0}; 
    p1.m = 1.;
    reb_add(r, p1);  
    
    struct reb_particle p2 = {0};
    p2.x = 1;
    p2.vy = 1;
    p2.m = 1e-3;
    reb_add(r, p2); 

    reb_integrate(r,10.);
}

