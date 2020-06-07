#include "rebound.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXSHELLS 100

extern unsigned long rebd_drift[MAXSHELLS];
extern unsigned long rebd_viol1[MAXSHELLS];
extern unsigned long rebd_viol2[MAXSHELLS];
double E0;

void heartbeat(struct reb_simulation* r){
    for(int i=0;i<r->ri_mercurana.Nmaxshells;i++){
        if (r->ri_mercurana.shellN_encounter){
        printf("%2d dom=%4d sub=%4d enc=%4d ",i, r->ri_mercurana.shellN_dominant[i], r->ri_mercurana.shellN_subdominant[i], r->ri_mercurana.shellN_encounter[i]);
        printf(" sub-passive=%4d enc-passive=%4d ", r->ri_mercurana.shellN_subdominant_passive[i], r->ri_mercurana.shellN_encounter_passive[i]);
        printf("%12lu ", rebd_drift[i]);
        printf("%12lu ", rebd_viol1[i]);
        printf("%12lu ", rebd_viol2[i]);
        printf("\n");
        }
    }
    printf("-------------\n");
    printf("maxshells %d\n", r->ri_mercurana.Nmaxshellsused);
    double E1 = reb_tools_energy(r);
    printf("dE/E= %e (offset=%e)  N= %d  N_active= %d  moved= %d\n",fabs((E0-E1)/E0), r->energy_offset, r->N, r->N_active, r->ri_mercurana.moved_particles);
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation_from_binary("out.tmp");
    r->exact_finish_time = 0;
    r->dt = 0.040;
    r->t = 0.0;
    r->heartbeat = heartbeat;
    r->testparticle_type = 1;
    r->integrator = REB_INTEGRATOR_MERCURANA;
    r->ri_mercurana.kappa = 1e-4;
    r->ri_mercurana.N_dominant = 1;
    r->ri_mercurana.Nmaxshells = 30;
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    for (int i=18;i<27;i++){
        reb_remove(r,i,1);
    }
    for (int i=3;i<27;i++){
        reb_remove(r,i,1);
    }
    for (int i=3;i<4;i++){
    }

    E0 = reb_tools_energy(r);
    printf("moved= %d\n", r->ri_mercurana.moved_particles);
    reb_step(r);
    double E1 = reb_tools_energy(r);
    printf("dE/E = %e\n",fabs((E0-E1)/E0));
    printf("moved= %d\n", r->ri_mercurana.moved_particles);
}

