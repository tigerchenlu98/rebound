/**
 * @file    integrator_mercurana.c
 * @brief   MERCURANA is a symplectic multi-step method.
 *          It is adaptive, can handle close encounters and works in complex
 *          hierarchies.  
 * @author  Hanno Rein, Daniel Tamayo
 * 
 * @section LICENSE
 * Copyright (c) 2020 Hanno Rein, Daniel Tamayo
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"
#include "integrator.h"
#include "gravity.h"
#include "integrator_mercurana.h"
#include "integrator_eos.h"
#include "collision.h"
#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b


// DEBUG Variables:

#define MAXSHELLS 100

unsigned long rebd_drift[MAXSHELLS];
unsigned long rebd_viol1[MAXSHELLS];
unsigned long rebd_viol2[MAXSHELLS];



// Machine independent implementation of pow(*,1./3.) using Newton's method.
// Speed is not an issue. Only used to calculate dcrit.
static double sqrt3(double a){
    double x = 1.;
    for (int k=0; k<200;k++){  // A smaller number should be ok too.
        double x2 = x*x;
        x += (a/x2-x)/3.;
    }
    return x;
}

// Helper functions for L_infinity
static double f(double x){
    if (x<0) return 0;
    return exp(-1./x);
}

static double dfdy(double x){
    if (x<0) return 0;
    return exp(-1./x)/(x*x);
}

// Infinitely differentiable switching function.
static double reb_integrator_mercurana_L_infinity(const struct reb_simulation* const r, double d, double ri, double ro){
    double y = (d-ri)/(ro-ri);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return f(y) /(f(y) + f(1.-y));
    }
}

// First derivative of the infinitely differentiable switching function.
static double reb_integrator_mercurana_dLdr_infinity(const struct reb_simulation* const r, double d, double ri, double ro){
    double y = (d-ri)/(ro-ri);
    double dydr = 1./(ro-ri);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 0.;
    }else{
        return dydr*(
                dfdy(y) /(f(y) + f(1.-y))
                -f(y) /(f(y) + f(1.-y))/(f(y) + f(1.-y)) * (dfdy(y) - dfdy(1.-y))
                );
    }
}

// This function returns the closest approach between particles p1 and p2 during a drift step with length dt
static double reb_mercurana_predict_rmin2(struct reb_particle p1, struct reb_particle p2, double dt){ 
    double dts = copysign(1.,dt); 
    dt = fabs(dt);
    const double dx1 = p1.x - p2.x; // distance at beginning
    const double dy1 = p1.y - p2.y;
    const double dz1 = p1.z - p2.z;
    const double r1 = (dx1*dx1 + dy1*dy1 + dz1*dz1);
    const double dvx1 = dts*(p1.vx - p2.vx); 
    const double dvy1 = dts*(p1.vy - p2.vy);
    const double dvz1 = dts*(p1.vz - p2.vz);
    const double dx2 = dx1 +dt*dvx1; // distance at end
    const double dy2 = dy1 +dt*dvy1;
    const double dz2 = dz1 +dt*dvz1;
    const double r2 = (dx2*dx2 + dy2*dy2 + dz2*dz2);
    const double t_closest = -(dx1*dvx1 + dy1*dvy1 + dz1*dvz1)/(dvx1*dvx1 + dvy1*dvy1 + dvz1*dvz1);
    const double dx3 = dx1+t_closest*dvx1; // closest approach
    const double dy3 = dy1+t_closest*dvy1;
    const double dz3 = dz1+t_closest*dvz1;
    const double r3 = (dx3*dx3 + dy3*dy3 + dz3*dz3);

    double rmin2 = MIN(r1,r2);
    if (t_closest/dt>=0. && t_closest/dt<=1.){
        rmin2 = MIN(rmin2, r3);
    }
    return rmin2;
}


// This functions records a physical collision between particles i and j, to be resolved later.
static void reb_mercurana_record_collision(struct reb_simulation* const r, unsigned int i, unsigned int j){
    struct reb_simulation_integrator_mercurana* rim = &(r->ri_mercurana);
    if (r->collisions_allocatedN<=rim->collisions_N){
        // Allocate memory if there is no space in array.
        r->collisions_allocatedN = r->collisions_allocatedN ? r->collisions_allocatedN * 2 : 32;
        r->collisions = realloc(r->collisions,sizeof(struct reb_collision)*r->collisions_allocatedN);
    }
    r->collisions[rim->collisions_N].p1 = i;
    r->collisions[rim->collisions_N].p2 = j;
    struct reb_ghostbox gb = {0};
    gb.shiftx = r->particles[i].x;
    gb.shifty = r->particles[i].y;
    gb.shiftz = r->particles[i].z;
    gb.shiftvx = r->particles[i].vx;
    gb.shiftvy = r->particles[i].vy;
    gb.shiftvz = r->particles[i].vz;
    r->collisions[rim->collisions_N].gb = gb;
    rim->collisions_N++;
}

static double reb_drift_from_straight_line(struct reb_particle p0, double dt0, struct reb_particle p1, double dt1){
    double dx = p0.x + dt0*p0.vx - (p1.x + dt1*p1.vx); 
    double dy = p0.y + dt0*p0.vy - (p1.y + dt1*p1.vy); 
    double dz = p0.z + dt0*p0.vz - (p1.z + dt1*p1.vz); 
    return sqrt(dx*dx + dy*dy + dz*dz);
}



static inline enum REB_PTYPE get_ptype_drifting(enum REB_ITYPE itype){
    switch (itype){
        case REB_ITYPE_DOM_DOM:
        case REB_ITYPE_DOM_SUB:
        case REB_ITYPE_DOM_SUBP:
            return REB_PTYPE_DOM; 
        case REB_ITYPE_SUB_DOM:
            return REB_PTYPE_SUB; 
        case REB_ITYPE_ENC_ENC:
        case REB_ITYPE_ENC_ENCP:
            return REB_PTYPE_ENC; 
        case REB_ITYPE_SUBP_DOM:
            return REB_PTYPE_SUBP; 
        case REB_ITYPE_ENCP_ENC:
            return REB_PTYPE_ENCP; 
        default:
            return REB_PTYPE_NONE;
    }
}
static inline enum REB_PTYPE get_ptype_interacting(enum REB_ITYPE itype){
    switch (itype){
        case REB_ITYPE_DOM_DOM:
        case REB_ITYPE_SUB_DOM:
        case REB_ITYPE_SUBP_DOM:
            return REB_PTYPE_DOM; 
        case REB_ITYPE_DOM_SUB:
            return REB_PTYPE_SUB; 
        case REB_ITYPE_ENC_ENC:
        case REB_ITYPE_ENCP_ENC:
            return REB_PTYPE_ENC; 
        case REB_ITYPE_DOM_SUBP:
            return REB_PTYPE_SUBP; 
        case REB_ITYPE_ENC_ENCP:
            return REB_PTYPE_ENCP; 
        default:
            return REB_PTYPE_NONE;
    }
}


// Checks for if particle mi is getting close to any particle in pisd_interacting. 
// If so, move interacting particle one shell deeper.
// If not, update maxdrift of particle mi.
void check_one(struct reb_simulation* r, double dt, unsigned int shell, int mi, 
        struct reb_mdd* mdd, struct reb_pisd pisd_interacting){
    struct reb_particle* const particles = r->particles;
    const double* const dcrit = r->ri_mercurana.dcrit[shell];
    const double* const p_t = r->ri_mercurana.p_t;
    mdd->p0 = particles[mi];
    mdd->t0 = r->ri_mercurana.t_now;
    mdd->maxdrift = 1e300;
    for (int j=0; j<pisd_interacting.shellN[shell]; j++){
        int mj = pisd_interacting.map[shell][j]; 
        if (mi==mj) continue; // no self-interactions
        // Drift particle mj to current time for rmin calculation
        struct reb_particle pj = particles[mj];
        double pjdrift = r->ri_mercurana.t_now - p_t[mj];
        pj.x += pjdrift * pj.vx;
        pj.y += pjdrift * pj.vy;
        pj.z += pjdrift * pj.vz;
        double rmin2 = reb_mercurana_predict_rmin2(particles[mi],pj,dt);
        double rsum = particles[mi].r+particles[mj].r;
        if (rmin2< rsum*rsum && r->collision==REB_COLLISION_DIRECT){
            // TODO: Think about right location for collision check
            reb_mercurana_record_collision(r,mi,mj);
        }
        double dcritsum = dcrit[mi]+dcrit[mj];
        if (rmin2< dcritsum*dcritsum){ 
            if (pisd_interacting.inshell[mj] == shell){
                // TODO: check if particle needs to be drifted!
                pisd_interacting.inshell[mj] = shell+1;
                pisd_interacting.map[shell+1][pisd_interacting.shellN[shell+1]] = mj;
                pisd_interacting.shellN[shell+1]++;
            }
        }else{
            // Interaction potentially not resolved in subshells. Need maxdrift.
            const double maxdrift = (sqrt(rmin2) - dcritsum)/2.;
            mdd->maxdrift = MIN(mdd->maxdrift,maxdrift);
        }
    }
}

static void check_maxdrift_violation(struct reb_simulation* r, double dt, unsigned int shell, enum REB_ITYPE itype){
    enum REB_PTYPE drifting = get_ptype_drifting(itype);
    enum REB_PTYPE interacting = get_ptype_interacting(itype);

    struct reb_pisd pisd_drifting = r->ri_mercurana.pisd[drifting];
    struct reb_pisd pisd_interacting = r->ri_mercurana.pisd[interacting];

    struct reb_particle* const particles = r->particles;
    
    for (int i=0; i<pisd_drifting.shellN[shell]; i++){
        int mi = pisd_drifting.map[shell][i]; 
        for (int s=0;s<shell;s++){
            struct reb_mdd* mdd = &(r->ri_mercurana.mdd[itype][s][mi]);
            double dt0 = r->ri_mercurana.t_now - mdd->t0 + dt;
            // TODO Check what function does
            double drift = reb_drift_from_straight_line(mdd->p0,dt0,particles[mi],dt);
            if (drift > mdd->maxdrift){
                check_one(r, dt, s, mi, mdd, pisd_interacting);
            }
        }
    }
}


static void check_this_shell(struct reb_simulation* r, double dt, unsigned int shell, enum REB_ITYPE itype){
    enum REB_PTYPE drifting = get_ptype_drifting(itype);
    enum REB_PTYPE interacting = get_ptype_interacting(itype);

    struct reb_pisd pisd_drifting = r->ri_mercurana.pisd[drifting];
    struct reb_pisd pisd_interacting = r->ri_mercurana.pisd[interacting];

    for (int i=0; i<pisd_drifting.shellN[shell]; i++){
        int mi = pisd_drifting.map[shell][i]; 
        struct reb_mdd* mdd = &(r->ri_mercurana.mdd[itype][shell][mi]);
        check_one(r, dt, shell, mi, mdd, pisd_interacting);
    }
}


// This function checks if there are any close encounters or physical collisions between particles in a given shell during a drift step of length dt. If a close encounter occures, particles are placed in deeper shells.
static void reb_mercurana_encounter_predict(struct reb_simulation* const r, double dt, int shell){
    struct reb_simulation_integrator_mercurana* rim = &(r->ri_mercurana);
    
    const int N_active = r->N_active==-1?r->N:r->N_active;
    
    if (shell==0){
        // Setup maps in outermost shell 
        rim->shellN_dominant[0] = rim->N_dominant;
        rim->shellN_subdominant[0] = N_active - rim->N_dominant;
        rim->shellN_subdominant_passive[0] = r->N - N_active;
        rim->shellN_encounter[0] = N_active - rim->N_dominant;
        rim->shellN_encounter_passive[0] = r->N - N_active;
        for (int i=0; i<rim->N_dominant; i++){
            map_dominant[i] = i; 
        }
        for (int i=rim->N_dominant; i<N_active; i++){
            map_subdominant[i-rim->N_dominant] = i; 
            map_encounter[i-rim->N_dominant] = i; 
        }
        for (int i=N_active; i<r->N; i++){
            map_subdominant_passive[i-N_active] = i; 
            map_encounter_passive[i-N_active] = i; 
        }
    }
    // First, setup arrays for current shell.
    for (int i=0; i<rim->shellN_encounter[shell]; i++){
        int mi = map_encounter[i]; 
        inshell_encounter[mi] = shell;
    }
    for (int i=0; i<rim->shellN_dominant[shell]; i++){
        int mi = map_dominant[i]; 
        inshell_dominant[mi] = shell;
    }
    for (int i=0; i<rim->shellN_subdominant[shell]; i++){
        int mi = map_subdominant[i]; 
        inshell_subdominant[mi] = shell;
    }
    for (int i=0; i<rim->shellN_encounter_passive[shell]; i++){
        int mi = map_encounter_passive[i]; 
        inshell_encounter[mi] = shell;
    }
    for (int i=0; i<rim->shellN_subdominant_passive[shell]; i++){
        int mi = map_subdominant_passive[i]; 
        inshell_subdominant[mi] = shell;
    }
    
    rim->collisions_N = 0;
    
    if (shell+1>=rim->Nmaxshells){ // does sub-shell exist?
        return; // If not, no prediction needed.
    }

    rim->shellN_encounter[shell+1] = 0;
    rim->shellN_encounter_passive[shell+1] = 0;
    rim->shellN_dominant[shell+1] = 0;
    rim->shellN_subdominant[shell+1] = 0;
    rim->shellN_subdominant_passive[shell+1] = 0;

    if (shell!=0){
        for (int itype=0; itype<8; itype++){ 
            check_maxdrift_violation(r, dt, shell, itype);
        } 
    }

    for (int itype=0; itype<8; itype++){ 
        check_this_shell(r, dt, shell, itype);
    }
    
    if (rim->collisions_N){
        //printf("collision\n");
        unsigned int N_before = r->N;
        reb_collision_search(r); // will resolve collisions
        if (N_before!=r->N){
            // Need to redo predict step as particles changed.
            reb_mercurana_encounter_predict(r, dt, shell);
        }
        rim->collisions_N = 0;
    }
}

// Main Kernel Operator: Interaction step. 
// y = timestep for acceleration
// v = timestep for jerk (0 if not used)
static void reb_integrator_mercurana_interaction_step(struct reb_simulation* r, double y, double v, unsigned int shell){
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    struct reb_particle* const particles = r->particles;
    r->gravity = REB_GRAVITY_MERCURANA; // needed here again for SimulationArchive
    rim->current_shell = shell;
    reb_update_acceleration(r);
    if (v!=0.){
        reb_calculate_and_apply_jerk(r,v);
    }
    unsigned int* map_encounter = rim->map_encounter[shell];
    unsigned int* map_encounter_passive = rim->map_encounter_passive[shell];
    unsigned int* map_dominant = rim->map_dominant[shell];
    unsigned int* map_subdominant = rim->map_subdominant[shell];
    unsigned int* map_subdominant_passive = rim->map_subdominant_passive[shell];
    unsigned int* inshell_encounter = rim->inshell_encounter;

    for (int i=0;i<rim->shellN_dominant[shell];i++){ // Apply acceleration. Jerk already applied.
        const int mi = map_dominant[i];
        particles[mi].vx += y*particles[mi].ax;
        particles[mi].vy += y*particles[mi].ay;
        particles[mi].vz += y*particles[mi].az;
    }
    for (int i=0;i<rim->shellN_encounter[shell];i++){ // Apply acceleration. Jerk already applied.
        const int mi = map_encounter[i];
        particles[mi].vx += y*particles[mi].ax;
        particles[mi].vy += y*particles[mi].ay;
        particles[mi].vz += y*particles[mi].az;
    }
    for (int i=0;i<rim->shellN_encounter_passive[shell];i++){ // Apply acceleration. Jerk already applied.
        const int mi = map_encounter_passive[i];
        particles[mi].vx += y*particles[mi].ax;
        particles[mi].vy += y*particles[mi].ay;
        particles[mi].vz += y*particles[mi].az;
    }
    if (shell>0){ // All particles are encounter particles in shell 0, no need for subdominant kick
        for (int i=0;i<rim->shellN_subdominant[shell];i++){ // Apply acceleration. Jerk already applied.
            const int mi = map_subdominant[i];
            if (inshell_encounter[mi]<shell){ // do not apply acceleration twice
                particles[mi].vx += y*particles[mi].ax;
                particles[mi].vy += y*particles[mi].ay;
                particles[mi].vz += y*particles[mi].az;
            }
        }
        for (int i=0;i<rim->shellN_subdominant_passive[shell];i++){ // Apply acceleration. Jerk already applied.
            const int mi = map_subdominant_passive[i];
            if (inshell_encounter[mi]<shell){ // do not apply acceleration twice // inshell_encounter = inshell_encounter_passive
                particles[mi].vx += y*particles[mi].ax;
                particles[mi].vy += y*particles[mi].ay;
                particles[mi].vz += y*particles[mi].az;
            }
        }
    }
}

// Main Kernel Operator: Drift step. 
static void reb_integrator_mercurana_drift_step(struct reb_simulation* const r, double a, unsigned int shell){
    rebd_drift[shell]++;
#ifndef OPENMP
    if (reb_sigint) return;
#endif
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    struct reb_particle* restrict const particles = r->particles;
    reb_mercurana_encounter_predict(r, a, shell);
    unsigned int* map_encounter = rim->map_encounter[shell];
    unsigned int* map_encounter_passive = rim->map_encounter_passive[shell];
    unsigned int* map_dominant = rim->map_dominant[shell];
    unsigned int* map_subdominant = rim->map_subdominant[shell];
    unsigned int* map_subdominant_passive = rim->map_subdominant_passive[shell];
    unsigned int* inshell_encounter = rim->inshell_encounter;
    unsigned int* inshell_dominant = rim->inshell_dominant;
    unsigned int* inshell_subdominant = rim->inshell_subdominant;
    
    if (shell+1<rim->Nmaxshells){ // does sub-shell exist? If so, do that first.
        // Are there particles in it?
        if (rim->shellN_encounter[shell+1]>0 || rim->shellN_encounter_passive[shell+1]>0 || rim->shellN_dominant[shell+1]>0){ // no need to check for subdominant, because dominant will also be filled
            rim->Nmaxshellsused = MAX(rim->Nmaxshellsused, shell+2);
            // advance all sub-shell particles
            unsigned int n = rim->n1?rim->n1:rim->n0;
            if (rim->n0>0 && shell==0){
                n = rim->n0; // use different number of substeps for first shell
            }
            rim->t_drifted[shell+1] = rim->t_drifted[shell]; // needed in case subshell was not active before
            const double as = a/n;
            reb_integrator_eos_preprocessor(r, as, shell+1, rim->phi1, reb_integrator_mercurana_drift_step, reb_integrator_mercurana_interaction_step);
            for (int i=0;i<n;i++){
                reb_integrator_eos_step(r, as, 1., 1., shell+1, rim->phi1, reb_integrator_mercurana_drift_step, reb_integrator_mercurana_interaction_step);
            }
            reb_integrator_eos_postprocessor(r, as, shell+1, rim->phi1, reb_integrator_mercurana_drift_step, reb_integrator_mercurana_interaction_step);
        }else{
            r->t += a;
        }
    }else{
        r->t += a;
    }
    
    for (int i=0;i<rim->shellN_dominant[shell];i++){  // loop over all particles in shell (includes subshells)
        int mi = map_dominant[i]; 
        if( inshell_dominant[mi]==shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
            p_t[mi] += a;
        }
    }
    for (int i=0;i<rim->shellN_subdominant[shell];i++){  // loop over all particles in shell (includes subshells)
        int mi = map_subdominant[i]; 
        if( inshell_subdominant[mi]==shell && inshell_encounter[mi]<=shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
            p_t[mi] += a;
        }
    }
    for (int i=0;i<rim->shellN_encounter[shell];i++){  // loop over all particles in shell (includes subshells)
        int mi = map_encounter[i]; 
        if( inshell_subdominant[mi]<shell && inshell_encounter[mi]==shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
            p_t[mi] += a;
        }
    }
    for (int i=0;i<rim->shellN_subdominant_passive[shell];i++){  // loop over all particles in shell (includes subshells)
        int mi = map_subdominant_passive[i]; 
        if( inshell_subdominant[mi]==shell && inshell_encounter[mi]<=shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
            p_t[mi] += a;
        }
    }
    for (int i=0;i<rim->shellN_encounter_passive[shell];i++){  // loop over all particles in shell (includes subshells)
        int mi = map_encounter_passive[i]; 
        if( inshell_subdominant[mi]<shell && inshell_encounter[mi]==shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
            p_t[mi] += a;
        }
    }
    rim->t_drifted[shell] += a;
}

// Part 1 only contains logic for setting up all the data structures. 
// The actual integration is done in part 2.
void reb_integrator_mercurana_part1(struct reb_simulation* r){
    if (r->var_config_N){
        reb_warning(r,"Mercurana does not work with variational equations.");
    }
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    if (rim->Nmaxshells<=0){
        reb_error(r,"Nmaxshells needs to be larger than 0.");
        return;
    }
    if (rim->Nmaxshells==1 && rim->n0 ){
        reb_error(r,"Nmaxshells>=2 is required if n0 is greater than 0.");
        return;
    }
    if (rim->Nmaxshells==2 && rim->n1){
        reb_error(r,"Nmaxshells>=3 is required if n1 is greater than 0.");
        return;
    }
    if (rim->Nmaxshells>1 && rim->kappa<=0.){
        reb_error(r,"kappa>0 is required if Nmaxshells>1.");
        return;
    }
    if (rim->Nmaxshells>1 && rim->n0==0){
        reb_error(r,"n0>0 is required if Nmaxshells>1.");
        return;
    }
    if (rim->N_dominant>r->N_active && r->N_active!=-1){
        reb_error(r,"The number of dominant particles N_dominant cannot be larger than the number of active particles N_active.");
        return;
    }
    
    const int N = r->N;
    
    if (rim->allocatedN<N){
        // dcrit
        if (rim->dcrit){
            for (int i=0;i<rim->Nmaxshells;i++){
                free(rim->dcrit[i]);
            }
        }
        rim->dt_shell = realloc(rim->dt_drift, sizeof(double)*(rim->Nmaxshells));
        rim->dcrit = realloc(rim->dcrit, sizeof(double*)*(rim->Nmaxshells));
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->dcrit[i] = malloc(sizeof(double)*N);
        }
        // map
        if (rim->map_encounter){
            for (int i=0;i<rim->Nmaxshells;i++){
                free(rim->map_encounter[i]);
                free(rim->map_encounter_passive[i]);
                free(rim->map_dominant[i]);
                free(rim->map_subdominant[i]);
                free(rim->map_subdominant_passive[i]);
                free(rim->maxdrift_encounter[i]);
                free(rim->maxdrift_dominant[i]);
                free(rim->p0[i]);
            }
        }
        rim->map_encounter = realloc(rim->map_encounter, sizeof(unsigned int*)*rim->Nmaxshells);
        rim->map_encounter_passive = realloc(rim->map_encounter_passive, sizeof(unsigned int*)*rim->Nmaxshells);
        rim->map_dominant = realloc(rim->map_dominant, sizeof(unsigned int*)*rim->Nmaxshells);
        rim->map_subdominant = realloc(rim->map_subdominant, sizeof(unsigned int*)*rim->Nmaxshells);
        rim->map_subdominant_passive = realloc(rim->map_subdominant_passive, sizeof(unsigned int*)*rim->Nmaxshells);
        rim->maxdrift_dominant = realloc(rim->maxdrift_dominant, sizeof(double*)*rim->Nmaxshells);
        rim->maxdrift_subdominant = realloc(rim->maxdrift_subdominant, sizeof(double*)*rim->Nmaxshells);
        rim->maxdrift_encounter = realloc(rim->maxdrift_encounter, sizeof(double*)*rim->Nmaxshells);
        rim->maxdrift_subdominant_passive = realloc(rim->maxdrift_subdominant_passive, sizeof(double*)*rim->Nmaxshells);
        rim->maxdrift_encounter_passive = realloc(rim->maxdrift_encounter_passive, sizeof(double*)*rim->Nmaxshells);
        rim->p0 = realloc(rim->p0, sizeof(struct reb_particle*)*rim->Nmaxshells);
        rim->p0_drifted = realloc(rim->t_drifted, sizeof(double*)*rim->Nmaxshells);
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->map_encounter[i] = malloc(sizeof(unsigned int)*N);
            rim->map_encounter_passive[i] = malloc(sizeof(unsigned int)*N);
            rim->map_dominant[i] = malloc(sizeof(unsigned int)*N);
            rim->map_subdominant[i] = malloc(sizeof(unsigned int)*N);
            rim->map_subdominant_passive[i] = malloc(sizeof(unsigned int)*N);
            rim->maxdrift_dominant[i] = malloc(sizeof(double)*N);
            rim->maxdrift_subdominant[i] = malloc(sizeof(double)*N);
            rim->maxdrift_encounter[i] = malloc(sizeof(double)*N);
            rim->maxdrift_subdominant_passive[i] = malloc(sizeof(double)*N);
            rim->maxdrift_encounter_passive[i] = malloc(sizeof(double)*N);
            rim->p0[i] = malloc(sizeof(struct reb_particle)*N);
            rim->p0_drifted[i] = malloc(sizeof(double)*N);
        }
        rim->p_drifted = realloc(rim->t_drifted, sizeof(double)*N);
        // inshell
        rim->inshell_encounter = realloc(rim->inshell_encounter, sizeof(unsigned int)*N);
        rim->inshell_dominant = realloc(rim->inshell_dominant, sizeof(unsigned int)*N);
        rim->inshell_subdominant = realloc(rim->inshell_subdominant, sizeof(unsigned int)*N);
        // shellN
        rim->shellN_encounter = realloc(rim->shellN_encounter, sizeof(unsigned int)*rim->Nmaxshells);
        rim->shellN_encounter_passive = realloc(rim->shellN_encounter_passive, sizeof(unsigned int)*rim->Nmaxshells);
        rim->shellN_dominant = realloc(rim->shellN_dominant, sizeof(unsigned int)*rim->Nmaxshells);
        rim->shellN_subdominant = realloc(rim->shellN_subdominant, sizeof(unsigned int)*rim->Nmaxshells);
        rim->shellN_subdominant_passive = realloc(rim->shellN_subdominant_passive, sizeof(unsigned int)*rim->Nmaxshells);
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->shellN_encounter[i] = 0;
            rim->shellN_encounter_passive[i] = 0;
            rim->shellN_dominant[i] = 0;
            rim->shellN_subdominant[i] = 0;
            rim->shellN_subdominant_passive[i] = 0;
        }
        
        // shellN_active
        //rim->shellN_active = realloc(rim->shellN_active, sizeof(unsigned int)*rim->Nmaxshells);

        rim->allocatedN = N;
        // If particle number increased (or this is the first step), need to calculate critical radii
        rim->recalculate_dcrit_this_timestep = 1;
    }

    if (rim->recalculate_dcrit_this_timestep){
        rim->recalculate_dcrit_this_timestep = 0;
        if (rim->is_synchronized==0){
            reb_integrator_mercurana_synchronize(r);
            reb_warning(r,"MERCURANA: Recalculating dcrit but pos/vel were not synchronized before.");
        }
        double dt0 = r->dt;
        double dt_shell = r->dt;
        for (int s=0;s<rim->Nmaxshells;s++){ // innermost shell has no dcrit
            for (int i=0;i<N;i++){
                double mi = r->particles[i].m;
                double dgrav = sqrt3(r->G*dt0*dt0*mi/rim->kappa);
                if (rim->Gm0r0){
                    double dgravrel = sqrt(sqrt(r->G*r->G*dt0*dt0*mi*mi/rim->Gm0r0/rim->kappa));
                    dgrav = MAX(dgrav,dgravrel);
                }
                if (rim->alpha!=0.5){
                    // might not machine independent!
                    rim->dcrit[s][i] = pow(dt_shell/dt0,rim->alpha) * dgrav;
                }else{
                    rim->dcrit[s][i] = sqrt(dt_shell/dt0) * dgrav;
                }
            }
            // Definition: ??
            // - n=0 is not allowed
            // - n=1 means the D in a DKD scheme is calculated using one sub-step with DKD (0.5*dt)
            // - n=2 means the D in a DKD scheme is calculated using two DKD sub steps (0.25*dt each)
            // - n=0 is not allowed
            // - n=1 means the D in a DKDKD scheme is calculates using one sub-step of DKDKD (0.33*dt)
            // - n=2 means the D in a DKDKD scheme is calculated using two DKDKD sub-step (0.1666*dt each)
            double longest_drift_step_in_shell = 0.5; 
            enum REB_EOS_TYPE phi  = s==0?rim->phi0:rim->phi1;
            switch(phi){
                case REB_EOS_LF:
                case REB_EOS_PMLF4:
                    longest_drift_step_in_shell = 0.5; 
                    break;
                case REB_EOS_LF4:
                    longest_drift_step_in_shell = reb_eos_lf4_a;
                    break;
                case REB_EOS_LF6:
                    longest_drift_step_in_shell = reb_eos_lf6_a[0]+reb_eos_lf6_a[1];
                    break; 
                case REB_EOS_LF8: 
                    longest_drift_step_in_shell = reb_eos_lf8_a[0]+reb_eos_lf8_a[1];
                    break;
                case REB_EOS_LF4_2: 
                    longest_drift_step_in_shell = 1.-2.*reb_eos_lf4_2_a;
                    break;
                case REB_EOS_LF8_6_4:
                    longest_drift_step_in_shell = reb_eos_lf8_6_4_a[2];   
                case REB_EOS_PMLF6:
                    longest_drift_step_in_shell = reb_eos_pmlf6_a[1]; 
                    break;
                case REB_EOS_PLF7_6_4:
                    longest_drift_step_in_shell = reb_eos_plf7_6_4_a[0];   
                    break;
            }
            dt_shell *= longest_drift_step_in_shell;
            rim->dt_shell[s] = dt_shell;
            unsigned int n = rim->n0;
            if (s>0 && rim->n1){
                n = rim->n1; // use different number of substeps for deep shells
                             // if n1 is not set, keep using n0
            }
            dt_shell /= n;
            // Initialize shell numbers to zero (not needed, but helps debugging)
            //rim->shellN[s] = 0;
            //rim->shellN_active[s] = 0;
        }

    }
    
    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_warning(r,"Mercurana only works with a direct collision search.");
    }
    
    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_MERCURANA){
        reb_warning(r,"Mercurana has it's own gravity routine. Gravity routine set by the user will be ignored.");
    }
    r->gravity = REB_GRAVITY_NONE; // Only temporary
    
    if (rim->L == NULL){
        // Setting default switching function
        rim->L = reb_integrator_mercurana_L_infinity;
        rim->dLdr = reb_integrator_mercurana_dLdr_infinity;
    }
}

// This routine performs one global timestep
void reb_integrator_mercurana_part2(struct reb_simulation* const r){
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    if (rim->allocatedN<r->N){ // Error occured earlier.
        return;
    }
   
    rim->t_now = 0; 
    for (int s=0;s<rim->Nmaxshells;s++){ // innermost shell has no dcrit
        for (int i=0;i<r->N;i++){
            rim->p0_t_dominant[s][i] =0;
            rim->p0_t_subdominant[s][i] =0;
            rim->p0_t_encounter[s][i] =0;
            rim->p0_t_subdominant_passive[s][i] =0;
            rim->p0_t_encounter_passive[s][i] =0;
        }
    }

    if (rim->is_synchronized){
        reb_integrator_eos_preprocessor(r, r->dt, 0, rim->phi0, reb_integrator_mercurana_drift_step, reb_integrator_mercurana_interaction_step);
    }
    reb_integrator_eos_step(r, r->dt, 1., 1., 0, rim->phi0, reb_integrator_mercurana_drift_step, reb_integrator_mercurana_interaction_step);

    rim->is_synchronized = 0;
    if (rim->safe_mode){
        reb_integrator_mercurana_synchronize(r);
    }

    //r->t+=r->dt;
    r->dt_last_done = r->dt;
}

// Apply post-processor to outermost splitting
void reb_integrator_mercurana_synchronize(struct reb_simulation* r){
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    if (rim->is_synchronized == 0){
        if (rim->L == NULL){
            // Setting default switching function
            rim->L = reb_integrator_mercurana_L_infinity;
            rim->dLdr = reb_integrator_mercurana_dLdr_infinity;
        }
        reb_integrator_eos_postprocessor(r, r->dt, 0, rim->phi0, reb_integrator_mercurana_drift_step, reb_integrator_mercurana_interaction_step);
        rim->is_synchronized = 1;
    }
}

void reb_integrator_mercurana_reset(struct reb_simulation* r){
    if (r->ri_mercurana.allocatedN){
        for (int i=0;i<r->ri_mercurana.Nmaxshells;i++){
            free(r->ri_mercurana.map_encounter[i]);
            free(r->ri_mercurana.map_encounter_passive[i]);
            free(r->ri_mercurana.map_dominant[i]);
            free(r->ri_mercurana.map_subdominant[i]);
            free(r->ri_mercurana.map_subdominant_passive[i]);
            free(r->ri_mercurana.dcrit[i]);
            free(r->ri_mercurana.maxdrift_dominant[i]);
            free(r->ri_mercurana.maxdrift_subdominant[i]);
            free(r->ri_mercurana.maxdrift_encounter[i]);
            free(r->ri_mercurana.maxdrift_subdominant_passive[i]);
            free(r->ri_mercurana.maxdrift_encounter_passive[i]);
            free(r->ri_mercurana.p0[i]);
            free(r->ri_mercurana.p0_drifted[i]);
        }
        free(r->ri_mercurana.map_encounter);
        free(r->ri_mercurana.map_encounter_passive);
        free(r->ri_mercurana.map_dominant);
        free(r->ri_mercurana.map_subdominant);
        free(r->ri_mercurana.map_subdominant_passive);
        free(r->ri_mercurana.dcrit);
        free(r->ri_mercurana.inshell_encounter);
        free(r->ri_mercurana.inshell_dominant);
        free(r->ri_mercurana.inshell_subdominant);
        free(r->ri_mercurana.shellN_encounter);
        free(r->ri_mercurana.shellN_encounter_passive);
        free(r->ri_mercurana.shellN_dominant);
        free(r->ri_mercurana.shellN_subdominant);
        free(r->ri_mercurana.shellN_subdominant_passive);
        free(r->ri_mercurana.p_drifted);
        free(r->ri_mercurana.p0_drifted);
        free(r->ri_mercurana.maxdrift_dominant);
        free(r->ri_mercurana.maxdrift_subdominant);
        free(r->ri_mercurana.maxdrift_encounter);
        free(r->ri_mercurana.maxdrift_subdominant_passive);
        free(r->ri_mercurana.maxdrift_encounter_passive);
        free(r->ri_mercurana.p0);
        free(r->ri_mercurana.dt_shell);
    }
    r->ri_mercurana.allocatedN = 0;
    r->ri_mercurana.map_encounter = NULL;
    r->ri_mercurana.map_encounter_passive = NULL;
    r->ri_mercurana.map_dominant = NULL;
    r->ri_mercurana.map_subdominant = NULL;
    r->ri_mercurana.map_subdominant_passive = NULL;
    r->ri_mercurana.dcrit = NULL;
    r->ri_mercurana.inshell_encounter = NULL;
    r->ri_mercurana.inshell_dominant = NULL;
    r->ri_mercurana.inshell_subdominant = NULL;
    r->ri_mercurana.shellN_encounter = NULL;
    r->ri_mercurana.shellN_encounter_passive = NULL;
    r->ri_mercurana.shellN_dominant = NULL;
    r->ri_mercurana.shellN_subdominant = NULL;
    r->ri_mercurana.shellN_subdominant_passive = NULL;
    r->ri_mercurana.maxdrift_dominant = NULL;
    r->ri_mercurana.maxdrift_subdominant = NULL;
    r->ri_mercurana.maxdrift_encounter = NULL;
    r->ri_mercurana.maxdrift_subdominant_passive = NULL;
    r->ri_mercurana.maxdrift_encounter_passive = NULL;
    r->ri_mercurana.p0 = NULL;
    
    r->ri_mercurana.phi0 = REB_EOS_LF;
    r->ri_mercurana.phi1 = REB_EOS_LF;
    r->ri_mercurana.n0 = 2;
    r->ri_mercurana.n1 = 0;
    r->ri_mercurana.kappa = 1e-3;
    r->ri_mercurana.Gm0r0 = 0.;
    r->ri_mercurana.alpha = 0.5;
    r->ri_mercurana.safe_mode = 1;
    r->ri_mercurana.Nmaxshells = 10;
    r->ri_mercurana.Nmaxshellsused = 1;
    r->ri_mercurana.recalculate_dcrit_this_timestep = 0;
    r->ri_mercurana.is_synchronized = 1;
    r->ri_mercurana.L = NULL;
    r->ri_mercurana.dLdr = NULL;
    r->ri_mercurana.collisions_N = 0;
    r->ri_mercurana.N_dominant = 0;
}

