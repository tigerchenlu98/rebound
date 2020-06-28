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
static inline void reb_mercurana_record_collision(struct reb_simulation* const r, unsigned int i, unsigned int j){
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

static inline double reb_drift_from_straight_line(struct reb_particle p0, double dt0, struct reb_particle p1, double dt1){
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
static inline unsigned int check_one(struct reb_simulation* r, double dt, unsigned int shell, int mi, 
        struct reb_mdd* mdd, struct reb_pisd pisd_interacting){
    struct reb_particle* const particles = r->particles;
    const double* const dcrit = r->ri_mercurana.dcrit[shell];
    double* const p_t = r->ri_mercurana.p_t;
    struct reb_particle pi = particles[mi];
    if (mdd){  // only calculate maxdrift if check_maxdrift==1
        mdd->p0 = particles[mi];
        mdd->t0 = r->ri_mercurana.t_now;
        mdd->maxdrift = 1e300;
    }
    unsigned int moved = 0;
    for (int j=0; j<pisd_interacting.shellN[shell]; j++){
        int mj = pisd_interacting.map[shell][j]; 
        if (mi==mj) continue; // no self-interactions
        // Drift particle mj to current time for rmin calculation
        struct reb_particle pj = particles[mj];
        double pjdrift = r->ri_mercurana.t_now - p_t[mj];
        if (pjdrift!=0){ printf("pjdrift: %e\n",pjdrift);};
        pj.x += pjdrift * pj.vx;
        pj.y += pjdrift * pj.vy;
        pj.z += pjdrift * pj.vz;
        double rmin2 = reb_mercurana_predict_rmin2(pi,pj,dt);
        double rsum = pi.r+pj.r;
        if (rmin2< rsum*rsum && r->collision==REB_COLLISION_DIRECT){
            // TODO: Think about right location for collision check
            reb_mercurana_record_collision(r,mi,mj);
        }
        double dcritsum = dcrit[mi]+dcrit[mj];
        if (rmin2< dcritsum*dcritsum){ 
            if (pisd_interacting.inshell[mj] == shell){
                moved++;
                // Drift particle to now
                particles[mj] = pj;
                p_t[mj] = r->ri_mercurana.t_now;
                // Move particle in subshell
                pisd_interacting.inshell[mj] = shell+1;
                pisd_interacting.map[shell+1][pisd_interacting.shellN[shell+1]] = mj;
                pisd_interacting.shellN[shell+1]++;
            }
        }else{
            if (mdd){  // only calculate maxdrift if check_maxdrift==1
                // Interaction potentially not resolved in subshells. Need maxdrift.
                const double maxdrift = (sqrt(rmin2) - dcritsum)/2.;
                mdd->maxdrift = MIN(mdd->maxdrift,maxdrift);
            }
        }
    }
    return moved;
}

static void reb_mercurana_encounter_predict(struct reb_simulation* const r, double dt, int shell){
    struct reb_simulation_integrator_mercurana* rim = &(r->ri_mercurana);
    
    rim->collisions_N = 0;
    
    if (shell+1>=rim->Nmaxshells){ // does sub-shell exist?
        return; // If not, no prediction needed.
    }

    for (int ptype=0; ptype<5; ptype++){ 
        rim->pisd[ptype].shellN[shell+1] = 0;
        for (int i=0; i<rim->pisd[ptype].shellN[shell]; i++){
            int mi = rim->pisd[ptype].map[shell][i]; 
            rim->pisd[ptype].inshell[mi] = shell;
        }
    }

    if (shell!=0 && rim->check_maxdrift){
        struct reb_particle* const particles = r->particles;
        for (int itype=0; itype<8; itype++){ 
            enum REB_PTYPE drifting = get_ptype_drifting(itype);
            enum REB_PTYPE interacting = get_ptype_interacting(itype);

            struct reb_pisd pisd_drifting = r->ri_mercurana.pisd[drifting];
            struct reb_pisd pisd_interacting = r->ri_mercurana.pisd[interacting];
            
            for (int i=0; i<pisd_drifting.shellN[shell]; i++){
                int mi = pisd_drifting.map[shell][i]; 
                for (int s=0;s<shell;s++){
                    struct reb_mdd* mdd = &(r->ri_mercurana.mdd[itype][s][mi]);
                    double dt0 = r->ri_mercurana.t_now - mdd->t0 + dt;
                    double drift = reb_drift_from_straight_line(mdd->p0,dt0,particles[mi],dt);
                    if (drift > mdd->maxdrift){
                        unsigned int Nmoved = check_one(r, dt, s, mi, mdd, pisd_interacting);
                        r->ri_mercurana.Nmoved += Nmoved;
                    }
                }
            }
        } 
    }

    for (int itype=0; itype<8; itype++){ 
        enum REB_PTYPE drifting = get_ptype_drifting(itype);
        enum REB_PTYPE interacting = get_ptype_interacting(itype);

        struct reb_pisd pisd_drifting = r->ri_mercurana.pisd[drifting];
        struct reb_pisd pisd_interacting = r->ri_mercurana.pisd[interacting];

        for (int i=0; i<pisd_drifting.shellN[shell]; i++){
            int mi = pisd_drifting.map[shell][i]; 
            struct reb_mdd* mdd = NULL;
            if (rim->check_maxdrift){
                mdd = &(r->ri_mercurana.mdd[itype][shell][mi]);
            }
            check_one(r, dt, shell, mi, mdd, pisd_interacting);
        }
    }
    
    if (rim->collisions_N){
        unsigned int N_before = r->N;
        reb_collision_search(r); 
        if (N_before!=r->N){
            reb_mercurana_encounter_predict(r, dt, shell);
        }
        rim->collisions_N = 0;
    }
}


static inline void reb_mercurana_grav_update_B(struct reb_simulation* const r, 
        const double* dcrit_i, const double* dcrit_c, const double* dcrit_o, const double G, const double softening2,
        int shell, double y, enum REB_PTYPE ptype_A, enum REB_PTYPE ptype_B){
    struct reb_pisd pisd_A = r->ri_mercurana.pisd[ptype_A];
    struct reb_pisd pisd_B = r->ri_mercurana.pisd[ptype_B];
    
    struct reb_particle* const particles = r->particles;
    double (*_L) (const struct reb_simulation* const r, double d, double dcrit, double fracin) = r->ri_mercurana.L;
#pragma omp parallel for
    for (int i=0; i<pisd_A.shellN[shell]; i++){
        const int mi = pisd_A.map[shell][i];
        for (int j=0; j<pisd_B.shellN[shell]; j++){
            const int mj = pisd_B.map[shell][j];
            if (mi==mj) continue; // avoid self interaction. Only needed for OPENMP
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dc_c = dcrit_c[mi]+dcrit_c[mj];
            double Lsum = 0.;
            if (dcrit_o){
                double dc_o = dcrit_o[mi]+dcrit_o[mj];
                Lsum -= _L(r,dr,dc_c,dc_o);
            }
            if (dcrit_i){
                double dc_i = dcrit_i[mi]+dcrit_i[mj];
                Lsum += _L(r,dr,dc_i,dc_c);
            }else{
                Lsum += 1; 
            }
            const double prefact = y*G*Lsum/(dr*dr*dr);
            const double prefacti = prefact*particles[mi].m;
            particles[mj].vx    += prefacti*dx;
            particles[mj].vy    += prefacti*dy;
            particles[mj].vz    += prefacti*dz;
        }
    }
}
static inline void reb_mercurana_grav_update_AB(struct reb_simulation* const r, 
        const double* dcrit_i, const double* dcrit_c, const double* dcrit_o, const double G, const double softening2,
        int shell, double y, enum REB_PTYPE ptype_A, enum REB_PTYPE ptype_B){
    
    struct reb_pisd pisd_A = r->ri_mercurana.pisd[ptype_A];
    struct reb_pisd pisd_B = r->ri_mercurana.pisd[ptype_B];

    struct reb_particle* const particles = r->particles;
    double (*_L) (const struct reb_simulation* const r, double d, double dcrit, double fracin) = r->ri_mercurana.L;
#ifndef OPENMP // OPENMP off, do O(1/2*N^2)
#pragma omp parallel for
    for (int i=0; i<pisd_A.shellN[shell]; i++){
        const int startj = (ptype_A==ptype_B)?(i+1):0;
        const int mi = pisd_A.map[shell][i];
        for (int j=startj; j<pisd_B.shellN[shell]; j++){
            const int mj = pisd_B.map[shell][j];
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double dr = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dc_c = dcrit_c[mi]+dcrit_c[mj];
            double Lsum = 0.;
            if (dcrit_o){
                double dc_o = dcrit_o[mi]+dcrit_o[mj];
                Lsum -= _L(r,dr,dc_c,dc_o);
            }
            if (dcrit_i){
                double dc_i = dcrit_i[mi]+dcrit_i[mj];
                Lsum += _L(r,dr,dc_i,dc_c);
            }else{
                Lsum += 1; 
            }
            const double prefact = y*G*Lsum/(dr*dr*dr);
            const double prefactj = -prefact*particles[mj].m;
            const double prefacti = prefact*particles[mi].m;
            particles[mi].vx    += prefactj*dx;
            particles[mi].vy    += prefactj*dy;
            particles[mi].vz    += prefactj*dz;
            particles[mj].vx    += prefacti*dx;
            particles[mj].vy    += prefacti*dy;
            particles[mj].vz    += prefacti*dz;
        }
    }
#else // OPENMP on, do O(N^2)
    reb_mercurana_grav_update_B(r, dcrit_i, dcrit_c, dcrit_o, G, softening2, shell, y, ptype_A, ptpye_B);
    reb_mercurana_grav_update_B(r, dcrit_i, dcrit_c, dcrit_o, G, softening2, shell, y, ptype_B, ptpye_A);
#endif
}

// Main Kernel Operator: Interaction step. 
// y = timestep for acceleration
// v = timestep for jerk (0 if not used)
static void reb_integrator_mercurana_interaction_step(struct reb_simulation* r, double y, double v, unsigned int shell){
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);

    const double* dcrit_i = NULL; // critical radius of inner shell
    const double* dcrit_c = NULL; // critical radius of current shell
    const double* dcrit_o = NULL; // critical radius of outer shell
    if (shell<rim->Nmaxshells-1){
        dcrit_i = r->ri_mercurana.dcrit[shell+1];
    }
    dcrit_c = r->ri_mercurana.dcrit[shell];
    if (shell>0){
        dcrit_o = r->ri_mercurana.dcrit[shell-1];
    }

    reb_mercurana_grav_update_AB(r, dcrit_i, dcrit_c, dcrit_o, G, softening2, shell, y, REB_PTYPE_ENC, REB_PTYPE_ENC);
    reb_mercurana_grav_update_AB(r, dcrit_i, dcrit_c, dcrit_o, G, softening2, shell, y, REB_PTYPE_DOM, REB_PTYPE_SUB);
    reb_mercurana_grav_update_AB(r, dcrit_i, dcrit_c, dcrit_o, G, softening2, shell, y, REB_PTYPE_DOM, REB_PTYPE_DOM);
    if (r->testparticle_type==1){
        reb_mercurana_grav_update_AB(r, dcrit_i, dcrit_c, dcrit_o, G, softening2, shell, y, REB_PTYPE_DOM, REB_PTYPE_SUBP);
        reb_mercurana_grav_update_AB(r, dcrit_i, dcrit_c, dcrit_o, G, softening2, shell, y, REB_PTYPE_ENC, REB_PTYPE_ENCP);
    }else{
        reb_mercurana_grav_update_B(r, dcrit_i, dcrit_c, dcrit_o, G, softening2, shell, y, REB_PTYPE_DOM, REB_PTYPE_SUBP);
        reb_mercurana_grav_update_B(r, dcrit_i, dcrit_c, dcrit_o, G, softening2, shell, y, REB_PTYPE_ENC, REB_PTYPE_ENCP);
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
    
    if (shell+1<rim->Nmaxshells){ // does sub-shell exist? If so, do that first.
        // Are there particles in it?
        int particle_in_subshell = 0;
        for (int ptype=0; ptype<5; ptype++){ 
            particle_in_subshell |= rim->pisd[ptype].shellN[shell+1];
        }
        if (particle_in_subshell){ 
            rim->Nmaxshellsused = MAX(rim->Nmaxshellsused, shell+2);
            // advance all sub-shell particles
            unsigned int n = rim->n1?rim->n1:rim->n0;
            if (rim->n0>0 && shell==0){
                n = rim->n0; // use different number of substeps for first shell
            }
            const double as = a/n;
            reb_integrator_eos_preprocessor(r, as, shell+1, rim->phi1, reb_integrator_mercurana_drift_step, reb_integrator_mercurana_interaction_step);
            for (int i=0;i<n;i++){
                reb_integrator_eos_step(r, as, 1., 1., shell+1, rim->phi1, reb_integrator_mercurana_drift_step, reb_integrator_mercurana_interaction_step);
            }
            reb_integrator_eos_postprocessor(r, as, shell+1, rim->phi1, reb_integrator_mercurana_drift_step, reb_integrator_mercurana_interaction_step);
        }else{
            r->t += a;
            rim->t_now += a; // measured from beginning of timestep to avoid roundoff issues
        }
    }else{
        r->t += a;
        rim->t_now += a;
    }
    
    for (int i=0;i<rim->pisd[REB_PTYPE_DOM].shellN[shell];i++){  
        int mi = rim->pisd[REB_PTYPE_DOM].map[shell][i]; 
        if( rim->pisd[REB_PTYPE_DOM].inshell[mi]==shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
            rim->p_t[mi] = rim->t_now;
        }
    }
    for (int i=0;i<rim->pisd[REB_PTYPE_SUB].shellN[shell];i++){  
        int mi = rim->pisd[REB_PTYPE_SUB].map[shell][i]; 
        if( rim->pisd[REB_PTYPE_SUB].inshell[mi]==shell && rim->pisd[REB_PTYPE_ENC].inshell[mi]<=shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
            rim->p_t[mi] = rim->t_now;
        }
    }
    for (int i=0;i<rim->pisd[REB_PTYPE_ENC].shellN[shell];i++){  
        int mi = rim->pisd[REB_PTYPE_ENC].map[shell][i]; 
        if( rim->pisd[REB_PTYPE_ENC].inshell[mi]==shell && rim->pisd[REB_PTYPE_SUB].inshell[mi]<shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
            rim->p_t[mi] = rim->t_now;
        }
    }
    for (int i=0;i<rim->pisd[REB_PTYPE_SUBP].shellN[shell];i++){  
        int mi = rim->pisd[REB_PTYPE_SUBP].map[shell][i]; 
        if( rim->pisd[REB_PTYPE_SUBP].inshell[mi]==shell && rim->pisd[REB_PTYPE_ENCP].inshell[mi]<=shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
            rim->p_t[mi] = rim->t_now;
        }
    }
    for (int i=0;i<rim->pisd[REB_PTYPE_ENCP].shellN[shell];i++){  
        int mi = rim->pisd[REB_PTYPE_ENCP].map[shell][i]; 
        if( rim->pisd[REB_PTYPE_ENCP].inshell[mi]==shell && rim->pisd[REB_PTYPE_SUBP].inshell[mi]<shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
            rim->p_t[mi] = rim->t_now;
        }
    }
    //if (r->heartbeat){ r->heartbeat(r); }               // Heartbeat
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
    if (rim->Nmaxshells>1 && rim->n0==0){
        reb_error(r,"n0>0 is required if Nmaxshells>1.");
        return;
    }
    if (rim->N_dominant>r->N_active && r->N_active!=-1){
        reb_error(r,"The number of dominant particles N_dominant cannot be larger than the number of active particles N_active.");
        return;
    }
    
    const int N = r->N;
    const int N_active = r->N_active==-1?N:r->N_active;
    const int N_dominant = rim->N_dominant;
    
    if (rim->allocatedN<N || (rim->mdd[0]==NULL && rim->check_maxdrift)){ // if check_maxdrift enabled during integration
        if (rim->dcrit){
            for (int i=0;i<rim->Nmaxshells;i++){
                free(rim->dcrit[i]);
            }
        }
        rim->dcrit = realloc(rim->dcrit, sizeof(double*)*(rim->Nmaxshells));
        rim->p_t = realloc(rim->p_t, sizeof(double*)*N);
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->dcrit[i] = malloc(sizeof(double)*N);
        }
        for (int ptype=0; ptype<5; ptype++){ 
            if (rim->pisd[ptype].map){
                for (int i=0;i<rim->Nmaxshells;i++){
                    free(rim->pisd[ptype].map[i]);
                }
            }
            rim->pisd[ptype].map = realloc(rim->pisd[ptype].map, sizeof(unsigned int*)*rim->Nmaxshells);
            for (int i=0;i<rim->Nmaxshells;i++){
                rim->pisd[ptype].map[i] = malloc(sizeof(unsigned int)*N);
            }
            rim->pisd[ptype].inshell = realloc(rim->pisd[ptype].inshell, sizeof(unsigned int*)*N);
            rim->pisd[ptype].shellN = realloc(rim->pisd[ptype].shellN, sizeof(unsigned int*)*rim->Nmaxshells);
        }
        if (rim->check_maxdrift){
            for (int itype=0; itype<8; itype++){ 
                if (rim->mdd[itype]){
                    for (int i=0;i<rim->Nmaxshells;i++){
                        free(rim->mdd[itype]);
                    }
                }
                rim->mdd[itype] = realloc(rim->mdd[itype], sizeof(unsigned int*)*rim->Nmaxshells);
                for (int i=0;i<rim->Nmaxshells;i++){
                    rim->mdd[itype][i] = malloc(sizeof(struct reb_mdd)*N);
                }
            }
        }


        rim->allocatedN = N;
        // If particle number increased (or this is the first step), need to calculate critical radii
        rim->recalculate_dcrit_this_timestep = 1;
    }

    // Mass ratio calculation
    if (rim->massratio==-1 || rim->recalculate_dcrit_this_timestep){
        if (N_dominant==0 || N==0 || N_dominant==N_active){
            rim->massratio = 1;
        }else{
            double Mmindom = r->particles[0].m;
            double Mmaxsub = 0;
            for (int i=1;i<N_dominant;i++){
                Mmindom = MIN(Mmindom,r->particles[i].m);
            }
            for (int i=N_dominant;i<N;i++){
                Mmaxsub = MAX(Mmaxsub,r->particles[i].m);
            }
            rim->massratio = Mmaxsub/Mmindom;
            if (!isnormal(rim->massratio)){
                reb_warning(r,"MERCURANA: Unable to automatically calculate massratio. Setting it to 1. Consider setting it manually.");
                rim->massratio = 1;
            }
        }
    }
    
    // Rmin calculation
    if (rim->rmin==-1 || rim->recalculate_dcrit_this_timestep){
        if (N>0){
            rim->rmin = r->particles[0].r;
            for (int i=1;i<N_active;i++){
                rim->rmin = MIN(rim->rmin,r->particles[i].r);
            }
            if (!isnormal(rim->rmin)){
                reb_warning(r,"MERCURANA: Unable to automatically calculate rmin. Setting it to 1. Consider setting it manually.");
                rim->rmin = 1;
            }
        }
    }
    
    // Rmax calculation
    if (rim->rref==-1 || rim->recalculate_dcrit_this_timestep){
        if (N>0){
            rim->rref = 0;
            for (int i=1;i<N_active;i++){
                double x = r->particles[i].x;
                double y = r->particles[i].y;
                double z = r->particles[i].z;
                rim->rref = MAX(rim->rref,sqrt(x*x+y*y+z*z));
            }
            if (!isnormal(rim->rref)){
                reb_warning(r,"MERCURANA: Unable to automatically calculate rref. Setting it to 1. Consider setting it manually.");
                rim->rref = 1;
            }
        }
    }
    // Epsilon calculation
    if (rim->epsilon==-1 || rim->recalculate_dcrit_this_timestep){
        if (N_dominant==0 || rim->massratio<=0. ){
            reb_warning(r,"MERCURANA: Unable to automatically calculate epsilon. Setting it to 1e-3. Consider setting it manually.");
            rim->epsilon = 1e-3;
        }else{
            double T0min = 1e300;
            for (int i=0;i<N_dominant;i++){
                for (int j=N_dominant;j<N;j++){
                    double x =  r->particles[i].x - r->particles[j].x;
                    double y =  r->particles[i].y - r->particles[j].y;
                    double z =  r->particles[i].z - r->particles[j].z;
                    double a = sqrt(x*x+y*y+z*z);
                    T0min = MIN(T0min,sqrt(a*a*a/(r->G*r->particles[i].m)));
                }
            }
            if (!isnormal(T0min) || T0min==1e300){
                reb_warning(r,"MERCURANA: Unable to automatically calculate epsilon. Setting it to 1e-3. Consider setting it manually.");
                rim->epsilon = 1e-3;
            }else{
                rim->epsilon = r->dt*r->dt/(T0min*T0min);
            }
        }
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
                // Constrain two body error
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
    r->gravity = REB_GRAVITY_MERCURANA; 
    
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
   
    rim->t_now = 0.;
    for (int i=0;i<r->N;i++){
        rim->p_t[i] = 0.;
    }
    for (int i=0;i<rim->Nmaxshells;i++){
        for (int ptype=0; ptype<5; ptype++){ 
            rim->pisd[ptype].shellN[i] = 0;
        }
    }
    
    // Setup maps in outermost shell 
    const int N_active = r->N_active==-1?r->N:r->N_active;
    
    rim->pisd[REB_PTYPE_DOM].shellN[0] = rim->N_dominant;
    rim->pisd[REB_PTYPE_SUB].shellN[0] = N_active - rim->N_dominant;
    rim->pisd[REB_PTYPE_ENC].shellN[0] = N_active - rim->N_dominant;
    rim->pisd[REB_PTYPE_SUBP].shellN[0] = r->N - N_active;
    rim->pisd[REB_PTYPE_ENCP].shellN[0] = r->N - N_active;
    for (int i=0; i<rim->pisd[REB_PTYPE_DOM].shellN[0]; i++){
        rim->pisd[REB_PTYPE_DOM].map[0][i] = i; 
        rim->pisd[REB_PTYPE_DOM].inshell[i] = 0; 
    }
    for (int i=0; i<rim->pisd[REB_PTYPE_SUB].shellN[0]; i++){
        rim->pisd[REB_PTYPE_SUB].map[0][i] = i + rim->N_dominant; 
        rim->pisd[REB_PTYPE_SUB].inshell[i + rim->N_dominant] = 0; 
    }
    for (int i=0; i<rim->pisd[REB_PTYPE_ENC].shellN[0]; i++){
        rim->pisd[REB_PTYPE_ENC].map[0][i] = i + rim->N_dominant; 
        rim->pisd[REB_PTYPE_ENC].inshell[i + rim->N_dominant] = 0; 
    }
    for (int i=0; i<rim->pisd[REB_PTYPE_SUBP].shellN[0]; i++){
        rim->pisd[REB_PTYPE_SUBP].map[0][i] = i + N_active; 
        rim->pisd[REB_PTYPE_SUBP].inshell[i + N_active] = 0; 
    }
    for (int i=0; i<rim->pisd[REB_PTYPE_ENCP].shellN[0]; i++){
        rim->pisd[REB_PTYPE_ENCP].map[0][i] = i + N_active; 
        rim->pisd[REB_PTYPE_ENCP].inshell[i + N_active] = 0;
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
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    if (r->ri_mercurana.allocatedN){
        for (int ptype=0; ptype<5; ptype++){ 
            if (rim->pisd[ptype].map){
                for (int i=0;i<rim->Nmaxshells;i++){
                    free(rim->pisd[ptype].map[i]);
                }
            }
            free(rim->pisd[ptype].map);
            free(rim->pisd[ptype].inshell);
            free(rim->pisd[ptype].shellN);
        }
        for (int itype=0; itype<8; itype++){ 
            if (rim->mdd[itype]){
                for (int i=0;i<rim->Nmaxshells;i++){
                    free(rim->mdd[itype][i]);
                }
            }
            free(rim->mdd[itype]);
        }
        for (int i=0;i<r->ri_mercurana.Nmaxshells;i++){
            free(r->ri_mercurana.dcrit[i]);
        }
        free(r->ri_mercurana.dcrit);
        free(r->ri_mercurana.p_t);
    }
    r->ri_mercurana.allocatedN = 0;
    r->ri_mercurana.dcrit = NULL;
    r->ri_mercurana.phi0 = REB_EOS_LF;
    r->ri_mercurana.phi1 = REB_EOS_LF;
    r->ri_mercurana.n0 = 2;
    r->ri_mercurana.n1 = 0;
    r->ri_mercurana.epsilon = -1;
    r->ri_mercurana.massratio = -1.;
    r->ri_mercurana.rmin = -1.;
    r->ri_mercurana.rref = -1.;
    r->ri_mercurana.safe_mode = 1;
    r->ri_mercurana.check_maxdrift = 1;
    r->ri_mercurana.Nmaxshells = 10;
    r->ri_mercurana.Nmaxshellsused = 1;
    r->ri_mercurana.Nmoved = 0;
    r->ri_mercurana.recalculate_dcrit_this_timestep = 0;
    r->ri_mercurana.is_synchronized = 1;
    r->ri_mercurana.L = NULL;
    r->ri_mercurana.dLdr = NULL;
    r->ri_mercurana.collisions_N = 0;
    r->ri_mercurana.N_dominant = 0;
}

