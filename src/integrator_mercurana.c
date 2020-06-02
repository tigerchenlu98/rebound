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
unsigned long rebd_kick[MAXSHELLS];
unsigned long rebd_md1a[MAXSHELLS];
unsigned long rebd_md1b[MAXSHELLS];
unsigned long rebd_md1c[MAXSHELLS];
unsigned long rebd_md2a[MAXSHELLS];
unsigned long rebd_md2b[MAXSHELLS];
unsigned long rebd_md2c[MAXSHELLS];
unsigned long rebd_md3a[MAXSHELLS];
unsigned long rebd_md3b[MAXSHELLS];
unsigned long rebd_md3c[MAXSHELLS];
unsigned long rebd_md4a[MAXSHELLS];
unsigned long rebd_md4b[MAXSHELLS];
unsigned long rebd_md4c[MAXSHELLS];



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
    const double t_closest = (dx1*dvx1 + dy1*dvy1 + dz1*dvz1)/(dvx1*dvx1 + dvy1*dvy1 + dvz1*dvz1);
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

static double reb_mercurana_predict_rmin2_drifted(struct reb_particle p1, struct reb_particle p2, double dt, double p2drift){ 
    struct reb_particle p2drifted = p2;
    p2drifted.x += p2drift * p2drifted.vx;
    p2drifted.y += p2drift * p2drifted.vy;
    p2drifted.z += p2drift * p2drifted.vz;
    return reb_mercurana_predict_rmin2(p1, p2drifted,dt);
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

// This function checks if there are any close encounters or physical collisions between particles in a given shell during a drift step of length dt. If a close encounter occures, particles are placed in deeper shells.
static void reb_mercurana_encounter_predict(struct reb_simulation* const r, double dt, int shell){
    struct reb_simulation_integrator_mercurana* rim = &(r->ri_mercurana);
    struct reb_particle* const particles = r->particles;
    const double* const dcrit = rim->dcrit[shell];
    int shellN_encounter = rim->shellN_encounter[shell];
    int shellN_dominant = rim->shellN_dominant[shell];
    int shellN_subdominant = rim->shellN_subdominant[shell];
    
    unsigned int* map_encounter = rim->map_encounter[shell];
    unsigned int* map_dominant = rim->map_dominant[shell];
    unsigned int* map_subdominant = rim->map_subdominant[shell];
    
    double** maxdrift_dominant = rim->maxdrift_dominant;
    double** maxdrift_encounter = rim->maxdrift_encounter;
    
    unsigned int* inshell_encounter = rim->inshell_encounter;
    unsigned int* inshell_dominant = rim->inshell_dominant;
    unsigned int* inshell_subdominant = rim->inshell_dominant;
    struct reb_particle** p0 = rim->p0;
    
    double* t_drifted = rim->t_drifted;
    double* dt_drift = rim->dt_drift;
    
    if (shell+1>=rim->Nmaxshells){ // does sub-shell exist?
        return;
    }
	
    rim->collisions_N = 0;

    rim->shellN_encounter[shell+1] = 0;
    rim->shellN_dominant[shell+1] = 0;
    rim->shellN_subdominant[shell+1] = 0;

    for (int i=0; i<shellN_encounter; i++){
        int mi = map_encounter[i]; 
        inshell_encounter[mi] = shell;
    }
    for (int i=0; i<shellN_dominant; i++){
        int mi = map_dominant[i]; 
        inshell_dominant[mi] = shell;
    }
    for (int i=0; i<shellN_subdominant; i++){
        int mi = map_subdominant[i]; 
        inshell_subdominant[mi] = shell;
    }
    if (shell==0){
        // Setup maps in outermost shell 
        rim->shellN_dominant[0] = rim->N_dominant;
        rim->shellN_subdominant[0] = r->N - rim->N_dominant;
        rim->shellN_encounter[0] = r->N - rim->N_dominant;
        shellN_encounter = rim->shellN_encounter[shell];
        shellN_dominant = rim->shellN_dominant[shell];
        shellN_subdominant = rim->shellN_subdominant[shell];
        for (int i=0; i<r->N; i++){
            p0[0][i] = particles[i]; 
            maxdrift_encounter[0][i] = 1e300; 
            maxdrift_dominant[0][i] = 1e300; 
        }
        for (int i=0; i<shellN_dominant; i++){
            map_dominant[i] = i; 
        }
        for (int i=0; i<shellN_subdominant; i++){
            map_subdominant[i] = shellN_dominant + i; 
            map_encounter[i] = shellN_dominant + i; 
        }
        // Repeat from above
        for (int i=0; i<shellN_encounter; i++){
            int mi = map_encounter[i]; 
            inshell_encounter[mi] = shell;
        }
        for (int i=0; i<shellN_dominant; i++){
            int mi = map_dominant[i]; 
            inshell_dominant[mi] = shell;
        }
        for (int i=0; i<shellN_subdominant; i++){
            int mi = map_subdominant[i]; 
            inshell_subdominant[mi] = shell;
        }
    }else{
        // Check for maxdrift violation of higher shells

        double buffer = 1.05; // prevent having to check the same particle over and over again in one timestep

        // Dominant - Dominant
        for (int i=0; i<shellN_dominant; i++){
            int mi = map_dominant[i]; 
            for (int s=0;s<shell;s++){
                double dt0 = t_drifted[shell] + dt - t_drifted[s];
                double dt1 = dt; 
                double drift = reb_drift_from_straight_line(p0[s][mi],dt0,particles[mi],dt1);
                if (drift>maxdrift_dominant[s][mi]){
                    rebd_md1a[shell]++;
                    //printf("MAXDRIFT VIOLATION DOM-DOM triggered shell. checking particles. %2d %2d  %3d   %e   %e\n",shell,s,mi,drift,maxdrift_dominant[s][mi]);
                    maxdrift_dominant[s][mi] = 1e300;
                    for (int j=0; j<rim->shellN_dominant[s]; j++){ 
                        int mj = rim->map_dominant[s][j]; 
                        if (mi==mj) continue;
                        // Not yet in current shell?
                        unsigned int mj_shell = inshell_dominant[mj];
                        if (mj_shell<shell){
                            rebd_md1b[shell]++;
                            double drift = t_drifted[shell] - t_drifted[mj_shell];
                            double rmin2 = reb_mercurana_predict_rmin2_drifted(particles[mi],particles[mj],dt_drift[s],drift);
                            double dcritsum = buffer*(dcrit[mi]+dcrit[mj]);
//TODO: compare with maxdrift which triggered search to avoid repeated searches
                            if (rmin2< dcritsum*dcritsum){ 
                                rebd_md1c[shell]++;
                                // Add mj to all higher shells as dominant:
                                inshell_dominant[mj] = shell;
                                printf("MAXDRIFT VIOLATION DOM-DOM triggered. moving particle.\n");
                                rim->moved_particles++;
                                for (int sa=mj_shell+1; sa<=shell; sa++){
                                    rim->map_dominant[sa][rim->shellN_dominant[sa]] = mj;
                                    rim->shellN_dominant[sa]++;
                                }
                                particles[mj].x += drift*particles[mj].vx;
                                particles[mj].y += drift*particles[mj].vy;
                                particles[mj].z += drift*particles[mj].vz;
                                p0[shell][mj] = particles[mj]; 
                            }else{
                                double dcritsum = dcrit[mi]+dcrit[mj];
                                const double maxdrift = (sqrt(rmin2) - dcritsum)/2.;
                                maxdrift_dominant[s][mi] = MIN(maxdrift_dominant[s][mi],maxdrift);
                            }
                        }
                    }
                }
            }
        }
        // Dominant - Subdominant
        for (int i=0; i<shellN_dominant; i++){
            int mi = map_dominant[i]; 
            for (int s=0;s<shell;s++){
                double dt0 = t_drifted[shell] + dt - t_drifted[s];
                double dt1 = dt; 
                double drift = reb_drift_from_straight_line(p0[s][mi],dt0,particles[mi],dt1);
                if (drift>maxdrift_encounter[s][mi]){
                    rebd_md2a[shell]++;
                    //printf("MAXDRIFT VIOLATION DOM-SUB triggered shell. checking particles. %2d %2d  %3d   %e   %e\n",shell,s,mi,drift,maxdrift_encounter[s][mi]);
                    maxdrift_encounter[s][mi] = 1e300;
                    for (int j=0; j<rim->shellN_subdominant[s]; j++){ 
                        int mj = rim->map_subdominant[s][j]; 
                        if (mi==mj) continue;
                        // Not yet in current shell?
                        unsigned int mj_shell = inshell_subdominant[mj];
                        if (mj_shell<shell){
                            rebd_md2b[shell]++;
                            double drift = t_drifted[shell] - t_drifted[mj_shell];
                            double rmin2 = reb_mercurana_predict_rmin2_drifted(particles[mi],particles[mj],dt_drift[s],drift);
                            double dcritsum = buffer*(dcrit[mi]+dcrit[mj]);
                            if (rmin2< dcritsum*dcritsum){ 
                                rebd_md2c[shell]++;
                                // Add mj to all higher shells as subdominant:
                                inshell_subdominant[mj] = shell;
                                printf("MAXDRIFT VIOLATION DOM-SUB triggered. moving particle.\n");
                                rim->moved_particles++;
                                for (int sa=mj_shell+1; sa<=shell; sa++){
                                    rim->map_subdominant[sa][rim->shellN_subdominant[sa]] = mj;
                                    rim->shellN_subdominant[sa]++;
                                }
                                particles[mj].x += drift*particles[mj].vx;
                                particles[mj].y += drift*particles[mj].vy;
                                particles[mj].z += drift*particles[mj].vz;
                                p0[shell][mj] = particles[mj]; 
                            }else{
                                double dcritsum = dcrit[mi]+dcrit[mj];
                                const double maxdrift = (sqrt(rmin2) - dcritsum)/2.;
                                maxdrift_encounter[s][mi] = MIN(maxdrift_encounter[s][mi],maxdrift);
                            }
                        }
                    }
                }
            }
        }
        // Subdominant - Dominant
        for (int i=0; i<shellN_subdominant; i++){
            int mi = map_subdominant[i]; 
            for (int s=0;s<shell;s++){
                double dt0 = t_drifted[shell] + dt - t_drifted[s];
                double dt1 = dt; 
                double drift = reb_drift_from_straight_line(p0[s][mi],dt0,particles[mi],dt1);
                if (drift>maxdrift_dominant[s][mi]){
                    rebd_md3a[shell]++;
                    //printf("MAXDRIFT VIOLATION SUB-DOM triggered shell. checking particles. %2d %2d  %3d   %e   %e\n",shell,s,mi,drift,maxdrift_dominant[s][mi]);
                    maxdrift_dominant[s][mi] = 1e300;
                    for (int j=0; j<rim->shellN_dominant[s]; j++){ 
                        int mj = rim->map_dominant[s][j]; 
                        if (mi==mj) continue; // should not be possible here, (dom always != sub)
                        // Not yet in current shell?
                        unsigned int mj_shell = inshell_dominant[mj];
                        if (mj_shell<shell){
                            rebd_md3b[shell]++;
                            double drift = t_drifted[shell] - t_drifted[mj_shell];
                            double rmin2 = reb_mercurana_predict_rmin2_drifted(particles[mi],particles[mj],dt_drift[s],drift);
                            double dcritsum = buffer*(dcrit[mi]+dcrit[mj]);
                            if (rmin2< dcritsum*dcritsum){ 
                                rebd_md3c[shell]++;
                                // Add mj to all higher shells as dominant:
                                inshell_dominant[mj] = shell;
                                printf("MAXDRIFT VIOLATION SUB-DOM triggered. moving particle.\n");
                                rim->moved_particles++;
                                for (int sa=mj_shell+1; sa<=shell; sa++){
                                    rim->map_dominant[sa][rim->shellN_dominant[sa]] = mj;
                                    rim->shellN_dominant[sa]++;
                                }
                                particles[mj].x += drift*particles[mj].vx;
                                particles[mj].y += drift*particles[mj].vy;
                                particles[mj].z += drift*particles[mj].vz;
                                p0[shell][mj] = particles[mj]; 
                            }else{
                                double dcritsum = dcrit[mi]+dcrit[mj];
                                const double maxdrift = (sqrt(rmin2) - dcritsum)/2.;
                                maxdrift_dominant[s][mi] = MIN(maxdrift_dominant[s][mi],maxdrift);
                            }
                        }
                    }
                }
            }
        }
        // Encounter - Encounter
        for (int i=0; i<shellN_encounter; i++){
            int mi = map_encounter[i]; 
            for (int s=0;s<shell;s++){
                double dt0 = t_drifted[shell] + dt - t_drifted[s];
                double dt1 = dt; 
                double drift = reb_drift_from_straight_line(p0[s][mi],dt0,particles[mi],dt1);
                if (drift>maxdrift_encounter[s][mi]){
                    rebd_md4a[shell]++;
                    //printf("MAXDRIFT VIOLATION ENC-ENC triggered shell. checking particles. %2d %2d  %3d   %e   %e\n",shell,s,mi,drift,maxdrift_encounter[s][mi]);
                    maxdrift_encounter[s][mi] = 1e300;
                    for (int j=0; j<rim->shellN_encounter[s]; j++){ 
                        int mj = rim->map_encounter[s][j]; 
                        if (mi==mj) continue;
                        // Not yet in current shell?
                        unsigned int mj_shell = inshell_encounter[mj];
                        if (mj_shell<shell){
                            rebd_md4b[shell]++;
                            double drift = t_drifted[shell] - t_drifted[mj_shell];
                            double rmin2 = reb_mercurana_predict_rmin2_drifted(particles[mi],particles[mj],dt_drift[s],drift);
                            double dcritsum = buffer*(dcrit[mi]+dcrit[mj]);
                            if (rmin2< dcritsum*dcritsum){ 
                                rebd_md4c[shell]++;
                                // Add mj to all higher shells as encounter:
                                inshell_encounter[mj] = shell;
                                printf("MAXDRIFT VIOLATION ENC-ENC triggered. moving particle.\n");
                                rim->moved_particles++;
                                for (int sa=mj_shell+1; sa<=shell; sa++){
                                    rim->map_encounter[sa][rim->shellN_encounter[sa]] = mj;
                                    rim->shellN_encounter[sa]++;
                                }
                                particles[mj].x += drift*particles[mj].vx;
                                particles[mj].y += drift*particles[mj].vy;
                                particles[mj].z += drift*particles[mj].vz;
                                p0[shell][mj] = particles[mj]; 
                            }else{
                                const double maxdrift = (sqrt(rmin2) - dcritsum)/2.;
                                maxdrift_encounter[s][mi] = MIN(maxdrift_encounter[s][mi],maxdrift);
                            }
                        }
                    }
                }
            }
        }
    }

    // Check interactions in shell
    // Dominant and dominant
    for (int i=0; i<shellN_dominant; i++){
        int mi = map_dominant[i]; 
        for (int j=i+1; j<shellN_dominant; j++){
            int mj = map_dominant[j]; 
            double rmin2 = reb_mercurana_predict_rmin2(particles[mi],particles[mj],dt);
            double rsum = r->particles[mi].r+r->particles[mj].r;
            if (rmin2< rsum*rsum && r->collision==REB_COLLISION_DIRECT){
                reb_mercurana_record_collision(r,mi,mj);
            }
            double dcritsum = dcrit[mi]+dcrit[mj];
            if (rmin2< dcritsum*dcritsum){ 
                if (inshell_dominant[mi] == shell){
                    inshell_dominant[mi] = shell+1;
                    rim->map_dominant[shell+1][rim->shellN_dominant[shell+1]] = mi;
                    rim->shellN_dominant[shell+1]++;
                    p0[shell+1][mi] = particles[mi]; 
                }
                if (inshell_dominant[mj] == shell){
                    inshell_dominant[mj] = shell+1;
                    rim->map_dominant[shell+1][rim->shellN_dominant[shell+1]] = mj;
                    rim->shellN_dominant[shell+1]++;
                    p0[shell+1][mj] = particles[mj]; 
                }
            }
        }
    }
    
    // Dominant and subdominant
    for (int i=0; i<shellN_dominant; i++){
        int mi = map_dominant[i]; 
        for (int j=0; j<shellN_subdominant; j++){
            int mj = map_subdominant[j]; 
            double rmin2 = reb_mercurana_predict_rmin2(particles[mi],particles[mj],dt);
            double rsum = r->particles[mi].r+r->particles[mj].r;
            if (rmin2< rsum*rsum && r->collision==REB_COLLISION_DIRECT){
                reb_mercurana_record_collision(r,mi,mj);
            }
            double dcritsum = dcrit[mi]+dcrit[mj];
            if (rmin2< dcritsum*dcritsum){ 
                if (inshell_dominant[mi] == shell){
                    inshell_dominant[mi] = shell+1;
                    rim->map_dominant[shell+1][rim->shellN_dominant[shell+1]] = mi;
                    rim->shellN_dominant[shell+1]++;
                    p0[shell+1][mi] = particles[mi]; 
                }
                if (inshell_subdominant[mj] == shell){
                    inshell_subdominant[mj] = shell+1;
                    rim->map_subdominant[shell+1][rim->shellN_subdominant[shell+1]] = mj;
                    rim->shellN_subdominant[shell+1]++;
                    p0[shell+1][mj] = particles[mj]; 
                }
            }
        }
    }
    
    // Encounter and encounter
    for (int i=0; i<shellN_encounter; i++){
        int mi = map_encounter[i]; 
        for (int j=i+1; j<shellN_encounter; j++){
            int mj = map_encounter[j]; 
            double rmin2 = reb_mercurana_predict_rmin2(particles[mi],particles[mj],dt);
            double rsum = r->particles[mi].r+r->particles[mj].r;
            if (rmin2< rsum*rsum && r->collision==REB_COLLISION_DIRECT){
                reb_mercurana_record_collision(r,mi,mj);
            }
            double dcritsum = dcrit[mi]+dcrit[mj];
            if (rmin2< dcritsum*dcritsum){ 
                if (inshell_encounter[mi] == shell){
                    inshell_encounter[mi] = shell+1;
                    rim->map_encounter[shell+1][rim->shellN_encounter[shell+1]] = mi;
                    rim->shellN_encounter[shell+1]++;
                    p0[shell+1][mi] = particles[mi]; 
                }
                if (inshell_encounter[mj] == shell){
                    inshell_encounter[mj] = shell+1;
                    rim->map_encounter[shell+1][rim->shellN_encounter[shell+1]] = mj;
                    rim->shellN_encounter[shell+1]++;
                    p0[shell+1][mj] = particles[mj]; 
                }
            }
        }
    }
    
    // Maxdrift calculation
    for (int i=0; i<shellN_encounter; i++){
        int mi = map_encounter[i]; 
        maxdrift_encounter[shell][mi] = 1e300; 
        maxdrift_dominant[shell][mi] = 1e300; 
        //TODO think about following
        p0[shell][mi] = particles[mi]; 
    }
    // Maxdrift calculation
    for (int i=0; i<shellN_dominant; i++){
        int mi = map_dominant[i]; 
        maxdrift_encounter[shell][mi] = 1e300; 
        maxdrift_dominant[shell][mi] = 1e300; 
        //TODO think about following
        p0[shell][mi] = particles[mi]; 
    }
    // Maxdrift calculation
    for (int i=0; i<shellN_subdominant; i++){
        int mi = map_subdominant[i]; 
        maxdrift_encounter[shell][mi] = 1e300; 
        maxdrift_dominant[shell][mi] = 1e300; 
        //TODO think about following
        p0[shell][mi] = particles[mi]; 
    }
    
    // Encounter and encounter
    for (int i=0; i<shellN_encounter; i++){
        int mi = map_encounter[i]; 
        for (int j=i+1; j<shellN_encounter; j++){
            int mj = map_encounter[j]; 
            // Interaction not resolved in subshells?
            if ((inshell_encounter[mi] == shell) || (inshell_encounter[mj] == shell)){
                double rmin2 = reb_mercurana_predict_rmin2(particles[mi],particles[mj],dt);
                double dcritsum = dcrit[mi]+dcrit[mj];
                const double maxdrift = (sqrt(rmin2) - dcritsum)/2.;
                maxdrift_encounter[shell][mi] = MIN(maxdrift_encounter[shell][mi],maxdrift);
                maxdrift_encounter[shell][mj] = MIN(maxdrift_encounter[shell][mj],maxdrift);
            }
        }
    }
    // Dominant and dominant
    for (int i=0; i<shellN_dominant; i++){
        int mi = map_dominant[i]; 
        for (int j=i+1; j<shellN_dominant; j++){
            int mj = map_dominant[j]; 
            // Interaction not resolved in subshells?
            if ((inshell_dominant[mi] == shell) || (inshell_dominant[mj] == shell)){
                double rmin2 = reb_mercurana_predict_rmin2(particles[mi],particles[mj],dt);
                double dcritsum = dcrit[mi]+dcrit[mj];
                const double maxdrift = (sqrt(rmin2) - dcritsum)/2.;
                maxdrift_dominant[shell][mi] = MIN(maxdrift_dominant[shell][mi],maxdrift);
                maxdrift_dominant[shell][mj] = MIN(maxdrift_dominant[shell][mj],maxdrift);
            }
        }
    }
    // Subdominant and dominant
    for (int i=0; i<shellN_dominant; i++){
        int mi = map_dominant[i]; 
        for (int j=0; j<shellN_subdominant; j++){
            int mj = map_subdominant[j]; 
            // Interaction not resolved in subshells?
            if ((inshell_dominant[mi] == shell) || (inshell_subdominant[mj] == shell)){
                double rmin2 = reb_mercurana_predict_rmin2(particles[mi],particles[mj],dt);
                double dcritsum = dcrit[mi]+dcrit[mj];
                const double maxdrift = (sqrt(rmin2) - dcritsum)/2.;
                maxdrift_encounter[shell][mi] = MIN(maxdrift_encounter[shell][mi],maxdrift);
                maxdrift_dominant[shell][mj] = MIN(maxdrift_dominant[shell][mj],maxdrift);
            }
        }
    }
    
    
    if (rim->collisions_N){
        printf("collision\n");
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
    rebd_kick[shell]++;
    struct reb_simulation_integrator_mercurana* const rim = &(r->ri_mercurana);
    struct reb_particle* const particles = r->particles;
    r->gravity = REB_GRAVITY_MERCURANA; // needed here again for SimulationArchive
    rim->current_shell = shell;
    reb_update_acceleration(r);
    if (v!=0.){
        reb_calculate_and_apply_jerk(r,v);
    }
    unsigned int* map_encounter = rim->map_encounter[shell];
    unsigned int* map_dominant = rim->map_dominant[shell];
    unsigned int* map_subdominant = rim->map_subdominant[shell];
    int shellN_encounter = rim->shellN_encounter[shell];
    int shellN_dominant = rim->shellN_dominant[shell];
    int shellN_subdominant = rim->shellN_subdominant[shell];
    unsigned int* inshell_encounter = rim->inshell_encounter;

    for (int i=0;i<shellN_dominant;i++){ // Apply acceleration. Jerk already applied.
        const int mi = map_dominant[i];
        particles[mi].vx += y*particles[mi].ax;
        particles[mi].vy += y*particles[mi].ay;
        particles[mi].vz += y*particles[mi].az;
    }
    for (int i=0;i<shellN_encounter;i++){ // Apply acceleration. Jerk already applied.
        const int mi = map_encounter[i];
        particles[mi].vx += y*particles[mi].ax;
        particles[mi].vy += y*particles[mi].ay;
        particles[mi].vz += y*particles[mi].az;
    }
    if (shell>0){ // All particles are encounter particles in shell 0, no need for subdominant kick
    for (int i=0;i<shellN_subdominant;i++){ // Apply acceleration. Jerk already applied.
        const int mi = map_subdominant[i];
        if (inshell_encounter[mi]<shell){ // do not apply acceleration twice
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
    unsigned int* map_dominant = rim->map_dominant[shell];
    unsigned int* map_subdominant = rim->map_subdominant[shell];
    int shellN_encounter = rim->shellN_encounter[shell];
    int shellN_dominant = rim->shellN_dominant[shell];
    int shellN_subdominant = rim->shellN_subdominant[shell];
    unsigned int* inshell_encounter = rim->inshell_encounter;
    unsigned int* inshell_dominant = rim->inshell_dominant;
    unsigned int* inshell_subdominant = rim->inshell_dominant;
    
    if (shell+1<rim->Nmaxshells){ // does sub-shell exist? If so, do that first.
        // Are there particles in it?
        if (rim->shellN_encounter[shell+1]>0 || rim->shellN_dominant[shell+1]>0){
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
    
    for (int i=0;i<shellN_dominant;i++){  // loop over all particles in shell (includes subshells)
        int mi = map_dominant[i]; 
        if( inshell_dominant[mi]==shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
        }
    }
    for (int i=0;i<shellN_subdominant;i++){  // loop over all particles in shell (includes subshells)
        int mi = map_subdominant[i]; 
        if( inshell_subdominant[mi]==shell && inshell_encounter[mi]<=shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
        }
    }
    for (int i=0;i<shellN_encounter;i++){  // loop over all particles in shell (includes subshells)
        int mi = map_encounter[i]; 
        if( inshell_subdominant[mi]<shell && inshell_encounter[mi]==shell){
            particles[mi].x += a*particles[mi].vx;
            particles[mi].y += a*particles[mi].vy;
            particles[mi].z += a*particles[mi].vz;
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
    if (rim->N_dominant>r->N_active && r->N_active!=-1){
        reb_error(r,"kappa>0 is required if Nmaxshells>1.");
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
        rim->t_drifted = realloc(rim->t_drifted, sizeof(double)*rim->Nmaxshells);
        rim->dt_drift = realloc(rim->dt_drift, sizeof(double)*(rim->Nmaxshells));
        rim->dcrit = realloc(rim->dcrit, sizeof(double*)*(rim->Nmaxshells));
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->dcrit[i] = malloc(sizeof(double)*N);
        }
        // map
        if (rim->map_encounter){
            for (int i=0;i<rim->Nmaxshells;i++){
                free(rim->map_encounter[i]);
                free(rim->map_dominant[i]);
                free(rim->map_subdominant[i]);
                free(rim->maxdrift_encounter[i]);
                free(rim->maxdrift_dominant[i]);
                free(rim->p0[i]);
            }
        }
        rim->map_encounter = realloc(rim->map_encounter, sizeof(unsigned int*)*rim->Nmaxshells);
        rim->map_dominant = realloc(rim->map_dominant, sizeof(unsigned int*)*rim->Nmaxshells);
        rim->map_subdominant = realloc(rim->map_subdominant, sizeof(unsigned int*)*rim->Nmaxshells);
        rim->maxdrift_encounter = realloc(rim->maxdrift_encounter, sizeof(double*)*rim->Nmaxshells);
        rim->maxdrift_dominant = realloc(rim->maxdrift_dominant, sizeof(double*)*rim->Nmaxshells);
        rim->p0 = realloc(rim->p0, sizeof(struct reb_particle*)*rim->Nmaxshells);
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->map_encounter[i] = malloc(sizeof(unsigned int)*N);
            rim->map_dominant[i] = malloc(sizeof(unsigned int)*N);
            rim->map_subdominant[i] = malloc(sizeof(unsigned int)*N);
            rim->maxdrift_encounter[i] = malloc(sizeof(double)*N);
            rim->maxdrift_dominant[i] = malloc(sizeof(double)*N);
            rim->p0[i] = malloc(sizeof(struct reb_particle)*N);
        }
        // inshell
        rim->inshell_encounter = realloc(rim->inshell_encounter, sizeof(unsigned int)*N);
        rim->inshell_dominant = realloc(rim->inshell_dominant, sizeof(unsigned int)*N);
        rim->inshell_subdominant = realloc(rim->inshell_subdominant, sizeof(unsigned int)*N);
        // shellN
        rim->shellN_encounter = realloc(rim->shellN_encounter, sizeof(unsigned int)*rim->Nmaxshells);
        rim->shellN_dominant = realloc(rim->shellN_dominant, sizeof(unsigned int)*rim->Nmaxshells);
        rim->shellN_subdominant = realloc(rim->shellN_subdominant, sizeof(unsigned int)*rim->Nmaxshells);
        for (int i=0;i<rim->Nmaxshells;i++){
            rim->shellN_encounter[i] = 0;
            rim->shellN_dominant[i] = 0;
            rim->shellN_subdominant[i] = 0;
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
            rim->dt_drift[s] = dt_shell;
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
    // done in predict now
    //rim->shellN[0] = r->N;
    //rim->shellN_active[0] = r->N_active==-1?r->N:r->N_active;
    
    for (int i=0;i<r->ri_mercurana.Nmaxshells;i++){
        rim->t_drifted[i] =0;
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
            free(r->ri_mercurana.map_dominant[i]);
            free(r->ri_mercurana.map_subdominant[i]);
            free(r->ri_mercurana.dcrit[i]);
            free(r->ri_mercurana.maxdrift_encounter[i]);
            free(r->ri_mercurana.maxdrift_dominant[i]);
            free(r->ri_mercurana.p0[i]);
        }
        free(r->ri_mercurana.map_encounter);
        free(r->ri_mercurana.map_dominant);
        free(r->ri_mercurana.map_subdominant);
        free(r->ri_mercurana.dcrit);
        free(r->ri_mercurana.inshell_encounter);
        free(r->ri_mercurana.inshell_dominant);
        free(r->ri_mercurana.inshell_subdominant);
        free(r->ri_mercurana.shellN_encounter);
        free(r->ri_mercurana.shellN_dominant);
        free(r->ri_mercurana.shellN_subdominant);
        free(r->ri_mercurana.t_drifted);
        free(r->ri_mercurana.maxdrift_encounter);
        free(r->ri_mercurana.maxdrift_dominant);
        free(r->ri_mercurana.p0);
    }
    r->ri_mercurana.allocatedN = 0;
    r->ri_mercurana.map_encounter = NULL;
    r->ri_mercurana.map_dominant = NULL;
    r->ri_mercurana.map_subdominant = NULL;
    r->ri_mercurana.dcrit = NULL;
    r->ri_mercurana.inshell_encounter = NULL;
    r->ri_mercurana.inshell_dominant = NULL;
    r->ri_mercurana.inshell_subdominant = NULL;
    r->ri_mercurana.shellN_encounter = NULL;
    r->ri_mercurana.shellN_dominant = NULL;
    r->ri_mercurana.shellN_subdominant = NULL;
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

