/**
 * @file    integrator_mercurius.c
 * @brief   MERCURIUS, a modified version of John Chambers' MERCURY algorithm
 *          using the IAS15 integrator and WHFast. It works with planet-planry
 *          collisions, test particles, and additional forces.
 * @author  Hanno Rein, Dan Tamayo
 * 
 * @section LICENSE
 * Copyright (c) 2019 Hanno Rein, Dan Tamayo
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
#include <math.h>
#include <string.h>
#include "rebound.h"
#include "integrator.h"
#include "gravity.h"
#include "integrator_mercurius.h"
#include "integrator_ias15.h"
#include "integrator_whfast.h"
#include "collision.h"
#include "input.h"
#include "output.h"

#define MIN(a, b) ((a) > (b) ? (b) : (a))    ///< Returns the minimum of a and b
#define MAX(a, b) ((a) > (b) ? (a) : (b))    ///< Returns the maximum of a and b

void reb_integator_mercurius_calculate_acceleration_ias15_part(struct reb_simulation* r);
void reb_integator_mercurius_calculate_acceleration_whfast_part(struct reb_simulation* r);

double reb_integrator_mercurius_L_mercury(const struct reb_simulation* const r, double d, double dcrit){
    // This is the changeover function used by the Mercury integrator.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return 10.*(y*y*y) - 15.*(y*y*y*y) + 6.*(y*y*y*y*y);
    }
}

static double f(double x){
    if (x<0) return 0;
    return exp(-1./x);
}

double reb_integrator_mercurius_L_infinity(const struct reb_simulation* const r, double d, double dcrit){
    // Infinitely differentiable function.
    double y = (d-0.1*dcrit)/(0.9*dcrit);
    if (y<0.){
        return 0.;
    }else if (y>1.){
        return 1.;
    }else{
        return f(y) /(f(y) + f(1.-y));
    }
}


void reb_integrator_mercurius_inertial_to_dh(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_vec3d com_pos = {0};
    struct reb_vec3d com_vel = {0};
    double mtot = 0.;
    const int N_active = (r->N_active==-1 || r->testparticle_type==1)?r->N:r->N_active;
    const int N = r->N;
    for (int i=0;i<N_active;i++){
        double m = particles[i].m;
        com_pos.x += m * particles[i].x;
        com_pos.y += m * particles[i].y;
        com_pos.z += m * particles[i].z;
        com_vel.x += m * particles[i].vx;
        com_vel.y += m * particles[i].vy;
        com_vel.z += m * particles[i].vz;
        mtot += m; 
    }
    com_pos.x /= mtot; com_pos.y /= mtot; com_pos.z /= mtot;
    com_vel.x /= mtot; com_vel.y /= mtot; com_vel.z /= mtot;
    // Particle 0 is also changed to allow for easy collision detection
    for (int i=N-1;i>=0;i--){ 
        particles[i].x -= particles[0].x;
        particles[i].y -= particles[0].y;
        particles[i].z -= particles[0].z;
        particles[i].vx -= com_vel.x;
        particles[i].vy -= com_vel.y;
        particles[i].vz -= com_vel.z;
    }
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    config->com_pos = com_pos;
    config->com_vel = com_vel;
}

void reb_integrator_mercurius_dh_to_inertial(struct reb_simulation* r){
    struct reb_particle* restrict const particles = r->particles;
    struct reb_particle temp = {0};
    const int N = r->N;
    const int N_active = (r->N_active==-1 || r->testparticle_type==1)?r->N:r->N_active;
    for (int i=1;i<N_active;i++){
        double m = particles[i].m;
        temp.x += m * particles[i].x;
        temp.y += m * particles[i].y;
        temp.z += m * particles[i].z;
        temp.vx += m * particles[i].vx;
        temp.vy += m * particles[i].vy;
        temp.vz += m * particles[i].vz;
        temp.m += m;
    }
    temp.m += r->particles[0].m;
    temp.x /= temp.m; 
    temp.y /= temp.m;
    temp.z /= temp.m;
    temp.vx /= particles[0].m; 
    temp.vy /= particles[0].m;
    temp.vz /= particles[0].m;
    // Use com to calculate central object's position.
    // This ignores previous values stored in particles[0].
    // Should not matter unless collisions occured.
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    particles[0].x = config->com_pos.x - temp.x; 
    particles[0].y = config->com_pos.y - temp.y; 
    particles[0].z = config->com_pos.z - temp.z; 

    for (int i=1;i<N;i++){
        particles[i].x += particles[0].x;
        particles[i].y += particles[0].y;
        particles[i].z += particles[0].z;
        particles[i].vx += config->com_vel.x;
        particles[i].vy += config->com_vel.y;
        particles[i].vz += config->com_vel.z;
    }
    particles[0].vx = config->com_vel.x - temp.vx; 
    particles[0].vy = config->com_vel.y - temp.vy; 
    particles[0].vz = config->com_vel.z - temp.vz; 
}


static void reb_mercurius_encounter_predict(struct reb_simulation* const r){
    // This function predicts close encounters during the timestep
    // It makes use of the old and new position and velocities obtained
    // after the Kepler step.
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    struct reb_particle* const particles = r->particles;
    struct reb_particle* const particles_backup = config->particles_backup;
    const double* const dcrit = config->dcrit;
    const int N = r->N;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    const double dt = r->dt;
    config->encounterN = 1;
    config->encounter_map[0] = 1;
    if (r->testparticle_type==1){
        config->tponly_encounter = 0; // testparticles affect massive particles
    }else{
        config->tponly_encounter = 1;
    }
    for (int i=1; i<N; i++){
        config->encounter_map[i] = 0;
    }
    for (int i=0; i<N_active; i++){
        for (int j=i+1; j<N; j++){
            const double dxn = particles[i].x - particles[j].x;
            const double dyn = particles[i].y - particles[j].y;
            const double dzn = particles[i].z - particles[j].z;
            const double dvxn = particles[i].vx - particles[j].vx;
            const double dvyn = particles[i].vy - particles[j].vy;
            const double dvzn = particles[i].vz - particles[j].vz;
            const double rn = (dxn*dxn + dyn*dyn + dzn*dzn);
            const double dxo = particles_backup[i].x - particles_backup[j].x;
            const double dyo = particles_backup[i].y - particles_backup[j].y;
            const double dzo = particles_backup[i].z - particles_backup[j].z;
            const double dvxo = particles_backup[i].vx - particles_backup[j].vx;
            const double dvyo = particles_backup[i].vy - particles_backup[j].vy;
            const double dvzo = particles_backup[i].vz - particles_backup[j].vz;
            const double ro = (dxo*dxo + dyo*dyo + dzo*dzo);

            const double drndt = (dxn*dvxn+dyn*dvyn+dzn*dvzn)*2.;
            const double drodt = (dxo*dvxo+dyo*dvyo+dzo*dvzo)*2.;

            const double a = 6.*(ro-rn)+3.*dt*(drodt+drndt); 
            const double b = 6.*(rn-ro)-2.*dt*(2.*drodt+drndt); 
            const double c = dt*drodt; 

            double rmin = MIN(rn,ro);

            const double s = b*b-4.*a*c;
            const double sr = sqrt(MAX(0.,s));
            const double tmin1 = (-b + sr)/(2.*a); 
            const double tmin2 = (-b - sr)/(2.*a); 
            if (tmin1>0. && tmin1<1.){
                const double rmin1 = (1.-tmin1)*(1.-tmin1)*(1.+2.*tmin1)*ro
                                     + tmin1*tmin1*(3.-2.*tmin1)*rn
                                     + tmin1*(1.-tmin1)*(1.-tmin1)*dt*drodt
                                     - tmin1*tmin1*(1.-tmin1)*dt*drndt;
                rmin = MIN(MAX(rmin1,0.),rmin);
            }
            if (tmin2>0. && tmin2<1.){
                const double rmin2 = (1.-tmin2)*(1.-tmin2)*(1.+2.*tmin2)*ro
                                     + tmin2*tmin2*(3.-2.*tmin2)*rn
                                     + tmin2*(1.-tmin2)*(1.-tmin2)*dt*drodt
                                     - tmin2*tmin2*(1.-tmin2)*dt*drndt;
                rmin = MIN(MAX(rmin2,0.),rmin);
            }

            double dcritmax2 = MAX(dcrit[i],dcrit[j]);
            dcritmax2 *= 1.21*dcritmax2;
            if (rmin < dcritmax2){
                if (config->encounter_map[i]==0){
                    config->encounter_map[i] = i;
                    config->encounterN++;
                }
                if (config->encounter_map[j]==0){
                    config->encounter_map[j] = j;
                    config->encounterN++;
                }
                if (j<N_active){ // Two massive particles have a close encounter
                    config->tponly_encounter = 0;
                }
            }
        }
    }
}
    
void reb_integrator_mercurius_interaction_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    for (int i=1;i<N;i++){
        particles[i].vx += dt*particles[i].ax;
        particles[i].vy += dt*particles[i].ay;
        particles[i].vz += dt*particles[i].az;
    }
}

void reb_integrator_mercurius_jump_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N_active = r->N_active==-1?r->N:r->N_active;
    const int N = r->testparticle_type==0 ? N_active: r->N;
    double px=0., py=0., pz=0.;
    for (int i=1;i<N;i++){
        px += r->particles[i].vx*r->particles[i].m; // in dh
        py += r->particles[i].vy*r->particles[i].m; 
        pz += r->particles[i].vz*r->particles[i].m;
    }
    px /= r->particles[0].m;
    py /= r->particles[0].m;
    pz /= r->particles[0].m;
    for (int i=1;i<N;i++){
        particles[i].x += dt*px;
        particles[i].y += dt*py;
        particles[i].z += dt*pz;
    }
}

void reb_integrator_mercurius_com_step(struct reb_simulation* const r, double dt){
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    config->com_pos.x += dt*config->com_vel.x;
    config->com_pos.y += dt*config->com_vel.y;
    config->com_pos.z += dt*config->com_vel.z;
}

void reb_integrator_mercurius_kepler_step(struct reb_simulation* const r, double dt){
    struct reb_particle* restrict const particles = r->particles;
    const int N = r->N;
    for (int i=1;i<N;i++){
        reb_whfast_kepler_solver(r,particles,r->G*particles[0].m,i,dt); // in dh
    }
}

static void reb_mercurius_encounter_step(struct reb_simulation* const r, const double _dt){
    // Only particles having a close encounter are integrated by IAS15.
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    if (config->encounterN<2){
        return; // If there are no particles (other than the star) having a close encounter, then there is nothing to do.
    }

    int i_enc = 0;
    config->encounterNactive = 0;
    for (unsigned int i=0; i<r->N; i++){
        if(config->encounter_map[i]){  
            struct reb_particle tmp = r->particles[i];      // Copy for potential use for tponly_encounter
            r->particles[i] = config->particles_backup[i];     // Use coordinates before whfast step
            config->encounter_map[i_enc] = i;
            i_enc++;
            if (r->N_active==-1 || i<r->N_active){
                config->encounterNactive++;
                if (config->tponly_encounter){
                    config->particles_backup[i] = tmp;         // Make copy of particles after the kepler step.
                                                            // used to restore the massive objects' states in the case
                                                            // of only massless test-particle encounters
                }
            }
        }
    }

    config->mode = 1;
    
    // run
    const double old_dt = r->dt;
    const double old_t = r->t;
    double t_needed = r->t + _dt; 
        
    printf("Need to do ias15 reset");
    exit(EXIT_FAILURE);
    // reb_integrator_ias15_reset(r); // TODO
    
    r->dt = 0.0001*_dt; // start with a small timestep.
    
    while(r->t < t_needed && fabs(r->dt/old_dt)>1e-14 ){
        struct reb_particle star = r->particles[0]; // backup velocity
        r->particles[0].vx = 0; // star does not move in dh 
        r->particles[0].vy = 0;
        r->particles[0].vz = 0;
        reb_integrator_ias15_step(NULL, r);  // TODO  need to do differently. first argumnet should not be null
        r->particles[0].vx = star.vx; // restore every timestep for collisions
        r->particles[0].vy = star.vy;
        r->particles[0].vz = star.vz;
        
        if (r->t+r->dt >  t_needed){
            r->dt = t_needed-r->t;
        }

        // Search and resolve collisions
        reb_collision_search(r);

        // Do any additional post_timestep_modifications.
        // Note: post_timestep_modifications is called here but also
        // at the end of the full timestep. The function thus needs
        // to be implemented with care as not to do the same 
        // modification multiple times. To do that, check the value of
        // r->ri_mercurius.mode
        if (r->post_timestep_modifications){
            r->post_timestep_modifications(r);
        }

        star.vx = r->particles[0].vx; // keep track of changed star velocity for later collisions
        star.vy = r->particles[0].vy;
        star.vz = r->particles[0].vz;
        if (r->particles[0].x !=0 || r->particles[0].y !=0 || r->particles[0].z !=0){
            // Collision with star occured
            // Shift all particles back to heliocentric coordinates
            // Ignore stars velocity:
            //   - will not be used after this
            //   - com velocity is unchained. this velocity will be used
            //     to reconstruct star's velocity later.
            for (int i=r->N-1; i>=0; i--){
                r->particles[i].x -= r->particles[0].x;
                r->particles[i].y -= r->particles[0].y;
                r->particles[i].z -= r->particles[0].z;
            }
        }
    }

    // if only test particles encountered massive bodies, reset the
    // massive body coordinates to their post Kepler step state
    if(config->tponly_encounter){
        for (int i=1;i<config->encounterNactive;i++){
            unsigned int mi = config->encounter_map[i];
            r->particles[mi] = config->particles_backup[mi];
        }
    }

    // Reset constant for global particles
    r->t = old_t;
    r->dt = old_dt;
    config->mode = 0;

}

void reb_integrator_mercurius_synchronize(struct reb_integrator* integrator, struct reb_simulation* r){
    struct reb_integrator_mercurius_config* const config = (struct reb_integrator_mercurius_config*) &(integrator->config);
    if (config->is_synchronized == 0){
        r->gravity = REB_GRAVITY_MERCURIUS; // needed here again for SimulationArchive
        config->mode = 0;
        if (config->L == NULL){
            // Setting default switching function
            config->L = reb_integrator_mercurius_L_mercury;
        }
        reb_integator_mercurius_calculate_acceleration_whfast_part(r);
        reb_integrator_mercurius_interaction_step(r,r->dt/2.);
        
        reb_integrator_mercurius_dh_to_inertial(r);

        config->recalculate_coordinates_this_timestep = 1; 
        config->is_synchronized = 1;
    }
}



double reb_integrator_mercurius_calculate_dcrit_for_particle(struct reb_simulation* r, unsigned int i){
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    const double m0 = r->particles[0].m;
    const double dx  = r->particles[i].x;  // in dh
    const double dy  = r->particles[i].y;
    const double dz  = r->particles[i].z;
    const double dvx = r->particles[i].vx - r->particles[0].vx; 
    const double dvy = r->particles[i].vy - r->particles[0].vy; 
    const double dvz = r->particles[i].vz - r->particles[0].vz; 
    const double _r = sqrt(dx*dx + dy*dy + dz*dz);
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;

    const double GM = r->G*(m0+r->particles[i].m);
    const double a = GM*_r / (2.*GM - _r*v2);
    const double vc = sqrt(GM/fabs(a));
    double dcrit = 0;
    // Criteria 1: average velocity
    dcrit = MAX(dcrit, vc*0.4*r->dt);
    // Criteria 2: current velocity
    dcrit = MAX(dcrit, sqrt(v2)*0.4*r->dt);
    // Criteria 3: Hill radius
    dcrit = MAX(dcrit, config->hillfac*a*cbrt(r->particles[i].m/(3.*r->particles[0].m)));
    // Criteria 4: physical radius
    dcrit = MAX(dcrit, 2.*r->particles[i].r);
    return dcrit;
}


void reb_integrator_mercurius_step(struct reb_integrator* integrator, struct reb_simulation* r){
    if (r->var_config_N){
        reb_warning(r,"Mercurius does not work with variational equations.");
    }
    
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    const int N = r->N;
    
    if (config->dcrit_allocatedN<N){
        // Need to safe these arrays in SimulationArchive
        config->dcrit              = realloc(config->dcrit, sizeof(double)*N);
        config->dcrit_allocatedN = N;
        // If particle number increased (or this is the first step), need to calculate critical radii
        config->recalculate_dcrit_this_timestep        = 1;
        // Heliocentric coordinates were never calculated.
        // This will get triggered on first step only (not when loaded from archive)
        config->recalculate_coordinates_this_timestep = 1;
    }
    if (config->allocatedN<N){
        // These arrays are only used within one timestep. 
        // Can be recreated without loosing bit-wise reproducibility
        config->particles_backup   = realloc(config->particles_backup,sizeof(struct reb_particle)*N);
        config->encounter_map      = realloc(config->encounter_map,sizeof(int)*N);
        config->allocatedN = N;
    }
    if (config->safe_mode || config->recalculate_coordinates_this_timestep){
        if (config->is_synchronized==0){
            reb_integrator_mercurius_synchronize(r->integrator_selected, r);
            reb_warning(r,"MERCURIUS: Recalculating heliocentric coordinates but coordinates were not synchronized before.");
        }
        reb_integrator_mercurius_inertial_to_dh(r);
        config->recalculate_coordinates_this_timestep = 0;
    }

    if (config->recalculate_dcrit_this_timestep){
        config->recalculate_dcrit_this_timestep = 0;
        if (config->is_synchronized==0){
            reb_integrator_mercurius_synchronize(r->integrator_selected, r);
            reb_integrator_mercurius_inertial_to_dh(r);
            config->recalculate_coordinates_this_timestep = 0;
            reb_warning(r,"MERCURIUS: Recalculating dcrit but pos/vel were not synchronized before.");
        }
        config->dcrit[0] = 2.*r->particles[0].r; // central object only uses physical radius
        for (int i=1;i<N;i++){
            config->dcrit[i] = reb_integrator_mercurius_calculate_dcrit_for_particle(r, i);
        }
    }
    
    // Calculate collisions only with DIRECT method
    if (r->collision != REB_COLLISION_NONE && r->collision != REB_COLLISION_DIRECT){
        reb_warning(r,"Mercurius only works with a direct collision search.");
    }
    
    // Calculate gravity with special function
    if (r->gravity != REB_GRAVITY_BASIC && r->gravity != REB_GRAVITY_MERCURIUS){
        reb_warning(r,"Mercurius has it's own gravity routine. Gravity routine set by the user will be ignored.");
    }
    r->gravity = REB_GRAVITY_MERCURIUS;
    config->mode = 0;
    
    if (config->L == NULL){
        // Setting default switching function
        config->L = reb_integrator_mercurius_L_mercury;
    }

    reb_integator_mercurius_calculate_acceleration_whfast_part(r);
   
    if (config->is_synchronized){
        reb_integrator_mercurius_interaction_step(r,r->dt/2.);
    }else{
        reb_integrator_mercurius_interaction_step(r,r->dt);
    }
    reb_integrator_mercurius_jump_step(r,r->dt/2.);
    reb_integrator_mercurius_com_step(r,r->dt); 
    
    // Make copy of particles before the kepler step.
    // Then evolve all particles in kepler step.
    // Result will be used in encounter prediction.
    // Particles having a close encounter will be overwritten 
    // later by encounter step.
    memcpy(config->particles_backup,r->particles,N*sizeof(struct reb_particle)); 
    reb_integrator_mercurius_kepler_step(r,r->dt);

    reb_mercurius_encounter_predict(r);
   
    reb_mercurius_encounter_step(r,r->dt);
    
    reb_integrator_mercurius_jump_step(r,r->dt/2.);
        
    config->is_synchronized = 0;
    if (config->safe_mode){
        reb_integrator_mercurius_synchronize(r->integrator_selected, r);
    }

    r->t+=r->dt;
    r->dt_last_done = r->dt;
}

    
enum MERCURIUS_CONFIG {
    HILLFAC = 1,
    SAFEMODE = 2,
    ISSYNCHRON = 3,
    DCRIT = 4,
    COMPOS = 5,
    COMVEL = 6,
};

size_t reb_integrator_mercurius_load(struct reb_integrator* integrator, struct reb_simulation* r, struct reb_input_stream* stream, struct reb_binary_field field){
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)integrator->config;
    switch (field.type){
        case REB_BF(MERCURIUS, HILLFAC):
            return reb_input_stream_fread(stream, &config->hillfac, field.size, 1);
        case REB_BF(MERCURIUS, SAFEMODE):
            return reb_input_stream_fread(stream, &config->safe_mode, field.size, 1);
        case REB_BF(MERCURIUS, ISSYNCHRON):
            return reb_input_stream_fread(stream, &config->is_synchronized, field.size, 1);
        case REB_BF(MERCURIUS, COMVEL):
            return reb_input_stream_fread(stream, &config->com_vel, field.size, 1);
        case REB_BF(MERCURIUS, COMPOS):
            return reb_input_stream_fread(stream, &config->com_pos, field.size, 1);
        case REB_BF(MERCURIUS, DCRIT):
            if(config->dcrit){
                free(config->dcrit);
            }
            config->dcrit_allocatedN = (int)(field.size/sizeof(double));
            if (field.size){
                config->dcrit = malloc(field.size);
                reb_input_stream_fread(stream, config->dcrit, field.size,1);
            }
            break;
        default:
            return 0;
    }
    return 0;
}

void reb_integrator_mercurius_save(struct reb_integrator* integrator, struct reb_simulation* r, struct reb_output_stream* stream){
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)integrator->config;
    reb_output_stream_write_field(stream, REB_BF(MERCURIUS, HILLFAC),    &config->hillfac,           sizeof(double));
    reb_output_stream_write_field(stream, REB_BF(MERCURIUS, SAFEMODE),   &config->safe_mode,         sizeof(unsigned int));
    reb_output_stream_write_field(stream, REB_BF(MERCURIUS, ISSYNCHRON), &config->is_synchronized,   sizeof(unsigned int));
    reb_output_stream_write_field(stream, REB_BF(MERCURIUS, DCRIT),      config->dcrit,              sizeof(double)*config->dcrit_allocatedN);
    reb_output_stream_write_field(stream, REB_BF(MERCURIUS, COMPOS),     &(config->com_pos),         sizeof(struct reb_vec3d));
    reb_output_stream_write_field(stream, REB_BF(MERCURIUS, COMVEL),     &(config->com_vel),         sizeof(struct reb_vec3d));

}

void* reb_integrator_mercurius_alloc(struct reb_integrator* integrator, struct reb_simulation* r){
    struct reb_integrator_mercurius_config* config = calloc(1, sizeof(struct reb_integrator_mercurius_config));
    config->hillfac = 3;
    return config;
}

void reb_integrator_mercurius_free(struct reb_integrator* integrator, struct reb_simulation* r){
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)integrator->config;
    if (config->dcrit){
        free(config->dcrit);
        config->dcrit = NULL;
    }
    if (config->particles_backup){
        free(config->particles_backup);
        config->particles_backup = NULL;
    }
    if (config->particles_backup_additionalforces){
        free(config->particles_backup_additionalforces);
        config->particles_backup_additionalforces = NULL;
    }
    if (config->particles_backup_additionalforces){
        free(config->particles_backup_additionalforces);
        config->particles_backup_additionalforces = NULL;
    }
    free(config);
    integrator->config = NULL;
}
    
void reb_integrator_mercurius_particle_add(struct reb_simulation* r, unsigned int i){
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    if (config->mode==0){ //WHFast part
        config->recalculate_dcrit_this_timestep       = 1;
        config->recalculate_coordinates_this_timestep = 1;
    }else{  // IAS15 part
        printf("Need to do ias15 reset");
        exit(EXIT_FAILURE);
        // reb_integrator_ias15_reset(r); // TODO
        if (config->dcrit_allocatedN<r->N){
            config->dcrit              = realloc(config->dcrit, sizeof(double)*r->N);
            config->dcrit_allocatedN = r->N;
        }
        config->dcrit[r->N-1] = reb_integrator_mercurius_calculate_dcrit_for_particle(r,r->N-1);
        if (config->allocatedN<r->N){
            config->particles_backup   = realloc(config->particles_backup,sizeof(struct reb_particle)*r->N);
            config->encounter_map      = realloc(config->encounter_map,sizeof(int)*r->N);
            config->allocatedN = r->N;
        }
        config->encounter_map[config->encounterN] = r->N-1;
        config->encounterN++;
        if (r->N_active==-1){ 
            // If global N_active is not set, then all particles are active, so the new one as well.
            // Otherwise, assume we're adding non active particle. 
            config->encounterNactive++;
        }
    }
}

void reb_integrator_mercurius_particle_remove(struct reb_simulation* r, unsigned int index, unsigned int keepSorted){
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    if (config->dcrit_allocatedN>0 && index<config->dcrit_allocatedN){
        for (int i=0;i<r->N-1;i++){
            if (i>=index){
                config->dcrit[i] = config->dcrit[i+1];
            }
        }
    }
    printf("Need to do ias15 reset");
    exit(EXIT_FAILURE);
    // reb_integrator_ias15_reset(r); // TODO
    if (config->mode==1){
        int after_to_be_removed_particle = 0;
        int encounter_index = -1;
        for (int i=0;i<config->encounterN;i++){
            if (after_to_be_removed_particle == 1){
                config->encounter_map[i-1] = config->encounter_map[i] - 1; 
            }
            if (config->encounter_map[i]==index){
                encounter_index = i;
                after_to_be_removed_particle = 1;
            }
        }
        if (encounter_index<config->encounterNactive){
            config->encounterNactive--;
        }
        config->encounterN--;
    }
}

void reb_integrator_mercurius_register(struct reb_simulation* r){
    struct reb_integrator* integrator = reb_simulation_register_integrator(r, "mercurius", 9);
    integrator->step        = reb_integrator_mercurius_step;
    integrator->synchronize = reb_integrator_mercurius_synchronize;
    integrator->alloc       = reb_integrator_mercurius_alloc;
    integrator->free        = reb_integrator_mercurius_free;
    integrator->load        = reb_integrator_mercurius_load;
    integrator->save        = reb_integrator_mercurius_save;
    integrator->particle_add= reb_integrator_mercurius_particle_add;
    integrator->particle_remove= reb_integrator_mercurius_particle_remove;
}

void reb_integator_mercurius_calculate_acceleration_whfast_part(struct reb_simulation* r){
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    struct reb_particle* const particles = r->particles;
    const int N = r->N;
    const int N_active = r->N_active;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const int _N_real   = N  - r->N_var;
    const int _N_active = ((N_active==-1)?_N_real:N_active);
    const int _testparticle_type   = r->testparticle_type;
    double (*_L) (const struct reb_simulation* const r, double d, double dcrit) = config->L;
    const double* const dcrit = config->dcrit;
#ifndef OPENMP
    for (int i=0; i<_N_real; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
    }
    for (int i=2; i<_N_active; i++){
        if (reb_sigint) return;
        for (int j=1; j<i; j++){
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[i],dcrit[j]);
            const double L = _L(r,_r,dcritmax);
            const double prefact = G*L/(_r*_r*_r);
            const double prefactj = -prefact*particles[j].m;
            const double prefacti = prefact*particles[i].m;
            particles[i].ax    += prefactj*dx;
            particles[i].ay    += prefactj*dy;
            particles[i].az    += prefactj*dz;
            particles[j].ax    += prefacti*dx;
            particles[j].ay    += prefacti*dy;
            particles[j].az    += prefacti*dz;
        }
    }
    const int startitestp = MAX(_N_active,2);
    for (int i=startitestp; i<_N_real; i++){
        if (reb_sigint) return;
        for (int j=1; j<_N_active; j++){
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[i],dcrit[j]);
            const double L = _L(r,_r,dcritmax);
            const double prefact = G*L/(_r*_r*_r);
            const double prefactj = -prefact*particles[j].m;
            particles[i].ax    += prefactj*dx;
            particles[i].ay    += prefactj*dy;
            particles[i].az    += prefactj*dz;
            if (_testparticle_type){
                const double prefacti = prefact*particles[i].m;
                particles[j].ax    += prefacti*dx;
                particles[j].ay    += prefacti*dy;
                particles[j].az    += prefacti*dz;
            }
        }
    }
#else // OPENMP
    particles[0].ax = 0; 
    particles[0].ay = 0; 
    particles[0].az = 0; 
#pragma omp parallel for schedule(guided)
    for (int i=1; i<_N_real; i++){
        particles[i].ax = 0; 
        particles[i].ay = 0; 
        particles[i].az = 0; 
        for (int j=1; j<_N_active; j++){
            if (i==j) continue;
            const double dx = particles[i].x - particles[j].x;
            const double dy = particles[i].y - particles[j].y;
            const double dz = particles[i].z - particles[j].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[i],dcrit[j]);
            const double L = _L(r,_r,dcritmax);
            const double prefact = -G*particles[j].m*L/(_r*_r*_r);
            particles[i].ax    += prefact*dx;
            particles[i].ay    += prefact*dy;
            particles[i].az    += prefact*dz;
        }
    }
    if (_testparticle_type){
        for (int i=1; i<_N_active; i++){
            for (int j=_N_active; j<_N_real; j++){
                const double dx = particles[i].x - particles[j].x;
                const double dy = particles[i].y - particles[j].y;
                const double dz = particles[i].z - particles[j].z;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                const double dcritmax = MAX(dcrit[i],dcrit[j]);
                const double L = _L(r,_r,dcritmax);
                const double prefact = -G*particles[j].m*L/(_r*_r*_r);
                particles[i].ax    += prefact*dx;
                particles[i].ay    += prefact*dy;
                particles[i].az    += prefact*dz;
            }
        }
    }
#endif // OPENMP
    
    if (r->additional_forces  && config->mode==0){
        // For Mercurius:
        // Additional forces are only calculated in the kick step, not during close encounter
        // shift pos and velocity so that external forces are calculated in inertial frame
        // Note: Copying avoids degrading floating point performance
        if(r->N>config->allocatedN_additionalforces){
            config->particles_backup_additionalforces = realloc(config->particles_backup_additionalforces, r->N*sizeof(struct reb_particle));
            config->allocatedN_additionalforces = r->N;
        }
        memcpy(config->particles_backup_additionalforces,r->particles,r->N*sizeof(struct reb_particle)); 
        reb_integrator_mercurius_dh_to_inertial(r);

        r->additional_forces(r);

        struct reb_particle* restrict const particles = r->particles;
        struct reb_particle* restrict const backup = config->particles_backup_additionalforces;
        for (int i=0;i<r->N;i++){
            particles[i].x = backup[i].x;
            particles[i].y = backup[i].y;
            particles[i].z = backup[i].z;
            particles[i].vx = backup[i].vx;
            particles[i].vy = backup[i].vy;
            particles[i].vz = backup[i].vz;
        }
    }
    printf("ERROR Need to call next funvtion from ias15 part");
    exit(EXIT_FAILURE);

}

void reb_integator_mercurius_calculate_acceleration_ias15_part(struct reb_simulation* r){
    // TODO: Need to call this from IAS15
    struct reb_integrator_mercurius_config* config = (struct reb_integrator_mercurius_config*)r->integrator_selected->config;
    struct reb_particle* const particles = r->particles;
    const double G = r->G;
    const double softening2 = r->softening*r->softening;
    const int _testparticle_type   = r->testparticle_type;
    double (*_L) (const struct reb_simulation* const r, double d, double dcrit) = config->L;
    const double m0 = r->particles[0].m;
    const double* const dcrit = config->dcrit;
    const int encounterN = config->encounterN;
    const int encounterNactive = config->encounterNactive;
    int* map = config->encounter_map;
#ifndef OPENMP
    particles[0].ax = 0; // map[0] is always 0 
    particles[0].ay = 0; 
    particles[0].az = 0; 
    // Acceleration due to star
    for (int i=1; i<encounterN; i++){
        int mi = map[i];
        const double x = particles[mi].x;
        const double y = particles[mi].y;
        const double z = particles[mi].z;
        const double _r = sqrt(x*x + y*y + z*z + softening2);
        double prefact = -G/(_r*_r*_r)*m0;
        particles[mi].ax    = prefact*x;
        particles[mi].ay    = prefact*y;
        particles[mi].az    = prefact*z;
    }
    // We're in a heliocentric coordinate system.
    // The star feels no acceleration
    // Interactions between active-active
    for (int i=2; i<encounterNactive; i++){
        int mi = map[i];
        for (int j=1; j<i; j++){
            int mj = map[j];
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
            const double L = _L(r,_r,dcritmax);
            double prefact = G*(1.-L)/(_r*_r*_r);
            double prefactj = -prefact*particles[mj].m;
            double prefacti = prefact*particles[mi].m;
            particles[mi].ax    += prefactj*dx;
            particles[mi].ay    += prefactj*dy;
            particles[mi].az    += prefactj*dz;
            particles[mj].ax    += prefacti*dx;
            particles[mj].ay    += prefacti*dy;
            particles[mj].az    += prefacti*dz;
        }
    }
    // Interactions between active-testparticle
    const int startitestp = MAX(encounterNactive,2);
    for (int i=startitestp; i<encounterN; i++){
        int mi = map[i];
        for (int j=1; j<encounterNactive; j++){
            int mj = map[j];
            const double dx = particles[mi].x - particles[mj].x;
            const double dy = particles[mi].y - particles[mj].y;
            const double dz = particles[mi].z - particles[mj].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
            const double L = _L(r,_r,dcritmax);
            double prefact = G*(1.-L)/(_r*_r*_r);
            double prefactj = -prefact*particles[mj].m;
            particles[mi].ax    += prefactj*dx;
            particles[mi].ay    += prefactj*dy;
            particles[mi].az    += prefactj*dz;
            if (_testparticle_type){
                double prefacti = prefact*particles[mi].m;
                particles[mj].ax    += prefacti*dx;
                particles[mj].ay    += prefacti*dy;
                particles[mj].az    += prefacti*dz;
            }
        }
    }
#else // OPENMP
    particles[0].ax = 0; // map[0] is always 0 
    particles[0].ay = 0; 
    particles[0].az = 0; 
    // We're in a heliocentric coordinate system.
    // The star feels no acceleration
#pragma omp parallel for schedule(guided)
    for (int i=1; i<encounterN; i++){
        int mi = map[i];
        particles[mi].ax = 0; 
        particles[mi].ay = 0; 
        particles[mi].az = 0; 
        // Acceleration due to star
        const double x = particles[mi].x;
        const double y = particles[mi].y;
        const double z = particles[mi].z;
        const double _r = sqrt(x*x + y*y + z*z + softening2);
        double prefact = -G/(_r*_r*_r)*m0;
        particles[mi].ax    += prefact*x;
        particles[mi].ay    += prefact*y;
        particles[mi].az    += prefact*z;
        for (int j=1; j<encounterNactive; j++){
            if (i==j) continue;
            int mj = map[j];
            const double dx = x - particles[mj].x;
            const double dy = y - particles[mj].y;
            const double dz = z - particles[mj].z;
            const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
            const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
            const double L = _L(r,_r,dcritmax);
            double prefact = -G*particles[mj].m*(1.-L)/(_r*_r*_r);
            particles[mi].ax    += prefact*dx;
            particles[mi].ay    += prefact*dy;
            particles[mi].az    += prefact*dz;
        }
    }
    if (_testparticle_type){
#pragma omp parallel for schedule(guided)
        for (int i=1; i<encounterNactive; i++){
            int mi = map[i];
            const double x = particles[mi].x;
            const double y = particles[mi].y;
            const double z = particles[mi].z;
            for (int j=encounterNactive; j<encounterN; j++){
                int mj = map[j];
                const double dx = x - particles[mj].x;
                const double dy = y - particles[mj].y;
                const double dz = z - particles[mj].z;
                const double _r = sqrt(dx*dx + dy*dy + dz*dz + softening2);
                const double dcritmax = MAX(dcrit[mi],dcrit[mj]);
                const double L = _L(r,_r,dcritmax);
                double prefact = -G*particles[mj].m*(1.-L)/(_r*_r*_r);
                particles[mi].ax    += prefact*dx;
                particles[mi].ay    += prefact*dy;
                particles[mi].az    += prefact*dz;
            }
        }
    }
#endif // OPENMP
}
