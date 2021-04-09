/**
 * @file 	integrator_sei.c
 * @brief 	Symplectic Epicycle Integrator (SEI).
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the Symplectic Epicycle Integrator 
 * (SEI). The integrator is described in detail in Rein & Tremaine 2011. 
 * It solves epicyclic motion exactly and is therefore exact up to machine
 * precision in the limit of no perturbing forces. When perturbing-forces
 * are of order eps, then the error of the scheme is O(eps dt^3). It also
 * makes use of two shear operators instead of a rotation to minimize 
 * systematic numerical round-off errors.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#include "rebound.h"
#include "particle.h"
#include "gravity.h"
#include "boundary.h"
#include "integrator.h"
#include "integrator_sei.h"
#include "output.h"
#include "input.h"


struct reb_integrator_sei_config {
    // Cached values
    double lastdt;      ///< Cached sin(), tan() for this value of dt.
    double sindt;       ///< Cached sin() 
    double tandt;       ///< Cached tan() 
    double sindtz;      ///< Cached sin(), z axis
    double tandtz;      ///< Cached tan(), z axis
};

enum SEI_CONFIG {
    LASTDT = 0,
    SINDT  = 1,
    TANDT  = 2,
    SINDTZ = 3,
    TANDTZ = 4,
};


size_t reb_integrator_sei_config_load(struct reb_integrator* integrator, struct reb_input_stream* stream, struct reb_binary_field field){
    struct reb_integrator_sei_config* config = (struct reb_integrator_sei_config*)integrator->config;
    switch (field.type){
        case REB_BF(SEI, LASTDT):
            return reb_input_stream_fread(stream, &config->lastdt, field.size, 1);
        case REB_BF(SEI, SINDT):
            return reb_input_stream_fread(stream, &config->sindt, field.size, 1);
        case REB_BF(SEI, TANDT):
            return reb_input_stream_fread(stream, &config->tandt, field.size, 1);
        case REB_BF(SEI, SINDTZ):
            return reb_input_stream_fread(stream, &config->sindtz, field.size, 1);
        case REB_BF(SEI, TANDTZ):
            return reb_input_stream_fread(stream, &config->tandtz, field.size, 1);
        default:
            return 0;
    }
}

void reb_integrator_sei_config_save(struct reb_integrator* integrator, struct reb_output_stream* stream){
    struct reb_integrator_sei_config* config = (struct reb_integrator_sei_config*)integrator->config;
    reb_output_stream_write_field(stream, REB_BF(SEI, LASTDT),   &(config->lastdt),   sizeof(double));
    reb_output_stream_write_field(stream, REB_BF(SEI, SINDT),    &(config->sindt),    sizeof(double));
    reb_output_stream_write_field(stream, REB_BF(SEI, TANDT),    &(config->tandt),    sizeof(double));
    reb_output_stream_write_field(stream, REB_BF(SEI, SINDTZ),   &(config->sindtz),   sizeof(double));
    reb_output_stream_write_field(stream, REB_BF(SEI, TANDTZ),   &(config->tandtz),   sizeof(double));
}

struct reb_integrator_sei_config* reb_integrator_sei_config_alloc(){
    struct reb_integrator_sei_config* config = calloc(1, sizeof(struct reb_integrator_sei_config)); // sets all variables to zero
    return config;
}
void reb_integrator_sei_config_free(struct reb_integrator_sei_config* config){
    free(config);
}

void reb_integrator_sei_register(struct reb_simulation* r){
    struct reb_integrator integrator = {
            .name = "sei",
               .id = 2,
               .config_load = reb_integrator_sei_config_load 
    };
    
    integrator.config = reb_integrator_sei_config_alloc();

    r->integrators_available = realloc(r->integrators_available, sizeof(struct reb_integrator) * r->integrators_available_N+1);
    r->integrators_available[r->integrators_available_N] = integrator;
    r->integrators_available_N++;

}

/**
 * @brief This function evolves a particle under the unperturbed
 * Hamiltonian H0 exactly up to machine precission.
 */
static void operator_H012(double dt, const double Omega, const double Omega_z, const struct reb_integrator_sei_config config, struct reb_particle* p){
		
	// Integrate vertical motion
	const double zx = p->z * Omega_z;
	const double zy = p->vz;
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double zt1 =  zx - config.tandtz*zy;			
	const double zyt =  config.sindtz*zt1 + zy;
	const double zxt =  zt1 - config.tandtz*zyt;	
	p->z  = zxt/Omega_z;
	p->vz = zyt;

	// Integrate motion in xy directions
	const double aO = 2.*p->vy + 4.*p->x*Omega;	// Center of epicyclic motion
	const double bO = p->y*Omega - 2.*p->vx;	

	const double ys = (p->y*Omega-bO)/2.; 		// Epicycle vector
	const double xs = (p->x*Omega-aO); 
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double xst1 =  xs - config.tandt*ys;			
	const double yst  =  config.sindt*xst1 + ys;
	const double xst  =  xst1 - config.tandt*yst;	

	p->x  = (xst+aO)    /Omega;			
	p->y  = (yst*2.+bO) /Omega - 3./4.*aO*dt;	
	p->vx = yst;
	p->vy = -xst*2. -3./2.*aO;
}

/**
 * @brief This function applies the acceleration due to the PHI1 term.
 * @details It is only exact if the forces are velocity independet (i.e. gravity).
 * If the forces are velocity dependent, it breaks the symmetry of the scheme,
 * making it firsr-order and non-symplectic. As long as these forces are small,
 * this should not be visible. However, it is worth keeping in mind. 
 * @param p reb_particle to evolve.
 * @param dt Timestep
 */
static void operator_phi1(double dt, struct reb_particle* p){
	// The force used here is for test cases 2 and 3 
	// in Rein & Tremaine 2011. 
	p->vx += p->ax * dt;
	p->vy += p->ay * dt;
	p->vz += p->az * dt;
}


void reb_integrator_sei_step(struct reb_simulation* const r, struct reb_integrator* integrator){
    struct reb_integrator_sei_config* config = (struct reb_integrator_sei_config*)integrator->config;
    r->gravity_ignore_terms = 0;
	const int N = r->N;
	struct reb_particle* const particles = r->particles;
        
    double const Omega = r->Omega;
    double Omega_z = r->Omega_z;
    if (isnan(Omega_z)){
        Omega_z = Omega;
    }
	
    if (config->lastdt!=r->dt){
        /**
         * Pre-calculates sin() and tan() needed for SEI. 
         */
        config->sindt = sin(Omega*(-r->dt/2.));
        config->tandt = tan(Omega*(-r->dt/4.));
        config->sindtz = sin(Omega_z*(-r->dt/2.));
        config->tandtz = tan(Omega_z*(-r->dt/4.));
        config->lastdt = r->dt;
	}

#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		operator_H012(r->dt, Omega, Omega_z, *config, &(particles[i]));
	}
	r->t+=r->dt/2.;

    reb_update_acceleration(r);

#pragma omp parallel for schedule(guided)
	for (int i=0;i<N;i++){
		operator_phi1(r->dt, &(particles[i]));
		operator_H012(r->dt, Omega, Omega_z, *config, &(particles[i]));
	}
	r->t+=r->dt/2.;
	r->dt_last_done = r->dt;
}

void reb_integrator_sei_synchronize(struct reb_simulation* r){
	// Do nothing.
}


