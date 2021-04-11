/**
 * @file 	simulation.c
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein
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

#include <stdlib.h>
#include <string.h>
#include "rebound.h"
#include "simulation.h"
   
   
struct reb_simulation* reb_simulation_init(){
    struct reb_simulation* r = calloc(1,sizeof(struct reb_simulation));
    reb_tools_init_srand(r);
    // Set variables to their default value.
    // Only non-zero ones need to be set explicitly. 
    r->G        = 1;
    r->Omega    = 1;
    r->Omega_z  = nan(""); // Will default to Omega
    r->dt       = 0.001;
    r->root_size    = -1;  // Will be set by user. If not, this generates an error
    r->root_nx  = 1;
    r->root_ny  = 1;
    r->root_nz  = 1;
    r->root_n   = 1;
    r->opening_angle2   = 0.25;
    
    r->N_active     = -1;   
    r->status       = REB_RUNNING;
    r->exact_finish_time    = 1;
    r->output_timing_last   = -1;
    
#ifdef OPENGL
    r->visualization= REB_VISUALIZATION_OPENGL;
#else // OPENGL
    r->visualization= REB_VISUALIZATION_NONE;
#endif // OPENGL
    r->boundary     = REB_BOUNDARY_NONE;
    r->gravity      = REB_GRAVITY_BASIC;
    r->collision    = REB_COLLISION_NONE;

    reb_integrator_ias15_register(r);
    reb_integrator_leapfrog_register(r);
    reb_integrator_sei_register(r);
    reb_integrator_janus_register(r);
    reb_integrator_whfast_register(r);
    reb_integrator_mercurius_register(r);
    reb_integrator_eos_register(r);
    reb_integrator_saba_register(r);

    for (int i=0; i<r->integrators_available_N; i++){
        if (r->integrators_available[i].init){
            r->integrators_available[i].config = r->integrators_available[i].init(&(r->integrators_available[i]), r);
        }
    }

    r->integrator_selected = r->integrators_available; // first integrator is the default (IAS15)


#ifndef LIBREBOUND
    printf("Process id: %d.\n", getpid());
#endif // LIBREBOUND
#ifdef OPENMP
    printf("Using OpenMP with %d threads per node.\n",omp_get_max_threads());
#endif // OPENMP
}

