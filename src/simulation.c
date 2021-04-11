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
#include <math.h>
#include "rebound.h"
#include "simulation.h"
#include "tools.h"
#include "output.h"
#include "input.h"
#include "tree.h"

#include "integrator_saba.h"
#include "integrator_whfast.h"
#include "integrator_ias15.h"
#include "integrator_sei.h"
#include "integrator_janus.h"
#include "integrator_eos.h"
#include "integrator_leapfrog.h"
#include "integrator_mercurius.h"
   
void reb_simulation_init(struct reb_simulation* r){
    // Set variables to their default value.
    // Only non-zero ones need to be set explicitly. 
    memset(r,0,sizeof(struct reb_simulation));
    reb_tools_init_srand(r);
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
        if (r->integrators_available[i].new){
            r->integrators_available[i].config = r->integrators_available[i].new(&(r->integrators_available[i]), r);
        }
    }
    r->integrator_selected = r->integrators_available; // first integrator is the default (IAS15)

#ifdef OPENGL
    r->visualization= REB_VISUALIZATION_OPENGL;
#else // OPENGL
    r->visualization= REB_VISUALIZATION_NONE;
#endif // OPENGL
#ifdef OPENMP
    printf("Using OpenMP with %d threads per node.\n",omp_get_max_threads());
#endif // OPENMP
}

struct reb_simulation* reb_simulation_new(){
    struct reb_simulation* r = malloc(sizeof(struct reb_simulation));
    reb_simulation_init(r); // includes memset(0)
    return r;
}

void reb_simulation_destroy(struct reb_simulation* const r){
    free(r->simulationarchive_filename);
    reb_tree_delete(r);
    if(r->display_data){
        pthread_mutex_destroy(&(r->display_data->mutex));
        free(r->display_data->r_copy);
        free(r->display_data->particles_copy);
        free(r->display_data->p_jh_copy);
        free(r->display_data->particle_data);
        free(r->display_data->orbit_data);
        free(r->display_data); // TODO: Free other pointers in display_data
    }
    free(r->gravity_cs);
    free(r->collisions);
    
    for (int i=0; i<r->integrators_available_N; i++){
        if (r->integrators_available[i].free){
            r->integrators_available[i].free(&(r->integrators_available[i]), r);
        }
    }

    if(r->free_particle_ap){
        for(int i=0; i<r->N; i++){
            r->free_particle_ap(&r->particles[i]);
        }
    }
    free(r->particles);
    free(r->particle_lookup_table);
    if (r->messages){
        for (int i=0;i<reb_max_messages_N;i++){
            free(r->messages[i]);
        }
    }
    free(r->messages);
    if (r->extras_cleanup){
        r->extras_cleanup(r);
    }
    free(r->var_config);
}

void reb_simulation_free(struct reb_simulation* const r){
    reb_simulation_destroy(r);
    free(r);
}


struct reb_simulation* reb_simulation_copy(struct reb_simulation* r){
    struct reb_simulation* r_copy = reb_simulation_new();
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    
    struct reb_output_stream stream;
    stream.buf = NULL;
    stream.size = 0;
    stream.allocated = 0;
    
    reb_output_stream_write_binary(&stream, r);
    
    r_copy->simulationarchive_filename = NULL; // Not sure about this one.
    
    struct reb_input_stream istream = {0};
    istream.mem_stream = stream.buf;
    while(reb_input_field(r_copy, &istream, &warnings)){ }
    free(stream.buf);
    r = reb_input_process_warnings(r, warnings);
    return r_copy;
}
