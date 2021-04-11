/**
 * @file    input.c
 * @brief   Parse command line options and read retart files.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section     LICENSE
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
#include <getopt.h>
#include <string.h>
#include "particle.h"
#include "rebound.h"
#include "collision.h"
#include "input.h"
#include "tree.h"
#include "integrator_sei.h"
#include "simulationarchive.h"

double reb_read_double(int argc, char** argv, const char* argument, double _default){
    char* value = reb_read_char(argc,argv,argument);
    if (value){
        return atof(value);
    }
    return _default;
}

int reb_read_int(int argc, char** argv, const char* argument, int _default){
    char* value = reb_read_char(argc,argv,argument);
    if (value){
        return atoi(value);
    }
    return _default;
}


char* reb_read_char(int argc, char** argv, const char* argument){
    opterr = 0;
    optind = 1;
    while (1) {
        struct option long_options[] = {
            {NULL, required_argument, 0, 'a'},
            {0,0,0,0}
        };

        long_options[0].name = argument;

        /* getopt_long stores the option index here.   */
        int option_index = 0;
        //              short options. format abc:d::
        int c = getopt_long (argc, argv, "", long_options, &option_index);

        /* Detect the end of the options.   */
        if (c == -1) break;

        switch (c){
            case 'a':
                return optarg;
                break;
            default:
                break;
        }
    }
    return NULL;
}

size_t reb_input_stream_fread(struct reb_input_stream* stream, void *restrict ptr, size_t size, size_t nitems){
    if (stream->mem_stream!=NULL){
        // read from memory
        memcpy(ptr,stream->mem_stream,size*nitems);
        stream->mem_stream = (char*)(stream->mem_stream)+ size*nitems;
        return size*nitems;
    }else if(stream->file_stream!=NULL){
        // read from file
        size_t s =  fread(ptr,size,nitems,stream->file_stream);
        return s;
    }
    return 0; 
}

static int reb_fseek(struct reb_input_stream* stream, long offset, int whence){
    if (stream->mem_stream!=NULL){
        // read from memory
        if (whence==SEEK_CUR){
            stream->mem_stream = (char*)(stream->mem_stream)+offset;
            return 0;
        }
        return -1;
    }else if(stream->file_stream!=NULL){
        // read from file
        return fseek(stream->file_stream,offset,whence);
    }
    return -1;
}

// Macro to read a single field from a binary file.
#define CASE(typename, value) case REB_BINARY_FIELD_TYPE_##typename: \
    {\
        reb_input_stream_fread(stream, value, field.size, 1);\
    }\
    break;

int reb_input_field(struct reb_simulation* r, struct reb_input_stream* stream, enum reb_input_binary_messages* warnings){
    struct reb_binary_field field;
    int numread = reb_input_stream_fread(stream, &field, sizeof(struct reb_binary_field), 1);
    if (numread<1){
        return 0; // End of file
    }
    switch (field.type){
        CASE(T,                  &r->t);
        CASE(G,                  &r->G);
        CASE(OMEGA,              &r->Omega);
        CASE(OMEGAZ,             &r->Omega_z);
        CASE(SOFTENING,          &r->softening);
        CASE(DT,                 &r->dt);
        CASE(DTLASTDONE,         &r->dt_last_done);
        CASE(N,                  &r->N);
        CASE(NVAR,               &r->N_var);
        CASE(VARCONFIGN,         &r->var_config_N);
        CASE(NACTIVE,            &r->N_active);
        CASE(TESTPARTICLETYPE,   &r->testparticle_type);
        CASE(TESTPARTICLEHIDEWARNINGS,   &r->testparticle_hidewarnings);
        CASE(HASHCTR,            &r->hash_ctr);
        CASE(OPENINGANGLE2,      &r->opening_angle2);
        CASE(STATUS,             &r->status);
        CASE(EXACTFINISHTIME,    &r->exact_finish_time);
        CASE(FORCEISVELOCITYDEP, &r->force_is_velocity_dependent);
        CASE(GRAVITYIGNORETERMS, &r->gravity_ignore_terms);
        CASE(OUTPUTTIMINGLAST,   &r->output_timing_last);
        CASE(SAVEMESSAGES,       &r->save_messages);
        CASE(EXITMAXDISTANCE,    &r->exit_max_distance);
        CASE(EXITMINDISTANCE,    &r->exit_min_distance);
        CASE(USLEEP,             &r->usleep);
        CASE(TRACKENERGYOFFSET,  &r->track_energy_offset);
        CASE(ENERGYOFFSET,       &r->energy_offset);
        CASE(BOXSIZE,            &r->boxsize);
        CASE(BOXSIZEMAX,         &r->boxsize_max);
        CASE(ROOTSIZE,           &r->root_size);
        CASE(ROOTN,              &r->root_n);
        CASE(ROOTNX,             &r->root_nx);
        CASE(ROOTNY,             &r->root_ny);
        CASE(ROOTNZ,             &r->root_nz);
        CASE(NGHOSTX,            &r->nghostx);
        CASE(NGHOSTY,            &r->nghosty);
        CASE(NGHOSTZ,            &r->nghostz);
        CASE(COLLISIONRESOLVEKEEPSORTED, &r->collision_resolve_keep_sorted);
        CASE(MINIMUMCOLLISIONVELOCITY, &r->minimum_collision_velocity);
        CASE(COLLISIONSPLOG,     &r->collisions_plog);
        CASE(MAXRADIUS,          &r->max_radius);
        CASE(COLLISIONSNLOG,     &r->collisions_Nlog);
        CASE(CALCULATEMEGNO,     &r->calculate_megno);
        CASE(MEGNOYS,            &r->megno_Ys);
        CASE(MEGNOYSS,           &r->megno_Yss);
        CASE(MEGNOCOVYT,         &r->megno_cov_Yt);
        CASE(MEGNOVART,          &r->megno_var_t);
        CASE(MEGNOMEANT,         &r->megno_mean_t);
        CASE(MEGNOMEANY,         &r->megno_mean_Y);
        CASE(MEGNON,             &r->megno_n);
        CASE(SAAUTOINTERVAL,     &r->simulationarchive_auto_interval);
        CASE(SAAUTOWALLTIME,     &r->simulationarchive_auto_walltime);
        CASE(SANEXT,             &r->simulationarchive_next);
        CASE(WALLTIME,           &r->walltime);
        CASE(COLLISION,          &r->collision);
        CASE(VISUALIZATION,      &r->visualization);
        CASE(BOUNDARY,           &r->boundary);
        CASE(GRAVITY,            &r->gravity);
        CASE(PYTHON_UNIT_L,      &r->python_unit_l);
        CASE(PYTHON_UNIT_M,      &r->python_unit_m);
        CASE(PYTHON_UNIT_T,      &r->python_unit_t);
        CASE(STEPSDONE,          &r->steps_done);
        CASE(SAAUTOSTEP,         &r->simulationarchive_auto_step);
        CASE(SANEXTSTEP,         &r->simulationarchive_next_step);
        CASE(RAND_SEED,          &r->rand_seed);
        case REB_BINARY_FIELD_TYPE_INTEGRATOR_SELECTED:
        {
            unsigned int id = 0;
            int found = 0;
            reb_input_stream_fread(stream, &id, field.size,1);
            for (int i=0; i<r->integrators_available_N; i++){
                if (r->integrators_available[i].id==id){
                    found = 1;
                    r->integrator_selected = &(r->integrators_available[i]);
                }
            }
            if (found==0){
                *warnings |= REB_INPUT_BINARY_WARNING_FIELD_UNKOWN;
            }
        }
        case REB_BINARY_FIELD_TYPE_PARTICLES:
            if(r->particles){
                free(r->particles);
            }
            r->allocatedN = (int)(field.size/sizeof(struct reb_particle));
            if (field.size){
                r->particles = malloc(field.size);
                reb_input_stream_fread(stream, r->particles, field.size,1);
            }
            if (r->allocatedN<r->N && warnings){
                *warnings |= REB_INPUT_BINARY_WARNING_PARTICLES;
            }
            for (int l=0;l<r->allocatedN;l++){
                r->particles[l].c = NULL;
                r->particles[l].ap = NULL;
                r->particles[l].sim = r;
            }
            if (r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE || r->collision==REB_COLLISION_LINETREE){
                for (int l=0;l<r->allocatedN;l++){
                    reb_tree_add_particle_to_tree(r, l);
                }
            }
            break;
        case REB_BINARY_FIELD_TYPE_VARCONFIG:
            if (r->var_config){
                free(r->var_config);
            }
            if (r->var_config_N>0){
                r->var_config = malloc(field.size);
                reb_input_stream_fread(stream, r->var_config, field.size,1);
                for (int l=0;l<r->var_config_N;l++){
                    r->var_config[l].sim = r;
                }
            }
            break;
        case REB_BINARY_FIELD_TYPE_END:
            return 0;
        case REB_BINARY_FIELD_TYPE_FUNCTIONPOINTERS:
            {
                int fpwarn;
                reb_input_stream_fread(stream, &fpwarn, field.size,1);
                if (fpwarn && warnings){
                    *warnings |= REB_INPUT_BINARY_WARNING_POINTERS;
                }
            }
            break;
        case REB_BINARY_FIELD_TYPE_HEADER:
            {
                long objects = 0;
                // Input header.
                const long bufsize = 64 - sizeof(struct reb_binary_field);
                char readbuf[bufsize], curvbuf[bufsize];
                const char* header = "REBOUND Binary File. Version: ";
                sprintf(curvbuf,"%s%s",header+sizeof(struct reb_binary_field), reb_version_str);
                
                objects += reb_input_stream_fread(stream, readbuf,sizeof(char),bufsize);
                // Note: following compares version, but ignores githash.
                if(strncmp(readbuf,curvbuf,bufsize)!=0){
                    *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
                }
            }
            break;
        default:
            {   
                // All others
                size_t found = 0;
                for (int i=0;i<r->integrators_available_N;i++){
                    struct reb_integrator* integrator = &(r->integrators_available[i]);
                    found |= integrator->load(integrator, r, stream, field);
                    if (found){
                        break;
                    }
                }
                // No match found. Unknown field.
                if (!found){
                    if (warnings){
                        *warnings |= REB_INPUT_BINARY_WARNING_FIELD_UNKOWN;
                    }
                    reb_fseek(stream, field.size, SEEK_CUR);
                }
            }
            break;
    }
    return 1;
} 

struct reb_simulation* reb_input_process_warnings(struct reb_simulation* r, enum reb_input_binary_messages warnings){
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        reb_error(r,"Cannot read binary file. Check filename and file contents.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_WARNING_VERSION){
        reb_warning(r,"Binary file was saved with a different version of REBOUND. Binary format might have changed.");
    }
    if (warnings & REB_INPUT_BINARY_WARNING_POINTERS){
        reb_warning(r,"You have to reset function pointers after creating a reb_simulation struct with a binary file.");
    }
    if (warnings & REB_INPUT_BINARY_WARNING_PARTICLES){
        reb_warning(r,"Binary file might be corrupted. Number of particles found does not match expected number.");
    }
    if (warnings & REB_INPUT_BINARY_ERROR_FILENOTOPEN){
        reb_error(r,"Error while reading binary file (file was not open).");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_ERROR_OUTOFRANGE){
        reb_error(r,"Index out of range.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_ERROR_SEEK){
        reb_error(r,"Error while trying to seek file.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_WARNING_FIELD_UNKOWN){
        reb_warning(r,"Unknown field found in binary file.");
    }
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        reb_error(r,"Cannot read binary file. Check filename and file contents.");
        if (r) free(r);
        return NULL;
    }
    if (warnings & REB_INPUT_BINARY_WARNING_CORRUPTFILE){
        reb_warning(r,"The binary file seems to be corrupted. An attempt has been made to recover parts of it. However, it might not be possible to append snapshots to the current file.");
    }
    return r;
}

struct reb_simulation* reb_create_simulation_from_binary(char* filename){
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    struct reb_simulation* r = reb_simulation_new();
    
    struct reb_simulationarchive* sa = malloc(sizeof(struct reb_simulationarchive)); 
    reb_read_simulationarchive_with_messages(sa, filename, NULL, &warnings);
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE){
        // Don't output an error if file does not exist, just return NULL.
        free(sa);
        return NULL;
    }else{
        reb_input_process_warnings(NULL, warnings);
    }
    reb_create_simulation_from_simulationarchive_with_messages(r, sa, -1, &warnings);
    reb_close_simulationarchive(sa);
    r = reb_input_process_warnings(r, warnings);
    return r;
}

