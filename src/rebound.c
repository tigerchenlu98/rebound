/**
 * @file    rebound.c
 * @brief   Main REBOUND control structures and routine, iteration loop.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section LICENSE
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
#include <sys/types.h>
#include <string.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <pthread.h>
#include <fcntl.h>
#include "rebound.h"
#include "boundary.h"
#include "gravity.h"
#include "collision.h"
#include "tree.h"
#include "output.h"
#include "tools.h"
#include "particle.h"
#include "input.h"
#include "binarydiff.h"
#include "simulationarchive.h"
#include "display.h"
#ifdef OPENMP
#include <omp.h>
#endif
#define MAX(a, b) ((a) < (b) ? (b) : (a))       ///< Returns the maximum of a and b
#define STRINGIFY(s) str(s)
#define str(s) #s

const int reb_max_messages_length = 1024;   // needs to be constant expression for array size
const int reb_max_messages_N = 10;
const char* reb_build_str = __DATE__ " " __TIME__;  // Date and time build string. 
const char* reb_version_str = "3.16.0";         // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* reb_githash_str = STRINGIFY(GITHASH);             // This line gets updated automatically. Do not edit manually.
    

static int reb_error_message_waiting(struct reb_simulation* const r);

void reb_simulation_set_integrator(struct reb_simulation* r, const char* name){
    for (int i=0; i<r->integrators_available_N; i++){
        if (strcmp(r->integrators_available[i].name,name)==0){
            r->integrator_selected = &(r->integrators_available[i]);
            return;
        }
    }
    printf("Error: Integrator not found."); // TODO
    exit(EXIT_FAILURE);
}


void reb_update_acceleration(struct reb_simulation* r){
    // Update and simplify tree. 
    // This function also creates the tree if called for the first time.
    if (r->tree_needs_update || r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE || r->collision==REB_COLLISION_LINETREE){
        // Check for root crossings.
        reb_boundary_check(r);     

        // Update tree (this will remove particles which left the box)
        reb_tree_update(r);          
    }
    if (r->tree_root!=NULL && r->gravity==REB_GRAVITY_TREE){
        // Update center of mass and quadrupole moments in tree in preparation of force calculation.
        reb_tree_update_gravity_data(r); 
    }

    // Main force calculation:
	reb_calculate_acceleration(r);
	if (r->N_var){
		reb_calculate_acceleration_var(r);
	}

	if (r->additional_forces){
        r->additional_forces(r);
    }
}

void* reb_simulation_get_integrator_config(struct reb_simulation* r, const char* name){
    for (int i=0; i<r->integrators_available_N; i++){
        if (strcmp(r->integrators_available[i].name,name)==0){
            return r->integrator_selected->config;
        }
    }
    printf("Error: Integrator not found."); // TODO
    exit(EXIT_FAILURE);
}

struct reb_integrator* reb_simulation_register_integrator(struct reb_simulation* r, const char* name, int id){
    for (int i=0; i<r->integrators_available_N; i++){
        if (r->integrators_available[i].id == id){
            printf("Error: Integrator ID already registered.");
            exit(EXIT_FAILURE);
        }
        if (strcmp(r->integrators_available[i].name,name)==0){
            printf("Error: Integrator name already registered.");
            exit(EXIT_FAILURE);
        }
    }

    struct reb_integrator integrator = {0};
    integrator.name = name;
    integrator.id   = id;

    r->integrators_available_N++;
    r->integrators_available = realloc(r->integrators_available, sizeof(struct reb_integrator) * r->integrators_available_N);
    r->integrators_available[r->integrators_available_N-1] = integrator;
    return &(r->integrators_available[r->integrators_available_N-1]);
}

void reb_steps(struct reb_simulation* const r, unsigned int N_steps){
    for (unsigned int i=0;i<N_steps;i++){
        reb_step(r);
    }
}
void reb_step(struct reb_simulation* const r){
    // Update walltime
    struct timeval time_beginning;
    gettimeofday(&time_beginning,NULL);

    if (r->pre_timestep_modifications){
        if (r->integrator_selected->synchronize){
            r->integrator_selected->synchronize(r->integrator_selected, r);
        }
        r->pre_timestep_modifications(r);
        //r->ri_whfast.recalculate_coordinates_this_timestep = 1; TODO Reimplement this
        //r->ri_mercurius.recalculate_coordinates_this_timestep = 1;
    }
   
    r->integrator_selected->step(r->integrator_selected, r);
    
    if (r->post_timestep_modifications){
        if (r->integrator_selected->synchronize){
            r->integrator_selected->synchronize(r->integrator_selected, r);
        }
        r->post_timestep_modifications(r);
        //r->ri_whfast.recalculate_coordinates_this_timestep = 1; TODO Reimplement this
        //r->ri_mercurius.recalculate_coordinates_this_timestep = 1;
    }

    // Do collisions here. We need both the positions and velocities at the same time.
    // Check for root crossings.
    reb_boundary_check(r);     
    if (r->tree_needs_update){
        // Update tree (this will remove particles which left the box)
        reb_tree_update(r);          
    }

    // Search for collisions using local and essential tree.
    reb_collision_search(r);
    
    // Update walltime
    struct timeval time_end;
    gettimeofday(&time_end,NULL);
    r->walltime += time_end.tv_sec-time_beginning.tv_sec+(time_end.tv_usec-time_beginning.tv_usec)/1e6;
    // Update step counter
    r->steps_done++; // This also counts failed IAS15 steps
}

void reb_exit(const char* const msg){
    // This function should also kill all children. 
    // Not implemented as pid is not easy to get to.
    // kill(pid, SIGKILL);
    fprintf(stderr,"\n\033[1mFatal error! Exiting now.\033[0m %s\n",msg);
    exit(EXIT_FAILURE);
}

void reb_message(struct reb_simulation* const r, char type, const char* const msg){
    int save_messages = 0;
    if (r != NULL){
        save_messages = r->save_messages;
    }
    if (!save_messages || strlen(msg)>=reb_max_messages_length){
        if (type=='w'){
            fprintf(stderr,"\n\033[1mWarning!\033[0m %s\n",msg);
        }else if (type=='e'){
            fprintf(stderr,"\n\033[1mError!\033[0m %s\n",msg);
        }
    }else{
        if (r->messages==NULL){
            r->messages = calloc(reb_max_messages_N,sizeof(char*));
        }
        int n = 0;
        for (;n<reb_max_messages_N;n++){
            if (r->messages[n]==NULL){
                break;
            }
        }
        if (n==reb_max_messages_N){
            free(r->messages[0]);
            for (int i=0;i<reb_max_messages_N-1;i++){
                r->messages[i] = r->messages[i+1];
            }
            r->messages[reb_max_messages_N-1] = NULL;
            n= reb_max_messages_N-1;
        }
        r->messages[n] = malloc(sizeof(char*)*reb_max_messages_length);
        r->messages[n][0] = type;
        strcpy(r->messages[n]+1, msg);
    }
}

void reb_warning(struct reb_simulation* const r, const char* const msg){
    reb_message(r, 'w', msg);
}

void reb_error(struct reb_simulation* const r, const char* const msg){
    reb_message(r, 'e', msg);
}

int reb_get_next_message(struct reb_simulation* const r, char* const buf){
    if (r->messages){
        char* w0 = r->messages[0];
        if (w0){
            for(int i=0;i<reb_max_messages_N-1;i++){
                r->messages[i] = r->messages[i+1];
            }
            r->messages[reb_max_messages_N-1] = NULL;
            strcpy(buf,w0);
            free(w0);
            return 1;
        }
    }
    return 0;
}

static int reb_error_message_waiting(struct reb_simulation* const r){
    if (r->messages){
        for (int i=0;i<reb_max_messages_N;i++){
            if (r->messages[i]!=NULL){
                if (r->messages[i][0]=='e'){
                    return 1;
                }
            }
        }
    }
    return 0;
}


void reb_configure_box(struct reb_simulation* const r, const double root_size, const int root_nx, const int root_ny, const int root_nz){
    r->root_size = root_size;
    r->root_nx = root_nx;
    r->root_ny = root_ny;
    r->root_nz = root_nz;
    // Setup box sizes
    r->boxsize.x = r->root_size *(double)r->root_nx;
    r->boxsize.y = r->root_size *(double)r->root_ny;
    r->boxsize.z = r->root_size *(double)r->root_nz;
    r->root_n = r->root_nx*r->root_ny*r->root_nz;
    r->boxsize_max = MAX(r->boxsize.x, MAX(r->boxsize.y, r->boxsize.z));
    if (r->root_nx <=0 || r->root_ny <=0 || r->root_nz <= 0){
        reb_exit("Number of root boxes must be greater or equal to 1 in each direction.");
    }
}


int reb_diff_simulations(struct reb_simulation* r1, struct reb_simulation* r2, int output_option){
    if (output_option!=1 && output_option!=2){
        // Not implemented
        return -1;
    }
    struct reb_output_stream stream1 = {0};
    reb_output_stream_write_binary(&stream1, r1);
    
    struct reb_output_stream stream2 = {0};
    reb_output_stream_write_binary(&stream2, r2);

    struct reb_input_stream istream1 = {.mem_stream = stream1.buf, .size = stream1.size, .file_stream = NULL};
    struct reb_input_stream istream2 = {.mem_stream = stream2.buf, .size = stream2.size, .file_stream = NULL};
    struct reb_output_stream ostream = {0};
    int ret = reb_binary_diff(&istream1, &istream2, &ostream, output_option);
    
    free(stream1.buf);
    free(stream2.buf);
    return ret;
}

    
int reb_check_exit(struct reb_simulation* const r, const double tmax, double* last_full_dt){
    while(r->status == REB_RUNNING_PAUSED){
        // Wait for user to disable paused simulation
        usleep(1000);
    }
    const double dtsign = copysign(1.,r->dt);   // Used to determine integration direction
    if (reb_error_message_waiting(r)){
        r->status = REB_EXIT_ERROR;
    }
    if (r->status>=0){
        // Exit now.
    }else if(tmax!=INFINITY){
        if(r->exact_finish_time==1){
            if ((r->t+r->dt)*dtsign>=tmax*dtsign){  // Next step would overshoot
                if (r->t==tmax){
                    r->status = REB_EXIT_SUCCESS;
                }else if(r->status == REB_RUNNING_LAST_STEP){
                    double tscale = 1e-12*fabs(tmax);   // Find order of magnitude for time
                    if (tscale<1e-200){     // Failsafe if tmax==0.
                        tscale = 1e-12;
                    }
                    if (fabs(r->t-tmax)<tscale){
                        r->status = REB_EXIT_SUCCESS;
                    }else{
                        // not there yet, do another step.
                        if (r->integrator_selected->synchronize){
                            r->integrator_selected->synchronize(r->integrator_selected, r);
                        }
                        r->dt = tmax-r->t;
                    }
                }else{
                    r->status = REB_RUNNING_LAST_STEP; // Do one small step, then exit.
                    if (r->integrator_selected->synchronize){
                        r->integrator_selected->synchronize(r->integrator_selected, r);
                    }
                    if (r->dt_last_done!=0.){   // If first timestep is also last, do not use dt_last_done (which would be 0.)
                        *last_full_dt = r->dt_last_done; // store last full dt before decreasing the timestep to match finish time
                    }
                    r->dt = tmax-r->t;
                }
            }else{
                if (r->status == REB_RUNNING_LAST_STEP){
                    // This will get executed if an adaptive integrator reduces
                    // the timestep in what was supposed to be the last timestep.
                    r->status = REB_RUNNING;
                }
            }
        }else{
            if (r->t*dtsign>=tmax*dtsign){  // Past the integration time
                r->status = REB_EXIT_SUCCESS; // Exit now.
            }
        }
    }
    if (r->N<=0){
        reb_warning(r,"No particles found. Will exit.");
        r->status = REB_EXIT_NOPARTICLES; // Exit now.
    }
    return r->status;
}


void reb_run_heartbeat(struct reb_simulation* const r){
    if (r->heartbeat){ r->heartbeat(r); }               // Heartbeat
    if (r->display_heartbeat){ reb_check_for_display_heartbeat(r); } 
    if (r->exit_max_distance){
        // Check for escaping particles
        const double max2 = r->exit_max_distance * r->exit_max_distance;
        const struct reb_particle* const particles = r->particles;
        const int N = r->N - r->N_var;
        for (int i=0;i<N;i++){
            struct reb_particle p = particles[i];
            double r2 = p.x*p.x + p.y*p.y + p.z*p.z;
            if (r2>max2){
                r->status = REB_EXIT_ESCAPE;
            }
        }
    }
    if (r->exit_min_distance){
        // Check for close encounters
        const double min2 = r->exit_min_distance * r->exit_min_distance;
        const struct reb_particle* const particles = r->particles;
        const int N = r->N - r->N_var;
        for (int i=0;i<N;i++){
            struct reb_particle pi = particles[i];
            for (int j=0;j<i;j++){
                struct reb_particle pj = particles[j];
                const double x = pi.x-pj.x;
                const double y = pi.y-pj.y;
                const double z = pi.z-pj.z;
                const double r2 = x*x + y*y + z*z;
                if (r2<min2){
                    r->status = REB_EXIT_ENCOUNTER;
                }
            }
        }
    }
    if (r->usleep > 0){
        usleep(r->usleep);
    }
}

////////////////////////////////////////////////////
///  Integrate functions and visualization stuff

struct reb_thread_info {
    struct reb_simulation* r;
    double tmax;
};

volatile sig_atomic_t reb_sigint;

void reb_sigint_handler(int signum) {
    // Handles graceful shutdown for interrupts
    if (signum == SIGINT){
        reb_sigint = 1;
    }
}

static void* reb_integrate_raw(void* args){
    reb_sigint = 0;
    signal(SIGINT, reb_sigint_handler);
    struct reb_thread_info* thread_info = (struct reb_thread_info*)args;
	struct reb_simulation* const r = thread_info->r;

    double last_full_dt = r->dt; // need to store r->dt in case timestep gets artificially shrunk to meet exact_finish_time=1
    r->dt_last_done = 0.; // Reset in case first timestep attempt will fail

    if (r->testparticle_hidewarnings==0 && reb_particle_check_testparticles(r)){
        reb_warning(r,"At least one test particle (type 0) has finite mass. This might lead to unexpected behaviour. Set testparticle_hidewarnings=1 to hide this warning.");
    }

    r->status = REB_RUNNING;
    reb_run_heartbeat(r);
    while(reb_check_exit(r,thread_info->tmax,&last_full_dt)<0){
#ifdef OPENGL
        if (r->display_data){
            if (r->display_data->opengl_enabled){ pthread_mutex_lock(&(r->display_data->mutex)); }
        }
#endif // OPENGL
        if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}
        reb_step(r); 
        reb_run_heartbeat(r);
        if (reb_sigint== 1){
            r->status = REB_EXIT_SIGINT;
        }
#ifdef OPENGL
        if (r->display_data){
            if (r->display_data->opengl_enabled){ pthread_mutex_unlock(&(r->display_data->mutex)); }
        }
#endif // OPENGL
    }

    if (r->integrator_selected->synchronize){
        r->integrator_selected->synchronize(r->integrator_selected, r);
    }
    if (r->display_heartbeat){                          // Display Heartbeat
        r->display_heartbeat(r); 
    }
    if(r->exact_finish_time==1){ // if finish_time = 1, r->dt could have been shrunk, so set to the last full timestep
        r->dt = last_full_dt; 
    }
    if (r->simulationarchive_filename){ reb_simulationarchive_heartbeat(r);}

    return NULL;
}

enum REB_STATUS reb_integrate(struct reb_simulation* const r, double tmax){
    struct reb_thread_info thread_info = {
        .r = r,
        .tmax = tmax, 
    };
    switch (r->visualization){
        case REB_VISUALIZATION_NONE:
            {
                if (r->display_data){
                    r->display_data->opengl_enabled = 0;
                }
                reb_integrate_raw(&thread_info);
            }
            break;
        case REB_VISUALIZATION_OPENGL:
            {
#ifdef OPENGL
                reb_display_init_data(r);
                r->display_data->opengl_enabled = 1;

                pthread_t compute_thread;
                if (pthread_create(&compute_thread,NULL,reb_integrate_raw,&thread_info)){
                    reb_error(r, "Error creating display thread.");
                }
                
                reb_display_init(r); // Display routines running on main thread.

                if (pthread_join(compute_thread,NULL)){
                    reb_error(r, "Error joining display thread.");
                }
#else // OPENGL
                reb_error(r,"REBOUND was not compiled/linked with OPENGL libraries.");
                return REB_EXIT_ERROR; 
#endif // OPENGL
            }
            break;
        case REB_VISUALIZATION_WEBGL:
            {
                reb_display_init_data(r);
                reb_integrate_raw(&thread_info);
            }
            break;
    }
    return r->status;
}

  
#ifdef OPENMP
void reb_omp_set_num_threads(int num_threads){
    omp_set_num_threads(num_threads);
}
#endif // OPENMP

const char* reb_logo[26] = {
"          _                           _  ",
"         | |                         | | ",
" _ __ ___| |__   ___  _   _ _ __   __| | ",
"| '__/ _ \\ '_ \\ / _ \\| | | | '_ \\ / _` | ",
"| | |  __/ |_) | (_) | |_| | | | | (_| | ",
"|_|  \\___|_.__/ \\___/ \\__,_|_| |_|\\__,_| ",
"                                         ",
"              `-:://::.`                 ",
"          `/oshhoo+++oossso+:`           ",
"       `/ssooys++++++ossssssyyo:`        ",
"     `+do++oho+++osssso++++++++sy/`      ",
"    :yoh+++ho++oys+++++++++++++++ss.     ",
"   /y++hooyyooshooo+++++++++++++++oh-    ",
"  -dsssdssdsssdssssssssssooo+++++++oh`   ",
"  ho++ys+oy+++ho++++++++oosssssooo++so   ",
" .d++oy++ys+++oh+++++++++++++++oosssod   ",
" -h+oh+++yo++++oyo+++++++++++++++++oom   ",
" `d+ho+++ys+++++oys++++++++++++++++++d   ",
"  yys++++oy+++++++oys+++++++++++++++s+   ",
"  .m++++++h+++++++++oys++++++++++++oy`   ",
"   -yo++++ss++++++++++oyso++++++++oy.    ",
"    .ss++++ho+++++++++++osys+++++yo`     ",
"      :ss+++ho+++++++++++++osssss-       ",
"        -ossoys++++++++++++osso.         ",
"          `-/oyyyssosssyso+/.            ",
"                ``....`                  ",
};
