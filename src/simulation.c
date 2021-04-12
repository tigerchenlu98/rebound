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

#include "simulation.h"
#include "rebound.h"

#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "boundary.h"
#include "collision.h"
#include "input.h"
#include "integrator_eos.h"
#include "integrator_ias15.h"
#include "integrator_janus.h"
#include "integrator_leapfrog.h"
#include "integrator_mercurius.h"
#include "integrator_saba.h"
#include "integrator_sei.h"
#include "integrator_whfast.h"
#include "output.h"
#include "tools.h"
#include "tree.h"
#include "particle.h"
#include "display.h"
#include "simulationarchive.h"

/////////////////////////////////////////
// Memory management

struct reb_simulation* reb_simulation_new() {
    struct reb_simulation* r = malloc(sizeof(struct reb_simulation));
    reb_simulation_init(r); // includes memset(0)
    return r;
}

void reb_simulation_free(struct reb_simulation* const r) {
    reb_simulation_destroy(r);
    free(r);
}

void reb_simulation_init(struct reb_simulation* r) {
    // Set variables to their default value.
    // Only non-zero ones need to be set explicitly.
    memset(r, 0, sizeof(struct reb_simulation));
    reb_tools_init_srand(r);
    r->G       = 1;
    r->Omega   = 1;
    r->Omega_z = nan(""); // Will default to Omega
    r->dt      = 0.001;

    // Tree related parameters
    r->root_size      = -1; // Will be set by user. If not, this generates an error
    r->root_nx        = 1;
    r->root_ny        = 1;
    r->root_nz        = 1;
    r->root_n         = 1;
    r->opening_angle2 = 0.25;

    r->N_active           = -1;
    r->status             = REB_RUNNING;
    r->exact_finish_time  = 1;
    r->output_timing_last = -1;

    r->boundary  = REB_BOUNDARY_NONE;
    r->gravity   = REB_GRAVITY_BASIC;
    r->collision = REB_COLLISION_NONE;

    // Integrators
    reb_integrator_ias15_register(r);
    reb_integrator_leapfrog_register(r);
    reb_integrator_sei_register(r);
    reb_integrator_janus_register(r);
    reb_integrator_whfast_register(r);
    reb_integrator_mercurius_register(r);
    reb_integrator_eos_register(r);
    reb_integrator_saba_register(r);

    for (int i = 0; i < r->integrators_available_N; i++) {
        if (r->integrators_available[i].new) {
            r->integrators_available[i].config =
                r->integrators_available[i].new(&(r->integrators_available[i]),
                                                r);
        }
    }
    r->integrator_selected =
        r->integrators_available; // first integrator is the default (IAS15)

#ifdef OPENGL
    r->visualization = REB_VISUALIZATION_OPENGL;
#else  // OPENGL
    r->visualization = REB_VISUALIZATION_NONE;
#endif // OPENGL
#ifdef OPENMP
    printf("Using OpenMP with %d threads per node.\n", omp_get_max_threads());
#endif // OPENMP
}

void reb_simulation_destroy(struct reb_simulation* const r) {
    free(r->simulationarchive_filename);
    reb_tree_delete(r);
    if (r->display_data) {
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

    for (int i = 0; i < r->integrators_available_N; i++) {
        if (r->integrators_available[i].free) {
            r->integrators_available[i].free(&(r->integrators_available[i]), r);
        }
    }

    if (r->free_particle_ap) {
        for (int i = 0; i < r->N; i++) {
            r->free_particle_ap(&r->particles[i]);
        }
    }
    free(r->particles);
    free(r->particle_lookup_table);
    if (r->messages) {
        for (int i = 0; i < reb_max_messages_N; i++) {
            free(r->messages[i]);
        }
    }
    free(r->messages);
    if (r->extras_cleanup) {
        r->extras_cleanup(r);
    }
    free(r->var_config);
}

struct reb_simulation* reb_simulation_copy(struct reb_simulation* r) {
    struct reb_simulation* r_copy           = reb_simulation_new();
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;

    struct reb_output_stream stream;
    stream.buf       = NULL;
    stream.size      = 0;
    stream.allocated = 0;

    reb_output_stream_write_binary(&stream, r);

    r_copy->simulationarchive_filename = NULL; // Not sure about this one.

    struct reb_input_stream istream = {0};
    istream.mem_stream              = stream.buf;
    while (reb_input_field(r_copy, &istream, &warnings)) {
    }
    free(stream.buf);
    r = reb_input_process_warnings(r, warnings);
    return r_copy;
}

int reb_simulation_diff(struct reb_simulation* r1, struct reb_simulation* r2, int output_option) {
    if (output_option != 1 && output_option != 2) {
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
    int ret                          = reb_binary_diff(&istream1, &istream2, &ostream, output_option);

    free(stream1.buf);
    free(stream2.buf);
    return ret;
}

/////////////////////////////////////////
// Time stepping

void reb_simulation_steps(struct reb_simulation* const r, unsigned int N_steps) {
    for (unsigned int i = 0; i < N_steps; i++) {
        reb_simulation_step(r);
    }
}

void reb_simulation_step(struct reb_simulation* const r) {
    // Update walltime
    struct timeval time_beginning;
    gettimeofday(&time_beginning, NULL);

    if (r->pre_timestep_modifications) {
        if (r->integrator_selected->synchronize) {
            r->integrator_selected->synchronize(r->integrator_selected, r);
        }
        r->pre_timestep_modifications(r);
        // r->ri_whfast.recalculate_coordinates_this_timestep = 1; TODO
        // Reimplement this
        // r->ri_mercurius.recalculate_coordinates_this_timestep = 1;
    }

    r->integrator_selected->step(r->integrator_selected, r);

    if (r->post_timestep_modifications) {
        if (r->integrator_selected->synchronize) {
            r->integrator_selected->synchronize(r->integrator_selected, r);
        }
        r->post_timestep_modifications(r);
        // r->ri_whfast.recalculate_coordinates_this_timestep = 1; TODO
        // Reimplement this
        // r->ri_mercurius.recalculate_coordinates_this_timestep = 1;
    }

    // Do collisions here. We need both the positions and velocities at the same
    // time. Check for root crossings.
    reb_boundary_check(r);
    if (r->tree_needs_update) {
        // Update tree (this will remove particles which left the box)
        reb_tree_update(r);
    }

    // Search for collisions using local and essential tree.
    reb_collision_search(r);

    // Update walltime
    struct timeval time_end;
    gettimeofday(&time_end, NULL);
    r->walltime += time_end.tv_sec - time_beginning.tv_sec +
                   (time_end.tv_usec - time_beginning.tv_usec) / 1e6;
    // Update step counter
    r->steps_done++; // This also counts failed IAS15 steps
}


/////////////////////////////////////////
///  Integrate functions 

static int reb_error_message_waiting(struct reb_simulation* const r) {
    if (r->messages) {
        for (int i = 0; i < reb_max_messages_N; i++) {
            if (r->messages[i] != NULL) {
                if (r->messages[i][0] == 'e') {
                    return 1;
                }
            }
        }
    }
    return 0;
}
void reb_run_heartbeat(struct reb_simulation* const r) {
    if (r->heartbeat) {
        r->heartbeat(r);
    } // Heartbeat
    if (r->display_heartbeat) {
        reb_check_for_display_heartbeat(r);
    }
    if (r->exit_max_distance) {
        // Check for escaping particles
        const double max2                          = r->exit_max_distance * r->exit_max_distance;
        const struct reb_particle* const particles = r->particles;
        const int N                                = r->N - r->N_var;
        for (int i = 0; i < N; i++) {
            struct reb_particle p = particles[i];
            double r2             = p.x * p.x + p.y * p.y + p.z * p.z;
            if (r2 > max2) {
                r->status = REB_EXIT_ESCAPE;
            }
        }
    }
    if (r->exit_min_distance) {
        // Check for close encounters
        const double min2                          = r->exit_min_distance * r->exit_min_distance;
        const struct reb_particle* const particles = r->particles;
        const int N                                = r->N - r->N_var;
        for (int i = 0; i < N; i++) {
            struct reb_particle pi = particles[i];
            for (int j = 0; j < i; j++) {
                struct reb_particle pj = particles[j];
                const double x         = pi.x - pj.x;
                const double y         = pi.y - pj.y;
                const double z         = pi.z - pj.z;
                const double r2        = x * x + y * y + z * z;
                if (r2 < min2) {
                    r->status = REB_EXIT_ENCOUNTER;
                }
            }
        }
    }
    if (r->usleep > 0) {
        usleep(r->usleep);
    }
}

int reb_check_exit(struct reb_simulation* const r, const double tmax, double* last_full_dt) {
    while (r->status == REB_RUNNING_PAUSED) {
        // Wait for user to disable paused simulation
        usleep(1000);
    }
    const double dtsign = copysign(1., r->dt); // Used to determine integration direction
    if (reb_error_message_waiting(r)) {
        r->status = REB_EXIT_ERROR;
    }
    if (r->status >= 0) {
        // Exit now.
    } else if (tmax != INFINITY) {
        if (r->exact_finish_time == 1) {
            if ((r->t + r->dt) * dtsign >= tmax * dtsign) { // Next step would overshoot
                if (r->t == tmax) {
                    r->status = REB_EXIT_SUCCESS;
                } else if (r->status == REB_RUNNING_LAST_STEP) {
                    double tscale = 1e-12 * fabs(tmax); // Find order of magnitude for time
                    if (tscale < 1e-200) {              // Failsafe if tmax==0.
                        tscale = 1e-12;
                    }
                    if (fabs(r->t - tmax) < tscale) {
                        r->status = REB_EXIT_SUCCESS;
                    } else {
                        // not there yet, do another step.
                        if (r->integrator_selected->synchronize) {
                            r->integrator_selected->synchronize(r->integrator_selected, r);
                        }
                        r->dt = tmax - r->t;
                    }
                } else {
                    r->status = REB_RUNNING_LAST_STEP; // Do one small step, then exit.
                    if (r->integrator_selected->synchronize) {
                        r->integrator_selected->synchronize(r->integrator_selected, r);
                    }
                    if (r->dt_last_done != 0.) {         // If first timestep is also last, do not use dt_last_done (which would be 0.)
                        *last_full_dt = r->dt_last_done; // store last full dt before decreasing the timestep to match finish time
                    }
                    r->dt = tmax - r->t;
                }
            } else {
                if (r->status == REB_RUNNING_LAST_STEP) {
                    // This will get executed if an adaptive integrator reduces
                    // the timestep in what was supposed to be the last timestep.
                    r->status = REB_RUNNING;
                }
            }
        } else {
            if (r->t * dtsign >= tmax * dtsign) { // Past the integration time
                r->status = REB_EXIT_SUCCESS;     // Exit now.
            }
        }
    }
    if (r->N <= 0) {
        reb_warning(r, "No particles found. Will exit.");
        r->status = REB_EXIT_NOPARTICLES; // Exit now.
    }
    return r->status;
}

struct reb_thread_info {
    struct reb_simulation* r;
    double tmax;
};

volatile sig_atomic_t reb_sigint;

void reb_sigint_handler(int signum) {
    // Handles graceful shutdown for interrupts
    if (signum == SIGINT) {
        reb_sigint = 1;
    }
}

static void* reb_integrate_raw(void* args) {
    reb_sigint = 0;
    signal(SIGINT, reb_sigint_handler);
    struct reb_thread_info* thread_info = (struct reb_thread_info*)args;
    struct reb_simulation* const r      = thread_info->r;

    double last_full_dt = r->dt; // need to store r->dt in case timestep gets artificially shrunk to meet exact_finish_time=1
    r->dt_last_done     = 0.;    // Reset in case first timestep attempt will fail

    if (r->testparticle_hidewarnings == 0 && reb_particle_check_testparticles(r)) {
        reb_warning(r, "At least one test particle (type 0) has finite mass. This might lead to unexpected behaviour. Set testparticle_hidewarnings=1 to hide this warning.");
    }

    r->status = REB_RUNNING;
    reb_run_heartbeat(r);
    while (reb_check_exit(r, thread_info->tmax, &last_full_dt) < 0) {
#ifdef OPENGL
        if (r->display_data) {
            if (r->display_data->opengl_enabled) {
                pthread_mutex_lock(&(r->display_data->mutex));
            }
        }
#endif // OPENGL
        if (r->simulationarchive_filename) {
            reb_simulationarchive_heartbeat(r);
        }
        reb_simulation_step(r);
        reb_run_heartbeat(r);
        if (reb_sigint == 1) {
            r->status = REB_EXIT_SIGINT;
        }
#ifdef OPENGL
        if (r->display_data) {
            if (r->display_data->opengl_enabled) {
                pthread_mutex_unlock(&(r->display_data->mutex));
            }
        }
#endif // OPENGL
    }

    if (r->integrator_selected->synchronize) {
        r->integrator_selected->synchronize(r->integrator_selected, r);
    }
    if (r->display_heartbeat) { // Display Heartbeat
        r->display_heartbeat(r);
    }
    if (r->exact_finish_time == 1) { // if finish_time = 1, r->dt could have been shrunk, so set to the last full timestep
        r->dt = last_full_dt;
    }
    if (r->simulationarchive_filename) {
        reb_simulationarchive_heartbeat(r);
    }

    return NULL;
}

enum REB_STATUS reb_integrate(struct reb_simulation* const r, double tmax) {
    struct reb_thread_info thread_info = {
        .r    = r,
        .tmax = tmax,
    };
    switch (r->visualization) {
    case REB_VISUALIZATION_NONE: {
        if (r->display_data) {
            r->display_data->opengl_enabled = 0;
        }
        reb_integrate_raw(&thread_info);
    } break;
    case REB_VISUALIZATION_OPENGL: {
#ifdef OPENGL
        reb_display_init_data(r);
        r->display_data->opengl_enabled = 1;

        pthread_t compute_thread;
        if (pthread_create(&compute_thread, NULL, reb_integrate_raw, &thread_info)) {
            reb_error(r, "Error creating display thread.");
        }

        reb_display_init(r); // Display routines running on main thread.

        if (pthread_join(compute_thread, NULL)) {
            reb_error(r, "Error joining display thread.");
        }
#else  // OPENGL
        reb_error(r, "REBOUND was not compiled/linked with OPENGL libraries.");
        return REB_EXIT_ERROR;
#endif // OPENGL
    } break;
    case REB_VISUALIZATION_WEBGL: {
        reb_display_init_data(r);
        reb_integrate_raw(&thread_info);
    } break;
    }
    return r->status;
}
