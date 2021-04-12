/**
 * @file    output.c
 * @brief   Output routines.
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "rebound.h"

#include "input.h"
#include "output.h"
#include "particle.h"
#include "tools.h"

void reb_output_stream_write(struct reb_output_stream* stream, void* restrict data, size_t size) {
    // Increase size
    int increased = 0;
    while (stream->allocated == 0 || stream->size + size > stream->allocated) {
        increased         = 1;
        stream->allocated = stream->allocated ? stream->allocated * 2 : 32;
    }
    if (increased) {
        stream->buf = realloc(stream->buf, stream->allocated);
    }
    // Copy data to buffer
    memcpy(stream->buf + stream->size, data, size);
    stream->size += size;
}

/**
 * @brief Same as reb_output_check but with a phase argument
 */
int reb_output_check_phase(struct reb_simulation* r, double interval, double phase) {
    double shift = r->t + interval * phase;
    if (floor(shift / interval) != floor((shift - r->dt) / interval)) {
        return 1;
    }
    // Output at beginning
    if (r->t == 0) {
        return 1;
    }
    return 0;
}

int reb_output_check(struct reb_simulation* r, double interval) {
    return reb_output_check_phase(r, interval, 0);
}

void reb_output_timing(struct reb_simulation* r, const double tmax) {
    const int N = r->N;
    int N_tot   = N;
    struct timeval tim;
    gettimeofday(&tim, NULL);
    double temp = tim.tv_sec + (tim.tv_usec / 1000000.0);
    if (r->output_timing_last == -1) {
        r->output_timing_last = temp;
    } else {
        printf("\r");
    }
    printf("N_tot= %- 9d  ", N_tot);
    if (strcmp(r->integrator_selected->name, "sei") == 0) {
        printf("t= %- 9f [orb]  ", r->t * r->Omega / 2. / M_PI);
    } else {
        printf("t= %- 9f  ", r->t);
    }
    printf("dt= %- 9f  ", r->dt);
    printf("cpu= %- 9f [s]  ", temp - r->output_timing_last);
    if (tmax > 0) {
        printf("t/tmax= %5.2f%%", r->t / tmax * 100.0);
    }
    fflush(stdout);
    r->output_timing_last = temp;
}

void reb_output_ascii(struct reb_simulation* r, char* filename) {
    const int N = r->N;
    FILE* of    = fopen(filename, "a");
    if (of == NULL) {
        reb_error(r, "Can not open file.");
        return;
    }
    for (int i = 0; i < N; i++) {
        struct reb_particle p = r->particles[i];
        fprintf(of, "%e\t%e\t%e\t%e\t%e\t%e\n", p.x, p.y, p.z, p.vx, p.vy, p.vz);
    }
    fclose(of);
}

void reb_output_orbits(struct reb_simulation* r, char* filename) {
    const int N = r->N;
    FILE* of    = fopen(filename, "a");
    if (of == NULL) {
        reb_error(r, "Can not open file.");
        return;
    }
    struct reb_particle com = r->particles[0];
    for (int i = 1; i < N; i++) {
        struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i], com);
        fprintf(of, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", r->t, o.a, o.e, o.inc, o.Omega, o.omega, o.l, o.P, o.f);
        com = reb_get_com_of_pair(com, r->particles[i]);
    }
    fclose(of);
}

void reb_output_stream_write_field(struct reb_output_stream* stream, enum REB_BINARY_FIELD_TYPE type, void* data, size_t size) {
    // Write a single field to a binary file.
    struct reb_binary_field field;
    // Memset forces padding to be set to 0 (not necessary but
    memset(&field, 0, sizeof(struct reb_binary_field));
    field.type = type;
    field.size = size;
    reb_output_stream_write(stream, &field, sizeof(struct reb_binary_field));
    reb_output_stream_write(stream, data, size);
}

void reb_output_stream_write_binary(struct reb_output_stream* stream, struct reb_simulation* r) {
    // Output header.
    char header[64] = "\0";
    int cwritten    = sprintf(header, "REBOUND Binary File. Version: %s", reb_version_str);
    snprintf(header + cwritten + 1, 64 - cwritten - 1, "%s", reb_githash_str);
    reb_output_stream_write(stream, header, sizeof(char) * 64);

    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_INTEGRATOR_SELECTED, &r->integrator_selected->id, sizeof(unsigned int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_T, &r->t, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_G, &r->G, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_OMEGA, &r->Omega, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_OMEGAZ, &r->Omega_z, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_SOFTENING, &r->softening, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_DT, &r->dt, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_DTLASTDONE, &r->dt_last_done, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_N, &r->N, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_NVAR, &r->N_var, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_VARCONFIGN, &r->var_config_N, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_NACTIVE, &r->N_active, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_TESTPARTICLETYPE, &r->testparticle_type, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_TESTPARTICLEHIDEWARNINGS, &r->testparticle_hidewarnings, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_HASHCTR, &r->hash_ctr, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_OPENINGANGLE2, &r->opening_angle2, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_STATUS, &r->status, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_EXACTFINISHTIME, &r->exact_finish_time, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_FORCEISVELOCITYDEP, &r->force_is_velocity_dependent, sizeof(unsigned int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_GRAVITYIGNORETERMS, &r->gravity_ignore_terms, sizeof(unsigned int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_OUTPUTTIMINGLAST, &r->output_timing_last, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_SAVEMESSAGES, &r->save_messages, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_EXITMAXDISTANCE, &r->exit_max_distance, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_EXITMINDISTANCE, &r->exit_min_distance, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_USLEEP, &r->usleep, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_TRACKENERGYOFFSET, &r->track_energy_offset, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_ENERGYOFFSET, &r->energy_offset, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_BOXSIZE, &r->boxsize, sizeof(struct reb_vec3d));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_BOXSIZEMAX, &r->boxsize_max, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_ROOTSIZE, &r->root_size, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_ROOTN, &r->root_n, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_ROOTNX, &r->root_nx, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_ROOTNY, &r->root_ny, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_ROOTNZ, &r->root_nz, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_NGHOSTX, &r->nghostx, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_NGHOSTY, &r->nghosty, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_NGHOSTZ, &r->nghostz, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_COLLISIONRESOLVEKEEPSORTED, &r->collision_resolve_keep_sorted, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_MINIMUMCOLLISIONVELOCITY, &r->minimum_collision_velocity, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_COLLISIONSPLOG, &r->collisions_plog, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_MAXRADIUS, &r->max_radius, 2 * sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_COLLISIONSNLOG, &r->collisions_Nlog, sizeof(long));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_CALCULATEMEGNO, &r->calculate_megno, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_MEGNOYS, &r->megno_Ys, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_MEGNOYSS, &r->megno_Yss, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_MEGNOCOVYT, &r->megno_cov_Yt, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_MEGNOVART, &r->megno_var_t, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_MEGNOMEANT, &r->megno_mean_t, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_MEGNOMEANY, &r->megno_mean_Y, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_MEGNON, &r->megno_n, sizeof(long));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_SAAUTOINTERVAL, &r->simulationarchive_auto_interval, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_SAAUTOWALLTIME, &r->simulationarchive_auto_walltime, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_SANEXT, &r->simulationarchive_next, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_WALLTIME, &r->walltime, sizeof(double));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_COLLISION, &r->collision, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_VISUALIZATION, &r->visualization, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_BOUNDARY, &r->boundary, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_GRAVITY, &r->gravity, sizeof(int));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_PYTHON_UNIT_L, &r->python_unit_l, sizeof(uint32_t));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_PYTHON_UNIT_M, &r->python_unit_m, sizeof(uint32_t));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_PYTHON_UNIT_T, &r->python_unit_t, sizeof(uint32_t));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_STEPSDONE, &r->steps_done, sizeof(unsigned long long));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_SAAUTOSTEP, &r->simulationarchive_auto_step, sizeof(unsigned long long));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_SANEXTSTEP, &r->simulationarchive_next_step, sizeof(unsigned long long));
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_RAND_SEED, &r->rand_seed, sizeof(unsigned int));
    int functionpointersused = 0;
    if (r->coefficient_of_restitution ||
        r->collision_resolve ||
        r->additional_forces ||
        r->heartbeat ||
        r->post_timestep_modifications ||
        r->free_particle_ap) {
        functionpointersused = 1;
    }
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_FUNCTIONPOINTERS, &functionpointersused, sizeof(int));
    {
        struct reb_binary_field field;
        memset(&field, 0, sizeof(struct reb_binary_field));
        field.type = REB_BINARY_FIELD_TYPE_PARTICLES;
        field.size = sizeof(struct reb_particle) * r->N;
        reb_output_stream_write(stream, &field, sizeof(struct reb_binary_field));
        // output one particle at a time to sanitize pointers.
        for (int l = 0; l < r->N; l++) {
            struct reb_particle op = r->particles[l];
            op.c                   = NULL;
            op.ap                  = NULL;
            op.sim                 = NULL;
            reb_output_stream_write(stream, &op, sizeof(struct reb_particle));
        }
    }
    if (r->var_config) {
        reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_VARCONFIG, r->var_config, sizeof(struct reb_variational_configuration) * r->var_config_N);
    }
    for (int i = 0; i < r->integrators_available_N; i++) {
        struct reb_integrator* integrator = &(r->integrators_available[i]);
        integrator->save(integrator, r, stream);
    }
    int end_null = 0;
    reb_output_stream_write_field(stream, REB_BINARY_FIELD_TYPE_END, &end_null, 0);
    struct reb_simulationarchive_blob blob = {0};
    reb_output_stream_write(stream, &blob, sizeof(struct reb_simulationarchive_blob));
}

void reb_output_binary(struct reb_simulation* r, const char* filename) {
    FILE* of = fopen(filename, "wb");
    if (of == NULL) {
        reb_error(r, "Can not open file.");
        return;
    }

    struct reb_output_stream stream;
    stream.buf       = NULL;
    stream.size      = 0;
    stream.allocated = 0;

    reb_output_stream_write_binary(&stream, r);
    fwrite(stream.buf, stream.size, 1, of);

    free(stream.buf);
    fclose(of);
}

void reb_output_binary_positions(struct reb_simulation* r, const char* filename) {
    const int N = r->N;
    FILE* of    = fopen(filename, "wb");
    if (of == NULL) {
        reb_error(r, "Can not open file.");
        return;
    }
    for (int i = 0; i < N; i++) {
        struct reb_vec3d v;
        v.x = r->particles[i].x;
        v.y = r->particles[i].y;
        v.z = r->particles[i].z;
        fwrite(&(v), sizeof(struct reb_vec3d), 1, of);
    }
    fclose(of);
}

void reb_output_velocity_dispersion(struct reb_simulation* r, char* filename) {
    const int N = r->N;
    // Algorithm with reduced roundoff errors (see wikipedia)
    struct reb_vec3d A = {.x = 0, .y = 0, .z = 0};
    struct reb_vec3d Q = {.x = 0, .y = 0, .z = 0};
    for (int i = 0; i < N; i++) {
        struct reb_vec3d Aim1 = A;
        struct reb_particle p = r->particles[i];
        A.x                   = A.x + (p.vx - A.x) / (double)(i + 1);
        if (strcmp(r->integrator_selected->name, "sei") == 0) {
            A.y = A.y + (p.vy + 1.5 * r->Omega * p.x - A.y) / (double)(i + 1);
        } else {
            A.y = A.y + (p.vy - A.y) / (double)(i + 1);
        }
        A.z = A.z + (p.vz - A.z) / (double)(i + 1);
        Q.x = Q.x + (p.vx - Aim1.x) * (p.vx - A.x);
        if (strcmp(r->integrator_selected->name, "sei") == 0) {
            Q.y = Q.y + (p.vy + 1.5 * r->Omega * p.x - Aim1.y) * (p.vy + 1.5 * r->Omega * p.x - A.y);
        } else {
            Q.y = Q.y + (p.vy - Aim1.y) * (p.vy - A.y);
        }
        Q.z = Q.z + (p.vz - Aim1.z) * (p.vz - A.z);
    }
    int N_tot              = N;
    struct reb_vec3d A_tot = A;
    struct reb_vec3d Q_tot = Q;
    Q_tot.x                = sqrt(Q_tot.x / (double)N_tot);
    Q_tot.y                = sqrt(Q_tot.y / (double)N_tot);
    Q_tot.z                = sqrt(Q_tot.z / (double)N_tot);
    FILE* of               = fopen(filename, "a");
    if (of == NULL) {
        reb_error(r, "Can not open file.");
        return;
    }
    fprintf(of, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", r->t, A_tot.x, A_tot.y, A_tot.z, Q_tot.x, Q_tot.y, Q_tot.z);
    fclose(of);
}
