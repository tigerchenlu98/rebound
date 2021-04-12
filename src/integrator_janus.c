/**
 * @file    integrator_janus.c
 * @brief   Janus integration scheme.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file implements the Janus integration scheme.  
 * Described in Rein & Tamayo 2017.
 * 
 * @section LICENSE
 * Copyright (c) 2017 Hanno Rein, Daniel Tamayo
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
#include <stdlib.h>

#include "rebound.h"

#include "boundary.h"
#include "gravity.h"
#include "input.h"
#include "integrator_janus.h"
#include "output.h"
#include "particle.h"
#include "tools.h"

#define REB_PARTICLE_INT_TYPE int64_t
/**
 * @brief Integer positions and velocities for particles. Used in JANUS integrator. 
 */
struct reb_particle_int {
    REB_PARTICLE_INT_TYPE x;
    REB_PARTICLE_INT_TYPE y;
    REB_PARTICLE_INT_TYPE z;
    REB_PARTICLE_INT_TYPE vx;
    REB_PARTICLE_INT_TYPE vy;
    REB_PARTICLE_INT_TYPE vz;
};

struct reb_integrator_janus_config {
    /**
     * @brief Scale of the problem. Positions get divided by this number before the conversion to an integer. 
     */
    double scale_pos;
    /**
     * @brief Scale of the problem. Velocities get divided by this number before the conversion to an integer. 
     */
    double scale_vel;
    /**
     * @brief Order of the scheme. Default is 6. 
     */
    unsigned int order; //TODO needs input/output
    /**
     * @brief If this flag is set, then janus will recalculate integer coordinates at
     * the next timestep.
     */
    unsigned int recalculate_integer_coordinates_this_timestep;
    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */
    struct reb_particle_int* p_int; ///< Integer particle pos/vel
    unsigned int allocated_N;       ///< Space allocated in arrays
    /**
     * @endcond
     */
};
void reb_integrator_janus_synchronize(struct reb_integrator* integrator, struct reb_simulation* r);

/**
 * Stucture derscribing one specific JANUS scheme
 **/
struct reb_janus_scheme {
    unsigned int order;  ///< Order of the scheme
    unsigned int stages; ///< Number of stages
    double gamma[17];    ///< Coefficients (padded with 0 if not used)
};

static struct reb_janus_scheme s1odr2 = {
    .order  = 2,
    .stages = 1,
    .gamma  = {1.,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

static struct reb_janus_scheme s5odr4 = {
    .order  = 4,
    .stages = 5,
    .gamma  = {0.41449077179437573714,
              0.41449077179437573714,
              -0.65796308717750294857,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

static struct reb_janus_scheme s9odr6a = {
    .order  = 6,
    .stages = 9,
    .gamma  = {0.39216144400731413928,
              0.33259913678935943860,
              -0.70624617255763935981,
              0.082213596293550800230,
              0.79854399093482996340,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

static struct reb_janus_scheme s15odr8 = {
    .order  = 8,
    .stages = 15,
    .gamma  = {.74167036435061295345,
              -.40910082580003159400,
              .19075471029623837995,
              -.57386247111608226666,
              .29906418130365592384,
              .33462491824529818378,
              .31529309239676659663,
              -.79688793935291635402,
              0, 0, 0, 0, 0, 0, 0, 0, 0}};

static struct reb_janus_scheme s33odr10c = {
    .order  = 10,
    .stages = 33,
    .gamma  = {0.12313526870982994083,
              0.77644981696937310520,
              0.14905490079567045613,
              -0.17250761219393744420,
              -0.54871240818800177942,
              0.14289765421841842100,
              -0.31419193263986861997,
              0.12670943739561041022,
              0.17444734584181312998,
              0.44318544665428572929,
              -0.81948900568299084419,
              0.13382545738489583020,
              0.64509023524410605020,
              -0.71936337169922060719,
              0.20951381813463649682,
              -0.26828113140636051966,
              0.83647216092348048955}};

static double gg(struct reb_janus_scheme s, unsigned int stage) {
    if (stage < (s.stages + 1) / 2) {
        return s.gamma[stage];
    } else {
        return s.gamma[s.stages - 1 - stage];
    }
}

static void to_int(struct reb_particle_int* psi, struct reb_particle* ps, unsigned int N, double scale_pos, double scale_vel) {
    for (unsigned int i = 0; i < N; i++) {
        psi[i].x  = ps[i].x / scale_pos;
        psi[i].y  = ps[i].y / scale_pos;
        psi[i].z  = ps[i].z / scale_pos;
        psi[i].vx = ps[i].vx / scale_vel;
        psi[i].vy = ps[i].vy / scale_vel;
        psi[i].vz = ps[i].vz / scale_vel;
    }
}
static void to_double(struct reb_particle* ps, struct reb_particle_int* psi, unsigned int N, double scale_pos, double scale_vel) {
    for (unsigned int i = 0; i < N; i++) {
        ps[i].x  = ((double)psi[i].x) * scale_pos;
        ps[i].y  = ((double)psi[i].y) * scale_pos;
        ps[i].z  = ((double)psi[i].z) * scale_pos;
        ps[i].vx = ((double)psi[i].vx) * scale_vel;
        ps[i].vy = ((double)psi[i].vy) * scale_vel;
        ps[i].vz = ((double)psi[i].vz) * scale_vel;
    }
}

static void drift(struct reb_simulation* r, double dt, double scale_pos, double scale_vel) {
    struct reb_integrator_janus_config* ri_janus = (struct reb_integrator_janus_config*)&(r->integrator_selected->config);
    const unsigned int N                         = r->N;
    for (unsigned int i = 0; i < N; i++) {
        ri_janus->p_int[i].x += (REB_PARTICLE_INT_TYPE)(dt * (double)ri_janus->p_int[i].vx * scale_vel / scale_pos);
        ri_janus->p_int[i].y += (REB_PARTICLE_INT_TYPE)(dt * (double)ri_janus->p_int[i].vy * scale_vel / scale_pos);
        ri_janus->p_int[i].z += (REB_PARTICLE_INT_TYPE)(dt * (double)ri_janus->p_int[i].vz * scale_vel / scale_pos);
    }
}

static void kick(struct reb_simulation* r, double dt, double scale_vel) {
    struct reb_integrator_janus_config* ri_janus = (struct reb_integrator_janus_config*)&(r->integrator_selected->config);
    const unsigned int N                         = r->N;
    for (unsigned int i = 0; i < N; i++) {
        ri_janus->p_int[i].vx += (REB_PARTICLE_INT_TYPE)(dt * r->particles[i].ax / scale_vel);
        ri_janus->p_int[i].vy += (REB_PARTICLE_INT_TYPE)(dt * r->particles[i].ay / scale_vel);
        ri_janus->p_int[i].vz += (REB_PARTICLE_INT_TYPE)(dt * r->particles[i].az / scale_vel);
    }
}

void reb_integrator_janus_step(struct reb_integrator* integrator, struct reb_simulation* r) {
    r->gravity_ignore_terms                      = 0;
    struct reb_integrator_janus_config* ri_janus = (struct reb_integrator_janus_config*)&(r->integrator_selected->config);
    const unsigned int N                         = r->N;
    const double dt                              = r->dt;
    const double scale_vel                       = ri_janus->scale_vel;
    const double scale_pos                       = ri_janus->scale_pos;
    if (ri_janus->allocated_N != N) {
        ri_janus->allocated_N                                   = N;
        ri_janus->p_int                                         = realloc(ri_janus->p_int, sizeof(struct reb_particle_int) * N);
        ri_janus->recalculate_integer_coordinates_this_timestep = 1;
    }

    if (ri_janus->recalculate_integer_coordinates_this_timestep == 1) {
        to_int(ri_janus->p_int, r->particles, N, scale_pos, scale_vel);
        ri_janus->recalculate_integer_coordinates_this_timestep = 0;
    }

    struct reb_janus_scheme s;
    switch (ri_janus->order) {
    case 2:
        s = s1odr2;
        break;
    case 4:
        s = s5odr4;
        break;
    case 6:
        s = s9odr6a;
        break;
    case 8:
        s = s15odr8;
        break;
    case 10:
        s = s33odr10c;
        break;
    default:
        s = s1odr2;
        reb_error(r, "Order not supported in JANUS.");
    }

    drift(r, gg(s, 0) * dt / 2., scale_pos, scale_vel);
    to_double(r->particles, ri_janus->p_int, r->N, scale_pos, scale_vel);

    reb_update_acceleration(r);

    switch (ri_janus->order) {
    case 2:
        s = s1odr2;
        break;
    case 4:
        s = s5odr4;
        break;
    case 6:
        s = s9odr6a;
        break;
    case 8:
        s = s15odr8;
        break;
    case 10:
        s = s33odr10c;
        break;
    default:
        s = s1odr2;
        reb_error(r, "Order not supported in JANUS.");
    }

    kick(r, gg(s, 0) * dt, scale_vel);
    for (unsigned int i = 1; i < s.stages; i++) {
        drift(r, (gg(s, i - 1) + gg(s, i)) * dt / 2., scale_pos, scale_vel);
        to_double(r->particles, ri_janus->p_int, N, scale_pos, scale_vel);
        reb_update_acceleration(r);
        kick(r, gg(s, i) * dt, scale_vel);
    }
    drift(r, gg(s, s.stages - 1) * dt / 2., scale_pos, scale_vel);

    // Small overhead here: Always get positions and velocities in floating point at
    // the end of the timestep.
    reb_integrator_janus_synchronize(integrator, r);

    r->t += r->dt;
}

void reb_integrator_janus_synchronize(struct reb_integrator* integrator, struct reb_simulation* r) {
    struct reb_integrator_janus_config* ri_janus = (struct reb_integrator_janus_config*)&(r->integrator_selected->config);
    if (ri_janus->allocated_N == r->N) {
        to_double(r->particles, ri_janus->p_int, r->N, ri_janus->scale_pos, ri_janus->scale_vel);
    }
}

enum JANUS_CONFIG {
    PINT     = 1,
    SCALEPOS = 2,
    SCALEVEL = 3,
    ORDER    = 4,
    RECALC   = 5,
};

size_t reb_integrator_janus_load(struct reb_integrator* integrator, struct reb_simulation* r, struct reb_input_stream* stream, struct reb_binary_field field) {
    struct reb_integrator_janus_config* config = (struct reb_integrator_janus_config*)integrator->config;
    REB_READ_FIELD(SCALEPOS, scale_pos);
    REB_READ_FIELD(SCALEVEL, scale_vel);
    REB_READ_FIELD(ORDER, order);
    REB_READ_FIELD(RECALC, recalculate_integer_coordinates_this_timestep);
    REB_IF_FIELD(PINT) {
        if (config->p_int) {
            free(config->p_int);
        }
        config->allocated_N = (int)(field.size / sizeof(struct reb_particle_int));
        if (field.size) {
            config->p_int = malloc(field.size);
            return reb_input_stream_fread(stream, config->p_int, field.size, 1);
        }
        return 0;
    }
    return 0;
}

void reb_integrator_janus_save(struct reb_integrator* integrator, struct reb_simulation* r, struct reb_output_stream* stream) {
    struct reb_integrator_janus_config* config = (struct reb_integrator_janus_config*)integrator->config;
    REB_WRITE_FIELD(SCALEPOS, scale_pos);
    REB_WRITE_FIELD(SCALEVEL, scale_vel);
    REB_WRITE_FIELD(ORDER, order);
    REB_WRITE_FIELD(RECALC, recalculate_integer_coordinates_this_timestep);
    REB_WRITE_FIELD_WITH_SIZE(PINT, config->p_int, sizeof(struct reb_particle_int) * config->allocated_N);
}

void* reb_integrator_janus_new(struct reb_integrator* integrator, struct reb_simulation* r) {
    struct reb_integrator_janus_config* config = calloc(1, sizeof(struct reb_integrator_janus_config));

    config->order     = 2;
    config->scale_pos = 1e-16;
    config->scale_vel = 1e-16;
    return config;
}

void reb_integrator_janus_free(struct reb_integrator* integrator, struct reb_simulation* r) {
    struct reb_integrator_janus_config* config = (struct reb_integrator_janus_config*)integrator->config;
    if (config->p_int) {
        free(config->p_int);
        config->p_int = NULL;
    }
    free(config);
    integrator->config = NULL;
}

void reb_integrator_janus_register(struct reb_simulation* r) {
    struct reb_integrator* integrator = reb_simulation_register_integrator(r, "janus", 8);

    integrator->step        = reb_integrator_janus_step;
    integrator->synchronize = reb_integrator_janus_synchronize;
    integrator->new         = reb_integrator_janus_new;
    integrator->free        = reb_integrator_janus_free;
    integrator->load        = reb_integrator_janus_load;
    integrator->save        = reb_integrator_janus_save;
}
