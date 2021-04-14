/**
 * @file    integrator_saba.c
 * @brief   SABA integrator family (Laskar and Robutel 2001, Blanes et al 2013).
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details This file implements the family of symplectic integrators
 * of Laskar and Robutel (2001), Blanes et al (2013), Farres et al (2013).
 * 
 * @section LICENSE
 * Copyright (c) 2019 Hanno Rein
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

#include "integrator_saba.h"
#include "input.h"
#include "integrator_whfast.h"
#include "output.h"
#include "particle.h"
#include "rebound.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAX(a, b) ((a) < (b) ? (b) : (a)) ///< Returns the maximum of a and b
#define MIN(a, b) ((a) > (b) ? (b) : (a)) ///< Returns the minimum of a and b
enum REB_SABA_TYPE {
    REB_SABA_1        = 0x0,   // WH
    REB_SABA_2        = 0x1,   // SABA2
    REB_SABA_3        = 0x2,   // SABA3
    REB_SABA_4        = 0x3,   // SABA4
    REB_SABA_CM_1     = 0x100, // SABACM1 (Modified kick corrector)
    REB_SABA_CM_2     = 0x101, // SABACM2 (Modified kick corrector)
    REB_SABA_CM_3     = 0x102, // SABACM3 (Modified kick corrector)
    REB_SABA_CM_4     = 0x103, // SABACM4 (Modified kick corrector)
    REB_SABA_CL_1     = 0x200, // SABACL1 (lazy corrector)
    REB_SABA_CL_2     = 0x201, // SABACL2 (lazy corrector)
    REB_SABA_CL_3     = 0x202, // SABACL3 (lazy corrector)
    REB_SABA_CL_4     = 0x203, // SABACL4 (lazy corrector)
    REB_SABA_10_4     = 0x4,   // SABA(10,4), 7 stages
    REB_SABA_8_6_4    = 0x5,   // SABA(8,6,4), 7 stages
    REB_SABA_10_6_4   = 0x6,   // SABA(10,6,4), 8 stages, default
    REB_SABA_H_8_4_4  = 0x7,   // SABAH(8,4,4), 6 stages
    REB_SABA_H_8_6_4  = 0x8,   // SABAH(8,6,4), 8 stages
    REB_SABA_H_10_6_4 = 0x9,   // SABAH(10,6,4), 9 stages
};

struct reb_integrator_saba_config {
    /**
     * @brief SABA type.
     * @details Available types include: SABA1, SABA2, SABA3, SABA4, SABACM1, 
     * SABACM2, SABACM3, SABACM4, SABACL1, SABACL2, SABACL3, SABACL4, 
     * SABA(10,4), SABA(8,6,4), SABA(10,6,4), SABAH(8,4,4), SABAH(8,6,4), 
     * and SABAH(10,6,4).
     */
    enum REB_SABA_TYPE type;
    unsigned int safe_mode;       ///< Safe_mode has the same functionality as in WHFast.
    unsigned int is_synchronized; ///< Flag to determine if current particle structure is synchronized
    /**
     * @brief Flaf that determines if the inertial coordinates generated are discarded in subsequent timesteps (Jacobi coordinates are used instead).
     * @details Danger zone! Only use this flag if you are absolutely sure
     * what you are doing. This is intended for
     * simulation which have to be reproducible on a bit by bit basis.
     */
    enum REB_WHFAST_COORDINATES coordinates;
    unsigned int keep_unsynchronized;
    struct reb_particle* p_jh;
    unsigned int allocated_N;
    struct reb_particle* p_temp;
    unsigned int allocated_Ntemp;
};

// Returns the number of stages for a given type of integrator (not including corrector)
static int reb_saba_stages(const int type) {
    switch (type) {
    case REB_SABA_1:
    case REB_SABA_CM_1:
    case REB_SABA_CL_1:
        return 1;
    case REB_SABA_2:
    case REB_SABA_CM_2:
    case REB_SABA_CL_2:
        return 2;
    case REB_SABA_3:
    case REB_SABA_CM_3:
    case REB_SABA_CL_3:
        return 3;
    case REB_SABA_4:
    case REB_SABA_CM_4:
    case REB_SABA_CL_4:
        return 4;
    case REB_SABA_H_8_4_4:
        return 6;
    case REB_SABA_10_4:
    case REB_SABA_8_6_4:
        return 7;
    case REB_SABA_10_6_4:
    case REB_SABA_H_8_6_4:
        return 8;
    case REB_SABA_H_10_6_4:
        return 9;
    default:
        return 0;
    }
}

// Some coefficients appear multiple times to simplify the loop structures.
const static double reb_saba_c[10][5] = {
    {
        0.5,
    }, // SABA1
    {
        0.2113248654051871177454256097490212721762,
        0.5773502691896257645091487805019574556476,
    }, // SABA2
    {
        0.1127016653792583114820734600217600389167,
        0.3872983346207416885179265399782399610833,
    }, // SABA3
    {
        0.06943184420297371238802675555359524745214,
        0.2605776340045981552106403648947824089476,
        0.3399810435848562648026657591032446872006,
    }, // SABA4
    {
        0.04706710064597250612947887637243678556564,
        0.1847569354170881069247376193702560968574,
        0.2827060056798362053243616565541452479160,
        -0.01453004174289681837857815229683813033908,
    }, // ABA(10,4)
    {
        0.0711334264982231177779387300061549964174,
        0.241153427956640098736487795326289649618,
        0.521411761772814789212136078067994229991,
        -0.333698616227678005726562603400438876027,
    }, // ABA(8,6,4)
    {
        0.03809449742241219545697532230863756534060,
        0.1452987161169137492940200726606637497442,
        0.2076276957255412507162056113249882065158,
        0.4359097036515261592231548624010651844006,
        -0.6538612258327867093807117373907094120024,
    }, // ABA(10,6,4)
    {0.2741402689434018761640565440378637101205, -0.1075684384401642306251105297063236526845,
     -0.0480185025906016926911954171508475065370, 0.7628933441747280943044988056386148982021}, // ABAH(8,4,4)
    {0.06810235651658372084723976682061164571212, 0.2511360387221033233072829580455350680082,
     -0.07507264957216562516006821767601620052338, -0.009544719701745007811488218957217113269121,
     0.5307579480704471776340674235341732001443}, // ABAH(8,6,4)
    {0.04731908697653382270404371796320813250988, 0.2651105235748785159539480036185693201078,
     -0.009976522883811240843267468164812380613143, -0.05992919973494155126395247987729676004016,
     0.2574761120673404534492282264603316880356}, // ABAH(10,6,4)
};
const static double reb_saba_d[10][5] = {
    {
        1.,
    },
    {
        0.5,
    },
    {
        0.2777777777777777777777777777777777777778,
        0.4444444444444444444444444444444444444444,
    },
    {0.1739274225687269286865319746109997036177, 0.3260725774312730713134680253890002963823},
    {
        0.1188819173681970199453503950853885936957,
        0.2410504605515015657441667865901651105675,
        -0.2732866667053238060543113981664559460630,
        0.8267085775712504407295884329818044835997,
    }, // ABA(10,4)
    {
        0.183083687472197221961703757166430291072,
        0.310782859898574869507522291054262796375,
        -0.0265646185119588006972121379164987592663,
        0.0653961422823734184559721793911134363710,
    }, // ABA(8,6,4)
    {
        0.09585888083707521061077150377145884776921,
        0.2044461531429987806805077839164344779763,
        0.2170703479789911017143385924306336714532,
        -0.01737538195906509300561788011852699719871,
    }, // ABA(10,6,4)
    {0.6408857951625127177322491164716010349386, -0.8585754489567828565881283246356000103664,
     0.7176896537942701388558792081639989754277}, // ABAH(8,4,4)
    {0.1684432593618954534310382697756917558148, 0.4243177173742677224300351657407231801453,
     -0.5858109694681756812309015355404036521923, 0.4930499927320125053698281000239887162321}, // ABAH(8,6,4)
    {0.1196884624585322035312864297489892143852, 0.3752955855379374250420128537687503199451,
     -0.4684593418325993783650820409805381740605, 0.3351397342755897010393098942949569049275,
     0.2766711191210800975049457263356834696055}, // ABAH(10,6,4)
};
const static double reb_saba_cc[4] = {
    0.08333333333333333333333333333333333333333,  // SABAC1
    0.01116454968463011276968973577058865137738,  // SABAC2
    0.005634593363122809402267823769797538671562, // SABAC3
    0.003396775048208601331532157783492144,       // SABAC4
};

static void reb_saba_corrector_step(struct reb_simulation* r, double cc) {
    struct reb_integrator_saba_config* const config = (struct reb_integrator_saba_config*)&(r->integrator_selected->config);
    struct reb_particle* const p_j                  = config->p_jh; // TODO Need to allocat
    struct reb_particle* const particles            = r->particles;
    const int N                                     = r->N;
    switch (config->type / 0x100) {
    case 1: // modified kick
        // Calculate normal kick
        reb_transformations_jacobi_to_inertial_pos(particles, p_j, particles, N, N);
        reb_simulation_calculate_acceleration(r);
        // Calculate jerk
        reb_whfast_calculate_jerk(r, p_j);

        for (int i = 0; i < N; i++) {
            const double prefact = r->dt * r->dt;
            particles[i].ax      = prefact * p_j[i].ax;
            particles[i].ay      = prefact * p_j[i].ay;
            particles[i].az      = prefact * p_j[i].az;
        }
        reb_whfast_interaction_step(r, cc * r->dt, config->coordinates, p_j);
        break;
    case 2: // lazy corrector
    {
        // Need temporary array to store old positions
        if (config->allocated_Ntemp != N) {
            config->allocated_Ntemp = N;
            config->p_temp          = realloc(config->p_temp, sizeof(struct reb_particle) * N);
        }
        struct reb_particle* p_temp = config->p_temp;

        // Calculate normal kick
        reb_transformations_jacobi_to_inertial_pos(particles, p_j, particles, N, N);
        reb_simulation_calculate_acceleration(r);
        reb_transformations_inertial_to_jacobi_acc(particles, p_j, particles, N, N);

        // make copy of original positions and accelerations
        memcpy(p_temp, p_j, r->N * sizeof(struct reb_particle));

        // WHT96 Eq 10.6
        const double prefac1 = r->dt * r->dt / 12.;
        for (unsigned int i = 1; i < N; i++) {
            p_j[i].x += prefac1 * p_temp[i].ax;
            p_j[i].y += prefac1 * p_temp[i].ay;
            p_j[i].z += prefac1 * p_temp[i].az;
        }

        // recalculate kick
        reb_transformations_jacobi_to_inertial_pos(particles, p_j, particles, N, N);
        reb_simulation_calculate_acceleration(r);
        reb_transformations_inertial_to_jacobi_acc(particles, p_j, particles, N, N);

        const double prefact = cc * r->dt * 12.;
        for (unsigned int i = 1; i < N; i++) {
            // Lazy implementer's commutator
            p_j[i].vx += prefact * (p_j[i].ax - p_temp[i].ax);
            p_j[i].vy += prefact * (p_j[i].ay - p_temp[i].ay);
            p_j[i].vz += prefact * (p_j[i].az - p_temp[i].az);
            // reset positions
            p_j[i].x = p_temp[i].x;
            p_j[i].y = p_temp[i].y;
            p_j[i].z = p_temp[i].z;
        }
    } break;
    default:
        return;
    };
}

void reb_integrator_saba_synchronize(struct reb_integrator* integrator, struct reb_simulation* const r) {
    struct reb_integrator_saba_config* const config = (struct reb_integrator_saba_config*)&(r->integrator_selected->config);
    struct reb_particle* const p_j                  = config->p_jh; // TODO Need to allocat
    int type                                        = config->type;
    struct reb_particle* sync_pj                    = NULL;
    if (config->keep_unsynchronized) {
        sync_pj = malloc(sizeof(struct reb_particle) * r->N);
        memcpy(sync_pj, config->p_jh, r->N * sizeof(struct reb_particle));
    }
    if (config->is_synchronized == 0) {
        const int N = r->N;
        if (type >= 0x100) { // correctors on
            // Drift already done, just need corrector
            reb_saba_corrector_step(r, reb_saba_cc[type % 0x100]);
        } else {
            reb_whfast_kepler_step(r, reb_saba_c[type % 0x100][0] * r->dt, config->coordinates, p_j);
            reb_whfast_com_step(r, reb_saba_c[type % 0x100][0] * r->dt, p_j);
        }
        reb_transformations_jacobi_to_inertial_posvel(r->particles, config->p_jh, r->particles, N, N);
        if (config->keep_unsynchronized) {
            memcpy(config->p_jh, sync_pj, r->N * sizeof(struct reb_particle));
            free(sync_pj);
        } else {
            config->is_synchronized = 1;
        }
    }
}

void reb_integrator_saba_step(struct reb_integrator* integrator, struct reb_simulation* const r) {
    struct reb_integrator_saba_config* const config = (struct reb_integrator_saba_config*)&(r->integrator_selected->config);
    struct reb_particle* const p_j                  = config->p_jh; // TODO Need to allocat
    const int type                                  = config->type;
    if (r->var_config_N > 0) {
        reb_error(r, "Variational particles are not supported in the SABA integrator.");
        return;
    }
    if (config->coordinates != REB_WHFAST_COORDINATES_JACOBI) {
        reb_error(r, "SABA integrator requires ri_whfast.coordinates to be set to Jacobi coordinates.");
        return;
    }
    if (config->keep_unsynchronized == 1 && config->safe_mode == 1) {
        reb_error(r, "config->keep_unsynchronized == 1 is not compatible with safe_mode. Must set config->safe_mode = 0.");
    }
    if (type != 0x0 && type != 0x1 && type != 0x2 && type != 0x3 &&
        type != 0x100 && type != 0x101 && type != 0x102 && type != 0x103 &&
        type != 0x200 && type != 0x201 && type != 0x202 && type != 0x203 &&
        type != 0x4 && type != 0x5 && type != 0x6 &&
        type != 0x7 && type != 0x8 && type != 0x9) {
        reb_error(r, "Invalid SABA integrator type used.");
        return;
    }
    if (type >= 0x100) {
        // Force Jacobi terms to be calculated in reb_simulation_calculate_acceleration if corrector is used
        r->gravity = REB_GRAVITY_JACOBI;
    } else {
        // Otherwise can do either way
        r->gravity_ignore_terms = 1;
    }
    // TODO REIMPLEMENT
    //if (reb_integrator_whfast_init(r)){
    //    // Non recoverable error occured.
    //    return;
    // }

    // TODO REIMPLEMENT
    // Only recalculate Jacobi coordinates if needed
    //if (config->safe_mode || ri_whfast->recalculate_coordinates_this_timestep){
    //    reb_integrator_whfast_from_inertial(r);
    //    ri_whfast->recalculate_coordinates_this_timestep = 0;
    //}
    if (type >= 0x100) { // Correctors on
        if (config->is_synchronized) {
            reb_saba_corrector_step(r, reb_saba_cc[type % 0x100]);
        } else {
            reb_saba_corrector_step(r, 2. * reb_saba_cc[type % 0x100]);
        }
        // First half DRIFT step
        reb_whfast_kepler_step(r, reb_saba_c[type % 0x100][0] * r->dt, config->coordinates, p_j);
        reb_whfast_com_step(r, reb_saba_c[type % 0x100][0] * r->dt, p_j);
    } else { // Correctors off
        if (config->is_synchronized) {
            // First half DRIFT step
            reb_whfast_kepler_step(r, reb_saba_c[type % 0x100][0] * r->dt, config->coordinates, p_j);
            reb_whfast_com_step(r, reb_saba_c[type % 0x100][0] * r->dt, p_j);
        } else {
            // Combined DRIFT step
            reb_whfast_kepler_step(r, 2. * reb_saba_c[type % 0x100][0] * r->dt, config->coordinates, p_j);
            reb_whfast_com_step(r, 2. * reb_saba_c[type % 0x100][0] * r->dt, p_j);
        }
    }

    reb_integrator_whfast_to_inertial(r, config->coordinates, p_j);

    reb_simulation_calculate_acceleration(r);

    const int stages                              = reb_saba_stages(type);
    const int N                                   = r->N;
    struct reb_particle* restrict const particles = r->particles;

    reb_whfast_interaction_step(r, reb_saba_d[type % 0x100][0] * r->dt, config->coordinates, p_j);

    for (int j = 1; j < stages; j++) {
        {
            int i = j;
            if (j > stages / 2) {
                i = stages - j;
            }
            reb_whfast_kepler_step(r, reb_saba_c[type % 0x100][i] * r->dt, config->coordinates, p_j);
            reb_whfast_com_step(r, reb_saba_c[type % 0x100][i] * r->dt, p_j);
        }
        {
            int i = j;
            if (j > (stages - 1) / 2) {
                i = stages - j - 1;
            }
            reb_transformations_jacobi_to_inertial_pos(particles, p_j, particles, N, N);
            reb_simulation_calculate_acceleration(r);
            reb_whfast_interaction_step(r, reb_saba_d[type % 0x100][i] * r->dt, config->coordinates, p_j);
        }
    }

    if (config->type >= 0x100) { // correctors on
        // Always need to do drift step if correctors are turned on
        reb_whfast_kepler_step(r, reb_saba_c[type % 0x100][0] * r->dt, config->coordinates, p_j);
        reb_whfast_com_step(r, reb_saba_c[type % 0x100][0] * r->dt, p_j);
    }

    config->is_synchronized = 0;
    if (config->safe_mode) {
        reb_integrator_saba_synchronize(integrator, r);
    }

    r->t += r->dt;
    r->dt_last_done = r->dt;
}

enum SABA_FIELD_TYPES {
    SAFEMODE    = 1,
    ISSYNCHRON  = 2,
    TYPE        = 3,
    KEEPUNSYNC  = 4,
    PJ          = 5,
    COORDINATES = 6,
};

size_t reb_integrator_saba_load(struct reb_integrator* integrator, struct reb_simulation* r, struct reb_input_stream* stream, struct reb_binary_field field) {
    struct reb_integrator_saba_config* config = (struct reb_integrator_saba_config*)integrator->config;
    REB_READ_FIELD(COORDINATES, coordinates);
    REB_READ_FIELD(SAFEMODE, safe_mode);
    REB_READ_FIELD(ISSYNCHRON, is_synchronized);
    REB_READ_FIELD(TYPE, type);
    REB_READ_FIELD(KEEPUNSYNC, keep_unsynchronized);
    REB_IF_FIELD(PJ) {
        if (config->p_jh) {
            free(config->p_jh);
        }
        config->allocated_N = (int)(field.size / sizeof(struct reb_particle));
        if (field.size) {
            config->p_jh = malloc(field.size);
            return reb_input_stream_fread(stream, config->p_jh, field.size, 1);
        }
    }
    return 0;
}

void reb_integrator_saba_save(struct reb_integrator* integrator, struct reb_simulation* r, struct reb_output_stream* stream) {
    struct reb_integrator_saba_config* config = (struct reb_integrator_saba_config*)integrator->config;
    REB_WRITE_FIELD(COORDINATES, coordinates);
    REB_WRITE_FIELD(SAFEMODE, safe_mode);
    REB_WRITE_FIELD(ISSYNCHRON, is_synchronized);
    REB_WRITE_FIELD(TYPE, type);
    REB_WRITE_FIELD(KEEPUNSYNC, keep_unsynchronized);
    REB_WRITE_FIELD_WITH_SIZE(PJ, config->p_jh, sizeof(struct reb_particle) * config->allocated_N);
}

void* reb_integrator_saba_new(struct reb_integrator* integrator, struct reb_simulation* r) {
    struct reb_integrator_saba_config* config = calloc(1, sizeof(struct reb_integrator_saba_config));
    config->type                              = REB_SABA_10_6_4;
    config->safe_mode                         = 1;
    config->is_synchronized                   = 1;
    config->keep_unsynchronized               = 0;
    return config;
}

void reb_integrator_saba_free(struct reb_integrator* integrator, struct reb_simulation* r) {
    struct reb_integrator_saba_config* config = (struct reb_integrator_saba_config*)integrator->config;
    if (config->p_jh) {
        free(config->p_jh);
        config->p_jh = NULL;
    }
    free(config);
    integrator->config = NULL;
}

void reb_integrator_saba_register(struct reb_simulation* r) {
    struct reb_integrator* integrator = reb_simulation_register_integrator(r, "saba", 10);
    integrator->step                  = reb_integrator_saba_step;
    integrator->synchronize           = reb_integrator_saba_synchronize;
    integrator->new                   = reb_integrator_saba_new;
    integrator->free                  = reb_integrator_saba_free;
    integrator->load                  = reb_integrator_saba_load;
    integrator->save                  = reb_integrator_saba_save;
}
