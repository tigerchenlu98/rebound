/**
 * @file 	integrator_leapfrog.c
 * @brief 	Leap-frog integration scheme.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the leap-frog integration scheme.
 * This scheme is second order accurate, symplectic and well suited for
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
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

#include "integrator_leapfrog.h"
#include "rebound.h"

void reb_integrator_leapfrog_step(struct reb_integrator* integrator,
                                  struct reb_simulation* r) {
    r->gravity_ignore_terms = 0;
    const int N             = r->N;
    const double dt         = r->dt;

    struct reb_particle* restrict const particles = r->particles;
#pragma omp parallel for schedule(guided)
    for (int i = 0; i < N; i++) {
        particles[i].x += 0.5 * dt * particles[i].vx;
        particles[i].y += 0.5 * dt * particles[i].vy;
        particles[i].z += 0.5 * dt * particles[i].vz;
    }
    r->t += dt / 2.;

    reb_update_acceleration(r);

#pragma omp parallel for schedule(guided)
    for (int i = 0; i < N; i++) {
        particles[i].vx += dt * particles[i].ax;
        particles[i].vy += dt * particles[i].ay;
        particles[i].vz += dt * particles[i].az;
        particles[i].x += 0.5 * dt * particles[i].vx;
        particles[i].y += 0.5 * dt * particles[i].vy;
        particles[i].z += 0.5 * dt * particles[i].vz;
    }

    r->t += dt / 2.;
    r->dt_last_done = r->dt;
}

void reb_integrator_leapfrog_register(struct reb_simulation* r) {
    struct reb_integrator* integrator =
        reb_simulation_register_integrator(r, "leapfrog", 4);
    integrator->step = reb_integrator_leapfrog_step;
    // All other function pointers are NULL by default
}
