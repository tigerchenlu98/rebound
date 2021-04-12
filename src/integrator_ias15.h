/**
 * @file    integrator_ias15.h
 * @brief   Interface for numerical particle integrator
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section LICENSE
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
#ifndef _INTEGRATOR_IAS15_H
#define _INTEGRATOR_IAS15_H
/**
 * @brief Generic 7d pointer, for internal use only (IAS15).
 */
struct reb_dp7 {
    double* p0; ///< 0 substep
    double* p1; ///< 1 substep
    double* p2; ///< 2 substep
    double* p3; ///< 3 substep
    double* p4; ///< 4 substep
    double* p5; ///< 5 substep
    double* p6; ///< 6 substep
};

/**
 * @brief This structure contains variables and pointer used by the IAS15 integrator.
 */
struct reb_integrator_ias15_config {
    /**
     * @brief This parameter controls the accuracy of the integrator.
     * @details Set to 0 to make IAS15 a non-adaptive integrator.
     * The default value is: 1e-9.
     **/
    double epsilon;

    /**
     * @brief The minimum allowed timestep.
     * @details The default value is 0 (no minimal timestep).
     * Set a finite value to this variable if the IAS15 integrator has problems
     * and the timestep becomes excessively small.
     **/
    double min_dt;

    /** 
     * @brief Flag that determines how relative acceleration error is estimated.
     * @details If set to 1, estimate the fractional error by max(acceleration_error)/max(acceleration), 
     * where max is take over all particles. If set to 0, estimate the fractional error by 
     * max(acceleration_error/acceleration).
     **/
    unsigned int epsilon_global;

    /**
     * @brief Counter how many times the iteration did not converge. 
     */
    unsigned long iterations_max_exceeded;

    int allocated_N; ///< Size of allocated arrays.

    double* at;   ///< Temporary buffer for acceleration
    double* x0;   ///<                      position (used for initial values at h=0)
    double* v0;   ///<                      velocity
    double* a0;   ///<                      acceleration
    double* csx;  ///<                      compensated summation for x
    double* csv;  ///<                      compensated summation for v
    double* csa0; ///<                      compensated summation for a

    struct reb_dp7 g;
    struct reb_dp7 b;
    struct reb_dp7 csb; ///< Compensated summation for b
    struct reb_dp7 e;

    // The following values are used for resetting the b and e coefficients if a timestep gets rejected
    struct reb_dp7 br;
    struct reb_dp7 er;

    int* map;            // map to particles (identity map for non-mercurius simulations)
    int map_allocated_N; // allocated size for map
};

void reb_integrator_ias15_register(struct reb_simulation* r);

void reb_integrator_ias15_step(struct reb_integrator* ingtegrator, struct reb_simulation* r); // TODO. Needed for mercurius. needs to be done differently
#endif
