/**
 * @file 	integrator_whfast.h
 * @brief 	Interface for numerical particle integrator
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2015 Hanno Rein, Daniel Tamayo
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
#ifndef _INTEGRATOR_WHFAST_H
#define _INTEGRATOR_WHFAST_H

#include "rebound.h"
enum REB_WHFAST_COORDINATES {
    REB_WHFAST_COORDINATES_JACOBI                 = 0, ///< Jacobi coordinates (default)
    REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC = 1, ///< Democratic Heliocentric coordinates
    REB_WHFAST_COORDINATES_WHDS                   = 2, ///< WHDS coordinates (Hernandez and Dehnen, 2017)
};
enum REB_WHFAST_KERNEL {
    REB_WHFAST_KERNEL_DEFAULT      = 0,
    REB_WHFAST_KERNEL_MODIFIEDKICK = 1,
    REB_WHFAST_KERNEL_COMPOSITION  = 2,
    REB_WHFAST_KERNEL_LAZY         = 3,
};

/**
 * @brief This structure contains variables used by the WHFast integrator.
 */
struct reb_integrator_whfast_config {
    /**
     * @brief This variable turns on/off different first symplectic correctors for WHFast.
     * @details These correctors remove terms of order O(eps*dt^2) 
     * - 0 (default): turns off all first correctors
     * - 3: uses third order (two-stage) first corrector 
     * - 5: uses fifth order (four-stage) first corrector 
     * - 7: uses seventh order (six-stage) first corrector 
     * - 11: uses eleventh order (ten-stage) first corrector 
     * - 17: uses 17th order (16-stage) first corrector 
     */
    unsigned int corrector;

    /**
     * @brief This variable turns on/off the second symplectic correctors for WHFast.
     * @details 
     * - 0 (default): turns off second correctors
     * - 1: uses second corrector 
     */
    unsigned int corrector2;

    /**
     * @brief This variable determines the kernel of the WHFast integrator.
     * @details 
     * - 0 (default): Uses a standard WH kick step 
     * - 1: uses the exact modified kick (for Newtonian gravity) 
     * - 2: uses the composition kernel  
     * - 3: uses the lazy implementer's modified kick   
     */
    enum REB_WHFAST_KERNEL kernel;

    /**
     * @brief Chooses the coordinate system for the WHFast algorithm. Default is Jacobi Coordinates.
     */
    enum REB_WHFAST_COORDINATES coordinates;

    /** 
     * @brief Setting this flag to one will recalculate Jacobi/heliocentric coordinates from the particle structure in the next timestep. 
     * @details After the timestep, the flag gets set back to 0. 
     * If you want to change particles after every timestep, you 
     * also need to set this flag to 1 before every timestep.
     * Default is 0.
     */
    unsigned int recalculate_coordinates_this_timestep;

    /**
     * @brief If this flag is set (the default), whfast will recalculate 
     * jacobi/heliocentric coordinates and synchronize
     * every timestep, to avoid problems with outputs or particle modifications
     * between timesteps. 
     * @details Setting it to 0 will result in a speedup, but care
     * must be taken to synchronize and recalculate jacobi coordinates when needed.
     * See AdvWHFast.ipynb in the python_tutorials folder (navigate to it on github
     * if you don't have ipython notebook installed).  The explanation is general, and
     * the python and C flags have the same names.
     */
    unsigned int safe_mode;

    /**
     * @brief Jacobi/heliocentric coordinates
     * @details This array contains the Jacobi/heliocentric
     * coordinates of all particles.
     * It is automatically filled and updated by WHfast.
     * Access this array with caution.
     */
    struct reb_particle* p_jh;

    /**
     * @brief Internal temporary array used for lazy implementer's kernel method
     */
    struct reb_particle* p_temp;

    /**
     * @brief Generate inertial coordinates at the end of the integration, but do not change the Jacobi/heliocentric coordinates
     * @details Danger zone! Only use this flag if you are absolutely sure
     * what you are doing. This is intended for
     * simulation which have to be reproducible on a bit by bit basis.
     */
    unsigned int keep_unsynchronized;

    /**
     * @cond PRIVATE
     * Internal data structures below. Nothing to be changed by the user.
     */

    unsigned int is_synchronized;                                      ///< Flag to determine if current particle structure is synchronized
    unsigned int allocated_N;                                          ///< Space allocated in p_jh array
    unsigned int allocated_Ntemp;                                      ///< Space allocated in p_temp array
    unsigned int timestep_warning;                                     ///< Counter of timestep warnings
    unsigned int recalculate_coordinates_but_not_synchronized_warning; ///< Counter of Jacobi synchronization errors
    /**
     * @endcond
     */
};

void reb_integrator_whfast_register(struct reb_simulation* r);                                                                                            ///< Internal function used to call a specific integrator
void reb_whfast_kepler_solver(const struct reb_simulation* const r, struct reb_particle* const restrict p_j, const double M, unsigned int i, double _dt); ///< Internal function (Main WHFast Kepler Solver)
void reb_whfast_calculate_jerk(struct reb_simulation* r, struct reb_particle* const jerk);                                                                ///< Calculates "jerk" term

void reb_integrator_whfast_from_inertial(struct reb_simulation* const r, enum REB_WHFAST_COORDINATES coordinates, struct reb_particle* const p_j);      ///< Internal function to the appropriate WHFast coordinates from inertial
void reb_integrator_whfast_to_inertial(struct reb_simulation* const r, enum REB_WHFAST_COORDINATES coordinates, struct reb_particle* const p_j);        ///< Internal function to move back from particular WHFast coordinates to inertial
void reb_integrator_whfast_reset(struct reb_simulation* r);                                                                                             ///< Internal function used to call a specific integrator
void reb_whfast_interaction_step(struct reb_simulation* const r, const double _dt, enum REB_WHFAST_COORDINATES coordinates, struct reb_particle* p_j);  ///< Internal function
void reb_whfast_jump_step(const struct reb_simulation* const r, const double _dt, enum REB_WHFAST_COORDINATES coordinates, struct reb_particle* p_j);   ///< Internal function
void reb_whfast_kepler_step(const struct reb_simulation* const r, const double _dt, enum REB_WHFAST_COORDINATES coordinates, struct reb_particle* p_j); ///< Internal function
void reb_whfast_com_step(const struct reb_simulation* const r, const double _dt, struct reb_particle* p_j);                                             ///< Internal function

#endif
