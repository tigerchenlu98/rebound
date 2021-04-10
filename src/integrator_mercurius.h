/**
 * @file 	integrator_mercurius.h
 * @brief 	Interface for numerical particle integrator
 * @author 	Hanno Rein 
 * 
 * @section 	LICENSE
 * Copyright (c) 2017 Hanno Rein
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
#ifndef _INTEGRATOR_MERCURIUS_H
#define _INTEGRATOR_MERCURIUS_H

/**
 * @brief This structure contains variables and pointer used by the MERCURIUS integrator.
 */
struct reb_integrator_mercurius_config {
   /**
    * @brief This is a function pointer to the force switching function used.
    * @details If NULL (the default), the MERCURY switching function will be used.
    * The argument d is the distance between two particles.
    * The argument dcrit is the maximum of the critical distances of the two particles.
    * The return value is a scalar between 0 and 1. If it always returns 1, then the
    * integrator becomes the standard Wisdom-Holman integrator.
    */
    double (*L) (const struct reb_simulation* const r, double d, double dcrit);  
    
    /** 
     * @brief Switching distance in units of the Hill radius 
     * @brief The switching distances for particles are calculated automatically
     * based on multiple criteria. One criterion calculates the Hill radius of 
     * particles and then multiplies it with the hillfac variable. 
     */ 
    double hillfac;        
    
    /** 
     * @brief Setting this flag to one will recalculate heliocentric coordinates from the particle structure at the beginning of the next timestep. 
     * @details After one timestep, the flag gets set back to 0. 
     * If you want to change particles after every timestep, you 
     * also need to set this flag to 1 before every timestep.
     * Default is 0.
     */ 
    unsigned int recalculate_coordinates_this_timestep;

    /** 
     * @brief Setting this flag to one will recalculate the critical switchover 
     * distances dcrit at the beginning of the next timestep. 
     * @details After one timestep, the flag gets set back to 0. 
     * If you want to recalculate dcrit at every timestep, you 
     * also need to set this flag to 1 before every timestep.
     * Default is 0.
     */ 
    unsigned int recalculate_dcrit_this_timestep;

    /**
     * @brief If this flag is set (the default), the integrator will 
     * recalculate heliocentric coordinates and synchronize after
     * every timestep, to avoid problems with outputs or particle modifications
     * between timesteps. 
     * @details Setting it to 0 will result in a speedup, but care
     * must be taken to synchronize and recalculate coordinates when needed.
     */
    unsigned int safe_mode;
    
    unsigned int is_synchronized;   ///< Flag to determine if current particle structure is synchronized
    unsigned int mode;              ///< Internal. 0 if WH is operating, 1 if IAS15 is operating.
    unsigned int encounterN;        ///< Number of particles currently having an encounter
    unsigned int encounterNactive;  ///< Number of particles currently having an encounter
    unsigned int tponly_encounter;  ///< Flag to determine if any of the encounters are between two massive bodies (0) or only involve test particles (1). Internal use only.
    unsigned int allocatedN;        ///< Current size of allocated internal arrays
    unsigned int allocatedN_additionalforces;        ///< Current size of allocated internal particles_backup_additionalforces array
    unsigned int dcrit_allocatedN;  ///< Current size of dcrit arrays
    double* dcrit;                  ///< Switching radii for particles
    struct reb_particle* particles_backup;     ///< Internal array, contains coordinates before Kepler step for encounter prediction
    struct reb_particle* particles_backup_additionalforces;     ///< Internal array, contains coordinates before Kepler step for encounter prediction
    int* encounter_map;             ///< Map to represent which particles are integrated with ias15
    struct reb_vec3d com_pos;       ///< Used internally to keep track of the centre of mass during the timestep
    struct reb_vec3d com_vel;       ///< Used internally to keep track of the centre of mass during the timestep
};


/**
 * @brief A force switching function for the MERCURIUS integrator. This function implements 
 * the same polynomial switching function used in MERCURY. 
 */
double reb_integrator_mercurius_L_mercury(const struct reb_simulation* const r, double d, double dcrit);           

/**
 * @brief A force switching function for the MERCURIUS integrator. This function implements 
 * an infinitely differentiable switching function. 
 */
double reb_integrator_mercurius_L_infinite(const struct reb_simulation* const r, double d, double dcrit);           

void reb_integrator_mercurius_register(struct reb_simulation* r);
#endif
