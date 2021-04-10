/**
 * @file 	integrator_janus.h
 * @brief 	Interface for numerical particle integrator
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
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
#ifndef _INTEGRATOR_JANUS_H
#define _INTEGRATOR_JANUS_H
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
    struct reb_particle_int* p_int;    ///< Integer particle pos/vel
    unsigned int allocated_N;                   ///< Space allocated in arrays
    /**
     * @endcond
     */
};
void reb_integrator_janus_register(struct reb_simulation* r);

#endif
