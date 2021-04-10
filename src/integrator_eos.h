/**
 * @file 	integrator_eos.h
 * @brief 	Interface for numerical particle integrator
 * @author 	Hanno Rein 
 * 
 * @section 	LICENSE
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
#ifndef _INTEGRATOR_EOS_H
#define _INTEGRATOR_EOS_H

/**
 * @brief Available operator splitting methods for phi0 and phi1 in EOS integrators.
 */
enum REB_EOS_TYPE {
    REB_EOS_LF = 0x00,      // 2nd order, standard leap-frog
    REB_EOS_LF4 = 0x01,     // 4th order, three function evaluations
    REB_EOS_LF6 = 0x02,     // 6th order, nine function evaluations
    REB_EOS_LF8 = 0x03,     // 8th order, seventeen funtion evaluations, see Blanes & Casa (2016), p91
    REB_EOS_LF4_2 = 0x04,   // generalized order (4,2), two force evaluations, McLachlan 1995
    REB_EOS_LF8_6_4= 0x05,  // generalized order (8,6,4), seven force evaluations
    REB_EOS_PLF7_6_4= 0x06, // generalized order (7,6,4), three force evaluations, pre- and post-processors
    REB_EOS_PMLF4 = 0x07,   // 4th order, one modified force evaluation, pre- and post-processors, Blanes et al. (1999)
    REB_EOS_PMLF6 = 0x08,   // 6th order, three modified force evaluations, pre- and post-processors, Blanes et al. (1999)
};

struct reb_integrator_eos_config {
    enum REB_EOS_TYPE phi0;         ///< Outer operator splitting scheme
    enum REB_EOS_TYPE phi1;         ///< Inner operator splitting scheme
    unsigned int n;                 ///

    unsigned int safe_mode;         ///< If set to 0, always combine drift steps at the beginning and end of phi0. If set to 1, n needs to be bigger than 1.
    unsigned int is_synchronized;   ///< Flag to indicate if the drift step at the end of the last timestep has been taken.
};

void reb_integrator_eos_register(struct reb_simulation* r);
#endif
