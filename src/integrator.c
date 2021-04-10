/**
 * @file 	integrator.c
 * @brief 	Integration schemes.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @details	This file implements the leap-frog integration scheme.  
 * This scheme is second order accurate, symplectic and well suited for 
 * non-rotating coordinate systems. Note that the scheme is formally only
 * first order accurate when velocity dependent forces are present.
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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "rebound.h"
#include "gravity.h"
#include "boundary.h"
#include "tree.h"
#include "output.h"
#include "integrator.h"
   
void reb_integrator_step(struct reb_simulation* r){
	switch(r->integrator){
		case REB_INTEGRATOR_IAS15:
		case REB_INTEGRATOR_LEAPFROG:
		case REB_INTEGRATOR_SEI:
		case REB_INTEGRATOR_SABA:
		case REB_INTEGRATOR_MERCURIUS:
		case REB_INTEGRATOR_WHFAST:
		case REB_INTEGRATOR_JANUS:
		case REB_INTEGRATOR_EOS:
            r->integrator_selected->step(r->integrator_selected, r);
			break;
		default:
			break;
	}
}
    
void reb_integrator_synchronize(struct reb_simulation* r){
	switch(r->integrator){
		case REB_INTEGRATOR_IAS15:
		case REB_INTEGRATOR_LEAPFROG:
		case REB_INTEGRATOR_SEI:
		case REB_INTEGRATOR_WHFAST:
		case REB_INTEGRATOR_MERCURIUS:
		case REB_INTEGRATOR_SABA:
		case REB_INTEGRATOR_JANUS:
		case REB_INTEGRATOR_EOS:
            if (r->integrator_selected->synchronize){
                r->integrator_selected->synchronize(r->integrator_selected, r);
            }
			break;
		default:
			break;
	}
}

void reb_integrator_reset(struct reb_simulation* r){
	r->integrator = REB_INTEGRATOR_IAS15;
	r->gravity_ignore_terms = 0;
	//reb_integrator_ias15_reset(r);
	//reb_integrator_mercurius_reset(r);
	//reb_integrator_sei_config_free(r->sei_config);   // TODO!
	//r->sei_config = reb_integrator_sei_config_alloc();
	//reb_integrator_janus_reset(r);
	//reb_integrator_eos_reset(r);
}

void reb_update_acceleration(struct reb_simulation* r){
    // Update and simplify tree. 
    // This function also creates the tree if called for the first time.
    if (r->tree_needs_update || r->gravity==REB_GRAVITY_TREE || r->collision==REB_COLLISION_TREE || r->collision==REB_COLLISION_LINETREE){
        // Check for root crossings.
        reb_boundary_check(r);     

        // Update tree (this will remove particles which left the box)
        reb_tree_update(r);          
    }
    if (r->tree_root!=NULL && r->gravity==REB_GRAVITY_TREE){
        // Update center of mass and quadrupole moments in tree in preparation of force calculation.
        reb_tree_update_gravity_data(r); 
    }

    // Main force calculation:
	reb_calculate_acceleration(r);
	if (r->N_var){
		reb_calculate_acceleration_var(r);
	}

	if (r->additional_forces){
        r->additional_forces(r);
    }
}

