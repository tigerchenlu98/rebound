/**
 * @file    rebound.c
 * @brief   Main REBOUND control structures and routine, iteration loop.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section LICENSE
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
#include "rebound.h"
#include "binarydiff.h"
#include "boundary.h"
#include "collision.h"
#include "display.h"
#include "gravity.h"
#include "input.h"
#include "output.h"
#include "particle.h"
#include "simulationarchive.h"
#include "tools.h"
#include "tree.h"
#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#ifdef OPENMP
#include <omp.h>
#endif
#define MAX(a, b) ((a) < (b) ? (b) : (a)) ///< Returns the maximum of a and b
#define STRINGIFY(s) str(s)
#define str(s) #s

const int reb_max_messages_length = 1024; // needs to be constant expression for array size
const int reb_max_messages_N      = 10;
const char* reb_build_str         = __DATE__ " " __TIME__; // Date and time build string.
const char* reb_version_str       = "3.16.0";              // **VERSIONLINE** This line gets updated automatically. Do not edit manually.
const char* reb_githash_str       = STRINGIFY(GITHASH);    // This line gets updated automatically. Do not edit manually.

void reb_simulation_set_integrator(struct reb_simulation* r, const char* name) {
    for (int i = 0; i < r->integrators_available_N; i++) {
        if (strcmp(r->integrators_available[i].name, name) == 0) {
            r->integrator_selected = &(r->integrators_available[i]);
            return;
        }
    }
    printf("Error: Integrator not found."); // TODO
    exit(EXIT_FAILURE);
}

void reb_update_acceleration(struct reb_simulation* r) {
    // Update and simplify tree.
    // This function also creates the tree if called for the first time.
    if (r->tree_needs_update || r->gravity == REB_GRAVITY_TREE || r->collision == REB_COLLISION_TREE || r->collision == REB_COLLISION_LINETREE) {
        // Check for root crossings.
        reb_boundary_check(r);

        // Update tree (this will remove particles which left the box)
        reb_tree_update(r);
    }
    if (r->tree_root != NULL && r->gravity == REB_GRAVITY_TREE) {
        // Update center of mass and quadrupole moments in tree in preparation of force calculation.
        reb_tree_update_gravity_data(r);
    }

    // Main force calculation:
    reb_calculate_acceleration(r);
    if (r->N_var) {
        reb_calculate_acceleration_var(r);
    }

    if (r->additional_forces) {
        r->additional_forces(r);
    }
}

void* reb_simulation_get_integrator_config(struct reb_simulation* r, const char* name) {
    for (int i = 0; i < r->integrators_available_N; i++) {
        if (strcmp(r->integrators_available[i].name, name) == 0) {
            return r->integrator_selected->config;
        }
    }
    printf("Error: Integrator not found."); // TODO
    exit(EXIT_FAILURE);
}

struct reb_integrator* reb_simulation_register_integrator(struct reb_simulation* r, const char* name, int id) {
    for (int i = 0; i < r->integrators_available_N; i++) {
        if (r->integrators_available[i].id == id) {
            printf("Error: Integrator ID already registered.");
            exit(EXIT_FAILURE);
        }
        if (strcmp(r->integrators_available[i].name, name) == 0) {
            printf("Error: Integrator name already registered.");
            exit(EXIT_FAILURE);
        }
    }

    struct reb_integrator integrator = {0};
    integrator.name                  = name;
    integrator.id                    = id;

    r->integrators_available_N++;
    r->integrators_available                                 = realloc(r->integrators_available, sizeof(struct reb_integrator) * r->integrators_available_N);
    r->integrators_available[r->integrators_available_N - 1] = integrator;
    return &(r->integrators_available[r->integrators_available_N - 1]);
}

void reb_exit(const char* const msg) {
    // This function should also kill all children.
    // Not implemented as pid is not easy to get to.
    // kill(pid, SIGKILL);
    fprintf(stderr, "\n\033[1mFatal error! Exiting now.\033[0m %s\n", msg);
    exit(EXIT_FAILURE);
}

void reb_message(struct reb_simulation* const r, char type, const char* const msg) {
    int save_messages = 0;
    if (r != NULL) {
        save_messages = r->save_messages;
    }
    if (!save_messages || strlen(msg) >= reb_max_messages_length) {
        if (type == 'w') {
            fprintf(stderr, "\n\033[1mWarning!\033[0m %s\n", msg);
        } else if (type == 'e') {
            fprintf(stderr, "\n\033[1mError!\033[0m %s\n", msg);
        }
    } else {
        if (r->messages == NULL) {
            r->messages = calloc(reb_max_messages_N, sizeof(char*));
        }
        int n = 0;
        for (; n < reb_max_messages_N; n++) {
            if (r->messages[n] == NULL) {
                break;
            }
        }
        if (n == reb_max_messages_N) {
            free(r->messages[0]);
            for (int i = 0; i < reb_max_messages_N - 1; i++) {
                r->messages[i] = r->messages[i + 1];
            }
            r->messages[reb_max_messages_N - 1] = NULL;
            n                                   = reb_max_messages_N - 1;
        }
        r->messages[n]    = malloc(sizeof(char*) * reb_max_messages_length);
        r->messages[n][0] = type;
        strcpy(r->messages[n] + 1, msg);
    }
}

void reb_warning(struct reb_simulation* const r, const char* const msg) {
    reb_message(r, 'w', msg);
}

void reb_error(struct reb_simulation* const r, const char* const msg) {
    reb_message(r, 'e', msg);
}

int reb_get_next_message(struct reb_simulation* const r, char* const buf) {
    if (r->messages) {
        char* w0 = r->messages[0];
        if (w0) {
            for (int i = 0; i < reb_max_messages_N - 1; i++) {
                r->messages[i] = r->messages[i + 1];
            }
            r->messages[reb_max_messages_N - 1] = NULL;
            strcpy(buf, w0);
            free(w0);
            return 1;
        }
    }
    return 0;
}


void reb_configure_box(struct reb_simulation* const r, const double root_size, const int root_nx, const int root_ny, const int root_nz) {
    r->root_size = root_size;
    r->root_nx   = root_nx;
    r->root_ny   = root_ny;
    r->root_nz   = root_nz;
    // Setup box sizes
    r->boxsize.x   = r->root_size * (double)r->root_nx;
    r->boxsize.y   = r->root_size * (double)r->root_ny;
    r->boxsize.z   = r->root_size * (double)r->root_nz;
    r->root_n      = r->root_nx * r->root_ny * r->root_nz;
    r->boxsize_max = MAX(r->boxsize.x, MAX(r->boxsize.y, r->boxsize.z));
    if (r->root_nx <= 0 || r->root_ny <= 0 || r->root_nz <= 0) {
        reb_exit("Number of root boxes must be greater or equal to 1 in each direction.");
    }
}




#ifdef OPENMP
void reb_omp_set_num_threads(int num_threads) {
    omp_set_num_threads(num_threads);
}
#endif // OPENMP

const char* reb_logo[26] = {
    "          _                           _  ",
    "         | |                         | | ",
    " _ __ ___| |__   ___  _   _ _ __   __| | ",
    "| '__/ _ \\ '_ \\ / _ \\| | | | '_ \\ / _` | ",
    "| | |  __/ |_) | (_) | |_| | | | | (_| | ",
    "|_|  \\___|_.__/ \\___/ \\__,_|_| |_|\\__,_| ",
    "                                         ",
    "              `-:://::.`                 ",
    "          `/oshhoo+++oossso+:`           ",
    "       `/ssooys++++++ossssssyyo:`        ",
    "     `+do++oho+++osssso++++++++sy/`      ",
    "    :yoh+++ho++oys+++++++++++++++ss.     ",
    "   /y++hooyyooshooo+++++++++++++++oh-    ",
    "  -dsssdssdsssdssssssssssooo+++++++oh`   ",
    "  ho++ys+oy+++ho++++++++oosssssooo++so   ",
    " .d++oy++ys+++oh+++++++++++++++oosssod   ",
    " -h+oh+++yo++++oyo+++++++++++++++++oom   ",
    " `d+ho+++ys+++++oys++++++++++++++++++d   ",
    "  yys++++oy+++++++oys+++++++++++++++s+   ",
    "  .m++++++h+++++++++oys++++++++++++oy`   ",
    "   -yo++++ss++++++++++oyso++++++++oy.    ",
    "    .ss++++ho+++++++++++osys+++++yo`     ",
    "      :ss+++ho+++++++++++++osssss-       ",
    "        -ossoys++++++++++++osso.         ",
    "          `-/oyyyssosssyso+/.            ",
    "                ``....`                  ",
};
