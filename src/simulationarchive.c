/**
 * @file 	simulationarchive.c
 * @brief 	Tools for creating and reading Simulation Archive binary files.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2016 Hanno Rein
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
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include "rebound.h"

#include "binarydiff.h"
#include "input.h"
#include "integrator_ias15.h"
#include "output.h"
#include "particle.h"
#include "tools.h"

// TODO: Naming very inconsistent!
void reb_create_simulation_from_simulationarchive_with_messages(struct reb_simulation* r, struct reb_simulationarchive* sa, long snapshot, enum reb_input_binary_messages* warnings) {
    FILE* inf = sa->inf;
    if (inf == NULL) {
        *warnings |= REB_INPUT_BINARY_ERROR_FILENOTOPEN;
        return;
    }
    if (snapshot < 0)
        snapshot += sa->nblobs;
    if (snapshot >= sa->nblobs || snapshot < 0) {
        *warnings |= REB_INPUT_BINARY_ERROR_OUTOFRANGE;
        return;
    }

    // load original binary file
    reb_simulation_destroy(r); // This might not be necessary
    reb_simulation_init(r);    // This neither
    r->simulationarchive_filename = NULL;

    fseek(inf, 0, SEEK_SET);
    struct reb_input_stream stream = {0};
    stream.file_stream             = inf;
    while (reb_input_field(r, &stream, warnings)) {
    }

    // Done?
    if (snapshot == 0)
        return;

    // Read SA snapshot
    if (fseek(inf, sa->offset[snapshot], SEEK_SET)) {
        *warnings |= REB_INPUT_BINARY_ERROR_SEEK;
        reb_simulation_free(r);
        return;
    }

    // Read in all fields
    while (reb_input_field(r, &stream, warnings)) {
    }

    return;
}

struct reb_simulation* reb_create_simulation_from_simulationarchive(struct reb_simulationarchive* sa, long snapshot) {
    if (sa == NULL)
        return NULL;
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    struct reb_simulation* r                = reb_simulation_new();
    reb_create_simulation_from_simulationarchive_with_messages(r, sa, snapshot, &warnings);
    r = reb_input_process_warnings(r, warnings);
    return r; // might be null if error occured
}

void reb_read_simulationarchive_with_messages(struct reb_simulationarchive* sa, const char* filename, struct reb_simulationarchive* sa_index, enum reb_input_binary_messages* warnings) {
    sa->inf = fopen(filename, "r");
    if (sa->inf == NULL) {
        *warnings |= REB_INPUT_BINARY_ERROR_NOFILE;
        return;
    }
    sa->filename = malloc(strlen(filename) + 1);
    strcpy(sa->filename, filename);

    // Get version
    fseek(sa->inf, 0, SEEK_SET);
    struct reb_binary_field field = {0};
    double t0                     = 0;
    do {
        fread(&field, sizeof(struct reb_binary_field), 1, sa->inf);
        switch (field.type) {
        case REB_BINARY_FIELD_TYPE_HEADER:
            //fseek(sa->inf,64 - sizeof(struct reb_binary_field),SEEK_CUR);
            {
                long objects = 0;
                // Input header.
                const long bufsize = 64 - sizeof(struct reb_binary_field);
                char readbuf[bufsize], curvbuf[bufsize];
                const char* header = "REBOUND Binary File. Version: ";
                sprintf(curvbuf, "%s%s", header + sizeof(struct reb_binary_field), reb_version_str);

                objects += fread(readbuf, sizeof(char), bufsize, sa->inf);
                // Note: following compares version, but ignores githash.
                if (strncmp(readbuf, curvbuf, bufsize) != 0) {
                    *warnings |= REB_INPUT_BINARY_WARNING_VERSION;
                }
            }
            break;
        case REB_BINARY_FIELD_TYPE_T:
            fread(&t0, sizeof(double), 1, sa->inf);
            break;
        case REB_BINARY_FIELD_TYPE_SAAUTOWALLTIME:
            fread(&(sa->auto_walltime), sizeof(double), 1, sa->inf);
            break;
        case REB_BINARY_FIELD_TYPE_SAAUTOINTERVAL:
            fread(&(sa->auto_interval), sizeof(double), 1, sa->inf);
            break;
        case REB_BINARY_FIELD_TYPE_SAAUTOSTEP:
            fread(&(sa->auto_step), sizeof(unsigned long long), 1, sa->inf);
            break;
        default:
            fseek(sa->inf, field.size, SEEK_CUR);
            break;
        }
    } while (field.type != REB_BINARY_FIELD_TYPE_END);

    // Make index
    if (sa_index == NULL) { // Need to construct offset index from file.
        long nblobsmax = 1024;
        sa->t          = malloc(sizeof(double) * nblobsmax);
        sa->offset     = malloc(sizeof(uint32_t) * nblobsmax);
        fseek(sa->inf, 0, SEEK_SET);
        int failed    = 0;
        int endoffile = 0;
        sa->nblobs    = 0;
        for (long i = 0; i < nblobsmax; i++) {
            struct reb_binary_field field = {0};
            sa->offset[i]                 = ftell(sa->inf);
            do {
                size_t r1 = fread(&field, sizeof(struct reb_binary_field), 1, sa->inf);
                if (r1 == 1) {
                    switch (field.type) {
                    case REB_BINARY_FIELD_TYPE_HEADER: {
                        int s1 = fseek(sa->inf, 64 - sizeof(struct reb_binary_field), SEEK_CUR);
                        if (s1) {
                            failed = 1;
                        }
                    } break;
                    case REB_BINARY_FIELD_TYPE_T: {
                        size_t r2 = fread(&(sa->t[i]), sizeof(double), 1, sa->inf);
                        if (r2 != 1) {
                            failed = 1;
                        }
                    } break;
                    default: {
                        int s2 = fseek(sa->inf, field.size, SEEK_CUR);
                        if (s2) {
                            failed = 1;
                        }
                    } break;
                    }
                } else {
                    endoffile = 1;
                }
            } while (endoffile == 0 && failed == 0 && field.type != REB_BINARY_FIELD_TYPE_END);
            if (endoffile == 0 && failed == 0) {
                struct reb_simulationarchive_blob blob = {0};
                size_t r3                              = fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, sa->inf);
                if (r3 != 1) {
                    failed = 1;
                }
            }
            if (failed) {
                break;
            }
            sa->nblobs = i;
            if (endoffile == 1) {
                break; // Normal exit
            }
            if (i == nblobsmax - 1) {
                nblobsmax += 1024;
                sa->t      = realloc(sa->t, sizeof(double) * nblobsmax);
                sa->offset = realloc(sa->offset, sizeof(uint32_t) * nblobsmax);
            }
        }
        if (failed) {
            if (sa->nblobs > 0) {
                *warnings |= REB_INPUT_BINARY_WARNING_CORRUPTFILE;
            } else {
                fclose(sa->inf);
                free(sa->filename);
                free(sa->t);
                free(sa->offset);
                free(sa);
                *warnings |= REB_INPUT_BINARY_ERROR_SEEK;
                return;
            }
        }

    } else { // reuse index from other SA
        // This is an optimzation for loading many large SAs.
        // It assumes the structure of this SA is *exactly* the same as in sa_index.
        // Unexpected behaviour if the shape is not the same.
        sa->nblobs = sa_index->nblobs;
        sa->t      = malloc(sizeof(double) * sa->nblobs);
        sa->offset = malloc(sizeof(uint32_t) * sa->nblobs);
        fseek(sa->inf, 0, SEEK_SET);
        // No need to read the large file, just copying the index.
        memcpy(sa->offset, sa_index->offset, sizeof(uint32_t) * sa->nblobs);
        memcpy(sa->t, sa_index->t, sizeof(double) * sa->nblobs);
    }
}

struct reb_simulationarchive* reb_open_simulationarchive(const char* filename) {
    struct reb_simulationarchive* sa        = malloc(sizeof(struct reb_simulationarchive));
    enum reb_input_binary_messages warnings = REB_INPUT_BINARY_WARNING_NONE;
    reb_read_simulationarchive_with_messages(sa, filename, NULL, &warnings);
    if (warnings & REB_INPUT_BINARY_ERROR_NOFILE) {
        // Don't output an error if file does not exist, just return NULL.
        free(sa);
        sa = NULL;
    } else {
        reb_input_process_warnings(NULL, warnings);
    }
    return sa;
}

void reb_close_simulationarchive(struct reb_simulationarchive* sa) {
    reb_free_simulationarchive_pointers(sa);
    free(sa);
}

void reb_free_simulationarchive_pointers(struct reb_simulationarchive* sa) {
    if (sa == NULL)
        return;
    if (sa->inf) {
        fclose(sa->inf);
    }
    free(sa->filename);
    free(sa->t);
    free(sa->offset);
}

void reb_simulationarchive_heartbeat(struct reb_simulation* const r) {
    if (r->simulationarchive_filename != NULL) {
        int modes = 0;
        if (r->simulationarchive_auto_interval != 0)
            modes++;
        if (r->simulationarchive_auto_walltime != 0.)
            modes++;
        if (r->simulationarchive_auto_step != 0)
            modes++;
        if (modes > 1) {
            reb_error(r, "Only use one of simulationarchive_auto_interval, simulationarchive_auto_walltime, or simulationarchive_auto_step");
        }
        if (r->simulationarchive_auto_interval != 0.) {
            const double sign = r->dt > 0. ? 1. : -1;
            if (sign * r->simulationarchive_next <= sign * r->t) {
                r->simulationarchive_next += sign * r->simulationarchive_auto_interval;
                //Snap
                reb_simulationarchive_snapshot(r, NULL);
            }
        }
        if (r->simulationarchive_auto_step != 0.) {
            if (r->simulationarchive_next_step <= r->steps_done) {
                r->simulationarchive_next_step += r->simulationarchive_auto_step;
                //Snap
                reb_simulationarchive_snapshot(r, NULL);
            }
        }
        if (r->simulationarchive_auto_walltime != 0.) {
            if (r->simulationarchive_next <= r->walltime) {
                r->simulationarchive_next += r->simulationarchive_auto_walltime;
                //Snap
                reb_simulationarchive_snapshot(r, NULL);
            }
        }
    }
}

void reb_simulationarchive_snapshot(struct reb_simulation* const r, const char* filename) {
    if (filename == NULL)
        filename = r->simulationarchive_filename;
    struct stat buffer;
    if (stat(filename, &buffer) < 0) {
        reb_output_binary(r, filename);
    } else {
        // File exists, append snapshot.
        // Create buffer containing original binary file
        FILE* of = fopen(filename, "r+b");
        fseek(of, 64, SEEK_SET); // Header
        struct reb_binary_field field;
        int bytesread;
        do {
            bytesread = fread(&field, sizeof(struct reb_binary_field), 1, of);
            fseek(of, field.size, SEEK_CUR);
        } while (field.type != REB_BINARY_FIELD_TYPE_END && bytesread);
        long size_old = ftell(of);
        char* buf_old = malloc(size_old);
        fseek(of, 0, SEEK_SET);
        fread(buf_old, size_old, 1, of);

        // Create buffer containing current binary file
        struct reb_output_stream new_stream;
        new_stream.buf       = NULL;
        new_stream.size      = 0;
        new_stream.allocated = 0;

        reb_output_stream_write_binary(&new_stream, r);

        // Create buffer containing diff
        struct reb_input_stream new_istream = {.mem_stream = new_stream.buf, .size = new_stream.size, .file_stream = NULL};
        struct reb_input_stream old_istream = {.mem_stream = buf_old, .size = size_old, .file_stream = NULL};
        struct reb_output_stream ostream    = {0};
        reb_binary_diff(&old_istream, &new_istream, &ostream, 0);

        // Update blob info and Write diff to binary file
        struct reb_simulationarchive_blob blob = {0};
        fseek(of, -sizeof(struct reb_simulationarchive_blob), SEEK_END);
        fread(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
        blob.offset_next = ostream.size + sizeof(struct reb_binary_field);
        fseek(of, -sizeof(struct reb_simulationarchive_blob), SEEK_END);
        fwrite(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);
        fwrite(ostream.buf, ostream.size, 1, of);
        field.type = REB_BINARY_FIELD_TYPE_END;
        field.size = 0;
        fwrite(&field, sizeof(struct reb_binary_field), 1, of);
        blob.index++;
        blob.offset_prev = blob.offset_next;
        blob.offset_next = 0;
        fwrite(&blob, sizeof(struct reb_simulationarchive_blob), 1, of);

        fclose(of);
        free(new_stream.buf);
        free(buf_old);
        free(ostream.buf);
    }
}

static int _reb_simulationarchive_automate_set_filename(struct reb_simulation* const r, const char* filename) {
    if (r == NULL)
        return -1;
    if (filename == NULL) {
        reb_error(r, "Filename missing.");
        return -1;
    }
    struct stat buffer;
    if (stat(filename, &buffer) == 0) {
        reb_warning(r, "File in use for SimulationArchive already exists. Snapshots will be appended.");
    }
    free(r->simulationarchive_filename);
    r->simulationarchive_filename = malloc((strlen(filename) + 1) * sizeof(char));
    strcpy(r->simulationarchive_filename, filename);
    return 0;
}

void reb_simulationarchive_automate_interval(struct reb_simulation* const r, const char* filename, double interval) {
    if (_reb_simulationarchive_automate_set_filename(r, filename) < 0)
        return;
    if (r->simulationarchive_auto_interval != interval) {
        // Only update simulationarchive_next if interval changed.
        // This ensures that interrupted simulations will continue
        // after being restarted from a simulationarchive
        r->simulationarchive_auto_interval = interval;
        r->simulationarchive_next          = r->t;
    }
}

void reb_simulationarchive_automate_walltime(struct reb_simulation* const r, const char* filename, double walltime) {
    if (_reb_simulationarchive_automate_set_filename(r, filename) < 0)
        return;
    // Note that this will create two snapshots if restarted.
    r->simulationarchive_auto_walltime = walltime;
    r->simulationarchive_next          = r->walltime;
}

void reb_simulationarchive_automate_step(struct reb_simulation* const r, const char* filename, unsigned long long step) {
    if (_reb_simulationarchive_automate_set_filename(r, filename) < 0)
        return;
    if (r->simulationarchive_auto_step != step) {
        // Only update simulationarchive_next if interval changed.
        // This ensures that interrupted simulations will continue
        // after being restarted from a simulationarchive
        r->simulationarchive_auto_step = step;
        r->simulationarchive_next_step = r->steps_done;
    }
}
