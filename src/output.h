/**
 * @file 	output.h
 * @brief 	Output routines.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
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
#ifndef _OUTPUT_H
#define _OUTPUT_H

struct reb_simulation;
struct reb_output_stream;

struct reb_output_stream {
    char* buf;
    size_t allocated;
    size_t size;
};

void reb_output_stream_write_binary(struct reb_output_stream*, struct reb_simulation* r);
void reb_output_stream_write_field(struct reb_output_stream* stream, enum REB_BINARY_FIELD_TYPE type, void* restrict data, size_t size);
void reb_output_stream_write(struct reb_output_stream* stream, void* restrict data, size_t size);

#endif
