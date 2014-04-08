/*
SWOCS: a software to check the numerical performance of different models in
	channel or furrow flows

Copyright 2011-2014, Javier Burguete Tolosa.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

	1. Redistributions of source code must retain the above copyright notice,
		this list of conditions and the following disclaimer.

	2. Redistributions in binary form must reproduce the above copyright notice,
		this list of conditions and the following disclaimer in the
		documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Javier Burguete Tolosa ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL Javier Burguete Tolosa OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * \file mesh.h
 * \brief Header file to define a mesh.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2014, Javier Burguete Tolosa.
 */

// in order to prevent multiple definitions
#ifndef MESH__H
#define MESH__H 1

/**
 * \struct _Mesh
 * \brief Struct to define a mesh.
 */
struct _Mesh
{
/**
 * \var node
 * \brief array of node structs.
 * \var n
 * \brief number of nodes.
 * \var type
 * \brief initial conditions type (1 dry, 2 longitudinal profile).
 */
	Node *node;
	int n, type;
};

/**
 * \typedef Mesh
 */
typedef struct _Mesh Mesh;

// member functions

int mesh_open(Mesh *mesh, Channel *channel);
int mesh_read(Mesh *mesh, Channel *channel, FILE *file);
void mesh_write_variables(Mesh *mesh, FILE *file);
void mesh_write_flows(Mesh *mesh, FILE *file);
double mesh_water_mass(Mesh *mesh);
double mesh_solute_mass(Mesh *mesh);

#endif
