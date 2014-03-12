/*
SWOCS: a software to check the numerical performance of different models in
	channel or furrow flows

Copyright 2011-2012, Javier Burguete Tolosa.

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
 * \file mesh.c
 * \brief Source file to define a mesh.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2012, Javier Burguete Tolosa.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "config.h"
#include "channel.h"
#include "node.h"
#include "mesh.h"

/**
 * \define DEBUG_MESH
 * \brief Macro to debug the mesh functions.
 */
#define DEBUG_MESH 0

/**
 * \fn int mesh_open(Mesh *mesh, Channel *channel)
 * \brief Function to open a mesh.
 * \param mesh
 * \brief mesh struct.
 * \param channel
 * \brief channel struct.
 * \return 0 on error, 1 on success.
 */
int mesh_open(Mesh *mesh, Channel *channel)
{
	unsigned int i;
	double ix;
	Node *node;
	mesh->node = node = (Node*)malloc(mesh->n * sizeof(Node));
	if (!mesh->node)
	{
		print_error("mesh: not enough memory");
		return 0;
	}
	ix = channel->length / (mesh->n - 1);
	for (i = 0; i < mesh->n; ++i)
	{
		node[i].ix = ix;
		node[i].x = i * ix;
		node_init(node + i, channel->geometry);
		memcpy(node[i].friction_coefficient, channel->friction_coefficient,
			3 * sizeof(double));
		memcpy(node[i].infiltration_coefficient,
			channel->infiltration_coefficient, 4 * sizeof(double));
		memcpy(node[i].diffusion_coefficient, channel->diffusion_coefficient,
			sizeof(double));
	}
	node[0].dx = node[mesh->n - 1].dx = 0.5 * ix;
	for (i = 0; ++i < mesh->n - 1;) node[i].dx = ix;
#if DEBUG_MESH
	for (i=0; i < mesh->n; ++i)
		printf("node %u:\nx=%lg ix=%lg dx=%lg\nzb=%lg B0=%lg Z=%lg zmax=%lg\n",
			i,
			node[i].x,
			node[i].ix,
			node[i].dx,
			node[i].zb,
			node[i].B0,
			node[i].Z,
			node[i].zmax);
#endif
	return 1;
}

/**
 * \fn void mesh_initial_conditions_dry(Mesh *mesh) 
 * \brief Function to read dry initial conditions.
 * \param mesh
 * \brief mesh struct.
 */
void mesh_initial_conditions_dry(Mesh *mesh)
{
	unsigned int i;
	Node *node = mesh->node;
	for (i = 0; i < mesh->n; ++i)
		node[i].U[0] = node[i].U[1] = node[i].U[2] = node[i].U[3] = node[i].U[4]
			= 0.; 
}

/**
 * \fn int mesh_initial_conditions_profile(Mesh *mesh, FILE *file)
 * \brief Function to read an initial conditions profile.
 * \param mesh
 * \brief mesh struct.
 * \param file
 * \brief input file.
 */
int mesh_initial_conditions_profile(Mesh *mesh, FILE *file)
{
	unsigned int i, j, n;
	double dx, *x, *A, *Q, *s;
	char *msg;
	Node *node = mesh->node;
	if (fscanf(file, "%u", &n) != 1 || n < 1)
	{
		msg = "mesh initial conditions profile: bad points number";
		goto bad2;
	}
	x = (double*)malloc(n * sizeof(double));
	A = (double*)malloc(n * sizeof(double));
	Q = (double*)malloc(n * sizeof(double));
	s = (double*)malloc(n * sizeof(double));
	if (!x || !A || !Q || !s)
	{
		msg = "mesh initial conditions profile: not enough memory";
		goto bad;
	}
	for (i = 0; i < n; ++i)
	{
		if (fscanf(file, "%lf%lf%lf%lf", x + i, A + i, Q + i, s + i) != 4 ||
			A[i] < 0. || s[i] < 0.)
		{
			msg = "mesh initial conditions profile: bad defined";
			goto bad;
		}
		if (i > 0 && x[i] < x[i - 1])
		{
			msg = "mesh initial conditions profile: bad order";
			goto bad;
		}
	} 
#if DEBUG_MESH
	for (i=0; i < n; ++i)
		printf("i=%u x=%lg A=%lg Q=%lg s=%lg\n", i, x[i], A[i], Q[i], s[i]);
#endif
	--n;
	for (i = j = 0; i < mesh->n; ++i)
	{
		while (node[i].x > x[j])
		{
			if (j < n) ++j; else break;
		}
		if (node[i].x <= x[0])
		{
			node[i].U[0] = A[0];
			node[i].U[1] = Q[0];
			node[i].U[2] = A[0] * s[0];
		}
		else if (node[i].x >= x[n])
		{
			node[i].U[0] = A[n];
			node[i].U[1] = Q[n];
			node[i].U[2] = A[n] * s[n];
		}
		else
		{
			dx = (node[i].x - x[j]) / (x[j - 1] - x[j]);
			node[i].U[0] = A[j] + dx * (A[j - 1] - A[j]);
			node[i].U[1] = Q[j] + dx * (Q[j - 1] - Q[j]);
			node[i].U[2] = node[i].U[0] * (s[j] + dx * (s[j - 1] - s[j]));
		}
		node[i].U[3] = node[i].U[4] = 0.;
	}
#if DEBUG_MESH
	for (i=0; i < mesh->n; ++i)
		printf("node:\nx=%lg A=%lg Q=%lg As=%lg\n",
			node[i].x,
			node[i].U[0],
			node[i].U[1],
			node[i].U[2]);
#endif
	return 1;

bad:
	free(x), free(A), free(Q), free(s);

bad2:
	print_error(msg);
	return 0;
}

/**
 * \fn int mesh_read(Mesh *mesh, Channel *channel, FILE *file)
 * \brief Function to read a mesh.
 * \param mesh
 * \brief mesh struct.
 * \param channel
 * \brief channel struct.
 * \param file
 * \brief input file.
 * \return 0 on error, 1 on succes.
 */
int mesh_read(Mesh *mesh, Channel *channel, FILE *file)
{
	char *msg;
	if (fscanf(file, "%u%u", &mesh->n, &mesh->type) != 2)
	{
		msg = "mesh: bad defined";
		goto bad;
	}
	if (mesh->n < 3)
	{
		msg = "mesh: bad nodes number";
		goto bad;
	}
#if DEBUG_MESH
	printf("mesh: n=%u type=%u\n", mesh->n, mesh->type);
#endif
	if (!mesh_open(mesh, channel)) return 0;
	switch (mesh->type)
	{
	case 1:
		mesh_initial_conditions_dry(mesh);
		break;
	case 2:
		if (!mesh_initial_conditions_profile(mesh, file)) return 0;
		break;
	default:
		msg = "mesh: bad type";
		goto bad;
	}
#if DEBUG_MESH
	printf("mass: water=%lg solute=%lg\n",
		mesh_water_mass(mesh), mesh_solute_mass(mesh));
#endif
	return 1;

bad:
	print_error(msg);
	return 0;
}

/**
 * \fn int mesh_write_variables(Mesh *mesh, FILE *file)
 * \brief Function to write the variables of a mesh in a file.
 * \param mesh
 * \brief mesh struct.
 * \param file
 * \brief input file.
 */
void mesh_write_variables(Mesh *mesh, FILE *file)
{
	unsigned int i;
	Node *node;
	for (i = 0; i < mesh->n; ++i)
	{
		node = mesh->node + i;
		fprintf(file, "%.14le %.14le %.14le %.14le %.14le %.14le %.14le\n",
			node->x,
			node->U[0],
			node->U[1],
			node->U[2],
			node->U[3],
			node->U[4],
			node->zb);
	}
}

/**
 * \fn int mesh_write_flows(Mesh *mesh, FILE *file)
 * \brief Function to write the flows of a mesh in a file.
 * \param mesh
 * \brief mesh struct.
 * \param file
 * \brief input file.
 */
void mesh_write_flows(Mesh *mesh, FILE *file)
{
	unsigned int i, n1;
	Node *node = mesh->node;
	n1 = mesh->n - 1;
	for (i = 0; i < n1; ++i)
	{
		fprintf(file, "%.14le %.14le %.14le %.14le %.14le\n",
			0.5 * (node[i].x + node[i + 1].x),
			(node[i + 1].U[1] * node[i + 1].u - node[i].U[1] * node[i].u)
				/ node[i].ix,
			0.5 * G * (node[i + 1].U[0] + node[i].U[0])
				* (node[i + 1].zb - node[i].zb) / node[i].ix,
			0.5 * G * (node[i + 1].U[0] + node[i].U[0])
				* (node[i + 1].h - node[i].h) / node[i].ix,
			0.25 * G * (node[i + 1].U[0] + node[i].U[0])
				* (node[i + 1].Sf + node[i].Sf));
	}
}

/**
 * \fn double mesh_water_mass(Mesh *mesh)
 * \brief Function to calculate the water mass in a mesh.
 * \param mesh
 * \brief mesh struct.
 * \return water mass
 */
double mesh_water_mass(Mesh *mesh)
{
	unsigned int i;
	double mass = 0.;
	Node *node = mesh->node;
	for (i = 0; i < mesh->n; ++i)
		mass += node[i].dx * (node[i].U[0] + node[i].U[3]);
	return mass;
}

/**
 * \fn double mesh_solute_mass(Mesh *mesh)
 * \brief Function to calculate the solute mass in a mesh.
 * \param mesh
 * \brief mesh struct.
 * \return solute mass.
 */
double mesh_solute_mass(Mesh *mesh)
{
	unsigned int i;
	double mass = 0.;
	Node *node = mesh->node;
	for (i = 0; i < mesh->n; ++i)
		mass += node[i].dx * (node[i].U[2] + node[i].U[4]);
	return mass;
}
