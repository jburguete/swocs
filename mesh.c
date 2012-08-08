/*
SWOCS: a software to check the numerical performance of different models in
	channel or furrow flows

Copyright 2011, Javier Burguete Tolosa.

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
 * \copyright Copyright 2011, Javier Burguete Tolosa.
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
	int i;
	double ix, Z;
	Node *node;
	mesh->node = node = (Node*)malloc(mesh->n * sizeof(Node));
	if (!mesh->node)
	{
		printf("mesh: not enough memory\n");
		return 0;
	}
	ix = channel->length / (mesh->n - 1);
	Z = channel->wall_slope;
	for (i = 0; i < mesh->n; ++i)
	{
		node[i].ix = ix;
		node[i].x = i * ix;
		node[i].zb = geometry_level(channel->geometry, node[i].x);
		node[i].B0 = channel->bottom_width;
		node[i].Z = Z;
		memcpy(node[i].friction_coefficient, channel->friction_coefficient,
			3 * sizeof(double));
		memcpy(node[i].infiltration_coefficient,
			channel->infiltration_coefficient, 4 * sizeof(double));
		memcpy(node[i].diffusion_coefficient, channel->diffusion_coefficient,
			sizeof(double));
	}
	node[0].dx = node[mesh->n - 1].dx = 0.5 * ix;
	for (i = 0; ++i < mesh->n - 1;) node[i].dx = ix;
#if DEBUG_MESH_OPEN
	for (i=0; i < mesh->n; ++i)
		printf("node:\nx=%lf ix=%lf dx=%lf\nzb=%lf B0=%lf Z=%lf\n",
			node[i].x,
			node[i].ix,
			node[i].dx,
			node[i].zb,
			node[i].B0,
			node[i].Z);
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
	int i;
	Node *node = mesh->node;
	for (i = 0; i < mesh->n; ++i)
		node[i].A = node[i].Q = node[i].As = node[i].Ai = node[i].Asi = 0.; 
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
	int i, j, n;
	double dx, *x, *A, *Q, *s;
	char *msg;
	Node *node = mesh->node;
	if (fscanf(file, "%d", &n) != 1 || n < 1)
	{
		msg = "mesh initial conditions profile: bad points number\n";
		goto bad2;
	}
	x = (double*)malloc(n * sizeof(double));
	A = (double*)malloc(n * sizeof(double));
	Q = (double*)malloc(n * sizeof(double));
	s = (double*)malloc(n * sizeof(double));
	if (!x || !A || !Q || !s)
	{
		msg = "mesh initial conditions profile: not enough memory\n";
		goto bad;
	}
	for (i = 0; i < n; ++i)
	{
		if (fscanf(file, "%lf%lf%lf%lf", x + i, A + i, Q + i, s + i) != 4 ||
			A[i] < 0. || s[i] < 0.)
		{
			msg = "mesh initial conditions profile: bad defined\n";
			goto bad;
		}
		if (i > 0 && x[i] < x[i - 1])
		{
			msg = "mesh initial conditions profile: bad order\n";
			goto bad;
		}
	} 
	--n;
	for (i = j = 0; i < mesh->n; ++i)
	{
		while (node[i].x > x[j])
		{
			if (j < n) ++j; else break;
		}
		if (node[i].x <= x[0])
		{
			node[i].A = A[0];
			node[i].Q = Q[0];
			node[i].As = A[0] * s[0];
		}
		else if (node[i].x >= x[n])
		{
			node[i].A = A[n];
			node[i].Q = Q[n];
			node[i].As = A[n] * s[n];
		}
		else
		{
			dx = (node[i].x - x[j]) / (x[j + 1] - x[j]);
			node[i].A = A[j] + dx * (A[j+1] - A[j]);
			node[i].Q = Q[j] + dx * (Q[j+1] - Q[j]);
			node[i].As = node[i].A * (s[j] + dx * (s[j+1] - s[j]));
		}
		node[i].Ai = node[i].Asi = 0.;
	}
	return 1;

bad:
	free(x), free(A), free(Q), free(s);

bad2:
	printf(msg);
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
	if (fscanf(file, "%d%d", &mesh->n, &mesh->type) != 2)
	{
		msg = "mesh: bad defined\n";
		goto bad;
	}
	if (mesh->n < 3)
	{
		msg = "mesh: bad nodes number\n";
		goto bad;
	}
#if DEBUG_MODEL_READ
	printf("mesh: n=%d type=%d\n", mesh->n, mesh->type);
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
		msg = "mesh: bad type\n";
		goto bad;
	}
	return 1;

bad:
	printf(msg);
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
	int i;
	Node *node;
	for (i = 0; i < mesh->n; ++i)
	{
		node = mesh->node + i;
printf("i=%d A=%.14le\n", i, node->A);
		fprintf(file, "%.14le %.14le %.14le %.14le %.14le %.14le\n",
			node->x,
			node->A,
			node->Q,
			node->As,
			node->Ai,
			node->Asi);
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
	int i, n1;
	Node *node = mesh->node;
	n1 = mesh->n - 1;
	for (i = 0; i < n1; ++i)
	{
		fprintf(file, "%.14le %.14le %.14le %.14le %.14le\n",
			0.5 * (node[i].x + node[i + 1].x),
			(node[i + 1].Q * node[i + 1].u - node[i].Q * node[i].u)
				/ node[i].ix,
			0.5 * G * (node[i + 1].A + node[i].A)
				* (node[i + 1].zb - node[i].zb) / node[i].ix,
			0.5 * G * (node[i + 1].A + node[i].A)
				* (node[i + 1].h - node[i].h) / node[i].ix,
			0.25 * G * (node[i + 1].A + node[i].A)
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
	int i;
	double mass = 0.;
	Node *node = mesh->node;
	for (i = 0; i < mesh->n; ++i)
		mass += node[i].dx * (node[i].A + node[i].Ai);
for (i=0; i<mesh->n; ++i) printf("i=%d A=%.14le Ai=%.14le\n", i, node[i].A, node[i].Ai);
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
	int i;
	double mass = 0.;
	Node *node = mesh->node;
	for (i = 0; i < mesh->n; ++i)
		mass += node[i].dx * (node[i].As + node[i].Asi);
	return mass;
}

