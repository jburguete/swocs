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
 * \file node.c
 * \brief Source file to define a mesh node.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2014, Javier Burguete Tolosa.
 */
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "channel.h"
#include "node.h"

double (*node_normal_discharge)(Node *node, double S);

/**
 * \fn void node_init(Node *node, Geometry *geometry)
 * \brief Function to calculate the initial parameters of a mesh node.
 * \param node
 * \brief node struct.
 * \param geometry
 * \brief geometry struct.
 */
void node_init(Node *node, Geometry *geometry)
{
	unsigned int i, i2, n1;
	double x, dx;
	x = node->x;
	n1 = geometry->n - 1;
	if (x <= geometry->x[0])
	{
		node->zb = geometry->zb[0];
		node->B0 = geometry->B0[0];
		node->Z = geometry->Z[0];
		node->zmax = geometry->zmax[0];
		return;
	}
	if (x >= geometry->x[n1])
	{
		node->zb = geometry->zb[n1];
		node->B0 = geometry->B0[n1];
		node->Z = geometry->Z[n1];
		node->zmax = geometry->zmax[n1];
		return;
	}
	for (i = 0; i < n1 - 1;)
	{
		i2 = (i + n1) / 2;
		if (x > geometry->x[i2]) i = i2; else n1 = i2;
	}
	dx = (x - geometry->x[i]) / (geometry->x[i + 1] - geometry->x[i]);
	node->zb = geometry->zb[i] + dx * (geometry->zb[i + 1] - geometry->zb[i]);
	node->B0 = geometry->B0[i] + dx * (geometry->B0[i + 1] - geometry->B0[i]);
	node->Z = geometry->Z[i] + dx * (geometry->Z[i + 1] - geometry->Z[i]);
	node->zmax = geometry->zmax[i]
		+ dx * (geometry->zmax[i + 1] - geometry->zmax[i]);
}

/**
 * \fn void node_depth(Node *node)
 * \brief Function to calculate the depth in a mesh node.
 * \param node
 * \brief node struct.
 */
void node_depth(Node *node)
{
	if (node->Z == 0.)
		node->h = node->U[0] / node->B0;
	else
		node->h = (sqrt(node->B0 * node->B0 + 4. * node->U[0] * node->Z)
			- node->B0) / (2 * node->Z);
	node->zs = node->zb + node->h;
}

/**
 * \fn void node_width(Node *node)
 * \brief Function to calculate the width in a mesh node.
 * \param node
 * \brief node struct.
 */
void node_width(Node *node)
{
	node->B = node->B0 + 2 * node->Z * node->h;
}

/**
 * \fn void node_perimeter(Node *node)
 * \brief Function to calculate the wetted perimeter in a mesh node.
 * \param node
 * \brief node struct.
 */
void node_perimeter(Node *node)
{
	node->P = node->B0 + 2 * sqrt(1 + node->Z * node->Z) * node->h;
}

/**
 * \fn void node_critical_velocity(Node *node)
 * \brief  * \brief Function to calculate the critical velocity in a mesh node
 * \param node
 * \brief node struct.
 */

void node_critical_velocity(Node *node)
{
	if (node->B > 0.) node->c = sqrt(G * node->U[0] / node->B); else node->c = 0.;
}

/**
 * \fn void node_subcritical_discharge(Node *node)
 * \brief Function to force a subcritical discharge in a mesh node.
 * \param node
 * \brief node struct.
 */
void node_subcritical_discharge(Node *node)
{
	node_depth(node);
	node_width(node);
	node_critical_velocity(node);
	node->U[1] = fmin(node->U[1], 0.99 * node->U[0] * node->c);
}

/**
 * \fn double node_critical_depth(Node *node, double Q)
 * \brief Function to calculate the critical depth in a mesh node.
 * \param node
 * \brief node struct.
 * \param Q
 * \brief discharge.
 * \return critical depth.
 */
double node_critical_depth(Node *node, double Q)
{
	double h[3], A[3], B[3], u[3], c[3];

	// searching a subcritical depth

	h[0] = 1.;
	do
	{
		h[0] *= 2;
		A[0] = h[0] * (node->B0 + h[0] * node->Z);
		B[0] = node->B0 + 2 * h[0] * node->Z;
		c[0] = G * A[0] / B[0];
		u[0] = Q / A[0];
		u[0] = u[0] * u[0];
	}
	while (u[0] > c[0]);

	// searching a supercritical depth

	h[1] = h[0];
	do
	{
		h[1] *= 0.5;
		A[1] = h[1] * (node->B0 + h[1] * node->Z);
		B[1] = node->B0 + 2 * h[1] * node->Z;
		c[1] = G * A[1] / B[1];
		u[1] = Q / A[1];
		u[1] = u[1] * u[1];
	}
	while (u[1] < c[1]);

	// searching the critical depth by a middle point algorithm

	do
	{
		h[2] = 0.5 * (h[0] + h[1]);
		A[2] = h[2] * (node->B0 + h[2] * node->Z);
		B[2] = node->B0 + 2 * h[2] * node->Z;
		c[2] = G * A[2] / B[2];
		u[2] = Q / A[2];
		u[2] = u[2] * u[2];
		if (u[2] < c[1]) h[0] = h[2]; else h[1] = h[2];

	}
	while (h[0] - h[1] > critical_depth_tolerance);
	return 0.5 * (h[0] + h[1]);
}

/**
 * \fn void node_friction_Manning(Node *node)
 * \brief Function to calculate the friction slope with the Manning model.
 * \param node
 * \brief node struct.
 */
void node_friction_Manning(Node *node)
{
	node->f = node->friction_coefficient[0] * node->friction_coefficient[0]
		* pow(node->P / node->U[0], 4./3.) / (node->U[0] * node->U[0]);
	node->Sf = node->f * node->U[1] * fabs(node->U[1]);
	node->beta = 1.;
}

double node_normal_discharge_Manning(Node *node, double S)
{
	return sqrt(S) * node->U[0] * pow(node->U[0] / node->P, 2./3.)
		/ node->friction_coefficient[0];
}

/**
 * \fn void node_friction_Manning_minimizing_losses(Node *node)
 * \brief Function to calculate the friction slope with the Manning model and
 *   minimizing energy losses velocity distribution in the cross section.
 * \param node
 * \brief node struct.
 */
void node_friction_Manning_minimizing_losses(Node *node)
{
	double k;
	node->f1 = node->B0 + node->Z * node->h;
	node->f2 = node->B0 + 0.75 * node->Z * node->h;
	node->f3 = node->B0 + 0.6 * node->Z * node->h;
	k = node->friction_coefficient[0] * node->f1 / (node->f2 * node->U[0]);
	node->f = pow(node->h, -4./3.) * k * k;
	node->Sf = node->f * node->U[1] * fabs(node->U[1]);
	node->beta = 49/48. * node->f1 * node->f3 / (node->f2 * node->f2);
}

double node_normal_discharge_Manning_minimizing_losses(Node *node, double S)
{
	double f1, f2;
	f1 = node->B0 + node->Z * node->h;
	f2 = node->B0 + 0.75 * node->Z * node->h;
	return sqrt(S) * node->U[0] * pow(node->h, 2./3.) * f2
		/ (node->friction_coefficient[0] * f1);
}

/**
 * \fn void node_infiltration_KostiakovLewis(Node *node)
 * \brief Function to calculate the infiltration with the Kostiakov-Lewis model.
 * \param node
 * \brief node struct.
 */

void node_infiltration_KostiakovLewis(Node *node)
{
	node->i = node->infiltration_coefficient[2];
	if (node->infiltration_coefficient[0] == 0.) return;
	node->i += node->infiltration_coefficient[0] *
		node->infiltration_coefficient[1]
		* pow(node->U[3] / (node->infiltration_coefficient[0]
		* node->infiltration_coefficient[3]),
		1. - 1. / node->infiltration_coefficient[1]);
}

/**
 * \fn void node_diffusion_Rutherford(Node *node)
 * \brief Function to calculate the diffusion coefficient with the Rutherford
 *   model.
 * \param node
 * \brief node struct.
 */
void node_diffusion_Rutherford(Node *node)
{
	node->Kx = node->diffusion_coefficient[0]
		* sqrt(G * node->P * node->U[0] * fabs(node->Sf));
}
