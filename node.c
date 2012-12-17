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
 * \file node.c
 * \brief Source file to define a mesh node.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2012, Javier Burguete Tolosa.
 */
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "channel.h"
#include "node.h"

/**
 * \fn void node_depth(Node *node)
 * \brief Function to calculate the depth in a mesh node.
 * \param node
 * \brief node struct.
 */
void node_depth(Node *node)
{
	if (node->Z == 0.)
		node->h = node->A / node->B0;
	else
		node->h = (sqrt(node->B0 * node->B0 + 4. * node->A * node->Z)
			- node->B0) / (2 * node->Z);
if (isnan(node->h)) printf("A=%lg B0=%lg Z=%lg\n", node->A, node->B0, node->Z);
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
	if (node->B > 0.) node->c = sqrt(G * node->A / node->B); else node->c = 0.;
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
	node->Q = fmin(node->Q, 0.99 * node->A * node->c);
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
	while (h[0]-h[1] > critical_depth_tolerance);
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
	node->Sf = node->friction_coefficient[0] * node->friction_coefficient[0]
		* node->u * fabs(node->u) * pow(node->P / node->A, 4./3.);
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
		* pow(node->Ai / (node->infiltration_coefficient[0]
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
		* sqrt(G * node->P * node->A * fabs(node->Sf));
}

/**
 * \fn void node_inlet(Node *node, Hydrogram *water, Hydrogram *solute, double t, double t2)
 * \brief Function to calculate the inlet boundary condition.
 * \param node
 * \brief node struct.
 * \param water
 * \brief water inlet hydrogram.
 * \param solute
 * \brief solute inlet hydrogram.
 * \param t
 * \brief actual time.
 * \param t2
 * \brief next time
 */
void node_inlet
	(Node *node, Hydrogram *water, Hydrogram *solute, double t, double t2)
{
	node->A += hydrogram_integrate(water, t, t2) / node->dx;
	node->As += hydrogram_integrate(solute, t, t2) / node->dx;
	node_subcritical_discharge(node);
}

/**
 * \fn void node_outlet_closed(Node *node)
 * \brief Function to calculate a closed outlet boundary condition.
 * \param node
 * \brief node struct.
*/
void node_outlet_closed(Node *node)
{
	node->Q = 0.;
}

/**
 * \fn void node_outlet_open(Node *node)
 * \brief Function to calculate an open outlet boundary condition.
 * \param node
 * \brief node struct.
 */
void node_outlet_open(Node *node)
{
	node_depth(node);
	node_width(node);
	node_critical_velocity(node);
	node->Q = fmax(node->Q, 1.01 * node->A * node->c);
}

