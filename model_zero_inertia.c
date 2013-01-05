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
 * \file model_zero_inertia.c
 * \brief Source file to define the zero-inertia model.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2012, Javier Burguete Tolosa.
 */
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "channel.h"
#include "node.h"
#include "mesh.h"
#include "model.h"
#include "model_zero_inertia.h"

/**
 * \fn void node_discharge_centre_zero_inertia_Manning(Node *node)
 * \brief Function to calculate the zero-inertia discharge with the Manning model
 *   using centred derivatives.
 * \param node
 * \brief node struct.
 */
void node_discharge_centre_zero_inertia_Manning(Node *node)
{
	double dz;
	dz = (node - 1)->zs - (node + 1)->zs;
	if (dz <= 0.) node->U[1] = 0.; else
		node->U[1] = sqrt(dz / ((node - 1)->ix + node->ix)) * node->U[0]
			* pow(node->U[0] / node->P, 2./3.) / node->friction_coefficient[0];
}

/**
 * \fn void node_discharge_right_zero_inertia_Manning(Node *node)
 * \brief Function to calculate the zero-inertia discharge with the Manning model
 *   using right derivatives.
 * \param node
 * \brief node struct.
 */
void node_discharge_right_zero_inertia_Manning(Node *node)
{
	double dz;
	dz = node->zs - (node + 1)->zs;
	if (dz <= 0.) node->U[1] = 0.; else
		node->U[1] = sqrt(dz / node->ix) * node->U[0]
			* pow(node->U[0] / node->P, 2./3.) / node->friction_coefficient[0];
}

/**
 * \fn void node_discharge_left_zero_inertia_Manning(Node *node)
 * \brief Function to calculate the zero-inertia discharge with the Manning model
 *   using left derivatives.
 * \param node
 * \brief node struct.
 */
void node_discharge_left_zero_inertia_Manning(Node *node)
{
	double dz;
	dz = (node - 1)->zs - node->zs;
	if (dz <= 0.) node->U[1] = 0.; else
		node->U[1] = sqrt(dz / (node - 1)->ix) * node->U[0]
			* pow(node->U[0] / node->P, 2./3.) / node->friction_coefficient[0];
}

/**
 * \fn void model_node_parameters_centre_zero_inertia(Model *model, Node *node)
 * \brief Function to calculate the numerical parameters of a node with the
 *   zero-inertia model using centred derivatives.
 * \param model
 * \brief model struct.
 * \param node
 * \brief node struct.
 */
void model_node_parameters_centre_zero_inertia(Model *model, Node *node)
{
	node_depth(node);
	node_width(node);
	node_perimeter(node);
	if (node->U[0] <= 0.)
	{
		node->s = node->U[1] = node->u = node->T = node->Sf = node->Kx
			= node->KxA = 0.;
	}
	else if (node->h < model->minimum_depth)
	{
		node->s = node->U[2] / node->U[0];
		node->U[1] = node->u = node->T = node->Sf = node->Kx = node->KxA = 0.;
	}
	else
	{
		node->s = node->U[2] / node->U[0];
		model->node_discharge_centre(node);
		node->u = node->U[1] / node->U[0];
		node->T = node->U[1] * node->s;
		model->node_friction(node);
		model->node_diffusion(node);
		node->KxA = node->Kx * node->U[0];
	}
	node->zs = node->zb + node->h;
	model->node_infiltration(node);
	node->Pi = node->P * node->i;
}

/**
 * \fn void model_node_parameters_right_zero_inertia(Model *model, Node *node)
 * \brief Function to calculate the numerical parameters of a node with the
 *   zero-inertia model using right derivatives.
 * \param model
 * \brief model struct.
 * \param node
 * \brief node struct.
 */
void model_node_parameters_right_zero_inertia(Model *model, Node *node)
{
	node_depth(node);
	node_width(node);
	node_perimeter(node);
	if (node->U[0] <= 0.)
	{
		node->s = node->U[1] = node->u = node->T = node->Kx = node->KxA = 0.;
	}
	else if (node->h < model->minimum_depth)
	{
		node->s = node->U[2] / node->U[0];
		node->U[1] = node->u = node->T = node->Kx = node->KxA = 0.;
	}
	else
	{
		node->s = node->U[2] / node->U[0];
		model->node_discharge_right(node);
		node->u = node->U[1] / node->U[0];
		node->T = node->U[1] * node->s;
		model->node_diffusion(node);
		node->KxA = node->Kx * node->U[0];
	}
	node->zs = node->zb + node->h;
	model->node_infiltration(node);
	node->Pi = node->P * node->i;
}

/**
 * \fn void model_node_parameters_left_zero_inertia(Model *model, Node *node)
 * \brief Function to calculate the numerical parameters of a node with the
 *   zero-inertia model using left derivatives.
 * \param model
 * \brief model struct.
 * \param node
 * \brief node struct.
 */
void model_node_parameters_left_zero_inertia(Model *model, Node *node)
{
	node_depth(node);
	node_width(node);
	node_perimeter(node);
	if (node->U[0] <= 0.)
	{
		node->s = node->U[1] = node->u = node->T = node->Kx = node->KxA = 0.;
	}
	else if (node->h < model->minimum_depth)
	{
		node->s = node->U[2] / node->U[0];
		node->U[1] = node->u = node->T = node->Kx = node->KxA = 0.;
	}
	else
	{
		node->s = node->U[2] / node->U[0];
		node->u = node->U[1] / node->U[0];
		node->T = node->U[1] * node->s;
		model->node_diffusion(node);
		node->KxA = node->Kx * node->U[0];
	}
	node->zs = node->zb + node->h;
	model->node_infiltration(node);
	node->Pi = node->P * node->i;
}

/**
 * \fn double node_1dt_max_zero_inertia(Node *node)
 * \brief Function to calculate the allowed maximum time step size in a node
 *   with the zero-inertia model.
 * \param node
 * \brief node struct.
 * \return inverse of the allowed maximum time step size.
 */
double node_1dt_max_zero_inertia(Node *node)
{
	double u;
	u =  5./3. * node->u - 4./3. * node->U[1] * sqrt(1 + node->Z * node->Z)
		/ (node->B * node->P);
	if (node->u > 0.)
		u += node->U[0] * pow(node->U[0] / node->P, 4./3.)
			/ (node->friction_coefficient[0] * node->friction_coefficient[0]
			* node->u * node->dx);
	return u / node->dx;
}

/**
 * \fn void node_flows_zero_inertia(Node *node1)
 * \brief Function to calculate the flux differences in a node with the
 *   zero-inertia model.
 * \param node1
 * \brief node struct.
 */
void node_flows_zero_inertia(Node *node1)
{
	Node *node2 = node1 + 1;
	node1->dF[0] = node2->U[1] - node1->U[1];
	node1->dF[2] = node2->T - node1->T;
}

/**
 * \fn double model_inlet_dtmax_zero_inertia(Model *model)
 * \brief Function to calculate the allowed maximum time step size at the inlet
 *   with the zero-inertia model.
 * \param model
 * \brief model struct.
 * \return allowed maximum time step size.
 */
double model_inlet_dtmax_zero_inertia(Model *model)
{
	double A, Q, h, B, c;
	Node *node = model->mesh->node;
	Q = hydrogram_discharge(model->channel->water_inlet, model->t);
	h = node_critical_depth(node, Q);
	A = h * (node->B0 + h * node->Z);
	B = node->B0 + 2 * h * node->Z;
	c = sqrt(G * A / B);
	return node->ix / c;
}
