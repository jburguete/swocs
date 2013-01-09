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

THIS SOFTWARE IS PROVIDED BY Javier Burguete ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL Javier Burguete OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/**
 * \file model_zero_inertia_implicit.c
 * \brief Source file to define the first order upwind implicit numerical model 
 *   applied to the zero inertia model.
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
#include "model_zero_inertia_implicit.h"

/**
 * \fn double node_area_with_level(Node *node, double zs)
 * \brief Function to calculate the wetted area on a mesh node.
 * \param node
 * \brief mesh node struct.
 * \param zs
 * \brief surface level.
 * \return wetted area.
 */
double node_area_with_level(Node *node, double zs)
{
	double h;
	h = zs - node->zb;
	return h * (node->B0 + node->Z * h);
}

/**
 * \fn void model_surface_flow_zero_inertia_implicit_multiply\
 *   (double *m, double *v, double *r)
 * \brief Function to multiply a matrix by a vector of the zero inertia model.
 * \param m
 * \brief multiplying matrix.
 * \param v
 * \brief vector to multiply.
 * \param r
 * \brief resulting vector.
 */
void model_surface_flow_zero_inertia_implicit_multiply
	(double *m, double *v, double *r)
{
	r[0] = m[0] * v[0];
	r[2] = m[2] * v[0] + m[3] * v[2];
}

/**
 * \fn void model_surface_flow_zero_inertia_implicit_invert\
 *   (double *m, double *i)
 * \brief Function to invert an implicit operator of the zero inertia model.
 * \param m
 * \brief implicit operator to invert.
 * \param i
 * \brief invert operator.
 */
void model_surface_flow_zero_inertia_implicit_invert(double *m, double *i)
{
	double d;
	d = m[0] * m[3];
	i[0] = m[3] / d;
	i[2] = - m[2] / d;
	i[1] = 0.;
	i[3] = m[0] / d;
}

/**
 * \fn void model_surface_flow_zero_inertia_stabilize(Model *model)
 * \brief Function to stabilize the numerical results of the zero inertia\
 *   surface flow model.
 * \param model
 * \brief model struct.
 */
void model_surface_flow_zero_inertia_stabilize(Model *model)
{
	unsigned int i, k;
	double zmax, dA, dA2;
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;

	zmax = node[0].zs;
	for (i = k = 0; ++i < mesh->n;)
	{
		if (node[i].zs > node[i - 1].zs)
		{
			dA = node[i].U[0] - node_area_with_level(node + i, node[i - 1].zs);
			dA2 = node_area_with_level(node + i - 1, zmax) - node[i - 1].U[0];
			if (dA > dA2)
			{
				node[i - 1].U[0] += dA2;
				dA -= dA2;
				node[i].U[0] -= dA;
				if (i < mesh->n - 1)
				{
					node[i + 1].U[0] += dA - dA2;
					node_depth(node + i + 1);
				}
			}
			else
			{
				node[i - 1].U[0] += dA;
				node[i].U[0] -= dA;
			}
			node_depth(node + i - 1);
			node_depth(node + i);
			k = 1;
		}
		zmax = node[i - 1].zs;
	}
	if (k) model_parameters(model);
}

/**
 * \fn void model_surface_flow_zero_inertia_implicit(Model *model)
 * \brief Function to make the surface flow with the upwind implicit numerical
 *   scheme.
 * \param model
 * \brief model struct.
 */
void model_surface_flow_zero_inertia_implicit(Model *model)
{
	unsigned int i, j, n1, iteration;
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;
	double k, l1, l2, odt, A[9], B[9], C[9], D[3],
		inlet_contribution[3], outlet_contribution[3],
		CC[mesh->n], DD[mesh->n], EE[mesh->n];

	n1 = mesh->n - 1;

	// saving some former time step variables

	for (i = 0; i < mesh->n; ++i)
	{
		for (j = 0; j < 3; ++j) node[i].Un[j] = node[i].U[j];
		node[i].Sfn = node[i].Sf;
	}

	// explicit part

	inlet_contribution[0] = - model->dt * node[0].U[1];
	inlet_contribution[2] = - model->dt * node[0].T;
	outlet_contribution[0] = model->dt * node[n1].U[1];
	outlet_contribution[2] = model->dt * node[n1].T;

	for (i = 0; i < n1; ++i) node_flows_zero_inertia(node + i);

	// implicit part

	iteration = 0;
	odt = model->theta * model->dt;

	do
	{

		model->inlet_contribution[0] = inlet_contribution[0];
		model->inlet_contribution[2] = inlet_contribution[2];
		model->outlet_contribution[0] = outlet_contribution[0];
		model->outlet_contribution[2] = outlet_contribution[2];

		for (i = 0; i < mesh->n; ++i)
		{
			if (node[i].h <= model->minimum_depth || node[i].U[1] == 0.)
			{
				for (j = 0; j < 4; ++j) node[i].Jp[j] = 0.;
				node[i].Jn[0] = 0.;
				continue;
			}
			l1 = node[i].U[1] * (5./3. / node[i].U[0]
				- 4./3. * sqrt(1 + node[i].Z * node[i].Z)
				/ (node[i].B * node[i].P));
			l2 = 0.5 * node[i].U[0] * node[i].U[0]
				* pow(node[i].U[0] / node[i].P, 4./3.)
				/ (node[i].U[1] * node[i].friction_coefficient[0]
				* node[i].friction_coefficient[0] * node[i].B * node[i].dx);
			node[i].Jp[0] = l1;
			node[i].Jp[1] = 0.;
			node[i].Jp[2] = (l1 - node[i].u) * node[i].s;
			node[i].Jp[3] = node[i].u;
			node[i].Jn[0] = - l2;
//			node[i].Jn[0] = 0.;
		}

		// variables updating

		for (j = 0; j < 4; ++j) B[j] = odt * node[0].Jp[j];
		node[0].dU[0] = node[0].dU[2] = 0.;
		node[0].U[0] = node[0].Un[0];
		node[0].U[2] = node[0].Un[2];
		for (i = 0; ++i <= n1;)
		{
			model_surface_flow_zero_inertia_implicit_multiply
				(B, node[i - 1].dU, D);
			for (j = 0; j < 4; ++j) A[j] = B[j] = odt * node[i].Jp[j];
			A[0] += node[i].dx;
			A[3] += node[i].dx;
			model_surface_flow_zero_inertia_implicit_invert(A, C);
			D[0] -= model->dt * node[i - 1].dF[0];
			D[2] -= model->dt * node[i - 1].dF[2];
			model_surface_flow_zero_inertia_implicit_multiply(C, D, node[i].dU);
//			node[i].U[0] = node[i].Un[0] + node[i].dU[0];
			node[i].U[2] = node[i].Un[2] + node[i].dU[2];
		}
		i = n1;
		model_surface_flow_zero_inertia_implicit_multiply(B, node[i].dU, D);
		model->outlet_contribution[0] += D[0];
		model->outlet_contribution[2] += D[2];

		for (i = 0; i < mesh->n; ++i)
		{
			DD[i] = node[i].dx;
			node[i].dU[0] *= node[i].dx;
		}
		for (i = 0; i < n1; ++i)
		{
			k = odt * fmin(node[i + 1].Jn[0], node[i].Jn[0]);
			CC[i] = EE[i] = k;
			DD[i] -= k;
			DD[i + 1] -= k;
		}
		for (i = 0; i < n1; ++i)
		{
			if (DD[i] == 0.) k = 0; else k = CC[i] / DD[i];
			DD[i + 1] -= k * EE[i];
			node[i + 1].dU[0] -= k * node[i].dU[0];
		}
		if (DD[i] == 0.) node[i].dU[0] = 0; else node[i].dU[0] /= DD[i];
		node[i].U[0] = node[i].Un[0] + node[i].dU[0];
		do
		{
			--i;
			if (DD[i] == 0.) node[i].dU[0] = 0;
			else node[i].dU[0] = (node[i].dU[0] - EE[i] * node[i+1].dU[0])
				/ DD[i];
			node[i].U[0] = node[i].Un[0] + node[i].dU[0];
		}
		while (i > 0);

		// boundary conditions

		model->model_inlet(model);
		for (j = 0; j < 4; ++j) A[j] = B[j] = odt * node[0].Jp[j];
		A[0] += node[0].dx;
		A[3] += node[0].dx;
		model_surface_flow_zero_inertia_implicit_invert(A, C);
		model_surface_flow_zero_inertia_implicit_multiply
			(C, model->inlet_contribution, node[0].dU);
		node[0].U[0] += node[0].dU[0];
		node[0].U[2] += node[0].dU[2];
		for (i = 0; ++i <= n1;)
		{
			model_surface_flow_zero_inertia_implicit_multiply
				(B, node[i - 1].dU, D);
			for (j = 0; j < 4; ++j) A[j] = B[j] = odt * node[i].Jp[j];
			A[0] += node[i].dx;
			A[3] += node[i].dx;
			model_surface_flow_zero_inertia_implicit_invert(A, C);
			model_surface_flow_zero_inertia_implicit_multiply(C, D, node[i].dU);
			node[i].U[0] += node[i].dU[0];
			node[i].U[2] += node[i].dU[2];
		}
		if (model->channel->type_inlet == 1) node_subcritical_discharge(node);
		model->model_outlet(model);
		i = n1;
		node[i].U[0] += model->outlet_contribution[0] / node[i].dx;
		node[i].U[2] += model->outlet_contribution[2] / node[i].dx;

		model_parameters(model);

		model_surface_flow_zero_inertia_stabilize(model);
	}
	while (++iteration < 1);
}
