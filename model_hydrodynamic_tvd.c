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
 * \file model_hydrodynamic_upwind.c
 * \brief Source file to define the high order TVD explicit numerical model 
 *   applied to the hydrodynamic model.
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
#include "model_hydrodynamic.h"
#include "model_hydrodynamic_tvd.h"

/**
 * \fn void model_surface_flow_hydrodynamic_limiter(double dW1, double dW2)
 * \brief Function to make the high order TVD flux limiter.
 * \param dW1
 * \brief centred high order flux.
 * \param dW2
 * \brief upwind high order flux.
 * \return flux limiter value.
 */
double model_surface_flow_hydrodynamic_limiter(double dW1, double dW2)
{
	double r;
	if (dW1 * dW2 <= 0.) return 0.;
	r = dW1 / dW2;
	return fmax(fmin(2. * r, 1.), fmin(r, 2.));
}

/**
 * \fn void model_surface_flow_hydrodynamic_tvd(Model *model)
 * \brief Function to make the surface flow with the high order TVD numerical
 *   scheme.
 * \param model
 * \brief model struct.
 */
void model_surface_flow_hydrodynamic_tvd(Model *model)
{
	int i, j, n1;
	double c, u, s, l1, l2, sA1, sA2, k1, k2, dt2, lp[3], ln[3];
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;

	model->inlet_contribution[0] = - model->dt * node[0].U[1];
	model->inlet_contribution[2] = - model->dt * node[0].T;

	n1 = mesh->n - 1;
	for (i = 0; i < n1; ++i)
	{
		model->node_flows(node + i);
		node[i].dFr[0] = node[i].dFr[1] = node[i].dFr[2] = node[i].dFl[0]
			= node[i].dFl[1] = node[i].dFl[2] = node[i].dWr[0] = node[i].dWr[1]
		   	= node[i].dWr[2] = node[i].dWl[0] = node[i].dWl[1] = node[i].dWl[2]
			= 0;
		if (node[i].h <= model->minimum_depth &&
			node[i + 1].h <= model->minimum_depth)
				continue;

		// Roe's averages

		c = sqrt(G * (node[i + 1].U[0] + node[i].U[0]) /
			(node[i + 1].B + node[i].B));
		sA1 = sqrt(node[i].U[0]);
		sA2 = sqrt(node[i + 1].U[0]);
		k2 = sA1 + sA2;
		k1 = sA1 / k2;
		k2 = sA2 / k2;
		u = k1 * node[i].u + k2 * node[i + 1].u;
		l1 = u + c;
		l2 = u - c;
		s = k1 * node[i].s + k2 * node[i + 1].s;
		node[i].l[0] = l1;
		node[i].l[1] = l2;
		node[i].l[2] = s;

		// first order wave decomposition

		if (u >= c)
		{
			for (j = 0; j < 3; ++j) node[i].dFl[j] = node[i].dF[j];
		}
		else
		{
			node[i].dFr[0] = 0.5 * (l1 * node[i].dF[0] - node[i].dF[1]) / c;
			node[i].dFr[1] = l2 * node[i].dFr[0];
			node[i].dFr[2] = s * node[i].dFr[0];
			for (j = 0; j < 3; ++j)
				node[i].dFl[j] = node[i].dF[j] - node[i].dFr[j];
		}

		// high order TVD correction

		lp[0] = fmax(0., l1);
		lp[1] = fmax(0., l2);
		lp[2] = fmax(0., u);
		ln[0] = fmin(0., l1);
		ln[1] = fmin(0., l2);
		ln[2] = fmin(0., u);
		node[i].dWl[0] = 0.5 * (1. - lp[0] * model->dt / node[i].ix)
			* (node[i].dFl[1] - l2 * node[i].dFl[0]) / c;
		node[i].dWl[1] = 0.5 * (1. - lp[1] * model->dt / node[i].ix)
			* (l1 * node[i].dFl[0] - node[i].dFl[1]) / c;
		node[i].dWl[2] = (1. - lp[1] * model->dt / node[i].ix)
			* (node[i].dFl[2] - s * node[i].dFl[0]);
		node[i].dWr[0] = 0.5 * (1. + ln[0] * model->dt / node[i].ix)
			* (node[i].dFr[1] - l2 * node[i].dFr[0]) / c;
		node[i].dWr[1] = 0.5 * (1. + ln[1] * model->dt / node[i].ix)
			* (l1 * node[i].dFr[0] - node[i].dFr[1]) / c;
		node[i].dWr[2] = (1. + ln[1] * model->dt / node[i].ix)
			* (node[i].dFr[2] - s * node[i].dFr[0]);

		// entropy correction

		if (node[i].l2 < 0. && node[i + 1].l2 > 0.)
			k1 = 0.25 * (node[i + 1].l2 - node[i].l2 - 2 * fabs(l2));
		else if (node[i].l1 < 0. && node[i + 1].l1 > 0.)
			k1 = 0.25 * (node[i + 1].l1 - node[i].l1 - 2 * fabs(l1));
		else continue;
		for (j = 0; j < 3; ++j)
		{
			k2 = k1 * (node[i + 1].U[j] - node[i].U[j]);
			node[i].dFl[j] += k2;
			node[i].dFr[j] -= k2;
		}
	}

	// variables updating

	for (i = 0; i < n1; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			node[i].U[j] -= model->dt * node[i].dFr[j] / node[i].dx;
			node[i + 1].U[j] -= model->dt * node[i].dFl[j] / node[i + 1].dx;

		}
	}

	// high order variables correction

	dt2 = 0.5 * model->dt;
	for (i = 1; i < n1 - 1; ++i)
	{
		for (j = 0; j < 3; ++j)
		{
			lp[j] = dt2 * node[i - 1].dWl[j]
				* model_surface_flow_hydrodynamic_limiter
					(node[i].dWl[j], node[i - 1].dWl[j]);
			ln[j] = dt2 * node[i].dWr[j]
				* model_surface_flow_hydrodynamic_limiter
					(node[i - 1].dWr[j], node[i].dWr[j]);
		}
		k1 = lp[0] + lp[1];
		k2 = ln[0] + ln[1];
		node[i + 1].U[0] += k1 / node[i + 1].dx;
		node[i - 1].U[0] += k2 / node[i - 1].dx;
		node[i].U[0] -= (k1 + k2) / node[i].dx;
		k1 = k1 * node[i].l[2] + lp[2];
		k2 = k2 * node[i].l[2] + ln[2];
		node[i + 1].U[2] += k1 / node[i + 1].dx;
		node[i - 1].U[2] += k2 / node[i - 1].dx;
		node[i].U[2] -= (k1 + k2) / node[i].dx;
		k1 = node[i].l[0] * lp[0] + node[i].l[1] * lp[1];
		k2 = node[i].l[0] * ln[0] + node[i].l[1] * ln[1];
		node[i + 1].U[1] += k1 / node[i + 1].dx;
		node[i - 1].U[1] += k2 / node[i - 1].dx;
		node[i].U[1] -= (k1 + k2) / node[i].dx;
	}

	// boundary correction

	model->model_inlet(model);
	node[0].U[0] += model->inlet_contribution[0] / node[0].dx;
	node[0].U[2] += model->inlet_contribution[2] / node[0].dx;
	node_subcritical_discharge(node);
	model->model_outlet(model);
}
