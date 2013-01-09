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
 * \file model_hydrodynamic_LaxFriedrichs.c
 * \brief Source file to define the Lax-Friedrichs numerical model applied to
 *   the hydrodynamic model.
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
#include "model_hydrodynamic_LaxFriedrichs.h"

/**
 * \fn void model_surface_flow_hydrodynamic_LaxFriedrichs(Model *model)
 * \brief Function to make the surface flow with the Lax-Friedrichs numerical
 *   scheme.
 * \param model
 * \brief model struct.
 */
void model_surface_flow_hydrodynamic_LaxFriedrichs(Model *model)
{
	unsigned int i, j, n1;
	double k1, k2;
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;

	model->inlet_contribution[0] = - model->dt * node[0].U[1];
	model->inlet_contribution[2] = - model->dt * node[0].T;

	n1 = mesh->n - 1;
	for (i = 0; i < n1; ++i)
	{
		node_flows_hydrodynamic(node + i);
		node[i].dFr[0] = node[i].dFr[1] = node[i].dFr[2] = node[i].dFl[0]
			= node[i].dFl[1] = node[i].dFl[2] = 0;
		if (node[i].h <= model->minimum_depth &&
			node[i + 1].h <= model->minimum_depth)
				continue;

		// artificial viscosity

		k1 = 0.5 * fmax(node[i + 1].c + fabs(node[i + 1].u),
			node[i].c + fabs(node[i].u));

		// wave decomposition

		for (j = 0; j < 3; ++j)
		{
			node[i].dFr[j] = node[i].dFl[j] = 0.5 * node[i].dF[j];
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

	// boundary correction

	model->model_inlet(model);
	node[0].U[0] += model->inlet_contribution[0] / node[0].dx;
	node[0].U[2] += model->inlet_contribution[2] / node[0].dx;
	if (model->channel->type_inlet == 1) node_subcritical_discharge(node);
	model->model_outlet(model);
}
