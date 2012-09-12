/*

model_zero_inertia_upwind.c: source file to define the upwind numerical model
	applied to the zero-inertia model

Copyright 2011, Javier Burguete Tolosa.

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

#include <stdio.h>
#include <math.h>
#include "config.h"
#include "channel.h"
#include "node.h"
#include "mesh.h"
#include "model.h"
#include "model_zero_inertia.h"
#include "model_zero_inertia_upwind.h"

/**
 * \fn void model_surface_flow_zero_inertia_upwind(Model *model)
 * \brief Function to make the surface flow with the upwind numerical scheme.
 * \param model
 * \brief model struct.
 */

void model_surface_flow_zero_inertia_upwind(Model *model)
{
	int i, n1;
	double c, s, l1, l2, sA1, sA2, k1, k2;
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;
	double inlet_water_contribution, inlet_solute_contribution;
	inlet_water_contribution = model->dt * node[0].Q;
	inlet_solute_contribution = model->dt * node[0].T;
	n1 = mesh->n - 1;
	for (i = 0; i < n1; ++i)
	{
		model->node_flows(node + i);
		node[i].dQr = node[i].dFr = node[i].dTr = node[i].dQl = node[i].dFl
			= node[i].dTl = 0;
		if (node[i].h <= model->minimum_depth &&
			node[i + 1].h <= model->minimum_depth)
				continue;
		c = sqrt(G * (node[i + 1].A + node[i].A) / (node[i + 1].B + node[i].B));
		sA1 = sqrt(node[i].A);
		sA2 = sqrt(node[i + 1].A);
		k2 = sA1 + sA2;
		k1 = sA1 / k2;
		k2 = sA2 / k2;
		s = k1 * node[i].s + k2 * node[i + 1].s;
		l1 = c;
		l2 = -c;
		node[i].dQr = 0.5 * (l1 * node[i].dQ - node[i].dF) / c;
		node[i].dFr = l2 * node[i].dQr;
		node[i].dTr = s * node[i].dQr;
		node[i].dQl = node[i].dQ - node[i].dQr;
		node[i].dFl = node[i].dF - node[i].dFr;
		node[i].dTl = node[i].dT - node[i].dTr;
	}
	for (i = 0; i < n1; ++i)
	{
		node[i].A -= model->dt * node[i].dQr / node[i].dx;
		node[i].Q -= model->dt * node[i].dFr / node[i].dx;
		node[i].As -= model->dt * node[i].dTr / node[i].dx;
		node[i + 1].A -= model->dt * node[i].dQl / node[i + 1].dx;
		node[i + 1].Q -= model->dt * node[i].dFl / node[i + 1].dx;
		node[i + 1].As -= model->dt * node[i].dTl / node[i + 1].dx;
	}
	node[0].A -= inlet_water_contribution / node[0].dx;
	node[0].As -= inlet_solute_contribution / node[0].dx;
}

