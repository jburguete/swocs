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
 * \file model_diffusive_upwind.c
 * \brief Source file to define the upwind numerical model applied to the
 *   diffusive model.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011, Javier Burguete Tolosa.
 */
#include <stdio.h>
#include <math.h>
#include "config.h"
#include "channel.h"
#include "node.h"
#include "mesh.h"
#include "model.h"
#include "model_diffusive_upwind.h"

/**
 * \fn void model_surface_flow_diffusive_upwind(Model *model)
 * \brief Function to make the surface flow with the upwind numerical scheme.
 * \param model
 * \brief model struct.
 */
void model_surface_flow_diffusive_upwind(Model *model)
{
	int i;
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;
	double inlet_water_contribution, inlet_solute_contribution;
	inlet_water_contribution = model->dt * node[0].Q;
	inlet_solute_contribution = model->dt * node[0].T;
	for (i = 0; ++i < mesh->n;)
	{
		model->node_flows(node + i - 1);
		node[i].A -= model->dt * node[i - 1].dQ / node[i].dx;
		node[i].As -= model->dt * node[i - 1].dT / node[i].dx;
	}
	node[0].A -= inlet_water_contribution / node[0].dx;
	node[0].As -= inlet_solute_contribution / node[0].dx;
}
