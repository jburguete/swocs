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
 * \file model.c
 * \brief Source file to define the numerical model.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2012, Javier Burguete Tolosa.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "config.h"
#include "channel.h"
#include "node.h"
#include "mesh.h"
#include "model.h"

/**
 * \fn void model_parameters(Model *model)
 * \brief Function to calculate the model parameters.
 * \param model
 * \brief model struct.
 */
void model_parameters(Model *model)
{
	int i, n1;
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;
	model->model_node_parameters_right(model, node);
	n1 = mesh->n - 1;
	for (i = 0; ++i < n1;)
		model->model_node_parameters_centre(model, node + i);
	model->model_node_parameters_left(model, node + i);
}

/**
 * \fn void model_infiltration(Model *model)
 * \brief Function to make the infiltration model.
 * \param model
 * \brief model struct.
 */
void model_infiltration(Model *model)
{
	int i;
	double Pidt;
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;
	for (i = 0; i < mesh->n; ++i)
	{
		Pidt = fmin(node[i].Pi * model->dt, node[i].U[0]);
		node[i].U[0] -= Pidt;
		node[i].U[3] += Pidt;
		Pidt *= node[i].s;
		node[i].U[2] -= Pidt;
		node[i].U[4] += Pidt;
	}
}

/**
 * \fn void model_diffusion_explicit(Model *model)
 * \brief Function to make the explicit diffusion model.
 * \param model
 * \brief model struct.
 */
void model_diffusion_explicit(Model *model)
{
	int i, n1;
	double dD;
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;
	n1 = mesh->n - 1;
	for (i = 0; i < n1; ++i)
	{
		dD = model->dt * fmin(node[i + 1].KxA, node[i].KxA)
			* (node[i + 1].s - node[i].s) / node[i].ix;
		node[i].U[2] += dD / node[i].dx;
		node[i + 1].U[2] -= dD / node[i + 1].dx;
	}
}

/**
 * \fn void model_diffusion_implicit(Model *model)
 * \brief Function to make the implicit diffusion model.
 * \param model
 * \brief model struct.
 */
void model_diffusion_implicit(Model *model)
{
	int i, n1;
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;
	double k, C[mesh->n], D[mesh->n], E[mesh->n], H[mesh->n];
	for (i = 0; i < mesh->n; ++i)
	{
		D[i] = node[i].U[0] * node[i].dx;
		H[i] = node[i].U[2] * node[i].dx;
	}
	n1 = mesh->n - 1;
	for (i = 0; i < n1; ++i)
	{
		k = model->dt * fmin(node[i + 1].KxA, node[i].KxA) / node[i].ix;
		C[i] = E[i] = -k;
		D[i] += k;
		D[i + 1] += k;
	}
	for (i = 0; i < n1; ++i)
	{
		if (D[i] == 0.) k = 0; else k = C[i] / D[i];
		D[i + 1] -= k * E[i];
		H[i + 1] -= k * H[i];
	}
	if (D[i] == 0.) H[i] = 0; else H[i] /= D[i];
	node[i].U[2] = H[i] * node[i].U[0];
	while (--i >= 0)
	{
		if (D[i] == 0.) H[i] = 0; else H[i] = (H[i] - E[i] * H[i+1]) / D[i];
		node[i].U[2] = H[i] * node[i].U[0];
	}
}

/**
 * \fn double model_node_diffusion_1dt_max(Node *node)
 * \brief Function to calculate the allowed maximum time step size in a node
 *   with the diffusion model.
 * \param node
 * \brief node struct.
 * \return inverse of the allowed maximum time step size.
 */
double model_node_diffusion_1dt_max(Node *node)
{
	return (2 * node->Kx + fabs(node->u) * node->dx) / (node->dx * node->dx);
}

/**
 * \fn void model_step(Model *model)
 * \brief Function to make a step of the numerical model.
 * \param model
 * \brief model struct.
 */
void model_step(Model *model)
{
	int i;
	double dtmax;
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;
	for (i = 0, dtmax = 0; i < mesh->n; ++i)
		dtmax = fmax(dtmax, model->node_1dt_max(node + i));
	if (model->type_diffusion == 1)
	{
		for (i = 0, dtmax = 0; i < mesh->n; ++i)
			dtmax = fmax(dtmax, model_node_diffusion_1dt_max(node + i));
	}
	dtmax = model->cfl / dtmax;
	dtmax = fmin(dtmax, model->model_inlet_dtmax(model));
	model->t2 = fmin(model->tfinal, model->t + dtmax);
	model->dt = model->t2 - model->t;
	model->model_surface_flow(model);
	model->model_diffusion(model);
	model_infiltration(model);
	model_parameters(model);
	model->t = model->t2;
}

/**
 * \fn int model_read(Model *model, char *file_name)
 * \brief Function to read the numerical model.
 * \param model
 * \brief model struct.
 * \param file_name
 * \brief name of the input data file
 * \return 0 on error, 1 on success.
 */
int model_read(Model *model, char *file_name)
{
	char *msg;
	FILE *file;

	file = fopen(file_name, "r");
	if (!file)
	{
		msg = "model: unable to open the input file\n";
		goto bad;
	}

	if (!channel_read(model->channel, file)) goto bad;
	switch (model->channel->friction_model)
	{
	case 1:
		model->node_friction = node_friction_Manning;
	}
	switch (model->channel->infiltration_model)
	{
	case 1:
		model->node_infiltration = node_infiltration_KostiakovLewis;
	}
	switch (model->channel->diffusion_model)
	{
	case 1:
		model->node_diffusion = node_diffusion_Rutherford;
	}

	if (!mesh_read(model->mesh, model->channel, file)) goto bad;

	model->model_inlet = model_inlet;
	switch (model->channel->type_outlet)
	{
	case 1:
		model->model_outlet = model_outlet_closed;
		break;
	case 2:
		model->model_outlet = model_outlet_open;
	}

	if (fscanf(file, "%lf%lf%lf%lf%d%d%d",
		&model->tfinal,
		&model->interval,
		&model->cfl,
		&model->minimum_depth,
		&model->type_surface_flow,
		&model->type_diffusion,
		&model->type_model) != 7)
	{
		msg = "model: bad data\n";
		goto bad;
	}
#if DEBUG_MODEL_READ
	printf("model:\n"
		"tfinal=%lf interval=%lf cfl=%lf\n"
		"type_surface_flow=%d type_model=%d\n",
		model->tfinal,
		model->interval,
		model->cfl,
		model->type_surface_flow,
		model->type_model);
#endif

	fclose(file);
	return 1;

bad:
	if (file) fclose(file);
	printf(msg);
	return 0;
}

/**
 * \fn void model_print(Model *model, int nsteps)
 * \brief Function to print a model stat.
 * \param model
 * \brief model struct.
 * \param nsteps
 * \brief number of time steps.
 */
void model_print(Model *model, int nsteps)
{
	printf(
		"main: steps number=%d t=%.14lg water mass=%.14lg solute mass=%.14lg\n",
		nsteps,
		model->t,
		mesh_water_mass(model->mesh),
		mesh_solute_mass(model->mesh));
}

/**
 * \fn void model_write_advance(Model *model, FILE *file)
 * \brief Function to write in a file the channel water advance.
 * \param model
 * \brief model struct.
 * \param file
 * \brief output file.
 */
void model_write_advance(Model *model, FILE *file)
{
	int i;
	Mesh *mesh = model->mesh;
	Node *node;
	for (i = 0; i < mesh->n; ++i)
	{
		node = mesh->node + i;
		if (node->U[0] == 0) break;
	}
	if (i) --i;
	fprintf(file, "%lg %lg\n", model->t, mesh->node[i].x);
}

/**
 * \fn int model_probes_read(Model *model, char *name)
 * \brief Function to read the model probes in a file.
 * \param model
 * \brief model struct.
 * \param name
 * \brief input file name.
 * \return 0 on error, 1 on success.
 */
int model_probes_read(Model *model, char *name)
{
	int i, j, k;
	double d, dmin, *x;
	char *msg;
	FILE *file;
	Probes *probes = model->probes;
	Node *node = model->mesh->node;
	file = fopen(name, "r");
	if (!file)
	{
		msg = "probes: unable to open the input file\n";
		goto bad2;
	}
	if (fscanf(file, "%d", &probes->n) != 1) goto bad;
	probes->x = x = (double*)malloc(probes->n * sizeof(double));
	probes->node = (int*)malloc(probes->n * sizeof(int));
	for (i = 0; i < probes->n; ++i)
	{
		if (fscanf(file, "%lf", x + i) != 1) goto bad;
		k = 0;
		dmin = fabs(x[i] - node[0].x);
		for (j = 0; ++j < model->mesh->n;)
		{
			d = fabs(x[i] - node[j].x);
			if (d < dmin)
			{
				dmin = d;
				k = j;
			}
		}
		probes->node[i] = k;
	}
	return 1;

bad:
	msg = "probes: bad data\n";
	fclose(file);

bad2:
	printf(msg);
	return 0;
}

/**
 * \fn void model_write_probes(Model *model, FILE *file)
 * \brief Function to write the model probes in a file.
 * \param model
 * \brief model struct.
 * \param file
 * \brief output file.
 */
void model_write_probes(Model *model, FILE *file)
{
int i;
	Probes *probes = model->probes;
	Node *node;

	// writing the time
	fprintf(file, "%lg ", model->t);
	for (i = 0; i < probes->n; ++i)
	{
		//writing the depth and the concentration of the i-th probe
		node = model->mesh->node + probes->node[i];
		fprintf(file, "%lg %lg ", node->h, node->s);
	}
	// writing a new row
	fprintf(file, "\n");
}

/**
 * \fn void model_inlet(Model *model)
 * \brief Function to calculate the inlet boundary condition.
 * \param model
 * \brief model struct.
 */
void model_inlet(Model *model)
{
	double t, t2;
	t = model->t;
	t2 = model->t2;
	model->inlet_contribution[0] +=
		hydrogram_integrate(model->channel->water_inlet, t, t2);
	model->inlet_contribution[2] +=
		hydrogram_integrate(model->channel->solute_inlet, t, t2);
}

/**
 * \fn void model_outlet_closed(Model *model)
 * \brief Function to calculate a closed outlet boundary condition.
 * \param model
 * \brief model struct.
*/
void model_outlet_closed(Model *model)
{
	Node *node = model->mesh->node + model->mesh->n - 1;
	node->U[1] = 0.;
}

/**
 * \fn void model_outlet_open(Model *model)
 * \brief Function to calculate a open outlet boundary condition.
 * \param model
 * \brief model struct.
 */
void model_outlet_open(Model *model)
{
	Node *node = model->mesh->node + model->mesh->n - 1;
	node_depth(node);
	node_width(node);
	node_critical_velocity(node);
	node->U[1] = fmax(node->U[1], node->U[0] * node->c);
}
