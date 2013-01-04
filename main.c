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
 * \file main.c
 * \brief Main source code.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2012, Javier Burguete Tolosa.
 */
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "config.h"
#include "channel.h"
#include "node.h"
#include "mesh.h"
#include "model.h"
#include "model_hydrodynamic.h"
#include "model_zero_advection.h"
#include "model_zero_inertia.h"
#include "model_kinematic.h"
#include "model_hydrodynamic_upwind.h"
#include "model_zero_advection_upwind.h"
#include "model_zero_inertia_upwind.h"
#include "model_kinematic_upwind.h"
#include "model_hydrodynamic_LaxFriedrichs.h"
#include "model_zero_advection_LaxFriedrichs.h"
//#include "model_zero_inertia_LaxFriedrichs.h"
//#include "model_kinematic_LaxFriedrichs.h"
#include "model_hydrodynamic_implicit.h"

/**
 * \var critical_depth_tolerance
 * \brief Accuracy calculating the critical depth.
 */
double critical_depth_tolerance = 0.001;

/**
 * \fn int main(int argn, char **argc)
 * \brief Main function
 */
int main(int argn, char **argc)
{
	int i;
	FILE *file, *file_advance, *file_probes;
	clock_t cpu;
	Model model[1];
	if (argn < 3 || argn == 6 || argn > 7)
	{
		printf("the sintaxis is:\n./SWOCS input_file "
			"output_variables_file "
			"[output_flows_file] [output_advance_file]"
			"[input_probes_file output_probes_file]\n");
		return 1;
	}

	if (!model_read(model, argc[1])) return 2;

	switch (model->type_model)
	{
	case 1:
		model->model_node_parameters_centre
			= model->model_node_parameters_right
			= model->model_node_parameters_left
			= model_node_parameters_hydrodynamic;
		model->node_1dt_max = node_1dt_max_hydrodynamic;
		model->node_flows = node_flows_hydrodynamic;
		model->model_inlet_dtmax = model_inlet_dtmax_hydrodynamic;
		goto hydrodynamic;
	case 2:
		model->model_node_parameters_centre
			= model->model_node_parameters_right
			= model->model_node_parameters_left
			= model_node_parameters_zero_advection;
		model->node_1dt_max = node_1dt_max_zero_advection;
		model->node_flows = node_flows_zero_advection;
		model->model_inlet_dtmax = model_inlet_dtmax_zero_advection;
		goto zero_advection;
	case 3:
		switch (model->channel->friction_model)
		{
		case 1:
			model->node_discharge_centre
				= node_discharge_centre_zero_inertia_Manning;
			model->node_discharge_right
				= node_discharge_right_zero_inertia_Manning;
			model->node_discharge_left
				= node_discharge_left_zero_inertia_Manning;
			model->model_node_parameters_centre
				= model_node_parameters_centre_zero_inertia;
			model->model_node_parameters_right
				= model_node_parameters_right_zero_inertia;
			model->model_node_parameters_left
				= model_node_parameters_left_zero_inertia;
		}
		model->node_1dt_max = node_1dt_max_zero_inertia;
		model->node_flows = node_flows_zero_inertia;
		model->model_inlet_dtmax = model_inlet_dtmax_zero_inertia;
		goto zero_inertia;
	case 4:
		switch (model->channel->friction_model)
		{
		case 1:
			model->node_discharge_centre
				= node_discharge_centre_kinematic_Manning;
			model->node_discharge_right
				= node_discharge_right_kinematic_Manning;
			model->node_discharge_left
				= node_discharge_left_kinematic_Manning;
			model->model_node_parameters_centre
				= model_node_parameters_centre_kinematic;
			model->model_node_parameters_right
				= model_node_parameters_right_kinematic;
			model->model_node_parameters_left
				= model_node_parameters_left_kinematic;
		}
		model->node_1dt_max = node_1dt_max_kinematic;
		model->node_flows = node_flows_kinematic;
		model->model_inlet_dtmax = model_inlet_dtmax_kinematic;
		goto kinematic;
	default:
		printf("model: bad type (%u)\n", model->type_model);
		return 2;
	}

hydrodynamic:
	switch (model->type_surface_flow)
	{
	case 1:
		model->model_surface_flow = model_surface_flow_hydrodynamic_upwind;
		goto calculate;
	case 2:
		model->model_surface_flow =
			model_surface_flow_hydrodynamic_LaxFriedrichs;
		goto calculate;
	case 3:
		model->model_surface_flow = model_surface_flow_hydrodynamic_implicit;
		goto calculate;
	default:
		printf("model: bad surface flow type\n");
		return 2;
	}

zero_advection:
	switch (model->type_surface_flow)
	{
	case 1:
		model->model_surface_flow = model_surface_flow_zero_advection_upwind;
		goto calculate;
	case 2:
		model->model_surface_flow =
			model_surface_flow_zero_advection_LaxFriedrichs;
		goto calculate;
	default:
		printf("model: bad surface flow type\n");
		return 2;
	}

zero_inertia:
	switch (model->type_surface_flow)
	{
	case 1:
		model->model_surface_flow = model_surface_flow_zero_inertia_upwind;
		break;
//	case 2:
//		model->model_surface_flow = model_surface_flow_zero_inertia_LaxFriedrichs;
//		break;
	default:
		printf("model: bad surface flow type\n");
		return 2;
	}

kinematic:
	switch (model->type_surface_flow)
	{
	case 1:
		model->model_surface_flow = model_surface_flow_kinematic_upwind;
		break;
//	case 2:
//		model->model_surface_flow = model_surface_flow_kinematic_LaxFriedrichs;
//		break;
	default:
		printf("model: bad surface flow type\n");
		return 2;
	}

calculate:
	switch (model->type_diffusion)
	{
	case 1:
		model->model_diffusion = model_diffusion_explicit;
		break;
	case 2:
		model->model_diffusion = model_diffusion_implicit;
		break;
	default:
		printf("model: bad diffusion type\n");
		return 2;
	}

	if (argn > 4)
	{
		// opening the advance file
		file_advance = fopen(argc[4], "w");

		if (argn > 6)
		{
			// opening the probes files
			if (!model_probes_read(model, argc[5])) return 2;
			file_probes = fopen(argc[6], "w");
			if (!file_probes)
			{
				printf("model: unable to open the probes output file\n");
				return 2;
			}
		}
	}

	// reset the clock
	cpu = clock();

	// init model parameters
	model_parameters(model);

	// main calculation bucle
	for (model->t = 0, i = 0; model->t < model->tfinal; ++i)
	{
		if (argn > 4)
		{
			// writing the advance
			model_write_advance(model, file_advance);

			// writing the probes
			if (argn > 6) model_write_probes(model, file_probes);
		}

		// model step
		model_step(model);
//		model_print(model, i);
	}

	// printing main results
	printf("cpu=%lg ", (clock() - cpu) / ((double)CLOCKS_PER_SEC));
	model_print(model, i);

	// writing result variables
	file = fopen(argc[2], "w");
	mesh_write_variables(model->mesh, file);
	fclose(file);

	// writing result flows
	if (argn > 3)
	{
		file = fopen(argc[3], "w");
		mesh_write_flows(model->mesh, file);
		fclose(file);

		if (argn > 4)
		{
			// closing the advance
			model_write_advance(model, file_advance);
			fclose(file_advance);

			if (argn > 6)
			{
				// closing the probes
				model_write_probes(model, file_probes);
				fclose(file_probes);
			}
		}
	}

	return 0;
}
