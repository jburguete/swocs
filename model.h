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
 * \file model.h
 * \brief Header file to define the numerical model.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011, Javier Burguete Tolosa.
 */

// in order to prevent multiple definitions
#ifndef MODEL__H
#define MODEL__H 1

/**
 * \struct _Probes
 * \brief Struct to define probes to save the evolution of the variables at a
 *   mesh cell.
 */
struct _Probes
{
/**
 * \var x
 * \brief array of x-coordinates of the probes.
 * \var node
 * \brief array of positions of the probes in the mesh.
 * \var n
 * \brief number of probes.
 */
	double *x;
	int *node, n;
};

/**
 * \typedef Probes
 */
typedef struct _Probes Probes;

/**
 * \struct _Model
 * \brief Struct to define a numerical model.
 */
struct _Model
{
/**
 * \var mesh
 * \brief mesh struct.
 * \var channel
 * \brief channel struct.
 * \var probes
 * \brief probes struct.
 * \var t
 * \brief actual time.
 * \var t2 
 * \brief next time.
 * \var dt
 * \brief time step size.
 * \var tfinal
 * \brief final time.
 * \var cfl
 * \brief CFL number.
 * \var interval
 * \brief time interval to save the data.
 * \var minimum_depth
 * \brief minimum depth allowing the water movement.
 * \var model_node_parameters_centre
 * \brief pointer to the function calculating the node parameters in a centred
 *   form.
 * \var model_node_parameters_right
 * \brief pointer to the function calculating the node parameters in an right
 *   upwind form.
 * \var model_node_parameters_left
 * \brief pointer to the function calculating the node parameters in an left
 *   upwind form.
 * \var node_1dt_max
 * \brief pointer to the function calculating the maximum allowed time at a
 *   node.
 * \var model_inlet_dtmax
 * \brief pointer to the function calculating the maximum allowed time at the
 *   inlet
 * \var node_flows
 * \brief pointer to the function calculating the node flows.
 * \var node_discharge_centre
 * \brief pointer to the function calculating the node discharge in a centred
 *   form.
 * \var node_discharge_right
 * \brief pointer to the function calculating the node discharge in an right
 *   upwind form.
 * \var node_discharge_left
 * \brief pointer to the function calculating the node discharge in an left
 *   upwind form.
 * \var node_friction
 * \brief pointer to the function calculating the node friction.
 * \var node_infiltration
 * \brief pointer to the function calculating the node infiltration.
 * \var node_diffusion
 * \brief pointer to the function calculating the node diffusion.
 * \var node_inlet
 * \brief pointer to the function calculating the inlet.
 * \var node_outlet
 * \brief pointer to the function calculating the outlet.
 * \var model_surface_flow
 * \brief pointer to the function defining the numerical surface flow scheme.
 * \var model_diffusion
 * \brief pointer to the function defining the numerical diffusion scheme.
 * \var type_surface_flow
 * \brief type of numerical surface flow scheme (1 McCormack, 2 upwind).
 * \var type_diffusion
 * \brief type of numerical diffusion scheme (1 explicit, 2 implicit).
 * \var type_model
 * \brief type of model (1 complete, 2 zero-inertia, 3 diffusive, 4 kinematic).
*/
	Mesh mesh[1];
	Channel channel[1];
	Probes probes[1];
	double t, t2, dt, tfinal, cfl, interval, minimum_depth;
	void (*model_node_parameters_centre)(struct _Model *model, Node *node);
	void (*model_node_parameters_right)(struct _Model *model, Node *node);
	void (*model_node_parameters_left)(struct _Model *model, Node *node);
	double (*node_1dt_max)(Node *node);
	double (*model_inlet_dtmax)(struct _Model *model);
	void (*node_flows)(Node *node1);
	void (*node_discharge_centre)(Node *node);
	void (*node_discharge_right)(Node *node);
	void (*node_discharge_left)(Node *node);
	void (*node_friction)(Node *node);
	void (*node_infiltration)(Node *node);
	void (*node_diffusion)(Node *node);
	void (*node_inlet)
		(Node *node, Hydrogram *water, Hydrogram *solute, double t, double t2);
	void (*node_outlet)(Node *node);
	void (*model_surface_flow)(struct _Model *model);
	void (*model_diffusion)(struct _Model *model);
	int type_surface_flow, type_diffusion, type_model;
};

/**
 * \typedef Model
 */
typedef struct _Model Model;

// member functions

void model_parameters(Model *model);
void model_infiltration(Model *model);
void model_diffusion_explicit(Model *model);
void model_diffusion_implicit(Model *model);
double model_node_diffusion_1dt_max(Node *node);
void model_step(Model *model);
int model_read(Model *model, char *file_name);
void model_print(Model *model, int nsteps);
void model_write_advance(Model *model, FILE *file);
int model_probes_read(Model *model, char *name);
void model_write_probes(Model *model, FILE *file);

#endif
