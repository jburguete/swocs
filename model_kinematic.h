/*
SWOCS: a software to check the numerical performance of different models in
	channel or furrow flows

Copyright 2011-2014, Javier Burguete Tolosa.

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
 * \file model_kinematic.h
 * \brief Header file to define the kinematic model.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2014, Javier Burguete Tolosa.
 */

// in order to prevent multiple definitions
#ifndef MODEL_KINEMATIC__H
#define MODEL_KINEMATIC__H 1

// member functions

void node_discharge_centre_kinematic(Node *node);
void node_discharge_right_kinematic(Node *node);
void node_discharge_left_kinematic(Node *node);
void model_node_parameters_kinematic(Model *model, Node *node,
	void (*node_discharge)(Node*));
void model_node_parameters_centre_kinematic(Model *model, Node *node);
void model_node_parameters_right_kinematic(Model *model, Node *node);
void model_node_parameters_left_kinematic(Model *model, Node *node);
double node_1dt_max_kinematic(Node *node);
void node_flows_kinematic(Node *node1);
double model_inlet_dtmax_kinematic(Model *model);

#endif
