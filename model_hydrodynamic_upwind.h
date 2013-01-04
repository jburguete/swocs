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
 * \file model_hydrodynamic_upwind.h
 * \brief Header file to define the first order upwind explicit numerical model 
 *   applied to the hydrodynamic model.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2012, Javier Burguete Tolosa.
 */

// in order to prevent multiple definitions

#ifndef MODEL_HYDRODYNAMIC_UPWIND__H
#define MODEL_HYDRODYNAMIC_UPWIND__H 1

// member functions

void model_surface_flow_hydrodynamic_upwind_inlet(Model *model);
void model_surface_flow_hydrodynamic_upwind(Model *model);

#endif
