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
 * \file channel.h
 * \brief Header file to define a channel.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011, Javier Burguete Tolosa.
 */

// in order to prevent multiple definitions
#ifndef CHANNEL__H
#define CHANNEL__H 1

/**
 * \struct _Hydrogram
 * \brief Struct to define a hydrogram.
 */
struct _Hydrogram
{
/**
 * \var t
 * \brief array of times.
 * \var Q
 * \brief array of discharges.
 * \var n
 * \brief number of points defining the hydrogram.
 */
	double *t, *Q;
	int n;
};

/**
 * \typedef Hydrogram
 */
typedef struct _Hydrogram Hydrogram;

/**
 * \struct _Channel
 * \brief Struct to define a channel.
 */
struct _Channel
{
/**
 * \var water_inlet
 * \brief hydrogram of water inlet.
 * \var solute_inlet
 * \brief hydrogram of solute inlet.
 * \var friction_coefficient
 * \brief array of friction coefficients.
 * \var infiltration_coefficient
 * \brief array of infiltration coefficients.
 * \var diffusion_coefficient
 * \brief array of diffusion coefficients.
 * \var length
 * \brief channel length.
 * \var slope
 * \brief channel slope.
 * \var bottom_width
 * \brief bottom width.
 * \var wall_slope
 * \brief slope of the lateral walls.
 * \var height
 * \brief channel height.
 * \var type_outlet
 * \brief type of outlet (1 closed, 2 open).
 * \var friction_model
 * \brief type of friction model (1 Gauckler-Manning).
 * \var infiltration_model
 * \brief type of infiltration model (1 Kostiakov-Lewis).
 * \var diffusion_model
 * \brief type of diffusion model (1 Rutherford).
 */
	Hydrogram water_inlet[1], solute_inlet[1]; 
	double friction_coefficient[3], infiltration_coefficient[4],
		diffusion_coefficient[1], slope, length, bottom_width, wall_slope,
		height;
	int type_outlet, friction_model, infiltration_model, diffusion_model;
};

/**
 * \typedef Channel
 */
typedef struct _Channel Channel;

// member functions

double interpolate(double x, double x1, double x2, double y1, double y2);

int hydrogram_read(Hydrogram *hydrogram, FILE *file);
double hydrogram_discharge(Hydrogram *hydrogram, double t);
double hydrogram_integrate(Hydrogram *hydrogram, double t1, double t2);

int channel_friction_read_Manning(Channel *channel, FILE *file);
int channel_infiltration_read_KostiakovLewis(Channel *channel, FILE *file);
int channel_diffusion_read_Rutherford(Channel *channel, FILE *file);
int channel_read(Channel *channel,FILE *file);

#endif