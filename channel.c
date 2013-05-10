/*
SWOCS: a software to check the numerical performance of different models in
	channel or furrow flows

Copyright 2011-2013, Javier Burguete Tolosa.

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
 * \file channel.c
 * \brief Source file to define a channel.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2013, Javier Burguete Tolosa.
 */
#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "channel.h"

/**
 * \define DEBUG_CHANNEL
 * \brief Macro to debug the channel struct functions.
 */
#define DEBUG_CHANNEL 0

/**
 * void print_error(char *msg)
 * \brief Function to print an error message.
 * \param msg
 * \brief error message.
 */
void print_error(char *msg)
{
	printf("ERROR!\n%s\n", msg);
}

/**
 * \fn double interpolate(double x, double x1, double x2, double y1, double y2)
 * \brief Function to calculate an interpolation.
 * \param x
 * \brief x-coordinate of the interpolation point.
 * \param x1
 * \brief x-coordinate of the first point.
 * \param x2
 * \brief x-coordinate of the second point.
 * \param y1
 * \brief y-coordinate of the first point.
 * \param y2
 * \brief y-coordinate of the second point.
 * \return y-coordinate of the interpolation point.
 */
double interpolate(double x, double x1, double x2, double y1, double y2)
{
	return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

/**
 * \fn int hydrogram_read(Hydrogram *hydrogram, FILE *file)
 * \brief Function to read the data of a hydrogram.
 * \param hydrogram
 * \brief hydrogram struct.
 * \param file
 * \brief input file.
 * \return 0 on error, 1 on success.
 */
int hydrogram_read(Hydrogram *hydrogram, FILE *file)
{
	unsigned int i;
	char *msg;
	if (fscanf(file, "%u", &hydrogram->n) != 1 || hydrogram->n < 1)
	{
		msg = "hydrogram: bad points number";
		goto bad;
	}
#if DEBUG_CHANNEL
	printf("hydrogram: n=%u\n", hydrogram->n);
#endif
	hydrogram->t = (double*)malloc(hydrogram->n * sizeof(double));
	hydrogram->Q = (double*)malloc(hydrogram->n * sizeof(double));
	if (!hydrogram->t || !hydrogram->Q)
	{
		msg = "hydrogram: not enough memory";
		goto bad;
	}
	for (i = 0; i < hydrogram->n; ++i)
	{
		if (fscanf(file, "%lf%lf", hydrogram->t + i, hydrogram->Q + i) != 2)
		{
			msg = "hydrogram: bad defined";
			goto bad;
		}
#if DEBUG_CHANNEL
		printf("hydrogram: t=%lg Q=%lg\n", hydrogram->t[i], hydrogram->Q[i]);
#endif
	}
	return 1;

bad:
	print_error(msg);
	return 0;
}

/**
 * \fn double hydrogram_discharge(Hydrogram *hydrogram, double t)
 * \brief Function to calculate the discharge in a hydrogram.
 * \param hydrogram
 * \brief hydrogram struct.
 * \param t
 * \brief time.
 * \return discharge.
 */
double hydrogram_discharge(Hydrogram *hydrogram, double t)
{
	unsigned int i, n1;
	n1 = hydrogram->n - 1;
	if (n1 == 0 || t <= hydrogram->t[0]) return hydrogram->Q[0];
	if (t >= hydrogram->t[n1]) return hydrogram->Q[n1];
	for (i = 0; t > hydrogram->t[i];) ++i;
	return interpolate(t, hydrogram->t[i], hydrogram->t[i - 1],
		hydrogram->Q[i], hydrogram->Q[i - 1]);
}

/**
 * \fn double hydrogram_integrate(Hydrogram *hydrogram, double t1, double t2)
 * \brief Function to integrate the mass flux in a hydrogram.
 * \param hydrogram
 * \brief hydrogram struct.
 * \param t1
 * \brief initial time.
 * \param t2
 * \brief final time.
 * \return integral of the mass flux.
 */
double hydrogram_integrate(Hydrogram *hydrogram, double t1, double t2)
{
	unsigned int i, j, n1;
	double Q1, Q2, I;
	n1 = hydrogram->n - 1;
	if (n1 == 0 || t2 <= hydrogram->t[0]) return hydrogram->Q[0] * (t2 - t1);
	if (t1 >= hydrogram->t[n1]) return hydrogram->Q[n1] * (t2 - t1);
	for (i = 0; t1 > hydrogram->t[i];) ++i;
	for (j = i; j < hydrogram->n && t2 > hydrogram->t[j];) ++j;
	if (i == j)
	{
		Q1 = interpolate(t1, hydrogram->t[i], hydrogram->t[i - 1],
			hydrogram->Q[i], hydrogram->Q[i - 1]);
		Q2 = interpolate(t2, hydrogram->t[i], hydrogram->t[i - 1],
			hydrogram->Q[i], hydrogram->Q[i - 1]);
		return 0.5 * (Q1 + Q2) * (t2 - t1);
	}
	if (i == 0)
	{
		Q1 = hydrogram->Q[0];
	}
	else
	{
		Q1 = interpolate(t1, hydrogram->t[i], hydrogram->t[i - 1],
			hydrogram->Q[i], hydrogram->Q[i - 1]);
	}
	I = 0.5 * (Q1 + hydrogram->Q[i]) * (hydrogram->t[i] - t1);
	while (++i < j)
	{
		I += 0.5 * (hydrogram->Q[i] - hydrogram->Q[i - 1])
			* (hydrogram->t[i] - hydrogram->t[i - 1]);
	}
	if (i == hydrogram->n)
	{
		Q2 = hydrogram->Q[n1];
	}
	else
	{
		Q2 = interpolate(t2, hydrogram->t[i], hydrogram->t[i - 1],
			hydrogram->Q[i], hydrogram->Q[i - 1]);
	}
	return I + 0.5 * (Q2 + hydrogram->Q[i - 1]) * (t2 - hydrogram->t[i - 1]);
}

/**
 * \fn int geometry_read(Geometry *geometry, FILE *file)
 * \brief Function to read the data of a channel geometry.
 * \param geometry
 * \brief channel geometry struct.
 * \param file
 * \brief input file.
 * \return 0 on error, 1 on success.
 */
int geometry_read(Geometry *geometry, FILE *file)
{
	unsigned int i;
	char *msg;
	if (fscanf(file, "%u", &geometry->n) != 1 || geometry->n < 2)
	{
		msg = "geometry: bad points number";
		goto bad;
	}
#if DEBUG_CHANNEL
	printf("geometry: n=%u\n", geometry->n);
#endif
	geometry->x = (double*)malloc(geometry->n * sizeof(double));
	geometry->zb = (double*)malloc(geometry->n * sizeof(double));
	geometry->B0 = (double*)malloc(geometry->n * sizeof(double));
	geometry->Z = (double*)malloc(geometry->n * sizeof(double));
	geometry->zmax = (double*)malloc(geometry->n * sizeof(double));
	if (!geometry->x || !geometry->zb || !geometry->B0 || !geometry->Z
		|| ! geometry->zmax)
	{
		msg = "geometry: not enough memory";
		goto bad;
	}
	for (i = 0; i < geometry->n; ++i)
	{
		if (fscanf(file, "%lf%lf%lf%lf%lf",
			geometry->x + i,
			geometry->zb + i,
			geometry->B0 + i,
			geometry->Z + i,
			geometry->zmax + i) != 5)
		{
			msg = "geometry: bad defined";
			goto bad;
		}
#if DEBUG_CHANNEL
		printf("geometry: x=%lg zb=%lg B0=%lg Z=%lg zmax=%lg\n",
			geometry->x[i], geometry->zb[i], geometry->B0[i], geometry->Z[i],
			geometry->zmax[i]);
#endif
		if (geometry->zmax[i] <= geometry->zb[i])
		{
			msg = "geometry: bad zmax";
			goto bad;
		}
		if (i > 0 && geometry->x[i] < geometry->x[i - 1])
		{
			msg = "geometry: bad x";
			goto bad;
		}
	}
	return 1;

bad:
	print_error(msg);
	return 0;
}

/**
 * \fn int channel_friction_read_Manning(Channel *channel, FILE *file)
 * \brief Function to read the friction coefficient of the Manning model.
 * \param channel
 * \brief channel struct.
 * \param file
 * \brief input file.
 * \return 0 on error, 1 on success.
 */
int channel_friction_read_Manning(Channel *channel, FILE *file)
{
	if (fscanf(file, "%lf", channel->friction_coefficient) != 1
		|| channel->friction_coefficient[0] < 0.)
	{
		print_error("channel friction: bad defined");
		return 0;;
	}
#if DEBUG_CHANNEL
	printf("channel friction: coefficient1=%lg\n",
		channel->friction_coefficient[0]);
#endif
	return 1;
}

/**
 * \fn int channel_infiltration_read_KostiakovLewis(Channel *channel, FILE *file)
 * \brief function to read the infiltration coefficients of the Kostiakov-Lewis
 *   model.
 * \param channel
 * \brief channel struct.
 * \param file
 * \brief input file.
 * \return 0 on error, 1 on success.
 */
int channel_infiltration_read_KostiakovLewis(Channel *channel, FILE *file)
{
	if (fscanf(file, "%lf%lf%lf%lf",
		channel->infiltration_coefficient,
		channel->infiltration_coefficient + 1,
		channel->infiltration_coefficient + 2,
		channel->infiltration_coefficient + 3) != 4
		|| channel->infiltration_coefficient[0] < 0.
		|| channel->infiltration_coefficient[1] < 0.
		|| channel->infiltration_coefficient[3] <= 0.)
	{
		print_error("channel infiltration: bad defined");
		return 0;;
	}
#if DEBUG_CHANNEL
	printf("channel infiltration:\n"
		"coefficient1=%lg\n"
		"coefficient2=%lg\n"
		"coefficient3=%lg\n"
		"coefficient4=%lg\n",
		channel->infiltration_coefficient[0],
		channel->infiltration_coefficient[1],
		channel->infiltration_coefficient[2],
		channel->infiltration_coefficient[3]);
#endif
	return 1;
}

/**
 * \fn int channel_diffusion_read_Rutherford(Channel *channel, FILE *file)
 * \brief Function to read the diffusion coefficient of the Rutherford model.
 * \param channel
 * \brief channel struct.
 * \param file
 * \brief input file.
 * \return 0 on error, 1 on success.
 */
int channel_diffusion_read_Rutherford(Channel *channel, FILE *file)
{
	if (fscanf(file, "%lf", channel->diffusion_coefficient) != 1
		|| channel->diffusion_coefficient[0] < 0.)
	{
		print_error("channel diffusion: bad defined");
		return 0;;
	}
	return 1;
}

/**
 * \fn int channel_read(Channel *channel, FILE *file)
 * \brief function to read a channel.
 * \param channel
 * \brief channel struct.
 * \param file
 * \brief input file.
 * \return 0 on error, 1 on success
 */
int channel_read(Channel *channel, FILE *file)
{
	char *msg;
	if (fscanf(file, "%u%u%u%u%u",
		&channel->type_inlet,
		&channel->type_outlet,
		&channel->friction_model,
		&channel->infiltration_model,
		&channel->diffusion_model) != 5)
	{
		msg = "channel: bad defined\n";
		goto bad;
	}
#if DEBUG_CHANNEL
	printf("channel:\n"
		"type_inlet=%u type_outlet=%u\n"
		"friction_model=%u infiltration_model=%u diffusion_model=%u\n",
		channel->type_inlet,
		channel->type_outlet,
		channel->friction_model,
		channel->infiltration_model,
		channel->diffusion_model);
#endif
	if (!geometry_read(channel->geometry, file))
	{
		msg = "channel: geometry";
		goto bad;
	}
	channel->length = channel->geometry->x[channel->geometry->n - 1]
		- channel->geometry->x[0];
	if (channel->length <= 0.)
	{
		msg = "channel: bad length";
		goto bad;
	}
	switch (channel->type_outlet)
	{
	case 1:
	case 2:
		break;
	default:
		msg = "channel: bad outlet";
		goto bad;
	}
	switch (channel->friction_model)
	{
	case 1:
		if (!channel_friction_read_Manning(channel, file)) return 0;
		break;
	default:
		msg = "channel: bad friction model";
		goto bad;
	}
	switch (channel->infiltration_model)
	{
	case 1:
		if (!channel_infiltration_read_KostiakovLewis(channel, file)) return 0;
		break;
	default:
		msg = "channel: bad infiltration model";
		goto bad;
	}
	switch (channel->diffusion_model)
	{
	case 1:
		if (!channel_diffusion_read_Rutherford(channel, file)) return 0;
		break;
	default:
		msg = "channel: bad diffusion model";
		goto bad;
	}
	if (!hydrogram_read(channel->water_inlet, file))
	{
		msg = "channel: inlet";
		goto bad;
	}
	if (!hydrogram_read(channel->solute_inlet, file))
	{
		msg = "channel: outlet";
		goto bad;
	}
	return 1;

bad:
	print_error(msg);
	return 0;
}
