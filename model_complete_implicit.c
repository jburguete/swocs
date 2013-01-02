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

#include <stdio.h>
#include <math.h>
#include "config.h"
#include "channel.h"
#include "node.h"
#include "mesh.h"
#include "model.h"
#include "model_complete_implicit.h"

/**
 * \fn void model_surface_flow_complete_multiply\
 *   (double *m, double *v, double *r)
 * \brief Function to multiply a matrix by a vector of the complete model.
 * \var m
 * \brief multiplying matrix.
 * \var v
 * \brief vector to multiply.
 * \var r
 * \brief resulting vector.
 */
void model_surface_flow complete_multiply(double *m, double *v, double *r)
{
	r[0] = m[0] * v[0] + m[1] * v[1];
	r[1] = m[3] * v[0] + m[4] * v[1];
	r[2] = m[6] * v[0] + m[7] * v[1] + m[8] * v[2];
}

/**
 * \fn void model_surface_flow_complete_invert(double *m, double *i)
 * \brief Function to invert an implicit operator of the complete model.
 * \var m
 * \brief implicit operator to invert.
 * \var i
 * \brief invert operator.
 */
void model_surface_flow complete_invert(double *m, double *i)
{
	double d;
	d = m[8] * (m[0] * m[4] - m[1] * m[3]);
	i[0] = m[4] * m[8] / d;
	i[3] = - m[3] * m[8] / d;
	i[6] = (m[3] * m[7] - m[4] * m[6]) / d;
	i[1] = - m[1] * m[8] / d;
	i[4] = m[0] * m[8] / d;
	i[7] = (m[1] * m[6] - m[0] * m[7]) / d;
	i[2] = 0.;
	i[5] = 0.;
	i[8] = (m[0] * m[4] - m[1] * m[3]) / d;
}

/**
 * \fn void model_surface_flow_complete_implicit(Model *model)
 * \brief Function to make the surface flow with the upwind implicit numerical
 * 	scheme.
 * \param model
 * \brief model struct.
 */
void model_surface_flow_complete_implicit(Model *model)
{
	int i, j, n1;
	double c, u, s, l1, l2, l3, c2, sA1, sA2, k1, k2, odt, A[9], B[9], C[9],
		D[3];
	Mesh *mesh = model->mesh;
	Node *node = mesh->node;
	for (i = 0; i < mesh->n; ++i)
		for (j = 0; j < 3; ++j) node[i].Un[j] = node[i].U[j];
	for (i = 0; i < mesh->n; ++i)
	{
		if (node[i].h <= model->minimum_depth)
		{
			for (j = 0; j < 9; ++j) node[i].Jp[j] = node[i].Jm[j] = 0.;
			continue;
		}
		c2 = 2. * node[i].c;
		l1 = fmax(0., node[i].l1);
		l2 = fmax(0., node[i].l2);
		l3 = fmax(0., node[i].u);
		node[i].Jp[0] = (node[i].l1 * l2 - node[i].l2 * l1) / c2;
		node[i].Jp[1] = (l1 -l2) / c2;
		node[i].Jp[2] = 0.;
		node[i].Jp[3] = - node[i].l1 * node[i].l2 * node[i].Jp[1];
		node[i].Jp[4] = (node[i].l1 * l1 - node[i].l2 * l2) / c2;
		node[i].Jp[5] = 0.;
		node[i].Jp[6] = (node[i].Jp[0] - l3) * node[i].s;
		node[i].Jp[7] = node[i].Jp[1] * node[i].s;
		node[i].Jp[8] = l3;
		l1 = fmin(0., node[i].l1);
		l2 = fmin(0., node[i].l2);
		l3 = fmin(0., node[i].u);
		node[i].Jm[0] = (node[i].l1 * l2 - node[i].l2 * l1) / c2;
		node[i].Jm[1] = (l1 -l2) / c2;
		node[i].Jm[2] = 0.;
		node[i].Jm[3] = - node[i].l1 * node[i].l2 * node[i].Jm[1];
		node[i].Jm[4] = (node[i].l1 * l1 - node[i].l2 * l2) / c2;
		node[i].Jm[5] = 0.;
		node[i].Jm[6] = (node[i].Jm[0] - l3) * node[i].s;
		node[i].Jm[7] = node[i].Jm[1] * node[i].s;
		node[i].Jm[8] = l3;
	}
	n1 = mesh->n - 1;
	for (i = 0; i < n1; ++i)
	{
		model->node_flows(node + i);
		node[i].dFr[0] = node[i].dFr[1] = node[i].dFr[2] = node[i].dFl[0] =
			node[i].dFl[1] = node[i].dFl[2] = 0;
		if (node[i].h <= model->minimum_depth &&
			node[i + 1].h <= model->minimum_depth)
				continue;

		// Roe's averages

		c = sqrt(G * (node[i + 1].U[0] + node[i].U[0]) /
			(node[i + 1].B + node[i].B));
		sA1 = sqrt(node[i].U[0]);
		sA2 = sqrt(node[i + 1].U[0]);
		k2 = sA1 + sA2;
		k1 = sA1 / k2;
		k2 = sA2 / k2;
		u = k1 * node[i].u + k2 * node[i + 1].u;
		l1 = u + c;
		l2 = u - c;

		// wave decomposition

		if (u >= c)
		{
			for (j = 0; j < 3; ++j) node[i].dFl[j] = node[i].dF[j];
		}
		else
		{
			s = k1 * node[i].s + k2 * node[i + 1].s;
			node[i].dFr[0] = 0.5 * (l1 * node[i].dF[0] - node[i].dF[1]) / c;
			node[i].dFr[1] = l2 * node[i].dFr[0];
			node[i].dFr[2] = s * node[i].dFr[0];
			for (j = 0; j < 3; ++j)
				node[i].dFl[j] = node[i].dF[j] - node[i].dFr[j];
		}

		// entropy correction

		if (node[i].l2 < 0. && node[i + 1].l2 > 0.)
			k1 = 0.25 * (node[i + 1].l2 - node[i].l2 - 2 * fabs(l2));
		else if (node[i].l1 < 0. && node[i + 1].l1 > 0.)
			k1 = 0.25 * (node[i + 1].l1 - node[i].l1 - 2 * fabs(l1));
		else continue;

		for (j = 0; j < 3; ++j)
		{
			k2 = k1 * (node[i + 1].U[j] - node[i].U[j]);
			node[i].dFl[j] += k2;
			node[i].dFr[j] -= k2;
		}
	}

	// variables updating

	odt = model->theta * model->dt;
	for (j = 0; j < 9; ++j) B[j] = odt * node[0].JP[j];
	for (j = 0; j < 3; ++j)
	{
		node[0].dU[j] = 0.;
		node[0].U[j] = node[0].Un[j];
	}
	for (i = 0; ++i <= n1;)
	{
		model_surface_flow_complete_implicit_multiply(B, node[i - 1].dU, D);
		for (j = 0; j < 9; ++j) A[j] = B[j] = odt * node[i].JP[j];
		A[0] += node[i].dx;
		A[4] += node[i].dx;
		A[8] += node[i].dx;
		model_surface_flow_complete_implicit_invert(A, C);
		for (j = 0; j < 3; ++j)
			D[j] -= model->dt * node[i - 1].dFl[j] / node[i].dx;
		model_surface_flow_complete_implicit_multiply(C, D, node[i].dU);
		for (j = 0; j < 3; ++j) node[i].U[j] = node[i].Un[j] + node[i].dU[j];
	}
	i = n1;
	for (j = 0; j < 9; ++j) B[j] = odt * node[i].JM[j];
	for (j = 0; j < 3; ++j) node[i].dU[j] = 0.;
	while (--i >= 0)
	{
		model_surface_flow_complete_implicit_multiply(B, node[i + 1].dU, D);
		for (j = 0; j < 9; ++j) A[j] = B[j] = odt * node[i].JP[j];
		A[0] += node[i].dx;
		A[4] += node[i].dx;
		A[8] += node[i].dx;
		model_surface_flow_complete_implicit_invert(A, C);
		for (j = 0; j < 3; ++j) D[j] -= model->dt * node[i].dFr[j] / node[i].dx;
		model_surface_flow_complete_implicit_multiply(C, D, node[i].dU);
		for (j = 0; j < 3; ++j) node[i].U[j] += node[i].dU[j];
	}
}
