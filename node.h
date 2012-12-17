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
 * \file node.h
 * \brief Header file to define a mesh node.
 * \author Javier Burguete Tolosa.
 * \copyright Copyright 2011-2012, Javier Burguete Tolosa.
 */

// in order to prevent multiple definitions
#ifndef NODE__H
#define NODE__H 1

/**
 * \struct _Node
 * \brief Struct to define a mesh node.
 */
struct _Node
{
/**
 * \var friction_coefficient
 * \brief array of friction coefficients.
 * \var infiltration_coefficient
 * \brief array of infiltration coefficients.
 * \var diffusion_coefficient
 * \brief array of diffusion coefficients.
 * \var x
 * \brief position.
 * \var dx
 * \brief cell size.
 * \var ix
 * \brief cell distance.
 * \var A
 * \brief wetted cross sectional area.
 * \var Ai
 * \brief infiltrated cross sectional area.
 * \var Q
 * \brief discharge.
 * \var s
 * \brief solute concentration.
 * \var si
 * \brief infiltrated solute concentration.
 * \var As
 * \brief A * s.
 * \var Asi
 * \brief A * si.
 * \var h
 * \brief depth.
 * \var Sf
 * \brief friction slope.
 * \var zb
 * \brief bottom level.
 * \var zs
 * \brief surface level.
 * \var P
 * \brief wetted perimeter.
 * \var B
 * \brief surface width.
 * \var u
 * \brief velocity.
 * \var c
 * \brief critical velocity.
 * \var l1
 * \brief first eigenvalue.
 * \var l2
 * \brief second eigenvalue.
 * \var i
 * \brief infiltration velocity.
 * \var Pi
 * \brief P * i.
 * \var Z
 * \brief lateral wall slope.
 * \var B0
 * \brief bottom width.
 * \var F
 * \brief A * u * u.
 * \var T
 * \brief Q * s.
 * \var Kx
 * \brief diffusion coefficient.
 * \var KxA
 * \brief Kx * A.
 * \var Kxi
 * \brief soil diffusion coefficient.
 * \var KxiA
 * \brief Kxi * A.
 * \var dQ
 * \brief mass flux difference.
 * \var dF
 * \brief momentum flux difference.
 * \var dT
 * \brief solute mass flux difference.
 * \var dQl
 * \brief left numerical mass flux difference.
 * \var dFl
 * \brief left numerical momentum flux difference.
 * \var dTl
 * \brief left numerical solute mass flux difference.
 * \var dQr
 * \brief right numerical mass flux difference.
 * \var dFr
 * \brief right numerical momentum flux difference.
 * \var dTr
 * \brief right numerical solute mass flux difference.
 * \var nu
 * \brief artificial viscosity coefficient.
 */
	double friction_coefficient[3], infiltration_coefficient[4],
		diffusion_coefficient[1], x, dx, ix, A, Ai, Q, s, si, As, Asi, h, Sf,
		zb, zs, P, B, u, c, l1, l2, i, Pi, Z, B0, F, T, Kx, KxA, Kxi, KxiA,
		dQ, dF, dT, dQl, dFl, dTl, dQr, dFr, dTr, nu;
};

/**
 * \typedef Node
 */
typedef struct _Node Node;

// global variables

extern double critical_depth_tolerance;

// member functions

void node_depth(Node *node);
void node_width(Node *node);
void node_perimeter(Node *node);
void node_critical_velocity(Node *node);
void node_subcritical_discharge(Node *node);
double node_critical_depth(Node *node, double Q);
void node_friction_Manning(Node *node);
void node_infiltration_KostiakovLewis(Node *node);
void node_diffusion_Rutherford(Node *node);
void node_inlet
	(Node *node, Hydrogram *water, Hydrogram *solute, double t, double t2);
void node_outlet_closed(Node *node);
void node_outlet_open(Node *node);

#endif
