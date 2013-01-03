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
 * \var U
 * \brief conserved variables vector.
 * \var s
 * \brief surface solute concentration.
 * \var si
 * \brief infiltrated solute concentration.
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
 * \var dF
 * \brief flux difference vector.
 * \var dFl
 * \brief left numerical flux difference vector.
 * \var dFr
 * \brief right numerical flux difference vector.
 * \var nu
 * \brief artificial viscosity coefficient.
 * \var Jp
 * \brief positive implicit operator.
 * \var Jn
 * \brief negative implicit operator.
 * \var Un
 * \brief former time step conserved variables vector.
 * \var dU
 * \brief conserved variables increment vector.
 */
	double friction_coefficient[3], infiltration_coefficient[4],
		diffusion_coefficient[1], x, dx, ix, U[5], s, si, h, Sf, zb, zs, P, B,
		u, c, l1, l2, i, Pi, Z, B0, F, T, Kx, KxA, Kxi, KxiA,
		dF[3], dFl[3], dFr[3], nu, Jp[9], Jn[9], Un[3], dU[3];
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

#endif
