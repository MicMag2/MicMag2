/*
 * Copyright 2012, 2013 by the Micromagnum authors.
 *
 * This file is part of MicroMagnum.
 * 
 * MicroMagnum is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * MicroMagnum is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with MicroMagnum.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"
#include "fs_fdm_slonchewski_cpu.h"
#include <iostream>

#include "mmm/constants.h"

void fs_fdm_slonchewski_cpu(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	const Matrix &Ms,
	const VectorMatrix &J,
	const Matrix &alpha,
	const Matrix &xi,
	const Matrix &alpha_hall,
	const VectorMatrix &M,
	VectorMatrix &dM, const Matrix &mu)
{
	// Calculate: 
	//   c1*(M x (M x p)) + c2*(M x p)
	//
	//   c1(theta): damping factor
	//   c2(theta): precession factor
	//   Ms*cos(theta) = M*p

	Matrix::ro_accessor Ms_acc(Ms), alpha_acc(alpha), xi_acc(xi), alpha_hall_acc(alpha_hall), mu_acc(mu);
	VectorMatrix::const_accessor J_acc(J);
	VectorMatrix::const_accessor M_acc(M);
	VectorMatrix::accessor dM_acc(dM);

	const int N = dim_x * dim_y * dim_z;
	for (int n=0; n<N; ++n) {
		const double alpha_n           = alpha_acc.at(n);
		const double xi_n              = xi_acc.at(n);
		const double Ms_n              = Ms_acc.at(n);
		const double mu_n              = mu_acc.at(n);
		const double alpha_hall_n      = alpha_hall_acc.at(n);

		const Vector3d J_n = J_acc.get(n);
		const Vector3d M_n = M_acc.get(n);
		const Vector3d e_z = Vector3d(0.0, 0.0, 1.0);
		
		if (J_n == Vector3d(0.0, 0.0, 0.0)) continue;
		if (Ms_n == 0.0) continue;

		// Calculate precession and damping terms
		const Vector3d p_n   = cross(J_n, e_z); // spin polarization J x e_z with J = technical current direction
		const Vector3d Mxp   = cross(M_n, p_n); // precession: u=mxp
		const Vector3d MxMxp = cross(M_n, Mxp); // damping:  t=mxu=mx(mxp)

		// add both terms to dm/dt in LLGE
		const double gamma_pr = GYROMAGNETIC_RATIO / (1.0 + alpha_n*alpha_n);

		Vector3d dM_n = dM_acc.get(n);
		//std::cout << "xi_n"<< xi_n << "alpha_hall_n" << alpha_hall_n << std::endl;
		//std::cout << "p_n" << p_n << std::endl;
		//const double Jabs = pow(J_n.x*J_n.x + J_n.y*J_n.y + J_n.z*J_n.z, 0.5);		// absolute value of current in p.
		const double a_j = (H_BAR * alpha_hall_n) / (2.* ELECTRON_CHARGE * Ms_n * delta_z * MU0);

		dM_n.x += (gamma_pr * a_j * ((1. + xi_n*alpha_n) * MxMxp.x/mu_n + (xi_n - alpha_n) * Mxp.x));
		dM_n.y += (gamma_pr * a_j * ((1. + xi_n*alpha_n) * MxMxp.y/mu_n + (xi_n - alpha_n) * Mxp.y));
		dM_n.z += (gamma_pr * a_j * ((1. + xi_n*alpha_n) * MxMxp.z/mu_n + (xi_n - alpha_n) * Mxp.z));
		dM_acc.set(n, dM_n);
		//std::cout << "MxMxp = " << MxMxp << "Mxp = "<< Mxp << std::endl;
		//std::cout << "n = " << n << "dM_n = "<< dM_n << std::endl;
	}
}

/*
void fdm_slonchewski(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	double a_j,
	const VectorMatrix &p, // spin polarization
	const Matrix &Ms,
	const Matrix &alpha,
	VectorMatrix &dM)
{
	// Calculate: 
	//   c1*(M x (M x p)) + c2*(M x p)
	//
	//   c1(theta): damping factor
	//   c2(theta): precession factor
	//   Ms*cos(theta) = M*p

	Matrix::ro_accessor Ms_acc(Ms), alpha_acc(alpha);
	VectorMatrix::const_accessor p_acc(p);
	VectorMatrix::const_accessor M_acc(M);
	VectorMatrix::accessor dM_acc(dM);

	const double Lambda2 = params.Lambda * params.Lambda;

	const int N = model()->mesh->totalNodes();
	for (int n=0; n<N; ++n) {
		const double alpha = alpha_acc.at(n);
		const double Ms    = Ms_acc.at(n);

		const Vector3d p = vector_get(p_acc, n);
		const Vector3d M = vector_get(M_acc, n);
		
		if (p == Vector3d(0.0, 0.0, 0.0)) continue;
		if (Ms == 0.0) continue;

		double a_j;
		switch (mode) {
			case MODE_FIXED_AJ: {
				a_j = params.fixed_a_j;
				break;
			}

			case MODE_VARIABLE_AJ: {
				const double cos_theta = dot(M, p) / Ms;
				const double theta = std::acos(cos_theta);
				const double sin_theta = std::sin(theta);
				const double f1 = H_BAR / (2*MU0*ELECTRON_CHARGE) / params.thickness;
				const double f2 = params.P * Lambda2 * sin_theta / ((Lambda2+1)+(Lambda2-1)*cos_theta);
				a_j = f1*f2 * params.J / Ms;
			}

			default: assert(0);
		}

		// Calculate precession and damping terms
		const Vector3d u = cross(M, p); // precession: u=mxp
		const Vector3d t = cross(M, u); // damping:  t=mxu=mx(mxp)

		// add both terms to dm/dt in LLGE
		const double gamma_pr = GYROMAGNETIC_RATIO / (1.0 + alpha*alpha);

		Vector3d dM_n = vector_get(dM_acc, n);
		dM_n.x += gamma_pr * a_j * (-t.x/Ms + u.x*alpha);
		dM_n.y += gamma_pr * a_j * (-t.y/Ms + u.y*alpha);
		dM_n.z += gamma_pr * a_j * (-t.z/Ms + u.z*alpha);
		dM_acc.set(n, dM_n);
	}
}
*/
