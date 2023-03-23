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
#include "fs_anisotropy_cpu.h"
#include "mmm/constants.h"
#include <iostream>
#include <cstddef>

// see comments at the end of the file
//#define USE_CROSS_PRODUCT_LIKE_OOMMF 1

double fs_uniaxial_anisotropy_cpu(
	const VectorMatrix &axis,
	const Matrix &k,
	const Matrix &mu,
	const VectorMatrix &M,
	VectorMatrix &H, const Matrix &Ms)
{

	VectorMatrix::const_accessor M_acc(M);
	VectorMatrix::accessor H_acc(H);
	VectorMatrix::const_accessor axis_acc(axis);
	Matrix::ro_accessor mu_acc(mu), k_acc(k), Ms_acc(Ms);

		double energy = 0.0;
	const size_t num_nodes = M.size();
	for (size_t i=0; i<num_nodes; ++i)  // linear index of (x,y,z)
	 {
				const double spin = mu_acc.at(i); // remember to replace with: spin_acc.at(i);
				const double k = k_acc.at(i);
                //const double Ms = Ms_acc.at(i);

				if (spin == 0.0 || k == 0) {
					H_acc.set(i, Vector3d(0.0, 0.0, 0.0));
					continue;
				}
				Vector3d M_i = M_acc.get(i)/(spin);
				const Vector3d axis = axis_acc.get(i);
        const double d = dot(axis, M_i);
				const Vector3d H_i = 2.0*k * d * axis/(MU0);//*Ms/spin
				H_acc.set(i, H_i);
				//std::cout <<M_i << H_i << std::endl;

				//Uniaxial Anisotropy energy sum
				energy += -k*d*d;//dot(M_i, H_i);
				}
			return energy;			
		}
double fs_cubic_anisotropy_cpu(
	const VectorMatrix &axis1,
	const VectorMatrix &axis2,
	const Matrix &k,
	const Matrix &mu,
	const VectorMatrix &M,
	VectorMatrix &H, const Matrix &Ms)
{
	VectorMatrix::const_accessor M_acc(M);
	VectorMatrix::accessor H_acc(H);
	VectorMatrix::const_accessor axis1_acc(axis1);
	VectorMatrix::const_accessor axis2_acc(axis2);
	Matrix::ro_accessor mu_acc(mu), k_acc(k);

	double energy_sum = 0.0;

	// Compute field
	const size_t num_nodes = M.size();
	for (size_t i=0; i<num_nodes; ++i) {
		const double spin = mu_acc.at(i);
		if (spin == 0.0) {
			H_acc.set(i, Vector3d(0.0, 0.0, 0.0));
		} else {
			const double k = k_acc.at(i);
			const Vector3d m = normalize(M_acc.get(i), 1.0);
			const Vector3d axis1 = axis1_acc.get(i);
			const Vector3d axis2 = axis2_acc.get(i);
			const Vector3d axis3 = cross(axis1, axis2);

			const double a1 = dot(axis1, m), a1sq = a1*a1;
			const double a2 = dot(axis2, m), a2sq = a2*a2;
			const double a3 = dot(axis3, m), a3sq = a3*a3;

			const Vector3d m1 = a1*axis1;
			const Vector3d m2 = a2*axis2;
			const Vector3d m3 = a3*axis3;
			const Vector3d H = (-2/MU0)*k/spin * ((a2sq+a3sq)*m1 + (a1sq+a3sq)*m2 + (a1sq+a2sq)*m3);

			H_acc.set(i, H);

			energy_sum += k * (a1sq*a2sq+a1sq*a3sq+a2sq*a3sq);
		}
	}

	return energy_sum;
}
  