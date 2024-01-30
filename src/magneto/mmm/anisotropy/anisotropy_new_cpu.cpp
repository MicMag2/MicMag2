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
#include <cstddef>
#include "config.h"
#include "anisotropy_new_cpu.h"
#include "mmm/constants.h"


/*
 * This function calculates the anisotropy field H_aniso and the anisotropy energy E_aniso based on a model
 * with a uniaxial and a cubic anisotropy. The uniaxial anisotropy is defined by the anisotropy constant
 * ku and the axis axisu, whereas the cubic anisotropy is defined by the anisotropy constant k1 and the axes
 * axis1 and axis2.
 * 
 * The anisotropy energy E_aniso is returned as a return value and the anisotropy field is stored in
 * the given VectorField H_aniso.
 */
double anisotropy_new_cpu(
	const VectorMatrix &axisu,
	const Matrix &ku1,
	const Matrix &ku2,
	const VectorMatrix &axis1,
	const VectorMatrix &axis2,
	const Matrix &k1,
	const Matrix &k2,
	const Matrix &Ms,
	const VectorMatrix &M,
	VectorMatrix &H_aniso)
{
	//
	// Calculate:
	//   H(x,y,z) = sum_i [ 2k_i(x,y,z)/(mu0*Ms(x,y,z)^2) * (M(x,y,z)*axis_i) * axis_i ]
	//
	VectorMatrix::const_accessor axisu_acc(axisu);
	Matrix::ro_accessor ku1_acc(ku1);
	Matrix::ro_accessor ku2_acc(ku2);
	VectorMatrix::const_accessor axis1_acc(axis1);
	VectorMatrix::const_accessor axis2_acc(axis2);
	Matrix::ro_accessor k1_acc(k1);
	Matrix::ro_accessor k2_acc(k2);
	Matrix::ro_accessor Ms_acc(Ms);
	VectorMatrix::const_accessor M_acc(M);
	VectorMatrix::accessor H_acc(H_aniso);
	
	double energy_sum = 0.0;
	
	// For performance reasons, we distinguish between three cases:
	// 1. only uniaxial anisotropy present
	// 2. only cubic anisotropy present
	// 3. both anisotropies present

	// Case 1
	if((k1.absMax()==0) && (k2.absMax()==0))
	{
		const size_t num_nodes = M.size();
		for (size_t i=0; i<num_nodes; ++i) {
			const double Ms = Ms_acc.at(i);
			const double ku1 = ku1_acc.at(i);
			const double ku2 = ku2_acc.at(i);
			if (Ms == 0.0 || (ku1 == 0.0 && ku2 == 0.0)) {
				H_acc.set(i, Vector3d(0.0, 0.0, 0.0));
			} else {
                const Vector3d m = M_acc.get(i) / Ms;
                
				const Vector3d axisu = axisu_acc.get(i);
				const double d = dot(m, axisu);
				const Vector3d Hu = (2.0 * ku1 / Ms / MU0) * d * axisu       +       (4.0 * ku2 / Ms / MU0) * d * d * d * axisu;
                const double Eu = ku1 * (1.0 - d*d)       +       ku2 * (1.0 - d*d)*(1.0 - d*d);
                
				H_acc.set(i, Hu);
				energy_sum += Eu;
			}
		}
		return energy_sum;
	}


	// Case 2
	if((ku1.absMax()==0) && (ku2.absMax()==0))
	{
		const size_t num_nodes = M.size();
		for (size_t i=0; i<num_nodes; ++i) {
			const double Ms = Ms_acc.at(i);
            const double k1 = k1_acc.at(i);
			const double k2 = k2_acc.at(i);
			if (Ms == 0.0 || (k1 == 0.0 && k2 == 0.0)) {
				H_acc.set(i, Vector3d(0.0, 0.0, 0.0));
			} else {
				const Vector3d m = M_acc.get(i)/Ms;
                
				const Vector3d axis1 = axis1_acc.get(i);
				const Vector3d axis2 = axis2_acc.get(i);
				const Vector3d axis3 = cross(axis1, axis2);

				const double a1 = dot(axis1, m), a1sq = a1*a1;
				const double a2 = dot(axis2, m), a2sq = a2*a2;
				const double a3 = dot(axis3, m), a3sq = a3*a3;

				const Vector3d m1 = a1*axis1;
				const Vector3d m2 = a2*axis2;
				const Vector3d m3 = a3*axis3;
				
				const Vector3d H1 = (-2.0/MU0)*k1/Ms * ((a2sq+a3sq)*m1 + (a1sq+a3sq)*m2 + (a1sq+a2sq)*m3)       +       (-2.0/MU0)*k2/Ms * ((a2sq*a3sq)*m1 + (a1sq*a3sq)*m2 + (a1sq*a2sq)*m3);
				const double E1 = k1 * (a1sq*a2sq+a1sq*a3sq+a2sq*a3sq)       +       k2 * (a1sq*a2sq*a3sq);
				
                H_acc.set(i, H1);
				energy_sum += E1;
			}
		}
		return energy_sum;
	}
	
	// Case 3
	const size_t num_nodes = M.size();
	for (size_t i=0; i<num_nodes; ++i) {
		const double Ms = Ms_acc.at(i);
        const double ku1 = ku1_acc.at(i);
		const double ku2 = ku2_acc.at(i);
        const double k1 = k1_acc.at(i);
		const double k2 = k2_acc.at(i);
		if (Ms == 0.0 || (ku1 == 0.0 && k1 == 0.0 && ku2 == 0.0 && k2 == 0.0)) {
			H_acc.set(i, Vector3d(0.0, 0.0, 0.0));
		} else {
			const Vector3d m = M_acc.get(i) / Ms;

            // uniaxial contribution
            const Vector3d axisu = axisu_acc.get(i);
			const double d = dot(m, axisu);
			const Vector3d Hu = (2.0 * ku1 / Ms / MU0) * d * axisu       +       (4.0 * ku2 / Ms / MU0) * d * d * d * axisu;
            const double Eu = ku1 * (1.0 - d*d)       +       ku2 * (1.0 - d*d)*(1.0 - d*d);
			
            // cubic contribution
			const Vector3d axis1 = axis1_acc.get(i);
			const Vector3d axis2 = axis2_acc.get(i);
			const Vector3d axis3 = cross(axis1, axis2);

			const double a1 = dot(axis1, m), a1sq = a1*a1;
			const double a2 = dot(axis2, m), a2sq = a2*a2;
			const double a3 = dot(axis3, m), a3sq = a3*a3;

			const Vector3d m1 = a1*axis1;
			const Vector3d m2 = a2*axis2;
			const Vector3d m3 = a3*axis3;
			
			const Vector3d H1 = (-2.0/MU0)*k1/Ms * ((a2sq+a3sq)*m1 + (a1sq+a3sq)*m2 + (a1sq+a2sq)*m3)       +       (-2.0/MU0)*k2/Ms * ((a2sq*a3sq)*m1 + (a1sq*a3sq)*m2 + (a1sq*a2sq)*m3);
			const double E1 = k1 * (a1sq*a2sq+a1sq*a3sq+a2sq*a3sq)       +       k2 * (a1sq*a2sq*a3sq);
            
            // both together
			H_acc.set(i, Hu+H1);
			energy_sum += Eu + E1;
		}
	}
	
	return energy_sum;
}

/*
double uniaxial_anisotropy_cpu(
	const VectorMatrix &axis,
	const Matrix &k,
	const Matrix &Ms,
	const VectorMatrix &M,
	VectorMatrix &H)
{
	//
	// Calculate:
	//   H(i) = 2(i)k/(mu0*Ms^2) * (M(i)*axis) * axis
	//
	VectorMatrix::const_accessor M_acc(M);
	VectorMatrix::accessor H_acc(H);
	VectorMatrix::const_accessor axis_acc(axis);
	Matrix::ro_accessor Ms_acc(Ms), k_acc(k);

	double energy_sum = 0.0;

	// Compute field
	const size_t num_nodes = M.size();
	for (size_t i=0; i<num_nodes; ++i) {
		const double     Ms = Ms_acc.at(i);
		const double      k = k_acc.at(i);

		if (Ms == 0.0 || k == 0.0) {
			H_acc.set(i, Vector3d(0.0, 0.0, 0.0));
		} else {
			const Vector3d axis = axis_acc.get(i);
			const Vector3d spin = M_acc.get(i) / Ms;

			const double d = dot(spin, axis);

			const Vector3d H = (2.0 * k * d / Ms / MU0) * axis;
			H_acc.set(i, H);

			energy_sum += k * (1.0 - d*d);
		}
	}

	return energy_sum;
}

double cubic_anisotropy_cpu(
	const VectorMatrix &axis1,
	const VectorMatrix &axis2,
	const Matrix &k,
	const Matrix &Ms,
	const VectorMatrix &M,
	VectorMatrix &H)
{
	VectorMatrix::const_accessor M_acc(M);
	VectorMatrix::accessor H_acc(H);
	VectorMatrix::const_accessor axis1_acc(axis1);
	VectorMatrix::const_accessor axis2_acc(axis2);
	Matrix::ro_accessor Ms_acc(Ms), k_acc(k);

	double energy_sum = 0.0;

	// Compute field
	const size_t num_nodes = M.size();
	for (size_t i=0; i<num_nodes; ++i) {
		const double Ms = Ms_acc.at(i);
		if (Ms == 0.0) {
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
			const Vector3d H = (-2/MU0)*k/Ms * ((a2sq+a3sq)*m1 + (a1sq+a3sq)*m2 + (a1sq+a2sq)*m3);
			H_acc.set(i, H);

			energy_sum += k * (a1sq*a2sq+a1sq*a3sq+a2sq*a3sq);
		}
	}

	return energy_sum;
}
*/
