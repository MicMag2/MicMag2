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
#include "dmi_cpu.h"
#include "mmm/constants.h"
#include <iostream>

double fs_dmi_cpu(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	bool periodic_x, bool periodic_y, bool periodic_z,
	const Matrix &Ms,
	const VectorMatrix &Dx,
	const VectorMatrix &Dy,
	const VectorMatrix &Dz,
	const VectorMatrix &M,
	VectorMatrix &H)
{
	const int dim_xy = dim_x * dim_y;

	VectorMatrix::const_accessor M_acc(M);
	VectorMatrix::accessor H_acc(H);
	Matrix::ro_accessor Ms_acc(Ms);
	VectorMatrix::const_accessor Dx_acc(Dx), Dy_acc(Dy), Dz_acc(Dz);

	double energy = 0.0;
	for (int z=0; z<dim_z; ++z) {
		for (int y=0; y<dim_y; ++y) {	
			for (int x=0; x<dim_x; ++x) {
				const int i = z*dim_xy + y*dim_x + x; // linear index of (x,y,z)
				const double Ms = Ms_acc.at(i);
				if (Ms == 0.0) {
					H_acc.set(i, Vector3d(0.0, 0.0, 0.0));
					continue;
				}
				
				const Vector3d Dx_delta = Dx_acc.get(i)*delta_x*delta_x;
				const Vector3d Dy_delta = Dy_acc.get(i)*delta_y*delta_y;
				const Vector3d Dz_delta = Dz_acc.get(i)*delta_z*delta_z;

				int idx_l = i -      1;
				int idx_r = i +      1;
				int idy_d = i -  dim_x;
				int idy_u = i +  dim_x;
				int idz_b = i - dim_xy;
				int idz_f = i + dim_xy;

				// wrap-around for periodic boundary conditions
				if (periodic_x) {
					if (x ==       0) idx_l += dim_x;
					if (x == dim_x-1) idx_r -= dim_x;
				}
				if (periodic_y) {
					if (y ==       0) idy_d += dim_xy;
					if (y == dim_y-1) idy_u -= dim_xy;
				}
				if (periodic_z) {
					if (z ==       0) idz_b += dim_xy*dim_z;
					if (z == dim_z-1) idz_f -= dim_xy*dim_z;
				}
				Vector3d sum(0.0, 0.0, 0.0);

				// left / right (X)
				if (x >       0 || periodic_x) {
					const double Ms_l = Ms_acc.at(idx_l);
					if (Ms_l != 0.0) sum -= cross(Dx_delta, M_acc.get(idx_l)/Ms_l);
					//if (M_acc.get(idx_l)[0]!=0.0) std::cout << "left_ghost" << M_acc.get(idx_l)/Ms_l << "; "<< "left_DMI" << cross(Dx_delta, M_acc.get(idx_l)) << std::endl;

				}
				if (x < dim_x-1 || periodic_x) {
					const double Ms_r = Ms_acc.at(idx_r);	
					if (Ms_r != 0.0) sum += cross(Dx_delta, M_acc.get(idx_r)/Ms_r);
					//if (M_acc.get(idx_r)[0]!=0.0) std::cout << "right_ghost" << M_acc.get(idx_r)/Ms_r <<  "; "<< "right_DMI" << cross(Dx_delta, M_acc.get(idx_r)) << std::endl;

				}
				// up / down (Y)
				if (y >       0 || periodic_y) {
					const double Ms_d = Ms_acc.at(idy_d);
					if (Ms_d != 0.0) sum -= cross(Dy_delta, M_acc.get(idy_d)/Ms_d);
				}
				if (y < dim_y-1 || periodic_y) {
					const double Ms_u = Ms_acc.at(idy_u);
					if (Ms_u != 0.0) sum += cross(Dy_delta, M_acc.get(idy_u)/Ms_u);
				}
				// forward / backward (Z)
				if (z >       0 || periodic_z) {
					const double Ms_b = Ms_acc.at(idz_b);
					if (Ms_b != 0.0) sum -= cross(Dz_delta, M_acc.get(idz_b)/Ms_b);
				}
				if (z < dim_z-1 || periodic_z) {
					const double Ms_f = Ms_acc.at(idz_f);
					if (Ms_f != 0.0) sum += cross(Dz_delta, M_acc.get(idz_f)/Ms_f);
				}

				// DMI field at (x,y,z)
				const Vector3d H_i = (1/MU0)*sum/Ms_acc.at(i);
				H_acc.set(i, H_i);
//std::cout << H_i << ',' << M_acc.get(i) << std::endl;
				// DMI energy sum
				energy += dot(M_acc.get(i), H_i);
			}
		}
	}

	energy *= -MU0/2.0*delta_x*delta_y*delta_z;
	return energy;
}


