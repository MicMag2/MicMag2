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
#include "interlayerExchange_cpu.h"
#include "mmm/constants.h"
#include <iostream>

static double interlayerExchange_cpu_nonperiodic(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	const Matrix &Ms,
	const VectorMatrix &intExchMatrix,
	const VectorMatrix &M,
	VectorMatrix &H
);

double interlayerExchange_cpu(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	bool periodic_x, bool periodic_y, bool periodic_z,
	const Matrix &Ms,
	const VectorMatrix &intExchMatrix,
	const VectorMatrix &M,
	VectorMatrix &H)
{
	//const bool periodic = periodic_x || periodic_y || periodic_z;
	return interlayerExchange_cpu_nonperiodic(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Ms, intExchMatrix, M, H);
}

static double interlayerExchange_cpu_nonperiodic(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	const Matrix &Ms,
	const VectorMatrix &intExchMatrix,
	const VectorMatrix &M,
	VectorMatrix &H)
{
	const int dim_xy = dim_x * dim_y;
	//const double wx = 1.0 / (delta_x * delta_x);
	//const double wy = 1.0 / (delta_y * delta_y);
	//const double wz = 1.0 / (delta_z * delta_z);

	VectorMatrix::const_accessor M_acc(M), pattern_acc(intExchMatrix);
	VectorMatrix::accessor H_acc(H);
	Matrix::ro_accessor Ms_acc(Ms);

	double energy = 0.0;
	double scale = delta_z * delta_z;
	for (int z=0; z<dim_z; ++z) {
		for (int y=0; y<dim_y; ++y) {	
			for (int x=0; x<dim_x; ++x) {
				const int i = z*dim_xy + y*dim_x + x; // linear index of (x,y,z)
				const double Ms = Ms_acc.at(i);
				if (Ms == 0.0) {
					H_acc.set(i, Vector3d(0.0, 0.0, 0.0));
					continue;
				}
				const Vector3d M_i = M_acc.get(i); // magnetization at (x,y,z)
				//std::cout<< "test" << std::endl;
				Vector3d sum(0.0, 0.0, 0.0);			
				
				int layer1 = pattern_acc.get(i)[0];
				int layer2 = pattern_acc.get(i)[1];
				if (layer1 != -1 && layer1 != z){
					int interact_linear = layer1*dim_xy + y*dim_x + x;

					// calculate interlayer exchange
					const double Ms_inter = Ms_acc.at(interact_linear);
					if (Ms_inter != 0.0) sum += pattern_acc.get(i)[2]*(M_acc.get(interact_linear) / Ms_inter);
				}
				if (layer2 != -1 && layer2 != z){
					int interact_linear = layer2*dim_xy + y*dim_x + x;

					// calculate interlayer exchange
					const double Ms_inter = Ms_acc.at(interact_linear);
					if (Ms_inter != 0.0) sum += pattern_acc.get(i)[2]*(M_acc.get(interact_linear) / Ms_inter);
				}

				// Exchange field at (x,y,z)
				const Vector3d H_i = (1.0/MU0) * (1.0/scale) * sum / Ms;
				H_acc.set(i, H_i);

				// Exchange energy sum
				energy += dot(M_i , H_i);
			}
		}
	}

	energy *= -MU0 * delta_x * delta_y * delta_z;
	return energy;
}
