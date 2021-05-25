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
#include "fs_fdm_temperature_cpu.h"
#include <iostream>
#include "mmm/constants.h"

#include <random>

void fs_fdm_temperature_cpu(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	const Matrix &Ms,
	const Matrix &alpha,
	const Matrix &kelv,
	const double dtime,
	const double step,
	const double seed,
	VectorMatrix &Hth)
{

	Matrix::ro_accessor Ms_acc(Ms), alpha_acc(alpha), kelv_acc(kelv);
	VectorMatrix::accessor Hth_acc(Hth);

	const int N = dim_x * dim_y * dim_z;
	const double volume = delta_x * delta_y * delta_z;

	std::default_random_engine generator;
	generator.seed( seed + step);
	std::normal_distribution<double> distribution(0,1);

	for (int n=0; n<N; ++n) {
		const double alpha_n           = alpha_acc.at(n);
		const double kelv_n            = kelv_acc.at(n);
		const double Ms_n              = Ms_acc.at(n);

		Vector3d H = Hth_acc.get(n);
	
		// calc random field
		if (kelv_n !=0 && Ms_n !=0){
			// Calculate sigma
			const double sigma = sqrt(2.*BOLTZMANN_CONSTANT*kelv_n*alpha_n/(GYROMAGNETIC_RATIO*Ms_n*dtime*MU0));
			
			// Create random vector
			H.x = sigma*distribution(generator);
			H.y = sigma*distribution(generator);
			H.z = sigma*distribution(generator);
		}
		else {
		        H.x = 0;
                        H.y = 0;
                        H.z = 0;
}
		
		Hth_acc.set(n, H);
	}
}

