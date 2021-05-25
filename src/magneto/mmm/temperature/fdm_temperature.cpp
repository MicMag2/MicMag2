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
#include <iostream>

#include "config.h"
#include "fdm_temperature.h"
#include "fdm_temperature_cpu.h"
#ifdef HAVE_CUDA
#include "fdm_temperature_cuda.h"
#include <cuda_runtime.h>
#endif

#include "Magneto.h"
#include "Benchmark.h"

void fdm_temperature(
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
	const bool use_cuda = isCudaEnabled();

	if (use_cuda) {
#ifdef HAVE_CUDA
		//std::cout << "test\n";
		CUTIC("fdm_temperature");
		fdm_temperature_cuda(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Ms, alpha, kelv,  dtime, step, seed, Hth, isCuda64Enabled());
		CUTOC("fdm_temperature");
#else
		assert(0);
#endif
	} else {
		TIC("fdm_temperature");
		fdm_temperature_cpu(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Ms, alpha, kelv, dtime, step, seed, Hth);
		TOC("fdm_temperature");
	}
}

