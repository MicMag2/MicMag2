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
#include "fs_fdm_slonchewski.h"
#include "fs_fdm_slonchewski_cpu.h"
#ifdef HAVE_CUDA
#include "fs_fdm_slonchewski_cuda.h"
#include <cuda_runtime.h>
#endif

#include "Magneto.h"
#include "Benchmark.h"

void fs_fdm_slonchewski(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	const Matrix &Ms,
	const VectorMatrix &j,
	const Matrix &alpha,
	const Matrix &xi,
	const Matrix &alpha_hall,
	const VectorMatrix &M,
	VectorMatrix &dM, const Matrix &mu)
{
	const bool use_cuda = isCudaEnabled();

	if (use_cuda) {
#ifdef HAVE_CUDA
		//std::cout << "test\n";
		CUTIC("fs_fdm_slonchewski");
		fs_fdm_slonchewski_cuda(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Ms, j, alpha, xi, alpha_hall, M, dM, isCuda64Enabled());
		CUTOC("fs_fdm_slonchewski");
#else
		assert(0);
#endif
	} else {
		TIC("fdm_slonchewski");
		fs_fdm_slonchewski_cpu(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Ms, j, alpha, xi, alpha_hall, M, dM, mu);
		TOC("fdm_slonchewski");
	}
}

