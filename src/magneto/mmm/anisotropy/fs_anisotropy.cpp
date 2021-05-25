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
#include "fs_anisotropy.h"
#include "fs_anisotropy_cpu.h"
#ifdef HAVE_CUDA
#include "fs_anisotropy_cuda.h"
#include <cuda_runtime.h>
#endif

#include "Magneto.h"
#include "Benchmark.h"

#include <cassert>

double fs_uniaxial_anisotropy(
	const VectorMatrix &axis,
	const       Matrix &k,
	const       Matrix &mu,
	const VectorMatrix &M,
	VectorMatrix &H, const Matrix &Ms)
{
	const bool use_cuda = isCudaEnabled();

	double energy_sum = 0.0;

	if (use_cuda) {
#ifdef HAVE_CUDA
		CUTIC("fs_uniaxial_anisotropy");
		energy_sum = fs_uniaxial_anisotropy_cuda(axis, k, mu, M, H, Ms, use_cuda);
		CUTOC("fs_uniaxial_anisotropy");
#else
		assert(0);
#endif
	} else {
		TIC("fs_uniaxial_anisotropy");
		energy_sum = fs_uniaxial_anisotropy_cpu(axis, k, mu, M, H, Ms);
		TOC("fs_uniaxial_anisotropy");
	}

	return energy_sum;
}

double fs_cubic_anisotropy(
	const VectorMatrix &axis1,
	const VectorMatrix &axis2,
	const       Matrix &k,
	const       Matrix &mu,
	const VectorMatrix &M,
	VectorMatrix &H, const Matrix &Ms)
{
	const bool use_cuda = isCudaEnabled();

	double energy_sum = 0.0;

	if (use_cuda) {
#ifdef HAVE_CUDA
		CUTIC("fs_cubic_anisotropy");
		energy_sum = fs_cubic_anisotropy_cuda(axis1, axis2, k, mu, M, H,Ms, use_cuda);
		CUTOC("fs_cubic_anisotropy");
#else
		assert(0);
#endif
	} else {
		TIC("fs_cubic_anisotropy");
		energy_sum = fs_cubic_anisotropy_cpu(axis1, axis2, k, mu, M, H, Ms);
		TOC("fs_cubic_anisotropy");
	}

	return energy_sum;
} 
