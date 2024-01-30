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
#include "anisotropy_new.h"
#include "anisotropy_new_cpu.h"
#ifdef HAVE_CUDA
#include "anisotropy_new_cuda.h"
#include <cuda_runtime.h>
#endif

#include "Magneto.h"
#include "Benchmark.h"

#include <cassert>

double anisotropy_new(
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
	const bool use_cuda = isCudaEnabled();

	double energy_sum = 0.0;

	if (use_cuda) {
#ifdef HAVE_CUDA
// 		CUTIC("uniaxial_anisotropy_new"); --> TODO:What does this do????
		energy_sum = anisotropy_new_cuda(axisu, ku1, ku2, axis1, axis2, k1, k2, Ms, M, H_aniso, isCuda64Enabled());
// 		CUTOC("uniaxial_anisotropy_new");
#else
		assert(0);
#endif
	} else {
// 		TIC("uniaxial_anisotropy_new");
		energy_sum = anisotropy_new_cpu(axisu, ku1, ku2, axis1, axis2, k1, k2, Ms, M, H_aniso);
// 		TOC("uniaxial_anisotropy_new");
	}

	return energy_sum;
}
