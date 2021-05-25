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
#include "dmi.h"
#include "dmi_cpu.h"
#ifdef HAVE_CUDA
#include "dmi_cuda.h"
#include <cuda_runtime.h>
#endif
#include "Magneto.h"
#include "Benchmark.h"

double dmi(
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
	double res = 0;
	const bool use_cuda = isCudaEnabled();

	if (use_cuda) {
#ifdef HAVE_CUDA
		//std::cout << "cuda\n";
		CUTIC("dmi");
		res = dmi_cuda(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, periodic_x, periodic_y, periodic_z, Ms, Dx, Dy, Dz, M, H, isCuda64Enabled());
		CUTOC("dmi");
#else
		assert(0);
#endif

	} else {
		TIC("dmi");
		res = dmi_cpu(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, periodic_x, periodic_y, periodic_z, Ms, Dx, Dy, Dz, M, H);
		TOC("dmi");
	}

	return res;
}

double dmi(
	const Field &Ms,
	const VectorMatrix &Dx,
	const VectorMatrix &Dy,
	const VectorMatrix &Dz,
	const VectorField &M,
	VectorField &H)
{
	const RectangularMesh &mesh = M.getMesh();

	int nx, ny, nz; mesh.getNumNodes(nx, ny, nz);
	double dx, dy, dz; mesh.getDelta(dx, dy, dz);
	std::string pbc; int pbc_reps; mesh.getPeriodicBC(pbc, pbc_reps);

	const bool px = pbc.find("x") != std::string::npos;
	const bool py = pbc.find("y") != std::string::npos;
	const bool pz = pbc.find("z") != std::string::npos;

	return dmi(nx, ny, nz, dx, dy, dz, px, py, pz, Ms, Dx, Dy, Dz, M, H);
}

