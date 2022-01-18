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
#include "interlayerExchange_cuda.h"

#include <cuda.h>

#include "config.h"
#include "mmm/constants.h"
#include <iostream>

static const int BLOCK_3D_SIZE_X = 8;
static const int BLOCK_3D_SIZE_Y = 8;
static const int BLOCK_3D_SIZE_Z = 8;

///////////////////////////////////////////////////////////////////
// KERNEL FOR and 3D meshes                                  //
///////////////////////////////////////////////////////////////////

template <typename real, bool periodic_x, bool periodic_y, bool periodic_z>
__global__ 
void kernel_interlayerExchange_3d(
	const real *Mx, const real *My, const real *Mz, 
	real *Hx, real *Hy, real *Hz, 
	const real *Ms, const real *patx, const real *paty, const real *patz,
	int dim_x, int dim_y, int dim_z, int logical_grid_dim_y, double scale)
{
	// Thread index (inside block)
	const int tx = threadIdx.x;
	const int ty = threadIdx.y;
	const int tz = threadIdx.z;

	// Cell index
	const int sx =  blockIdx.x                       * BLOCK_3D_SIZE_X + tx;
	const int sy = (blockIdx.y % logical_grid_dim_y) * BLOCK_3D_SIZE_Y + ty;
	const int sz = (blockIdx.y / logical_grid_dim_y) * BLOCK_3D_SIZE_Z + tz;

	// Copy into shared mem
	//__shared__ real shared[4][BLOCK_3D_SIZE_Z+2][BLOCK_3D_SIZE_Y+2][BLOCK_3D_SIZE_X+2];

	if (sx < dim_x && sy < dim_y && sz < dim_z) {
		
		// I. Prepare indices
		const int i = sx + dim_x*sy + dim_x*dim_y*sz; // center pos
		
		if (Ms[i] == 0.0) {
			Hx[i] = 0.0;
			Hy[i] = 0.0;
			Hz[i] = 0.0;
		}
		else {

			real sumx = 0.0;
			real sumy = 0.0;
			real sumz = 0.0;


			int layer1 = patx[i];
			int layer2 = paty[i];
			if (layer1 != -1 && layer1 != sz){
				int interact_linear =  sx + dim_x*sy + dim_x*dim_y*layer1; // cell to interact with

				// calculate interlayer exchange
				if (Ms[interact_linear] != 0.0) {
					sumx += patz[i]*(Mx[interact_linear] / Ms[interact_linear]);
					sumy += patz[i]*(My[interact_linear] / Ms[interact_linear]);
					sumz += patz[i]*(Mz[interact_linear] / Ms[interact_linear]);
				}
			}
			if (layer2 != -1 && layer2 != sz){
				int interact_linear =  sx + dim_x*sy + dim_x*dim_y*layer2; // cell to interact with

				// calculate interlayer exchange
				if (Ms[interact_linear] != 0.0) {
					sumx += patz[i]*(Mx[interact_linear] / Ms[interact_linear]);
					sumy += patz[i]*(My[interact_linear] / Ms[interact_linear]);
					sumz += patz[i]*(Mz[interact_linear] / Ms[interact_linear]);
				}
			}

			// Exchange field at (x,y,z)
			Hx[i] = (1.0/MU0) * (1.0/scale) * sumx / Ms[i];
			Hy[i] = (1.0/MU0) * (1.0/scale) * sumy / Ms[i];
			Hz[i] = (1.0/MU0) * (1.0/scale) * sumz / Ms[i];
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// MAIN ROUTINE THAT CALLS THE KERNELS                                      //
//////////////////////////////////////////////////////////////////////////////

template <typename real>
double interlayerExchange_cuda_impl(
	int dim_x, int dim_y, int dim_z,		
	double delta_x, double delta_y, double delta_z,
	bool periodic_x, bool periodic_y, bool periodic_z,
	const Matrix &Ms,
	const VectorMatrix &intExchMatrix,
	const VectorMatrix &M,
	VectorMatrix &H)
{
	{
		typename Matrix_const_cuda_accessor<real>::t Ms_acc(Ms); 
		typename VectorMatrix_const_cuda_accessor<real>::t M_acc(M);
		typename VectorMatrix_const_cuda_accessor<real>::t pattern_acc(intExchMatrix);
		typename VectorMatrix_cuda_accessor<real>::t H_acc(H);

		const real *Mx = M_acc.ptr_x(); real *Hx = H_acc.ptr_x();
		const real *My = M_acc.ptr_y(); real *Hy = H_acc.ptr_y();
		const real *Mz = M_acc.ptr_z(); real *Hz = H_acc.ptr_z();

		const real *patx = pattern_acc.ptr_x();
		const real *paty = pattern_acc.ptr_y();
		const real *patz = pattern_acc.ptr_z();

		// Precalculate weights that are used in the kernels.
		//const real wx = static_cast<real>(2.0/MU0) / (delta_x * delta_x);
		//const real wy = static_cast<real>(2.0/MU0) / (delta_y * delta_y);
		//const real wz = static_cast<real>(2.0/MU0) / (delta_z * delta_z);


		dim3 block_dim(BLOCK_3D_SIZE_X, BLOCK_3D_SIZE_Y, BLOCK_3D_SIZE_Z);
		dim3 grid_dim(
			(dim_x + BLOCK_3D_SIZE_X-1) / BLOCK_3D_SIZE_X, 
			(dim_y + BLOCK_3D_SIZE_Y-1) / BLOCK_3D_SIZE_Y,
			(dim_z + BLOCK_3D_SIZE_Z-1) / BLOCK_3D_SIZE_Z
		);

		// Only 2-dimensional grids are supported, so ...
		const int logical_grid_dim_y = grid_dim.y;
		grid_dim.y *= grid_dim.z;
		grid_dim.z = 1;

		#define INTEXCH_3D(bx,by,bz) if (periodic_x == bx && periodic_y == by && periodic_z == bz) kernel_interlayerExchange_3d<real, bx, by, bz><<<grid_dim, block_dim>>>(Mx, My, Mz, Hx, Hy, Hz, Ms_acc.ptr(), patx, paty, patz, dim_x, dim_y, dim_z, logical_grid_dim_y, delta_z * delta_z);
		INTEXCH_3D(false, false, false)
		INTEXCH_3D(false, false,  true)
		INTEXCH_3D(false,  true, false)
		INTEXCH_3D(false,  true,  true)
		INTEXCH_3D( true, false, false)
		INTEXCH_3D( true, false,  true)
		INTEXCH_3D( true,  true, false)
		INTEXCH_3D( true,  true,  true)
		#undef INTEXCH_3D

		checkCudaLastError("gpu_interlayerExchange(): kernel_interlayerExchange_3d execution failed!");

		CUDA_THREAD_SYNCHRONIZE();
	}

	// and calculate exchange energy
	//std::cout << M.dotSum(H) <<std::endl;
	return -MU0 * M.dotSum(H) * delta_x * delta_y * delta_z;
}

double interlayerExchange_cuda(
	int dim_x, int dim_y, int dim_z,		
	double delta_x, double delta_y, double delta_z,
	bool periodic_x, bool periodic_y, bool periodic_z,
	const Matrix &Ms,
	const VectorMatrix &intExchMatrix,
	const VectorMatrix &M,
	VectorMatrix &H,
	bool cuda64)
{

#ifdef HAVE_CUDA_64
	if (cuda64)
	return interlayerExchange_cuda_impl<double>(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, periodic_x, periodic_y, periodic_z, Ms, intExchMatrix, M, H);
	else
#endif
	return interlayerExchange_cuda_impl<float>(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, periodic_x, periodic_y, periodic_z, Ms, intExchMatrix, M, H);
}
