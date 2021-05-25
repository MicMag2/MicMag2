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
#include "fdm_temperature_cuda.h"
#include "curand.h"
#include "curand_kernel.h"

#include <cuda.h>

#include "config.h"
#include "mmm/constants.h"

#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16

///////////////////////////////////////////////////////////////////
// KERNELS FOR 2D and 3D meshes                                  //
///////////////////////////////////////////////////////////////////


template <typename real>
__global__ 
void kernel_fdm_temperature_2d(
	const real dim_x, const real dim_y,
	const real delta_x, const real delta_y, const real delta_z, 
        real *Hthx, real *Hthy, real *Hthz,
        const real *Ms,
        const real *alpha,
        const real *kelv,
        const real dtime,
	const real step,
	const real seed)
{
	// Cell index
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	const int y = blockIdx.y * blockDim.y + threadIdx.y;

	const real volume = delta_x*delta_y*delta_z;

	curandState_t state;

	if (x < dim_x && y < dim_y) {
		const int k = y*dim_x + x; // k: linear index of (x,y)

		/* Each thread gets same seed, a different sequence number, no offset */
		curand_init(seed+step, k, 0, &state);

		const real Ms_k    = Ms[k];
		const real alpha_k = alpha[k];
		const real kelv_k  = kelv[k];

		if (Ms_k > 0.0 && (kelv_k > 0)) {
			real sigma = sqrt(2.*BOLTZMANN_CONSTANT*kelv_k*alpha_k/(GYROMAGNETIC_RATIO*volume*Ms_k*dtime*MU0));
			Hthx[k] = sigma * curand_normal(&state);
			Hthy[k] = sigma * curand_normal(&state);
			Hthz[k] = sigma * curand_normal(&state);
		} else {
			Hthx[k] = 0;
			Hthy[k] = 0;
			Hthz[k] = 0;
		}
		
	}
}

template <typename real>
__global__ 
void kernel_fdm_temperature_3d(
	const real dim_x, const real dim_y, const real dim_z, 
	const real delta_x, const real delta_y, const real delta_z, 
	real *Hthx, real *Hthy, real *Hthz, 
	const real *Ms, 
	const real *alpha,
	const real *kelv,
	const real dtime,
	const real step,
	const real seed)
{

	// Cell index
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	const int y = blockIdx.y * blockDim.y + threadIdx.y;
	const int dim_xy = dim_x * dim_y;

	const real volume = delta_x*delta_y*delta_z;

	curandState_t state;

	for (int z=0; z<dim_z; ++z) {
		if (x < dim_x && y < dim_y) {
			const int k = z*dim_xy + y*dim_x + x; // k: linear index of (x,y)

			/* Each thread gets same seed, a different sequence number, no offset */
			curand_init(seed+step, k, 0, &state);

			const real Ms_k    = Ms[k];
			const real alpha_k = alpha[k];	
			const real kelv_k  = kelv[k];	
			
			if (Ms_k > 0.0 && (kelv_k > 0)) {

				real sigma = sqrt(2.*BOLTZMANN_CONSTANT*kelv_k*alpha_k/(GYROMAGNETIC_RATIO*volume*Ms_k*dtime*MU0));

				Hthx[k] = sigma * curand_normal(&state);//*curand_normal(curandState_t *state); 
				Hthy[k] = sigma * curand_normal(&state); 
				Hthz[k] = sigma * curand_normal(&state);

			} else {
				Hthx[k] = 0;
				Hthy[k] = 0;
				Hthz[k] = 0;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// MAIN ROUTINE THAT CALLS THE KERNELS                                      //
//////////////////////////////////////////////////////////////////////////////

template <typename real>
void fdm_temperature_cuda_impl(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	const Matrix &Ms,
	const Matrix &alpha,
	const Matrix &kelv,
	const double &dtime,
	const double &step,
	const double &seed,
	VectorMatrix &Hth)
{
	{
		typename Matrix_const_cuda_accessor<real>::t Ms_acc(Ms), kelv_acc(kelv), alpha_acc(alpha); 
		typename VectorMatrix_cuda_accessor<real>::t Hth_acc(Hth);

		const bool is_2d = (dim_z == 1);
		if (is_2d) { // call 2d kernel
			const dim3 grid_dim((dim_x+BLOCK_SIZE_X-1) / BLOCK_SIZE_X, (dim_y+BLOCK_SIZE_Y-1) / BLOCK_SIZE_Y);
			const dim3 block_dim(BLOCK_SIZE_X, BLOCK_SIZE_Y);

			kernel_fdm_temperature_2d<real><<<grid_dim, block_dim>>>(dim_x, dim_y, delta_x, delta_y, delta_z, Hth_acc.ptr_x(), Hth_acc.ptr_y(), Hth_acc.ptr_z(), Ms_acc.ptr(), alpha_acc.ptr(), kelv_acc.ptr(), dtime, step, seed);

			checkCudaLastError("gpu_SLONCHEWSKI(): kernel_slonchewski_2d execution failed!");

			CUDA_THREAD_SYNCHRONIZE();

		} else { // call 3d kernel
			const dim3 grid_dim((dim_x+BLOCK_SIZE_X-1) / BLOCK_SIZE_X, (dim_y+BLOCK_SIZE_Y-1) / BLOCK_SIZE_Y);
			const dim3 block_dim(BLOCK_SIZE_X, BLOCK_SIZE_Y);
	

			kernel_fdm_temperature_3d<real><<<grid_dim, block_dim>>>(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Hth_acc.ptr_x(), Hth_acc.ptr_y(), Hth_acc.ptr_z(), Ms_acc.ptr(), alpha_acc.ptr(), kelv_acc.ptr(), dtime, step, seed);

			checkCudaLastError("gpu_SLONCHEWSKI(): kernel_slonchewski_3d execution failed!");

			CUDA_THREAD_SYNCHRONIZE();
		}
	}
}

void fdm_temperature_cuda(
        int dim_x, int dim_y, int dim_z,
        double delta_x, double delta_y, double delta_z,
        const Matrix &Ms,
        const Matrix &alpha,
        const Matrix &kelv,
        const double dtime,
        const double step,
	const double seed,
        VectorMatrix &Hth,
	bool cuda64)
{

#ifdef HAVE_CUDA_64
	if (cuda64)
	return fdm_temperature_cuda_impl<double>(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Ms, alpha, kelv, dtime, step, seed, Hth);
	else
#endif
	return fdm_temperature_cuda_impl<float>(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Ms, alpha, kelv, dtime, step, seed, Hth);
}
