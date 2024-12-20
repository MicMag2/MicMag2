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
#include "fs_dmi_cuda.h"

#include <cuda.h>

#include "config.h"
#include "mmm/constants.h"

static const int BLOCK_2D_SIZE_X = 16;
static const int BLOCK_2D_SIZE_Y = 16;

static const int BLOCK_3D_SIZE_X = 8;
static const int BLOCK_3D_SIZE_Y = 8;
static const int BLOCK_3D_SIZE_Z = 8;

///////////////////////////////////////////////////////////////////
// KERNELS FOR 2D and 3D meshes                                  //
///////////////////////////////////////////////////////////////////

template <typename real, bool periodic_x, bool periodic_y>
__global__ 
void fs_kernel_dmi_2d(
	const real *Mx, const real *My, const real *Mz, 
	real *Hx, real *Hy, real *Hz,
	const real delta_x, const real delta_y, const real delta_z, 
	const real *Ms,
	const real *Dx_x, const real *Dx_y, const real *Dx_z, 
	const real *Dy_x, const real *Dy_y, const real *Dy_z, 
	int dim_x, int dim_y)
{
	// Thread index (inside block)
	const int tx = threadIdx.x;
	const int ty = threadIdx.y;

	// Cell index
	const int sx = blockIdx.x * BLOCK_2D_SIZE_X + tx;
	const int sy = blockIdx.y * BLOCK_2D_SIZE_Y + ty;

	// Shared mem
	__shared__ real shared[4][2+BLOCK_2D_SIZE_X][2+BLOCK_2D_SIZE_Y]; // mx,my,mz,Ms

	// I. Prepare indices
	const int i = sx + dim_x*sy; // center pos

	int idx_l = i-     1;
	int idx_r = i+     1;
	int idx_u = i- dim_x;
	int idx_d = i+ dim_x;

	if (periodic_x) {
		if (sx ==       0) idx_l += dim_x;
		if (sx == dim_x-1) idx_r -= dim_x;
	}
	if (periodic_y) {
		if (sy ==       0) idx_u += dim_x*dim_y;
		if (sy == dim_y-1) idx_d -= dim_x*dim_y;
	}

	// II. Prepare shared memory

	// these are indices into shared memory: shared[..][TZ][TY][TX]
	const int TX = tx+1; // own pos
	const int TY = ty+1;
	real Ms_i;

	if (sx < dim_x && sy < dim_y) {
		// Regular, non-ghost cells
		Ms_i = Ms[i];
		if (Ms_i != 0.0) { 
			shared[0][TY][TX] = Mx[i] / Ms_i;
			shared[1][TY][TX] = My[i] / Ms_i;
			shared[2][TY][TX] = Mz[i] / Ms_i;
		}
		shared[3][TY][TX] = Ms_i;

		// Copy ghost cells
		if (threadIdx.x == 0) { // left ghost cells
			real Ms_l = 0;
			if (periodic_x || (sx > 0)) {
				Ms_l = Ms[idx_l];
				if (Ms_l != 0) {
					shared[0][TY][TX-1] = Mx[idx_l] / Ms_l; 
					shared[1][TY][TX-1] = My[idx_l] / Ms_l; 
					shared[2][TY][TX-1] = Mz[idx_l] / Ms_l;
				}
			}
			shared[3][TY][TX-1] = Ms_l;
		} 
		if ((threadIdx.x == BLOCK_2D_SIZE_X-1) || (sx == dim_x-1)) { // right ghost cells
			real Ms_r = 0;
			if (periodic_x || (sx < dim_x-1)) {
				Ms_r = Ms[idx_r];
				if (Ms_r != 0) { 
					shared[0][TY][TX+1] = Mx[idx_r] / Ms_r; 
					shared[1][TY][TX+1] = My[idx_r] / Ms_r; 
					shared[2][TY][TX+1] = Mz[idx_r] / Ms_r; 
				}
			}
			shared[3][TY][TX+1] = Ms_r; 
		}
		if (threadIdx.y == 0) { // top ghost cells
			real Ms_u = 0;
			if (periodic_y || (sy > 0)) {
				Ms_u = Ms[idx_u];
				if (Ms_u != 0) { 
					shared[0][TY-1][TX] = Mx[idx_u] / Ms_u; 
					shared[1][TY-1][TX] = My[idx_u] / Ms_u; 
					shared[2][TY-1][TX] = Mz[idx_u] / Ms_u;
				}
			}
			shared[3][TY-1][TX] = Ms_u;
		} 
		if ((threadIdx.y == BLOCK_2D_SIZE_Y-1) || (sy == dim_y-1)) { // bottom ghost cells
			real Ms_d = 0;
			if (periodic_y || (sy < dim_y-1)) {
				Ms_d = Ms[idx_d];
				if (Ms_d != 0) { 
					shared[0][TY+1][TX] = Mx[idx_d] / Ms_d; 
					shared[1][TY+1][TX] = My[idx_d] / Ms_d; 
					shared[2][TY+1][TX] = Mz[idx_d] / Ms_d; 
				}
			}
			shared[3][TY+1][TX] = Ms_d; 
		}
	}

	__syncthreads();

	if (sx < dim_x && sy < dim_y) {
		const real Dxx = Dx_x[i];
        	const real Dxy = Dx_y[i];
        	const real Dxz = Dx_z[i];

        	const real Dyx = Dy_x[i];
        	const real Dyy = Dy_y[i];
        	const real Dyz = Dy_z[i];

		// III. Compute the finite differences
		if (Ms_i > 0) {
			real sum[3] = {0,0,0};

			if (shared[3][TY][TX-1] != 0) {
				sum[0] -= (Dxy*shared[2][TY][TX-1] - Dxz*shared[1][TY][TX-1]); 
				sum[1] -= (Dxz*shared[0][TY][TX-1] - Dxx*shared[2][TY][TX-1]);  
				sum[2] -= (Dxx*shared[1][TY][TX-1] - Dxy*shared[0][TY][TX-1]);  
			}

			if (shared[3][TY][TX+1] != 0) {
				sum[0] += (Dxy*shared[2][TY][TX+1] - Dxz*shared[1][TY][TX+1]); 
				sum[1] += (Dxz*shared[0][TY][TX+1] - Dxx*shared[2][TY][TX+1]);  
				sum[2] += (Dxx*shared[1][TY][TX+1] - Dxy*shared[0][TY][TX+1]); 
			}

			if (shared[3][TY-1][TX] != 0) {
				sum[0] -= (Dyy*shared[2][TY-1][TX] - Dyz*shared[1][TY-1][TX]); 
				sum[1] -= (Dyz*shared[0][TY-1][TX] - Dyx*shared[2][TY-1][TX]);  
				sum[2] -= (Dyx*shared[1][TY-1][TX] - Dyy*shared[0][TY-1][TX]);
			}

			if (shared[3][TY+1][TX] != 0) {
				sum[0] += (Dyy*shared[2][TY+1][TX] - Dyz*shared[1][TY+1][TX]); 
				sum[1] += (Dyz*shared[0][TY+1][TX] - Dyx*shared[2][TY+1][TX]);  
				sum[2] += (Dyx*shared[1][TY+1][TX] - Dyy*shared[0][TY+1][TX]);
			}

			const real factor = 1.0 / (MU0*Ms_i);
			Hx[i] = sum[0] * factor;
			Hy[i] = sum[1] * factor;
			Hz[i] = sum[2] * factor;
		} else {
			Hx[i] = 0;
			Hy[i] = 0;
			Hz[i] = 0;
		}
	} 
}

template <typename real, bool periodic_x, bool periodic_y, bool periodic_z>
__global__ 
void fs_kernel_dmi_3d(
	const real *Mx, const real *My, const real *Mz, 
	real *Hx, real *Hy, real *Hz, 
	const real delta_x, const real delta_y, const real delta_z, 
	const real *Ms, 
	const real *Dx_x, const real *Dx_y, const real *Dx_z, 
	const real *Dy_x, const real *Dy_y, const real *Dy_z, 
	const real *Dz_x, const real *Dz_y, const real *Dz_z, 
	int dim_x, int dim_y, int dim_z, int logical_grid_dim_y)
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
	__shared__ real shared[4][BLOCK_3D_SIZE_Z+2][BLOCK_3D_SIZE_Y+2][BLOCK_3D_SIZE_X+2];

	// I. Prepare indices
	const int dim_xy = dim_x*dim_y;
	const int i = sx + dim_x*sy + dim_xy*sz; // center pos

	int idx_l = i -      1;
	int idx_r = i +      1;
	int idx_u = i -  dim_x;
	int idx_d = i +  dim_x;
	int idx_f = i - dim_xy;
	int idx_b = i + dim_xy;

	if (periodic_x) {
		if (sx ==       0) idx_l += dim_x;
		if (sx == dim_x-1) idx_r -= dim_x;
	}
	if (periodic_y) {
		if (sy ==       0) idx_u += dim_xy;
		if (sy == dim_y-1) idx_d -= dim_xy;
	}
	if (periodic_z) {
		if (sz ==       0) idx_f += dim_xy*dim_z;
		if (sz == dim_z-1) idx_d -= dim_xy*dim_z;
	}

	// II. Prepare shared memory

	// these are indices into shared memory: shared[..][TZ][TY][TX]
	const int TX = tx+1; // own pos in shared mem
	const int TY = ty+1;
	const int TZ = tz+1;
	real Ms_i;

	if (sx < dim_x && sy < dim_y && sz < dim_z) {
		// Regular, non-ghost cells
		Ms_i = Ms[i];
		if (Ms_i > 0.0) { 
			shared[0][TZ][TY][TX] = Mx[i] / Ms_i; 
			shared[1][TZ][TY][TX] = My[i] / Ms_i; 
			shared[2][TZ][TY][TX] = Mz[i] / Ms_i; 
		}
		shared[3][TZ][TY][TX] = Ms_i;

		// Copy Ghost cells
		if (threadIdx.x == 0) { // left ghost cells
			real Ms_l = 0.0;
			if (periodic_x || (sx > 0)) {
				Ms_l = Ms[idx_l];
				if (Ms_l > 0.0) {
					shared[0][TZ][TY][TX-1] = Mx[idx_l] / Ms_l; 
					shared[1][TZ][TY][TX-1] = My[idx_l] / Ms_l; 
					shared[2][TZ][TY][TX-1] = Mz[idx_l] / Ms_l;
				}
			}
			shared[3][TZ][TY][TX-1] = Ms_l;
		} 
		if ((threadIdx.x == BLOCK_3D_SIZE_X-1) || (sx == dim_x-1)) { // right ghost cells
			real Ms_r = 0.0;
			if (periodic_x || sx < dim_x-1) {
				Ms_r = Ms[idx_r];
				if (Ms_r > 0.0) { 
					shared[0][TZ][TY][TX+1] = Mx[idx_r] / Ms_r; 
					shared[1][TZ][TY][TX+1] = My[idx_r] / Ms_r; 
					shared[2][TZ][TY][TX+1] = Mz[idx_r] / Ms_r; 
				}
			}
			shared[3][TZ][TY][TX+1] = Ms_r; 
		}
		if (threadIdx.y == 0) { // top ghost cells
			real Ms_u = 0.0;
			if (periodic_y || (sy > 0)) {
				Ms_u = Ms[idx_u];
				if (Ms_u > 0.0) { 
					shared[0][TZ][TY-1][TX] = Mx[idx_u] / Ms_u; 
					shared[1][TZ][TY-1][TX] = My[idx_u] / Ms_u; 
					shared[2][TZ][TY-1][TX] = Mz[idx_u] / Ms_u;
				}
			}
			shared[3][TZ][TY-1][TX] = Ms_u;
		} 
		if ((threadIdx.y == BLOCK_3D_SIZE_Y-1) || (sy == dim_y-1)) { // bottom ghost cells
			real Ms_d = 0.0;
			if (periodic_y || (sy < dim_y-1)) {
				Ms_d = Ms[idx_d];
				if (Ms_d > 0.0f) { 
					shared[0][TZ][TY+1][TX] = Mx[idx_d] / Ms_d; 
					shared[1][TZ][TY+1][TX] = My[idx_d] / Ms_d; 
					shared[2][TZ][TY+1][TX] = Mz[idx_d] / Ms_d; 
				}
			}
			shared[3][TZ][TY+1][TX] = Ms_d; 
		}
		if (threadIdx.z == 0) { // front ghost cells
			real Ms_f = 0.0;
			if (periodic_z || (sz > 0)) {
				Ms_f = Ms[idx_f];
				if (Ms_f > 0.0) { 
					shared[0][TZ-1][TY][TX] = Mx[idx_f] / Ms_f; 
					shared[1][TZ-1][TY][TX] = My[idx_f] / Ms_f; 
					shared[2][TZ-1][TY][TX] = Mz[idx_f] / Ms_f;
				}
			}
			shared[3][TZ-1][TY][TX] = Ms_f;
		} 
		if ((threadIdx.z == BLOCK_3D_SIZE_Z-1) || (sz == dim_z-1)) { // back ghost cells
			real Ms_b = 0.0;
			if (periodic_z || (sz < dim_z-1)) {
				Ms_b = Ms[idx_b];
				if (Ms_b > 0.0) { 
					shared[0][TZ+1][TY][TX] = Mx[idx_b] / Ms_b; 
					shared[1][TZ+1][TY][TX] = My[idx_b] / Ms_b; 
					shared[2][TZ+1][TY][TX] = Mz[idx_b] / Ms_b;
				}
			}
			shared[3][TZ+1][TY][TX] = Ms_b;
		}
	}
	__syncthreads();
	if (sx < dim_x && sy < dim_y && sz < dim_z) {
		const real Dxx = Dx_x[i];
  		const real Dxy = Dx_y[i];
	        const real Dxz = Dx_z[i];

        	const real Dyx = Dy_x[i];
	        const real Dyy = Dy_y[i];
	        const real Dyz = Dy_z[i];

        	const real Dzx = Dz_x[i];
        	const real Dzy = Dz_y[i];
        	const real Dzz = Dz_z[i];

		// III. Compute the finite differences
		if (Ms_i > 0.0) {
			real sum[3] = {0,0,0};

			if (shared[3][TZ][TY][TX-1] != 0) {
				sum[0] -= (Dxy*shared[2][TZ][TY][TX-1] - Dxz*shared[1][TZ][TY][TX-1]); 
				sum[1] -= (Dxz*shared[0][TZ][TY][TX-1] - Dxx*shared[2][TZ][TY][TX-1]);  
				sum[2] -= (Dxx*shared[1][TZ][TY][TX-1] - Dxy*shared[0][TZ][TY][TX-1]);  
			}

			if (shared[3][TZ][TY][TX+1] != 0) {
				sum[0] += (Dxy*shared[2][TZ][TY][TX+1] - Dxz*shared[1][TZ][TY][TX+1]); 
				sum[1] += (Dxz*shared[0][TZ][TY][TX+1] - Dxx*shared[2][TZ][TY][TX+1]);  
				sum[2] += (Dxx*shared[1][TZ][TY][TX+1] - Dxy*shared[0][TZ][TY][TX+1]); 
			}

			if (shared[3][TZ][TY-1][TX] != 0) {
				sum[0] -= (Dyy*shared[2][TZ][TY-1][TX] - Dyz*shared[1][TZ][TY-1][TX]); 
				sum[1] -= (Dyz*shared[0][TZ][TY-1][TX] - Dyx*shared[2][TZ][TY-1][TX]);  
				sum[2] -= (Dyx*shared[1][TZ][TY-1][TX] - Dyy*shared[0][TZ][TY-1][TX]);
			}

			if (shared[3][TZ][TY+1][TX] != 0) {
				sum[0] += (Dyy*shared[2][TZ][TY+1][TX] - Dyz*shared[1][TZ][TY+1][TX]); 
				sum[1] += (Dyz*shared[0][TZ][TY+1][TX] - Dyx*shared[2][TZ][TY+1][TX]);  
				sum[2] += (Dyx*shared[1][TZ][TY+1][TX] - Dyy*shared[0][TZ][TY+1][TX]);
			}

			if (shared[3][TZ-1][TY][TX] != 0) {
				sum[0] -= (Dzy*shared[2][TZ-1][TY][TX] - Dzz*shared[1][TZ-1][TY][TX]); 
				sum[1] -= (Dzz*shared[0][TZ-1][TY][TX] - Dzx*shared[2][TZ-1][TY][TX]);  
				sum[2] -= (Dzx*shared[1][TZ-1][TY][TX] - Dzy*shared[0][TZ-1][TY][TX]);
			}

			if (shared[3][TZ+1][TY][TX] != 0) {
				sum[0] += (Dzy*shared[2][TZ+1][TY][TX] - Dzz*shared[1][TZ+1][TY][TX]); 
				sum[1] += (Dzz*shared[0][TZ+1][TY][TX] - Dzx*shared[2][TZ+1][TY][TX]);  
				sum[2] += (Dzx*shared[1][TZ+1][TY][TX] - Dzy*shared[0][TZ+1][TY][TX]);
			}

			const real factor = 1.0 / (MU0*Ms_i);
			Hx[i] = sum[0] * factor;
			Hy[i] = sum[1] * factor;
			Hz[i] = sum[2] * factor;

		} else {
			Hx[i] = 0.0;
			Hy[i] = 0.0;
			Hz[i] = 0.0;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// MAIN ROUTINE THAT CALLS THE KERNELS                                      //
//////////////////////////////////////////////////////////////////////////////

template <typename real>
double fs_dmi_cuda_impl(
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
	{
		typename Matrix_const_cuda_accessor<real>::t Ms_acc(Ms); 
		typename VectorMatrix_const_cuda_accessor<real>::t M_acc(M), Dx_acc(Dx), Dy_acc(Dy), Dz_acc(Dz);
		typename VectorMatrix_cuda_accessor<real>::t H_acc(H);

		const real *Mx = M_acc.ptr_x(); real *Hx = H_acc.ptr_x();
		const real *My = M_acc.ptr_y(); real *Hy = H_acc.ptr_y();
		const real *Mz = M_acc.ptr_z(); real *Hz = H_acc.ptr_z();
		const real *Dx_x = Dx_acc.ptr_x(); const real *Dx_y = Dx_acc.ptr_y(); const real *Dx_z = Dx_acc.ptr_z(); 
		const real *Dy_x = Dy_acc.ptr_x(); const real *Dy_y = Dy_acc.ptr_y(); const real *Dy_z = Dy_acc.ptr_z(); 
		const real *Dz_x = Dz_acc.ptr_x(); const real *Dz_y = Dz_acc.ptr_y(); const real *Dz_z = Dz_acc.ptr_z(); 



		const bool is_2d = (dim_z == 1);
		if (is_2d) { // call 2d kernel
			const dim3 grid_dim(
				(dim_x + BLOCK_2D_SIZE_X-1) / BLOCK_2D_SIZE_X, 
				(dim_y + BLOCK_2D_SIZE_Y-1) / BLOCK_2D_SIZE_Y,
				1
			);
			const dim3 block_dim(BLOCK_2D_SIZE_X, BLOCK_2D_SIZE_Y, 1);

			#define DMI_2D(bx,by) if (periodic_x == bx && periodic_y == by) fs_kernel_dmi_2d<real, bx, by><<<grid_dim, block_dim>>>(Mx, My, Mz, Hx, Hy, Hz, delta_x, delta_y, delta_z, Ms_acc.ptr(), Dx_x, Dx_y, Dx_z, Dy_x, Dy_y, Dy_z, dim_x, dim_y);
			DMI_2D(false, false)
			DMI_2D(false,  true)
			DMI_2D( true, false)
			DMI_2D( true,  true)
			#undef DMI_2D

			checkCudaLastError("gpu_dmi(): fs_kernel_dmi_2d execution failed!");

			CUDA_THREAD_SYNCHRONIZE();

		} else { // call 3d kernel
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

			#define DMI_3D(bx,by,bz) if (periodic_x == bx && periodic_y == by && periodic_z == bz) fs_kernel_dmi_3d<real, bx, by, bz><<<grid_dim, block_dim>>>(Mx, My, Mz, Hx, Hy, Hz, delta_x, delta_y, delta_z, Ms_acc.ptr(), Dx_x, Dx_y, Dx_z, Dy_x, Dy_y, Dy_z, Dz_x, Dz_y, Dz_z, dim_x, dim_y, dim_z, logical_grid_dim_y);
			DMI_3D(false, false, false)
			DMI_3D(false, false,  true)
			DMI_3D(false,  true, false)
			DMI_3D(false,  true,  true)
			DMI_3D( true, false, false)
			DMI_3D( true, false,  true)
			DMI_3D( true,  true, false)
			DMI_3D( true,  true,  true)
			#undef DMI_3D

			checkCudaLastError("gpu_DMI(): fs_kernel_dmi_3d execution failed!");

			CUDA_THREAD_SYNCHRONIZE();
		}
	}

	// and calculate dmi energy
	return -MU0/2.0 * M.dotSum(H) * delta_x * delta_y * delta_z;
}

double fs_dmi_cuda(
	int dim_x, int dim_y, int dim_z,		
	double delta_x, double delta_y, double delta_z,
	bool periodic_x, bool periodic_y, bool periodic_z,
	const Matrix &Ms,
	const VectorMatrix &Dx,
	const VectorMatrix &Dy,
	const VectorMatrix &Dz,
	const VectorMatrix &M,
	VectorMatrix &H,
	bool cuda64)
{
#ifdef HAVE_CUDA_64
	if (cuda64)
	return fs_dmi_cuda_impl<double>(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, periodic_x, periodic_y, periodic_z, Ms, Dx, Dy, Dz, M, H);
	else
#endif
	return fs_dmi_cuda_impl<float>(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, periodic_x, periodic_y, periodic_z, Ms, Dx, Dy, Dz, M, H);
}
