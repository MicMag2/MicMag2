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
#include "fdm_slonchewski_cuda.h"

#include <cuda.h>

#include "config.h"
#include "mmm/constants.h"

#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16

// Cross product: A = B x C
template <typename real>
static __inline__ __device__ void gpu_cross(
	real *ax, real *ay, real *az,
	const real  bx, const real  by, const real  bz,
	const real  cx, const real  cy, const real  cz)
{
	*ax = by*cz - bz*cy;
	*ay = bz*cx - bx*cz;
	*az = bx*cy - by*cx;
}

///////////////////////////////////////////////////////////////////
// KERNELS FOR 2D and 3D meshes                                  //
///////////////////////////////////////////////////////////////////

template <typename real>
__global__ 
void kernel_fdm_slonchewski_2d(
	const real dim_x, const real dim_y, 
	const real delta_x, const real delta_y, const real delta_z, 
	const real *Mx, const real *My, const real *Mz, 
	real *dMx, real *dMy, real *dMz, 
	const real *Ms, 
	const real *jx, const real *jy, const real *jz, 
	const real *alpha,
	const real *xi,
	const real *alpha_hall)
{
	// Cell index
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	const int y = blockIdx.y * blockDim.y + threadIdx.y;


	if (x < dim_x && y < dim_y) {
		const int k = y*dim_x + x; // k: linear index of (x,y)

		const real Mx_k = Mx[k], My_k = My[k], Mz_k = Mz[k];
		const real jx_k = jx[k], jy_k = jy[k], jz_k = jz[k];
		const real Ms_k    = Ms[k];
		const real alpha_k = alpha[k];	
		const real xi_k = xi[k];	
		const real alpha_hall_k = alpha_hall[k];
		real p_x, p_y, p_z;
		real Mxp_x, Mxp_y, Mxp_z;
		real MxMxp_x, MxMxp_y, MxMxp_z;


		if (Ms_k > 0.0 && (jx_k != 0 || jy_k != 0 || jz_k != 0)) {

			gpu_cross(&p_x, &p_y, &p_z,
				  (real)(0.0), (real)(0.0), (real)(1.0),
				  jx_k, jy_k, jz_k); // spin polarization e_z x J with J = technical current direction
			gpu_cross(&Mxp_x, &Mxp_y, &Mxp_z,
				  Mx_k, My_k, Mz_k,
				  p_x, p_y, p_z); // precession: u=mxp
			gpu_cross(&MxMxp_x, &MxMxp_y, &MxMxp_z,
				  Mx_k, My_k, Mz_k,
				  Mxp_x, Mxp_y, Mxp_z); // damping:  t=mxu=mx(mxp)

			const real gamma_pr = GYROMAGNETIC_RATIO / (1.0 + alpha_k*alpha_k);
			//const double Jabs = - pow(jx_k*jx_k + jy_k*jy_k + jz_k*jz_k, 0.5);  // absolute value of current in p.
			const real a_j = (H_BAR * alpha_hall_k) / (2.* ELECTRON_CHARGE * Ms_k * delta_z * MU0);

			dMx[k] += gamma_pr * a_j * ((1.0 + xi_k*alpha_k) * MxMxp_x/Ms_k + (xi_k - alpha_k) *  Mxp_x);
			dMy[k] += gamma_pr * a_j * ((1.0 + xi_k*alpha_k) * MxMxp_y/Ms_k + (xi_k - alpha_k) *  Mxp_y);
			dMz[k] += gamma_pr * a_j * ((1.0 + xi_k*alpha_k) * MxMxp_z/Ms_k + (xi_k - alpha_k) *  Mxp_z);


		} else {
			dMx[k] += 0.0;
			dMy[k] += 0.0;
			dMz[k] += 0.0;
		}
	}
}

template <typename real>
__global__ 
void kernel_fdm_slonchewski_3d(
	const real dim_x, const real dim_y, const real dim_z, 
	const real delta_x, const real delta_y, const real delta_z, 
	const real *Mx, const real *My, const real *Mz, 
	real *dMx, real *dMy, real *dMz, 
	const real *Ms, 
	const real *jx,	const real *jy,	const real *jz,
	const real *alpha,
	const real *xi,
	const real *alpha_hall)
{

	// Cell index
	const int x = blockIdx.x * blockDim.x + threadIdx.x;
	const int y = blockIdx.y * blockDim.y + threadIdx.y;
	const int dim_xy = dim_x * dim_y;

	for (int z=0; z<dim_z; ++z) {
		if (x < dim_x && y < dim_y) {
			const int k = z*dim_xy + y*dim_x + x; // k: linear index of (x,y)

			const real Mx_k = Mx[k], My_k = My[k], Mz_k = Mz[k];
			const real jx_k = jx[k], jy_k = jy[k], jz_k = jz[k];
			const real Ms_k    = Ms[k];
			const real alpha_k = alpha[k];	
			const real xi_k = xi[k];	
			const real alpha_hall_k = alpha_hall[k];
			real p_x, p_y, p_z;
			real Mxp_x, Mxp_y, Mxp_z;
			real MxMxp_x, MxMxp_y, MxMxp_z;

			if (Ms_k > 0.0 && (jx_k != 0 || jy_k != 0 || jz_k != 0)) {

				gpu_cross<real>(&p_x, &p_y, &p_z,
				  	  (real)(0.0), (real)(0.0), (real)(1.0),
				  	  jx_k, jy_k, jz_k); // spin polarization e_z x J with J = technical current direction
				gpu_cross<real>(&Mxp_x, &Mxp_y, &Mxp_z,
					  Mx_k, My_k, Mz_k,
					  p_x, p_y, p_z); // precession: u=mxp
				gpu_cross<real>(&MxMxp_x, &MxMxp_y, &MxMxp_z,
					  Mx_k, My_k, Mz_k,
					  Mxp_x, Mxp_y, Mxp_z); // damping:  t=mxu=mx(mxp)

				const real gamma_pr = GYROMAGNETIC_RATIO / (1.0 + alpha_k*alpha_k);
				//const double Jabs = - pow(jx_k*jx_k + jy_k*jy_k + jz_k*jz_k, 0.5);  // absolute value of current in p.
				const real a_j = (H_BAR * alpha_hall_k) / (2.* ELECTRON_CHARGE * Ms_k * delta_z * MU0);

				dMx[k] += gamma_pr * a_j * (-(1.0 + xi_k*alpha_k) * MxMxp_x/Ms_k - (xi_k - alpha_k) *  Mxp_x);
				dMy[k] += gamma_pr * a_j * (-(1.0 + xi_k*alpha_k) * MxMxp_y/Ms_k - (xi_k - alpha_k) *  Mxp_y);
				dMz[k] += gamma_pr * a_j * (-(1.0 + xi_k*alpha_k) * MxMxp_z/Ms_k - (xi_k - alpha_k) *  Mxp_z);

			} else {
				dMx[k] += 0.0;
				dMy[k] += 0.0;
				dMz[k] += 0.0;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// MAIN ROUTINE THAT CALLS THE KERNELS                                      //
//////////////////////////////////////////////////////////////////////////////

template <typename real>
void fdm_slonchewski_cuda_impl(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	const Matrix &Ms,
	const VectorMatrix &J,
	const Matrix &alpha,
	const Matrix &xi,
	const Matrix &alpha_hall,
	const VectorMatrix &M,
	VectorMatrix &dM)
{
	{
		typename Matrix_const_cuda_accessor<real>::t Ms_acc(Ms), alpha_hall_acc(alpha_hall), alpha_acc(alpha), xi_acc(xi); 
		typename VectorMatrix_const_cuda_accessor<real>::t M_acc(M), J_acc(J);
		typename VectorMatrix_cuda_accessor<real>::t dM_acc(dM);

		//const real *Mx = M_acc.ptr_x(); real *dMx = dM_acc.ptr_x();
		//const real *My = M_acc.ptr_y(); real *dMy = dM_acc.ptr_y();
		//const real *Mz = M_acc.ptr_z(); real *dMz = dM_acc.ptr_z();
		//const real *px = p_acc.ptr_x(); const real *py = p_acc.ptr_y(); const real *pz = p_acc.ptr_z(); 


		const bool is_2d = (dim_z == 1);
		if (is_2d) { // call 2d kernel
			const dim3 grid_dim((dim_x+BLOCK_SIZE_X-1) / BLOCK_SIZE_X, (dim_y+BLOCK_SIZE_Y-1) / BLOCK_SIZE_Y);
			const dim3 block_dim(BLOCK_SIZE_X, BLOCK_SIZE_Y);

			kernel_fdm_slonchewski_2d<real><<<grid_dim, block_dim>>>(dim_x, dim_y, delta_x, delta_y, delta_z, M_acc.ptr_x(), M_acc.ptr_y(), M_acc.ptr_z(), dM_acc.ptr_x(), dM_acc.ptr_y(), dM_acc.ptr_z(), Ms_acc.ptr(), J_acc.ptr_x(), J_acc.ptr_y(), J_acc.ptr_z(), alpha_acc.ptr(), xi_acc.ptr(), alpha_hall_acc.ptr());

			checkCudaLastError("gpu_SLONCHEWSKI(): kernel_slonchewski_2d execution failed!");

			CUDA_THREAD_SYNCHRONIZE();

		} else { // call 3d kernel
			const dim3 grid_dim((dim_x+BLOCK_SIZE_X-1) / BLOCK_SIZE_X, (dim_y+BLOCK_SIZE_Y-1) / BLOCK_SIZE_Y);
			const dim3 block_dim(BLOCK_SIZE_X, BLOCK_SIZE_Y);
	

			kernel_fdm_slonchewski_3d<real><<<grid_dim, block_dim>>>(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, M_acc.ptr_x(), M_acc.ptr_y(), M_acc.ptr_z(), dM_acc.ptr_x(), dM_acc.ptr_y(), dM_acc.ptr_z(), Ms_acc.ptr(), J_acc.ptr_x(), J_acc.ptr_y(), J_acc.ptr_z(), alpha_acc.ptr(), xi_acc.ptr(), alpha_hall_acc.ptr());

			checkCudaLastError("gpu_SLONCHEWSKI(): kernel_slonchewski_3d execution failed!");

			CUDA_THREAD_SYNCHRONIZE();
		}
	}
}

void fdm_slonchewski_cuda(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	const Matrix &Ms,
	const VectorMatrix &J, // technical current
	const Matrix &alpha,
	const Matrix &xi,
	const Matrix &alpha_hall,
	const VectorMatrix &M,
	VectorMatrix &dM,
	bool cuda64)
{
#ifdef HAVE_CUDA_64
	if (cuda64)
	return fdm_slonchewski_cuda_impl<double>(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Ms, J, alpha, xi, alpha_hall, M, dM);
	else
#endif
	return fdm_slonchewski_cuda_impl<float>(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Ms, J, alpha, xi, alpha_hall, M, dM);
}
