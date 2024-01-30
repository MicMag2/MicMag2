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
#include "anisotropy_new_cuda.h"
#include "mmm/constants.h"

#include <stdexcept>

static const int GRID_SIZE = 32; // TODO: use device properties
static const int BLOCK_SIZE = 128;



// hack: use sum routines from matrix library..
template <typename real> 
double cuda_compute_sum(const real *src, int N)
{
        extern double cuda_sum(const real *src, int N);
        return cuda_sum(src, N);
}


template <typename real>
__global__ void kernel_anisotropy_new(
        const real *axu_ptr, const real *ayu_ptr, const real *azu_ptr,
        const real *ku1_ptr, const real *ku2_ptr,
        const real *ax1_ptr, const real *ay1_ptr, const real *az1_ptr,
        const real *ax2_ptr, const real *ay2_ptr, const real *az2_ptr,
        const real *k1_ptr,  const real *k2_ptr,
        const real *Ms_ptr,
        const real *Mx_ptr, const real *My_ptr, const real *Mz_ptr,
        real *Hx_ptr, real *Hy_ptr, real *Hz_ptr,
        real *E_ptr,
        real ku1absmax, real ku2absmax,
        real k1absmax, real k2absmax,
        int N)
{
        const int     tid = blockDim.x * blockIdx.x + threadIdx.x;
        const int threadN = blockDim.x * gridDim.x;
        
        // For performance reasons, we distinguish between three cases:
        // 1. only uniaxial anisotropy present
        // 2. only cubic anisotropy present
        // 3. both anisotropies present

        // Case 1
        if(k1absmax==0 && k2absmax==0)
        {
                for (int i=tid; i<N; i+=threadN) {
                        // Parameters at cell i
                        const real Ms = Ms_ptr[i]; 
                        const real ku1 = ku1_ptr[i];
						const real ku2 = ku2_ptr[i];
                        if (Ms == 0.0 || (ku1 == 0.0 && ku2 == 0.0)) {
                              Hx_ptr[i] = 0.0;
                              Hy_ptr[i] = 0.0;
                              Hz_ptr[i] = 0.0;
                              E_ptr[i] = 0.0;
                        } else {
                                // compute d = cos(angle between M and unit axis)
                                const real Mx = Mx_ptr[i], My = My_ptr[i], Mz = Mz_ptr[i];
                                const real axu = axu_ptr[i], ayu = ayu_ptr[i], azu = azu_ptr[i];
                                const real  d = (Mx*axu + My*ayu + Mz*azu) / Ms;
                                
                                const real fu = 2.0 * ku1 * d / real(MU0) / Ms + 4.0 * ku2 * d * d * d / real(MU0) / Ms;
                                const real Hxu = fu * axu;
                                const real Hyu = fu * ayu;
                                const real Hzu = fu * azu;
                                const real Eu = ku1 * (1.0 - d*d) + ku2 * (1.0 - d*d)*(1.0 - d*d);
                                
                                Hx_ptr[i] = Hxu;
                                Hy_ptr[i] = Hyu;
                                Hz_ptr[i] = Hzu;
                                E_ptr[i] = Eu;
                        }
                }
                return;
        }
        
        // Case 2
        if(ku1absmax==0 && ku2absmax==0 )
        {
                for (int i=tid; i<N; i+=threadN) {
                        // Parameters at cell i
                        const real Ms = Ms_ptr[i];
                        const real k1 = k1_ptr[i];
						const real k2 = k2_ptr[i];
                        if (Ms == 0.0 || (k1 == 0.0 && k2 == 0.0)) {
                              Hx_ptr[i] = 0.0;
                              Hy_ptr[i] = 0.0;
                              Hz_ptr[i] = 0.0;
                              E_ptr[i] = 0.0;
                        } else {
                              // Parameters for cell i
                              const real Mx = Mx_ptr[i], My = My_ptr[i], Mz = Mz_ptr[i];
                              const real ax1 = ax1_ptr[i], ay1 = ay1_ptr[i], az1 = az1_ptr[i];
                              const real ax2 = ax2_ptr[i], ay2 = ay2_ptr[i], az2 = az2_ptr[i];
                              // Third axis: axis3 = cross(axis1, axis2)
                              const real ax3 = ay1*az2 - az1*ay2;
                              const real ay3 = az1*ax2 - ax1*az2;
                              const real az3 = ax1*ay2 - ay1*ax2;
                              // a1, a2, a3: coordinates of unit magnetization in coordinate system with base (axis1, axis2, axis3).
                              const real a1 = (ax1*Mx + ay1*My + az1*Mz) / Ms;
                              const real a2 = (ax2*Mx + ay2*My + az2*Mz) / Ms;
                              const real a3 = (ax3*Mx + ay3*My + az3*Mz) / Ms;

                              const real a1sq = a1*a1;
                              const real a2sq = a2*a2;
                              const real a3sq = a3*a3;

                              const real f1 = -2.0 * k1 / real(MU0) / Ms;
							  const real f2 = -2.0 * k2 / real(MU0) / Ms;
                              const real Hx1 = f1 * ((a2sq+a3sq)*a1*ax1 + (a1sq+a3sq)*a2*ax2 + (a1sq+a2sq)*a3*ax3);
                              const real Hy1 = f1 * ((a2sq+a3sq)*a1*ay1 + (a1sq+a3sq)*a2*ay2 + (a1sq+a2sq)*a3*ay3);
                              const real Hz1 = f1 * ((a2sq+a3sq)*a1*az1 + (a1sq+a3sq)*a2*az2 + (a1sq+a2sq)*a3*az3);
                              const real Hx2 = f2 * ((a2sq*a3sq)*a1*ax1 + (a1sq*a3sq)*a2*ax2 + (a1sq*a2sq)*a3*ax3);
                              const real Hy2 = f2 * ((a2sq*a3sq)*a1*ay1 + (a1sq*a3sq)*a2*ay2 + (a1sq*a2sq)*a3*ay3);
                              const real Hz2 = f2 * ((a2sq*a3sq)*a1*az1 + (a1sq*a3sq)*a2*az2 + (a1sq*a2sq)*a3*az3);
                              const real E1 = k1 * (a1sq*a2sq+a1sq*a3sq+a2sq*a3sq);
							  const real E2 = k2 * (a1sq*a2sq*a3sq);
                              
                              Hx_ptr[i] = Hx1 + Hx2;
                              Hy_ptr[i] = Hy1 + Hy2;
                              Hz_ptr[i] = Hz1 + Hz2;
                              E_ptr[i] = E1 + E2;
                        }
                }
                return;
        }
        
        // Case 3
        for (int i=tid; i<N; i+=threadN) {
                // Parameters at cell i
                const real Ms = Ms_ptr[i];
                const real ku1 = ku1_ptr[i];
                const real ku2 = ku2_ptr[i];
				const real k1 = k1_ptr[i];
                const real k2 = k2_ptr[i];
                if (Ms == 0.0 || (ku1 == 0.0 && k1 == 0.0 && ku2 == 0.0 && k2 == 0.0)) {
                      Hx_ptr[i] = 0.0;
                      Hy_ptr[i] = 0.0;
                      Hz_ptr[i] = 0.0;
                      E_ptr[i] = 0.0;
                } else {
                      // Parameters for cell i
                      const real Mx = Mx_ptr[i], My = My_ptr[i], Mz = Mz_ptr[i];
                      
                      // uniaxial contribution
                      const real axu = axu_ptr[i], ayu = ayu_ptr[i], azu = azu_ptr[i];
                      const real  d = (Mx*axu + My*ayu + Mz*azu) / Ms;
                      
					  const real fu = 2.0 * ku1 * d / real(MU0) / Ms + 4.0 * ku2 * d * d * d / real(MU0) / Ms;
                      const real Hxu = fu * axu;
                      const real Hyu = fu * ayu;
                      const real Hzu = fu * azu;
                      const real Eu = ku1 * (1.0 - d*d) + ku2 * (1.0 - d*d)*(1.0 - d*d);
                      
                      // cubic contribution
                      const real ax1 = ax1_ptr[i], ay1 = ay1_ptr[i], az1 = az1_ptr[i];
                      const real ax2 = ax2_ptr[i], ay2 = ay2_ptr[i], az2 = az2_ptr[i];
                      // Third axis: axis3 = cross(axis1, axis2)
                      const real ax3 = ay1*az2 - az1*ay2;
                      const real ay3 = az1*ax2 - ax1*az2;
                      const real az3 = ax1*ay2 - ay1*ax2;
                      // a1, a2, a3: coordinates of unit magnetization in coordinate system with base (axis1, axis2, axis3).
                      const real a1 = (ax1*Mx + ay1*My + az1*Mz) / Ms;
                      const real a2 = (ax2*Mx + ay2*My + az2*Mz) / Ms;
                      const real a3 = (ax3*Mx + ay3*My + az3*Mz) / Ms;

                      const real a1sq = a1*a1;
                      const real a2sq = a2*a2;
                      const real a3sq = a3*a3;
                  
					  const real f1 = -2.0 * k1 / real(MU0) / Ms;
					  const real f2 = -2.0 * k2 / real(MU0) / Ms;
                      const real Hx1 = f1 * ((a2sq+a3sq)*a1*ax1 + (a1sq+a3sq)*a2*ax2 + (a1sq+a2sq)*a3*ax3);
                      const real Hy1 = f1 * ((a2sq+a3sq)*a1*ay1 + (a1sq+a3sq)*a2*ay2 + (a1sq+a2sq)*a3*ay3);
                      const real Hz1 = f1 * ((a2sq+a3sq)*a1*az1 + (a1sq+a3sq)*a2*az2 + (a1sq+a2sq)*a3*az3);
                      const real Hx2 = f2 * ((a2sq*a3sq)*a1*ax1 + (a1sq*a3sq)*a2*ax2 + (a1sq*a2sq)*a3*ax3);
                      const real Hy2 = f2 * ((a2sq*a3sq)*a1*ay1 + (a1sq*a3sq)*a2*ay2 + (a1sq*a2sq)*a3*ay3);
                      const real Hz2 = f2 * ((a2sq*a3sq)*a1*az1 + (a1sq*a3sq)*a2*az2 + (a1sq*a2sq)*a3*az3);
                      const real E1 = k1 * (a1sq*a2sq+a1sq*a3sq+a2sq*a3sq);
					  const real E2 = k2 * (a1sq*a2sq*a3sq);
                              
                      // both together

                      Hx_ptr[i] = Hxu + Hx1 + Hx2;
                      Hy_ptr[i] = Hyu + Hy1 + Hy2;
                      Hz_ptr[i] = Hzu + Hz1 + Hz2;
                      E_ptr[i] = Eu + E1 + E2;
                }
        }
        return;
}

//////////////////////////////////////////////////////////////////////////////
// MAIN ROUTINE THAT CALLS THE KERNELS                                      //
//////////////////////////////////////////////////////////////////////////////


template <typename real>
double anisotropy_new_cuda_impl(
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
        typename VectorMatrix_cuda_accessor<real>::t H_acc(H_aniso);
        typename VectorMatrix_const_cuda_accessor<real>::t axisu_acc(axisu), axis1_acc(axis1), axis2_acc(axis2), M_acc(M);
        typename Matrix_const_cuda_accessor<real>::t ku1_acc(ku1), ku2_acc(ku2), k1_acc(k1), k2_acc(k2), Ms_acc(Ms);

        const real *axu_ptr = axisu_acc.ptr_x();
        const real *ayu_ptr = axisu_acc.ptr_y();
        const real *azu_ptr = axisu_acc.ptr_z();
        const real *ax1_ptr = axis1_acc.ptr_x();
        const real *ay1_ptr = axis1_acc.ptr_y();
        const real *az1_ptr = axis1_acc.ptr_z();
        const real *ax2_ptr = axis2_acc.ptr_x();
        const real *ay2_ptr = axis2_acc.ptr_y();
        const real *az2_ptr = axis2_acc.ptr_z();
        const real *ku1_ptr  = ku1_acc.ptr();
		const real *ku2_ptr  = ku2_acc.ptr();
        const real *k1_ptr  = k1_acc.ptr();
		const real *k2_ptr  = k2_acc.ptr();
        const real *Ms_ptr = Ms_acc.ptr();
        const real *Mx_ptr = M_acc.ptr_x();
        const real *My_ptr = M_acc.ptr_y();
        const real *Mz_ptr = M_acc.ptr_z();
              real *Hx_ptr = H_acc.ptr_x();
              real *Hy_ptr = H_acc.ptr_y();
              real *Hz_ptr = H_acc.ptr_z();

        real *E_ptr = 0;
        const int N = M.size();
        const cudaError_t err = cudaMalloc((void**)&E_ptr, N * sizeof(real));
        if (err != cudaSuccess) throw std::runtime_error("cubic_anisotropy_new_cuda_impl: could not allocate memory on GPU");

        real ku1absmax = ku1.absMax();
		real ku2absmax = ku2.absMax();
        real k1absmax = k1.absMax();
        real k2absmax = k2.absMax();
        
        kernel_anisotropy_new<<<GRID_SIZE, BLOCK_SIZE>>>(
                axu_ptr, ayu_ptr, azu_ptr,
                ku1_ptr, ku2_ptr, 
                ax1_ptr, ay1_ptr, az1_ptr,
                ax2_ptr, ay2_ptr, az2_ptr, 
                k1_ptr,  k2_ptr,
                Ms_ptr,
                Mx_ptr, My_ptr, Mz_ptr,
                Hx_ptr, Hy_ptr, Hz_ptr,
                E_ptr,
                ku1absmax, ku2absmax,
                k1absmax, k2absmax,
                N
        );
        checkCudaLastError("kernel_cubic_anisotropy_new() execution failed");
        CUDA_THREAD_SYNCHRONIZE();

        const double energy_sum = cuda_compute_sum<real>(E_ptr, N);
        cudaFree(E_ptr);
        return energy_sum;
}

double anisotropy_new_cuda(
	const VectorMatrix &axisu,
    const Matrix &ku1,
	const Matrix &ku2,
    const VectorMatrix &axis1,
    const VectorMatrix &axis2,
    const Matrix &k1,
	const Matrix &k2,
    const Matrix &Ms,
    const VectorMatrix &M,
    VectorMatrix &H_aniso,
	bool cuda64)
{
#ifdef HAVE_CUDA_64
	if (cuda64) 
	return anisotropy_new_cuda_impl<double>(axisu, ku1, ku2, axis1, axis2, k1, k2, Ms, M, H_aniso); 
	else
#endif
	return anisotropy_new_cuda_impl<float>(axisu, ku1, ku2, axis1, axis2, k1, k2, Ms, M, H_aniso);
}