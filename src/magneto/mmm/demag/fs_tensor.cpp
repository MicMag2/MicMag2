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

/*
 * This file was developed by Benjamin Krueger.
 */

#include "config.h"
#include "fs_tensor.h"
#include <cfloat>
#include <cstdlib>
#include "../mmm/constants.h"
#include "Logger.h"

double I_T(int row, int column, double dx, double dy, double dz,  int i, int j,  int k)
{
	double rx = dx*i;
	double ry = dy*j;
	double rz = dz*k;
	double r = sqrt(rx*rx +ry*ry + rz*rz);

	if (i == 0 && j == 0 && k == 0) {return 0 ;}


	if (row == 0 && column == 0)
	{ return (1.0/(4.0*PI))*(3*rx*rx - r*r)/pow(r, 5);}

	if (row == 1 && column == 1)
	{ return (1.0/(4.0*PI))*(3*ry*ry - r*r)/pow(r, 5);}

	if (row == 2 && column == 2)
	{ return (1.0/(4.0*PI))*(3*rz*rz - r*r)/pow(r, 5);}
	
	if (row == 0 && column == 1)
	{ return (1.0/(4.0*PI))*(3*rx*ry)/pow(r, 5);}

	if (row == 1 && column == 2)
	{ return (1.0/(4.0*PI))*(3*rz*ry)/pow(r, 5);}

	if (row == 0 && column == 2)
	{ return (1.0/(4.0*PI))*(3*rx*rz)/pow(r, 5);}

 return 0 ;}

Matrix fs_calculateDemagTensor(long double lx, long double ly, long double lz, int nx, int ny, int nz, int ex, int ey, int ez, int repeat_x, int repeat_y, int repeat_z)
{
	LOG_DEBUG << "calculateDemagTensor: calculating with " << LDBL_MANT_DIG << " bit working precision.";
	if (LDBL_MANT_DIG < 64) LOG_WARN << "calculateDemagTensor: Your working precision is below 64 bit. This might result in too low accuracy.";

	Matrix N_matrix = zeros(Shape(6, ex, ey, ez)); 
	{
		Matrix::rw_accessor N_acc(N_matrix);

		const bool no_infinity_correction = std::getenv("MAGNUM_DEMAG_NO_INFINITY_CORRECTION");

		int percent = 0;

		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++)
			for (int k = 0; k < nz; k++)
			{
				double Nxx = 0, Nyy = 0, Nzz = 0, Nxy = 0, Nyz = 0, Nxz = 0;
				for (int rx = 0; rx < repeat_x; ++rx)
				for (int ry = 0; ry < repeat_y; ++ry)
				for (int rz = 0; rz < repeat_z; ++rz) {
					int mcs = 1;
					if (i+rx*nx == 0) mcs *= 2;
					if (j+ry*ny == 0) mcs *= 2;
					if (k+rz*nz == 0) mcs *= 2;
					Nxx += -I_T(0, 0, lx, ly, lz, i+rx*nx, j+ry*ny, k+rz*nz)/mcs;
					Nyy += -I_T(1, 1, lx, ly, lz, i+rx*nx, j+ry*ny, k+rz*nz)/mcs;
					Nzz += -I_T(2, 2, lx, ly, lz, i+rx*nx, j+ry*ny, k+rz*nz)/mcs;
					Nxy += -I_T(0, 1, lx, ly, lz, i+rx*nx, j+ry*ny, k+rz*nz)/mcs;
		 			Nyz += -I_T(1, 2, lx, ly, lz, i+rx*nx, j+ry*ny, k+rz*nz)/mcs;
					Nxz += -I_T(0, 2, lx, ly, lz, i+rx*nx, j+ry*ny, k+rz*nz)/mcs;
				}

				for (int o = -1; o <= 1; o += 2)
				for (int p = -1; p <= 1; p += 2)
				for (int q = -1; q <= 1; q += 2) {
					// (I,J,K): Coordinates in K matrix for octant (o,p,q).
					const int I = (o*i+ex) % ex;
					const int J = (p*j+ey) % ey;
					const int K = (q*k+ez) % ez;
					N_acc.at(0, I, J, K) +=     Nxx;
					N_acc.at(1, I, J, K) += o*p*Nxy;
					N_acc.at(2, I, J, K) += o*q*Nxz;
					N_acc.at(3, I, J, K) +=     Nyy;
					N_acc.at(4, I, J, K) += p*q*Nyz;
					N_acc.at(5, I, J, K) +=     Nzz;
				}
			}

			while (100.0 * (i+1) / nx >= percent) {
				LOG_INFO << "  " << percent << "%";
				percent += 5;
			}
		}
	}

	return N_matrix;
}
