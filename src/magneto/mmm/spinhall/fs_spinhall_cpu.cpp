#include "config.h"
#include "fs_spinhall_cpu.h"
#include "mmm/constants.h"
#include <iostream>

double fs_spinhall_cpu(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	const Matrix &mu,
	const Matrix &J,
	const double &t,
	const double &Ms,
	const VectorMatrix &M,
	const double &Theta,
	VectorMatrix &H
)
{
		VectorMatrix::const_accessor M_acc(M);
		VectorMatrix::accessor H_acc(H);
		Matrix::ro_accessor spin_acc(mu), J_acc(J);
		Vector3d prod;

        int dim_xy= dim_x*dim_y;
		double energy = 0.0;
	for (int z=0; z<dim_z; ++z) {
		for (int y=0; y<dim_y; ++y) {	
			for (int x=0; x<dim_x; ++x) {
				const int i = z*dim_xy + y*dim_x + x; // linear index of (x,y,z)
				const double spin = spin_acc.at(i);
				double J_i = J_acc.at(i);
				const Vector3d M_i = M_acc.get(i);
				if (spin == 0.0) {
					H_acc.set(i, Vector3d(0.0, 0.0, 0.0));
					continue;
				}			
				prod = Vector3d(-M_i[2], 0.0, M_i[0]);
				Vector3d H_i = Theta*J_i*2*PI*H_BAR/(2*ELECTRON_CHARGE*Ms*t)*prod/sqrt(dot(prod, prod));
				energy += dot(M_i, H_i);
}}}
return energy; }
