#ifndef CPU_FDM_SPINHALL_H
#define CPU_FDM_SPINHALL_H

#include "config.h"
#include "matrix/matty.h"

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
);


#endif
