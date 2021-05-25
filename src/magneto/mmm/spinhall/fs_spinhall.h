#ifndef SH_H
#define SH_H

#include "config.h"
#include "matrix/matty.h"
#include "mesh/VectorField.h"
#include "mesh/Field.h"

double fs_spinhall(
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

double fs_spinhall(
	const Field &mu,
	const Field &J,
	const double &t,
	const double &Ms,
	const VectorField &M,
	const double &Theta,
	VectorField &H
);

#endif
