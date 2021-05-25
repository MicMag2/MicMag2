#ifndef COSINECHECK_H
#define COSINECHECK_H

#include "Magneto.h"
#include "config.h"
#include "matrix/matty.h"
#include "mesh/VectorField.h"
#include "mesh/Field.h"
#include "Benchmark.h"

class cosinecheck{
	public:
		int idx_x0;
		int idx_y0;
		int idx_z0;
		Vector3d r_0;
		Vector3d r_1;
		int idx_x1;
		int idx_y1;
		int idx_z1;
		float cos_r;
		float cos_u;
		float cos_f;


	cosinecheck(
	int dim_x, int dim_y, 
	int dim_z, const VectorMatrix M,
	const Matrix Ms, float treshold);
	};
	
#endif
