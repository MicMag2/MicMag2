#include "config.h"
#include "matrix/matty.h"
#include "mesh/VectorField.h"
#include "mesh/Field.h"
#include "Magneto.h"
#include "Benchmark.h"
#include "cosinecheck.h"
#include <iostream>


	cosinecheck::cosinecheck(
	int dim_x, int dim_y, 
	int dim_z, const VectorMatrix M,
	const Matrix Ms, float treshold)

{
		int idx_x0 = -1;
		int idx_y0 = -1;
		int idx_z0 = -1;
		int idx_x1 = -2;
		int idx_y1 = -2;
		int idx_z1 = -2;
		const float costresh = cos(treshold);
		VectorMatrix::const_accessor M_acc(M);
		Matrix::ro_accessor Ms_acc(Ms);

        int dim_xy= (dim_x)*(dim_y);
	for (int z=0; z<dim_z; z++) {
		for (int y=0; y<dim_y; y++) {	
			for (int x=0; x<dim_x; x++) {
				const int i = z*(dim_xy) + y*(dim_x) + x; // linear index of (x,y,z)
				double Ms = Ms_acc.at(i);
				Vector3d M_i = M_acc.get(i);


				int idx_r = i + 1;
				int idx_u = i + dim_x;
				int idx_f = i + dim_xy;
				float cos_r = 2;
				float cos_u = 2;
				float cos_f = 2;

				
				// right (X)								
				if (x < dim_x-1){
					double Ms_r = Ms_acc.at(idx_r);	
					Vector3d M_r = M_acc.get(idx_r);
					cos_r = dot(M_i, M_r)/(Ms_r*Ms);
				}
				// up (Y)
				if (y < dim_y-1) {
					double Ms_u = Ms_acc.at(idx_u);
					Vector3d M_u = M_acc.get(idx_u);
					cos_u = dot(M_i, M_u)/(Ms_u*Ms);
				};
				// forward (Z)
				if (z < dim_z-1) {
					double Ms_f = Ms_acc.at(idx_f);
					Vector3d M_f = M_acc.get(idx_f);
					cos_f = dot(M_i, M_f)/(Ms_f*Ms);
				};

				/*std::cout << x << "\t" << y << "\t" << cos_r << "\t" << cos_u << std::endl;*/

				if (cos_r < costresh)
					{
					 if (idx_x1 == -2) {
					 idx_x0 = x;
					 idx_y0 = y;
					 idx_z0 = z;
					 idx_y1 = idx_y0;
					 idx_z1 = idx_z0;
					 idx_x1 = idx_x0 + 1;} 
					 else {	if(x == idx_x1){idx_x1 += 1;};
					 		if(idx_y1 < y){idx_y1 = y;};
					 		if(idx_z1 < z){idx_z1 = z;};};};		
				
				if (cos_u < costresh)
					{
					 if (idx_y1 == -2) {
					 idx_x0 = x;
					 idx_y0 = y;
					 idx_z0 = z;
					 idx_x1 = idx_x0;
					 idx_z1 = idx_z0;
					 idx_y1 = idx_y0 + 1;} 
					 else {	if(y == idx_y1){idx_y1 += 1;};
					 		if(idx_x1 < x){idx_x1 = x;};
					 		if(idx_z1 < z){idx_z1 = z;};};};		
					/*std::cout << "coordinates are: " << x << "," << y << "," << z << "," << "cos_u is : " << cos_u << "\t" << costresh <<std::endl;*/
					
				if (cos_f < costresh)
					{
					 if (idx_z1 == -2) {
					 idx_x0 = x;
					 idx_y0 = y;
					 idx_z0 = z;
					 idx_y1 = idx_y0;
					 idx_x1 = idx_x0;
					 idx_z1 = idx_z0 + 1;} 
					 else {	if(z == idx_z1){idx_z1 += 1;};
					 		if(idx_y1 < y){idx_y1 = y;};
					 		if(idx_x1 < x){idx_x1 = x;};};};		

					 
					 };};};


	r_0 = Vector3d(idx_x0, idx_y0, idx_z0);
	/*std::cout << r_0 << std::endl;*/
	r_1 = Vector3d(idx_x1, idx_y1, idx_z1); 
	

					
					 }
