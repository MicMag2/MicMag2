#ifndef AMRESISTANCE_H
#define AMRESISTANCE_H

#include "config.h"

#include "CurrentCalculator.h"
#include "matrix/matty.h"

#include <memory>

// Wrapper class for CurrentCalculator that uses the matrix classes. The AMR class is exposed
// via the Python interface.

class AMR
{
public:
	AMR(
		int dim_x, int dim_y, int dim_z, double delta_x, double delta_y, double delta_z, 
		std::vector<int> contact1, std::vector<int> contact2, std::vector<int> conductive_cells
	);
	virtual ~AMR();

	void calculate(const VectorMatrix &M);
	void get_j(VectorMatrix &j, double U);
	void get_phi(Matrix &phi, double U);
	double get_resistance() { return cc->get_R(); }

	void setSigma(const Matrix &sigma);
	void setAMR(const Matrix &amr, int amr_dimension);

private:
	const int dim_x, dim_y, dim_z;
	const double delta_x, delta_y, delta_z;
	std::unique_ptr<CurrentCalculator> cc;
};

#endif
