%{
#include "mmm/amr/AMR.h"
%}

class AMR
{
public:
        AMR(int dim_x, int dim_y, int dim_z, double delta_x, double delta_y, double delta_z, std::vector<int>, std::vector<int>, std::vector<int>);
	virtual ~AMR();

	void calculate(const VectorMatrix &M);
	void get_j(VectorMatrix &j, double U);
	void get_phi(Matrix &phi, double U);
	double get_resistance();

	void setSigma(const Matrix &sigma);
	void setAMR(const Matrix &amr, int amr_dimension);
};
