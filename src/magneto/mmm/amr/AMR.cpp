#include "config.h"
#include "AMR.h"

#include <stdexcept>
#include <cassert>
#include <iostream>
using namespace std;

AMR::AMR(
	int dim_x, int dim_y, int dim_z, double delta_x, double delta_y, double delta_z,
	std::vector<int> conductive_cells, std::vector<int> contact1_cells, std::vector<int> contact2_cells
)
	: dim_x(dim_x), dim_y(dim_y), dim_z(dim_z), delta_x(delta_x), delta_y(delta_y), delta_z(delta_z)
{
	// Setup cell type array.
	std::vector<CurrentCalculator::CellType> cell_types(dim_x*dim_y*dim_z, CurrentCalculator::NONCONDUCTIVE);
	for (unsigned i=0; i<conductive_cells.size(); ++i) cell_types[conductive_cells[i]] = CurrentCalculator::CONDUCTIVE;
	for (unsigned i=0; i<  contact1_cells.size(); ++i) cell_types[  contact1_cells[i]] = CurrentCalculator::CONTACT1;
	for (unsigned i=0; i<  contact2_cells.size(); ++i) cell_types[  contact2_cells[i]] = CurrentCalculator::CONTACT2;

	// Create actual calculator object
	cc.reset(new CurrentCalculator(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, cell_types, CurrentCalculator::CONTACT2CONTACT_ORDERING));
}

AMR::~AMR()
{
}

void AMR::calculate(const VectorMatrix &M)
{
	if (int(M.dimX()) != dim_x || int(M.dimY()) != dim_y || int(M.dimZ()) != dim_z) throw std::runtime_error("AMR::calculate: Invalid matrix size for 'M' parameter.");

	VectorMatrix::const_accessor M_acc(M);
	cc->set_m((double*)M_acc.ptr_x(), (double*)M_acc.ptr_y(), (double*)M_acc.ptr_z());
	cc->solve(1e-8);
}

void AMR::get_j(VectorMatrix &j, double U)
{
	if (int(j.dimX()) != dim_x || int(j.dimY()) != dim_y || int(j.dimZ()) != dim_z) throw std::runtime_error("AMR::get_j: Invalid matrix size for 'j' parameter.");

	VectorMatrix::accessor j_acc(j);
	cc->get_j((double*)j_acc.ptr_x(), (double*)j_acc.ptr_y(), (double*)j_acc.ptr_z(), U);
}

void AMR::get_phi(Matrix &phi, double U)
{
	if (int(phi.dimX()) != dim_x || int(phi.dimY()) != dim_y || int(phi.dimZ()) != dim_z) throw std::runtime_error("AMR::get_phi: Invalid matrix size for 'phi' parameter.");

	Matrix::wo_accessor phi_acc(phi);
	cc->get_Phi(phi_acc.ptr(), U);
}

void AMR::setSigma(const Matrix &sigma)
{
	if (int(sigma.dimX()) != dim_x || int(sigma.dimY()) != dim_y || int(sigma.dimZ()) != dim_z) throw std::runtime_error("AMR::setSigma: Invalid matrix size for 'sigma' parameter.");

	Matrix::ro_accessor sigma_acc(sigma);
	cc->set_sigma(const_cast<double*>(sigma_acc.ptr()));
}

void AMR::setAMR(const Matrix &amr, int amr_dimension)
{
	if (int(amr.dimX()) != dim_x || int(amr.dimY()) != dim_y || int(amr.dimZ()) != dim_z) throw std::runtime_error("AMR::setAMR: Invalid matrix size for 'amr' parameter.");

	Matrix::ro_accessor amr_acc(amr);
	cc->set_amr(const_cast<double*>(amr_acc.ptr()), amr_dimension);
}
