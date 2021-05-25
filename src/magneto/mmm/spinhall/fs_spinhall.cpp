#include "config.h"
#include "fs_spinhall.h"
#include "fs_spinhall_cpu.h"

#include "Magneto.h"
#include "Benchmark.h"

double fs_spinhall(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	const Matrix &mu,
	const Matrix &J,
	const double &t,
	const double &Ms,
	const VectorMatrix &M,
	const double &Theta,
	VectorMatrix &H)
	{ 
	double res = 0;
	TIC("spinhall");
	res = fs_spinhall_cpu(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, mu, J, t, Ms, M, Theta, H);
	TOC("spinhall");
		
		return res;
}

double fs_spinhall(
	const Field &mu,
	const Field &J,
	const double &t,
	const double &Ms,
	const VectorField &M,
	const double &Theta,
	VectorField &H)
	{
	const RectangularMesh &mesh = M.getMesh();

	int nx, ny, nz; mesh.getNumNodes(nx, ny, nz);
	double dx, dy, dz; mesh.getDelta(dx, dy, dz);

	return fs_spinhall(nx, ny, nz, dx, dy, dz, mu, J, t, Ms, M, Theta, H);
}

