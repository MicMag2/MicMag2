#ifndef CURRENTCALCULATOR_H
#define CURRENTCALCULATOR_H

#include "config.h"
#include "LinearSolver.h"
#include <vector>
#include <memory>

class CurrentCalculator
{
	public:
	enum CellType {
		NONCONDUCTIVE = 0,
		CONDUCTIVE = 1,
		CONTACT1 = 2,
		CONTACT2 = 3
	};
	enum OrderingTypes {
		NATURAL_ORDERING = 0,
		CONTACT2CONTACT_ORDERING = 1,
		TILE_ORDERING = 2,
	};

	private:
	int iteration;
        int nnf;
        bool AIsUpToDate;

	int nx;
	int ny;
	int nz;
	double dx;
	double dy;
	double dz;

	double tickTime[5];

        std::unique_ptr<LinearSolver> solver;
	std::vector<int> contacts;
	std::vector<double> B;
	std::vector<double> Phi;
	std::vector<int> nf;
	std::vector<int> nf2;
	std::vector<double> L;
	std::vector<int> Lc;
	std::vector<int> Ln;
	std::vector<double> U;
	std::vector<int> Uc;
	std::vector<int> Un;
	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<double> m_z;
	std::vector<double> sigma_0;
	std::vector<double> amr;
	std::vector<int> parallel;
	std::vector<CellType> cell_types;

        void init(int nx, int ny, int nz, double dx, double dy, double dz, std::vector<CellType>& cell_types, OrderingTypes ordering);
	void sort(std::vector<int>& list1, std::vector<int>& list2);
	double get_sigma(int kx, int ky, int kz);
	double get_sigma_xx(int kx, int ky, int kz);
	double get_sigma_yy(int kx, int ky, int kz);
	double get_sigma_zz(int kx, int ky, int kz);
	double get_sigma_xy(int kx, int ky, int kz);
	double get_sigma_xz(int kx, int ky, int kz);
	double get_sigma_yz(int kx, int ky, int kz);
	double get_sigma_yx(int kx, int ky, int kz);
	double get_sigma_zx(int kx, int ky, int kz);
	double get_sigma_zy(int kx, int ky, int kz);
        int index_of_neighbor(int kx, int ky, int kz, int lx, int ly, int lz);
        int index(int kx, int ky, int kz);
	void calculateRowOfA(std::vector<double>& A, int kx, int ky, int kz);
        void calculate_j(double* j_x, double* j_y, double* j_z, int kx, int ky, int kz);
	void prepareLinearSystem();
	void tick(int i);
	void tack(int i);

	public:
        CurrentCalculator(int nx, int ny, int nz, double dx, double dy, double dz, std::vector<CellType>& cell_types, OrderingTypes ordering);
	~CurrentCalculator();

	void set_m(double* m_x, double* m_y, double* m_z);
	void set_sigma(double* sigma);
	void set_amr(double* amr, int dim);
	void set_methodGMRES(int m);
	void set_methodAGMRES(int m_max, int m_min, int d, double cr_max, double cr_min);
	void set_preconditioning(int k, double omega);
	void solve(double rel);
        void get_j(double* j_x, double* j_y, double* j_z, double deltaU);
	void get_Phi(double* Phi, double deltaU);
        double get_R();
};

#endif

