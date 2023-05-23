#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "config.h"

#include <vector>

class LinearSolver
{
	private:
	const static int METHOD_GMRES = 1;
	const static int METHOD_AGMRES = 2;

        int size;
	int iteration;
        int method;
	double omega;
	int kssor;
	int mKrylov;
	int max_m_agmres;
	int min_m_agmres;
	double max_cr_agmres;
	double min_cr_agmres;
	double d_agmres;
	int num_last_steps;

	double tickTime[6];

	std::vector<double> B;
	std::vector<double> Phi;
	double* L;
	int* Lc;
	int* Ln;
	double* U;
	int* Uc;
	int* Un;
	std::vector<int>* parallel;
	std::vector<double> S;
	std::vector<int> Sc;
	std::vector<int> Sn;
	std::vector< std::vector<double> > p;

	std::vector< std::vector<double> > v;

        double iterate_GMRES(double rel, int type, int m);
        void multA(std::vector<double>& x, std::vector<double>& y);
        void multPrecond(std::vector<double>& x, std::vector<double>& y);
	void ssor(std::vector<double>& x, std::vector<double>& y, int kssor);
	double scalar(std::vector<double>& x, std::vector<double>& y);
	void givens(std::vector< std::vector<double> >& C, std::vector<double>& b, int x, int y);
	void tick(int i);
	void tack(int i);

	public:
        LinearSolver(int size);
        ~LinearSolver();

	void solve(double rel, std::vector<double>& L, std::vector<int>& Lc, std::vector<int>& Ln, std::vector<double>& U, std::vector<int>& Uc, std::vector<int>& Un, std::vector<double>& pB, std::vector<double>& pPhi, std::vector<int>& parallel);
	void set_methodGMRES(int m);
	void set_methodAGMRES(int m_max, int m_min, int d, double cr_max, double cr_min);
	void set_preconditioning(int k, double omega);
}
;

#endif

