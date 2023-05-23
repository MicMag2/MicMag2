#include "config.h"
#include "LinearSolver.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

#include <sys/time.h>
inline void LinearSolver::tick(int i)
{
	timeval ts;
	gettimeofday(&ts, 0);
	tickTime[i] -= (double)(ts.tv_sec) + (double)(ts.tv_usec) / 1000000;
	return;
}
inline void LinearSolver::tack(int i)
{
	timeval ts;
	gettimeofday(&ts, 0);
	tickTime[i] += (double)(ts.tv_sec) + (double)(ts.tv_usec) / 1000000;
	return;
}

inline void LinearSolver::givens(vector< vector<double> >& A, vector<double>& b, int x, int y)
{
  double r = sqrt(A[x][y]*A[x][y]+A[y][y]*A[y][y]);
  if (r != 0)
  {
    double s = A[x][y]/r;
    double c = A[y][y]/r;

    double tmp;
    for (unsigned int i = y; i < b.size(); i++)
    {
      tmp = A[y][i];
      A[y][i] = c*A[y][i]+s*A[x][i];
      A[x][i] = c*A[x][i]-s*tmp;
    }

    tmp = b[y];
    b[y] = c*b[y]+s*b[x];
    b[x] = c*b[x]-s*tmp;
  }

  return;
}

void LinearSolver::solve(double rel, vector<double>& L, vector<int>& Lc, vector<int>& Ln, vector<double>& U, vector<int>& Uc, vector<int>& Un, vector<double>& pB, vector<double>& pPhi, std::vector<int>& parallel)
{tick(0);
  this->L = &(L[0]);
  this->Lc = &(Lc[0]);
  this->Ln = &(Ln[0]);
  this->U = &(U[0]);
  this->Uc = &(Uc[0]);
  this->Un = &(Un[0]);
  this->parallel = &(parallel);

  double start_residual;
  double residual;

  iteration++;
  //take ssor initial value as initial guess for solver, not parallelized
//  if (iteration == 1) ssor(pB, Phi, 100);  //successive over-relaxation
  if (iteration == 1) for (int k = 0; k < size; k++) Phi[k] = 0;

  bool identical = true;
  for (int k = 0; k < size; k++)
    if (p[num_last_steps-1][k] < (1-1e-10)*Phi[k] || p[num_last_steps-1][k] > (1+1e-10)*Phi[k]) identical = false;

  if (!identical)
  {
    for (int i = 0; i < num_last_steps-1; i++)
      for (int k = 0; k < size; k++) p[i][k] = p[i+1][k];

    for (int k = 0; k < size; k++) p[num_last_steps-1][k] = Phi[k];

    if (iteration > num_last_steps)
    {
tick(3);
      vector< vector<double> > C(num_last_steps-1);
      for (int k = 0; k < num_last_steps-1; k++) C[k].resize(size);
      vector< vector<double> > C2(num_last_steps-1);
      for (int k = 0; k < num_last_steps-1; k++) C2[k].resize(num_last_steps-1);
      vector<double> t(size);
      vector<double> b(num_last_steps-1);

      for (int q = 0; q < num_last_steps-1; q++)
      {
        for (int l = 0; l < size; l++) t[l] = p[q][l]-p[num_last_steps-1][l];
        multA(t, C[q]);
      }

      for (int k = 0; k < num_last_steps-1; k++)
        for (int q = 0; q <= k; q++)
        {
          double c = scalar(C[k], C[q]);
          C2[k][q] = c;
          C2[q][k] = c;
        }

      multA(p[num_last_steps-1], t);
      for (int k = 0; k < size; k++) t[k] = pB[k]-t[k];
      for (int k = 0; k < num_last_steps-1; k++)
        b[k] = scalar(C[k], t);

      for (int l = 0; l < num_last_steps-1-1; l++)
        for (int k = num_last_steps-1-1; k > l; k--)
          givens(C2, b, k, l);

      bool singular = false;
      for (int k = 0; k < num_last_steps-1; k++)
        if (C2[k][k] == 0) singular = true;

      if (!singular)
      {
        for (int l = num_last_steps-1-1; l >= 0; l--)
        {
          for (int k = num_last_steps-1-1; k > l; k--) b[l] -= C2[l][k]*b[k];
          b[l] /= C2[l][l];
        }

        for (int i = 0; i < num_last_steps-1; i++)
          for (int k = 0; k < size; k++) Phi[k] += b[i]*(p[i][k]-p[num_last_steps-1][k]);
      }
tack(3);
    }
  }

  multPrecond(pB, B);
  start_residual = sqrt(scalar(B, B));
  residual = start_residual;



  if (method == METHOD_GMRES)
  {
    while (residual > rel*start_residual)
    {
      residual = iterate_GMRES(rel*start_residual, METHOD_GMRES, mKrylov);
    }
  }

  if (method == METHOD_AGMRES)
  {
    int cm = max_m_agmres;
    double cr = 1;

    while (residual > rel*start_residual)
    {
      cr = residual;
      residual = iterate_GMRES(rel*start_residual, METHOD_GMRES, cm);
      cr = residual/cr;

      if (cr > max_cr_agmres) cm = max_m_agmres;
      else if (cr > min_cr_agmres) cm -= d_agmres;
      if (cm < min_m_agmres) cm = max_m_agmres;
    }
  }

  for (int k = 0; k < size; k++) pPhi[k] = Phi[k];

tack(0);
  return;
}

inline void LinearSolver::multPrecond(vector<double>& x, vector<double>& y)
{
tick(2);
//  ssor(x, y, kssor);
  for (int k = 0; k < size; k++) y[k] = x[k];
tack(2);
  return;
}

inline void LinearSolver::ssor(vector<double>& x, vector<double>& y, int kssor)
{
  vector<double> t(size);
  vector<double> l(size);

  for (int k = size; k >= 0; k--)
  {
    t[k] = (2-omega)*omega*x[k];
    for (int n = Un[k]; n < Un[k+1]; n++) t[k] -= omega*U[n]*t[Uc[n]];
  }

  for (int k = 0; k < size; k++)
  {
    l[k] = 0;
    for (int n = Ln[k]; n < Ln[k+1]; n++) l[k] += omega*L[n]*y[Lc[n]];
    y[k] = t[k]-l[k];
  }


  for (int i = 0; i < kssor-1; i++)
  {
    for (int k = 0; k < size; k++) l[k] += (omega-1)*y[k];

    for (int k = size-1; k >= 0; k--)
    {
      t[k] = (2-omega)*(omega*x[k]-l[k]);
      for (int n = Un[k]; n < Un[k+1]; n++) t[k] -= omega*U[n]*t[Uc[n]];
    }

    for (int k = 0; k < size; k++)
    {
      double ls = l[k];
      l[k] = 0;
      for (int n = Ln[k]; n < Ln[k+1]; n++) l[k] += omega*L[n]*y[Lc[n]];
      y[k] = t[k]+ls-l[k];
    }
  }

  return;
}

inline void LinearSolver::multA(vector<double>& x, vector<double>& y)
{
tick(1);
  #ifdef USE_MULTI_THREAD
  #pragma omp parallel for
  #endif
  for (int k = 0; k < size; k++)
  {
    y[k] = x[k];
    for (int n = Ln[k]; n < Ln[k+1]; n++) y[k] += L[n]*x[Lc[n]];
    for (int n = Un[k]; n < Un[k+1]; n++) y[k] += U[n]*x[Uc[n]];
  }

tack(1);
  return;
}

inline double LinearSolver::scalar(vector<double>& x, vector<double>& y)
{
tick(4);
  double sum = 0;
  #ifdef USE_MULTI_THREAD
  #pragma omp parallel for reduction(+ : sum)
  #endif
  for (int k = 0; k < size; k++) sum += x[k]*y[k];
tack(4);
  return sum;
}

inline double LinearSolver::iterate_GMRES(double rel, int type, int m)
{
tick(5);
  vector<double> h((m+1)*(m+1));
  vector<double> d(m+1);
  vector<double> y(m+1);
  v.resize(m+1);
  for (int k = 0; k < m+1; k++) v[k].resize(size);
  vector<double> s(m+1);
  vector<double> c(m+1);
//  vector<double> t(size);

  int m2 = m;

    multA(Phi, v[0]);
//    multA(Phi, t);
//    multPrecond(t, v[0]);

    for (int k = 0; k < size; k++) v[0][k] = B[k]-v[0][k];

    d[0] = sqrt(scalar(v[0] ,v[0]));
    if (d[0] < rel)
    {
      tack(5);
      return d[0];
    }
    double residual = d[0];
    for (int k = 0; k < size; k++) v[0][k] /= d[0];

    for (int j = 0; j < (m+1)*(m+1); j++) h[j] = 0;

    for (int j = 0; j < m; j++)
    {
      multA(v[j], v[j+1]);
//      multA(v[j], t);
//      multPrecond(t, v[j+1]);

      for (int i = 0; i <= j; i++) h[(m+1)*i+j] = scalar(v[i], v[j+1]);

      #ifdef USE_MULTI_THREAD
      #pragma omp parallel for
      #endif
      for (int k = 0; k < size; k++)
        for (int i = 0; i <= j; i++)
          v[j+1][k] -= h[(m+1)*i+j]*v[i][k];

      h[(m+1)*(j+1)+j] = sqrt(scalar(v[j+1], v[j+1]));

      #ifdef USE_MULTI_THREAD
      #pragma omp parallel for
      #endif
      for (int k = 0; k < size; k++)
        v[j+1][k] /= h[(m+1)*(j+1)+j];

      for (int l = 0; l < j; l++)
      {
        double hu = h[(m+1)*l+j];
        double hl = h[(m+1)*(l+1)+j];
        h[(m+1)*l+j] = c[l]*hu+s[l]*hl;
        h[(m+1)*(l+1)+j] = c[l]*hl-s[l]*hu;
      }

      s[j] = h[(m+1)*(j+1)+j]/sqrt(h[(m+1)*j+j]*h[(m+1)*j+j]+h[(m+1)*(j+1)+j]*h[(m+1)*(j+1)+j]);
      c[j] = h[(m+1)*j+j]/sqrt(h[(m+1)*j+j]*h[(m+1)*j+j]+h[(m+1)*(j+1)+j]*h[(m+1)*(j+1)+j]);

      double hu = h[(m+1)*j+j];
      double hl = h[(m+1)*(j+1)+j];
      h[(m+1)*j+j] = c[j]*hu+s[j]*hl;
      h[(m+1)*(j+1)+j] = c[j]*hl-s[j]*hu;

      d[j+1] = -s[j]*d[j];
      d[j] = c[j]*d[j];

      if (fabs(d[j+1]) < rel)
      {
        m2 = j+1;
        break;
      }

    }

    //Minimize H y = d
    for (int i = m2-1; i >= 0; i--)
    {
      y[i] = d[i];
      for (int j = m2-1; j > i; j--) y[i] -= y[j]*h[(m+1)*i+j];
      y[i] /= h[(m+1)*i+i];
    }

    //Calculate new Phi
    #ifdef USE_MULTI_THREAD
    #pragma omp parallel for
    #endif
    for (int k = 0; k < size; k++)
      for (int i = 0; i < m2; i++)
        Phi[k] += v[i][k]*y[i];

  residual = fabs(d[m2]);
tack(5);
  return residual;
}

void LinearSolver::set_preconditioning(int k, double omega)
{
  this->kssor = k;
  this->omega = omega;
  return;
}

void LinearSolver::set_methodGMRES(int m)
{
  this->method = METHOD_GMRES;
  this->mKrylov = m;
  return;
}

void LinearSolver::set_methodAGMRES(int m_max, int m_min, int d, double cr_max, double cr_min)
{
  this->method = METHOD_AGMRES;
  this->max_m_agmres = m_max;
  this->min_m_agmres = m_min;
  this->max_cr_agmres = cr_max;
  this->min_cr_agmres = cr_min;
  this->d_agmres = d;
  return;
}

LinearSolver::LinearSolver(int size)
{
  this->method = METHOD_AGMRES;
  this->omega = 1.85;
  this->kssor = 1;
  this->mKrylov = 25;
  this->max_cr_agmres = 0.99;
  this->min_cr_agmres = 0.3;
  this->d_agmres = 3;
  this->max_m_agmres = 25;
  this->min_m_agmres = 1;
  this->iteration = 0;

  this->size = size;
  B.resize(size);

  num_last_steps = 4;
  p.resize(num_last_steps);	
  for (int k = 0; k < num_last_steps; k++) p[k].resize(size);
  Phi.resize(size);
  for (int k = 0; k < size; k++)
  {
    this->Phi[k] = 0;
  }

  for (int k = 0; k < 6; k++)
  {
    tickTime[k] = 0;
  }
}

LinearSolver::~LinearSolver()
{
  cerr << "LinearSolver: The calculation of " << iteration << " solutions took " << tickTime[0] << " s (average per solution: " << tickTime[0]/iteration << " s)\n";
  cerr << "LinearSolver: Time for matrix-vector multiplications: " << tickTime[1] << " s (average per solution: " << tickTime[1]/iteration << " s)\n";
  cerr << "LinearSolver: Time for applying the preconditioner: " << tickTime[2] << " s (average per solution: " << tickTime[2]/iteration << " s)\n";
  cerr << "LinearSolver: Time for performing the extrapolations: " << tickTime[3] << " s (average per solution: " << tickTime[3]/iteration << " s)\n";
  cerr << "LinearSolver: Time for calculating the scalar products: " << tickTime[4] << " s (average per solution: " << tickTime[4]/iteration << " s)\n";
  cerr << "LinearSolver: Time for performing the GMRES iterations: " << tickTime[5] << " s (average per solution: " << tickTime[5]/iteration << " s)\n";
}

