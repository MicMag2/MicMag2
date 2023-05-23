#include "config.h"
#include "CurrentCalculator.h"
#include "LinearSolver.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;


#include <sys/time.h>
inline void CurrentCalculator::tick(int i)
{
	timeval ts;
	gettimeofday(&ts, 0);
	tickTime[i] -= (double)(ts.tv_sec) + (double)(ts.tv_usec) / 1000000;
	return;
}
inline void CurrentCalculator::tack(int i)
{
	timeval ts;
	gettimeofday(&ts, 0);
	tickTime[i] += (double)(ts.tv_sec) + (double)(ts.tv_usec) / 1000000;
	return;
}

inline void CurrentCalculator::sort(std::vector<int>& list1, std::vector<int>& list2)
{
  int size = list1.size();

  if (size <= 1) return;

  int size1 = size/2;
  int size2 = size-size1;

  std::vector<int> list11(size1);
  std::vector<int> list12(size2);
  std::vector<int> list21(size1);
  std::vector<int> list22(size2);
  for (int i = 0; i < size1; i++)
  {
    list11[i] = list1[i];
    list21[i] = list2[i];
  }
  for (int i = 0; i < size2; i++)
  {
    list12[i] = list1[i+size1];
    list22[i] = list2[i+size1];
  }

  sort(list11, list21);
  sort(list12, list22);

  int i = 0;
  int j = 0;
  while (i < size1 || j < size2)
  {
    if (!(i == size1) && (j == size2 || list11[i] < list12[j]))
    {
      list1[i+j] = list11[i];
      list2[i+j] = list21[i];
      i++;
    }
    else
    {
      list1[i+j] = list12[j];
      list2[i+j] = list22[j];
      j++;
    }
  }

  return;
}

inline double CurrentCalculator::get_sigma(int kx, int ky, int kz)
{
  if (kx < 0 || ky < 0 || kz < 0 || kx >= nx || ky >= ny || kz >= nz) return 0;
  return sigma_0[ny*nx*kz+nx*ky+kx];
}

inline double CurrentCalculator::get_sigma_xx(int kx, int ky, int kz)
{
  return get_sigma(kx, ky, kz)*(1+amr[ny*nx*kz+nx*ky+kx]*(m_x[ny*nx*kz+nx*ky+kx]*m_x[ny*nx*kz+nx*ky+kx]-1.0/3.0));
}

inline double CurrentCalculator::get_sigma_yy(int kx, int ky, int kz)
{
  return get_sigma(kx, ky, kz)*(1+amr[ny*nx*kz+nx*ky+kx]*(m_y[ny*nx*kz+nx*ky+kx]*m_y[ny*nx*kz+nx*ky+kx]-1.0/3.0));
}

inline double CurrentCalculator::get_sigma_zz(int kx, int ky, int kz)
{
  return get_sigma(kx, ky, kz)*(1+amr[ny*nx*kz+nx*ky+kx]*(m_z[ny*nx*kz+nx*ky+kx]*m_z[ny*nx*kz+nx*ky+kx]-1.0/3.0));
}

inline double CurrentCalculator::get_sigma_xy(int kx, int ky, int kz)
{
  return get_sigma(kx, ky, kz)*amr[ny*nx*kz+nx*ky+kx]*m_x[ny*nx*kz+nx*ky+kx]*m_y[ny*nx*kz+nx*ky+kx];
}

inline double CurrentCalculator::get_sigma_yx(int kx, int ky, int kz)
{
  return get_sigma(kx, ky, kz)*amr[ny*nx*kz+nx*ky+kx]*m_x[ny*nx*kz+nx*ky+kx]*m_y[ny*nx*kz+nx*ky+kx];
}

inline double CurrentCalculator::get_sigma_yz(int kx, int ky, int kz)
{
  return get_sigma(kx, ky, kz)*amr[ny*nx*kz+nx*ky+kx]*m_y[ny*nx*kz+nx*ky+kx]*m_z[ny*nx*kz+nx*ky+kx];
}

inline double CurrentCalculator::get_sigma_zy(int kx, int ky, int kz)
{
  return get_sigma(kx, ky, kz)*amr[ny*nx*kz+nx*ky+kx]*m_y[ny*nx*kz+nx*ky+kx]*m_z[ny*nx*kz+nx*ky+kx];
}

inline double CurrentCalculator::get_sigma_xz(int kx, int ky, int kz)
{
  return get_sigma(kx, ky, kz)*amr[ny*nx*kz+nx*ky+kx]*m_x[ny*nx*kz+nx*ky+kx]*m_z[ny*nx*kz+nx*ky+kx];
}

inline double CurrentCalculator::get_sigma_zx(int kx, int ky, int kz)
{
  return get_sigma(kx, ky, kz)*amr[ny*nx*kz+nx*ky+kx]*m_x[ny*nx*kz+nx*ky+kx]*m_z[ny*nx*kz+nx*ky+kx];
}

inline int CurrentCalculator::index_of_neighbor(int kx, int ky, int kz, int lx, int ly, int lz)
{
  return 27*(nf2[ny*nx*lz+nx*ly+lx])+13+(kx-lx)+3*(ky-ly)+9*(kz-lz);
}

inline int CurrentCalculator::index(int kx, int ky, int kz)
{
  return ny*nx*kz+nx*ky+kx;
}

inline void CurrentCalculator::calculateRowOfA(vector<double>& A, int kx, int ky, int kz)
{
  double sigma_bar;
  int k = ny*nx*kz+nx*ky+kx;

  if (nf2[k] == -1) return;

  int c = 27*nf2[k]+13;

  for (int i = -1; i <= 1; i += 2)
  {
    if (get_sigma(kx+i, ky, kz) != 0)
    {
      sigma_bar = 2*(get_sigma_xx(kx, ky, kz)*get_sigma_xx(kx+i, ky, kz))/(get_sigma_xx(kx, ky, kz)+get_sigma_xx(kx+i, ky, kz));

      for (int p = 0; p <= 1; p++)
      {
        A[c+1*i*p] += (2*p-1)*sigma_bar*dy*dz/dx;

        for (int l = -1; l <= 1; l += 2)
        {
          double f = i*l*sigma_bar/(2*get_sigma(kx+i*p, ky, kz));

          for (int q = 0; q <= 1; q++)
          {
            A[c+1*i*p+3*l*q] += f*dz*get_sigma_xy(kx+i*p, ky, kz)*get_sigma(kx+i*p, ky+l*q, kz)/(get_sigma(kx+i*p, ky+l, kz)+get_sigma(kx+i*p, ky, kz));

            A[c+1*i*p+9*l*q] += f*dy*get_sigma_xz(kx+i*p, ky, kz)*get_sigma(kx+i*p, ky, kz+l*q)/(get_sigma(kx+i*p, ky, kz+l)+get_sigma(kx+i*p, ky, kz));
          }
        }
      }
    }

    if (get_sigma(kx, ky+i, kz) != 0)
    {
      sigma_bar = 2*(get_sigma_yy(kx, ky, kz)*get_sigma_yy(kx, ky+i, kz))/(get_sigma_yy(kx, ky, kz)+get_sigma_yy(kx, ky+i, kz));

      for (int p = 0; p <= 1; p++)
      {
        A[c+3*i*p] += (2*p-1)*sigma_bar*dx*dz/dy;

        for (int l = -1; l <= 1; l += 2)
        {
          double f = i*l*sigma_bar/(2*get_sigma(kx, ky+i*p, kz));

          for (int q = 0; q <= 1; q++)
          {
            A[c+3*i*p+1*l*q] += f*dz*get_sigma_yx(kx, ky+i*p, kz)*get_sigma(kx+l*q, ky+i*p, kz)/(get_sigma(kx+l, ky+i*p, kz)+get_sigma(kx, ky+i*p, kz));

            A[c+3*i*p+9*l*q] += f*dx*get_sigma_yz(kx, ky+i*p, kz)*get_sigma(kx, ky+i*p, kz+l*q)/(get_sigma(kx, ky+i*p, kz+l)+get_sigma(kx, ky+i*p, kz));
          }
        }
      }
    }

    if (get_sigma(kx, ky, kz+i) != 0)
    {
      sigma_bar = 2*(get_sigma_zz(kx, ky, kz)*get_sigma_zz(kx, ky, kz+i))/(get_sigma_zz(kx, ky, kz)+get_sigma_zz(kx, ky, kz+i));

      for (int p = 0; p <= 1; p++)
      {
        A[c+9*i*p] += (2*p-1)*sigma_bar*dx*dy/dz;

        for (int l = -1; l <= 1; l += 2)
        {
          double f = i*l*sigma_bar/(2*get_sigma(kx, ky, kz+i*p));

          for (int q = 0; q <= 1; q++)
          {
            A[c+9*i*p+1*l*q] += f*dy*get_sigma_zx(kx, ky, kz+i*p)*get_sigma(kx+l*q, ky, kz+i*p)/(get_sigma(kx+l, ky, kz+i*p)+get_sigma(kx, ky, kz+i*p));

            A[c+9*i*p+3*l*q] += f*dx*get_sigma_zy(kx, ky, kz+i*p)*get_sigma(kx, ky+l*q, kz+i*p)/(get_sigma(kx, ky+l, kz+i*p)+get_sigma(kx, ky, kz+i*p));
          }
        }
      }
    }
  }
 
  return;
}

inline void CurrentCalculator::calculate_j(double* j_x, double* j_y, double* j_z, int kx, int ky, int kz)
{
  if (get_sigma(kx, ky, kz) == 0)
  {
    j_x = 0;
    j_y = 0;
    j_z = 0;
    return;
  }	

  double sigma_bar;
  int k = ny*nx*kz+nx*ky+kx;

  vector<double> Ax(27, 0);
  vector<double> Ay(27, 0);
  vector<double> Az(27, 0);
  int c = 13;

  for (int i = -1; i <= 1; i += 2)
  {
    if (get_sigma(kx+i, ky, kz) != 0)
    {
      sigma_bar = 2*(get_sigma_xx(kx, ky, kz)*get_sigma_xx(kx+i, ky, kz))/(get_sigma_xx(kx, ky, kz)+get_sigma_xx(kx+i, ky, kz));

      for (int p = 0; p <= 1; p++)
      {
        Ax[c+1*i*p] += i*(2*p-1)*sigma_bar*dy*dz/dx;

        for (int l = -1; l <= 1; l += 2)
        {
          double f = l*sigma_bar/(2*get_sigma(kx+i*p, ky, kz));

          for (int q = 0; q <= 1; q++)
          {
            Ax[c+1*i*p+3*l*q] += f*dz*get_sigma_xy(kx+i*p, ky, kz)*get_sigma(kx+i*p, ky+l*q, kz)/(get_sigma(kx+i*p, ky+l, kz)+get_sigma(kx+i*p, ky, kz));

            Ax[c+1*i*p+9*l*q] += f*dy*get_sigma_xz(kx+i*p, ky, kz)*get_sigma(kx+i*p, ky, kz+l*q)/(get_sigma(kx+i*p, ky, kz+l)+get_sigma(kx+i*p, ky, kz));
          }
        }
      }
    }

    if (get_sigma(kx, ky+i, kz) != 0)
    {
      sigma_bar = 2*(get_sigma_yy(kx, ky, kz)*get_sigma_yy(kx, ky+i, kz))/(get_sigma_yy(kx, ky, kz)+get_sigma_yy(kx, ky+i, kz));

      for (int p = 0; p <= 1; p++)
      {
        Ay[c+3*i*p] += i*(2*p-1)*sigma_bar*dx*dz/dy;

        for (int l = -1; l <= 1; l += 2)
        {
          double f = l*sigma_bar/(2*get_sigma(kx, ky+i*p, kz));

          for (int q = 0; q <= 1; q++)
          {
            Ay[c+3*i*p+1*l*q] += f*dz*get_sigma_yx(kx, ky+i*p, kz)*get_sigma(kx+l*q, ky+i*p, kz)/(get_sigma(kx+l, ky+i*p, kz)+get_sigma(kx, ky+i*p, kz));

            Ay[c+3*i*p+9*l*q] += f*dx*get_sigma_yz(kx, ky+i*p, kz)*get_sigma(kx, ky+i*p, kz+l*q)/(get_sigma(kx, ky+i*p, kz+l)+get_sigma(kx, ky+i*p, kz));
          }
        }
      }
    }

    if (get_sigma(kx, ky, kz+i) != 0)
    {
      sigma_bar = 2*(get_sigma_zz(kx, ky, kz)*get_sigma_zz(kx, ky, kz+i))/(get_sigma_zz(kx, ky, kz)+get_sigma_zz(kx, ky, kz+i));

      for (int p = 0; p <= 1; p++)
      {
        Az[c+9*i*p] += i*(2*p-1)*sigma_bar*dx*dy/dz;

        for (int l = -1; l <= 1; l += 2)
        {
          double f = l*sigma_bar/(2*get_sigma(kx, ky, kz+i*p));

          for (int q = 0; q <= 1; q++)
          {
            Az[c+9*i*p+1*l*q] += f*dy*get_sigma_zx(kx, ky, kz+i*p)*get_sigma(kx+l*q, ky, kz+i*p)/(get_sigma(kx+l, ky, kz+i*p)+get_sigma(kx, ky, kz+i*p));

            Az[c+9*i*p+3*l*q] += f*dx*get_sigma_zy(kx, ky, kz+i*p)*get_sigma(kx, ky+l*q, kz+i*p)/(get_sigma(kx, ky+l, kz+i*p)+get_sigma(kx, ky, kz+i*p));
          }
        }
      }
    }
  }

  int i = 0;
  for (int dz = -1; dz <= 1; dz++)
  {
    for (int dy = -1; dy <= 1; dy++)
    {
      for (int dx = -1; dx <= 1; dx++)
      {
        if (Ax[i] != 0) j_x[k] += Ax[i]*Phi[nx*ny*(kz+dz)+nx*(ky+dy)+kx+dx];
        if (Ay[i] != 0) j_y[k] += Ay[i]*Phi[nx*ny*(kz+dz)+nx*(ky+dy)+kx+dx];
        if (Az[i] != 0) j_z[k] += Az[i]*Phi[nx*ny*(kz+dz)+nx*(ky+dy)+kx+dx];
        i++;
      }
    }
  }

  return;
}





inline void CurrentCalculator::prepareLinearSystem()
{
  vector<double> A(27*nnf);
  for (int k = 0; k < 27*nnf; k++) A[k] = 0;

  for (int kz = 0; kz < nz; kz++)
    for (int ky = 0; ky < ny; ky++)
      #ifdef USE_MULTI_THREAD
      #pragma omp parallel for
      #endif
      for (int kx = 0; kx < nx; kx++)
        calculateRowOfA(A, kx, ky, kz);

  #ifdef USE_MULTI_THREAD
  #pragma omp parallel for
  #endif
  for (int k = 0; k < nnf; k++)
  {
    double diag = A[27*k+13];
    for (int n = 0; n < 27; n++) A[27*k+n] /= diag;

    B[k] = 0;

    Ln[k+1] = 0;
    Un[k+1] = 0;

    int n = 0;
    for (int nkz = -1; nkz <= 1; nkz++)
    {
      for (int nky = -1; nky <= 1; nky++)
      {
        for (int nkx = -1; nkx <= 1; nkx++)
        {
          int c = nf[k]+ny*nx*nkz+nx*nky+nkx;
          if (A[27*k+n] != 0 && nf2[c] == -1)
          {
            B[k] -= A[27*k+n]*contacts[c];
            A[27*k+n] = 0;
          }

          if (A[27*k+n] != 0 && k > nf2[c]) Ln[k+1]++;
          if (A[27*k+n] != 0 && k < nf2[c]) Un[k+1]++;

          n++;
        }
      }
    }
  }

  for (int k = 0; k < nnf; k++)
  {
    Ln[k+1] += Ln[k];
    Un[k+1] += Un[k];
  }

  L.resize(Ln[nnf]);
  Lc.resize(Ln[nnf]);
  U.resize(Un[nnf]);
  Uc.resize(Un[nnf]);

  #ifdef USE_MULTI_THREAD
  #pragma omp parallel for
  #endif
  for (int k = 0; k < nnf; k++)
  {
    int nL = Ln[k];
    int nU = Un[k];
    int n = 0;
    for (int nkz = -1; nkz <= 1; nkz++)
    {
      for (int nky = -1; nky <= 1; nky++)
      {
        for (int nkx = -1; nkx <= 1; nkx++)
        {
          int c = nf[k]+ny*nx*nkz+nx*nky+nkx;

          if (A[27*k+n] != 0 && k > nf2[c])
          {
            L[nL] = A[27*k+n];
            Lc[nL] = nf2[c];
            nL++;
          }
          if (A[27*k+n] != 0 && k < nf2[c])
          {
            U[nU] = A[27*k+n];
            Uc[nU] = nf2[c];
            nU++;
          }
          n++;

        }
      }
    }
  }

  return;
}







void CurrentCalculator::set_m(double* m_x, double* m_y, double* m_z)
{
  for (int k = 0; k < nx*ny*nz; k++)
  {
    double norm = sqrt(m_x[k]*m_x[k]+m_y[k]*m_y[k]+m_z[k]*m_z[k]);
    if (norm != 0)
    {
      this->m_x[k] = m_x[k]/norm;
      this->m_y[k] = m_y[k]/norm;
      this->m_z[k] = m_z[k]/norm;
    }
    else
    {
      this->m_x[k] = 0;
      this->m_y[k] = 0;
      this->m_z[k] = 0;
    }
  }
  AIsUpToDate = false;
  return;
}

void CurrentCalculator::set_sigma(double* sigma)
{
  for (int kz = 0; kz < nz; kz++)
  {
    for (int ky = 0; ky < ny; ky++)
    {
      for (int kx = 0; kx < nx; kx++)
      {
        if (sigma[ny*nx*kz+nx*ky+kx] == 0 && cell_types[ny*nx*kz+nx*ky+kx] != NONCONDUCTIVE)
        {
          cerr << "sigma = 0 found where conductive cell was expected\n";
          cerr << "position: " << kx << " " << ky << " " << kz << "\n";
          exit(0);
        }
        if (sigma[ny*nx*kz+nx*ky+kx] != 0 && cell_types[ny*nx*kz+nx*ky+kx] == NONCONDUCTIVE)
        {
          cerr << "sigma > 0 found where non-conductive cell was expected\n";
          cerr << "position: " << kx << " " << ky << " " << kz << "\n";
          exit(0);
        }
        this->sigma_0[ny*nx*kz+nx*ky+kx] = sigma[ny*nx*kz+nx*ky+kx];
      }
    }
  }
  AIsUpToDate = false;
  return;
}

void CurrentCalculator::set_amr(double* amr, int dim)
{
  for (int k = 0; k < nx*ny*nz; k++) this->amr[k] = -3*dim*amr[k]/(3*dim+(2*dim-3)*amr[k]);
  AIsUpToDate = false;
  return;
}

void CurrentCalculator::set_preconditioning(int k, double omega)
{
  solver->set_preconditioning(k, omega);
  return;
}

void CurrentCalculator::set_methodGMRES(int m)
{
  solver->set_methodGMRES(m);
  return;
}

void CurrentCalculator::set_methodAGMRES(int m_max, int m_min, int d, double cr_max, double cr_min)
{
  solver->set_methodAGMRES(m_max, m_min, d, cr_max, cr_min);
  return;
}

void CurrentCalculator::solve(double rel)
{
tick(0);
  iteration++;
  if (!AIsUpToDate)
  {
tick(1);
    prepareLinearSystem();
    AIsUpToDate = true;
tack(1);
  }

  vector<double> Phi2(nnf);
  for (int l = 0; l < nnf; l++) Phi2[l] = Phi[nf[l]];

tick(2);
  solver->solve(rel, L, Lc, Ln, U, Uc, Un, B, Phi2, parallel);
tack(2);

  for (int k = 0; k < nx*ny*nz; k++) Phi[k] = contacts[k];
  for (int l = 0; l < nnf; l++) Phi[nf[l]] = Phi2[l];
tack(0);
//cout << iteration << " " << tickTime[0] << "\n";
  return;
}

double CurrentCalculator::get_R()
{
tick(0);
tick(3);
  double integral = 0;

  //x surfaces
  for (int kz = 0; kz < nz; kz++)
  {
    for (int ky = 0; ky < ny; ky++)
    {
      for (int kx = 1; kx < nx; kx++)
      {
        if (get_sigma(kx-1, ky, kz) != 0 && get_sigma(kx, ky, kz) != 0 && contacts[index(kx-1, ky, kz)] != contacts[index(kx, ky, kz)])
        {
          //calculate conductivity
          double sigma_bar = 2*(get_sigma_xx(kx-1, ky, kz)*get_sigma_xx(kx, ky, kz))/
                           (get_sigma_xx(kx-1, ky, kz)+get_sigma_xx(kx, ky, kz));

          //add current through the surface to neigbouring cells
          for (int i = -1; i <= 0; i++)//add to left (i=-1) or right (i=0) cell
          {
            for (int p = -1; p <= 0; p++)
              integral += contacts[index(kx+i, ky, kz)]*Phi[index(kx+p, ky, kz)]*(2*p+1)*(2*i+1)*sigma_bar*dy*dz/dx;

            for (int l = -1; l <= 0; l++)
              for (int p = -1; p <= 1; p += 2)
                for (int q = 0; q <= 1; q++)
                  if (get_sigma(kx+l, ky+p*q, kz) != 0)
                    integral += contacts[index(kx+i, ky, kz)]*Phi[index(kx+l, ky+p*q, kz)]*p*(2*i+1)*sigma_bar*get_sigma_xy(kx+l, ky, kz)/get_sigma_xx(kx+l, ky, kz)*dz/2*get_sigma(kx+l, ky+p*q, kz)/(get_sigma(kx+l, ky, kz)+get_sigma(kx+l, ky+p, kz));

            for (int l = -1; l <= 0; l++)
              for (int p = -1; p <= 1; p += 2)
                for (int q = 0; q <= 1; q++)
                  if (get_sigma(kx+l, ky, kz+p*q) != 0)
                    integral += contacts[index(kx+i, ky, kz)]*Phi[index(kx+l, ky, kz+p*q)]*p*(2*i+1)*sigma_bar*get_sigma_xz(kx+l, ky, kz)/get_sigma_xx(kx+l, ky, kz)*dy/2*get_sigma(kx+l, ky, kz+p*q)/(get_sigma(kx+l, ky, kz)+get_sigma(kx+l, ky, kz+p));
          }
        }
      }
    }
  }

  //y surfaces
  for (int kz = 0; kz < nz; kz++)
  {
    for (int ky = 1; ky < ny; ky++)
    {
      for (int kx = 0; kx < nx; kx++)
      {
        if (get_sigma(kx, ky-1, kz) != 0 && get_sigma(kx, ky, kz) != 0 && contacts[index(kx, ky-1, kz)] != contacts[index(kx, ky, kz)])
        {
          //calculate conductivity
          double sigma_bar = 2*(get_sigma_yy(kx, ky-1, kz)*get_sigma_yy(kx, ky, kz))/
                           (get_sigma_yy(kx, ky-1, kz)+get_sigma_yy(kx, ky, kz));

          //add current through the surface to neigbouring cells
          for (int i = -1; i <= 0; i++)//add to left (i=-1) or right (i=0) cell
          {
            for (int p = -1; p <= 0; p++)
              integral += contacts[index(kx, ky+i, kz)]*Phi[index(kx, ky+p, kz)]*(2*p+1)*(2*i+1)*sigma_bar*dx*dz/dy;

            for (int l = -1; l <= 0; l++)
              for (int p = -1; p <= 1; p += 2)
                for (int q = 0; q <= 1; q++)
                  if (get_sigma(kx+p*q, ky+l, kz) != 0)
                    integral += contacts[index(kx, ky+i, kz)]*Phi[index(kx+p*q, ky+l, kz)]*p*(2*i+1)*sigma_bar*get_sigma_yx(kx, ky+l, kz)/get_sigma_yy(kx, ky+l, kz)*dz/2*get_sigma(kx+p*q, ky+l, kz)/(get_sigma(kx, ky+l, kz)+get_sigma(kx+p, ky+l, kz));

            for (int l = -1; l <= 0; l++)
              for (int p = -1; p <= 1; p += 2)
                for (int q = 0; q <= 1; q++)
                  if (get_sigma(kx, ky+l, kz+p*q) != 0)
                    integral += contacts[index(kx, ky+i, kz)]*Phi[index(kx, ky+l, kz+p*q)]*p*(2*i+1)*sigma_bar*get_sigma_yz(kx, ky+l, kz)/get_sigma_yy(kx, ky+l, kz)*dx/2*get_sigma(kx, ky+l, kz+p*q)/(get_sigma(kx, ky+l, kz)+get_sigma(kx, ky+l, kz+p));
          }
        }
      }
    }
  }

  //z surfaces
  for (int kz = 1; kz < nz; kz++)
  {
    for (int ky = 0; ky < ny; ky++)
    {
      for (int kx = 0; kx < nx; kx++)
      {
        if (get_sigma(kx, ky, kz-1) != 0 && get_sigma(kx, ky, kz) != 0 && contacts[index(kx, ky, kz-1)] != contacts[index(kx, ky, kz)])
        {
          //calculate conductivity
          double sigma_bar = 2*(get_sigma_zz(kx, ky, kz-1)*get_sigma_zz(kx, ky, kz))/
                           (get_sigma_zz(kx, ky, kz-1)+get_sigma_zz(kx, ky, kz));

          //add current through the surface to neigbouring cells
          for (int i = -1; i <= 0; i++)//add to left (i=-1) or right (i=0) cell
          {
            for (int p = -1; p <= 0; p++)
              integral += contacts[index(kx, ky, kz+i)]*Phi[index(kx, ky, kz+p)]*(2*p+1)*(2*i+1)*sigma_bar*dx*dy/dz;

            for (int l = -1; l <= 0; l++)
              for (int p = -1; p <= 1; p += 2)
                for (int q = 0; q <= 1; q++)
                  if (get_sigma(kx+p*q, ky, kz+l) != 0)
                    integral += contacts[index(kx, ky, kz+i)]*Phi[index(kx+p*q, ky, kz+l)]*p*(2*i+1)*sigma_bar*get_sigma_zx(kx, ky, kz+l)/get_sigma_zz(kx, ky, kz+l)*dy/2*get_sigma(kx+p*q, ky, kz+l)/(get_sigma(kx, ky, kz+l)+get_sigma(kx+p, ky, kz+l));

            for (int l = -1; l <= 0; l++)
              for (int p = -1; p <= 1; p += 2)
                for (int q = 0; q <= 1; q++)
                  if (get_sigma(kx, ky+p*q, kz+l) != 0)
                    integral += contacts[index(kx, ky, kz+i)]*Phi[index(kx, ky+p*q, kz+l)]*p*(2*i+1)*sigma_bar*get_sigma_zy(kx, ky, kz+l)/get_sigma_zz(kx, ky, kz+l)*dx/2*get_sigma(kx, ky+p*q, kz+l)/(get_sigma(kx, ky, kz+l)+get_sigma(kx, ky+p, kz+l));
          }
        }
      }
    }
  }
tack(0);
tack(3);
  return 4/integral;
}

void CurrentCalculator::get_j(double* j_x, double* j_y, double* j_z, double deltaU)
{
tick(0);
tick(4);
  for (int k = 0; k < nx*ny*nz; k++)
  {
    j_x[k] = 0;
    j_y[k] = 0;
    j_z[k] = 0;
  }

  for (int kz = 0; kz < nz; kz++)
    for (int ky = 0; ky < ny; ky++)
      #ifdef USE_MULTI_THREAD
      #pragma omp parallel for
      #endif
      for (int kx = 0; kx < nx; kx++)
        calculate_j(j_x, j_y, j_z, kx, ky, kz);

  for (int k = 0; k < nx*ny*nz; k++)
  {
    j_x[k] *= deltaU/(4*dy*dz);
    j_y[k] *= deltaU/(4*dx*dz);
    j_z[k] *= deltaU/(4*dx*dy);
  }

tack(0);
tack(4);
  return;
}

void CurrentCalculator::get_Phi(double* Phi, double deltaU)
{
  for (int k = 0; k < nx*ny*nz; k++) Phi[k] = this->Phi[k]*deltaU/2;
  return;
}





CurrentCalculator::CurrentCalculator(int nx, int ny, int nz, double dx, double dy, double dz, std::vector<CellType>& cell_types, OrderingTypes ordering)
{
  init(nx, ny, nz, dx, dy, dz, cell_types, ordering);
}

void CurrentCalculator::init(int nx, int ny, int nz, double dx, double dy, double dz, std::vector<CellType>& cell_types, OrderingTypes ordering)
{
  this->cell_types.resize(nx*ny*nz);
  for (int k = 0; k < nx*ny*nz; k++) this->cell_types[k] = cell_types[k];

  this->contacts.resize(nx*ny*nz);
  for (int k = 0; k < nx*ny*nz; k++)
  {
    if (cell_types[k] == CONTACT1)
      this->contacts[k] = 1;
    else if (cell_types[k] == CONTACT2)
      this->contacts[k] = -1;
    else
      this->contacts[k] = 0;
  }

  parallel.resize(1);
  parallel[0] = 0;
  nf.resize(nx*ny*nz);
  nf2.resize(nx*ny*nz);
  nnf = 0;

  if (ordering == NATURAL_ORDERING)
  {
    for (int k = 0; k < nx*ny*nz; k++)
    {
      if (cell_types[k] == CONDUCTIVE)
      {
        nf[nnf] = k;
        nf2[k] = nnf;
        nnf++;
      }
      else
      {
        nf2[k] = -1;
      }
    }
    parallel.push_back(nnf);
    parallel.push_back(nnf);
    parallel.push_back(nnf);
    parallel.push_back(nnf);
    parallel.push_back(nnf);
  }
  else if (ordering == CONTACT2CONTACT_ORDERING)
  {
    vector<int> node_array1(nx*ny*nz, 0);
    vector<int> node_array2(nx*ny*nz, 0);

    {
    vector<int> old_nodes_x(0);
    vector<int> old_nodes_y(0);
    vector<int> old_nodes_z(0);
    vector<int> new_nodes_x(0);
    vector<int> new_nodes_y(0);
    vector<int> new_nodes_z(0);

    for (int kz = 0; kz < nz; kz++)
      for (int ky = 0; ky < ny; ky++)
        for (int kx = 0; kx < nx; kx++)
        {
          node_array1[nx*ny*kz+nx*ky+kx] = -1;
          if (contacts[nx*ny*kz+nx*ky+kx] == 1)
          {
            old_nodes_x.push_back(kx);
            old_nodes_y.push_back(ky);
            old_nodes_z.push_back(kz);
          }
        }

    for (int i = 0; i < nx*ny*nz; i++)
    {
      for (std::vector<int>::size_type l = 0; l < old_nodes_x.size(); l++)
      {
        int kx = old_nodes_x[l];
        int ky = old_nodes_y[l];
        int kz = old_nodes_z[l];

        if (kx < nx-1 && node_array1[nx*ny*kz+nx*ky+(kx+1)] == -1 && cell_types[nx*ny*kz+nx*ky+(kx+1)] != NONCONDUCTIVE)
        {
          node_array1[nx*ny*kz+nx*ky+(kx+1)] = i;
          new_nodes_x.push_back(kx+1);
          new_nodes_y.push_back(ky);
          new_nodes_z.push_back(kz);
        }

        if (kx > 0 && node_array1[nx*ny*kz+nx*ky+(kx-1)] == -1 && cell_types[nx*ny*kz+nx*ky+(kx-1)] != NONCONDUCTIVE)
        {
          node_array1[nx*ny*kz+nx*ky+(kx-1)] = i;
          new_nodes_x.push_back(kx-1);
          new_nodes_y.push_back(ky);
          new_nodes_z.push_back(kz);
        }

        if (ky < ny-1 && node_array1[nx*ny*kz+nx*(ky+1)+kx] == -1 && cell_types[nx*ny*kz+nx*(ky+1)+kx] != NONCONDUCTIVE)
        {
          node_array1[nx*ny*kz+nx*(ky+1)+kx] = i;
          new_nodes_x.push_back(kx);
          new_nodes_y.push_back(ky+1);
          new_nodes_z.push_back(kz);
        }

        if (ky > 0 && node_array1[nx*ny*kz+nx*(ky-1)+kx] == -1 && cell_types[nx*ny*kz+nx*(ky-1)+kx] != NONCONDUCTIVE)
        {
          node_array1[nx*ny*kz+nx*(ky-1)+kx] = i;
          new_nodes_x.push_back(kx);
          new_nodes_y.push_back(ky-1);
          new_nodes_z.push_back(kz);
        }
      }
      old_nodes_x.resize(0);
      old_nodes_y.resize(0);
      old_nodes_z.resize(0);
      old_nodes_x.swap(new_nodes_x);
      old_nodes_y.swap(new_nodes_y);
      old_nodes_z.swap(new_nodes_z);
    }

    }

    {
    vector<int> old_nodes_x(0);
    vector<int> old_nodes_y(0);
    vector<int> old_nodes_z(0);
    vector<int> new_nodes_x(0);
    vector<int> new_nodes_y(0);
    vector<int> new_nodes_z(0);

    for (int kz = 0; kz < nz; kz++)
      for (int ky = 0; ky < ny; ky++)
        for (int kx = 0; kx < nx; kx++)
        {
          node_array2[nx*ny*kz+nx*ky+kx] = -1;
          if (contacts[nx*ny*kz+nx*ky+kx] == -1)
          {
            old_nodes_x.push_back(kx);
            old_nodes_y.push_back(ky);
            old_nodes_z.push_back(kz);
          }
        }

    for (int i = 0; i < nx*ny*nz; i++)
    {
      for (std::vector<int>::size_type l = 0; l < old_nodes_x.size(); l++)
      {
        int kx = old_nodes_x[l];
        int ky = old_nodes_y[l];
        int kz = old_nodes_z[l];

        if (kx < nx-1 && node_array2[nx*ny*kz+nx*ky+(kx+1)] == -1 && cell_types[nx*ny*kz+nx*ky+(kx+1)] != NONCONDUCTIVE)
        {
          node_array2[nx*ny*kz+nx*ky+(kx+1)] = i;
          new_nodes_x.push_back(kx+1);
          new_nodes_y.push_back(ky);
          new_nodes_z.push_back(kz);
        }

        if (kx > 0 && node_array2[nx*ny*kz+nx*ky+(kx-1)] == -1 && cell_types[nx*ny*kz+nx*ky+(kx-1)] != NONCONDUCTIVE)
        {
          node_array2[nx*ny*kz+nx*ky+(kx-1)] = i;
          new_nodes_x.push_back(kx-1);
          new_nodes_y.push_back(ky);
          new_nodes_z.push_back(kz);
        }

        if (ky < ny-1 && node_array2[nx*ny*kz+nx*(ky+1)+kx] == -1 && cell_types[nx*ny*kz+nx*(ky+1)+kx] != NONCONDUCTIVE)
        {
          node_array2[nx*ny*kz+nx*(ky+1)+kx] = i;
          new_nodes_x.push_back(kx);
          new_nodes_y.push_back(ky+1);
          new_nodes_z.push_back(kz);
        }

        if (ky > 0 && node_array2[nx*ny*kz+nx*(ky-1)+kx] == -1 && cell_types[nx*ny*kz+nx*(ky-1)+kx] != NONCONDUCTIVE)
        {
          node_array2[nx*ny*kz+nx*(ky-1)+kx] = i;
          new_nodes_x.push_back(kx);
          new_nodes_y.push_back(ky-1);
          new_nodes_z.push_back(kz);
        }
    }
      old_nodes_x.resize(0);
      old_nodes_y.resize(0);
      old_nodes_z.resize(0);
      old_nodes_x.swap(new_nodes_x);
      old_nodes_y.swap(new_nodes_y);
      old_nodes_z.swap(new_nodes_z);
    }

    }


    for (int i = 0; i < nx*ny*nz; i++)
    {
      nf2[i] = -1;
      node_array1[i] -= node_array2[i];
    }

    for (int i = 0; i < nx*ny*nz; i++) node_array2[i] = i;

    sort(node_array1, node_array2);

    for (int i = 0; i < nx*ny*nz; i++)
    {
      if (cell_types[node_array2[i]] == CONDUCTIVE)
      {
        nf[nnf] = node_array2[i];
        nf2[node_array2[i]] = nnf;
        nnf++;
      }
    }

    parallel.push_back(nnf);
    parallel.push_back(nnf);
    parallel.push_back(nnf);
    parallel.push_back(nnf);
    parallel.push_back(nnf);
  }
  else if (ordering == TILE_ORDERING)
  {
    int npx = 2;
    int npy = 2;

    for (int py = 0; py < npy; py++)
    {
      for (int px = 0; px < npx; px++)
      {
        for (int ky = 0; ky < ny/npy; ky++)
        {
          for (int kx = 0; kx < nx/npx; kx++)
          {
            if (kx == 0 || ky == 0)
            {
              for (int kz = 0; kz < nz; kz++)
              {
                int k = nx*ny*kz+nx*(ky+py*ny/npy)+(kx+px*nx/npx);
                if (cell_types[k] == CONDUCTIVE)
                {
                  nf[nnf] = k;
                  nf2[k] = nnf;
                  nnf++;
                }
                else
                {
                  nf2[k] = -1;
                }
              }
            }
          }
        }
      }
    }
    parallel.push_back(nnf);

    for (int py = 0; py < npy; py++)
    {
      for (int px = 0; px < npx; px++)
      {
        for (int ky = 1; ky < ny/npy; ky++)
        {
          for (int kx = 1; kx < nx/npx; kx++)
          {
            for (int kz = 0; kz < nz; kz++)
            {
              int k = nx*ny*kz;
              if (py % 2 == 0) k += nx*(ny/npy-ky+py*ny/npy);
              else k += nx*(ky+py*ny/npy);
              if (px % 2 == 0) k += nx/npx-kx+px*nx/npx;
              else k += kx+px*nx/npx;
              if (cell_types[k] == CONDUCTIVE)
              {
                nf[nnf] = k;
                nf2[k] = nnf;
                nnf++;
              }
              else
              {
                nf2[k] = -1;
              }
            }
          }
        }
        parallel.push_back(nnf);
      }
    }
  }
  else
  {
    cerr << "No valid ordering selected!\n";
  }


  solver.reset(new LinearSolver(nnf));
  Ln.resize(nnf+1);
  Un.resize(nnf+1);
  B.resize(nnf);

  this->nx = nx;
  this->ny = ny;
  this->nz = nz;
  this->dx = dx;
  this->dy = dy;
  this->dz = dz;

  this->AIsUpToDate = false;
  this->iteration = 0;

  Phi.resize(nx*ny*nz);
  m_x.resize(nx*ny*nz);
  m_y.resize(nx*ny*nz);
  m_z.resize(nx*ny*nz);
  sigma_0.resize(nx*ny*nz);
  amr.resize(nx*ny*nz);

  for (int k = 0; k < nx*ny*nz; k++)
  {
    Phi[k] = 0;
    sigma_0[k] = 0;
    amr[k] = 0;
    m_x[k] = 0;
    m_y[k] = 0;
    m_z[k] = 0;
  }

  for (int k = 0; k < 5; k++)
  {
    tickTime[k] = 0;
  }
}

CurrentCalculator::~CurrentCalculator()
{
  cerr << "CurrentCalculator: The calculation of " << iteration << " solutions took " << tickTime[0] << " s (average per solution: " << tickTime[0]/iteration << " s)\n";
  cerr << "CurrentCalculator: Time for preparation of the linear system: " << tickTime[1] << " s (average per solution: " << tickTime[1]/iteration << " s)\n";
  cerr << "CurrentCalculator: Time for solving the linear system: " << tickTime[2] << " s (average per solution: " << tickTime[2]/iteration << " s)\n";
  cerr << "CurrentCalculator: Time for calculation of the samples resistivity: " << tickTime[3] << " s (average per solution: " << tickTime[3]/iteration << " s)\n";
  cerr << "CurrentCalculator: Time for calculating the current: " << tickTime[4] << " s (average per solution: " << tickTime[4]/iteration << " s)\n";
}

