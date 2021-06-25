#ifndef CHARGE
#define CHARGE

#include "config.h"
#include "matrix/matty.h"
#include "mesh/Field.h"
#include "mesh/VectorField.h"



double topology_charge_continuous(int zi, const VectorField &M, const Field &Ms);
Field topology_charge_density_continuous(const VectorField &M, const Field &Ms);

inline double sigmaA(const Vector3d s1,const Vector3d s2,const Vector3d s3);

double topology_charge_berg_luescher_dual_lattice(int zi, const VectorField &M, const Field &Ms);
Field topology_charge_berg_luescher_density_dual_lattice(const VectorField &M, const Field &Ms);

double topology_charge_berg_luescher(int zi, const VectorField &M, const Field &Ms);
Field topology_charge_berg_luescher_density(const VectorField &M, const Field &Ms);

#endif
