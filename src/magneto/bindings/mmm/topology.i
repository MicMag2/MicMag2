%{
#include "mmm/topology/charge.h"
%}

double topology_charge_continuous(int zi, const VectorField &M, const Field &Ms);
Field topology_charge_density_continuous(const VectorField &M, const Field &Ms);

double topology_charge_berg_luescher_dual_lattice(int zi, const VectorField &M, const Field &Ms);
Field topology_charge_berg_luescher_density_dual_lattice(const VectorField &M, const Field &Ms);

double topology_charge_berg_luescher(int zi, const VectorField &M, const Field &Ms);
Field topology_charge_berg_luescher_density(const VectorField &M, const Field &Ms);