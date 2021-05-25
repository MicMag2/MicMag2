#ifndef MULTISCALE_H
#define MULTISCALE_H

#include "config.h"
#include "matrix/matty.h"

#include <math.h>
#include <iostream>


double Hhb(double x, double y, double z, double Lx, double Ly, double Lz, double Mx, double My, double Mz);
double Hmm(double x, double y, double z, double Lx, double Ly, double Lz, double mx, double my, double mz);

#endif
