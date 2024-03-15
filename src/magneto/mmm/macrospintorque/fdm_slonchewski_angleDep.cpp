/*
 * Copyright 2012, 2013 by the Micromagnum authors.
 *
 * This file is part of MicroMagnum.
 * 
 * MicroMagnum is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * MicroMagnum is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with MicroMagnum.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include "config.h"
#include <string>
#include "fdm_slonchewski_angleDep.h"
#include "fdm_slonchewski_cpu_angleDep.h"
#ifdef HAVE_CUDA
#include "fdm_slonchewski_cuda_angleDep.h"
#include <cuda_runtime.h>
#endif

#include "Magneto.h"
#include "Benchmark.h"

#include "../external_code/exprtk.hpp"

void fdm_slonchewski_angleDep(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	double time,
	const Matrix &Ms,
	const Matrix &alpha,
	const VectorMatrix &j,
	const Matrix &xi,
	const Matrix &alpha_hall,
	const std::string &fn_DL,
	const std::string &fn_FL,
	const VectorMatrix &M,
	VectorMatrix &dM)
{
	const bool use_cuda = isCudaEnabled();

	if (use_cuda) {
#ifdef HAVE_CUDA


		//calculate function values in advance (not possible on gpu)
		VectorMatrix::const_accessor M_acc(M);
		Matrix parse_a_DL = xi;   //quick 'n' dirty initialization
		Matrix parse_a_FL = xi;
		Matrix::rw_accessor parse_a_DL_acc(parse_a_DL), parse_a_FL_acc(parse_a_FL);

		typedef exprtk::symbol_table<double> symbol_table_t;
		typedef exprtk::expression<double>     expression_t;
		typedef exprtk::parser<double>             parser_t;

		double var_theta;
		double var_phi;

		symbol_table_t symbol_table;
		symbol_table.add_variable("x", var_theta);
		symbol_table.add_variable("y", var_phi);
		symbol_table.add_variable("t", time);
		symbol_table.add_constants();

		expression_t expression_DL;
		expression_t expression_FL;
		expression_DL.register_symbol_table(symbol_table);
		expression_FL.register_symbol_table(symbol_table);

		parser_t parser_DL;
		parser_t parser_FL;
		parser_DL.compile(fn_DL,expression_DL);
		parser_FL.compile(fn_FL,expression_FL);

		const int N = dim_x * dim_y * dim_z;
		for (int n=0; n<N; ++n) {
			const Vector3d M_n = M_acc.get(n);	
			if (M_n  == Vector3d(0.0, 0.0, 0.0)) continue;

			// evaluate angular dependence
			const double norm       = sqrt(M_n.x*M_n.x + M_n.y*M_n.y + M_n.z*M_n.z);
			
			const double sigy		= M_n.y/sqrt(M_n.y*M_n.y);
			const double xynorm		= sqrt(M_n.x*M_n.x + M_n.y*M_n.y);

			var_theta = acos(M_n.z/norm);
			var_phi   = sigy * acos(M_n.x/xynorm);

			//const double f_DL = expression_DL.value();
			//const double f_FL = expression_FL.value();
			//std::cout << var_x << " " << f_DL << " " << f_FL << std::endl;

			parse_a_DL_acc.at(n)  = expression_DL.value();	//overwrite parser object with evaluated function call
			parse_a_FL_acc.at(n)  = expression_FL.value();
		}

		//temporarily unlock accessors for calling them in cuda code
		parse_a_DL.writeUnlock(0);
		parse_a_FL.writeUnlock(0);
		
		CUTIC("fdm_slonchewski_ANGLEDEP");
		fdm_slonchewski_cuda_angleDep(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, Ms, j, alpha, xi, alpha_hall, parse_a_DL, parse_a_FL, M, dM, isCuda64Enabled());
		CUTOC("fdm_slonchewski_ANGLEDEP");

		//lock accessors again, otherwise destructor crashes
		parse_a_DL.writeLock(0);
		parse_a_FL.writeLock(0);
#else
		assert(0);
#endif
	} else {
		TIC("fdm_slonchewski_ANGLEDEP");
		fdm_slonchewski_cpu_angleDep(dim_x, dim_y, dim_z, delta_x, delta_y, delta_z, time, Ms, alpha, j, xi, alpha_hall, fn_DL, fn_FL, M, dM);
		TOC("fdm_slonchewski_ANGLEDEP");
	}
}
