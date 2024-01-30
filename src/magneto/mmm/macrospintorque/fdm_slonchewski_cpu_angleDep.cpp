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

#include "config.h"
#include "fdm_slonchewski_cpu_angleDep.h"
#include <iostream>
#include <string>
#include <math.h>

#include "mmm/constants.h"


#include "../external_code/cpp_string2function_parser/exprtk/exprtk.hpp"

void fdm_slonchewski_cpu_angleDep(
	int dim_x, int dim_y, int dim_z,
	double delta_x, double delta_y, double delta_z,
	double r_x, double r_y, double r_z,
	double time,
	const Matrix &Ms,
	const Matrix &alpha,
	const VectorMatrix &j,
	const Matrix &a_DL,
	const Matrix &a_FL,
	const std::string &fn_DL,
	const std::string &fn_FL,
	const VectorMatrix &M,
	VectorMatrix &dM)
{
	// Calculate: 
	//   c1*(M x (M x p)) + c2*(M x p)
	//
	//   c1(theta): damping factor
	//   c2(theta): precession factor
	//   Ms*cos(theta) = M*p

	Matrix::ro_accessor Ms_acc(Ms), alpha_acc(alpha), a_DL_acc(a_DL), a_FL_acc(a_FL);
	VectorMatrix::const_accessor j_acc(j);
	VectorMatrix::const_accessor M_acc(M);
	VectorMatrix::accessor dM_acc(dM);

	typedef exprtk::symbol_table<double> symbol_table_t;
	typedef exprtk::expression<double>     expression_t;
	typedef exprtk::parser<double>             parser_t;

	double var_x;

	symbol_table_t symbol_table;
	symbol_table.add_variable("x", var_x);
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
        //double vals[];
	const Vector3d ref = Vector3d(r_x, r_y, r_z);
	for (int n=0; n<N; ++n) {
		const double alpha_n = alpha_acc.at(n);
		const double Ms_n    = Ms_acc.at(n);

		const Vector3d j_n = j_acc.get(n);
		const Vector3d M_n = M_acc.get(n);
		
		if (j_n  == Vector3d(0.0, 0.0, 0.0)) continue;
		if (Ms_n == 0.0) continue;
		const Vector3d e_z = Vector3d(0.0, 0.0, 1.0);

		// Calculate precession and damping terms
		const Vector3d p_n   = cross(e_z, j_n); // spin polarization e_z x J with J = technical current direction
		const Vector3d Mxp   = cross(M_n, p_n); // precession: u=mxp
		const Vector3d MxMxp = cross(M_n, Mxp); // damping:  t=mxu=mx(mxp)

		// add both terms to dm/dt in LLGE
		const double gamma_pr = GYROMAGNETIC_RATIO / (1.0 + alpha_n*alpha_n);


		// evaluate angular dependence
		const double dotProduct = ref.x*M_n.x + ref.y*M_n.y + ref.z*M_n.z;
		const double norm       = sqrt((ref.x*ref.x + ref.y*ref.y + ref.z*ref.z)*(M_n.x*M_n.x + M_n.y*M_n.y + M_n.z*M_n.z));

		var_x = acos(dotProduct/norm);

		//const double f_DL = expression_DL.value();
		//const double f_FL = expression_FL.value();
		//std::cout << var_x << " " << f_DL << " " << f_FL << std::endl;


		const double a_DL_n  = a_DL_acc.at(n)*expression_DL.value();
		const double a_FL_n  = a_FL_acc.at(n)*expression_FL.value();
		
		Vector3d dM_n = dM_acc.get(n);
		dM_n.x += gamma_pr * ((alpha_n * a_FL_n + a_DL_n) * MxMxp.x/Ms_n + (a_FL_n - alpha_n * a_DL_n) * Mxp.x);
		dM_n.y += gamma_pr * ((alpha_n * a_FL_n + a_DL_n) * MxMxp.y/Ms_n + (a_FL_n - alpha_n * a_DL_n) * Mxp.y);
		dM_n.z += gamma_pr * ((alpha_n * a_FL_n + a_DL_n) * MxMxp.z/Ms_n + (a_FL_n - alpha_n * a_DL_n) * Mxp.z);
		dM_acc.set(n, dM_n);
	}
}

