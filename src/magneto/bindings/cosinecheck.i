//%module cosinecheck

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

%{
#include "../cosinecheck.h"
%}


class cosinecheck{
	public:
		int idx_x0;
		int idx_y0;
		int idx_z0;
		Vector3d r_0;
		Vector3d r_1;
		int idx_x1;
		int idx_y1;
		int idx_z1;
		float cos_r;
		float cos_u;
		float cos_f;

	cosinecheck(
	int dim_x, int dim_y, 
	int dim_z, const VectorMatrix M,
	const Matrix Ms, float treshold);
	};


