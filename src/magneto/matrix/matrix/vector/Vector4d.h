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

#ifndef MATTY_VECTOR4D_H
#define MATTY_VECTOR4D_H

#include "config.h"
#include <ostream>
#include <cmath>
#include <cassert>
#include <immintrin.h>
namespace matty {

/** @file */ 

/**
 * Vector with three elements.
 */
struct Vector4d
{
	/**
	 * Constructor. Initializes all elements to zero.
	 */
	inline Vector4d();
	
	/**
	 * Constructor with element initialization.
	 */
	inline Vector4d(double x, double y, double z, double w);

	/**
	 * Constructor with equal x,y,z element initialization.
	 */
	inline explicit Vector4d(double xyzw);

	/**
	 * Default destructor.
	 */
	inline ~Vector4d();

	/**
	 * Default copy constructor.
	 */
	inline Vector4d(const Vector4d &other);

	/**
	 * Get component, c = 0, 1, or 2.
	 */
	inline double &operator[](int c)
	{
		switch (c) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
            case 3: return w;        
		}

		assert(0);
		static double foo = 0.0; return foo; // supress compiler warning...
	}

	inline const double &operator[](int c) const
	{
		switch (c) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
            case 3: return w;
		}

		assert(0);
		static double foo = 0.0; return foo; // supress compiler warning...
	}

	/**
	 * Default assignment operator.
	 */
	inline Vector4d &operator=(const Vector4d &other);

	/**
	 * Adds the vector 'other' to this vector.
	 */
	inline Vector4d &operator+=(const Vector4d &other);

	/**
	 * Substracts the vector 'other' from this vector.
	 */
	inline Vector4d &operator-=(const Vector4d &other);

	/**
	 * Copy contents of the vector 'other' to this vector.
	 */
	inline void assign(Vector4d &other);

	/**
	 * Explicitly assign the elements of this vector.
	 */
	inline void assign(double x, double y, double z, double w);

	/**
	 * Return the length of this vector.
	 */
	inline double abs() const;

	/**
	 * Return the squared length of this vector, faster than abs().
	 */
	inline double abs_squared() const;

	/**
	 * Normalize this vector to length 'norm'.
	 */
	inline void normalize(double norm = 1.0);

	double x, y, z, w; /** the elements of the vector */
};

/**
 * Returns the dot/scalar product of the vectors 'lhs' and 'rhs'.
 */
inline double dot(const Vector4d &lhs, const Vector4d &rhs)
{
	return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
}

/**
 * Returns a vector pointing into the same direction with length 'len'.
 */
inline Vector4d normalize(const Vector4d &lhs, double len = 1.0)
{
	Vector4d res(lhs);
	res.normalize(len);
	return res;
}

/**
 * Returns the cross product of the vectors lhs and rhs.
 */
inline Vector4d cross(const Vector4d &lhs, const Vector4d &rhs)
{
	return Vector4d(
		lhs.y*rhs.z - lhs.z*rhs.y,
		lhs.z*rhs.x - lhs.x*rhs.z,
		lhs.x*rhs.y - lhs.y*rhs.x,
        0
	);
}





/**
 * Returns the result of scaling the vector lhs with the scalar rhs.
 */
inline Vector4d operator*(const Vector4d &lhs, const double &rhs)
{
	return Vector4d(
		lhs.x * rhs,
		lhs.y * rhs,
		lhs.z * rhs,
        0
	);
}

/**
 * Returns the result of scaling the vector rhs with the scalar lhs.
 */
inline Vector4d operator*(const double lhs, const Vector4d &rhs)
{
	return Vector4d(
		rhs.x * lhs,
		rhs.y * lhs,
		rhs.z * lhs,
        0
	);
}

/**
 * Returns the result of dividing the elements of the vector lhs by the scalar rhs.
 */
inline Vector4d operator/(const Vector4d &lhs, const double &rhs)
{
	return Vector4d(
		lhs.x / rhs,
		lhs.y / rhs,
		lhs.z / rhs,
        0
	);
}

/**
 * Returns the result of summing the vectors lhs and rhs.
 */
inline Vector4d operator+(const Vector4d &lhs, const Vector4d &rhs)
{
	return Vector4d(
		lhs.x + rhs.x,
		lhs.y + rhs.y,
		lhs.z + rhs.z,
        0
	);
}

/**
 * Returns the result of subtracting the vectors lhs from rhs.
 */
inline Vector4d operator-(const Vector4d &lhs, const Vector4d &rhs)
{
	return Vector4d(
		lhs.x - rhs.x,
		lhs.y - rhs.y,
		lhs.z - rhs.z,
        0
	);
}

/**
 * Compares the vector lhs and rhs element-wise.
 */
inline bool operator==(const Vector4d &lhs, const Vector4d &rhs)
{
	return (lhs.x == rhs.x) && (lhs.y == rhs.y) && (lhs.z == rhs.z);
}

/**
 * Compares the vector lhs and rhs element-wise.
 */
inline bool operator!=(const Vector4d &lhs, const Vector4d &rhs)
{
	return !(lhs == rhs);
}

/**
 * Output to stream.
 */
std::ostream &operator<<(std::ostream &out, const Vector4d &vec);

////////////////////////////////////////////////////////////////////////////////

Vector4d::Vector4d()
	: x(0), y(0), z(0), w(0)
{
}

Vector4d::Vector4d(double x, double y, double z, double w)
	: x(x), y(y), z(z), w(w)
{
}

Vector4d::Vector4d(double xyzw)
	: x(xyzw), y(xyzw), z(xyzw), w(xyzw)
{
}


Vector4d::~Vector4d()
{
}

Vector4d::Vector4d(const Vector4d &other)
{
	assign(other.x, other.y, other.z, other.w);
}

Vector4d &Vector4d::operator=(const Vector4d &other)
{
	assign(other.x, other.y, other.z, other.w);
	return *this;
}

//#Vector4d &Vector4d::operator=(const float number)
//#//{
//	assign(number, number, number, 0);
//	return *this;
//}
Vector4d &Vector4d::operator+=(const Vector4d &other)
{
	x += other.x; y += other.y; z += other.z;
	return *this;
}

Vector4d &Vector4d::operator-=(const Vector4d &other)
{
	x -= other.x; y -= other.y; z -= other.z;
	return *this;
}

double Vector4d::abs() const 
{ 
	return std::sqrt(abs_squared()); 
}

double Vector4d::abs_squared() const 
{ 
	return dot(*this, *this);
}

void Vector4d::normalize(double norm)
{
	const double len = abs();
	if (len == 0.0) {
		return;
	}

	const double factor = norm / len;
	x *= factor;
	y *= factor;
	z *= factor;
}

void Vector4d::assign(Vector4d &other)
{
	assign(other.x, other.y, other.z, other.w);
}

void Vector4d::assign(double x, double y, double z,double w )
{
	this->x = x; this->y = y; this->z = z, this->w = w;
}
}

#endif

