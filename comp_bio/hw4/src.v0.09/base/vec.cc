
/*
 * Code to accompany the PHYS 462/562 Computational Biophysics course.
 * 
 * Copyright (c) 2009-2010 Luis R Cruz Cruz, Philadelphia, PA
 *
 * Please address questions to:
 *  Luis R Cruz Cruz
 *  Physics Department
 *  Drexel University
 *  3141 Chestnut St.
 *  Philadelphia, PA 19104
 */

#include "vec.h"
#include "macros.h"

#include <math.h>

/* ============================================= */

double ScalarProd(const vec &u, const vec &v)
{
	return ( u.x*v.x + u.y*v.y + u.z*v.z );
}

void VectorProd(vec &w, const vec &u, const vec &v)
{
	w.x = u.y * v.z - u.z * v.y;
	w.y = u.z * v.x - u.x * v.z;
	w.z = u.x * v.y - u.y * v.x;
}

double Angle(const vec &w, const vec &u, const vec &v)
{
	vec a = w - u;
	vec b = v - u;
	a /= a.norm();
	b /= b.norm();

	vec s;
	VectorProd(s,a,b);

	double t = ScalarProd(a,b);

	//double ang1 = asin( s.norm() );
	double ang2 = acos( t        );
	/*
	cout << "asin: " << s.norm() << " " << DEGREE(ang1)
	  << " acos: " << t << " " << DEGREE(ang2) << "\n";
	*/

	// ALL ANGLES WILL BE < 180, RETURN THEN THE ACOS RESULT.
	return DEGREE(ang2);
}

/* ============================================= */

ostream &operator<<(ostream &os, const vec &v)
{
	os << v.x << " " 
		<< v.y << " " 
		<< v.z << " " 
		;
	return os;
}

istream &operator>>(istream &os, vec &v)
{
	os >> v.x 
		>> v.y 
		>> v.z 
		;

	return os;
}                           

vec operator+(const vec &u, const vec &v)
{
  vec w;
  w.x = u.x + v.x;
  w.y = u.y + v.y;
  w.z = u.z + v.z;
  return w;
}

vec operator-(const vec &u, const vec &v)
{
  vec w;
  w.x = u.x - v.x;
  w.y = u.y - v.y;
  w.z = u.z - v.z;
  return w;
}

vec operator-(const vec &u)
{
  vec w;
  w.x = - u.x;
  w.y = - u.y;
  w.z = - u.z;
  return w;
}

vec operator/(const vec &u, const double &v)
{
  vec w;
  w.x = u.x / v;
  w.y = u.y / v;
  w.z = u.z / v;
  return w;
}

vec operator*(const vec &u, const int &v)
{
  vec w;
  w.x = u.x * double(v);
  w.y = u.y * double(v);
  w.z = u.z * double(v);
  return w;
}

vec operator*(const vec &u, const double &v)
{
  vec w;
  w.x = u.x * v;
  w.y = u.y * v;
  w.z = u.z * v;
  return w;
}

vec operator*(const double &v, const vec &u)
{
  vec w;
  w.x = u.x * v;
  w.y = u.y * v;
  w.z = u.z * v;
  return w;
}

vec operator/=(vec &v, const double &u)
{
  v.x /= u;
  v.y /= u;
  v.z /= u;
  return v;
}

vec operator*=(vec &v, const double &u)
{
  v.x *= u;
  v.y *= u;
  v.z *= u;
  return v;
}

vec operator+=(vec &v, const vec &u)
{
  v.x += u.x;
  v.y += u.y;
  v.z += u.z;
  return v;
}

vec operator-=(vec &v, const vec &u)
{
  v.x -= u.x;
  v.y -= u.y;
  v.z -= u.z;
  return v;
}

int operator==(const vec &v, const vec &u)
{
  if ( (v.x == u.x) && (v.y == u.y) && (v.z == u.z) )
  	return 1;
  else
	return 0;
}

int operator!=(const vec &v, const vec &u)
{
  if ( (v.x != u.x) || (v.y != u.y) || (v.z != u.z) )
  	return 1;
  else
	return 0;
}

int operator<(vec &v, vec &u)
{
  if ( (v.x<u.x) && (v.y<u.y) && (v.z<u.z) )
	return 1;
  else
	return 0;
}

int operator>(vec &v, vec &u)
{
  if ( (v.x>u.x) && (v.y>u.y) && (v.z>u.z) )
	return 1;
  else
	return 0;
}

/* ================================================ */

vec::vec(void)
{
  x = 0.0;
  y = 0.0;
  z = 0.0;
}

vec::vec(double ox, double oy, double oz)
{
  x = ox;
  y = oy;
  z = oz;
}

vec::vec(const vec &v)
{
  x = v.x;
  y = v.y;
  z = v.z;
}

vec &vec::operator=(const vec &v)
{
  x = v.x;
  y = v.y;
  z = v.z;

  return *this;
}

vec::~vec(void)
{
}

double vec::norm(void)
{
  return ( sqrt(x*x + y*y + z*z) );
}

double vec::XYnorm(void)
{
  return ( sqrt(x*x + y*y) );
}

double vec::Znorm(void)
{
  return ( ABS(z) );
}

double vec::norm2(void)
{
  return ( x*x + y*y + z*z );
}

double vec::XYnorm2(void)
{
  return ( (x*x + y*y) );
}

double vec::Znorm2(void)
{
  return ( (z*z) );
}


vec &vec::abs(void)
{
  x = ((x)<0) ? (-(x)) : (x);
  y = ((y)<0) ? (-(y)) : (y);
  z = ((z)<0) ? (-(z)) : (z);

  return *this;
}

void vec::reset(void)
{
	x = 0.;
	y = 0.;
	z = 0.;
}

void vec::clear(void)
{
	x = 0.;
	y = 0.;
	z = 0.;
}

vec imageDistance(const vec &u, const vec &v, double boxinv, double box)
{
	// RETURNS A VECTOR THAT IS THE IMAGE DISTANCE, AS IN SIMULATIONS
	// WITH PERIODIC BOUNDARY CONDITIONS, THAT IS THE SHORTEST DISTANCE.
	double xi = IMAGE(u.x-v.x, boxinv, box);
	double yi = IMAGE(u.y-v.y, boxinv, box);
	double zi = IMAGE(u.z-v.z, boxinv, box);
	/*
	double xi = 0.;
	double yi = 0.;
	double zi = 0.;
	*/

	vec w(xi, yi, zi);
	return w;
}

/* ================================================ */
