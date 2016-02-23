
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

#ifndef _VEC_H_
#define _VEC_H_

#include <iostream>
using namespace std;

class vec
{
  friend ostream &operator<<(ostream &os, const vec &v);
  friend istream &operator>>(istream &os, vec &v);
  friend vec operator-(const vec &u, const vec &v);
  friend vec operator-(const vec &u);
  friend vec operator+(const vec &u, const vec &v);
  friend vec operator/(const vec &u, const double &v);
  friend vec operator*(const vec &u, const int &v);
  friend vec operator*(const vec &u, const double &v);
  friend vec operator*(const double &v, const vec &u);
  friend vec operator/=(vec &v, const double &u);
  friend vec operator*=(vec &v, const double &u);
  friend vec operator+=(vec &v, const vec &u);
  friend vec operator-=(vec &v, const vec &u);
  friend int operator==(const vec &v, const vec &u);
  friend int operator!=(const vec &v, const vec &u);
  friend int operator<(vec &v, vec &u);
  friend int operator>(vec &v, vec &u);

  public: 
  	vec(void);
  	vec(double ox, double oy, double oz);
  	~vec(void);

	double norm(void);
	double XYnorm(void);
	double Znorm(void);
	double norm2(void);
	double XYnorm2(void);
	double Znorm2(void);
	vec &abs(void);

	// PUT ALL ELEMENTS TO ZERO
	void reset(void);
	void clear(void);

  // ===============================================
  // DEFAULT: DECLARE PRIVATE THE COPY AND EQUALITY OPERATORS
  private:
  //vec(const vec &z);
  //vec &operator=(const vec &v);
  // IMPLEMENTATIONS OF THE COPY AND EQUALITY OPERATORS
  public:
  	vec(const vec &v);
    vec &operator=(const vec &v);
  // ===============================================

  public:
	double x, y, z;

};

double ScalarProd(const vec &u, const vec &v);
void VectorProd(vec &w, const vec &u, const vec &v);
double Angle(const vec &w, const vec &u, const vec &v);
vec imageDistance(const vec &u, const vec &v, double boxinv, double box);

#endif // _VEC_H_
