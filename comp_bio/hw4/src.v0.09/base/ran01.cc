
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

#include <math.h>
#include <iostream>

#include "macros.h"

#include "ran01.h"

using namespace std;

double rand01(void)
{
  // THIS IS AN ALTERNATIVE FUNCTION FOR RANDOM NUMBERS.
  // CHECK THE SIZE OF 'RAND_MAX'. LAST TIME THAT CHECKED
  // WAS RAND_MAX = 2147483647

  double a;

  a = (double) rand() / ( (double) RAND_MAX );

  return a;
}

void ran_gaussian_01(double* nums)
{
  // METHOD 1.
  // FUNCTION THAT RETURNS RANDOM NUMBERS DRAWN FROM A
  // GASSIAN (NORMAL) DISTRIBUTION. 
  // IT USES THE rand01 FUNCTION

  int i;

  double xi1 = rand01();
  double xi2 = rand01();

  nums[0] = sqrt( -2.0 * log(xi1) ) * cos( 2.0*PI * xi2 );
  nums[1] = sqrt( -2.0 * log(xi1) ) * sin( 2.0*PI * xi2 );
}

//double ran_gaussian_02(void)
double ran_gaussian_02(int a)
{
  // METHOD 2.
  // FUNCTION THAT RETURNS RANDOM NUMBERS DRAWN FROM A
  // GASSIAN (NORMAL) DISTRIBUTION. 
  // IT USES THE rand01 FUNCTION

  int i;
  double x = 0.;

  for (i=0; i<a; i++)
  {
	x += rand01();
  }
  //x -= 6.;
  x -= a/2.;

  return x;
}

double ran_gaussian_02_mu_sigma(double mu, double sigma)
{
  // METHOD 2.
  // FUNCTION THAT RETURNS RANDOM NUMBERS DRAWN FROM A
  // GASSIAN (NORMAL) DISTRIBUTION WITH 'mu' MEAN AND 'sigma' VARIANCE.
  // IT USES THE rand01 FUNCTION

  int i;
  double x = 0.;

  for (i=0; i<12; i++)
  {
	x += rand01();
  }
  x -= 6.;

  return (mu + sigma*x);
}

// =============================
