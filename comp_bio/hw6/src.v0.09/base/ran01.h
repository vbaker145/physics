
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

#include  <stdlib.h>

#ifndef RAN01_H
#define RAN01_H

// THIS IS AN ALTERNATIVE FUNCTION FOR RANDOM NUMBERS.
// CHECK THE SIZE OF 'RAND_MAX'. LAST TIME THAT CHECKED
// WAS RAND_MAX = 2147483647
// TO IMPLEMENT, JUST REPLACE ran3 FOR ran01
double rand01(void);

// FUNCTION THAT RETURNS RANDOM NUMBERS DRAWN FROM A
// GAUSSIAN DISTRIBUTION WITH 0 MEAN AND UNIT VARIANCE
// METHOD 1
void ran_gaussian_01(double* nums);

// FUNCTION THAT RETURNS RANDOM NUMBERS DRAWN FROM A
// GAUSSIAN DISTRIBUTION WITH 0 MEAN AND UNIT VARIANCE
// METHOD 2
double ran_gaussian_02(int);

// FUNCTION THAT RETURNS RANDOM NUMBERS DRAWN FROM A
// GAUSSIAN DISTRIBUTION WITH 'mu' MEAN AND 'sigma' VARIANCE
// METHOD 2
double ran_gaussian_02_mu_sigma(double mu, double sigma);

#endif // RAN01_H
