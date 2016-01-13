
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

// THIS IS AN ALTERNATIVE FUNCTION FOR RANDOM NUMBERS.
// CHECK THE SIZE OF 'RAND_MAX'. LAST TIME THAT CHECKED
// WAS RAND_MAX = 2147483647
// TO IMPLEMENT, JUST REPLACE ran3 FOR ran01

#include "ran01.h"

double rand01(void)
{
  double a;

  a = (double) rand() / ( (double) RAND_MAX );

  return a;
}

// =============================
