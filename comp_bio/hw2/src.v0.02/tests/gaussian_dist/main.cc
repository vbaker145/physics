
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


#include <iostream>
#include <stdlib.h>

#include "ran01.h"

using namespace std;

int main(int argc, char* argv[])
{
  // PROGRAM TO TEST THE ran_gaussian_01() SUBROUTINE BY
  // GENERATING AN ARRAY OF RANDOM NUMBERS DRAWN FROM
  // A GAUSSIAN DISTRIBUTION

  // =============================================
  // SET SEED FOR RANDOM GENERATOR
  long seed;
  long time_sub = 1234234670;
  // USE TIME AS SEED
  seed = - long( time(NULL) - time_sub );
  // SEED THE GENERATOR
  srand(seed);

  // =============================================
  // GENERATE THE ARRAY
  double *nums = new double[2];
  for (int i=0; i<100000; i++)
  {
	// PLAIN UNIFORM RANDOM DIST
	cout << rand01() << "\n";

	// METH0D 1
	//ran_gaussian_01(nums);
	//cout << nums[0] << "\n"
	//     << nums[1] << "\n";

	// METHOD 2
	//cout << ran_gaussian_02() << "\n";
  }

  delete [] nums;

  return 1;
}
