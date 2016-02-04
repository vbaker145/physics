
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

#include <fstream>
#include <math.h>

#include "defs.h"
#include "macros.h"

#include "verletlist.h"

#include "ensemble.h"
#include "forces.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const verletlist &z)
{
  // FRIEND FUNCTION TO WRITE ELEMENTS OF DATA
//os << z.N  << " "
//<< "\n"
//;

  return os;
}


/*
 * ----------------------------------------------------------------------
 *                 INTEGRATOR CLASS
 */

verletlist::verletlist(void)
{
}


verletlist::~verletlist(void) 
{ 
  int i;
  for (i=0; i<N; i++)
	delete [] list[i];
  delete [] list;

  delete [] number_list;

  delete [] r_orig;
}


verletlist::verletlist(long Np, double rc, double rv)
{
  // CREATE A VERLET LIST WITH rv AS THE SIZE OF THE SKIN
  // AND rc THE SIZE OF THE POTENTIAL CUTOFF
  int i,j;

  // NUMBER OF PARTICLES
  N = Np;
  // SIZE OF POTENTIAL CUTOFF
  r_c = rc;
  // SIZE OF OUTER SKIN OF NEIGHBORHOOD
  r_skin = rv;

  // GET SOME SPACE FOR:
  // - VERLET LIST
  list = new long*[N]; 
  for (i=0; i<N; i++)
  	list[i] = new long[N];
  // - LIST OF NUMBER OF PARTICLES
  number_list = new long[N];
  // - LIST OF ORIGINAL POSITIONS
  r_orig = new vec[N];
}

int verletlist::check_verletlist(ensemble &sys)
{
  int i,j;
  int update = 0;
  double r2;
  double ave_r = 0.;
  vec dr;

  // DISTANCE TO CHECK	
  //double crit_dist2 = CUAD(r_c - r_skin);
  // SHOULD RATHER CHECK AT HALF OF THE DISTANCE INSTEAD
  double crit_dist2 = CUAD( (r_c - r_skin)/2. );

  // CHECK PARTICLE CURRENT POSITIONS WITH
  // ORIGINAL POSITIONS
  for (i=0; i<N; i++)
  {
	// FIND OUT DISPLACEMENT SINCE LAST UPDATE
	dr = sys.atoms[i].r - r_orig[i];

	// * DO NOT NEED IMAGE SINCE DISTANCES ARE
	// * FROM THE SAME PARTCLE.
	//dr.x = IMAGE(dr.x, sys.inv_Lx, sys.Lx);
	//dr.y = IMAGE(dr.y, sys.inv_Ly, sys.Ly);
	//dr.z = IMAGE(dr.z, sys.inv_Lz, sys.Lz);

	r2 = dr.norm2();

	// NOT USED IN THE CALCULATION. JUST TO
	// CHECK ON AVERAGE DISPLACEMENT.
	ave_r += r2;

	if (r2>=crit_dist2)
	{
	  update = 1;
	  break;
	}
  }
  //PRINT2("Ave r2",ave_r/N);

  return update;
}

void verletlist::update_verletlist(ensemble &sys)
{
  // BRING UP TO DATE THE VERLET LIST
  int i,j;
  double r_skin2 = CUAD(r_skin);
  double r2;
  vec dr;

  // UPDATE STORED POSITIONS
  for (i=0; i<N; i++)
  {
	r_orig[i] = sys.atoms[i].r;
  }

  // INITIALIZE LIST OF # OF NEIGHBORING PARTICLES
  for (i=0; i<N; i++)
  	number_list[i] = 0;

  // ADD PARTICLES TO LISTS DEPENDING
  // ON DISTANCE
  for (i=0; i<N-1; i++)
	for (j=i+1; j<N; j++)
	{
	  dr = sys.atoms[i].r - sys.atoms[j].r;
	  // USE IMAGE DISTANCES
	  dr.x = IMAGE(dr.x, sys.inv_Lx, sys.Lx);
	  dr.y = IMAGE(dr.y, sys.inv_Ly, sys.Ly);
	  dr.z = IMAGE(dr.z, sys.inv_Lz, sys.Lz);

	  r2 = dr.norm2();

	  // CHECK WHETHER IS A NEIGHBOR
	  if (r2 <= r_skin2)
	  {
		// ADD IDs OF NEIGHBORING PARTICLES
		// TO EACH OTHERS LISTS
		list[i][number_list[i]] = j;
		list[j][number_list[j]] = i;

		// UPDATE THE NUMBER OF NEIGHBOR LISTS
		number_list[i]++;
		number_list[j]++;
	  }
	}
}

int verletlist::is_in_list(int ii, int jj)
{
  // RETURNS 1 IF PARTICLE jj IS IN THE LIST
  // OF PARTICLE ii
  int i,j;
  int is_in = 0;
  for (j=0; j<number_list[ii]; j++)
  {
	if ( jj == list[ii][j] ) is_in = 1;
  }
  return is_in;
}

/*
 * ----------------------------------------------------------------------
 */
