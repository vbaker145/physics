
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

#include "defs.h"
#include "macros.h"

#include "atom.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const atom &z)
{
  // FRIEND FUNCTION TO WRITE ELEMENTS OF DATA
  os << z.m  << " "
     << z.d  << " "
     << z.r  << " "
     << z.r_old  << " "
     << z.v  << " "
     << z.f  << " "
	 << "\n"
	;

  return os;
}


/*
 * ----------------------------------------------------------------------
 *                 DATA CLASS
 */

atom::atom(void)
{
  // DEFAULT INITIALIZATIONS.

  m = -1.;
  d = -1.;
}


atom::~atom(void) { }

void atom::output_position(void)
{
	cout << r << "\n";
}


void atom::output_velocity(void)
{
	cout << v << "\n";
}


/*
 * ----------------------------------------------------------------------
 */
