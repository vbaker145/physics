
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
     << z.v_new  << " "
     << z.f  << " "
     << z.f_old  << " "
     << z.fixed  << " "
     << z.theta  << " "
     << z.phi  << " "
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

  m = 1.;
  d = -1.;
  fixed = 0; // MOBILE ATOM

  r = vec(0.,0.,0.);
  r_old = vec(0.,0.,0.);
  v = vec(0.,0.,0.);
  v_new = vec(0.,0.,0.);
  f = vec(0.,0.,0.);
  f_old = vec(0.,0.,0.);

  theta = 0.;
  phi   = 0.;
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

short atom::is_fixed(void)
{
  return fixed;
}

short atom::is_not_fixed(void)
{
  return !fixed;
}

void atom::set_fixed(void)
{
  fixed = 1;
}


/*
 * ----------------------------------------------------------------------
 */
