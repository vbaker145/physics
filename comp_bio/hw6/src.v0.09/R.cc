
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
#include <iomanip>

#include "defs.h"
#include "macros.h"

#include "R.h"

#include "ensemble.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const R_meas &z)
{
  // FRIEND FUNCTION TO WRITE ELEMENTS OF DATA
//os << z.N  << " "
//<< "\n"
//;

  return os;
}


/*
 * ----------------------------------------------------------------------
 *                 END-TO-END AND RADIUS OF GYRATION CLASS
 */

R_meas::R_meas(void) { }

R_meas::~R_meas(void) 
{ 
  delete [] Rg ;
  delete [] Rg2;
  delete [] Eoe;
  delete [] Eoe_vec;
}

R_meas::R_meas(ensemble& sys) 
{ 
  Rg      = new double[sys.atoms_per_chain];
  Rg2     = new double[sys.atoms_per_chain];
  Eoe     = new double[sys.atoms_per_chain];
  Eoe_vec = new    vec[sys.atoms_per_chain];
}

void R_meas::calc(ensemble& sys, ofstream& out_rg, ofstream& out_rg2, ofstream& out_eoe, ofstream& out_eoe_vec, ofstream& out_eoe2)
{
  // FUNCTION THAT CALCULATES THE END-TO-END DISTANCE AND RADIUS OF GYRATION
  // PER CHAIN.
  int i,j,k;
  int atomi;
  double tot_mass;
  vec r;

  vec Rcm; 	// CENTER OF MASS

  // LOOP OVER CHAINS
  for (j=0; j<sys.num_chains; j++)
  {
	// LABEL FOR THE FIRST ATOM IN THIS CHAIN
	atomi = j * sys.atoms_per_chain;

	// INITS FOR THIS CHAIN
	Rcm.clear();
	tot_mass = 0.;
	Rg[j]    = 0.;
	Rg2[j]   = 0.;
	Eoe[j]   = 0.;

	// CALCULATE CENTER OF MASS OF THIS PARTICULAR CHAIN.
	for (i=0; i<sys.atoms_per_chain; i++)
	{
	  r   = sys.atoms[i+atomi].r;
	  // GET IMAGE
	  //r.x = IMAGE(r.x, sys.inv_Lx, sys.Lx);
	  //r.y = IMAGE(r.y, sys.inv_Ly, sys.Ly);
	  //r.z = IMAGE(r.z, sys.inv_Lz, sys.Lz);

	  Rcm      += sys.atoms[i+atomi].m * r;
	  tot_mass += sys.atoms[i+atomi].m;
	}
	Rcm /= tot_mass;

	// CALCULATE RADIUS OF GYRATION
	for (i=0; i<sys.atoms_per_chain; i++)
	{ 
	  r   = sys.atoms[i+atomi].r;
	  // GET IMAGE
	  //r.x = IMAGE(r.x, sys.inv_Lx, sys.Lx);
	  //r.y = IMAGE(r.y, sys.inv_Ly, sys.Ly);
	  //r.z = IMAGE(r.z, sys.inv_Lz, sys.Lz);

	  Rg2[j] += (r - Rcm).norm2();
	}
	Rg2[j] /= sys.atoms_per_chain;
	Rg[j]   = sqrt( Rg2[j] );

	// CALCULATE END-TO-END DISTANCE
	Eoe_vec[j] = sys.atoms[atomi].r - sys.atoms[sys.atoms_per_chain+atomi-1].r;
	Eoe[j]     = Eoe_vec[j].norm();
  }

  // WRITE TO FILE
  out_rg << setprecision(14) << sys.simulation_time << " ";
  out_rg2 << setprecision(14) << sys.simulation_time << " ";
  out_eoe << setprecision(14) << sys.simulation_time << " ";
  out_eoe_vec << setprecision(14) << sys.simulation_time << " ";
  out_eoe2 << setprecision(14) << sys.simulation_time << " ";
  for (j=0; j<sys.num_chains; j++)
  {
	out_rg  << Rg[j] << " ";
	out_rg2 << Rg2[j] << " ";
	out_eoe << Eoe[j] << " ";
	out_eoe_vec << Eoe_vec[j] << " ";
	out_eoe2 << Eoe[j]*Eoe[j] << " ";
  }
  out_rg  << "\n";
  out_rg2  << "\n";
  out_eoe << "\n";
  out_eoe_vec << "\n";
  out_eoe2 << "\n";
}

/*
 * ----------------------------------------------------------------------
 */
