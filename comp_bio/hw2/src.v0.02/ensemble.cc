
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

#include "ensemble.h"

#include "ran01.h"

/*
 * ----------------------------------------------------------------------
 *                 FRIEND MEMBERS
 */

using namespace std;

ostream &operator<<(ostream &os, const ensemble &z)
{
  // FRIEND FUNCTION TO WRITE ELEMENTS OF DATA
  os << z.N  << " "
     << z.delta_t  << " "
     << z.density  << " "
     << z.Lx  << " "
     << z.Ly  << " "
#if NDIM == 3
     << z.Lz  << " "
#endif
     << z.inv_Lx  << " "
     << z.inv_Ly  << " "
#if NDIM == 3
     << z.inv_Lz  << " "
#endif
     << z.total_time  << " "
     << z.total_num_steps  << " "
     << z.temperature  << " "
     << z.pot_energy  << " "
     << z.kin_energy  << " "
     << z.energy  << " "
     << z.virSum  << " "
     << z.pressure  << " "
	<< "\n"
	;
/*
  // WRITE ALL THE ATOMS
  for (int i=0; i<z.N; i++)
  {
	os <<        i  << " "
	   << z.atoms[i]  << " "
	  //<< "\n"
	  ;
  }
*/
  return os;
}


/*
 * ----------------------------------------------------------------------
 *                 DATA CLASS
 */

ensemble::ensemble(void)
{
  // DEFAULT INITIALIZATIONS.
  // THESE SHOULD BE OVERRIDEN BY VALUES
  // READ FROM FILE.
  N = -1;
  delta_t = -1.;
  density = -1.;
  Lx = -1.;
  Ly = -1.;
  Lz = -1.;
  inv_Lx = -1.;
  inv_Ly = -1.;
  inv_Lz = -1.;
  total_time = -1.;
  total_num_steps = -1;
  temperature = -1.;
  requested_temperature = -1.;
  pot_energy = -1.;
  kin_energy = -1.;
  energy = -1.;
  virSum = -1.;
  pressure = -1.;
}

ensemble::~ensemble(void) 
{ 
	delete [] atoms;
}

void ensemble::readFile(char* the_file, ofstream& warn)
{
  // THIS FUNCTION READS INITIAL VALUES FROM FILE

  ifstream inp(the_file);

  char a[30];

  inp >> density >> a
      >> delta_t >> a
      >> Lx >> a
      >> Ly >> a
      >> Lz >> a
      >> total_time >> a
      >> temperature >> a
	;

  // -------------------------------------------------------
  // CALCULATE INVERSE OF SIZES OF THE SYSTEM
  inv_Lx = 1.0/Lx;
  inv_Ly = 1.0/Ly;
  inv_Lz = 1.0/Lz;
  
  // -------------------------------------------------------
  // CALCULATE THE TOTAL NUMBER OF PARTICLES ONCE HAVE THE
  // DENSITY AND SIZE OF THE BOX.
  // BECAUSE N IS AN INTEGER AND DENSITY IS A DOUBLE,
  // N HAS TO BE ROUNDED, AND DENSITY HAS TO BE RECALCULATED
  // FOR THE REAL VALUE.
  double density_orig = density;
#if NDIM == 2
  N 	  = ANINT(density * Lx * Ly);
  density = double(N) / (Lx*Ly);
#else
  N 	  = ANINT(density * Lx * Ly * Lz);
  density = double(N) / (Lx*Ly*Lz);
#endif

  // REPORT IF COULD NOT MEET THE REQUIRED DENSITY
  if (density_orig != density)
	warn << " * * WARNING: Requested density could not be met. New: " 
	  << density
	  << " requested: "
	  << density_orig
	  << " where abs(new-old) = "
	  << ABS(density_orig - density)
	  << "\n";

  // CREATE ATOMS
  atoms = new atom[N];

  // -------------------------------------------------------
  // CALCULATE THE TOTAL NUMBER OF STEPS USING THE TOTAL TIME AND
  // TIME STEP
  total_num_steps = ANINT( total_time / delta_t );

  // -------------------------------------------------------
  // SAVE THE REQUESTED TEMPERATURE
  requested_temperature = temperature;
}

void ensemble::init_positions(int a)
{
  // INITIALIZE POSITIONS
  switch(a)
  {
	case 1:
	  init_lattice();
	  break;
	case 2:
	  init_random_pos_check_energy();
	  break;
	default:
	  cout << "WARNING: no position initialization given " << a << "\n";
  }
}

void ensemble::init_velocities(int b)
{
  // INITIALIZE VELOCITIES
  switch(b)
  {
	case 1:
	  init_random_vel();
	  break;
	default:
	  cout << "WARNING: no velocity initialization given " << b << "\n";
  }
}

void ensemble::init_lattice(void)
{
  // PUT ATOMS IN A LATTICE
  
  // HOW MANY ATOMS ON THE SIDE
#if NDIM == 2
  int num_atoms_side = int( sqrt( N ) + 1.0 );
#else
  int num_atoms_side = int( pow( N, 1./3. ) + 1.0 );
#endif
  // DISTANCE INCREMENTS
  double delx = Lx / num_atoms_side * 1.;
  double dely = Ly / num_atoms_side * 1.;
  double delz = Lz / num_atoms_side * 1.;

  for (int i=0; i<N; i++)
  {
	// GET LATTICE COORDINATES
	int coor_k = int( double(i) / (num_atoms_side*num_atoms_side) ) + 1 ;
	int rest   = i - (coor_k - 1) * (num_atoms_side*num_atoms_side);
#if NDIM == 2
	int coor_j = int( double(i) / num_atoms_side ) + 1;
	int rest2  = i - (coor_j - 1) * num_atoms_side;
#else
	int coor_j = int( double(rest) / num_atoms_side ) + 1;
	int rest2  = rest - (coor_j - 1) * num_atoms_side;
#endif
	int coor_i = rest2;

	// ASSIGN DISTANCES
#if NDIM == 2
	atoms[i].r = vec( coor_i*delx, coor_j*dely, 0. );
#else
	atoms[i].r = vec( coor_i*delx, coor_j*dely, coor_k*delz );
#endif
  }
}

void ensemble::init_random_pos_check_energy(void)
{
  // PUT PARTICLES AT RANDOM POSITIONS, BUT ONLY
  // ACCEPT POSITIONS THAT CONTIBUTE TO THE ENERGY
  // BASED ON A THRESHOLD VALUE

  // CRITERION FOR ACCEPTANCE
  // GOOD FOR 2D
  double threshold_energy = 5.5;
  // GOOD FOR 3D
  //double threshold_energy = 3.3;

  int i,j;
  vec trial_position;
  double trial_energy;
  double running_energy;
  double potential_energy = 0.;
  int accept;
  
  vec dr;
  double r2, inv_r2, inv_r6;
  double force_val;

  // GENERATE RANDOM POSITION FOR THE FIRST PARTICLE
#if NDIM == 2
  trial_position = vec( rand01()*Lx, rand01()*Ly, 0. );
#else
  trial_position = vec( rand01()*Lx, rand01()*Ly, rand01()*Lz );
#endif
  // ACCEPT THAT POSITION
  atoms[0].r = trial_position;

  // GENERATE FOR THE REST OF THE PARTICLES
  i=1;
  while (i<N)
  {
  	running_energy = 0.;

  	trial_position.reset();
#if NDIM == 2
  	trial_position = vec( rand01()*Lx, rand01()*Ly, 0. );
#else
  	trial_position = vec( rand01()*Lx, rand01()*Ly, rand01()*Lz );
#endif
	//PRINT3(i,"TRIAL POSITION",trial_position);

	// CHECK LJ POTENTIAL ENERGY CONTRIBUTION DUE TO THIS NEW
	// PARTICLE POSITION. REJECT IF LARGER THAN THRESHOLD.
	j=0;
	accept = 1;
	while (j<i)
	{
	  // CALCULATE DISTANCE BETWEEN ATOMS
	  dr = trial_position - atoms[j].r;
	  // USE IMAGE DISTANCES
	  dr.x = IMAGE(dr.x, inv_Lx, Lx);
	  dr.y = IMAGE(dr.y, inv_Ly, Ly);
	  dr.z = IMAGE(dr.z, inv_Lz, Lz);

	  r2 = dr.norm2();
	  inv_r2 = 1.0 / r2;
	  inv_r6 = inv_r2 * inv_r2 * inv_r2;
	  // LJ ENERGY IN REDUCED UNITS
	  trial_energy    = 4.*inv_r6*(inv_r6 - 1.);
	  running_energy += trial_energy;

	  if (trial_energy > threshold_energy)
	  {
		// REJECT POSITION DUE TO TOO HIGH ENERGY CONTRIBUTION.
		// SET CONDITION TO GET OUT OF THE LOOP EARLY.
		j = i;
		accept = 0;
	  }
	  j++;
	}

	// ACCEPT OR NOT THE NEW POSITION
	if (accept)
	{
	  // UPDATE POSITION IF ACCEPTED
  	  atoms[i].r = trial_position;
	  // ADVANCE TO THE NEXT PARTICLE
	  i++;

	  // ACCUMULATE CURRENT CONTRINBUTION TO THE
	  // TOTAL POTENTIAL ENERGY
  	  potential_energy += running_energy;

	  // REPORT ON PROGRESS
	  PRINT4("position progress = ",double(i)*100./double(N)," % with ave pot energy ",potential_energy/double(i));
	}
  }
}

void ensemble::init_random_vel(void)
{
  int i;
  // GIVE RANDOM INITIAL VELOCITIES

  vec vel_cm;
  double vel_cm2 = 0.; 

  for (i=0; i<N; i++)
  {
#if NDIM == 2
	atoms[i].v = vec( (rand01()-0.5), (rand01()-0.5), 0. );
#else
	atoms[i].v = vec( (rand01()-0.5), (rand01()-0.5), (rand01()-0.5) );
#endif
	vel_cm = vel_cm + atoms[i].v;
	vel_cm2 += atoms[i].v.norm2();
  }

  vel_cm  /= double(N);
  vel_cm2 /= double(N);

  // SCALE FACTOR
  double fs = sqrt( 3.0 * temperature / vel_cm2 ); 

  for (i=0; i<N; i++)
  {
	atoms[i].v = ( atoms[i].v - vel_cm ) * fs;
	// UPDATE THE OLD POSITION
	atoms[i].r_old = atoms[i].r - atoms[i].v * delta_t;
  }

  // TEST TEMPERATURE
  /*
  double test_T = 0.;
  for (i=0; i<N; i++)
  {
	test_T += atoms[i].v.norm2();
  }
  test_T /= (3.*N-3);
  PRINT2("hey check on T ",test_T);
  */
}

void ensemble::output_properties(int counter, double real_time, ofstream& out_f)
{
  //cout << "* * TEMP = " << temperature << "\n";
  out_f << setprecision(14)
	<< counter << " "
	<< real_time << " "
	<< temperature  << " "
	<< kin_energy / double(N) << " "
	<< pot_energy / double(N) << " "
	<< energy / double(N) << " "
	<< pressure << " "
	<< "\n";
}

void ensemble::output_positions(void)
{
  // WRITE ALL THE ATOMS
  for (int i=0; i<N; i++)
  {
	   atoms[i].output_position();
  }
}

void ensemble::output_positions(char* outfile)
{
  ofstream out_f(outfile);

  // WRITE ALL THE ATOMS
  for (int i=0; i<N; i++)
  {
	   out_f << atoms[i].r << "\n";
  }
  out_f.close();
}

void ensemble::output_position(int a)
{
  // WRITE ONLY FOR GIVEN ATOM
  cout << atoms[a].r << "\n";
}

void ensemble::output_velocities(void)
{
  // WRITE ALL THE ATOMS
  for (int i=0; i<N; i++)
  {
	   atoms[i].output_velocity();
  }
}

void ensemble::output_velocities(char* outfile)
{
  ofstream out_f(outfile);

  // WRITE ALL THE ATOMS
  for (int i=0; i<N; i++)
  {
	   out_f << atoms[i].v << "\n";
  }
  out_f.close();
}

void ensemble::output_speeds(char* outfile)
{
  ofstream out_f(outfile);

  // WRITE ALL THE ATOMS
  for (int i=0; i<N; i++)
  {
	   out_f << atoms[i].v.norm() << "\n";
  }
  out_f.close();
}

void ensemble::output_total_linear_momentum(void)
{
  // WRITE TOTAL LINEAR MOMENTUM
  vec big_P;
  for (int i=0; i<N; i++)
  {
	big_P = big_P + atoms[i].v;
  }
  PRINT2("TOTAL MOMENTUM = ",big_P);
}

void ensemble::output_total_linear_momentum(ofstream& out_p)
{
  // WRITE TOTAL LINEAR MOMENTUM
  vec big_P;
  for (int i=0; i<N; i++)
  {
	big_P = big_P + atoms[i].v;
  }
  out_p << big_P << "\n";
}

void ensemble::clear_forces(void)
{
  for (int i=0; i<N; i++)
	atoms[i].f.reset();
}

void ensemble::clear_forces(int ii)
{
	atoms[ii].f.reset();
}

double ensemble::calc_temperature(void)
{
  int i;
  // CALCULATE INSTANTANEOUS TEMPERATURE

  temperature = 0.;
  for (i=0; i<N; i++)
  {
	temperature += atoms[i].v.norm2();
  }
#if NDIM==2
  return temperature / (2.*N-2.);
#else
  return temperature / (3.*N-3.);
#endif
}

void ensemble::calc_pressure(void)
{
#if NDIM == 2
  pressure = density * temperature + virSum / (Lx*Ly*Lz*2);
#else
  pressure = density * temperature + virSum / (Lx*Ly*Lz*3);
#endif
}

/*
 * ----------------------------------------------------------------------
 */
