
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
     << z.do_polymers  << " "
     << z.num_chains  << " "
     << z.atoms_per_chain  << " "
     << z.length_bonds  << " "
     << z.total_time  << " "
     << z.total_num_steps  << " "
     << z.temperature  << " "
     << z.pot_energy  << " "
     << z.kin_energy  << " "
     << z.energy  << " "
     << z.virSum  << " "
     << z.pressure  << " "
     << z.num_not_fixed  << " "
     << z.init_pos  << " "
     << z.threshold_energy  << " "
     << z.init_vel  << " "
     << z.dynamics  << " "
     << z.damping_constant  << " "
     << z.time_meas  << " "
     << z.nhis_pcf  << " "
     << z.do_vacf  << " "
     << z.do_pcf  << " "
     << z.do_R  << " "
     << z.do_tag  << " "
     << z.tag_particle  << " "
     << z.mc_dr  << " "
     << z.do_mc_adjust  << " "
     << z.rc  << " "
     << z.rl  << " "
     << z.XY_J  << " "
     << z.theta_o  << " "
     << z.do_minimize  << " "
     << z.two_part_vel  << " "
	<< "\n"
	;

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
  run_pressure = -1.;
  num_not_fixed = -1;
  init_pos = -1;
  threshold_energy = -1;
  init_vel = -1;
  dynamics = -1;
  damping_constant = -1;
  time_meas = -1;
  nhis_pcf = -1;
  do_vacf = 0;
  do_pcf = 0;
  do_R = 0;

  do_tag = 0;
  tag_particle = -1;
  reset_tag_pos();

  mc_dr = 0.;
  do_mc_adjust = 0;
  rc = -1.;
  rl = -1.;
  XY_J = -1.;
  theta_o = 1.;
  do_minimize = 0;
  two_part_vel = .4941;

  do_polymers = 0;
  num_chains = 2;
  atoms_per_chain = 6;
  length_bonds = 1.;
}

ensemble::~ensemble(void) 
{ 
	delete [] atoms;
}

void ensemble::readFile(char* the_file)
{
  // THIS FUNCTION READS INITIAL VALUES FROM FILE
  ifstream inp(the_file);

  char a[30];

  inp >> density >> a
      >> delta_t >> a
      >> Lx >> a
      >> Ly >> a
      >> Lz >> a
      >> do_polymers >> a
      >> num_chains >> a
      >> atoms_per_chain >> a
      >> length_bonds >> a
      >> total_time >> a
      >> temperature >> a
      >> pressure >> a
      >> init_pos >> a
      >> threshold_energy >> a
      >> init_vel >> a
      >> dynamics >> a
      >> damping_constant >> a
      >> time_meas >> a
      >> nhis_pcf >> a
      >> do_vacf >> a
      >> do_pcf >> a
      >> do_R >> a
      >> do_tag >> a
      >> tag_particle >> a
      >> mc_dr >> a
      >> do_mc_adjust >> a
      >> rc >> a
      >> rl >> a
      >> XY_J >> a
      >> theta_o >> a
      >> do_minimize >> a
      >> two_part_vel >> a
	;
/*
  cout << density << " "
      << delta_t << " "
      << Lx << " "
      << Ly << " "
      << Lz << " "
      << total_time << " "
      << temperature << " "
      << init_pos << " "
      << threshold_energy << " "
      << init_vel << " "
      << dynamics << " "
      << damping_constant << " "
      << time_meas << " "
      << nhis_pcf << " "
      << do_vacf << " "
      << do_pcf << " "
      << do_R << " "
      << do_tag << " "
      << tag_particle << " "
      << mc_dr << " "
      << do_mc_adjust << " "
      << rc << " "
      << rl << " "
      << XY_J << " "
      << theta_o << " "
      << do_minimize << " "
	  << "\n"
	;
*/
  inp.close();
  //cout << "out of read\n";
}

void ensemble::writeFile(char* the_file)
{
  // THIS FUNCTION WRITES THE INITIAL VALUES TO FILE
  ofstream inp(the_file);

  inp << density << "\trequested_particle_density\n"
      << delta_t << "\tdelta_t\n"
      << Lx << "\tLx\n"
      << Ly << "\tLy\n"
      << Lz << "\tLz\n"
      << do_polymers << "\tdo_polymers\n"
      << num_chains << "\tnum_chains\n"
      << atoms_per_chain << "\tatoms_per_chain\n"
      << length_bonds << "\tlength_bonds\n"
      << total_time << "\ttotal_time\n"
      << temperature << "\ttemperature\n"
      << pressure << "\tpressure\n"
      << init_pos << "\tinit_pos\n"
      << threshold_energy << "\tthresh_energ\n"
      << init_vel << "\tinit_vel\n"
      << dynamics << "\tdynamics\n"
      << damping_constant << "\tdamping_constant\n"
      << time_meas << "\ttime_measurement\n"
      << nhis_pcf << "\tnhis_pcf_bins\n"
      << do_vacf << "\tdo_vacf\n"
      << do_pcf << "\tdo_pcf\n"
      << do_R << "\tdo_R\n"
      << do_tag << "\tdo_tag\n"
      << tag_particle << "\ttag_particle\n"
      << mc_dr << "\tmc_dr\n"
      << do_mc_adjust << "\tdo_mc_adjust\n"
      << rc << "\trc_cutoff\n"
      << rl << "\trl_cutoff\n"
      << XY_J << "\tXY_J\n"
      << theta_o << "\tdelta_theta\n"
      << do_minimize << "\tdo_minimize\n"
      << two_part_vel << "\ttwo_part_vel\n"
	;

  inp.close();
}

void ensemble::apply_inits(ofstream& warn)
{
  // USE THE PARAMETERS JUST READ FROM FILE
  // TO INITIALIZE THE SYSTEM

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
  // LIQUIDS
  N = ANINT(density * Lx * Ly);
  // CHECK IF THIS IS A TWO-PARTICLE SIMULATION
  if (init_pos==2 || init_pos==3) N = 2;  
  // POLYMERS
  if (do_polymers) 
	N = num_chains * atoms_per_chain;
  density = double(N) / (Lx*Ly);
#else
  // LIQUIDS
  N = ANINT(density * Lx * Ly * Lz);
  // CHECK IF THIS IS A TWO-PARTICLE SIMULATION
  if (init_pos==2 || init_pos==3) N = 2;  
  // POLYMERS
  if (do_polymers) 
	N = num_chains * atoms_per_chain;
  density = double(N) / (Lx*Ly*Lz);
#endif
  // AT THIS POINT ALL ATOMS ARE MOBILE
  num_not_fixed = N;

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

  // -------------------------------------------------------
  // INITIALIZE INSTANTANEOUS PRESSURE
  run_pressure = pressure;
}

void ensemble::init_positions(int a)
{
  // INITIALIZE POSITIONS
  switch(a)
  {
	case 0:
	  init_lattice();
	  break;
	case 1:
	  init_random_pos_check_energy();
	  break;
	case 2:
	  init_random_pos_two_particles();
	  break;
	case 3:
	  init_xtal_pos_two_particles();
	  break;
	case 4:
	  init_xtal_polymer();
	  break;
	case 5:
	  init_from_file();
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
	case 0:
	  init_zero_vel();
	  break;
	case 1:
	  init_random_vel();
	  break;
	case 2:
	  init_two_part_vel();
	  break;
	case 3:
	  init_polymers_vel();
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
  //double threshold_energy = 7.5;
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

void ensemble::init_random_pos_two_particles(void)
{
  // PUT ONLY TWO PARTICLES AT RANDOM POSITIONS.

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

  //float eps = 1.1;
  float eps = 5.;

  // PUT THE OTHER ONE
  i=1;
  while (i<2)
  {
  	running_energy = 0.;

  	trial_position.reset();
#if NDIM == 2
  	//trial_position = vec( rand01()*Lx, rand01()*Ly, 0. );
  	trial_position = vec( atoms[0].r.x+eps, atoms[0].r.y+eps, 0. );
#else
  	//trial_position = vec( rand01()*Lx, rand01()*Ly, rand01()*Lz );
  	trial_position = vec( atoms[0].r.x+eps, atoms[0].r.y+eps,  atoms[0].r.z+eps);
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

void ensemble::init_xtal_pos_two_particles(void)
{
  // PUT ONLY TWO PARTICLES AT RANDOM POSITIONS.

  int i,j;
  vec trial_position;
  
  // GENERATE RANDOM POSITION FOR THE FIRST PARTICLE
  //trial_position = vec( 0., 0., 0. );
  trial_position = vec( 0.0 , 1., 0. );
  // ACCEPT THAT POSITION
  atoms[0].r = trial_position;

  trial_position.reset();
  trial_position = vec( 1.9 , 1., 0. );
  atoms[1].r = trial_position;

}

void ensemble::init_xtal_polymer(void)
{
  // PUT ATOMS IN A POLYMER LATTICE
  
  // HOW MANY ATOMS ON THE SIDE
  int num_atoms_side = atoms_per_chain;
  // DISTANCE INCREMENTS
  double delx = Lx / num_atoms_side * 1.;
  double dely = length_bonds;
  double delz = Lz / num_atoms_side * 1.;

  double init_x = .1;
  double init_y = 2.;
  double init_z = .1;
  double pos_x, pos_y, pos_z;
  int id;

  double shiftx = Lx / 2. - 1.;
  //double shifty = Ly / 2. - 1.;
  double shifty = 1.;
  double shiftz = Lz / 2. - 1.;

  for (int i=0; i<num_chains; i++)
  {
	pos_y = i * init_y;
#if NDIM == 2
	pos_z = 0.;
#else
	pos_z = init_z;
#endif
	for (int j=0; j<atoms_per_chain; j++)
	{
	  pos_x = j * length_bonds;
	// ASSIGN DISTANCES
	id = i*atoms_per_chain + j;
#if NDIM == 2
	atoms[id].r = vec( pos_x-shiftx, pos_y-shifty, 0. );
#else
	atoms[id].r = vec( pos_x-shiftx, pos_y-shifty, pos_z-shiftz );
#endif

	}
  }
}

void ensemble::init_from_file(void)
{
  // READ ATOM POSITIONS FROM FILE
  char file_in[30] = "inputs/positions.xyz";
  input_positions(file_in);
}

void ensemble::init_zero_vel(void)
{
  int i;

  // ZERO VELOCITIES
  for (i=0; i<N; i++)
  {
	atoms[i].v.reset();
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
	// NO VELOCITY FOR FIXED PARTICLES
	if (atoms[i].is_not_fixed())
	{
#if NDIM == 2
	  atoms[i].v = vec( (rand01()-0.5), (rand01()-0.5), 0. );
#else
	  atoms[i].v = vec( (rand01()-0.5), (rand01()-0.5), (rand01()-0.5) );
#endif
	  vel_cm = vel_cm + atoms[i].v;
	  vel_cm2 += atoms[i].v.norm2();
	}
	else
	  atoms[i].v.reset();
  }

  //vel_cm  /= double(N);
  vel_cm  /= double(num_not_fixed);
  //vel_cm2 /= double(N);
  vel_cm2 /= double(num_not_fixed);

  // SCALE FACTOR
  double fs = sqrt( 3.0 * temperature / vel_cm2 ); 

  for (i=0; i<N; i++)
  {
	if (atoms[i].is_not_fixed())
	{
	  atoms[i].v = ( atoms[i].v - vel_cm ) * fs;
	  // UPDATE THE OLD POSITION
	  atoms[i].r_old = atoms[i].r - atoms[i].v * delta_t;
	  /*
	  cout << " first " << atoms[i].r_old 
		<< " * " << atoms[i].r 
		<< " * " << atoms[i].v 
		<< " * " <<  delta_t << "\n";
	  */
	}
	/*
	  cout << "* " << i << " " 
		<< atoms[i] << "\n";
	cout << " RANDOM VEL POS " << atoms[i].r  << " ****** " << atoms[i].r_old << "\n";
	*/
  }

}

void ensemble::init_two_part_vel(void)
{
  int i;
  // GIVE INITIAL VELOCITIES
  // TO A TWO-PARTICLE SYSTEM

	  atoms[0].v.reset();
	  atoms[1].v.reset();

  float base = 0.4;
  // 1 - JUST ESCAPING
  //float eps = .1;
  // 2 - ORBIT THEN ESCAPE
  //float eps = .095;
  //float eps = .09417;
  // 3 - ORBIT - BUMP
  float eps = .0941;
  //float eps = .094;
  //float eps = .090;

  //atoms[0].v = vec( 0., -base-eps, 0. );
  atoms[0].v = vec( 0., -this->two_part_vel, 0. );
  //atoms[1].v = vec( 0., base+eps, 0. );
  atoms[1].v = vec( 0., this->two_part_vel, 0. );

  // SETUP OLD POSITIONS BASED ON VELOCITIES
  atoms[0].r_old = atoms[0].r - atoms[0].v * delta_t;
  atoms[1].r_old = atoms[1].r - atoms[1].v * delta_t;
}

void ensemble::init_polymers_vel(void)
{
  int i;
  // GIVE INITIAL VELOCITIES
  // TO MONOMERS IN A POLYMER

  vec vel_cm;
  double vel_cm2 = 0.; 

  for (i=0; i<N; i++)
  {
	// NO VELOCITY FOR FIXED PARTICLES
	if (atoms[i].is_not_fixed())
	{
#if NDIM == 2
	  if (i%2==0)
	  	atoms[i].v = vec( -0.5, -0.5, 0. );
	  else
	  	atoms[i].v = vec( +0.5, +0.5, 0. );
#else
	  if (i%2==0)
	  	atoms[i].v = vec( -0.5, -0.5, -0.5 );
	  else
	  	atoms[i].v = vec( +0.5, +0.5, +0.5 );
#endif
	  vel_cm = vel_cm + atoms[i].v;
	  vel_cm2 += atoms[i].v.norm2();
	}
	else
	  atoms[i].v.reset();
  }

  //vel_cm  /= double(N);
  vel_cm  /= double(num_not_fixed);
  //vel_cm2 /= double(N);
  vel_cm2 /= double(num_not_fixed);

  // SCALE FACTOR
  double fs = sqrt( 3.0 * temperature / vel_cm2 ); 

  for (i=0; i<N; i++)
  {
	if (atoms[i].is_not_fixed())
	{
	  atoms[i].v = ( atoms[i].v - vel_cm ) * fs;
	  // UPDATE THE OLD POSITION
	  atoms[i].r_old = atoms[i].r - atoms[i].v * delta_t;
	  /*
	  cout << " first " << atoms[i].r_old 
		<< " * " << atoms[i].r 
		<< " * " << atoms[i].v 
		<< " * " <<  delta_t << "\n";
	  */
	}
	/*
	  cout << "* " << i << " " 
		<< atoms[i] << "\n";
	cout << " POLYMER VEL POS " << atoms[i].r  << " ****** " << atoms[i].r_old << "\n";
	*/
  }

}

void ensemble::scale_vel(void)
{
  int i;
  // SCALE ALL VELOCITIES BY SQRT(TO/T)

  //double factor = 0.1;
  double factor = 1.0;
  double scale_v = factor * sqrt(requested_temperature/temperature);

  cout << "Scaling:\n";
  for (i=0; i<N; i++)
  {
	// NO VELOCITY FOR FIXED PARTICLES
	if (atoms[i].is_not_fixed())
	{
	  atoms[i].v *= scale_v;
	}
	cout << i << " ";
  }
  cout << "\n";
}

void ensemble::init_random_angles(void)
{
  int i;
  // GIVE RANDOM INITIAL ANGLES
  for (i=0; i<N; i++)
  {
	  atoms[i].theta = 2.*PI*rand01();
	  atoms[i].phi   =    PI*rand01();
  }
}

void ensemble::output_properties(int counter, double real_time, ofstream& out_f)
{
  //cout << "* * TEMP = " << temperature << "\n";
  out_f << setprecision(14)
	<< counter << " "
	<< real_time << " "
	<< temperature  << " "
	//<< kin_energy / double(N) << " "
	<< kin_energy / double(num_not_fixed) << " "
	//<< pot_energy / double(N) << " "
	<< pot_energy / double(num_not_fixed) << " "
	//<< energy / double(N) << " "
	<< energy / double(num_not_fixed) << " ";

	// DONT KNOW HOW TO CALCULATE CORRECT VIRIAL PRESSURE
	// IN NPT
	if (dynamics!=9)
	out_f << run_pressure << " ";
    else
	out_f << 0 << " ";

  // MC WITH ANGLES
  if ( dynamics == 7 )
	out_f << order_param_LJXY() << " ";
  else if ( dynamics == 8 )
	out_f << order_param_LJLL() << " ";

  // MC AT NPT
  if ( dynamics == 9 )
	out_f << density << " ";

	  out_f << "\n";

	  /*
  cout << setprecision(14)
	<< counter << " "
	<< real_time << " "
	<< temperature  << " "
	//<< kin_energy / double(N) << " "
	<< kin_energy / double(num_not_fixed) << " "
	//<< pot_energy / double(N) << " "
	<< pot_energy / double(num_not_fixed) << " "
	//<< energy / double(N) << " "
	<< energy / double(num_not_fixed) << " "
	<< run_pressure << " ";

  // MC WITH ANGLES
  if ( dynamics == 7 || dynamics == 8 )
	cout << order_param_LJXY() << " ";

  // MC AT NPT
  if ( dynamics == 9 )
	cout << density << " ";

	  cout << "\n";

	   */

}

void ensemble::output_positions(void)
{
  vec r;

  // WRITE ALL THE ATOMS
  for (int i=0; i<N; i++)
  {
	//atoms[i].output_position();
	r = atoms[i].r;
	// CALCULATE IMAGE POSITION
	r.x = IMAGE(r.x, inv_Lx, Lx);
	r.y = IMAGE(r.y, inv_Ly, Ly);
	r.z = IMAGE(r.z, inv_Lz, Lz);

	cout << r << "\n";
  }
}

void ensemble::output_positions(int write_image, char* outfile)
{
  // WRITE POSITIONS OF ALL ATOMS TO FILE.
  // SPECIFY IF WANT REAL write_image = 0
  // OR WRAPPED DISTANCES write_image = 1
  vec r;
  ofstream out_f(outfile);

  // WRITE NUMBER OF PARTICLES
  //out_f << N << "\n";

  // WRITE ALL THE ATOMS
  for (int i=0; i<N; i++)
  {
	if (write_image==0)
	{
	  out_f << atoms[i].r << "\n";
	}
	else
	{
	  r = atoms[i].r;
	  // CALCULATE IMAGE POSITION
	  r.x = IMAGE(r.x, inv_Lx, Lx);
	  r.y = IMAGE(r.y, inv_Ly, Ly);
	  r.z = IMAGE(r.z, inv_Lz, Lz);

	  out_f << r << "\n";
	}
  }
  out_f.close();
}

void ensemble::input_positions(char* infile)
{
  // READ POSITIONS FROM FILE
  int num_N;
  ifstream in_f(infile);

  // CHECK THAT FILE HAS EQUAL NUMBER
  // OF PARTICLES AS THE CURRENT SIMULATION
  in_f >> num_N;

  if (num_N!=N)
  {
	cout << "ERROR: Particles in file do not match simulation "
	  << num_N << " vs " << N << "\n";
	return;
  }

  // READ FOR ALL THE ATOMS
  double dum_rx, dum_ry, dum_rz;
  for (int i=0; i<N; i++)
  {
	in_f >> dum_rx;
	in_f >> dum_ry;
	in_f >> dum_rz;

	atoms[i].r = vec(dum_rx, dum_ry, dum_rz);
  }
  in_f.close();
}

void ensemble::output_position(int a)
{
  vec r;
  // WRITE ONLY FOR GIVEN ATOM
  //cout << atoms[a].r << "\n";

  r = atoms[a].r;
  // CALCULATE IMAGE POSITION
  r.x = IMAGE(r.x, inv_Lx, Lx);
  r.y = IMAGE(r.y, inv_Ly, Ly);
  r.z = IMAGE(r.z, inv_Lz, Lz);

  cout << r << "\n";
}

vec ensemble::get_image_position(int a)
{
  vec r;
  // WRITE ONLY FOR GIVEN ATOM
  //cout << atoms[a].r << "\n";

  r = atoms[a].r;
  // CALCULATE IMAGE POSITION
  r.x = IMAGE(r.x, inv_Lx, Lx);
  r.y = IMAGE(r.y, inv_Ly, Ly);
  r.z = IMAGE(r.z, inv_Lz, Lz);

  return r;
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

  // WRITE V FOR ALL THE ATOMS
  // INCLUDE ALSO ABS VAL
  for (int i=0; i<N; i++)
  {
	   //out_f << atoms[i].v << "\n";
	   out_f << atoms[i].v << " " << atoms[i].v.norm() << "\n";
  }
  out_f.close();
}

void ensemble::output_total_linear_momentum(void)
{
  // WRITE TOTAL LINEAR MOMENTUM
  vec big_P;
  for (int i=0; i<N; i++)
  {
	// EXCLUDE FIXED ATOMS
	if (atoms[i].is_not_fixed())
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
	// EXCLUDE FIXED ATOMS
	if (atoms[i].is_not_fixed())
		big_P = big_P + atoms[i].v;
  }
  out_p << big_P << "\n";
}

double ensemble::output_pressure(void)
{
  return run_pressure;
}

double ensemble::output_temperature(void)
{
  return temperature;
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
	// EXCLUDE FIXED ATOMS
	if (atoms[i].is_not_fixed())
		temperature += atoms[i].v.norm2();
  }
#if NDIM==2
  //return temperature / (2.*N-2.);
  return temperature / (2.*num_not_fixed-2.);
#else
  //return temperature / (3.*N-3.);
  return temperature / (3.*num_not_fixed-3.);
#endif
}

void ensemble::calc_pressure(void)
{
#if NDIM == 2
  run_pressure = density * temperature + virSum / (Lx*Ly*Lz*2);
  //PRINT7("*",density,temperature,virSum,Lx,Ly,Lz);
#else
  run_pressure = density * temperature + virSum / (Lx*Ly*Lz*3);
#endif
}

void ensemble::fix_atoms_beyond_radius(double radius)
{
  // FIXES ATOMS THAT ARE OUTSIDE A CIRCLE OF RADIUS radius
  int i,j;
  
  vec r_im;
  double radius2 = radius * radius;
  double r2;

  for (i=0; i<N; i++)
  {
	  r_im = atoms[i].r;
	  // USE IMAGE POSITION
	  r_im.x = IMAGE(r_im.x, inv_Lx, Lx);
	  r_im.y = IMAGE(r_im.y, inv_Ly, Ly);
	  r_im.z = IMAGE(r_im.z, inv_Lz, Lz);

	  r2 = r_im.norm2();

	  if (r2>=radius2)
	  {
		atoms[i].set_fixed();
		// UPDATE NUMBER OF MOBILE ATOMS
		num_not_fixed--;
	  }
  }
}

void ensemble::reset_tag_pos(void)
{
  tag_particle_pos.clear();
}

void ensemble::add_tag_pos(vec &r)
{
  tag_particle_pos.add_one(r);
}

void ensemble::record_tag_particles_position(void)
{
  vec im_r = get_image_position(tag_particle);
  //add_tag_pos( atoms[tag_particle].r );
  add_tag_pos( im_r );
}

void ensemble::write_tag_pos(void)
{
  // APPEND RESULTS TO FILE
  ofstream out_tag("results/tag.dat",ios_base::app);

  out_tag << "\n";
  for (int i=0; i<=tag_particle_pos.last_one; i++)
  {
	out_tag << tag_particle_pos[i] << "\n";
  }
  out_tag.close();
}

double ensemble::order_param_LJXY(void)
{
  int i;
  double order_param = 0.;
  double dum_x, dum_y;

#if NDIM==3
  PRINT1("*** NO IMPLEMENTATION OF ORDER PARAM IN 3D!");
#else
  for (i=0; i<N; i++)
  {
	if ( atoms[i].is_not_fixed() )
	{
	  dum_x += cos( atoms[i].theta );
	  //dum_x += CUAD( cos( atoms[i].theta ) );
	  //dum_x += cos(2.* atoms[i].theta );

	  dum_y += sin( atoms[i].theta );
	  //dum_y += CUAD( sin( atoms[i].theta ) );
	  //dum_y += sin(2.* atoms[i].theta );
	}
  }
  dum_x /= num_not_fixed;
  dum_y /= num_not_fixed;

  order_param = sqrt( CUAD(dum_x) + CUAD(dum_y) );
#endif
  return order_param;
}

double ensemble::order_param_LJLL(void)
{
  int i;
  double order_param = 0.;
  double dum_x, dum_y;

#if NDIM==3
  PRINT1("*** NO IMPLEMENTATION OF ORDER PARAM IN 3D!");
#else
  for (i=0; i<N; i++)
  {
	if ( atoms[i].is_not_fixed() )
	{
	  //dum_x += cos( atoms[i].theta );
	  //dum_x += CUAD( cos( atoms[i].theta ) );
	  dum_x += cos(2.* atoms[i].theta );

	  //dum_y += sin( atoms[i].theta );
	  //dum_y += CUAD( sin( atoms[i].theta ) );
	  dum_y += sin(2.* atoms[i].theta );
	}
  }
  dum_x /= num_not_fixed;
  dum_y /= num_not_fixed;

  order_param = sqrt( CUAD(dum_x) + CUAD(dum_y) );
#endif
  return order_param;
}

void ensemble::adjust_box_sizes(double lx, double ly, double lz)
{
  // -------------------------------------------------------
  // READJUST BOX LENGTHS
  Lx = lx;
  Ly = ly;
  Lz = lz;

  // -------------------------------------------------------
  // RECALCULATE INVERSE OF SIZES OF THE SYSTEM
  inv_Lx = 1.0/Lx;
  inv_Ly = 1.0/Ly;
  inv_Lz = 1.0/Lz;
}

void ensemble::calc_density(void)
{
#if NDIM == 2
  density = double(N) / (Lx*Ly);
#else
  density = double(N) / (Lx*Ly*Lz);
#endif
}

void ensemble::scale_positions(double sx, double sy, double sz)
{
  for (int i=0; i<N; i++)
  {
	atoms[i].r.x *= sx;
	atoms[i].r.y *= sy;
	atoms[i].r.z *= sz;
  }
}

/*
 * ----------------------------------------------------------------------
 */
