
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

#ifndef MACROS_H
#define MACROS_H

// PRINTING MACROS
#define PRINT1(x) \
  cout << (x) << "\n";
#define PRINT2(x,y) \
  cout << (x) << " " << (y) << "\n" ;
#define PRINT3(x,y,z) \
  cout << (x) << " " << (y) << " " << (z) << "\n" ;
#define PRINT4(x,y,z,w) \
  cout << (x) << " " << (y) << " " << (z) << " " << (w) << "\n" ;
#define PRINT5(x,y,z,w,v) \
  cout << (x) << " " << (y) << " " << (z) << " " << (w) << " " << (v) << "\n" ;
#define PRINT6(x,y,z,w,v,u) \
  cout << (x) << " " << (y) << " " << (z) << " " << (w) << " " << (v) << " " << (u) << "\n" ;
#define PRINT7(x,y,z,w,v,u,h) \
  cout << (x) << " " << (y) << " " << (z) << " " << (w) << " " << (v) << " " << (u) << " " << (h) << "\n" ;



// USEFUL CONSTANTS

#ifndef PI
#define PI       3.14159265358979323846
#endif

// USEFUL UTILITIES

#define SIGN(x)  ( ((x)<0) ? (-1) : (1) )
#define ABS(x)   ( ((x)<0) ? (-(x)) : (x))
#define MAX(x,y) ( (((x) > (y)) ? (x) : (y)) )
#define MIN(x,y) ( (((x) < (y)) ? (x) : (y)) )
#define CUAD(x)  ( (x) * (x) )

// -- ROUNDING FUNCTION
#define ANINT(a) \
	( int( (a) + SIGN(a)*(.5) ) )

// -- MINIMUM IMAGE DISTANCE MACRO
#define IMAGE(r, boxinv, box) \
        ( (r) - ANINT( (r) * (boxinv) ) * (box) )

// -- RETURN ANGLE IN DEGREES
# define DEGREE(theta) \
        ( (theta) * 180. / PI )
// -- RETURN ANGLE IN RADIANS
# define RAD(theta) \
        ( (theta) * PI / 180. )


#endif // MACROS_H
