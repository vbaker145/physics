
/*
 * Luis Cruz
 * Center for Polymer Studies
 * Boston University
 * ccruz@bu.edu
 * April 2001 - July, 2001
 *
 */   

#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "disp_func.h"
//#include "macros.h"
//
using namespace std;

float *rainbow(float radix, float* value)
      // GIVEN A VALUE radix IT IS MAPPED INTO A
      // RAINBOW. RETURNS AN {R,G,B} VECTOR.
      // 0.0 =< radix <= 1.0
{
      float r, g, b;

	  // NORMALIZE AND SHIFT THE GIVEN VALUE SO THAT COLOR
	  // ARE SHIFTED UPWARDS (SO THAT RED IS NOT LOW AND HIGH
	  // AT THE SAME TIME)
	  radix = radix * float(0.899) + float(0.1);

      if ( (radix<0.) || (radix>1.0) )
      {
          cout << "out of bounds in rainbow -- " << radix << "\n";
          exit(1);
      }

      // CHECK IN WHAT INTERVAL OF THE RAINBOX radix IS.
      int interval  = int(radix * 6.0);
      float residue = radix * 6.0 - float(interval);

      switch(interval)
      {
          case 0:
          {
             r = 1.0;
             g = residue;
             b = 0.0;
             break;
          }
          case 1:
          {
             r = 1.0-residue;
             g = 1.0;
             b = 0.0;
             break;
          }
          case 2:
          {
             r = 0.0;
             g = 1.0;
             b = residue;
             break;
          }
          case 3:
          {
             r = 0.0;
             g = 1.0-residue;
             b = 1.0;
             break;
          }
          case 4:
          {
             r = residue;
             g = 0.0;
             b = 1.0;
             break;
          }
          case 5:
          {
             r = 1.0;
             g = 0.0;
             b = 1.0-residue;
             break;
          }
          case 6:
          {
             r = 1.0;
             g = 0.0;
             b = 0.0;
             break;
          }
          default:
          {
             cout << "there is no default in rainbow " << interval << "\n";
             exit(1);
             break;
          }
      }

      value[0] = r;
      value[1] = g;
      value[2] = b;

      return value;
}

float ramp_func(float x)
{
/*
      float func_x = 1./ (1. + exp( dq.beta *(dq.mu-x ) ));
      float func_1 = 1./ (1. + exp( dq.beta *(dq.mu-1.) ));
      float dum2 = 1. - func_1 + func_x;
*/
	float beta,mu;
	// DEFAULTS
        beta = 4.;
        mu   = .5;

#ifdef LINUX
    beta = 4.;
    mu   = .5;
#endif
#ifdef ALPHA4
    beta = 4.;
//  beta = 10.;
//  beta = 15.;
    mu   = .5;
#endif
#ifdef ONYX
    beta = 1.;
    mu   = .5;
#endif
//#if defined(wx_msw)
//    beta = 4.;
//    mu   = .5;
//#endif

      float func_x = 1./ (1. + exp( beta *(mu-x ) ));
      float func_1 = 1./ (1. + exp( beta *(mu-1.) ));
      float func_0 = 1./ (1. + exp( beta *(mu-0.) ));
      float A      = 1./(func_1-func_0);
      float dum2   = A * ( func_x - func_0 );

      return dum2;
}

/*
 * --------------------------------------------------------------
 */
