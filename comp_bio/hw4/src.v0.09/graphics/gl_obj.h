
/*
 * Luis Cruz
 * Center for Polymer Studies
 * Boston University
 * ccruz@bu.edu
 * April 2001 - July, 2001
 *
 */   

#ifndef GL_OBJ_H
#define GL_OBJ_H

#include "config_fltk.h"
#include <fstream>

#include "GL/gl.h"

#ifndef BACKGROUND
# include "Fl/Fl.H"
# include "Fl/Fl_Gl_Window.H"
#endif

#include "defs.h"

#include "atom.h"
#include "ensemble.h"
#include "integrator.h"
#include "minima.h"
#include "forces.h"
#include "diffusion.h"
#include "paircorrfunc.h"
#include "R.h"
#include "verletlist.h"
#include "stepwise.h"
#include "langevin.h"
#include "montecarlo.h"


class gl_obj : public Fl_Gl_Window
{
//		friend ostream &operator<<(ostream &os, gl_obj &z);
//		friend istream &operator>>(istream &os, gl_obj &z); 
public:
	gl_obj(int x,int y,int w,int h,const char *l=0);
	~gl_obj(void);

	void clean(void);

	void DoDisplayLists(void);
	void DisplayLists(void);

	void save_all(char *, gl_obj&);
	void read_all(char *, gl_obj&);

    GLuint blue_circle_DL;
    GLuint red_circle_DL;
    GLuint white_line_tag_DL;

    GLuint red_radius_circle_DL;
    GLuint green_diam_circle_DL;
    GLuint green_cross_circle_DL;

    GLuint magenta_full_circle_DL;

    // =====================

    void  v_angle(double angle){vAng=angle;};
    double v_angle(){return vAng;};

    void  h_angle(double angle){hAng=angle;};
    double h_angle(){return hAng;};

    void  panx(double x){xshift=x;};
    double panx(void){return xshift;};
    
    void  pany(double y){yshift=y;};
    double pany(void){return yshift;};

    /*The widget class draw() override.
     *
     * The draw() function initialize Gl for another round of drawing
     * then calls specialized functions for drawing each of the
     * entities displayed in the cube view.
     *
     */
    void draw();    
	int handle(int event);

  // =============================================
  // CREATE AND INITIALIZE:
  // - THE SYSTEM
  ensemble sys;
  // - ACTION OBJECT
  forces physics;
  // - FINITE DIFF. INTEGRATOR
  integrator integrate;
  // - LANGEVIN INTEGRATOR
  langevin brownian;
  // - MINIMIZATION OBJECT
  minima minimize;
  // - STEPWISE INTERACTION ALGORITHM
  stepwise hardsphere;
  // - MC ALGORITHM
  montecarlo mc;

  verletlist *neighborsverlet;

  diffusion *diff;
  paircorrfunc *pcf;
  R_meas *R;

  double time_meas;
  int nhis;
  long n;

private:
    void  drawCanvas();
	void  show_credits(void);
    
    double vAng,hAng;
    double xshift,yshift;
    double size;

	int   window_width, window_height;

    // =====================
};

#endif // GL_OBJ_H
