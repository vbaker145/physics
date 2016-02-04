
/*
 * Luis Cruz
 * Center for Polymer Studies
 * Boston University
 * ccruz@bu.edu
 * April 2001 - July, 2001
 *
 */   

#include "gl_obj.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iomanip>
# include <GL/glut.h>

#include "defs.h"
#include "macros.h"

#include "disp_func.h"
#include "colordefs.h"

using namespace std;

void quicksort(float *a, int *b, int low, int high)
{
  //  low IS THE LOWER INDEX, high IS THE UPPER INDEX
  //  OF THE REGION OF ARRAY a THAT IS TO BE SORTED.
  //  ARRAY b WILL GET ALSO SORTED ACCORDING TO ARRAY a.
  int i=low, j=high;
  float h;

  // COMPARISON ELEMENT x
  int  ind = int((low+high)/2); 
  float  x = a[ind];

  // PARTITION
  do
  {    
	while (a[i]<x) i++; 
	while (a[j]>x) j--;
	if (i<=j)
	{
	  h=a[i]; a[i]=a[j]; a[j]=h;
	  h=b[i]; b[i]=b[j]; b[j]=h;
	  i++; j--;
	}
  } while (i<=j);

  // USE RECURSION
  if (low<j)  quicksort(a, b, low, j);
  if (i<high) quicksort(a, b, i, high);
}
	
/*
======================================================================
		OPENGL FUNCTIONS
*/

/* 
static void printstring(void *font, char *string)
{
  int len,i;
 
  len=(int)strlen(string);
  for(i=0;i<len;i++)
    glutBitmapCharacter(font,string[i]);
}
*/

/*
======================================================================
		FRIENDS
*/

//ostream &operator<<(ostream &os, gl_obj &z)
void gl_obj::save_all(char *a, gl_obj &z)
{
		int i,j;

		ofstream os(a);

		// WRITE DATA
		//os << *z.dq << "\n";
		// WRITE GUIDATA
		//os << *z.guidq << "\n";
		// WRITE LATTICE
		//os << *z.p << "\n";

		//os << setprecision(PRECISION);

		/*
		os << z.do_measure   << " "       
				<< z.numVertColsOld << " "
				<< z.stepOnAverage << " "
				<< z.accum_ave_flag << " "
				<< z.particleFlux << " "
				<< z.particleFluxStdDev << " "
				<< z.pressureDer << " "
				<< z.pressureDerStdDev << " "
				<< z.momentumX << " "
				<< z.momentumXStdDev << " "
				<< z.momentumY << " "
				<< z.momentumYStdDev << " "
				<< z.pureBlueAveDensity << " "
				<< z.pureBlueAveDensityStdDev << " "
				<< z.pureRedAveDensity << " "
				<< z.pureRedAveDensityStdDev << " "
				<< z.pureGreenAveDensity << " "
				<< z.pureGreenAveDensityStdDev << " "
				<< z.stepOnDensityAve << " "
				<< z.perMixedSites << " "
				<< z.rgBlueBubble << " "
				<< z.prev_max_intensity_color << " "
				<< z.new_max_intensity_color << " "
				<< z.prev_max_density_color << " "
				<< z.new_max_density_color << " "
				<< z.vAng << " "
				<< z.hAng << " "
				<< z.xshift << " "
				<< z.yshift << " "
				<< z.size << " "
				<< z.window_width << " "
				<< z.window_height << " "
				<< z.num_types << " "
				<< z.whichParttype << " "
				<< z.whichPart0Interac << " "
				<< z.whichPart1Interac << " "
				<< z.whichConstForce << " "
				<< z.flag_collapse_density_y << " "
				<< "\n";
				*/

		// --------------------------------------------------
		// DUMMY PARAMETERS. CREATED FOR FUTURE COMPATIBILITY.
		// --------------------------------------------------
		/*
		os << " ";
		for (i=0; i<z.num_dummy_params; i++)
				os << z.dummy_params[i] << " ";
		os << z.num_dummy_params << " ";
		*/

		// -------------------------------
		// NO MORE WRITING AFTER THIS LINE.
		// -------------------------------

//		return os;
  		os.close();
}

//istream &operator>>(istream &os, gl_obj &z)
void gl_obj::read_all(char *a, gl_obj &z)
{
		int i,j;

		ifstream os(a);

		/*
		// READ DATA
		os >> *z.dq ;
		// READ GUIDATA
		os >> *z.guidq ;
		// READ LATTICE
		os >> *z.p ;

		os >> setprecision(PRECISION);

		os >> z.do_measure   
				>> z.numVertColsOld 
				>> z.stepOnAverage 
				>> z.accum_ave_flag 
				>> z.particleFlux 
				>> z.particleFluxStdDev 
				>> z.pressureDer 
				>> z.pressureDerStdDev 
				>> z.momentumX 
				>> z.momentumXStdDev 
				>> z.momentumY 
				>> z.momentumYStdDev 
				>> z.pureBlueAveDensity 
				>> z.pureBlueAveDensityStdDev 
				>> z.pureRedAveDensity 
				>> z.pureRedAveDensityStdDev 
				>> z.pureGreenAveDensity
				>> z.pureGreenAveDensityStdDev
				>> z.stepOnDensityAve 
				>> z.perMixedSites 
				>> z.rgBlueBubble 
				>> z.prev_max_intensity_color 
				>> z.new_max_intensity_color 
				>> z.prev_max_density_color 
				>> z.new_max_density_color 
				>> z.vAng 
				>> z.hAng 
				>> z.xshift 
				>> z.yshift 
				>> z.size 
				>> z.window_width 
				>> z.window_height 
				>> z.num_types
				>> z.whichParttype 
				>> z.whichPart0Interac
				>> z.whichPart1Interac 
				>> z.whichConstForce 
				>> z.flag_collapse_density_y
				;
		*/

		// --------------------------------------------------
		// DUMMY PARAMETERS. CREATED FOR FUTURE COMPATIBILITY.
		// --------------------------------------------------
		/*
		for (int i=0; i<z.num_dummy_params; i++)
				os >> z.dummy_params[i];
		int new_num_dummy_params;
		os >> new_num_dummy_params;
		cout << "# gl_obj.cc: DATA FILE: " << new_num_dummy_params 
				<< " DUMMY PARAMS, CURRENT: " << z.num_dummy_params
				<< "\n";
				*/

		// -------------------------------
		// NO MORE READING AFTER THIS LINE.
		// -------------------------------

//		return os;
  		os.close();
}

/*
======================================================================
		MEMBERS
*/

gl_obj::gl_obj(int x,int y,int w,int h,const char *l)
            : Fl_Gl_Window(x,y,w,h,l) 
{
  int i,j;

	window_width  = w;
	window_height = h;

//	show_credits();
}

gl_obj::~gl_obj(void)
{
}


void gl_obj::clean(void)
{
	cout << "no cleaning of lattice yet\n";
}

void gl_obj::drawCanvas() 
{
//  PRINT1("in draw cube");
    int i,j,l;

	double low_limit = -1.0;
	double hig_limit = +1.0;

    int size_x = 10;
    int size_y = 10;
	float x_left  = -size_x/2.;
	float x_right =  size_x/2.;
	float y_up    =  size_y/2.;
	float y_down  = -size_y/2.;
	float general_zoom = 1.;

    // REFERENCE LINES
    glBegin(GL_LINES);

        glColor3f(1.,0.,0.);
       glVertex3f(x_left,y_down,0.);
       glVertex3f(x_left,y_up,0.);

        glColor3f(0.,1.,0.);
       glVertex3f(x_left,y_up,0.);
       glVertex3f(x_right,y_up,0.);

        glColor3f(0.,0.,1.);
       glVertex3f(x_right,y_up,0.);
       glVertex3f(x_right,y_down,0.);

        glColor3f(1.,0.,1.);
       glVertex3f(x_right,y_down,0.);
       glVertex3f(x_left,y_down,0.);

    glEnd();
}

void gl_obj::draw() 
{
  int i,j,k;

  // ARRAYS FOR 3D SORTING AND DRAWING
  // ACCORDING TO Z POSITION.
  //   a ARRAY OF Z POSITIONS
  float *a = NULL;
  //   b ARRAY OF INDICES
  int   *b = NULL;
  if (sys.N>-1)
  {
	a = new float[sys.N];
	b = new int[sys.N];
  }

  //PRINT1("in draw");
  int size_x = 10;
  int size_y = 10;
  int size_z = 1;
  float general_zoom = 1.;

  if (!valid()) 
  {
	glLoadIdentity();
	glViewport(0,0,w(),h());
	//glOrtho(-box_size_x,box_size_x,-box_size_y,box_size_y,-20000,10000);
	//glOrtho(-box_size_x,box_size_x,-box_size_y,box_size_y,-1,1);
	//glOrtho(-.3,size_x+.3,-.3,size_y+.3,-1,1);
	glOrtho(-.3-size_x/2.,size_x/2.+.3,-.3-size_y/2.,size_y/2.+.3,-1,1);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }

  //  glClearColor(1.0,1.0,1.0,0.);
  //  glClearColor(0.9,0.9,0.9,0.);
  glClearColor(0.8,0.8,0.8,0.);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glPushMatrix();

  //glScalef(general_zoom,general_zoom,general_zoom);
  //glTranslatef(-xshift, -yshift, 0);
  //glTranslatef(-size_x/2, -size_y/2, 0);

  //    glRotatef(hAng,0,1,0); 
  //    glRotatef(vAng,1,0,0);
  //    glRotatef(3.1415,1,0,0);

  // DRAW BOUNDARY AND OBSTACLES
  drawCanvas();

  // UPDATE MEASUREMENTS
  // DONE IN SIMULATION ANIM
  //	  if ( accum_ave_flag ) 
  //		measure();

  //float scale_fac = 10.;
  float scale_fac_x = size_x/sys.Lx;
  float scale_fac_y = size_y/sys.Ly;
  float scale_fac_z = size_z/sys.Lz;
  // DRAW PARTICLES
  vec r;
  double ang_theta;

#if NDIM == 3
  // SORT DRAWING BY Z POSITION
  for (k=0; k<sys.N; k++)
  {
	a[k] = (sys.get_image_position(k)).z;
	b[k] = k;
  }
  // SORT
  if (sys.N>-1)
  	quicksort(a, b, 0, sys.N-1);
  // DRAWING BY SORTING
  for (k=0; k<sys.N; k++)
  {
	i = b[k];
#else
  // NO SORTING, JUST DRAW IN ANY ORDER
  for (i=0; i<sys.N; i++)
  {
#endif

	// GET PARTICLE POSITION AND ANGLE
	r         = sys.get_image_position(i);
	ang_theta = sys.atoms[i].theta;

	glPushMatrix();
	glTranslatef(r.x*scale_fac_x, r.y*scale_fac_y, r.z*scale_fac_z);
	// DRAW THE DISPLAY LIST
	if ( sys.dynamics !=7 && sys.dynamics !=8 )
	{
	  // NO ROTATION
		glScalef(scale_fac_x,scale_fac_y,scale_fac_z);
#if NDIM == 3
		glCallList(magenta_full_circle_DL);
#else

#  if CONFINED == 0
		glCallList(blue_circle_DL);
#  else
		if (sys.atoms[i].is_not_fixed())
			glCallList(blue_circle_DL);
		else
			glCallList(red_circle_DL);
#  endif

# if DRAWVERLETNEIGH == 1
		// DRAW VERLET LIST NEIGHBORS OF ONE PARTICLE
		if ( neighborsverlet->is_in_list(2,i) )
			glCallList(red_circle_DL);
		if ( i == 2 )
			glCallList(green_cross_circle_DL);
# endif
#endif
	}
	else if ( sys.dynamics ==7 )
	{
	  // ROTATE
	    glRotatef(DEGREE(ang_theta), 0., 0., 1.);
		glScalef(scale_fac_x,scale_fac_y,scale_fac_z);
		glCallList(red_radius_circle_DL);
	}
	else if ( sys.dynamics ==8 )
	{
	  // ROTATE
	    glRotatef(DEGREE(ang_theta), 0., 0., 1.);
		glScalef(scale_fac_x,scale_fac_y,scale_fac_z);
		glCallList(green_diam_circle_DL);
	}

	glPopMatrix();
  }

  // DRAW TAGGED PARTICLE TRAJECTORY
  if (sys.do_tag)
  {
	glPushMatrix();
	//glColor4f(WHITE);
	glColor4f(BLACK);
	glBegin(GL_LINES);
	for (i=0; i<sys.tag_particle_pos.last_one; i++)
	{
	  glVertex2f(sys.tag_particle_pos[i].x*scale_fac_x, sys.tag_particle_pos[i].y*scale_fac_y);
	  glVertex2f(sys.tag_particle_pos[i+1].x*scale_fac_x, sys.tag_particle_pos[i+1].y*scale_fac_y);
	}
	glEnd();
	glPopMatrix();
  }

  glPopMatrix();

  if (a!=NULL)
  	delete [] a;
  if (b!=NULL)
  	delete [] b;
}


void gl_obj::DoDisplayLists(void)
{
    //PRINT1("In DoDisplayLists");
    // DISPLAY LIST THAT MAKES THE FRAME OF THE SIMULATION

    blue_circle_DL     = glGenLists(1);
    red_circle_DL      = glGenLists(1);
    white_line_tag_DL  = glGenLists(1);
    red_radius_circle_DL      = glGenLists(1);
    green_diam_circle_DL      = glGenLists(1);
    green_cross_circle_DL      = glGenLists(1);
	magenta_full_circle_DL	  =  glGenLists(1);
	DisplayLists();
}

void gl_obj::DisplayLists(void)
{
   //PRINT1("In DisplayLists");

  int i,j;
  //float radius = 0.003;
  float radius = .5;
  //float radius = 1.;
  float angle;
  int max_angle = 101;

  // DISPLAY LIST FOR BLUE CIRCLE
  glNewList(blue_circle_DL, GL_COMPILE);
    glColor4f(BLUE);
    //radius = 1.0;
    //radius = 0.5;
    glBegin(GL_LINE_LOOP);
    for (i=0; i<max_angle; i++)
    {
	  angle = 2.*PI*i/(max_angle - 1);
	  glVertex2f(radius*cos(angle), radius*sin(angle));
    }
    glEnd();
  glEndList();

  // DISPLAY LIST FOR RED CIRCLE
  glNewList(red_circle_DL, GL_COMPILE);
    glColor4f(RED);
    //radius = 1.0;
    glBegin(GL_LINE_LOOP);
    for (i=0; i<max_angle; i++)
    {
	  angle = 2.*PI*i/(max_angle - 1);
	  glVertex2f(radius*cos(angle), radius*sin(angle));
    }
    glEnd();
  glEndList();

  // DISPLAY LIST FOR RED+RADIUS CIRCLE 
  glNewList(red_radius_circle_DL, GL_COMPILE);
    glColor4f(RED);
    //radius = 1.0;
    //radius = 0.5;
    glBegin(GL_LINE_LOOP);
    for (i=0; i<max_angle; i++)
    {
	  angle = 2.*PI*i/(max_angle - 1);
	  glVertex2f(radius*cos(angle), radius*sin(angle));
    }
    glEnd();
    glBegin(GL_LINES);
	  angle = PI/2.;
	  glVertex2f(radius*cos(angle), radius*sin(angle));
	  glVertex2f(0.,0.);
    glEnd();
  glEndList();

  // DISPLAY LIST FOR GREEN+DIAMETER CIRCLE 
  glNewList(green_diam_circle_DL, GL_COMPILE);
    glColor4f(MAGENTA);
    //radius = 1.0;
    //radius = 0.5;
    glBegin(GL_LINE_LOOP);
    for (i=0; i<max_angle; i++)
    {
	  angle = 2.*PI*i/(max_angle - 1);
	  glVertex2f(radius*cos(angle), radius*sin(angle));
    }
    glEnd();
    glBegin(GL_LINES);
	  angle = PI/2.;
	  glVertex2f(radius*cos(angle), radius*sin(angle));
	  glVertex2f(-radius*cos(angle), -radius*sin(angle));
    glEnd();
  glEndList();

  // DISPLAY LIST FOR GREEN+CROSSHAIR CIRCLE 
  glNewList(green_cross_circle_DL, GL_COMPILE);
    glColor4f(GREEN);
    //radius = 1.0;
    //radius = 0.5;
    glBegin(GL_LINE_LOOP);
    for (i=0; i<max_angle; i++)
    {
	  angle = 2.*PI*i/(max_angle - 1);
	  glVertex2f(radius*cos(angle), radius*sin(angle));
    }
    glEnd();
    glBegin(GL_LINES);
	  angle = PI/2.;
	  glVertex2f(radius*cos(angle), radius*sin(angle));
	  glVertex2f(-radius*cos(angle), -radius*sin(angle));
    glEnd();
    glBegin(GL_LINES);
	  angle = 0.;
	  glVertex2f(radius*cos(angle), radius*sin(angle));
	  glVertex2f(-radius*cos(angle), -radius*sin(angle));
    glEnd();
  glEndList();

  // DISPLAY LIST FOR MAGENTA FILLED CIRCLE
  glNewList(magenta_full_circle_DL, GL_COMPILE);
    glColor4f(MAGENTA);
	glLineWidth(1.);
    //radius = 1.0;
    //radius = 0.5;
	int   num_inner_circ  = 8;
	//int   num_inner_circ  = 10;
	float incr_inner_circ = 1./num_inner_circ;
	float inner_circ = 0.;
	float radius_incr = radius / num_inner_circ; 
	for (j=0; j<num_inner_circ; j++)
	{
	  glBegin(GL_LINE_LOOP);
	  for (i=0; i<max_angle; i++)
	  {
		angle = 2.*PI*i/(max_angle - 1);
		glVertex2f(radius*cos(angle), radius*sin(angle));
	  }
	  glEnd();
	  radius     -= radius_incr;
	  inner_circ += incr_inner_circ;
	  glLineWidth(3.);
      glColor4f(1.0, inner_circ, 1.0, 1.0);
	}
	glLineWidth(1.);
  glEndList();
}

int gl_obj::handle(int event) 
{
  /*
    int box_size_x = 100;
    int box_size_y = 100;
    int size_x = 100;
    int size_y = 100;
	float general_zoom = 1.;

	   int i,j;
	   int x_ext = size_x;
	   int y_ext = size_y;
	   double x_trans = xshift+size_x/2.0;
	   double y_trans = yshift+size_y/2.0;
	   double mousex, mousey;
	   double scale = general_zoom;
	   char outputHandle[256];

       // MOUSE POSITION IN double COORDINATES
       guidq->mouse_x = (double(Fl::event_x())/double(window_width)) * 2.0 * guidq->box_size_x / scale - guidq->box_size_x/scale + xshift + guidq->size_x/2.0;
       guidq->mouse_y = (double(window_height-Fl::event_y())/double(window_height)) * 2.0 * guidq->box_size_y / scale - guidq->box_size_y/scale + yshift + guidq->size_y/2.0;

       switch(event) 
       {
        case FL_PUSH:
	    {
         // DETECT DOUBLE CLICK
         if (Fl::event_is_click())
         {
          // Fl::event_button() = 1     left mouse button
          // Fl::event_button() = 2     middle mouse button
          // Fl::event_button() = 3     right mouse button
          if (Fl::event_button() == 1)
          {
            // DETERMINE CLOSEST RULER TO SELECT
            double distanceToRuler = double(guidq->size_x);
            int   selected = -1;
            for (i=0; i<numVertCols; i++)
            {
                    if ( ABS(guidq->mouse_x-double(vertColsX[i])) < distanceToRuler )
                    {
                        distanceToRuler = ABS(guidq->mouse_x-double(vertColsX[i]));
                        selected = i;
                    }
            }
            // MARK THE SELECTED RULER
            if (selected != -1)
            {
                for (i=0; i<numVertCols; i++)
                {
                    vertColsSelected[i] = 0;
                }
                vertColsSelected[selected] = 1;
            }
            redraw();
          }
         }

         // FOR POPUP DIALOG WITH RIGHT MOUSE BUTTON
         if (Fl::event_button() == 3)
         {
            if ( (guidq->mouse_x>0) && (guidq->mouse_x<guidq->size_x)
              && (guidq->mouse_y>0) && (guidq->mouse_y<guidq->size_y)
               )
            {
		      // POSITION THE POPUP FOR SAVING TO FILE
   		      int popupx = mainWin->x() + winPopupSaveVars->w()/2 + Fl::event_x();
              int popupy = mainWin->y() + winPopupSaveVars->h()/2 + Fl::event_y();
		      winPopupSaveVars->position(popupx, popupy);
		      // SHOW THE POPUP
//			  winPopupSaveVars->show();
//			  sprintf(outputHandle, "x: %d y:%d",int(guidq->mouse_x),int(guidq->mouse_y));
		      PRINT5("( x:",int(guidq->mouse_x),") (y:",int(guidq->mouse_y)," )");
		      mouseOutpos->label(outputHandle);
            }
         }
         return 1;       	
	    }
        case FL_DRAG:
        {
          // DRAGGING WITH LEFT MOUSE BUTTON
          if (Fl::event_button() == 1)
          {
            // MOVE RULER THAT IS SELECTED ...
            int selected = oneVertColIsSelected();
            if (selected != -1)
            {
                // ... AND IS CLOSE TO THE POINTER
                int mouseRange = 500;
                if ( ABS(guidq->mouse_x-double(vertColsX[selected])) < mouseRange )
                {
                   if ( (guidq->mouse_x>0) && (guidq->mouse_x<guidq->size_x) )
                     vertColsX[selected] = int( guidq->mouse_x );
                }
            }
            redraw();
          }
          return 1;
        }
       	case FL_RELEASE:    
        {
         	//... MOUSE UP EVENT ...
            int i;
            for (i=0; i<numVertCols; i++)
            {
                    vertColsSelected[i] = 0;
            }
            redraw();
         	return 1;
        }
	default:
         	return 0;
       }
	   */
/*
 * THIS IS THE COMPLETE COLLECTION OF EVENTS TO HANDLE
       switch(event) 
       {
        case FL_PUSH:
		PRINT3("mouse at",Fl::event_x(),Fl::event_y());
         	//... mouse down event ...
         	//... position in Fl::event_x() and Fl::event_y()
         	return 1;
		case FL_KEYBOARD:
	    {
         	//... keypress, 
			//	key is in Fl::event_key(), ascii in Fl::event_text()
         	//... Return 1 if you understand/use the keyboard event, 
			//	0 otherwise...
			PRINT3("keyboard",Fl::event_key(),Fl::event_text());
         	return 1;
	    }
       	case FL_DRAG:
		PRINT1("dragging");
         	//... mouse moved while down event ...
         	return 1;
       	case FL_RELEASE:    
		PRINT1("release");
         	//... mouse up event ...
         	return 1;
       	case FL_FOCUS :
		PRINT1("focus ");
         	//... Return 1 if you want keyboard events, 0 otherwise
         	return 1;
       	case FL_UNFOCUS :
		PRINT1("unfocus ");
         	//... Return 1 if you want keyboard events, 0 otherwise
         	return 1;
       	case FL_SHORTCUT:
		PRINT1("shortcut");
         	//... shortcut, 
		//	key is in Fl::event_key(), ascii in Fl::event_text()
         	//... Return 1 if you understand/use the shortcut event, 
		//	0 otherwise...
         	return 1;
	case FL_NO_EVENT:
		PRINT1(" FL_NO_EVENT");
         	return 1;
	case FL_ENTER:
		PRINT1("FL_ENTER");
         	return 1;
	case FL_LEAVE:
		PRINT1("FL_LEAVE");
         	return 1;
	case FL_CLOSE:
		PRINT1("FL_CLOSE");
         	return 1;
	case FL_MOVE:
		PRINT1("FL_MOVE");
         	return 1;
	case FL_DEACTIVATE:
		PRINT1("FL_DEACTIVATE");
         	return 1;
	case FL_ACTIVATE:
		PRINT1("FL_ACTIVATE");
         	return 1;
	case FL_HIDE:
		PRINT1("FL_HIDE");
         	return 1;
	case FL_SHOW:
		PRINT1("FL_SHOW");
         	return 1;
	case FL_PASTE:
		PRINT1("FL_PASTE");
         	return 1;
	case FL_SELECTIONCLEAR:
		PRINT1("FL_SELECTIONCLEAR");
         	return 1;
       	default:
		cout << "event not found " << event << "\n";
         	return 0;
       }
*/
}

void gl_obj::show_credits(void)
{
  cout << "#------------------------------------------\n";
  cout << "#-- Many Body Particle Simulation Program --\n";
  cout << "#-- Version 0.06                          --\n";
  cout << "#-- Using the FLTK toolkit  vrs. 1.1.9    --\n";
  cout << "#-- Luis Cruz                             --\n";
  cout << "#-- ccruz@drexel.edu                      --\n";
  cout << "#-- Drexel University                     --\n";
  cout << "#-- 2009-2010                             --\n";
  cout << "#------------------------------------------\n";
  cout << "#                                          \n";
}

/* ------------------------------------------------- */

