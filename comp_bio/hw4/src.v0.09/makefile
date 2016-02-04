
OBJECTS		=${OBJDIR}/vec.o \
			${OBJDIR}/atom.o ${OBJDIR}/ensemble.o ${OBJDIR}/ran01.o \
			${OBJDIR}/forces.o ${OBJDIR}/integrator.o ${OBJDIR}/minima.o  \
			${OBJDIR}/diffusion.o ${OBJDIR}/paircorrfunc.o ${OBJDIR}/verletlist.o  \
			${OBJDIR}/stepwise.o ${OBJDIR}/langevin.o ${OBJDIR}/montecarlo.o \
			${OBJDIR}/wait.o ${OBJDIR}/R.o

OBJ_GRAPHICS = ${OBJDIR}/gl_obj.o ${OBJDIR}/liquids.o  ${OBJDIR}/disp_func.o

OBJ_MAIN_MD =  ${OBJDIR}/main.o 


OBJDIR 		= objects
BASE_DIR    = ${PWD}
CURR_DIR    = `basename ${PWD}`

DEFINES		= 
#OPT 		=
OPT 		= -O3
LIBS_GRAPHICS = -L/usr/lib/fltk-1.1 -lfltk -lfltk_gl -lGLU -lGL -L/usr/X11R6/lib -lXext -lX11 -lm
LIBS 		= -lm ${LIBS_GRAPHICS}
INCLUDES 	= -I. -Ibase -Igraphics
CPP			= g++
EXEC		= run_me_vince

# --------------------------------------------------


first_make:
	make -s -f makefile usage

usage:
	echo ""
	echo "TARGETS FOR THE LIQUIDS PROGRAM:"
	echo ""
	echo "   MAIN PROGRAM:"
	echo "          1. make md"
	echo "          2. make md_graphics"
	echo ""

1: do_md
2: do_md_graphics

# ------ TARGETS -----------------------------------

do_md:
	make -f makefile all OBJ_GRAPHICS='' LIBS_GRAPHICS=''
do_md_graphics:
	make -f makefile all OBJ_MAIN_MD='' EXEC='run_vis'

# --------------------------------------------------

all: ${OBJDIR} ${OBJECTS} ${OBJ_GRAPHICS} ${OBJ_MAIN_MD}
	${CPP} ${OPT} ${OBJECTS} ${OBJ_GRAPHICS} ${OBJ_MAIN_MD} -o ${EXEC} ${LIBS}

clean:
	rm -f ${OBJDIR}/*.o* ${EXEC}
	rm -f *~ base/*~ graphics/*~
	rm -f *rsyncBak*  base/*rsyncBak*  graphics/*rsyncBak* run_vis
resclean:
	cd results; \
	rm -f gr.dat r2t.dat thermo.dat vacf.dat params.txt tag.dat totlinmom.dat warnings.txt *rsyncBak*  *~; \
	cd ..

realclean: clean resclean
	rm -f out* *.dat

${OBJDIR}:
	mkdir ${OBJDIR}

tar:
	cd .. ; \
	tar cvzf CompBio_${CURR_DIR}.tgz ${CURR_DIR}

# --------------------------------------------------

${OBJDIR}/main.o:	main.cc atom.h ensemble.h integrator.h forces.h minima.h defs.h base/macros.h diffusion.h paircorrfunc.h R.h verletlist.h stepwise.h langevin.h montecarlo.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@

${OBJDIR}/ran01.o:	base/ran01.cc base/ran01.h  base/macros.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} base/${@F:.o=.cc} -o $@
${OBJDIR}/vec.o:	base/vec.cc base/vec.h base/macros.h defs.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} base/${@F:.o=.cc} -o $@
${OBJDIR}/wait.o:	base/wait.cc base/wait.h base/macros.h defs.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} base/${@F:.o=.cc} -o $@
${OBJDIR}/atom.o:	atom.cc atom.h defs.h base/macros.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/ensemble.o:	ensemble.cc ensemble.h defs.h base/macros.h base/ran01.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/forces.o:	forces.cc forces.h defs.h base/macros.h ensemble.h verletlist.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/integrator.o:	integrator.cc integrator.h defs.h base/macros.h  forces.h ensemble.h forces.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/minima.o:	minima.cc minima.h defs.h base/macros.h forces.h ensemble.h base/wait.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/diffusion.o:	diffusion.cc diffusion.h defs.h base/macros.h ensemble.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/paircorrfunc.o:	paircorrfunc.cc paircorrfunc.h defs.h base/macros.h ensemble.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/R.o:	R.cc R.h defs.h base/macros.h ensemble.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/verletlist.o:	verletlist.cc verletlist.h defs.h base/macros.h forces.h ensemble.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/stepwise.o:	stepwise.cc stepwise.h defs.h base/macros.h ensemble.h verletlist.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@

${OBJDIR}/gl_obj.o: graphics/gl_obj.cc graphics/gl_obj.h defs.h base/macros.h atom.h ensemble.h integrator.h forces.h minima.h defs.h diffusion.h paircorrfunc.h R.h verletlist.h stepwise.h langevin.h montecarlo.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} graphics/${@F:.o=.cc} -o $@
${OBJDIR}/disp_func.o: graphics/disp_func.cc graphics/disp_func.h defs.h base/macros.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} graphics/${@F:.o=.cc} -o $@
${OBJDIR}/liquids.o: graphics/liquids.h graphics/gl_obj.h graphics/liquids.cxx defs.h base/macros.h atom.h ensemble.h integrator.h forces.h minima.h defs.h diffusion.h paircorrfunc.h R.h verletlist.h stepwise.h langevin.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} graphics/${@F:.o=.cxx} -o $@
${OBJDIR}/langevin.o:	langevin.cc langevin.h defs.h base/macros.h  forces.h ensemble.h forces.h base/ran01.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/montecarlo.o:	montecarlo.cc montecarlo.h defs.h base/macros.h  forces.h ensemble.h forces.h base/ran01.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@


# -
