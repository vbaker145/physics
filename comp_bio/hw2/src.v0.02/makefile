
OBJECTS		= ${OBJDIR}/main.o ${OBJDIR}/vec.o \
			${OBJDIR}/atom.o ${OBJDIR}/ensemble.o ${OBJDIR}/ran01.o \
			${OBJDIR}/forces.o ${OBJDIR}/integrator.o ${OBJDIR}/minima.o 

OBJDIR 		= objects
BASE_DIR    = ${PWD}
CURR_DIR    = `basename ${PWD}`

DEFINES		= 
OPT 		= -O3
LIBS 		= -lm
INCLUDES 	= -I. -Ibase
CPP			= g++
EXEC		= run_me

# --------------------------------------------------

all: ${OBJDIR} ${OBJECTS}
	${CPP} ${OPT} ${OBJECTS} -o ${EXEC} ${LIBS}

clean:
	rm -f ${OBJDIR}/*.o* ${EXEC}
	rm -f *~ base/*~
	rm -f *rsyncBak*  base/*rsyncBak* 
clean_out:
	rm -f out*

${OBJDIR}:
	mkdir ${OBJDIR}

tar:
	cd .. ; \
	tar cvzf CompBio_${CURR_DIR}.tgz ${CURR_DIR}

# --------------------------------------------------

${OBJDIR}/main.o:	main.cc atom.h ensemble.h integrator.h forces.h minima.h defs.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@

${OBJDIR}/ran01.o:	base/ran01.cc base/ran01.h  base/macros.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} base/${@F:.o=.cc} -o $@
${OBJDIR}/vec.o:	base/vec.cc base/vec.h base/macros.h defs.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} base/${@F:.o=.cc} -o $@
${OBJDIR}/atom.o:	atom.cc atom.h defs.h base/macros.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/ensemble.o:	ensemble.cc ensemble.h defs.h base/macros.h base/ran01.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/forces.o:	forces.cc forces.h defs.h base/macros.h ensemble.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/integrator.o:	integrator.cc integrator.h defs.h base/macros.h  forces.h ensemble.h forces.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/minima.o:	minima.cc minima.h defs.h base/macros.h forces.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@

# -
