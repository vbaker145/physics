INCLUDES = -I../../ -I../../base
OBJECTS = ../../objects
CPP		= g++
OPT		= 
DEFINES = 

all:	main.o
	g++ main.o ${OBJECTS}/ran01.o -o gauss_rand -lm

clean:
	rm -f *.o *~ out* *rsyncBak*

main.o:	main.cc ../../base/ran01.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@

