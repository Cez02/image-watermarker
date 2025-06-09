#----------------------------------------------------------
CXX=g++
CXXFLAGS = -I. -std=c++20
CFLAGS= -I. -std=c++20
LFLAGS= 

default: iwm_creator iwm_checker

clean:
	rm iwm_creator iwm_checker

iwm_creator: iwm_creator.cpp
	g++ ${CFLAGS} -o iwm_creator iwm_creator.cpp  -Wall -O -I/usr/include/SDL2 -lSDL2_image -lSDL2 -lfftw3

iwm_checker: iwm_checker.cpp
	g++ ${CFLAGS} -o iwm_checker iwm_checker.cpp  -Wall -O -I/usr/include/SDL2 -lSDL2_image -lSDL2 -lfftw3


%.o:    %.C
	${CXX} -c ${CFLAGS} $<
