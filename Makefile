#----------------------------------------------------------
CXX=g++
CFLAGS= -I. 
LFLAGS= 

default: 

clean:
	rm imgwatermarker imgwatermark_checker

imgwatermarker: imgwatermarker.cpp
	g++ -o imgwatermarker imgwatermarker.cpp  -Wall -O -I/usr/include/SDL2 -lSDL2_image -lSDL2 -lfftw3

imgwatermark_checker: imgwatermarkchecker.cpp
	g++ -o imgwatermark_checker imgwatermarkchecker.cpp  -Wall -O -I/usr/include/SDL2 -lSDL2_image -lSDL2 -lfftw3


%.o:    %.C
	${CXX} -c ${CFLAGS} $<
