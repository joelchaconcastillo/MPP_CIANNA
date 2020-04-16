main: MA.o MPP.o utils.o main.cpp 
	g++ -std=c++11 -O2 -o $@ main.cpp MA.o MPP.o utils.o

MPP.o: MPP.h MPP.cpp
	g++ -std=c++11 -O2 -c -o $@ MPP.cpp

MA.o: MA.h MPP.h MA.cpp
	g++ -std=c++11 -O2 -c -o $@ MA.cpp

utils.o: utils.cpp utils.h
	g++ -O2 -c -o $@ utils.cpp

clean:
	rm -f main MPP.o MA.o utils.o
