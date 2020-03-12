main: MA.o MPP.o utils.o main.cpp 
	g++ -std=c++11 -O2 -o $@ main.cpp MA.o MPP.o utils.o -fopenmp

MPP.o: MPP.h MPP.cpp
	g++ -std=c++11 -O2 -c -o $@ MPP.cpp -fopenmp

MA.o: MA.h MPP.h MA.cpp
	g++ -std=c++11 -O2 -c -o $@ MA.cpp -fopenmp

utils.o: utils.cpp utils.h
	g++ -O2 -c -o $@ utils.cpp

clean:
	rm -f main MPP.o MA.o utils.o
