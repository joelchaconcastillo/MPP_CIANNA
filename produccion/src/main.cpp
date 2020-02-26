#include "MA.h"
/*
   Input arguments:
   1) file of Dishes information
   2) file of Constraintment nutrients
   3) Days
   4) Time
   5) out file
*/
int main(int argc, char **argv){

	srand(time(NULL)); //check this!!! ??
	int N = 50;
	double pc = 0.9;
	double pm = 0.01;
	double finalTime = atof(argv[4])*60;//25 * 60;

	//loading the input data...
	MPP_Problem STP;
        STP.load_data(argc, argv);

        //everithinh seem's to be OK, thus starting the algorithm...
	MA ma(N, pc, pm, finalTime);
	MPP::MPP_problem = &STP;
	ma.run();
}
