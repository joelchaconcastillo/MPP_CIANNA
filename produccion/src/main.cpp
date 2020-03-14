#include "MA.h"
/*
   Input arguments:
   1) file of Dishes information
   2) file of Constraintment nutrients
   3) Days
   4) Time
   5) out file
*/
//int crossoverType;
int main(int argc, char **argv){

	srand(1); //check this!!! ??
//	srand(time(NULL)); //check this!!! ??
	int N = 10;
	double pc = 1.0;
	double pm = 0.01;
	double finalTime = atof(argv[4])*60;//25 * 60;
	crossoverType= PAIR_BASED_CROSSOVER;//UNIFORM_CROSSOVER; //UNIFORM2_CROSSOVER, PAIR_BASED_CROSSOVER

	//loading the input data...
	MPP_Problem STP;
        STP.load_data(argc, argv);
	MPP::MPP_problem = &STP;
        //everithing seem's to be OK, thus starting the algorithm...
	MA ma(N, pc, pm, finalTime);
 	
	ma.run();

//	ExtendedIndividual *ei = new ExtendedIndividual();
//	ei->ind.init();
//	ei->ind.full_search();

}
