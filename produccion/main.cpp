#include "MA.h"

int main(int argc, char **argv){

	srand(time(NULL)); //check this!!! ??
	int N = 50;
	double pc = 0.9;
	double pm = 0.01;
	double finalTime = atof(argv[4])*60;//25 * 60;
	

	//getting the input data...
	MPP_Problem STP;
        STP.load_data(argc, argv);

        //everithinh seem's to be OK, thus starting the algorithm...
	MA ma(N, pc, pm, finalTime);
	//string file = string(argv[1]);
//	MPP::MPP_problem = &STP;
	ma.run(STP);
}
