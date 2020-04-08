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
	srand(time(NULL)); //check this!!! ??
	int N = 10;
	double pc = 1.0;
	double pm = 0.01;
	double finalTime = atof(argv[4])*60;//25 * 60;
	vector<vector<int>> conf_day, opt_conf(N_OPT_DAY), unique_opt_time(N_TIMES);
        vector<int> inv_unique_opt_time(N_OPT_DAY);
	conf_day.push_back({BREAKFAST, MORNING_SNACK, STARTER_1, MAIN_COURSE_1, EVENING_SNACK, DINNER});
	conf_day.push_back({BREAKFAST, MORNING_SNACK, STARTER_2, MAIN_COURSE_2, EVENING_SNACK, DINNER});
	for(int i = 0; i < conf_day.size(); i++)
	{
	   for(int j = 0; j  < conf_day[i].size(); j++)
	   {
	     opt_conf[conf_day[i][j]].push_back(i);//it maps a day configuration with the encoded individual day..
	     if( find(unique_opt_time[j].begin(), unique_opt_time[j].end(), conf_day[i][j]) == unique_opt_time[j].end() )
 	     unique_opt_time[j].push_back(conf_day[i][j]);
	     inv_unique_opt_time[conf_day[i][j]] = j;
	   }
	}

	crossoverType= PAIR_BASED_CROSSOVER;//UNIFORM_CROSSOVER; //UNIFORM2_CROSSOVER, PAIR_BASED_CROSSOVER

	//loading the input data...
	MPP_Problem STP;
        STP.load_data(argc, argv);
	STP.conf_day = conf_day;
	STP.opt_conf = opt_conf;
	STP.unique_opt_time = unique_opt_time;
	STP.inv_unique_opt_time = inv_unique_opt_time;
//	STP.weights =  {0.15, 0.3, 0.05, 0.05, 0.3, 0.15};
	STP.weights =  {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	STP.priority_time=  {inv_unique_opt_time[MAIN_COURSE_1], inv_unique_opt_time[STARTER_1], inv_unique_opt_time[BREAKFAST], inv_unique_opt_time[DINNER], inv_unique_opt_time[MORNING_SNACK], inv_unique_opt_time[EVENING_SNACK]};


	MPP::MPP_problem = &STP;


        //everithing seem's to be OK, thus starting the algorithm...
	MA ma(N, pc, pm, finalTime);
 	
	ma.run();

//	ExtendedIndividual *ei = new ExtendedIndividual();
//	ei->ind.init();
//	ei->ind.full_search();

}
