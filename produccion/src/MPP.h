#ifndef __MPP_H__
#define __MPP_H__

#include <bits/stdc++.h>
using namespace std;
#define FOREACH(i, v) for (__typeof((v).begin()) i = (v).begin(); i != (v).end(); i++)
//v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner, v_both_snack;
#define EPSILON 1e-50

#define CATEGORY_1 1
#define CATEGORY_2 2
#define CATEGORY_BOTH 0
#define GLOBAL 1
#define DIARIA 2

//encoded times by each day of the individual
#define N_OPT_DAY 8
#define BREAKFAST 0
#define MORNING_SNACK 1
#define STARTER_1 2
#define MAIN_COURSE_1 3
#define EVENING_SNACK 4
#define DINNER 5
#define STARTER_2 6
#define MAIN_COURSE_2 7


//#define BOTH_SNACK 8
////////crossover type......
#define PAIR_BASED_CROSSOVER 1
#define UNIFORM_CROSSOVER 2
#define UNIFORM2_CROSSOVER 3
#define WEIGHT_DAY 1.0e6
#define DAYS_FAVORITE 7*3
#define DAYS_NO_FAVORITE 7*4
#define ITERATIONS_LS 100000
//extern volatile bool finished;
extern int crossoverType;
extern int nDias;
void printBest();
struct Neighbor {
	int variable;
	int newValue;
};

struct Neighbor_swap{
	int day1, day2;
};


struct infoDishes {
        int description;	
	string time_day; //time related with the dish
	vector<double> v_nutrient_value;    //nutriments meta-data...those indexes should be in the same order than the v_contraints vector...
	int category; //category, at this point is 1 or 2
	bool favorite; //true if this is a favorite dish...
};
struct constraint_nutrient
{
   double min, max;
   string name;
   int type;
};
class MPP_Problem{
	public:
		MPP_Problem();
		~MPP_Problem(){
		}

 		void load_data(int argc, char **argv);
		void load_constraints(char *Plates_file);
		void load_dishes(char *Constraints_file);
		inline int random_dish(int time_dish){return rand()%((int)v_times_dishes[time_dish].size());}

		vector<vector<infoDishes> > v_times_dishes;  // the same as ->  vector<int> v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner;
		vector<constraint_nutrient> v_constraints;
		unordered_map<string, int> dic_nut_id;
		vector<int> v_constraint_global, v_constraint_day;
		string out_filename;
		vector<vector<int>> conf_day;
	        vector<vector<int>> time_conf;
		int max_description_id;
};
class MPP{
	public:
		MPP(){
		}
		~MPP(){
		}

//////menu planning
		void evaluate();
		void restart();
		//Individual *clone() const;
		void init (); //initialization 
		void dependentMutation(double pm);
		void dependentCrossover(MPP &i2);
		void uniformCrossover(MPP &i2);
		void uniform2Crossover(MPP &i2);
		void pairBasedCrossover(MPP &i2);
		void localSearch();
		pair<double, double> First_Improvement_Hill_Climbing(vector<Neighbor> &neighbors, vector<int> &current_sol);
		pair<double, double> First_Improvement_Hill_Climbing_swap(vector<Neighbor_swap> &neighbors, pair<double, double> &currentResult, vector<int> &x_var);
		int getDistance(MPP &ind2); 
		void full_search();
		void exportcsv();
		virtual void print(ostream &os) const;


		void calculateFeasibilityDegree2();
//!menu planning

		vector<int> x_var;
		double fitness;
		static MPP_Problem *MPP_problem;

	private:
		void calculateFeasibilityDegree();
		void my_next_permutation(vector<int> &perm, vector<int> &v_max_opt);

		double init_incremental_evaluation(vector< vector<double> > &globalPlan, vector< vector< vector<double> > >&nutriment_per_day, vector<int> &sol);
		double inc_eval_feas_time(vector< vector<double> > &globalPlan, vector< vector<vector<double> > > &nutriment_per_day, vector<int> &current_sol, Neighbor &new_neighbor, double current_infeasibility);

		void update_data_incremental_eval(vector< vector<double> > &globalPlan, vector< vector<vector<double> > > &nutriment_per_day, vector<int> &current_sol, Neighbor &new_neighbor);

		void swap_days(vector<int> &data, int day1, int day2);
		inline void perturb_day(vector<int> &data, int day){ for(int k = 0; k < N_OPT_DAY; k++) data[day*N_OPT_DAY + k] = MPP_problem->random_dish(k);}

		double calculateVariability(vector<int> &current_sol);

		double calculateVariability();
		bool day_constraint(infoDishes &dish1, infoDishes &dish2);
		int heaviestNut, heaviestType;
		double valorFac, variabilidadObj;//factibility and variability of the current solution..
		set<int> badDays;
		vector<bool> used_IDs, used_IDs_day;
};

#endif
