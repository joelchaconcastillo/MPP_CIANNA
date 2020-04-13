#ifndef __MPP_H__
#define __MPP_H__

#include <bits/stdc++.h>
using namespace std;
#define FOREACH(i, v) for (__typeof((v).begin()) i = (v).begin(); i != (v).end(); i++)
//v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner, v_both_snack;
#define EPSILON 1e-10

#define N_CATEGORIES 3
#define CATEGORY_1 1
#define CATEGORY_2 2
#define CATEGORY_BOTH 0
#define GLOBAL 1
#define DIARIA 2

//encoded times by each day of the individual
#define N_TIMES 6
#define N_OPT_DAY 8
#define BREAKFAST 0
#define MORNING_SNACK 1
#define STARTER_1 2
#define STARTER_2 3
#define MAIN_COURSE_1 4
#define MAIN_COURSE_2 5
#define EVENING_SNACK 6
#define DINNER 7



//#define BOTH_SNACK 8
////////crossover type......
#define PAIR_BASED_CROSSOVER 1
#define UNIFORM_CROSSOVER 2
#define UNIFORM2_CROSSOVER 3
#define WEIGHT_DAY 1.0e6
#define DAYS_FAVORITE 7*3
#define DAYS_NO_FAVORITE 7*4
#define ITERATIONS_LS 100000

#define W_VAR_DAY 1000
#define W_VAR_GLOBAL 100
#define W_VAR_GLOBAL_CAT 1
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
struct Solution_LS
{
   set<int> badDays;
   int heaviestNut, heaviestType;
   vector<double> obj_values;
   vector<int> x_var;
   vector< vector<double> > globalPlan;
   vector< vector< vector<double> > > nutriment_per_day;
   vector< vector< vector < int > > > time_id_day_table;
   vector< vector< int > > time_diff;
   vector<vector<int>> uniq_per_day;
};
class MPP_Problem{
	public:
		MPP_Problem();
		~MPP_Problem(){
		}

 		void load_data(int argc, char **argv);
		void load_constraints(char *Plates_file);
		void load_dishes(char *Constraints_file);
		inline int random_dish(int time_dish){return rand()%((int)v_opt_dishes[time_dish].size());}

		vector<vector<infoDishes> > v_opt_dishes;  // the same as ->  vector<int> v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner;
		vector<constraint_nutrient> v_constraints;
		unordered_map<string, int> dic_nut_id;
		vector<int> v_constraint_global, v_constraint_day;
		string out_filename;
		vector<vector<int>> conf_day, opt_conf, unique_opt_time;
		vector<int> inv_unique_opt_time;
		vector<double> weights;
	 	vector<int> priority_time;
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
		void evaluate(vector<int> &sol, vector<double> &objs);
		void restart();
		//Individual *clone() const;
		void init (); //initialization 
		void dependentMutation(double pm);
		void dependentCrossover(MPP &i2);
		void uniformCrossover(MPP &i2);
		void uniform2Crossover(MPP &i2);
		void pairBasedCrossover(MPP &i2);
		void localSearch();
		void First_Improvement_Hill_Climbing(vector<Neighbor> &neighbors, vector<int> &current_sol, vector<double> &objs);
		void First_Improvement_Hill_Climbing_swap(vector<Neighbor_swap> &neighbors, vector<int> &best_sol, vector<double> &best_objs);
		int getDistance(MPP &ind2); 
		void exportcsv();
		virtual void print(ostream &os) const;
//!menu planning

		vector<int> x_var;
		double fitness;
		static MPP_Problem *MPP_problem;

	private:
		void calculateFeasibilityDegree();
		void calculateFeasibilityDegree(vector<int> &sol, double &feas);

		void init_incremental_evaluation(struct Solution_LS &current);
		void init_incremental_evaluation2(struct Solution_LS &current);
		void inc_eval(struct Solution_LS &current, Neighbor &new_neighbor, vector<double> &new_objs);
		void update_inc(struct Solution_LS &current, Neighbor &new_neighbor, vector<double> &new_objs);
		void swap_days(vector<int> &data, int day1, int day2);
  	        inline void perturb_opt(vector<int> &data, int day, int which)
		{
	           set<int> id_Day;
		   for(int opt = 0; opt < N_OPT_DAY; opt++)
		     id_Day.insert(MPP_problem->v_opt_dishes[opt][data[day*N_OPT_DAY+opt]].description);
		   int rd = MPP_problem->random_dish(which);
	           int id = MPP_problem->v_opt_dishes[which][rd].description;
		   while(id_Day.find(id) != id_Day.end())//it forces to different dishes in a day
		   { 
		      rd = MPP_problem->random_dish(which);
		      id = MPP_problem->v_opt_dishes[which][rd].description;
	           }
		   data[day * N_OPT_DAY + which] = rd;
		}
		inline void perturb_day(vector<int> &data, int day){ 		
			set<int> dish;
			for(int k = 0; k < N_OPT_DAY; k++) 
		 	{
			   int rd = MPP_problem->random_dish(k);
			   int id = MPP_problem->v_opt_dishes[k][rd].description;
			   while(dish.find(id) != dish.end())
			   { 
			      rd = MPP_problem->random_dish(k);
			      id = MPP_problem->v_opt_dishes[k][rd].description;
			   }
			   dish.insert(id);
			   data[day*N_OPT_DAY + k] = rd;
			
			}
		}
 		inline double f(pair<int, int> data_dcn){ return ((double)data_dcn.first + ( 1.0 - ((double)data_dcn.second/(double)nDias)));}

 		inline void update_dcn_pair(int diff, pair<int, int> &p_dcn){ if(diff < p_dcn.first)p_dcn = make_pair(diff, 1);else if(diff == p_dcn.first)p_dcn.second++; }

		void calculateVariability(vector<int> &sol, vector<double> &objs);
		int heaviestNut, heaviestType;
		double valorFac, variabilidadObj;//factibility and variability of the current solution..
	        vector<double> obj_values;//{feasiblity, varibility x times}
		set<int> badDays;
		bool comp_objs(vector<double> &variability_v1, vector<double> &variability_v2);



		void  First_Improvement_Hill_Climbing2(vector<Neighbor> &neighbors, vector<int> &best_sol, vector<double> &best_objs);

};
#endif
