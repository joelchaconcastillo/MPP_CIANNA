#ifndef __MPP_H__
#define __MPP_H__

#include <bits/stdc++.h>
using namespace std;
#define FOREACH(i, v) for (__typeof((v).begin()) i = (v).begin(); i != (v).end(); i++)
//v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner, v_both_snack;

#define CATEGORY_1 1
#define CATEGORY_2 2
#define CATEGORY_BOTH 0

//encoded times by each day of the individual
#define N_OPT_DAY 8
#define BREAKFAST 0
#define MORNING_SNACK 1
#define STARTER_1 2
#define STARTER_2 3
#define MAIN_COURSE_1 4
#define MAIN_COURSE_2 5
#define EVENING_SNACK 6
#define DINNER 7

//#define N_FOODS_DAY 6
////Options of the v_times_dishes..
//#define BREAKFAST 0 
//#define MORNING_SNACK 1
//#define STARTER 2
//#define MAIN_COURSE 3
//#define EVENING_SNACK 4
//#define DINNER 5

//#define BOTH_SNACK 8
////////crossover type......
#define PAIR_BASED_CROSSOVER 1
#define UNIFORM_CROSSOVER 2
#define UNIFORM2_CROSSOVER 3

//extern volatile bool finished;
extern int crossoverType;
extern int nDias;
void printBest();
struct infoDishes {
        int description;	
	string time_day; //time related with the dish
	vector<double> v_nutrient_value;    //nutriments meta-data...
	int category; //category, at this point is 1 or 2
	bool favorite; //true if this is a favorite dish...
};
struct constraint_nutrient
{
   double min, max;
   string type, name;
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

//		vector<infoDishes> v_dishes;
		vector<vector<infoDishes> > v_times_dishes;  // the same as ->  vector<int> v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner;
		vector<constraint_nutrient> v_constraints;
		unordered_map<string, int> dic_nut_id;
		vector<int> v_constraint_global, v_constraint_day;
		string out_filename;
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
		int getDistance(MPP &ind2); 
		double inline getMaximum(const int i) const { cerr << "ErrorL llama a getMaximum" << endl; exit(-1); return 0; }
		double inline getMinimum(const int i) const { cerr << "Error: llama a getMinimum" << endl; exit(-1); return 0; }
		//unsigned int inline getOptDirection(const int i) const { return MAXIMIZE; }
		void exportcsv();
		virtual void print(ostream &os) const;

//!menu planning

		vector<int> x_var;
		double fitness;
		static MPP_Problem *MPP_problem;

	private:
		void calculateFeasibilityDegree();
		int heaviestNut, heaviestType;
		double valorFac, variabilidadObj;//factibility and variability of the current solution..
		set<int> badDays;

};

#endif
