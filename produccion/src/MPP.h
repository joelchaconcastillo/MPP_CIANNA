#ifndef __MPP_H__
#define __MPP_H__

#include <bits/stdc++.h>
using namespace std;
#define FOREACH(i, v) for (__typeof((v).begin()) i = (v).begin(); i != (v).end(); i++)
//v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner, v_both_snack;

#define N_TIMES_DAY 7
//the times indexes should start with zero
#define BREAKFAST 0 
#define MORNING_SNACK 1
#define STARTER 2
#define MAIN_COURSE 3
#define EVENING_SNACK 4
#define DINNER 5
#define BOTH_SNACK 6

//extern volatile bool finished;
void printBest();
struct infoPlates {
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

		vector<infoPlates> v_dishes;
		vector<vector<int> > v_times_dishes;  //vector<int> v_breakfast, v_morning_snack, v_starter, v_main_course, v_evening_snack, v_dinner, v_both_snack;
		vector<constraint_nutrient> v_constraints;
		vector<int> v_constraint_global, v_constraint_day;
		int nDias;
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
		virtual void print(ostream &os) const;

//!menu planning

//		vector<bool> I;
		vector<int> x_var;
		double fitness;
		static MPP_Problem *MPP_problem;
	private:
		void calculateFeasibilityDegree();
		//int heaviestNut, heaviestType;
		//double valorFac, precioObj;
		//set<int> badDays;

};

#endif
