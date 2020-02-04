#ifndef _MPP_INDIVIDUAL_H_
#define _MPP_INDIVIDUAL_H_

#include <set>
#include "Individual.h"
#include "MPP_Problem.h"

const int UNIFORM_CROSSOVER = 0;
const int PAIR_BASED_CROSSOVER = 1;
const int UNIFORM2_CROSSOVER = 2;
class MPP : public Individual {
	public:
		MPP();
		void evaluate();
		void restart();
		Individual *clone() const;
		bool init (const vector<string> &params);
		void dependentMutation(double pm);
		void dependentCrossover(Individual *i2);
		void uniformCrossover(Individual *i2);
		void uniform2Crossover(Individual *i2);
		void pairBasedCrossover(Individual *i2);
		void localSearch();
		double getDistance(Individual &ind2); 
		double inline getMaximum(const int i) const { cerr << "ErrorL llama a getMaximum" << endl; exit(-1); return 0; }
		double inline getMinimum(const int i) const { cerr << "Error: llama a getMinimum" << endl; exit(-1); return 0; }
		unsigned int inline getOptDirection(const int i) const { return MAXIMIZE; }
		virtual void print(ostream &os) const;

	private:
		void calculateFeasibilityDegree();
		int heaviestNut, heaviestType;
		double valorFac, precioObj;
		set<int> badDays;

};
#endif
