#ifndef __MPP_H__
#define __MPP_H__

#include <bits/stdc++.h>
using namespace std;
#define FOREACH(i, v) for (__typeof((v).begin()) i = (v).begin(); i != (v).end(); i++)

extern volatile bool finished;

void printBest();

struct infoPlatos {
	string nombre;				  //Nombre del plato
	float precio;           //Precio del plato
	float cantidad;					//Cantidad en gramos
	vector<float> infoN;    //Informacion nutricional del plato
	vector<string> alg;     //Alergenos del plato
	vector<string> inc;     //Incompatibilidades alimenticias del plato
	vector<int> gruposAl;
};

class MPP_Problem{
	public:
		MPP_Problem();
		~MPP_Problem(){
		
		}
		vector<infoPlatos> v_primerosPlatos;
		vector<infoPlatos> v_segundosPlatos;
		vector<infoPlatos> v_postres;
		int nDias;

};

class MPP{
	public:
		MPP(){
		
		}
		~MPP(){
		
		}

//////menu planning
		MPP();
		void evaluate();
		void restart();
		Individual *clone() const;
		bool init (const vector<string> &params);
		void dependentMutation(double pm);
		void dependentCrossover(MPP &i2);
		void uniformCrossover(MPP &i2);
		void uniform2Crossover(MPP &i2);
		void pairBasedCrossover(MPP &i2);
		void localSearch();
		double getDistance(MPP &ind2); 
		double inline getMaximum(const int i) const { cerr << "ErrorL llama a getMaximum" << endl; exit(-1); return 0; }
		double inline getMinimum(const int i) const { cerr << "Error: llama a getMinimum" << endl; exit(-1); return 0; }
		unsigned int inline getOptDirection(const int i) const { return MAXIMIZE; }
		virtual void print(ostream &os) const;

//!menu planning

		vector<bool> I;
		long long fitness;
		static MPP_Problem *MPP_problem;
	private:
		void calculateFeasibilityDegree();
		int heaviestNut, heaviestType;
		double valorFac, precioObj;
		set<int> badDays;

};

#endif
