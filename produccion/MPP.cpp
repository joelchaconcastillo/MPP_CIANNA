#include <signal.h>
#include "MPP.h"
#include "utils.h"
using namespace std;

int currentTol;
long long Globalbest = 1e16;
MPP bestI;
int generation = 0;
bool volatile finished = false;
MPP_Problem* MPP::MPP_problem;
long long maxThreshold;


//////////////////////////// Funciones del individuo ////////////////////////////////
MPP_Problem::MPP_Problem(){

}

void MPP::calculateFeasibilityDegree(){
	valorFac = 0;
	double infoNPlan[num_nutr];
	bzero(infoNPlan, sizeof(infoNPlan));
	badDays.clear();
	for(int i = 0; i < nDias; i++){
		int x = i*num_tipoPlato;
		double dayNutr[num_nutr];
		bzero(dayNutr, sizeof(dayNutr));
		for(unsigned int j = 0; j < num_nutr; j++){
			dayNutr[j] += v_primerosPlatos[round(getVar(x))].infoN[j];
			dayNutr[j] += v_segundosPlatos[round(getVar(x+1))].infoN[j];
			dayNutr[j] += v_postres[round(getVar(x+2))].infoN[j];
			infoNPlan[j] += dayNutr[j];
		}
		for (int j = 0; j < FORCED_INDEXES_SIZE; j++){
			int index = FORCED_INDEXES[j];
			if (dayNutr[index] < ingR[index] * FORCED_MIN[j]){
				valorFac += pow((ingR[index] * FORCED_MIN[j] - dayNutr[index]) / ingR[index], 2) * 1000000.0;
				badDays.insert(i);
			} else if (dayNutr[index] > ingR[index] * FORCED_MAX[j]){
				valorFac += pow((dayNutr[index] - ingR[index] * FORCED_MAX[j]) / ingR[index], 2) * 1000000.0;
				badDays.insert(i);
			}
		}
	}


	double heaviestValue = 0;
	heaviestNut = -1;
	for(unsigned int i = 0; i < num_nutr; i++){
		if ((i == CALCIUM_INDEX) || (i == POTASIUM_INDEX) || (i == IRON_INDEX)) continue;
		if (infoNPlan[i] < ingR[i] * minReq[i] * nDias){
			double v = pow((ingR[i] * minReq[i] * nDias - infoNPlan[i]) / (ingR[i] * nDias), 2);
			valorFac += v;
			if (v > heaviestValue){
				heaviestValue = v;
				heaviestNut = i;
				heaviestType = -1;
			}
		} 
		if (infoNPlan[i] > ingR[i] * maxReq[i] * nDias){
			double v = pow((infoNPlan[i] - ingR[i] * maxReq[i] * nDias) / (ingR[i] * nDias), 2);
			valorFac += v;
			if (v > heaviestValue){
				heaviestValue = v;
				heaviestNut = i;
				heaviestType = 1;
			}
		}
	}
}

void MPP::evaluate(){
	precioObj = 0;
	for(int i = 0; i < nDias; i++){
		int x = i*num_tipoPlato;
		precioObj += v_primerosPlatos[round(getVar(x))].precio + v_segundosPlatos[round(getVar(x+1))].precio + v_postres[round(getVar(x+2))].precio;
	}
	calculateFeasibilityDegree();
	setObj(0, 1.0 / (valorFac * 1e12 + precioObj));
}

Individual *MPP::clone() const {
	MPP *newInd = new MPP();
	newInd->badDays = badDays;
	newInd->heaviestNut = heaviestNut;
	newInd->heaviestType = heaviestType;
	newInd->valorFac = valorFac;
	newInd->precioObj = precioObj;
	return newInd;
}

bool MPP::init(const vector<string> &params){//Dias Cruce
	if (params.size() != 2) { cerr << "Error parametros de MPP incorrectos" << endl; return false;}
	nDias = atoi(params[0].c_str());
	string crossover = params[1];
	if (crossover == "Uniform"){
		crossoverType = UNIFORM_CROSSOVER;
	} else if (crossover == "PairBasedCrossover"){
		crossoverType = PAIR_BASED_CROSSOVER;
	} else if (crossover == "Uniform2"){
		crossoverType = UNIFORM2_CROSSOVER;
	} else {
		cerr << "Error en la especificacion del crossover" << endl;
		return false;
	}
	setNumberOfVar(nDias * 3);
	setNumberOfObj(1);
	set_VectoresPlatos("/home/carlos.segura/svn/oplink/algorithms/team/src/plugins/problems/MPP/instancias/postres.txt", v_postres);
	set_VectoresPlatos("/home/carlos.segura/svn/oplink/algorithms/team/src/plugins/problems/MPP/instancias/primerosplatos.txt", v_primerosPlatos); 
	set_VectoresPlatos("/home/carlos.segura/svn/oplink/algorithms/team/src/plugins/problems/MPP/instancias/segundosplatos.txt", v_segundosPlatos);
	return true;
}

void MPP::restart(){
	cout << "Llama a inicializar" << endl;
	for (int i = 0; i < nDias; i++){
		for (int j = 0; j < 3; j++){
			setVar(i * 3 + j, random() % NPLATOS[j]);
		}
	}
	evaluate();
	localSearch();
}

struct Food {
	int p1, p2, p3;
	bool operator<(const Food &i2) const {
		return (make_pair(p1, make_pair(p2, p3)) < make_pair(i2.p1, make_pair(i2.p2, i2.p3)));
	}
}; 

void MPP::dependentCrossover(Individual *i2){
	if (crossoverType == PAIR_BASED_CROSSOVER){
		pairBasedCrossover(i2);
	} else if (crossoverType == UNIFORM_CROSSOVER){
		uniformCrossover(i2);
	} else if (crossoverType == UNIFORM2_CROSSOVER){
		uniform2Crossover(i2);
	}
}

void MPP::uniformCrossover(Individual *i2){
	for (int i = 0; i < nDias; i++){
		if (rand() > (RAND_MAX / 2)){
			swap(var[i*3], i2->var[i*3]);
			swap(var[i*3+1], i2->var[i*3+1]);
			swap(var[i*3+2], i2->var[i*3+2]);
		} 
	}
}

void MPP::uniform2Crossover(Individual *i2){
	for (int i = 0; i < nDias * 3; i++){
		if (rand() > (RAND_MAX / 2)){
			swap(var[i], i2->var[i]);
		} 
	}
}

void MPP::pairBasedCrossover(Individual *i2){
	vector<Food> pendingI1, pendingI2;

	map<Food, int> f1;
	int dist = 0;
	for (int i = 0; i < nDias; i++){
		Food f;
		f.p1 = round(getVar(i*3));
		f.p2 = round(getVar(i*3 + 1));
		f.p3 = round(getVar(i*3 + 2));
		f1[f]++;
	}

	int fixed = 0;
	for (int i = 0; i < nDias; i++){
		Food f;
		f.p1 = round(i2->getVar(i*3));
		f.p2 = round(i2->getVar(i*3 + 1));
		f.p3 = round(i2->getVar(i*3 + 2));
		if (f1.count(f)){//Comida en ambos
			setVar(fixed*3, f.p1);
			setVar(fixed*3+1, f.p2);
			setVar(fixed*3+2, f.p3);
			i2->setVar(fixed*3, f.p1);
			i2->setVar(fixed*3+1, f.p2);
			i2->setVar(fixed*3+2, f.p3);
			fixed++;
			f1[f]--;
			if (f1[f] == 0){
				f1.erase(f);
			}
		} else {
			pendingI2.push_back(f);
			dist++;
		}
	}
	for (map<Food, int>::iterator it = f1.begin(); it != f1.end(); it++){
		for (int j = 0; j < it->second; j++){
			pendingI1.push_back(it->first);
		}
	}
	if (pendingI1.size() != pendingI2.size()){ cerr << "Error interno. PendingI1 != PendingI2" << endl; exit(-1); }
	random_shuffle(pendingI1.begin(), pendingI1.end());
	int next = nDias - pendingI1.size();
	for (int i = 0; i < pendingI1.size(); i++){
		Food f1 = pendingI1[i];
		Food f2 = pendingI2[i];
		if (rand() < RAND_MAX / 2.0){
			swap(f1, f2);
		} 
		setVar(next * 3, f1.p1);
		setVar(next * 3 + 1, f1.p2);
		setVar(next * 3 + 2, f1.p3);
		i2->setVar(next * 3, f2.p1);
		i2->setVar(next * 3 + 1, f2.p2);
		i2->setVar(next * 3 + 2, f2.p3);
		next++;
	}
}

void MPP::dependentMutation(double pm){
}
struct Neighbor {
	int variable;
	int newValue;
};

void MPP::localSearch( ) {
	vector<Neighbor> neighbors;
	for (int i = 0; i < nDias; i++){
		for (int j = 0; j < 3; j++){
			for (int k = 0; k < NPLATOS[j]; k++){
				Neighbor n;
				n.variable = i * 3 + j;
				n.newValue = k;
				neighbors.push_back(n);
			}
		}
	}
	vector<double> bestIndividual = var;
	evaluate();
	pair<double, double> bestResult = make_pair(valorFac, precioObj);
	for (int i = 0; i < 100; i++){
		evaluate();
		pair<double, double> currentResult = make_pair(valorFac, precioObj);
		bool improved = true;
		while(improved){
			improved = false;
			random_shuffle(neighbors.begin(), neighbors.end());
			for (int i = 0; i < neighbors.size(); i++){
				int currentValue = var[neighbors[i].variable];
				var[neighbors[i].variable] = neighbors[i].newValue;
				evaluate();
				pair<double, double> newResult = make_pair(valorFac, precioObj);
				if (newResult >= currentResult){
					var[neighbors[i].variable] = currentValue;
				} else {
					improved = true;
					currentResult = newResult;
				}
			}
		}

		if (currentResult >= bestResult){
			var = bestIndividual;
		} else {
			bestResult = currentResult;
			bestIndividual = var;
		}

		evaluate();
		if (badDays.size() == 0){
			if (heaviestNut != -1){
				vector< pair<double, int> > infoNut;
				for (int i = 0; i < nDias; i++){
					double total = v_primerosPlatos[var[i*3]].infoN[heaviestNut] + 
												 v_segundosPlatos[var[i*3 + 1]].infoN[heaviestNut] + 
											   v_postres[var[i*3 + 2]].infoN[heaviestNut];
					infoNut.push_back(make_pair(total, i));
				}
				sort(infoNut.begin(), infoNut.end());
				if (heaviestType == 1) reverse(infoNut.begin(), infoNut.end());
				int selectedDay = infoNut[random() % 6].second;
				var[selectedDay * 3] = (rand() % NPLATOS[0]);
				var[selectedDay * 3 + 1] = (rand() % NPLATOS[1]);
				var[selectedDay * 3 + 2] = (rand() % NPLATOS[2]);
			} else {
				int selectedDay = random() % nDias;
				//cout << "Selecciona dia " << selectedDay << endl;
				var[selectedDay * 3] = (rand() % NPLATOS[0]);
				var[selectedDay * 3 + 1] = (rand() % NPLATOS[1]);
				var[selectedDay * 3 + 2] = (rand() % NPLATOS[2]);
			}
		} else {
			for (auto it = badDays.begin(); it != badDays.end(); it++){
				int day = *it;
				int which = random() % 3;
				var[day * 3 + which] = (rand() % NPLATOS[which]); 
			}
		}
	}
	var = bestIndividual;
	evaluate();
	//cout << "Final " << valorFac << " " << precioObj << endl;
}


double MPP::getDistance(Individual &ind2) {
	map<Food, int> f1;
	int dist = 0;
	for (int i = 0; i < nDias; i++){
		Food f;
		f.p1 = round(getVar(i*3));
		f.p2 = round(getVar(i*3 + 1));
		f.p3 = round(getVar(i*3 + 2));
		f1[f]++;
	}
	for (int i = 0; i < nDias; i++){
		Food f;
		f.p1 = round(ind2.getVar(i*3));
		f.p2 = round(ind2.getVar(i*3 + 1));
		f.p3 = round(ind2.getVar(i*3 + 2));
		if (f1.count(f)){
			f1[f]--;
			if (f1[f] == 0){
				f1.erase(f);
			}
		} else {
			dist++;
		}
	}
	return dist;
}

void MPP::print(ostream &os) const {
	for (int i = 0; i < var.size(); i++){
		os << (int)round(var[i]) << " ";
	}
	os << valorFac << " " << precioObj << endl;
}

long long fitnessBest = 1e15;


