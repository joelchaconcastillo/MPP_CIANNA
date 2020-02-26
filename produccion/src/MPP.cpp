#include <signal.h>
#include "MPP.h"
#include "utils.h"
using namespace std;

const string WHITESPACE = " \n\r\t\f\v";

string ltrim(const string& s)
{
	size_t start = s.find_first_not_of(WHITESPACE);
	return (start == std::string::npos) ? "" : s.substr(start);
}

string rtrim(const string& s)
{
	size_t end = s.find_last_not_of(WHITESPACE);
	return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

string trim(const string& s)
{
	return rtrim(ltrim(s));
}

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

void MPP_Problem::load_data(int argc, char **argv)
{
    if(argc < 6)
    {
	 cout << "N\'umero de argumentos inv\'alidos" <<endl;
	 exit(EXIT_FAILURE);
    }
    this->nDias = atoi(argv[3]);
    ////reading the information.......
    v_times_dishes.resize(N_TIMES_DAY); //N_times_dishes...
    load_constraints(argv[2]); 
    load_dishes(argv[1]);
}
void MPP_Problem::load_dishes(char *c_filename)
{

   	ifstream ifs;
	struct infoPlates str_dish;
	ifs.open(c_filename, ifstream::in);
	if (ifs.is_open())
	{
	   string line, word;
	   getline(ifs, line);
	   stringstream line_commas(line);
	   vector<string> column_names;
	//Here is assumed that the first two columns are of the "DESCRIPCION" and "TIEMPO" respectively, thereafter are the meta-data and the last two columns are "CATEGORIA" and "FAVORITO" respectively.
	   while(getline(line_commas ,word, ',')) //first getting the tag-information of each column....
	   column_names.push_back(word);
	   while (ifs.good())
	   {
	       string cell;
	       getline(ifs, cell, ',');
	       if(trim(cell).empty()) break; //The file has an extra empty line
               str_dish.description = stoi(cell);
	       getline(ifs, cell, ',');
               str_dish.time_day = trim(cell);
               for(int i = 0; i < column_names.size()-4; i++)
               {
	          getline(ifs, cell, ',');
		  str_dish.v_nutrient_value.push_back(stod(cell));
               }
	       getline(ifs, cell, ',');
	       str_dish.category = stoi(cell);
	       getline(ifs, cell, '\n');
               str_dish.favorite = (bool)stoi(cell);
		
	       if(str_dish.time_day == "DESAYUNO") v_times_dishes[BREAKFAST].push_back(v_dishes.size());
	       else if(str_dish.time_day == "COLACION_MATUTINA") v_times_dishes[MORNING_SNACK].push_back(v_dishes.size());
	       else if(str_dish.time_day == "COMIDA_ENTRADA") v_times_dishes[STARTER].push_back(v_dishes.size());
	       else if(str_dish.time_day == "COMIDA_PRINCIPAL") v_times_dishes[MAIN_COURSE].push_back(v_dishes.size());
	       else if(str_dish.time_day == "COLACION_VESPERTINA") v_times_dishes[EVENING_SNACK].push_back(v_dishes.size());
	       else if(str_dish.time_day == "CENA") v_times_dishes[DINNER].push_back(v_dishes.size());
	       else if(str_dish.time_day == "COLACION_AMBAS") v_times_dishes[BOTH_SNACK].push_back(v_dishes.size());
	       else
		{
		   cout <<"Tiempo del día no desconocido"<<endl;
		   exit(EXIT_FAILURE);
		}
	       v_dishes.push_back(str_dish);
	    }
		ifs.close();
	} else {
		cout << "\n\nError. No se ha podido leer el archivo de platos."<<endl;
		exit(EXIT_FAILURE);
	}
}
void MPP_Problem::load_constraints(char *c_filename)
{
   	ifstream ifs;
	struct constraint_nutrient str_constraint_nutrient;
	ifs.open(c_filename, ifstream::in);
	if (ifs.is_open())
	{
	   string line, word;
	  //headers.....
	   getline(ifs, line);
	   stringstream line_commas(line);
	   vector<string> column_names;
	   while(getline(line_commas ,word, ',')) //first getting the tag-information of each column....
	   column_names.push_back(word);
	   while (ifs.good())
	   {
	       string cell;
	       getline(ifs, cell, ',');
	       if(trim(cell).empty()) break; //The file has an extra empty line
               str_constraint_nutrient.type = trim(cell);
	       getline(ifs, cell, ',');
	       str_constraint_nutrient.name = trim(cell);
	       getline(ifs, cell, ',');
	       str_constraint_nutrient.min = stod(cell);
	       getline(ifs, cell, '\n'); //last word...
               str_constraint_nutrient.max = stod(cell);

	       if( str_constraint_nutrient.type == "GLOBAL") v_constraint_global.push_back((int)v_constraints.size());
	       else if( str_constraint_nutrient.type == "DIA") v_constraint_day.push_back((int)v_constraints.size());
	       else{
	   cout << str_constraint_nutrient.type<<endl;
		   cout << "Se desconoce un tipo de restricci\'on, \'unicamente puede ser por d\'ia o global"<<endl;
		   exit(EXIT_FAILURE);
		}
	       v_constraints.push_back(str_constraint_nutrient);
	    }
		ifs.close();
	} else {
		cout << "\n\nError. No se ha podido leer el archivo de restricciones."<<endl;
		exit(EXIT_FAILURE);
	}
}
void MPP::calculateFeasibilityDegree(){
///	valorFac = 0;
///	double infoNPlan[num_nutr];
///	bzero(infoNPlan, sizeof(infoNPlan));
///	badDays.clear();
///	for(int i = 0; i < nDias; i++){
///		int x = i*num_tipoPlato;
///		double dayNutr[num_nutr];
///		bzero(dayNutr, sizeof(dayNutr));
///		for(unsigned int j = 0; j < num_nutr; j++){
///			dayNutr[j] += v_primerosPlatos[round(getVar(x))].infoN[j];
///			dayNutr[j] += v_segundosPlatos[round(getVar(x+1))].infoN[j];
///			dayNutr[j] += v_postres[round(getVar(x+2))].infoN[j];
///			infoNPlan[j] += dayNutr[j];
///		}
///		for (int j = 0; j < FORCED_INDEXES_SIZE; j++){
///			int index = FORCED_INDEXES[j];
///			if (dayNutr[index] < ingR[index] * FORCED_MIN[j]){
///				valorFac += pow((ingR[index] * FORCED_MIN[j] - dayNutr[index]) / ingR[index], 2) * 1000000.0;
///				badDays.insert(i);
///			} else if (dayNutr[index] > ingR[index] * FORCED_MAX[j]){
///				valorFac += pow((dayNutr[index] - ingR[index] * FORCED_MAX[j]) / ingR[index], 2) * 1000000.0;
///				badDays.insert(i);
///			}
///		}
///	}
///
///
///	double heaviestValue = 0;
///	heaviestNut = -1;
///	for(unsigned int i = 0; i < num_nutr; i++){
///		if ((i == CALCIUM_INDEX) || (i == POTASIUM_INDEX) || (i == IRON_INDEX)) continue;
///		if (infoNPlan[i] < ingR[i] * minReq[i] * nDias){
///			double v = pow((ingR[i] * minReq[i] * nDias - infoNPlan[i]) / (ingR[i] * nDias), 2);
///			valorFac += v;
///			if (v > heaviestValue){
///				heaviestValue = v;
///				heaviestNut = i;
///				heaviestType = -1;
///			}
///		} 
///		if (infoNPlan[i] > ingR[i] * maxReq[i] * nDias){
///			double v = pow((infoNPlan[i] - ingR[i] * maxReq[i] * nDias) / (ingR[i] * nDias), 2);
///			valorFac += v;
///			if (v > heaviestValue){
///				heaviestValue = v;
///				heaviestNut = i;
///				heaviestType = 1;
///			}
///		}
///	}
}

void MPP::evaluate(){
//	precioObj = 0;
//	for(int i = 0; i < nDias; i++){
//		int x = i*num_tipoPlato;
//		precioObj += v_primerosPlatos[round(getVar(x))].precio + v_segundosPlatos[round(getVar(x+1))].precio + v_postres[round(getVar(x+2))].precio;
//	}
//	calculateFeasibilityDegree();
//	setObj(0, 1.0 / (valorFac * 1e12 + precioObj));
}

//Individual *MPP::clone() const {
//	MPP *newInd = new MPP();
//	newInd->badDays = badDays;
//	newInd->heaviestNut = heaviestNut;
//	newInd->heaviestType = heaviestType;
//	newInd->valorFac = valorFac;
//	newInd->precioObj = precioObj;
//	return newInd;
//}

void MPP::init(){
   x_var.resize(N_TIMES_DAY*MPP_problem->nDias);
   for (int i = 0; i < MPP_problem->nDias; i++){
	for (int j = 0; j < N_TIMES_DAY; j++){
            x_var[i*N_TIMES_DAY+j] = rand()%(int)MPP_problem->v_times_dishes[j].size();
	   }
	}
}

void MPP::restart(){
	for (int i = 0; i < MPP_problem->nDias; i++){
		for (int j = 0; j < N_TIMES_DAY; j++){
		        x_var[i*N_TIMES_DAY+j] = rand()%(int)MPP_problem->v_times_dishes[j].size();
		}
	}
	evaluate();
	localSearch();
}

//struct Food {
//	int p1, p2, p3;
//	bool operator<(const Food &i2) const {
//		return (make_pair(p1, make_pair(p2, p3)) < make_pair(i2.p1, make_pair(i2.p2, i2.p3)));
//	}
//}; 

void MPP::dependentCrossover(MPP &i2){
//	if (crossoverType == PAIR_BASED_CROSSOVER){
//		pairBasedCrossover(i2);
//	} else if (crossoverType == UNIFORM_CROSSOVER){
//		uniformCrossover(i2);
//	} else if (crossoverType == UNIFORM2_CROSSOVER){
//		uniform2Crossover(i2);
//	}
}

void MPP::uniformCrossover(MPP &i2)
{
//	for (int i = 0; i < nDias; i++){
//		if (rand() > (RAND_MAX / 2)){
//			swap(var[i*3], i2->var[i*3]);
//			swap(var[i*3+1], i2->var[i*3+1]);
//			swap(var[i*3+2], i2->var[i*3+2]);
//		} 
//	}
}

void MPP::uniform2Crossover(MPP &i2){
//	for (int i = 0; i < nDias * 3; i++){
//		if (rand() > (RAND_MAX / 2)){
//			swap(var[i], i2->var[i]);
//		} 
//	}
}

void MPP::pairBasedCrossover(MPP &i2)
{
//	vector<Food> pendingI1, pendingI2;
//
//	map<Food, int> f1;
//	int dist = 0;
//	for (int i = 0; i < nDias; i++){
//		Food f;
//		f.p1 = round(getVar(i*3));
//		f.p2 = round(getVar(i*3 + 1));
//		f.p3 = round(getVar(i*3 + 2));
//		f1[f]++;
//	}
//
//	int fixed = 0;
//	for (int i = 0; i < nDias; i++){
//		Food f;
//		f.p1 = round(i2->getVar(i*3));
//		f.p2 = round(i2->getVar(i*3 + 1));
//		f.p3 = round(i2->getVar(i*3 + 2));
//		if (f1.count(f)){//Comida en ambos
//			setVar(fixed*3, f.p1);
//			setVar(fixed*3+1, f.p2);
//			setVar(fixed*3+2, f.p3);
//			i2->setVar(fixed*3, f.p1);
//			i2->setVar(fixed*3+1, f.p2);
//			i2->setVar(fixed*3+2, f.p3);
//			fixed++;
//			f1[f]--;
//			if (f1[f] == 0){
//				f1.erase(f);
//			}
//		} else {
//			pendingI2.push_back(f);
//			dist++;
//		}
//	}
//	for (map<Food, int>::iterator it = f1.begin(); it != f1.end(); it++){
//		for (int j = 0; j < it->second; j++){
//			pendingI1.push_back(it->first);
//		}
//	}
//	if (pendingI1.size() != pendingI2.size()){ cerr << "Error interno. PendingI1 != PendingI2" << endl; exit(-1); }
//	random_shuffle(pendingI1.begin(), pendingI1.end());
//	int next = nDias - pendingI1.size();
//	for (int i = 0; i < pendingI1.size(); i++){
//		Food f1 = pendingI1[i];
//		Food f2 = pendingI2[i];
//		if (rand() < RAND_MAX / 2.0){
//			swap(f1, f2);
//		} 
//		setVar(next * 3, f1.p1);
//		setVar(next * 3 + 1, f1.p2);
//		setVar(next * 3 + 2, f1.p3);
//		i2->setVar(next * 3, f2.p1);
//		i2->setVar(next * 3 + 1, f2.p2);
//		i2->setVar(next * 3 + 2, f2.p3);
//		next++;
//	}
}

void MPP::dependentMutation(double pm){
}
//struct Neighbor {
//	int variable;
//	int newValue;
//};

void MPP::localSearch( ) {
///	vector<Neighbor> neighbors;
///	for (int i = 0; i < nDias; i++){
///		for (int j = 0; j < 3; j++){
///			for (int k = 0; k < NPLATOS[j]; k++){
///				Neighbor n;
///				n.variable = i * 3 + j;
///				n.newValue = k;
///				neighbors.push_back(n);
///			}
///		}
///	}
///	vector<double> bestIndividual = var;
///	evaluate();
///	pair<double, double> bestResult = make_pair(valorFac, precioObj);
///	for (int i = 0; i < 100; i++){
///		evaluate();
///		pair<double, double> currentResult = make_pair(valorFac, precioObj);
///		bool improved = true;
///		while(improved){
///			improved = false;
///			random_shuffle(neighbors.begin(), neighbors.end());
///			for (int i = 0; i < neighbors.size(); i++){
///				int currentValue = var[neighbors[i].variable];
///				var[neighbors[i].variable] = neighbors[i].newValue;
///				evaluate();
///				pair<double, double> newResult = make_pair(valorFac, precioObj);
///				if (newResult >= currentResult){
///					var[neighbors[i].variable] = currentValue;
///				} else {
///					improved = true;
///					currentResult = newResult;
///				}
///			}
///		}
///
///		if (currentResult >= bestResult){
///			var = bestIndividual;
///		} else {
///			bestResult = currentResult;
///			bestIndividual = var;
///		}
///
///		evaluate();
///		if (badDays.size() == 0){
///			if (heaviestNut != -1){
///				vector< pair<double, int> > infoNut;
///				for (int i = 0; i < nDias; i++){
///					double total = v_primerosPlatos[var[i*3]].infoN[heaviestNut] + 
///												 v_segundosPlatos[var[i*3 + 1]].infoN[heaviestNut] + 
///											   v_postres[var[i*3 + 2]].infoN[heaviestNut];
///					infoNut.push_back(make_pair(total, i));
///				}
///				sort(infoNut.begin(), infoNut.end());
///				if (heaviestType == 1) reverse(infoNut.begin(), infoNut.end());
///				int selectedDay = infoNut[random() % 6].second;
///				var[selectedDay * 3] = (rand() % NPLATOS[0]);
///				var[selectedDay * 3 + 1] = (rand() % NPLATOS[1]);
///				var[selectedDay * 3 + 2] = (rand() % NPLATOS[2]);
///			} else {
///				int selectedDay = random() % nDias;
///				//cout << "Selecciona dia " << selectedDay << endl;
///				var[selectedDay * 3] = (rand() % NPLATOS[0]);
///				var[selectedDay * 3 + 1] = (rand() % NPLATOS[1]);
///				var[selectedDay * 3 + 2] = (rand() % NPLATOS[2]);
///			}
///		} else {
///			for (auto it = badDays.begin(); it != badDays.end(); it++){
///				int day = *it;
///				int which = random() % 3;
///				var[day * 3 + which] = (rand() % NPLATOS[which]); 
///			}
///		}
///	}
///	var = bestIndividual;
///	evaluate();
	//cout << "Final " << valorFac << " " << precioObj << endl;
}


int MPP::getDistance(MPP &ind2) {
return 0.0;
//	map<Food, int> f1;
//	int dist = 0;
//	for (int i = 0; i < nDias; i++){
//		Food f;
//		f.p1 = round(getVar(i*3));
//		f.p2 = round(getVar(i*3 + 1));
//		f.p3 = round(getVar(i*3 + 2));
//		f1[f]++;
//	}
//	for (int i = 0; i < nDias; i++){
//		Food f;
//		f.p1 = round(ind2.getVar(i*3));
//		f.p2 = round(ind2.getVar(i*3 + 1));
//		f.p3 = round(ind2.getVar(i*3 + 2));
//		if (f1.count(f)){
//			f1[f]--;
//			if (f1[f] == 0){
//				f1.erase(f);
//			}
//		} else {
//			dist++;
//		}
//	}
//	return dist;
}

void MPP::print(ostream &os) const {
//	for (int i = 0; i < var.size(); i++){
//		os << (int)round(var[i]) << " ";
//	}
//	os << valorFac << " " << precioObj << endl;
}

long long fitnessBest = 1e15;


