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

MPP_Problem* MPP::MPP_problem;

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
    out_filename = string(argv[5]);
    nDias = atoi(argv[3]);

    ////reading the information.......
    v_times_dishes.resize(N_OPT_DAY); //N_times_dishes...
    load_constraints(argv[2]); 
    load_dishes(argv[1]);
}
void MPP_Problem::load_dishes(char *c_filename)
{

   	ifstream ifs;
	struct infoDishes str_dish;
	ifs.open(c_filename, ifstream::in);
	if (ifs.is_open())
	{
	   string line, word;
	   getline(ifs, line);
	   stringstream line_commas(line);
	   vector<string> column_names;
	//Here is assumed that the first two columns are of the "DESCRIPCION" and "TIEMPO" respectively, thereafter are the meta-data and the last two columns are "CATEGORIA" and "FAVORITO" respectively.
	   while(getline(line_commas ,word, ',')) //first getting the tag-information of each column....
	   column_names.push_back(trim(word));
	   while (ifs.good())
	   {
	       string cell;

	       str_dish.v_nutrient_value.resize(v_constraints.size(), 0.0); 
               for(int i = 0; i < column_names.size(); i++)
               {
		  if(i==column_names.size()-1)
	             getline(ifs, cell, '\n');
		  else
	             getline(ifs, cell, ',');
		  cell = trim(cell);
	       if(cell.empty()) break; //The file has an extra empty line
		  if( column_names[i] == "DESCRIPCION") 
               		str_dish.description = stoi(cell);
		  else if( column_names[i] == "TIEMPO")
               		str_dish.time_day = trim(cell);
		  else if( column_names[i] == "CATEGORIA")
	       		str_dish.category = stoi(cell);
		  else if( column_names[i] == "FAVORITO")
               		str_dish.favorite = (bool)stoi(cell);
		  else{
			if( dic_nut_id.find(column_names[i]) != dic_nut_id.end())
		           str_dish.v_nutrient_value[dic_nut_id[column_names[i]]] = stod(cell); //the nutrient values need to be stored in the same order that the contraints..
			//else { cout << column_names[i]<<endl;}// cout << "error interno"<<endl; exit(EXIT_FAILURE);}
			}
               }
	       if(cell.empty()) break; //The file has an extra empty line

		
	       if(str_dish.time_day == "DESAYUNO") v_times_dishes[BREAKFAST].push_back(str_dish);
	       else if(str_dish.time_day == "COLACION_MATUTINA") v_times_dishes[MORNING_SNACK].push_back(str_dish);
	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_1) v_times_dishes[STARTER_1].push_back(str_dish);
	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_2) v_times_dishes[STARTER_2].push_back(str_dish);
	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_BOTH) //this dish blong to both categories...
		{
		 v_times_dishes[STARTER_1].push_back(str_dish);
		 v_times_dishes[STARTER_2].push_back(str_dish);
		}
	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_1) v_times_dishes[MAIN_COURSE_1].push_back(str_dish);
	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_2) v_times_dishes[MAIN_COURSE_2].push_back(str_dish);
	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_BOTH) //both categories..
		{
		  v_times_dishes[MAIN_COURSE_1].push_back(str_dish);
		  v_times_dishes[MAIN_COURSE_2].push_back(str_dish);
		}
	       else if(str_dish.time_day == "COLACION_VESPERTINA") v_times_dishes[EVENING_SNACK].push_back(str_dish);
	       else if(str_dish.time_day == "CENA") v_times_dishes[DINNER].push_back(str_dish);
	       else if(str_dish.time_day == "COLACION_AMBAS")
	       {
 		 v_times_dishes[MORNING_SNACK].push_back(str_dish);
 		 v_times_dishes[EVENING_SNACK].push_back(str_dish);
	       }
	       else
	       {
		   cout <<"Tiempo del dia desconocido"<<endl;
		   exit(EXIT_FAILURE);
	       }
//	      if( str_dish.category == CATEGORY_BOTH)
//	      {
//		str_dish.category = CATEGORY_1;
//	       v_dishes.push_back(str_dish);
//		str_dish.category = CATEGORY_2;
//	       v_dishes.push_back(str_dish);
//	      }
//	      else
//	       v_dishes.push_back(str_dish);
	    }
		ifs.close();
	} else {
		cout << "\n\nError. No se ha podido leer el archivo de platos."<<endl;
		exit(EXIT_FAILURE);
	}
	cout << "platillos... " <<endl;
       for(int i = 0; i < N_OPT_DAY; i++) cout << v_times_dishes[i].size() << " ";
       cout << endl;
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
	   column_names.push_back(trim(word));
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
	       dic_nut_id[str_constraint_nutrient.name] = (int) v_constraints.size();
	       if( str_constraint_nutrient.type == "GLOBAL") v_constraint_global.push_back((int)v_constraints.size());
	       else if( str_constraint_nutrient.type == "DIARIA") v_constraint_day.push_back((int)v_constraints.size());
	       else{
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
void MPP::calculateFeasibilityDegree2(){
	int num_nutr = (int)MPP_problem->v_constraints.size();
	vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
	double infoNPlan[num_nutr];
	bzero(infoNPlan, sizeof(infoNPlan));
	cout << "diario"<<endl;
	for(int i = 0; i < nDias; i++){
	cout << "Dia: " <<i<<endl;
		int x = i*N_OPT_DAY;
		double dayNutr[num_nutr];
		bzero(dayNutr, sizeof(dayNutr));
		for(unsigned int j = 0; j < num_nutr; j++){
			for(unsigned int k = 0; k < N_OPT_DAY; k++)
		 	   dayNutr[j] += MPP_problem->v_times_dishes[k][x_var[x+k]].v_nutrient_value[j];
			infoNPlan[j] += dayNutr[j]; //for global nutr..
		}
               //////////daily nutrients...
		for (int j = 0; j < (int)MPP_problem->v_constraint_day.size(); j++){
			int index = MPP_problem->v_constraint_day[j];
			double min = v_constraints[index].min;
			double max = v_constraints[index].max;
			if(dayNutr[index] < min){
				cout <<"<" <<v_constraints[index].name <<":"<< dayNutr[index]<<"--"<<min << endl;
			}
			 if (dayNutr[index] > max){
				cout <<">"<<v_constraints[index].name <<":"<< dayNutr[index]<<"--"<<max <<endl;
			}
		}
	cout <<endl<<endl;
	}
	cout << "global"<<endl;
	////////////////global nutriments...
	double heaviestValue = 0;
	heaviestNut = -1;
	for(unsigned int i = 0; i < (int) MPP_problem->v_constraint_global.size(); i++){
		int index = MPP_problem->v_constraint_global[i];
		double min = v_constraints[index].min*nDias;
		double max = v_constraints[index].max*nDias;
		if (infoNPlan[index] < min){
			cout <<"<"<<v_constraints[index].name <<":"<< infoNPlan[index]<<"--"<<min <<endl;
		} 
		if (infoNPlan[index] > max){
			cout <<">"<<v_constraints[index].name <<":"<< infoNPlan[index]<<"--"<<max <<endl;
		}
	}
}
void MPP::calculateFeasibilityDegree(){
	valorFac = 0.0;
	int num_nutr = (int)MPP_problem->v_constraints.size();
	vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
	double infoNPlan[num_nutr];
	bzero(infoNPlan, sizeof(infoNPlan));
	badDays.clear();
	for(int i = 0; i < nDias; i++){
		int x = i*N_OPT_DAY;
		double dayNutr[num_nutr];
		bzero(dayNutr, sizeof(dayNutr));
		for(unsigned int j = 0; j < num_nutr; j++){
			for(unsigned int k = 0; k < N_OPT_DAY; k++)
		 	   dayNutr[j] += MPP_problem->v_times_dishes[k][x_var[x+k]].v_nutrient_value[j];
			infoNPlan[j] += dayNutr[j]; //for global nutr..
		}
               //////////daily nutrients...
		for (int j = 0; j < (int)MPP_problem->v_constraint_day.size(); j++){
			int index = MPP_problem->v_constraint_day[j];
			double minv = v_constraints[index].min;
			double maxv = v_constraints[index].max;
			if(dayNutr[index] < minv){
				
				valorFac +=pow((minv - dayNutr[index])/minv, 2.0)*1.0e6;
				//valorFac +=pow((min - dayNutr[index])/((max+min)*0.5), 4.0)*1.0e6;//	 pow((ingR[index] * FORCED_MIN[j] - dayNutr[index]) / ingR[index], 2) * 1000000.0;
		//		cout << "-"<<v_constraints[index].name <<"-"<< dayNutr[index]<<"__"<<min <<" ";
				badDays.insert(i);
			}
			 if (dayNutr[index] > maxv){
				valorFac +=pow((dayNutr[index] - maxv)/maxv, 2.0)*1.0e6;
			//	valorFac +=pow((dayNutr[index]-max)/((max+min)*0.5), 4.0)*1.0e6;// pow((dayNutr[index] - ingR[index] * FORCED_MAX[j]) / ingR[index], 2) * 1000000.0;
	//			cout << "-"<<v_constraints[index].name <<"-"<< dayNutr[index]<<"__"<<max <<" ";
				badDays.insert(i);
			}
		}
	//cout <<endl<<endl;
	}
	////////////////global nutriments...
	double heaviestValue = 0;
	heaviestNut = -1;
	for(unsigned int i = 0; i < (int) MPP_problem->v_constraint_global.size(); i++){
		int index = MPP_problem->v_constraint_global[i];
		double min = v_constraints[index].min*nDias;
		double max = v_constraints[index].max*nDias;

		if (infoNPlan[index] < min){
			//double v = pow((min - infoNPlan[index])/((max+min)*0.5), 2.0);//pow((ingR[i] * minReq[i] * nDias - infoNPlan[i]) / (ingR[i] * nDias), 2);
			double v = pow((min - infoNPlan[index])/min, 2.0);//pow((ingR[i] * minReq[i] * nDias - infoNPlan[i]) / (ingR[i] * nDias), 2);
			valorFac += v;
			if (v > heaviestValue){
				heaviestValue = v;
				heaviestNut = index;
				heaviestType = -1;
			}
		} 
		if (infoNPlan[index] > max){
			double v = pow( (infoNPlan[index] - max)/max, 2.0);//pow((infoNPlan[i] - ingR[i] * maxReq[i] * nDias) / (ingR[i] * nDias), 2);
			valorFac += v;
			if (v > heaviestValue){
				heaviestValue = v;
				heaviestNut = index;
				heaviestType = 1;
			}
		}
	}
}
/*
 The objective is defined as follows:
  max f(x) + g(y)
  where 
     x is the varaibility by day 
     y is the global variability 
*/
void MPP::evaluate(){
    variabilidadObj = 0.0;
    vector<vector<infoDishes> > &v_times_dishes = MPP_problem->v_times_dishes;
    double variability_day= 0.0, variability_global = 0.0;
    //fitness by day...    
   int NN = v_times_dishes.size();
    vector<vector<int>> minv_d(NN), cont_d(NN), dish_d(NN), last_d(NN);
    for(int k = 0; k < N_OPT_DAY; k++)
    {
	minv_d[k] = vector<int>(v_times_dishes[k].size(), INT_MAX);
	cont_d[k] = vector<int>(v_times_dishes[k].size(), 0);
	dish_d[k] = vector<int>(v_times_dishes[k].size(), 0);
	last_d[k] = vector<int>(v_times_dishes[k].size(), -1);
    }

   for(int i = 0 ; i < nDias; i++)
   {
       if(v_times_dishes[STARTER_1][x_var[i*N_OPT_DAY + STARTER_1]].description != v_times_dishes[STARTER_2][x_var[i*N_OPT_DAY + STARTER_2]].description ) variability_day++;
       if(v_times_dishes[MAIN_COURSE_1][x_var[i*N_OPT_DAY + MAIN_COURSE_1]].description != v_times_dishes[MAIN_COURSE_2][x_var[i*N_OPT_DAY + MAIN_COURSE_2]].description ) variability_day++;
       if(v_times_dishes[MORNING_SNACK][x_var[i*N_OPT_DAY + MORNING_SNACK]].description != v_times_dishes[EVENING_SNACK][x_var[i*N_OPT_DAY + EVENING_SNACK]].description ) variability_day++;
     for(int k = 0; k < N_OPT_DAY; k++)
     {
       int id_dish = x_var[i*N_OPT_DAY + k];
       if(last_d[k][id_dish] != -1)
	{
	  int diff = i - last_d[k][id_dish];
	  if( minv_d[k][id_dish] > diff)
	  {
	   minv_d[k][id_dish] = diff;
	   cont_d[k][id_dish] = 0;
	  }
	  else if(minv_d[k][id_dish] == diff )
	    cont_d[k][id_dish]++;
	}
	last_d[k][id_dish] = i;
     }
   }
  for(int k = 0; k < N_OPT_DAY; k++)
   for(int i = 0; i < NN; i++)
   {
    if(last_d[k][i] == -1) continue;
     variability_global +=  minv_d[k][i] + cont_d[k][i];
   }
  variabilidadObj= variability_day + pow(variability_global,2);
  calculateFeasibilityDegree();
  //fitness = 1.0/(valorFac + 1.0e4/variabilidadObj);
  fitness =valorFac ;//- 1.0e3*variabilidadObj;// 1.0/(valorFac+0.0001);
//  cout << valorFac <<  " " <<variabilidadObj <<endl;

}
void MPP::init(){
   x_var.resize(N_OPT_DAY*nDias);
   for (int i = 0; i < nDias; i++){
	for (int j = 0; j < N_OPT_DAY; j++){
	      x_var[i*N_OPT_DAY+j] = MPP_problem->random_dish(j);
	   }
	}
 evaluate();
}
void MPP::restart(){
	for (int i = 0; i < nDias; i++){
		for (int j = 0; j < N_OPT_DAY; j++){
	      	  x_var[i*N_OPT_DAY+j] = MPP_problem->random_dish(j);
		}
	}
	evaluate();
	localSearch();
}

struct Food {
	int p[N_OPT_DAY];//check memmory usage, instead vector<int>..
//	bool operator<(const Food &i2) const {
//		return (make_pair(p1, make_pair(p2, p3)) < make_pair(i2.p1, make_pair(i2.p2, i2.p3)));
//	}
        bool operator<(const Food &i2) const {
		for(int i = 0; i < N_OPT_DAY; i++)
		{
		   if(p[i] < i2.p[i]) return true;
		   else if(p[i] > i2.p[i]) return false;
		   else return false;
		}
	}
}; 

void MPP::dependentCrossover(MPP &i2){
	if (crossoverType == PAIR_BASED_CROSSOVER){
		pairBasedCrossover(i2);
	} else if (crossoverType == UNIFORM_CROSSOVER){
		uniformCrossover(i2);
	} else if (crossoverType == UNIFORM2_CROSSOVER){
		uniform2Crossover(i2);
	}
	else
	{
	   cout << "Operador de cruce desconocido "<<endl;
	   exit(EXIT_FAILURE);
	}
}
void MPP::uniformCrossover(MPP &i2)
{
	for (int i = 0; i < nDias; i++){
	   if (rand() > (RAND_MAX / 2)){
		for(int j = 0; j < N_OPT_DAY; j++)
			swap(x_var[i*N_OPT_DAY+ j], i2.x_var[i*N_OPT_DAY + j]);
	    } 
	}
}

void MPP::uniform2Crossover(MPP &i2){

	for (int i = 0; i < (int)x_var.size(); i++){
		if (rand() > (RAND_MAX / 2)){
			swap(x_var[i], i2.x_var[i]);
		} 
	}
}

void MPP::pairBasedCrossover(MPP &i2)
{
	vector<Food> pendingI1, pendingI2;

	map<Food, int> f1;
	for (int i = 0; i < nDias; i++){
		Food f;
		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = x_var[i*N_OPT_DAY+j];
		f1[f]++;
	}

	for (int i = 0; i < nDias; i++){
		Food f;

		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = i2.x_var[i*N_OPT_DAY+j];
		if (f1.count(f)){//Comida en ambos
		   
 			for(int j = 0; j < N_OPT_DAY; j++)
			{
			  x_var[i*N_OPT_DAY + j] = f.p[j];
			  i2.x_var[i*N_OPT_DAY + j] = f.p[j];
			}
			f1[f]--;
			if (f1[f] == 0){
				f1.erase(f);
			}
		} else {
			pendingI2.push_back(f);
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
 		for(int j = 0; j < N_OPT_DAY; j++)
	        {
	           x_var[i*N_OPT_DAY + j] = f1.p[j];
	           i2.x_var[i*N_OPT_DAY + j] =f2.p[j];
		}
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
		for (int j = 0; j <N_OPT_DAY; j++){
			for (int k = 0; k < (int) MPP_problem->v_times_dishes[j].size(); k++){
				Neighbor n;
				n.variable = i * N_OPT_DAY + j;
				n.newValue = k;
				neighbors.push_back(n);
			}
		}
	}
	vector<int> bestIndividual = x_var;
	evaluate();
	//cout <<"entra..."<<valorFac<<endl;
	pair<double, double> bestResult = make_pair(valorFac, -variabilidadObj);
	for (int i = 0; i < 1000; i++){
		evaluate();
		pair<double, double> currentResult = make_pair(valorFac, -variabilidadObj);
		bool improved = true;
		while(improved){
			improved = false;
			random_shuffle(neighbors.begin(), neighbors.end());
			for (int i = 0; i < neighbors.size(); i++){
				int currentValue = x_var[neighbors[i].variable];
				x_var[neighbors[i].variable] = neighbors[i].newValue;
				evaluate();
				pair<double, double> newResult = make_pair(valorFac, -variabilidadObj); //note: variability is maximized..
				if (newResult >= currentResult){
					x_var[neighbors[i].variable] = currentValue;
				} else {
					improved = true;
					currentResult = newResult;
				}
			}
		}
		if (currentResult >= bestResult){
			x_var = bestIndividual;
		} else {
			bestResult = currentResult;
			bestIndividual = x_var;
			//cout << "--> "<< valorFac<<endl;
			cout << currentResult.first <<endl;
		}

		evaluate();
		if (badDays.empty()){
			int selectedDay = -1;
			if (heaviestNut != -1){
				vector< pair<double, int> > infoNut;
				for (int i = 0; i < nDias; i++){
					double total = 0.0;
					for(int k = 0; k < N_OPT_DAY; k++) total +=  MPP_problem->v_times_dishes[k][x_var[i*N_OPT_DAY + k]].v_nutrient_value[heaviestNut];
					infoNut.push_back(make_pair(total, i));
				}
				sort(infoNut.begin(), infoNut.end());
				if (heaviestType == 1) reverse(infoNut.begin(), infoNut.end());
				int nbestday = rand()%min(nDias ,  6);
				selectedDay = infoNut[nbestday].second;
				//cout << nDias<<endl;
			} else {
				selectedDay = rand() % nDias;
			}
//	cout << selectedDay<< " " << heaviestNut<<endl;
			for(int k = 0; k < N_OPT_DAY; k++) x_var[selectedDay*N_OPT_DAY + k] = MPP_problem->random_dish(k);

		} else {
			cout << "Dias: " << badDays.size() << endl;
			vector<int> v;
			for (auto it = badDays.begin(); it != badDays.end(); it++){
				v.push_back(*it);
			}
			random_shuffle(v.begin(), v.end());
			for (auto it = v.begin(); it != v.end(); it++){

				int day = *it;
				if (random()%2){
					for(int k = 0; k < N_OPT_DAY; k++) x_var[day*N_OPT_DAY + k] = MPP_problem->random_dish(k);
				} else {
					int which = rand() % N_OPT_DAY;
					x_var[day * N_OPT_DAY + which] = MPP_problem->random_dish(which);
				}
				//break;
				//cout << which << " "<<	x_var[day * N_OPT_DAY + which]<<endl;
			}
		}
	}
	x_var = bestIndividual;
	evaluate();
	cout <<"sale--- "<< valorFac<< " " <<variabilidadObj << endl;
	calculateFeasibilityDegree2();
	
	//cout << "Final " << valorFac << " " << precioObj << endl;
   exit(0);
}

/*

*/
int MPP::getDistance(MPP &ind2) {
	map<Food, int> f1;
	int dist = 0;
	for (int i = 0; i < nDias; i++){
		Food f;
		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = x_var[i*N_OPT_DAY+j];
		f1[f]++;
	}
	for (int i = 0; i < nDias; i++){
		Food f;
		for(int j = 0; j < N_OPT_DAY; j++) f.p[j] = ind2.x_var[i*N_OPT_DAY+j];
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
	for (int i = 0; i < x_var.size(); i++){
		os << x_var[i] << " ";
	}
	os << fitness <<endl;
//	os << valorFac << " " << precioObj << endl;
}

void MPP::exportcsv()
{
   ofstream ofs;
   ofs.open(MPP_problem->out_filename.c_str());
   ofs<<"DIA , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA \n";
    
   for(int i = 0; i < nDias; i++)
   {
	ofs << i+1;
	for(int j = 0; j < N_OPT_DAY; j++) ofs<< " , "<<MPP_problem->v_times_dishes[j][x_var[i*N_OPT_DAY + j]].description;
	ofs<<"\n";
   }
   ofs.close();
}

void MPP::full_search()
{
 vector< vector<int> > feasible_solutions;
 vector<pair<double, double > > fit_sol;
 fill(x_var.begin(), x_var.end(),0);
 long int cont = 0, max_perm =1;
   for(int max_opt = 0; max_opt < N_OPT_DAY; max_opt++)max_perm *= (int)MPP_problem->v_times_dishes[max_opt].size();
//  max_perm = pow(10, N_OPT_DAY);
 evaluate();
 pair< double, double> bestResult;
 bestResult.first = valorFac;
 vector<int> x_best = x_var;
 while(cont++ < max_perm)
 {

   evaluate();
   double current = valorFac;
   if( current < bestResult.first)
   {
	bestResult.first = current;
	bestResult.second = variabilidadObj;
	x_best = x_var;		
        cout << bestResult.first << " " <<bestResult.second<< " " <<cont<<endl;
   }
   else if( current <= 1e-9) //feasible solution
   {
	feasible_solutions.push_back(x_var);
        fit_sol.push_back(make_pair(valorFac, variabilidadObj));
   }

   x_var[0]++;
   for(int max_opt = 0; max_opt < N_OPT_DAY; max_opt++)
   {
      if(x_var[max_opt] >= MPP_problem->v_times_dishes[max_opt].size())
      {
	x_var[max_opt] = 0;
	if(max_opt+1 < N_OPT_DAY)
        x_var[max_opt+1]++;
      }
   }
 }
 if(feasible_solutions.empty()) 
  {
    feasible_solutions.push_back(x_best);
     fit_sol.push_back(bestResult);
  }
  
   cout << bestResult.first << " " <<bestResult.second<<endl;
   ofstream ofs;
   ofs.open(MPP_problem->out_filename.c_str());
   ofs<<"Factibilidad , Variabilidad , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA \n";
 ///exporting solutions....
 for(int i = 0; i < feasible_solutions.size(); i++)
 {
    ofs << fit_sol[i].first << " , "<< fit_sol[i].second<< " , ";
    for(int j = 0; j < feasible_solutions[i].size(); j++)
    {
	for(int j = 0; j < N_OPT_DAY; j++) ofs<< " , "<<MPP_problem->v_times_dishes[j][x_var[i*N_OPT_DAY + j]].description;
	ofs<<"\n";
    } 
    ofs.close();
 }
}
