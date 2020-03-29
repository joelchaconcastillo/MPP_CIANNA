#include <signal.h>
#include <omp.h>
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
        max_description_id = 0;
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
	       max_description_id = max(max_description_id, str_dish.description);
	       if(cell.empty()) break; //The file has an extra empty line
		
	       if(str_dish.time_day == "DESAYUNO") v_times_dishes[BREAKFAST].push_back(str_dish);
	       else if(str_dish.time_day == "COLACION_MATUTINA") v_times_dishes[MORNING_SNACK].push_back(str_dish);
	   //    else if(str_dish.time_day == "COMIDA_ENTRADA") v_times_dishes[STARTER].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_1) v_times_dishes[STARTER_1].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_2) v_times_dishes[STARTER_2].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_BOTH) //this dish blong to both categories...
		else if(str_dish.time_day == "COMIDA_ENTRADA" )
		{
		 v_times_dishes[STARTER_1].push_back(str_dish);
		 v_times_dishes[STARTER_2].push_back(str_dish);
		}
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL") v_times_dishes[MAIN_COURSE].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_1) v_times_dishes[MAIN_COURSE_1].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_2) v_times_dishes[MAIN_COURSE_2].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_BOTH) //both categories..
		else if(str_dish.time_day == "COMIDA_PRINCIPAL")
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
		if( trim(cell) == "GLOBAL") str_constraint_nutrient.type=GLOBAL;
		else if( trim(cell) == "DIARIA") str_constraint_nutrient.type=DIARIA;
               //str_constraint_nutrient.type = trim(cell);
	       
	       getline(ifs, cell, ',');
	       str_constraint_nutrient.name = trim(cell);
	       getline(ifs, cell, ',');
	       str_constraint_nutrient.min = stod(cell);
	       getline(ifs, cell, '\n'); //last word...
               str_constraint_nutrient.max = stod(cell);
	       dic_nut_id[str_constraint_nutrient.name] = (int) v_constraints.size();
	       if( str_constraint_nutrient.type == GLOBAL)
		{
		 str_constraint_nutrient.min *=nDias;
		 str_constraint_nutrient.max *=nDias;
		 v_constraint_global.push_back((int)v_constraints.size());
		}
	       else if( str_constraint_nutrient.type == DIARIA) v_constraint_day.push_back((int)v_constraints.size());
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
}

/*
 The objective is defined as follows:
  max f(x) + g(y)
  where 
     x is the varaibility by day 
     y is the global variability 
*/
void MPP::evaluate(){
  calculateFeasibilityDegree();
  variabilidadObj = calculateVariability(x_var);
  fitness =valorFac - variabilidadObj;
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

void MPP::localSearch( ) {
	vector<Neighbor> neighbors;
	for (int i = 0; i < nDias; i++)
        {
	  for (int j = 0; j <N_OPT_DAY; j++)
	  {
	     for (int k = 0; k < (int) MPP_problem->v_times_dishes[j].size(); k++)
	     {
		Neighbor n;
		n.variable = i * N_OPT_DAY + j;
		n.newValue = k;
		neighbors.push_back(n);
	     }
	   }
	 }

       vector<Neighbor_swap> neighbors_swap;
  for(int i = 0; i < nDias; i++)
     for(int j = i+1; j < nDias; j++) neighbors_swap.push_back({i, j});

	vector<int> bestIndividual = x_var;
	evaluate();
	cout <<"entra..."<<valorFac<<endl;
	pair<double, double> bestResult = make_pair(valorFac, -variabilidadObj);

	for (int i = 0; i < ITERATIONS_LS; i++){
		pair<double, double> currentResult = First_Improvement_Hill_Climbing(neighbors, x_var);
//		currentResult = First_Improvement_Hill_Climbing_swap(neighbors_swap, currentResult, x_var);
		
		if (currentResult >= bestResult){
			x_var = bestIndividual;
		} else {
			bestResult = currentResult;
			bestIndividual = x_var;
	evaluate();	
			cout << currentResult.first << " " << currentResult.second <<" " << badDays.size()<< " " << i<<endl;
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
			} else {
				selectedDay = rand() % nDias;
			}
//	cout << selectedDay<< " " << heaviestNut<<endl;
			for(int k = 0; k < N_OPT_DAY; k++) x_var[selectedDay*N_OPT_DAY + k] = MPP_problem->random_dish(k);

		} else {
	//		cout << "Dias: " << badDays.size() << endl;
			vector<int> v;
			for (auto it = badDays.begin(); it != badDays.end(); it++){
				v.push_back(*it);
			}
			random_shuffle(v.begin(), v.end());
			for (auto it = v.begin(); it != v.end(); it++){

				int day = *it;
		//		if (random()%2){
		//			for(int k = 0; k < N_OPT_DAY; k++) x_var[day*N_OPT_DAY + k] = MPP_problem->random_dish(k);
		//		} else {
					int which = rand() % N_OPT_DAY;
					x_var[day * N_OPT_DAY + which] = MPP_problem->random_dish(which);
		//		}
				break;
				//cout << which << " "<<	x_var[day * N_OPT_DAY + which]<<endl;
			}
		}
	}
 
	//First_Improvement_Hill_Climbing_swap(neighbors_swap, bestResult, bestIndividual);
	x_var = bestIndividual;
	evaluate();
	cout <<"sale--- "<< valorFac<< " " <<variabilidadObj << endl;
	exportcsv();
//	calculateFeasibilityDegree2();
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
   ofs << "DIA ";
   set<int> times_selected_per_day;
   for(int i = 0; i < MPP_problem->conf_day.size(); i++)
     for(auto t =MPP_problem->conf_day[i].begin(); t != MPP_problem->conf_day[i].end(); t++) times_selected_per_day.insert(*t);
   for(auto a = times_selected_per_day.begin(); a != times_selected_per_day.end(); a++)
   {
	   if(*a == BREAKFAST) ofs << " , DESAYUNO ";
	   if(*a == MORNING_SNACK) ofs << " , COLACION_MATUTINA ";	
	   if(*a == STARTER_1) ofs << " , COMIDA_ENTRADA ";	
	   if(*a == STARTER_2) ofs << " , COMIDA_ENTRADA ";	
	   if(*a == MAIN_COURSE_1) ofs << " , COMIDA_PRINCIPAL ";	
	   if(*a == MAIN_COURSE_2) ofs << " , COMIDA_PRINCIPAL ";	
	   if(*a == EVENING_SNACK) ofs << " , COLACION_VESPERTINA ";	
	   if(*a == DINNER) ofs << " , CENA ";	
   }
//   ofs<<"DIA , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA ";
   //ofs<<"DIA , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA ";
    for(auto i = MPP_problem->dic_nut_id.begin(); i !=  MPP_problem->dic_nut_id.end(); i++) ofs <<" , "<<i->first ;
	ofs<< "\n";
    
   for(int i = 0; i < nDias; i++)
   {
	ofs << i+1;
   	for(auto a = times_selected_per_day.begin(); a != times_selected_per_day.end(); a++)
	  ofs<< " , "<<MPP_problem->v_times_dishes[(*a)][x_var[i*N_OPT_DAY + (*a)]].description;
    	   for(auto ij = MPP_problem->dic_nut_id.begin(); ij !=  MPP_problem->dic_nut_id.end(); ij++)
	   {
		ofs <<" , \" (";
		for(int c = 0; c < MPP_problem->conf_day.size(); c++)
		{
		double sum_nut = 0.0;
     		  for(auto t =MPP_problem->conf_day[c].begin(); t != MPP_problem->conf_day[c].end(); t++) 
	              sum_nut +=MPP_problem->v_times_dishes[*t][x_var[i*N_OPT_DAY + (*t)]].v_nutrient_value[ij->second];
			if( c>0) ofs<<",";
			ofs <<sum_nut ;
		}
 		 ofs<<") ["<<MPP_problem->v_constraints[ij->second].min<<","<<MPP_problem->v_constraints[ij->second].max<<"] \"" ;
	   }
	ofs<<"\n";
   }
   ofs.close();
}

void MPP::full_search()
{
 vector<vector<infoDishes> > times = MPP_problem->v_times_dishes;
 vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
 vector< vector<int> > feasible_solutions;
 vector<pair<double, double > > fit_sol;
 long int cont = 0, max_perm =1;
 //information to get each permutation of the options-feasible space..
 fill(x_var.begin(), x_var.end(),0);
 vector<int> v_max_opt;
 for(int max_opt = 0; max_opt < MPP_problem->time_conf.size(); max_opt++)
 {
   if( MPP_problem->time_conf[max_opt].empty()) continue;
   int opt_s = (int)MPP_problem->v_times_dishes[max_opt].size();
   max_perm *= opt_s;
   v_max_opt.push_back(opt_s);
 }
 cout << max_perm<<endl;
 evaluate();
 pair< double, double> bestResult = make_pair(valorFac, -variabilidadObj);
 vector<int> x_best = x_var;
 int num_nutr = (int)MPP_problem->v_constraints.size();
 int day_constraints = (int)MPP_problem->v_constraint_day.size();
 for(long cont = 0; cont < max_perm; cont++)
 {
   if(cont > 0 ) my_next_permutation(x_var, v_max_opt);
   double in_valorFac = 0.0;
   for (int j = 0; j < day_constraints; j++)
   {
      double accum_nut = 0.0;
      int index = MPP_problem->v_constraint_day[j];
      for(unsigned int k = 0; k < MPP_problem->time_conf.size(); k++)
      {
         if(MPP_problem->time_conf[k].empty()) continue;
         accum_nut += times[k][x_var[k]].v_nutrient_value[index];
      }	   
      double minv = v_constraints[index].min;
      double maxv = v_constraints[index].max;
      double middle = (maxv+minv)*0.5;
           in_valorFac += (accum_nut < minv)?((minv - accum_nut)/middle)*((minv - accum_nut)/middle)*WEIGHT_DAY:0;
           in_valorFac += (accum_nut > maxv)?((accum_nut - maxv)/middle)*((accum_nut - maxv)/middle)*WEIGHT_DAY:0;
           if(in_valorFac > bestResult.first) break; //first optimization..
   }
   if( in_valorFac <= bestResult.first)
   {
      bestResult.first = in_valorFac;
      variabilidadObj = calculateVariability(x_var);
      x_best = x_var;		
      cout << bestResult.first << " " <<-variabilidadObj<< " " <<cont<<endl;
   }
   if( in_valorFac == 0.0) //feasible solution
   {
      variabilidadObj = calculateVariability(x_var);
      feasible_solutions.push_back(x_var);
      fit_sol.push_back(make_pair(in_valorFac, -variabilidadObj));
   }
   if( (cont % (long)1e8 )== 0)
   {
      cout << "========\n " << (double)cont/(double)max_perm<<endl;
      cout << bestResult.first << " " <<bestResult.second<< " " <<cont<<endl;
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
   ofs << "DIA ";
   set<int> times_selected_per_day;
   for(int i = 0; i < MPP_problem->conf_day.size(); i++)
     for(auto t =MPP_problem->conf_day[i].begin(); t != MPP_problem->conf_day[i].end(); t++) times_selected_per_day.insert(*t);
   for(auto a = times_selected_per_day.begin(); a != times_selected_per_day.end(); a++)
   {
	   if(*a == BREAKFAST) ofs << " , DESAYUNO ";
	   if(*a == MORNING_SNACK) ofs << " , COLACION_MATUTINA ";	
	   if(*a == STARTER_1) ofs << " , COMIDA_ENTRADA ";	
	   if(*a == STARTER_2) ofs << " , COMIDA_ENTRADA ";	
	   if(*a == MAIN_COURSE_1) ofs << " , COMIDA_PRINCIPAL ";	
	   if(*a == MAIN_COURSE_2) ofs << " , COMIDA_PRINCIPAL ";	
	   if(*a == EVENING_SNACK) ofs << " , COLACION_VESPERTINA ";	
	   if(*a == DINNER) ofs << " , CENA ";	
   }
//   ofs<<"DIA , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA ";
   //ofs<<"DIA , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA ";
    for(auto i = MPP_problem->dic_nut_id.begin(); i !=  MPP_problem->dic_nut_id.end(); i++) ofs <<" , "<<i->first ;
	ofs<< "\n";
  for(int j = 0; j <   feasible_solutions.size(); j++)
  {
   x_var = feasible_solutions[j];
   for(int i = 0; i < nDias; i++)
   {
	ofs << fit_sol[j].first<<"-"<<fit_sol[j].second;
   	for(auto a = times_selected_per_day.begin(); a != times_selected_per_day.end(); a++)
	  ofs<< " , "<<MPP_problem->v_times_dishes[(*a)][x_var[i*N_OPT_DAY + (*a)]].description;
    	   for(auto ij = MPP_problem->dic_nut_id.begin(); ij !=  MPP_problem->dic_nut_id.end(); ij++)
	   {
		ofs <<" , \" (";
		for(int c = 0; c < MPP_problem->conf_day.size(); c++)
		{
		double sum_nut = 0.0;
     		  for(auto t =MPP_problem->conf_day[c].begin(); t != MPP_problem->conf_day[c].end(); t++) 
	              sum_nut +=MPP_problem->v_times_dishes[*t][x_var[i*N_OPT_DAY + (*t)]].v_nutrient_value[ij->second];
			if( c>0) ofs<<",";
			ofs <<sum_nut ;
		}
 		 ofs<<") ["<<MPP_problem->v_constraints[ij->second].min<<","<<MPP_problem->v_constraints[ij->second].max<<"] \"" ;
	   }
	ofs<<"\n";
   }
  }
   ofs.close();
}
double MPP::init_incremental_evaluation(vector<vector< double > > &globalPlan, vector< vector< vector<double> > > &nutriment_per_day, vector<int> &sol)
{ 
	int num_nutr = (int)MPP_problem->v_constraints.size();
	vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
        vector<int> &v_constraint_global = MPP_problem->v_constraint_global, &v_constraint_day = MPP_problem->v_constraint_day;
	globalPlan.clear();
	nutriment_per_day.clear();
        globalPlan.assign((int)MPP_problem->conf_day.size(),vector<double> (num_nutr, 0.0 ));
	nutriment_per_day.assign( (int)MPP_problem->conf_day.size(), vector<vector<double> > (nDias, vector<double> (num_nutr, 0)));
        double unfeasibility_next  = 0.0;
        for(int a = 0; a < MPP_problem->conf_day.size(); a++)
	{
 	  for(int j = 0; j < num_nutr; j++)
	  {
	   for(int i = 0; i < nDias; i++)
	   {
		 for(int b = 0; b < MPP_problem->conf_day[a].size(); b++)
		 {
		   int k = MPP_problem->conf_day[a][b];
	    	   nutriment_per_day[a][i][j] += MPP_problem->v_times_dishes[k][sol[i*N_OPT_DAY + k]].v_nutrient_value[j];
		 }
	        globalPlan[a][j] += nutriment_per_day[a][i][j]; 	
	        if( v_constraints[j].type == DIARIA)
                {
                   double minv = v_constraints[j].min;
                   double maxv = v_constraints[j].max;
	  	   double middle = (maxv+minv)*0.5;
	           double nut = nutriment_per_day[a][i][j];
	           if( nut < minv) unfeasibility_next += ((minv - nut)/middle)*((minv - nut)/middle)*WEIGHT_DAY;
	           else if (nut > maxv) unfeasibility_next +=((nut - maxv)/middle)*((nut - maxv)/middle)*WEIGHT_DAY;
                }
	   }
	   if( v_constraints[j].type == GLOBAL)
           {
                   double minv = v_constraints[j].min;
                   double maxv = v_constraints[j].max;
	  	   double middle = (maxv+minv)*0.5;
	           double nut = globalPlan[a][j];
	           if( nut < minv) unfeasibility_next += ((middle - nut)/middle)*((middle - nut)/middle);
	           else if (nut > maxv) unfeasibility_next += ((nut - middle)/middle)*((nut - middle)/middle);
            }
	  }
	}
   return unfeasibility_next;
}
void MPP::calculateFeasibilityDegree(){
	valorFac = 0.0;
	int num_nutr = (int)MPP_problem->v_constraints.size();
	vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
	badDays.clear();

	double heaviestValue = 0;
	heaviestNut = -1;
        for(int a = 0; a < MPP_problem->conf_day.size(); a++)
        {
           for(int i = 0; i < num_nutr; i++)
           {
              double minv = v_constraints[i].min;
              double maxv = v_constraints[i].max;
	      double middle = (maxv+minv)*0.5;
              double globalNutr = 0.0;
	      for(int j = 0; j < nDias; j++)
	      {
	         double dayNutr = 0.0;
		 for(int b = 0; b < MPP_problem->conf_day[a].size(); b++)
		 {
		   int k = MPP_problem->conf_day[a][b];
	   	   dayNutr += MPP_problem->v_times_dishes[k][x_var[j*N_OPT_DAY+k]].v_nutrient_value[i];
		 }
	   	globalNutr += dayNutr;
	   	if(v_constraints[i].type == DIARIA)
             	{
	               if(dayNutr  < minv) valorFac += ((minv - dayNutr)/middle)*((minv - dayNutr)/middle)*WEIGHT_DAY, badDays.insert(j);
	               else if (dayNutr > maxv) valorFac +=((dayNutr - maxv)/middle)*((dayNutr - maxv)/middle)*WEIGHT_DAY, badDays.insert(j);
           	}
	      }
	      if(v_constraints[i].type == GLOBAL)
              {
	        if(globalNutr < minv)
	        {
	          double v = ((minv - globalNutr)/middle)*((minv - globalNutr)/middle);
	          valorFac += v;
	          if( v >  heaviestValue)
	          {
	            heaviestValue = v;
	            heaviestNut = i;
	            heaviestType = -1;
	          }
	        }
	        else if (globalNutr > maxv)
	        {
	          double v =((globalNutr - maxv)/middle)*((globalNutr - maxv)/middle);
	          valorFac +=v;
	          if( v >  heaviestValue)
	          {
	            heaviestValue = v;
	            heaviestNut = i;
	            heaviestType = 1;
	          }
	        }
              }
           }
	}
}
double MPP::inc_eval_feas_time( vector< vector<double> > &globalPlan, vector< vector<vector<double> > > &nutriment_per_day, vector<int> &current_sol, Neighbor &new_neighbor, double current_infeasibility)
{
    int num_nutr = (int)MPP_problem->v_constraints.size();
    vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
    vector<int> &v_constraint_global = MPP_problem->v_constraint_global, &v_constraint_day = MPP_problem->v_constraint_day;
    vector<vector<infoDishes> > &v_times_dishes = (MPP_problem->v_times_dishes);
    int day =  new_neighbor.variable/N_OPT_DAY;
    int time = new_neighbor.variable%N_OPT_DAY;
    double new_partial_infeasibility = 0.0, original_partial_infeasibility = 0.0 ;
   for(int b = 0; b < MPP_problem->time_conf[time].size(); b++)
   {	
	int a = MPP_problem->time_conf[time][b];
    for(unsigned int j = 0; j < num_nutr; j++)
    {
       	//update sumatory of nutriments....
	double new_nut_value = (-v_times_dishes[time][current_sol[new_neighbor.variable]].v_nutrient_value[j] + v_times_dishes[time][new_neighbor.newValue].v_nutrient_value[j]);
	new_nut_value = (fabs(new_nut_value)>EPSILON)?new_nut_value:0.0;
       if(v_constraints[j].type == DIARIA)
       {
          double minv = v_constraints[j].min;
          double maxv = v_constraints[j].max;
	  double middle = (maxv+minv)*0.5;
	  double nut = nutriment_per_day[a][day][j] + new_nut_value;
  	  double original_nut = nutriment_per_day[a][day][j];	
	  if( nut < minv)new_partial_infeasibility+= ((minv- nut)/middle)*((minv - nut)/middle)*WEIGHT_DAY;
	  else if (nut > maxv) new_partial_infeasibility+=((nut - maxv)/middle)*((nut - maxv)/middle)*WEIGHT_DAY;
	  if( original_nut  < minv)original_partial_infeasibility += ((minv - original_nut)/middle)*((minv - original_nut)/middle)*WEIGHT_DAY;
	  else if (original_nut > maxv) original_partial_infeasibility +=((original_nut - maxv)/middle)*((original_nut - maxv)/middle)*WEIGHT_DAY;
       }
       else if(v_constraints[j].type == GLOBAL)
       { 
          double minv = v_constraints[j].min;
          double maxv = v_constraints[j].max;
	  double middle = (maxv+minv)*0.5;
          double nut = globalPlan[a][j] + new_nut_value;
          double original_nut = globalPlan[a][j];
          if(nut < minv) new_partial_infeasibility+= ((minv - nut)/(middle))*((minv - nut)/(middle));
          else if (nut > maxv) new_partial_infeasibility+=((nut - maxv)/middle)*((nut - maxv)/middle);
          if(original_nut < minv) original_partial_infeasibility += ((minv - original_nut)/(middle))*((minv- original_nut)/(middle));
          else if ( original_nut > maxv) original_partial_infeasibility +=((original_nut - maxv)/middle)*((original_nut - maxv)/middle);
       }
     }
   }
    return  current_infeasibility - original_partial_infeasibility + new_partial_infeasibility;
}
void MPP::update_data_incremental_eval(vector< vector<double> > &globalPlan, vector< vector<vector<double> > > &nutriment_per_day, vector<int> &current_sol, Neighbor &new_neighbor)
{
    int num_nutr = (int)MPP_problem->v_constraints.size();
    vector<vector<infoDishes> > &v_times_dishes = (MPP_problem->v_times_dishes);
    int day =  new_neighbor.variable/N_OPT_DAY;
    int time = new_neighbor.variable%N_OPT_DAY;
   for(int b = 0; b < MPP_problem->time_conf[time].size(); b++)
   {	
    int a = MPP_problem->time_conf[time][b];
    for(unsigned int j = 0; j < num_nutr; j++)//check new neighbor..
    {
       	//update sumatory of nutriments....
	double diff = (-v_times_dishes[time][current_sol[new_neighbor.variable]].v_nutrient_value[j] + v_times_dishes[time][new_neighbor.newValue].v_nutrient_value[j]);
        nutriment_per_day[a][day][j] += diff;
	globalPlan[a][j] += diff;
    }
   }
    current_sol[new_neighbor.variable]=new_neighbor.newValue;
}
void MPP::my_next_permutation(vector<int> &perm, vector<int> &v_max_opt)
{
      x_var[0]++;
      for(int max_opt = 0; max_opt < (int)v_max_opt.size(); max_opt++)
      {
         if(x_var[max_opt] >= v_max_opt[max_opt])
         {
           x_var[max_opt] = 0;
           if(max_opt+1 < (int)v_max_opt.size())
           x_var[max_opt+1]++;
         }
      }
}
pair<double, double> MPP::First_Improvement_Hill_Climbing(vector<Neighbor> &neighbors, vector<int> &current_sol)
{
   vector< vector<double> > globalPlan;
   vector< vector< vector<double> > > nutriment_per_day;

   double current_infeasibility = init_incremental_evaluation(globalPlan, nutriment_per_day, current_sol);
   evaluate();
   double current_variability = calculateVariability(current_sol);
   bool improved = true;
   while(improved)
   {
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for (int i = 0; i < neighbors.size(); i++)
     {
        //incremental evaluation...
	double new_infeasibility =  inc_eval_feas_time(globalPlan, nutriment_per_day, current_sol, neighbors[i], current_infeasibility);

	if( new_infeasibility < current_infeasibility)
	{
	   improved = true;
	   update_data_incremental_eval(globalPlan, nutriment_per_day, current_sol, neighbors[i]);
	   current_infeasibility = new_infeasibility;
//	   evaluate();
	//   if(fabs(new_infeasibility - valorFac) >0.3) exit(0);
	}
	else if(fabs(new_infeasibility - current_infeasibility) < EPSILON) //to check: epsilon...
        { 
	    int tmp = current_sol[neighbors[i].variable];
	    current_sol[neighbors[i].variable] = neighbors[i].newValue;
	    double new_variability = calculateVariability(current_sol);
	    
	   if(current_variability < new_variability)
	   {	
	      improved = true;
	      current_variability = new_variability;
	   }
	   else 
	    current_sol[neighbors[i].variable] = tmp;
        }
      }
    }
	evaluate();
    return make_pair(valorFac ,-variabilidadObj);//this avoid a numerical error provoked by the incremental evaluation..
    //return make_pair(current_infeasibility, -current_variability);
}
double MPP::calculateVariability()
{
}
double MPP::calculateVariability(vector<int> &current_sol)
{
   vector<vector<infoDishes> > &v_times_dishes = MPP_problem->v_times_dishes;
   double variability_day = 0.0, variability_global = 0.0, variability_cat_day=0.0, variability_fav=0.0;
   int max_variability_day =0;
   unordered_set<int> gl_ids;
   unordered_map<int, int> min_dist_day_id, last_day;
   unordered_map<int, bool> favorite;
   
   for(int d = 0; d < nDias; d++)
   {
        unordered_set<int> day_ids, snacks_cat, starter_cat, main_course_cat;
	for(int i = 0; i < MPP_problem->time_conf.size(); i++)
	{
	   if(MPP_problem->time_conf[i].empty())continue;
	   max_variability_day++;
	   int id = v_times_dishes[i][current_sol[d*N_OPT_DAY + i]].description;
	   favorite[id] =v_times_dishes[i][current_sol[d*N_OPT_DAY + i]].favorite;
	   day_ids.insert(id);
	   gl_ids.insert(id);
	   if( last_day.find(id) != last_day.end() && i == MAIN_COURSE_1)
	      min_dist_day_id[id] = min(min_dist_day_id[id], last_day[id]);
	   else min_dist_day_id[id] = nDias;
	      last_day[id] = d;

	   //cat by snacks..
	   if( i == MORNING_SNACK || i == EVENING_SNACK) snacks_cat.insert(id);
        //cat by starter..
	   else if( i == STARTER_1 || i == STARTER_2) starter_cat.insert(id);
	//cat by main_course....
	   else if( i == MAIN_COURSE_1|| i == MAIN_COURSE_2) main_course_cat.insert(id);
	}
	variability_cat_day += snacks_cat.size() + starter_cat.size() + main_course_cat.size();
	variability_day +=day_ids.size();
   }
   variability_global = gl_ids.size();
   for(auto i = min_dist_day_id.begin(); i != min_dist_day_id.end(); i++)
    {
	if(favorite[i->first] && i->second > DAYS_FAVORITE )
	  variability_fav += nDias;
	else if(!favorite[i->first] && i->second > DAYS_NO_FAVORITE )
	  variability_fav += nDias;
	else variability_fav += i->second;
    }
   variability_day /= (max_variability_day);
   variability_global /=(MPP_problem->max_description_id);
   variability_fav /= (nDias*MPP_problem->max_description_id);
   variability_cat_day /=(nDias*6);
   return variability_day*10.0 + variability_global*0.5 + variability_fav*5.0 + variability_cat_day*0.5;
}
bool MPP::day_constraint(infoDishes &dish1, infoDishes &dish2)
{
   if(dish1.category == CATEGORY_BOTH || dish2.category == CATEGORY_BOTH)
	if(dish1.description != dish2.description) return true;
   else if(dish1.category != dish2.category)
	if(dish1.description != dish2.description) return true;
   return false;
}
pair<double, double> MPP::First_Improvement_Hill_Climbing_swap(vector<Neighbor_swap> &neighbors, pair<double, double> &bestResult, vector<int> &bestIndividual)
{
  bool improved= true;
  vector<int> current = bestIndividual;
  pair<double, double> currentResult = bestResult;
	cout <<"entra... " << currentResult.second<<endl;
  while(improved)
  {
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for(int i = 0; i < neighbors.size(); i++)
     {
	for(int ii = 0; ii < MPP_problem->time_conf.size(); ii++)
	{
	   if(MPP_problem->time_conf[ii].empty())continue;
	   swap(current[neighbors[i].day1*N_OPT_DAY + ii], current[neighbors[i].day2*N_OPT_DAY + ii] );
	}
	currentResult.second = -calculateVariability(current);
	if (currentResult >= bestResult)
        {
	   for(int ii = 0; ii < MPP_problem->time_conf.size(); ii++)
	   {
	      if(MPP_problem->time_conf[ii].empty())continue;
	      swap(current[neighbors[i].day1*N_OPT_DAY + ii], current[neighbors[i].day2*N_OPT_DAY + ii] );
	   }
//			current = bestIndividual;
	}
	else 
	{
            improved = true;
	    bestResult = currentResult;
	    for(int ii = 0; ii < MPP_problem->time_conf.size(); ii++)
	    {
	      if(MPP_problem->time_conf[ii].empty())continue;
	      swap(bestIndividual[neighbors[i].day1*N_OPT_DAY + ii], bestIndividual[neighbors[i].day2*N_OPT_DAY + ii] );
	   }
	}
     }
  }
	cout <<"sale... " << bestResult.second<<endl;
  return bestResult;
}
