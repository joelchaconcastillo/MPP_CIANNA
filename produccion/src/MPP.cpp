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
void MPP::calculateFeasibilityDegree(){
	valorFac = 0.0;
	int num_nutr = (int)MPP_problem->v_constraints.size();
	vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
	badDays.clear();

	double heaviestValue = 0;
	heaviestNut = -1;
        for(int a = 0; a < nopts; a++)
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
	         //for(int k = 0; k < N_OPT_DAY; k++)
		 for(int b = 0; b < N_OPT_DAY-2; b++)
		 {
		   int k = opts[a*(N_OPT_DAY-2)+b];
	   	   dayNutr += MPP_problem->v_times_dishes[k][x_var[j*N_OPT_DAY+k]].v_nutrient_value[i];
		 }
	   	globalNutr += dayNutr;
	   	if(v_constraints[i].type == DIARIA)
             	{
	               if(dayNutr  < minv) valorFac += ((middle - dayNutr)/middle)*((middle - dayNutr)/middle)*WEIGHT_DAY, badDays.insert(j);
	               else if (dayNutr > maxv) valorFac +=((dayNutr - middle)/middle)*((dayNutr - middle)/middle)*WEIGHT_DAY, badDays.insert(j);
           	}
	      }
	      if(v_constraints[i].type == GLOBAL)
              {
	        if(globalNutr < minv)
	        {
	          double v = ((middle - globalNutr)/middle)*((middle - globalNutr)/middle);
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
	          double v =((globalNutr - middle)/middle)*((globalNutr - middle)/middle);
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
	vector<int> bestIndividual = x_var;
	evaluate();
	cout <<"entra..."<<valorFac<<endl;
	pair<double, double> bestResult = make_pair(valorFac, -variabilidadObj);

	for (int i = 0; i < 1000; i++){
		pair<double, double> currentResult = First_Improvement_Hill_Climbing(neighbors, x_var);
		if (currentResult >= bestResult){
			x_var = bestIndividual;
		} else {
			bestResult = currentResult;
			bestIndividual = x_var;
	evaluate();	
			cout << currentResult.first << " " << currentResult.second <<" " << badDays.size()<<endl;
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
	x_var = bestIndividual;
	evaluate();
	cout <<"sale--- "<< valorFac<< " " <<variabilidadObj << endl;
	exportcsv();
	calculateFeasibilityDegree2();
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
   ofs<<"DIA , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA ";
   //ofs<<"DIA , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA ";
    for(auto i = MPP_problem->dic_nut_id.begin(); i !=  MPP_problem->dic_nut_id.end(); i++) ofs <<" , "<<i->first ;
	ofs<< "\n";
    
   for(int i = 0; i < nDias; i++)
   {
	ofs << i+1;
	for(int j = 0; j < N_OPT_DAY; j++)
	  ofs<< " , "<<MPP_problem->v_times_dishes[j][x_var[i*N_OPT_DAY + j]].description;

    	   for(auto ij = MPP_problem->dic_nut_id.begin(); ij !=  MPP_problem->dic_nut_id.end(); ij++)
	   {
		double sum_nut = 0.0;
		for(int j = 0; j < N_OPT_DAY; j++)
	          sum_nut +=MPP_problem->v_times_dishes[j][x_var[i*N_OPT_DAY + j]].v_nutrient_value[ij->second];
 		 ofs <<" , \" "<< sum_nut <<" ["<<MPP_problem->v_constraints[ij->second].min<<","<<MPP_problem->v_constraints[ij->second].max<<"] \"" ;
	   }
	ofs<<"\n";
   }
   ofs.close();
}

void MPP::full_search()
{
omp_set_num_threads(24);
 vector< vector<int> > feasible_solutions;
 vector<pair<double, double > > fit_sol;
 long int cont = 0, max_perm =1;
 //information to get each permutation of the options-feasible space..
 fill(x_var.begin(), x_var.end(),0);
 vector<int> v_max_opt;
 for(int max_opt = 0; max_opt < N_OPT_DAY; max_opt++)
 {
   int opt_s = (int)MPP_problem->v_times_dishes[max_opt].size();
   max_perm *= opt_s;
   v_max_opt.push_back(opt_s);
 }
 
 cout << max_perm<<endl;
 evaluate();
 pair< double, double> bestResult = make_pair(valorFac, -variabilidadObj);
 vector<int> x_best = x_var;
 vector<vector<infoDishes> > times = MPP_problem->v_times_dishes;

 int num_nutr = (int)MPP_problem->v_constraints.size();
 vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
 for(long cont = 0; cont < max_perm; cont++)
 {
       if(cont > 0 ) my_next_permutation(x_var, v_max_opt);
  // vector<vector<int>> parallel_x_var;
  // for(long i = 0; i < (long)1e8 && cont < max_perm; i++, cont++) 
  // {
  //      parallel_x_var.push_back(x_var);
  // }
   int day_constraints = (int)MPP_problem->v_constraint_day.size();
 //  #pragma omp parallel for shared(bestResult, x_best, feasible_solutions, fit_sol)
   //for(long i = 0; i < parallel_x_var.size(); i++)
   //{ 
       double valorFac = 0.0;
        for (int j = 0; j < day_constraints; j++)
        {
          double accum_nut = 0.0;
          int index = MPP_problem->v_constraint_day[j];
           for(unsigned int k = 0; k < N_OPT_DAY; k++)
                 accum_nut += MPP_problem->v_times_dishes[k][x_var[k]].v_nutrient_value[index];
           double minv = v_constraints[index].min, maxv = v_constraints[index].max;
           valorFac += (accum_nut < minv)?((minv - accum_nut)/minv)*((minv - accum_nut)/minv)*WEIGHT_DAY:0;
           valorFac += (accum_nut > maxv)?((accum_nut - maxv)/maxv)*((accum_nut - maxv)/maxv)*WEIGHT_DAY:0;
           if(valorFac > bestResult.first) break; //first optimization..
        }
      double current = valorFac;
      if( current < bestResult.first)
      {
//           evaluate();
           bestResult.first = current;
   //        bestResult.second = variabilidadObj;
           x_best = x_var;		
           cout << bestResult.first << " " <<bestResult.second<< " " <<cont<<endl;
      }
      if( current <= 1e-9) //feasible solution
      {
    //       evaluate();
           feasible_solutions.push_back(x_var);
           fit_sol.push_back(make_pair(valorFac, variabilidadObj));
      }
  // }
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
   ofs<<"Factibilidad , Variabilidad , DESAYUNO , COLACION_MATUTINA , COMIDA_ENTRADA , COMIDA_ENTRADA , COMIDA_PRINCIPAL , COMIDA_PRINCIPAL , COLACION_VESPERTINA , CENA \n";
 ///exporting solutions....
 for(int i = 0; i < feasible_solutions.size(); i++)
 {
    x_var=feasible_solutions[i];
    evaluate();
    ofs << this->valorFac << " , "<< this->variabilidadObj<< " , ";
    for(int j = 0; j < feasible_solutions[i].size(); j++)
    {
	for(int j = 0; j < N_OPT_DAY; j++) ofs<< " , "<<MPP_problem->v_times_dishes[j][x_var[i*N_OPT_DAY + j]].description;
	ofs<<"\n";
    } 
    ofs.close();
 }
}
double MPP::init_incremental_evaluation(vector<vector< double > > &globalPlan, vector< vector< vector<double> > > &nutriment_per_day, vector<int> &sol)
{ 
	int num_nutr = (int)MPP_problem->v_constraints.size();
	vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
        vector<int> &v_constraint_global = MPP_problem->v_constraint_global, &v_constraint_day = MPP_problem->v_constraint_day;
	globalPlan.clear();
	nutriment_per_day.clear();
        globalPlan.assign(nopts,vector<double> (num_nutr, 0.0 ));
	nutriment_per_day.assign( nopts, vector<vector<double> > (nDias, vector<double> (num_nutr, 0)));
        double unfeasibility_next  = 0.0;
        for(int a = 0; a < nopts; a++)
	{
 	  for(int j = 0; j < num_nutr; j++)
	  {
	   for(int i = 0; i < nDias; i++)
	   {
		 for(int b = 0; b < N_OPT_DAY-2; b++)
		 {
		   int k = opts[a*(N_OPT_DAY-2)+b];
	    	   nutriment_per_day[a][i][j] += MPP_problem->v_times_dishes[k][sol[i*N_OPT_DAY + k]].v_nutrient_value[j];
		 }
	        globalPlan[a][j] += nutriment_per_day[a][i][j]; 	
	        if( v_constraints[j].type == DIARIA)
                {
                   double minv = v_constraints[j].min;
                   double maxv = v_constraints[j].max;
	  	   double middle = (maxv+minv)*0.5;
	           double nut = nutriment_per_day[a][i][j];
	           if( nut < minv) unfeasibility_next+= ((middle - nut)/middle)*((middle - nut)/middle)*WEIGHT_DAY;
	           else if (nut > maxv) unfeasibility_next+=((nut - middle)/middle)*((nut - middle)/middle)*WEIGHT_DAY;
                }
	   }
	   if( v_constraints[j].type == GLOBAL)
           {
                   double minv = v_constraints[j].min;
                   double maxv = v_constraints[j].max;
	  	   double middle = (maxv+minv)*0.5;
	           double nut = globalPlan[a][j];
	           if( nut < minv) unfeasibility_next+= ((middle - nut)/middle)*((middle - nut)/middle);
	           else if (nut > maxv) unfeasibility_next+=((nut - middle)/middle)*((nut - middle)/middle);
            }
	  }
	}
   return unfeasibility_next;
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
  	cout  << time<<endl; 
   for(int a = 0; a < 2; a++)
   {	
    if(time == STARTER_1 && a != 0) continue;
    if(time == MAIN_COURSE_1 && a != 0) continue;
    if(time == STARTER_2 && a != 1) continue;
    if(time == MAIN_COURSE_2 && a != 1) continue;

    for(unsigned int j = 0; j < num_nutr; j++)
    {
       	//update sumatory of nutriments....
	double new_nut_value = (-v_times_dishes[time][current_sol[new_neighbor.variable]].v_nutrient_value[j] + v_times_dishes[time][new_neighbor.newValue].v_nutrient_value[j]);
       if(v_constraints[j].type == DIARIA)
       {
          double minv = v_constraints[j].min;
          double maxv = v_constraints[j].max;
	  double middle = (maxv+minv)*0.5;
	  double nut = nutriment_per_day[a][day][j] + new_nut_value;
  	  double original_nut = nutriment_per_day[a][day][j];	
	  if( nut < minv) new_partial_infeasibility+= ((middle - nut)/middle)*((middle - nut)/middle)*WEIGHT_DAY;
	  else if (nut > maxv) new_partial_infeasibility+=((nut - middle)/middle)*((nut - middle)/middle)*WEIGHT_DAY;
	  if( original_nut  < minv) original_partial_infeasibility += ((middle - original_nut)/middle)*((middle - original_nut)/middle)*WEIGHT_DAY;
	  else if (original_nut > maxv) original_partial_infeasibility +=((original_nut - middle)/middle)*((original_nut - middle)/middle)*WEIGHT_DAY;
       }
       else if(v_constraints[j].type == GLOBAL)
       { 
          double minv = v_constraints[j].min;
          double maxv = v_constraints[j].max;
	  double middle = (maxv+minv)*0.5;
          double nut = globalPlan[a][j] + new_nut_value;
          double original_nut = globalPlan[a][j];
          if(nut < minv) new_partial_infeasibility+= ((middle - nut)/(middle))*((middle - nut)/(middle));
          else if (nut > maxv) new_partial_infeasibility+=((nut - middle)/middle)*((nut - middle)/middle);
          if(original_nut < minv) original_partial_infeasibility += ((middle - original_nut)/(middle))*((middle- original_nut)/(middle));
          else if ( original_nut > maxv) original_partial_infeasibility +=((original_nut - middle)/middle)*((original_nut - middle)/middle);
       }
    }
   }
//    return  - original_partial_infeasibility + new_partial_infeasibility;
    return  current_infeasibility - original_partial_infeasibility + new_partial_infeasibility;
}
void MPP::update_data_incremental_eval(vector< vector<double> > &globalPlan, vector< vector<vector<double> > > &nutriment_per_day, vector<int> &current_sol, Neighbor &new_neighbor)
{
    int num_nutr = (int)MPP_problem->v_constraints.size();
    vector<vector<infoDishes> > &v_times_dishes = (MPP_problem->v_times_dishes);
    int day =  new_neighbor.variable/N_OPT_DAY;
    int time = new_neighbor.variable%N_OPT_DAY;
	cout <<time<<endl;
   for(int a = 0; a < nopts; a++)
   {
    if(time == STARTER_1 && a != 0) continue;
    if(time == MAIN_COURSE_1 && a != 0) continue;
    if(time == STARTER_2 && a != 1) continue;
    if(time == MAIN_COURSE_2 && a != 1) continue;
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
	int tmp2  = current_sol[neighbors[i].variable];
	current_sol[neighbors[i].variable] = neighbors[i].newValue;
//	double new_infeasibility = init_incremental_evaluation(globalPlan, nutriment_per_day, current_sol);
	evaluate();
    current_sol[neighbors[i].variable] = tmp2;

	double new_infeasibility =  inc_eval_feas_time(globalPlan, nutriment_per_day, current_sol, neighbors[i], current_infeasibility);
cout << new_infeasibility << " " <<valorFac<<endl;
cout << fabs(valorFac - new_infeasibility)<<endl;
	if(fabs(valorFac - new_infeasibility) > 0.5) exit(0);

	if( new_infeasibility < current_infeasibility)
//	if( (current_infeasibility - new_infeasibility) > EPSILON)
//	if(new_infeasibility < 0)
	{
	   improved = true;
	   update_data_incremental_eval(globalPlan, nutriment_per_day, current_sol, neighbors[i]);
	//   current_sol[neighbors[i].variable]=neighbors[i].newValue;
   //	   init_incremental_evaluation(globalPlan, nutriment_per_day, current_sol);
	   current_infeasibility = new_infeasibility;

	}
	else if(fabs(new_infeasibility - current_infeasibility) < EPSILON) //to check: epsilon...
//	else if( fabs(new_infeasibility) < EPSILON)
//	else if(new_infeasibility == current_infeasibility)
        { 
	  continue;
	    //double new_variability = calculateVariability(current_sol, neighbors[i]);
	    int tmp = current_sol[neighbors[i].variable];
	    current_sol[neighbors[i].variable] = neighbors[i].newValue;
	    double new_variability = calculateVariability(current_sol);
	    
	   if(current_variability < new_variability)
	   {	
	      improved = true;
	      current_variability = new_variability;
	   //   current_sol[neighbors[i].variable] = neighbors[i].newValue;
	   }
	   else 
	    current_sol[neighbors[i].variable] = tmp;
        }
      }
    }
	evaluate();
	return make_pair(valorFac ,-current_variability);
    return make_pair(current_infeasibility, -current_variability);
}
double MPP::calculateVariability()
{
}
double MPP::calculateVariability(vector<int> &current_sol)
{
 return 0;
   vector<vector<infoDishes> > &v_times_dishes = MPP_problem->v_times_dishes;
   double variability_day= 0.0, variability_global = 0.0;
       int min_diff = nDias;
    for(int k = 0; k < N_OPT_DAY; k++)
    {
       vector<int> cont_d(nDias, 0);
       vector<int> last_d((int)v_times_dishes[k].size(), -1);
       for(int i = 0 ; i < nDias; i++)
       {
	//Daily constraints.....
	   if(k == MORNING_SNACK && day_constraint(v_times_dishes[MORNING_SNACK][current_sol[i*N_OPT_DAY + MORNING_SNACK]], v_times_dishes[MORNING_SNACK][current_sol[i*N_OPT_DAY + EVENING_SNACK]]))
		variability_day++;
	   if(k == MAIN_COURSE_1 && day_constraint(v_times_dishes[MAIN_COURSE_1][current_sol[i*N_OPT_DAY + MAIN_COURSE_1]], v_times_dishes[MAIN_COURSE_2][current_sol[i*N_OPT_DAY + MAIN_COURSE_2]]))
		variability_day++;
	   if(k == STARTER_1 && day_constraint(v_times_dishes[STARTER_1][current_sol[i*N_OPT_DAY + STARTER_1]], v_times_dishes[STARTER_2][current_sol[i*N_OPT_DAY + STARTER_2]]))
		variability_day++;
	//PONER RESTRICCIÓN POR DÍA DE LA ENTRADAS Y DE LOS PLATILLOS FUERTE...!!!!!!!!!!!
	//if( k == MAIN_COURSE_1)
	{
	   int id_dish = current_sol[i*N_OPT_DAY + k];//v_times_dishes[k][current_sol[i*N_OPT_DAY + k]];
	   if(last_d[id_dish] != -1 )
	   {
	      int diff = i - last_d[id_dish];
	      cont_d[diff]++;
	      min_diff =  min(min_diff, diff);
	   }
	   last_d[id_dish] = i;
	}
       }
//		cout << min_diff<<endl;
//	if( k == MAIN_COURSE_1)
	{

	if(min_diff != nDias)
//	if( min_diff > 2) variability_global += 100;
//	variability_global += min_diff;//cont_d[min_diff];//cont_d[min_diff]*(min_diff+1);
	if( min_diff == 1) variability_global += 1.0/(cont_d[min_diff]+EPSILON);
	else variability_global +=min_diff;
	}
    }

// cout << variability_global << "____"<< endl;
 return 1e3*variability_day + variability_global;
}
bool MPP::day_constraint(infoDishes &dish1, infoDishes &dish2)
{
   if(dish1.category == CATEGORY_BOTH || dish2.category == CATEGORY_BOTH)
	if(dish1.description != dish2.description) return true;
   else if(dish1.category != dish2.category)
	if(dish1.description != dish2.description) return true;
   return false;
}
