#include <chrono> 
#include "MPP.h"
#include "utils.h"
using namespace std;
using namespace std::chrono;

const string WHITESPACE = " \n\r\t\f\v";
ostream & operator << (ostream &out, const vector<double> &data)
{
	for(int i = 0; i < data.size(); i++) out << data[i]<< " ";
  return out;
}
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
    v_opt_dishes.resize(N_OPT_DAY); //N_times_dishes...
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
		 	{
			   double value = stod(cell);
//			   value = floor(value * 10 + 0.5)/10;
		           str_dish.v_nutrient_value[dic_nut_id[column_names[i]]] = value; //the nutrient values need to be stored in the same order that the contraints..
			}
			//else { cout << column_names[i]<<endl;}// cout << "error interno"<<endl; exit(EXIT_FAILURE);}
			}
               }
	       max_description_id = max(max_description_id, str_dish.description);
	       if(cell.empty()) break; //The file has an extra empty line
		
	       if(str_dish.time_day == "DESAYUNO") v_opt_dishes[BREAKFAST].push_back(str_dish);
	       else if(str_dish.time_day == "COLACION_MATUTINA") v_opt_dishes[MORNING_SNACK].push_back(str_dish);
	   //    else if(str_dish.time_day == "COMIDA_ENTRADA") v_opt_dishes[STARTER].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_1) v_opt_dishes[STARTER_1].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_2) v_opt_dishes[STARTER_2].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_ENTRADA" && str_dish.category == CATEGORY_BOTH) //this dish blong to both categories...
		else if(str_dish.time_day == "COMIDA_ENTRADA" )
		{
		 v_opt_dishes[STARTER_1].push_back(str_dish);
		 v_opt_dishes[STARTER_2].push_back(str_dish);
		}
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL") v_opt_dishes[MAIN_COURSE].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_1) v_opt_dishes[MAIN_COURSE_1].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_2) v_opt_dishes[MAIN_COURSE_2].push_back(str_dish);
//	       else if(str_dish.time_day == "COMIDA_PRINCIPAL" && str_dish.category == CATEGORY_BOTH) //both categories..
		else if(str_dish.time_day == "COMIDA_PRINCIPAL")
		{
		  v_opt_dishes[MAIN_COURSE_1].push_back(str_dish);
		  v_opt_dishes[MAIN_COURSE_2].push_back(str_dish);
		}
	       else if(str_dish.time_day == "COLACION_VESPERTINA") v_opt_dishes[EVENING_SNACK].push_back(str_dish);
	       else if(str_dish.time_day == "CENA") v_opt_dishes[DINNER].push_back(str_dish);
	       else if(str_dish.time_day == "COLACION_AMBAS")
	       {
 		 v_opt_dishes[MORNING_SNACK].push_back(str_dish);
 		 v_opt_dishes[EVENING_SNACK].push_back(str_dish);
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
       for(int i = 0; i < N_OPT_DAY; i++) cout << v_opt_dishes[i].size() << " ";
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
void MPP::evaluate(){
  calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
  calculateVariability(x_var, obj_values);
  fitness = obj_values[0]+obj_values[1];
  //variability tchebycheff approach
   double ref = nDias*10;
   double max1 = -ref, sum1=0.0;
  for(int a = 2; a < obj_values.size(); a++)
  {
     double v1 = obj_values[a];
	sum1 += fabs(v1-ref);
	max1 = max(max1,fabs(v1-(ref))/MPP_problem->weights[a-2]);
  }
  max1 += 0.0001*sum1; 
  fitness = 1e4*(obj_values[0]+obj_values[1]) + max1;
}
void MPP::evaluate(vector<int> &sol, vector<double> &objs){
  calculateFeasibilityDegree(sol, objs[0], objs[1]);
  calculateVariability(sol, objs);
}
void MPP::init(){
   x_var.resize(N_OPT_DAY*nDias);
   obj_values.resize(N_TIMES+2, 0);
   for (int i = 0; i < nDias; i++) perturb_day(x_var, i);
   evaluate();
}
void MPP::restart(){
	for (int i = 0; i < nDias; i++) perturb_day(x_var, i);
	evaluate();
	localSearch();
}
struct Food {
	int p[N_OPT_DAY];//check memmory usage, instead vector<int>..
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
   for (int i = 0; i < nDias; i++)
   {
      if (rand() > (RAND_MAX / 2))
	for(int j = 0; j < N_OPT_DAY; j++)
	   swap(x_var[i*N_OPT_DAY+ j], i2.x_var[i*N_OPT_DAY + j]);
   }
}

void MPP::uniform2Crossover(MPP &i2){
   for (int i = 0; i < (int)x_var.size(); i++)
   {
	if (rand() > (RAND_MAX / 2))
		swap(x_var[i], i2.x_var[i]);
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
	   for (int k = 0; k < (int) MPP_problem->v_opt_dishes[j].size(); k++)
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
     {
      for(int j = i+1; j < nDias; j++)
      {
         neighbors_swap.push_back({i, j});
      }
     }

     calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
     calculateVariability(x_var, obj_values);

     vector<int> bestIndividual = x_var;
     vector<double> best_objs = obj_values;
     //load incremental evaluation values...
     cout << "entra--- " << obj_values<<endl;
     auto start = high_resolution_clock::now(); 
     for (int i = 0; i < ITERATIONS_LS; i++)
     {
        calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
        calculateVariability(x_var, obj_values);
        First_Improvement_Hill_Climbing(neighbors, x_var, obj_values);
        First_Improvement_Hill_Climbing_swap(neighbors_swap, x_var, obj_values);
        if(comp_objs(obj_values, best_objs))
        {
	  best_objs = obj_values;
          bestIndividual = x_var;
	  cout << obj_values<< " "<< i<< " -- " << badDaysFeas.size() << endl;
        }
        else x_var = bestIndividual;

        calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
        calculateVariability(x_var, obj_values);
        if(badDaysFeas.empty())
        {
	   int selectedDay = -1;
	   if (heaviestNut != -1)
           {
		vector< pair<double, int> > infoNut;
		for (int i = 0; i < nDias; i++)
		{
		   double total = 0.0;
		   for(int k = 0; k < N_OPT_DAY; k++) total +=  MPP_problem->v_opt_dishes[k][x_var[i*N_OPT_DAY + k]].v_nutrient_value[heaviestNut];
		   infoNut.push_back(make_pair(total, i));
		}
		sort(infoNut.begin(), infoNut.end());
		if (heaviestType == 1) reverse(infoNut.begin(), infoNut.end());
		int nbestday = rand()%min(nDias ,  6);
		selectedDay = infoNut[nbestday].second;
	   } 
	   else 
	   {
		if( rand()%2 && !badDaysVar.empty())
		{
	 	   vector<int> tmp;
	           for(auto it = badDaysVar.begin(); it != badDaysVar.end(); it++) tmp.push_back(*it);
		   random_shuffle(tmp.begin(), tmp.end());
		   selectedDay = tmp[0];	
		}
		else
		 selectedDay = rand() % nDias;
	   }
	   perturb_day(x_var, selectedDay);
	   oneDaylocalSearch(x_var, selectedDay);
	 } 
	 else 
	 {
	   for (auto it = badDaysFeas.begin(); it != badDaysFeas.end(); it++)
	   {
	      perturb_day(x_var, *it);
	      oneDaylocalSearch(x_var, *it);
	   }
 	 }
     }
     auto stop = high_resolution_clock::now(); 
     auto duration = duration_cast<microseconds>(stop - start); 
     cout << duration.count()/1.0e6 << endl; 
     x_var = bestIndividual;
     calculateFeasibilityDegree(x_var, obj_values[0], obj_values[1]);
     calculateVariability(x_var, obj_values);
     cout << "sale--- " << obj_values<<endl;
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
	  ofs<< " , "<<MPP_problem->v_opt_dishes[(*a)][x_var[i*N_OPT_DAY + (*a)]].description;
    	   for(auto ij = MPP_problem->dic_nut_id.begin(); ij !=  MPP_problem->dic_nut_id.end(); ij++)
	   {
		ofs <<" , \" (";
		for(int c = 0; c < MPP_problem->conf_day.size(); c++)
		{
		double sum_nut = 0.0;
     		  for(auto t =MPP_problem->conf_day[c].begin(); t != MPP_problem->conf_day[c].end(); t++) 
	              sum_nut +=MPP_problem->v_opt_dishes[*t][x_var[i*N_OPT_DAY + (*t)]].v_nutrient_value[ij->second];
			if( c>0) ofs<<",";
			ofs <<sum_nut ;
		}
 		 ofs<<") ["<<MPP_problem->v_constraints[ij->second].min<<","<<MPP_problem->v_constraints[ij->second].max<<"] \"" ;
	   }
	ofs<<"\n";
   }
   ofs.close();
}
void MPP::First_Improvement_Hill_Climbing_swap(vector<Neighbor_swap> &neighbors, vector<int> &best_sol, vector<double> &best_objs)
{
  bool improved= true;
  vector<int> current_sol = best_sol;
  vector<double> current_objs = best_objs;
  while(improved)
  {
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for(int i = 0; i < neighbors.size(); i++)
     {
	swap_days(current_sol, neighbors[i].day1, neighbors[i].day2);
	calculateVariability(current_sol, current_objs);
	if( comp_objs(current_objs, best_objs)) 
	{
            improved = true;
	    best_objs= current_objs;
	    swap_days(best_sol, neighbors[i].day1, neighbors[i].day2);
	}
	else 
	    swap_days(current_sol, neighbors[i].day1, neighbors[i].day2);
     }
  }
}
void MPP::swap_days(vector<int> &data, int day1, int day2)
{
   for(int ii = 0; ii < MPP_problem->opt_conf.size(); ii++)
   {
      swap(data[day1*N_OPT_DAY + ii], data[day2*N_OPT_DAY + ii] );
   }
}
void  MPP::First_Improvement_Hill_Climbing(vector<Neighbor> &neighbors, vector<int> &best_sol, vector<double> &best_objs)
{
   //incremental evaluation values...
   struct Solution_LS best;
   best.obj_values = best_objs;
   best.x_var = best_sol;
   vector<double> new_objs = best_objs;
   init_incremental_evaluation(best);
   bool improved = true;
   while(improved)
   {
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for (int i = 0; i < neighbors.size(); i++)
     {
	int day = neighbors[i].variable/N_OPT_DAY, opt =  neighbors[i].variable%N_OPT_DAY;
        struct infoDishes &dish_in = MPP_problem->v_opt_dishes[opt][neighbors[i].newValue], &dish_out = MPP_problem->v_opt_dishes[opt][best.x_var[day*N_OPT_DAY + opt]];
	if(best.uniq_per_day[day][dish_in.description]>0) continue;
        //incremental evaluation...
	inc_eval(best, neighbors[i], new_objs);
        if(comp_objs(new_objs, best.obj_values))
	{
	   improved = true;
	   update_inc(best, neighbors[i], new_objs); 
	}
      }
    }
    best_sol = best.x_var;
    best_objs = best.obj_values;
    evaluate(best_sol, best_objs); //check again.. the overall values....//avoid numerical error of the incremental evaluation sum..
}
////////////////////////////////////////////////////////////////////////////////////////////////
void MPP::inc_eval(struct Solution_LS &current, Neighbor &new_neighbor, vector<double> &new_objs)
{
    int num_nutr = (int)MPP_problem->v_constraints.size();
    vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
    vector<int> &v_constraint_global = MPP_problem->v_constraint_global, &v_constraint_day = MPP_problem->v_constraint_day;
    vector<vector<infoDishes> > &v_opt_dishes = (MPP_problem->v_opt_dishes);
    int day =  new_neighbor.variable/N_OPT_DAY;
    int opt = new_neighbor.variable%N_OPT_DAY;
    double new_partial_infeasibility_day = 0.0, original_partial_infeasibility_day = 0.0 ;
    double new_partial_infeasibility_global = 0.0, original_partial_infeasibility_global = 0.0 ;
   for(int b = 0; b < MPP_problem->opt_conf[opt].size(); b++)
   {	
    int a = MPP_problem->opt_conf[opt][b];
    for(unsigned int j = 0; j < num_nutr; j++)
    {
       	//update sumatory of nutriments....
	double new_nut_value = (-v_opt_dishes[opt][current.x_var[new_neighbor.variable]].v_nutrient_value[j] + v_opt_dishes[opt][new_neighbor.newValue].v_nutrient_value[j]);
          double minv = v_constraints[j].min;
          double maxv = v_constraints[j].max;
	  double middle = (maxv+minv)*0.5;
       if(v_constraints[j].type == DIARIA)
       {
	  double nut = current.nutriment_per_day[a][day][j] + new_nut_value;
  	  double original_nut = current.nutriment_per_day[a][day][j];	
	  if( nut < minv)new_partial_infeasibility_day+= ((minv- nut)/middle)*((minv - nut)/middle)*WEIGHT_DAY;
	  else if (nut > maxv) new_partial_infeasibility_day+=((nut - maxv)/middle)*((nut - maxv)/middle)*WEIGHT_DAY;
	  if( original_nut  < minv)original_partial_infeasibility_day += ((minv - original_nut)/middle)*((minv - original_nut)/middle)*WEIGHT_DAY;
	  else if (original_nut > maxv) original_partial_infeasibility_day +=((original_nut - maxv)/middle)*((original_nut - maxv)/middle)*WEIGHT_DAY;
       }
       else if(v_constraints[j].type == GLOBAL)
       { 
          double nut = current.globalPlan[a][j] + new_nut_value;
          double original_nut = current.globalPlan[a][j];
          if(nut < minv) new_partial_infeasibility_global+= ((minv - nut)/(middle))*((minv - nut)/(middle));
          else if (nut > maxv) new_partial_infeasibility_global+=((nut - maxv)/middle)*((nut - maxv)/middle);
          if(original_nut < minv) original_partial_infeasibility_global += ((minv - original_nut)/(middle))*((minv- original_nut)/(middle));
          else if ( original_nut > maxv) original_partial_infeasibility_global +=((original_nut - maxv)/middle)*((original_nut - maxv)/middle);
       }
     }
   }
   new_objs[0]  = current.obj_values[0] - original_partial_infeasibility_day + new_partial_infeasibility_day;
   new_objs[1]  = current.obj_values[1] - original_partial_infeasibility_global + new_partial_infeasibility_global;

   if(new_objs[0] != current.obj_values[0]) return; //kind of branch procedure....
   if(new_objs[1] != current.obj_values[1]) return;
   //variability... this code-part will be optimized...
   int tmp = current.x_var[new_neighbor.variable];
   current.x_var[new_neighbor.variable] = new_neighbor.newValue;
   calculateVariability(current.x_var, new_objs); 
   current.x_var[new_neighbor.variable] = tmp;
  return;
  int time_opt = MPP_problem->inv_unique_opt_time[opt];
  struct infoDishes &dish_out = v_opt_dishes[opt][current.x_var[day*N_OPT_DAY + opt]];
  struct infoDishes &dish_in = v_opt_dishes[opt][new_neighbor.newValue];
  int max_day_in = (dish_in.favorite)?DAYS_FAVORITE:DAYS_NO_FAVORITE;
  int max_day_out = (dish_out.favorite)?DAYS_FAVORITE:DAYS_NO_FAVORITE;
  int id_day_in = dish_in.description;
  int id_day_out = dish_out.description;
  vector<int> &old_id_days_in = current.time_id_day_table[time_opt][id_day_in], &old_id_days_out = current.time_id_day_table[time_opt][id_day_out]; 
  int size_in = old_id_days_in.size(), size_out = old_id_days_out.size();

  //cout << "###################"<<endl;  
 // for(int i = 0; i < old_id_days_out.size(); i++) cout << old_id_days_out[i]<< " "; cout <<endl;
 // for(int i = 0; i < old_id_days_in.size(); i++) cout << old_id_days_in[i]<< " "; cout <<endl;
 // cout << "time: "<< time_opt << " day: " << day<<endl;
 // for(int i = 0; i < current.time_diff[time_opt].size(); i++) cout << i << " "<<current.time_diff[time_opt][i]<< " "<<endl;
  ////removing an element..
   if(size_out > 1)
   {
      int idx_v = lower_bound(old_id_days_out.begin(), old_id_days_out.end(), day) - old_id_days_out.begin();
      if( idx_v == 0) current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v+1] - old_id_days_out[idx_v])]--;
      else if(idx_v == size_out-1) current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v] - old_id_days_out[idx_v-1])]--;
      else
      {
         current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v+1] - old_id_days_out[idx_v])]--;
         current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v] - old_id_days_out[idx_v-1])]--;
        if(size_out > 2)
         current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v+1] - old_id_days_out[idx_v-1])]++;
      }
   }
   if(size_in > 0)
   {
      int idx_u = lower_bound(old_id_days_in.begin(), old_id_days_in.end(), day) - old_id_days_in.begin();
      if(size_in > 1)
	if( idx_u > 0  && idx_u < size_in)
          current.time_diff[time_opt][min(max_day_in, old_id_days_in[idx_u] - old_id_days_in[idx_u-1])]--;

      if( idx_u == 0 ) current.time_diff[time_opt][min(max_day_in, old_id_days_in[idx_u] - day)]++;
      else if( idx_u == size_in ) current.time_diff[time_opt][min(max_day_in, day - old_id_days_in[idx_u-1])]++;
      else
      {
	current.time_diff[time_opt][min(max_day_in, old_id_days_in[idx_u] - day)]++;
	current.time_diff[time_opt][min(max_day_in, day-old_id_days_in[idx_u-1])]++;
      }
   }

//  for(int i = 0; i < current.time_diff[time_opt].size(); i++) cout << i << " "<<current.time_diff[time_opt][i]<< " "<<endl;
   double t;
   int i;
   for(i = 0; i < current.time_diff[time_opt].size(); i++) 
        if(current.time_diff[time_opt][i] > 0)
	{
         //new_objs[time_opt+1] = f(make_pair(i, current.time_diff[time_opt][i]));
         t  = f(make_pair(i, current.time_diff[time_opt][i]));
	  break;
	}
	if( i >= (int)current.time_diff[time_opt].size())
         t = f(make_pair(nDias+1, 0));
//	 new_objs[time_opt+2] = t;
         //new_objs[time_opt+1] = f(make_pair(nDias+1, 0));
	
	if( fabs(t - new_objs[time_opt+2]) > 0.001)
	{
	  cout <<"nope..."<<endl;
	cout << t << " " << new_objs[time_opt+2]<<endl;
	exit(0);
	}

///    //restore values...
    if(size_out > 1)
   {
      int idx_v = lower_bound(old_id_days_out.begin(), old_id_days_out.end(), day) - old_id_days_out.begin();
      if( idx_v == 0) current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v+1] - old_id_days_out[idx_v])]++;
      else if(idx_v == size_out-1) current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v] - old_id_days_out[idx_v-1])]++;
      else
      {
         current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v+1] - old_id_days_out[idx_v])]++;
         current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v] - old_id_days_out[idx_v-1])]++;
        if(size_out > 2)
         current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v+1] - old_id_days_out[idx_v-1])]--;
      }
   }
   if(size_in > 0)
   {
      int idx_u = lower_bound(old_id_days_in.begin(), old_id_days_in.end(), day) - old_id_days_in.begin();
      if(size_in > 1)
	if( idx_u > 0  && idx_u < size_in)
          current.time_diff[time_opt][min(max_day_in, old_id_days_in[idx_u] - old_id_days_in[idx_u-1])]++;

      if( idx_u == 0 ) current.time_diff[time_opt][min(max_day_in, old_id_days_in[idx_u] - day)]--;
      else if( idx_u == size_in ) current.time_diff[time_opt][min(max_day_in, day - old_id_days_in[idx_u-1])]--;
      else
      {
	current.time_diff[time_opt][min(max_day_in, old_id_days_in[idx_u] - day)]--;
	current.time_diff[time_opt][min(max_day_in, day-old_id_days_in[idx_u-1])]--;
      }
   }

}

bool MPP::comp_objs(vector<double> &variability_v1, vector<double> &variability_v2)
{
  //feasibility
  if( variability_v1[0] < variability_v2[0]) return true;
  else if( variability_v1[0] > variability_v2[0]) return false;
  if( variability_v1[1] < variability_v2[1]) return true;
  else if( variability_v1[1] > variability_v2[1]) return false;


  //variability tchebycheff approach
   double max1 = -nDias, max2=-nDias, sum1=0.0, sum2=0.0;
//  if( variability_v1[0] == 0) 
//	MPP_problem->weights =  {0.15, 0.3, 0.05, 0.05, 0.3, 0.15};
  for(int a = 0; a < MPP_problem->priority_time.size(); a++)
  {
     int time_opt =a;// MPP_problem->priority_time[a];
     double v1 = variability_v1[time_opt+2];
     double v2 = variability_v2[time_opt+2];
	double ref = nDias*10;
	sum1 += fabs(v1-ref);
	sum2 += fabs(v2-ref);
	max1 = max(max1,fabs(v1-(ref))/MPP_problem->weights[time_opt]);
	max2 = max(max2,fabs(v2-(ref))/MPP_problem->weights[time_opt]);
  }
  max1 += 0.0001*sum1; 
  max2 += 0.0001*sum2; 
 // cout << max1 << " "<<max2<<endl;
 if( max1 < max2) return true;
 return false;
}

void MPP::update_inc(struct Solution_LS &current, Neighbor &neighbor, vector<double> &new_objs)
{
    int num_nutr = (int)MPP_problem->v_constraints.size();
    vector<vector<infoDishes> > &v_opt_dishes = (MPP_problem->v_opt_dishes);
    int day =  neighbor.variable/N_OPT_DAY;
    int opt = neighbor.variable%N_OPT_DAY;
    struct infoDishes &dish_in = MPP_problem->v_opt_dishes[opt][neighbor.newValue], &dish_out = MPP_problem->v_opt_dishes[opt][current.x_var[day*N_OPT_DAY + opt]];
   for(int b = 0; b < MPP_problem->opt_conf[opt].size(); b++)
   {	
    int a = MPP_problem->opt_conf[opt][b];
    for(unsigned int j = 0; j < num_nutr; j++)//check new neighbor..
    {
       	//update sumatory of nutriments....
	double diff = (-v_opt_dishes[opt][current.x_var[neighbor.variable]].v_nutrient_value[j] + v_opt_dishes[opt][neighbor.newValue].v_nutrient_value[j]);
        current.nutriment_per_day[a][day][j] += diff;
	current.globalPlan[a][j] += diff;
    }
   }
  int time_opt = MPP_problem->inv_unique_opt_time[opt];
  int max_day_in = (dish_in.favorite)?DAYS_FAVORITE:DAYS_NO_FAVORITE;
  int max_day_out = (dish_out.favorite)?DAYS_FAVORITE:DAYS_NO_FAVORITE;
  int id_day_in = dish_in.description;
  int id_day_out = dish_out.description;
  vector<int> &old_id_days_in = current.time_id_day_table[time_opt][id_day_in], &old_id_days_out = current.time_id_day_table[time_opt][id_day_out]; 
  int size_in = old_id_days_in.size(), size_out = old_id_days_out.size();
  ////removing an element..

   int idx_v = lower_bound(old_id_days_out.begin(), old_id_days_out.end(), day) - old_id_days_out.begin();
   if(size_out > 1)
   {
      if( idx_v == 0) current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v+1] - old_id_days_out[idx_v])]--;
      else if(idx_v == size_out-1) current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v] - old_id_days_out[idx_v-1])]--;
      else
      {
         current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v+1] - old_id_days_out[idx_v])]--;
         current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v] - old_id_days_out[idx_v-1])]--;
        if(size_out > 2)
         current.time_diff[time_opt][min(max_day_out, old_id_days_out[idx_v+1] - old_id_days_out[idx_v-1])]++;
      }
   }
   int idx_u = lower_bound(old_id_days_in.begin(), old_id_days_in.end(), day) - old_id_days_in.begin();
   if(size_in > 0)
   {
      if(size_in > 1)
        if( idx_u > 0  && idx_u < size_in)
        {
          current.time_diff[time_opt][min(max_day_in, old_id_days_in[idx_u] - old_id_days_in[idx_u-1])]--;
        }

      if( idx_u == 0 ) current.time_diff[time_opt][min(max_day_in, old_id_days_in[idx_u] - day)]++;
      else if( idx_u == size_in ) current.time_diff[time_opt][min(max_day_in, day - old_id_days_in[idx_u-1])]++;
      else
      {
        current.time_diff[time_opt][min(max_day_in, old_id_days_in[idx_u] - day)]++;
        current.time_diff[time_opt][min(max_day_in, day-old_id_days_in[idx_u-1])]++;
      }
   }
//cout << "##################"<<endl;  
//for(int i = 0; i < old_id_days_out.size(); i++) cout << old_id_days_out[i]<< " "; cout <<endl;
//  for(int i = 0; i < old_id_days_in.size(); i++) cout << old_id_days_in[i]<< " "; cout <<endl;
   if(size_in == 0) old_id_days_in.push_back(day);
   else
   {
      old_id_days_in.insert(old_id_days_in.begin()+idx_u, day);
   }
   old_id_days_out.erase(old_id_days_out.begin()+idx_v);
    

//for(int i = 0; i < old_id_days_out.size(); i++) cout << old_id_days_out[i]<< " "; cout <<endl;
//  for(int i = 0; i < old_id_days_in.size(); i++) cout << old_id_days_in[i]<< " "; cout <<endl;

//cout << "##################"<<endl;  

   current.x_var[neighbor.variable]= neighbor.newValue;
   current.obj_values = new_objs;
   current.uniq_per_day[day][dish_in.description]++, current.uniq_per_day[day][dish_out.description]--;	
}


void MPP::init_incremental_evaluation(struct Solution_LS &current)
{ 
        //feasibility information
	int num_nutr = (int)MPP_problem->v_constraints.size();
	vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
        vector<int> &v_constraint_global = MPP_problem->v_constraint_global, &v_constraint_day = MPP_problem->v_constraint_day;
        int &max_description_id = MPP_problem->max_description_id;
        vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;
   	current.uniq_per_day.assign(nDias+1, vector<int> (MPP_problem->max_description_id+1, 0));
        current.globalPlan.assign((int)MPP_problem->conf_day.size(),vector<double> (num_nutr, 0.0 ));
	current.nutriment_per_day.assign( (int)MPP_problem->conf_day.size(), vector<vector<double> > (nDias, vector<double> (num_nutr, 0)));
	current.obj_values.assign(N_TIMES+2, 0.0);
        for(int a = 0; a < MPP_problem->conf_day.size(); a++)
	{
 	  for(int j = 0; j < num_nutr; j++)
	  {
	   for(int i = 0; i < nDias; i++)
	   {
		 for(int b = 0; b < MPP_problem->conf_day[a].size(); b++)
		 {
		   int k = MPP_problem->conf_day[a][b];
	    	   current.nutriment_per_day[a][i][j] += v_opt_dishes[k][current.x_var[i*N_OPT_DAY + k]].v_nutrient_value[j];
		 }
	        current.globalPlan[a][j] += current.nutriment_per_day[a][i][j]; 	
	        if( v_constraints[j].type == DIARIA)
                {
                   double minv = v_constraints[j].min;
                   double maxv = v_constraints[j].max;
	  	   double middle = (maxv+minv)*0.5;
	           double nut = current.nutriment_per_day[a][i][j];
	           if( nut < minv) current.obj_values[0] += ((minv - nut)/middle)*((minv - nut)/middle)*WEIGHT_DAY;
	           else if (nut > maxv) current.obj_values[0] +=((nut - maxv)/middle)*((nut - maxv)/middle)*WEIGHT_DAY;
                }
	   }
	   if( v_constraints[j].type == GLOBAL)
           {
                   double minv = v_constraints[j].min;
                   double maxv = v_constraints[j].max;
	  	   double middle = (maxv+minv)*0.5;
	           double nut = current.globalPlan[a][j];
	           if( nut < minv) current.obj_values[1] += ((minv - nut)/middle)*((minv - nut)/middle);
	           else if (nut > maxv) current.obj_values[1] += ((nut - maxv)/middle)*((nut - maxv)/middle);
            }
	  }
	}
    //variability information...
    vector<pair<int, int>> min_dcn(N_TIMES, make_pair(nDias+1, 0)), min_dcn_cat(N_TIMES, make_pair(nDias+1, 0));
    current.time_id_day_table.assign(N_TIMES, vector< vector<int> > (max_description_id+1));
    current.time_diff.assign(N_TIMES, vector< int > (nDias+1, 0)); 
    current.time_diff_cat.assign(N_TIMES, vector< int > (nDias+1, 0)); 
    vector<int> last_day_seen(N_TIMES*(max_description_id+1), -1), last_day_seen_cat(N_TIMES*(max_description_id+1), -1);
    double v_global_id = 0, v_global_cat = 0;
    for(int d = 0; d < nDias; d++)
    {
      for(int i = 0; i < N_TIMES; i++) 
      {
        for(int j = 0; j < MPP_problem->unique_opt_time[i].size(); j++)
        {
	   int opt = MPP_problem->unique_opt_time[i][j];
	   struct infoDishes &dish = v_opt_dishes[opt][current.x_var[d*N_OPT_DAY + opt]];
           current.time_id_day_table[i][dish.description].push_back(d);
	   current.uniq_per_day[d][dish.description]++;
           if(last_day_seen[i*N_TIMES+dish.description] != -1)
           {
              int diff = d-last_day_seen[i*N_TIMES+dish.description];
              diff = min((dish.favorite)?DAYS_FAVORITE:DAYS_NO_FAVORITE, diff);
              current.time_diff[i][diff]++;
              update_dcn_pair(diff, min_dcn[i]);
           }
           last_day_seen[i*N_TIMES+dish.description] = d;

	   if(last_day_seen_cat[i*N_TIMES+dish.category] != -1)
           {
              int diff = d-last_day_seen_cat[i*N_TIMES+dish.category];
              diff = min((dish.favorite)?DAYS_FAVORITE:DAYS_NO_FAVORITE, diff);
              current.time_diff_cat[i][diff]++;
              update_dcn_pair(diff, min_dcn_cat[i]);
           }
           last_day_seen_cat[i*N_TIMES+dish.category] = d;
	}
      }
  //  v_global_cat += (min_dcn_cat.first + ( 1.0 - ((double)min_dcn_cat.second/(double)nDias)));
    }
    for(int i = 0; i < N_TIMES; i++)
    {
      current.obj_values[i+2] = f(min_dcn[i]);
    }

//    for(int i = 0; i < current.time_id_day_table.size(); i++)
//    {
//	cout << "time.. " <<i <<endl;
//	for(int j = 0; j < current.time_id_day_table[i].size(); j++)
//	{
//	   if(current.time_id_day_table[i][j].empty()) continue;
//	   for(int k = 0; k < current.time_id_day_table[i][j].size(); k++)
//		cout << current.time_id_day_table[i][j][k]<< " ";
//	   cout <<endl;
//	}
//    }
//   exit(0);

//   double total_variability = W_VAR_GLOBAL*v_global_id;// + W_VAR_GLOBAL_CAT*v_global_cat;
}
void MPP::calculateFeasibilityDegree(vector<int> &sol, double &feas_day, double &feas_global){
	feas_day = 0.0;
	feas_global = 0.0;
	int num_nutr = (int)MPP_problem->v_constraints.size();
	vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
	badDaysFeas.clear();
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
	   	   dayNutr += MPP_problem->v_opt_dishes[k][sol[j*N_OPT_DAY+k]].v_nutrient_value[i];
		 }
	   	globalNutr += dayNutr;
	   	if(v_constraints[i].type == DIARIA)
             	{
	               if(dayNutr  < minv) feas_day += ((minv - dayNutr)/middle)*((minv - dayNutr)/middle)*WEIGHT_DAY, badDaysFeas.insert(j);
	               else if (dayNutr > maxv) feas_day +=((dayNutr - maxv)/middle)*((dayNutr - maxv)/middle)*WEIGHT_DAY, badDaysFeas.insert(j);
           	}
	      }
	      if(v_constraints[i].type == GLOBAL)
              {
	        if(globalNutr < minv)
	        {
	          double v = ((minv - globalNutr)/middle)*((minv - globalNutr)/middle);
	          feas_global += v;
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
	          feas_global +=v;
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
void MPP::oneDaylocalSearch(vector<int> &solution, int day) {
     vector<Neighbor> neighbors;
     int num_nutr = (int)MPP_problem->v_constraints.size();
     vector<int> best_solution(solution.begin()+day*N_OPT_DAY, solution.begin() + (day+1)*N_OPT_DAY );

     for (int j = 0; j <N_OPT_DAY; j++)
     {
       for (int k = 0; k < (int) MPP_problem->v_opt_dishes[j].size(); k++)
       {
	     Neighbor n;
	     n.variable = j;
	     n.newValue = k;
	     neighbors.push_back(n);
       }
     }
     double best_feasibility = 0.0; 
     feasibility_day(best_solution, best_feasibility);
//     cout << "entra day " << best_feasibility<<endl;
     double current_feasibility = best_feasibility;
     vector<int> current_solution = best_solution;
     auto start = high_resolution_clock::now(); 
     for (int i = 0; i < 2000; i++)
     {
	//first improvement hc
 	First_Improvement_Hill_Climbing_Day(neighbors, current_solution, current_feasibility);
	///save the best....
 	if( current_feasibility < best_feasibility)
	{
	//   cout << current_feasibility<<endl;
	   best_feasibility = current_feasibility;
	   best_solution = current_solution;
	}
	else
	{
	   current_solution = best_solution;
	   current_feasibility = best_feasibility;
	}
	//restart....
	int which = rand() % N_OPT_DAY;
	perturb_opt(current_solution, 0, which);
	feasibility_day(current_solution, current_feasibility);
	if(current_feasibility == 0.0) break;
     }
     auto stop = high_resolution_clock::now(); 
     auto duration = duration_cast<microseconds>(stop - start); 
//     cout << duration.count()/1.0e6 << endl; 
//     cout << "sale--- " << best_feasibility<<endl;
     //exit(0);
    for(int i = 0; i  < N_OPT_DAY; i++) solution[day*N_OPT_DAY + i] = best_solution[i];
}
void  MPP::First_Improvement_Hill_Climbing_Day(vector<Neighbor> &neighbors, vector<int> &best_solution, double &best_feasibility)
{
   //incremental evaluation values...
   double current_feasibility = best_feasibility; 
   vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;
   int &max_description_id = MPP_problem->max_description_id;
   vector<bool> uniq_per_day(max_description_id, false);

   int num_nutr = (int)MPP_problem->v_constraints.size();
   vector<vector<double>> nut_info;
   for(int i = 0 ; i  < N_OPT_DAY; i++) uniq_per_day[v_opt_dishes[i][best_solution[i]].description] = true;
   init_inc_eval_day(best_solution, nut_info);
   bool improved = true;
   while(improved)
   {
     improved = false;
     random_shuffle(neighbors.begin(), neighbors.end());
     for (int i = 0; i < neighbors.size(); i++)
     {
	int day = neighbors[i].variable/N_OPT_DAY, opt =  neighbors[i].variable;
        struct infoDishes &dish_in = MPP_problem->v_opt_dishes[opt][neighbors[i].newValue];
        struct infoDishes &dish_out = MPP_problem->v_opt_dishes[opt][best_solution[opt]];
	if(uniq_per_day[dish_in.description]) continue;
        //incremental evaluation...
	double new_feasibility = inc_eval_day(neighbors[i], nut_info , best_solution, best_feasibility);
	if( new_feasibility < best_feasibility)
	{
	   improved = true;
	   update_inc_day(nut_info, neighbors[i], best_solution);
	   best_feasibility = new_feasibility;
   	 //  init_inc_eval_day(best_solution, nut_info);
	   best_feasibility = max(0.0, best_feasibility);
	  // feasibility_day(best_solution, new_feasibility);
	  // cout << fabs(new_feasibility-best_feasibility)<<endl;
	   uniq_per_day[dish_in.description] = true;
	   uniq_per_day[dish_out.description] = false;
	}
      }
    }
}
void MPP::feasibility_day(vector<int> &best_solution, double &best_feasibility)
{
   vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
   vector<int> &v_constraint_day = MPP_problem->v_constraint_day;
   //get daily contribution
   best_feasibility = 0.0;
   for(int a = 0; a < MPP_problem->conf_day.size(); a++)
   {
   for(int ij = 0; ij < v_constraint_day.size(); ij++)
     {
        double dayNutr = 0;
        int i =  v_constraint_day[ij];
        double minv = v_constraints[i].min;
        double maxv = v_constraints[i].max;
	double middle = (maxv+minv)*0.5;
	for(int b = 0; b < MPP_problem->conf_day[a].size(); b++)
	{
	  int k = MPP_problem->conf_day[a][b]; 
	  dayNutr += MPP_problem->v_opt_dishes[k][best_solution[k]].v_nutrient_value[i];
        }
	if(dayNutr  < minv) best_feasibility += ((minv - dayNutr)/middle)*((minv - dayNutr)/middle);
	else if (dayNutr > maxv) best_feasibility +=((dayNutr - maxv)/middle)*((dayNutr - maxv)/middle);
     }
   }
}
void MPP::init_inc_eval_day(vector<int> &current_solution, vector<vector<double>> &nut_info)
{
   int num_nutr = (int)MPP_problem->v_constraints.size();
   vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
   vector<int> &v_constraint_day = MPP_problem->v_constraint_day;
   vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;

   nut_info.assign((int)MPP_problem->conf_day.size(), vector<double>(num_nutr, 0.0));
   for(int a = 0; a < MPP_problem->conf_day.size(); a++)
   {
     for(int i = 0; i < v_constraint_day.size(); i++)
     {
       int j =  v_constraint_day[i];
       for(int b = 0; b < MPP_problem->conf_day[a].size(); b++)
       {
         int k = MPP_problem->conf_day[a][b];
         nut_info[a][j] += v_opt_dishes[k][current_solution[k]].v_nutrient_value[j];
       }
     }
   }

}
double MPP::inc_eval_day(Neighbor &new_neighbor, vector<vector<double>> &nut_info, vector<int> &current_solution, double current_feasibility)
{
   vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;
   vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
   vector<int> &v_constraint_day = MPP_problem->v_constraint_day;
   double new_partial_infeasibility_day = 0.0, original_partial_infeasibility_day = 0.0 ;
   int opt = new_neighbor.variable;
   for(int b = 0; b < MPP_problem->opt_conf[opt].size(); b++)
   {	
       int a = MPP_problem->opt_conf[opt][b];
       for(int i = 0; i < v_constraint_day.size(); i++)
       {
          int j =  v_constraint_day[i];
         //update sumatory of nutriments....
         double new_nut_value = (-v_opt_dishes[opt][current_solution[opt]].v_nutrient_value[j] + v_opt_dishes[opt][new_neighbor.newValue].v_nutrient_value[j]);
         double minv = v_constraints[j].min;
         double maxv = v_constraints[j].max;
         double middle = (maxv+minv)*0.5;
         double nut = nut_info[a][j] + new_nut_value;
         double original_nut = nut_info[a][j];	
      
         if( nut < minv)new_partial_infeasibility_day+= ((minv- nut)/middle)*((minv - nut)/middle);
         else if (nut > maxv) new_partial_infeasibility_day+=((nut - maxv)/middle)*((nut - maxv)/middle);
         if( original_nut  < minv)original_partial_infeasibility_day += ((minv - original_nut)/middle)*((minv - original_nut)/middle);
         else if (original_nut > maxv) original_partial_infeasibility_day +=((original_nut - maxv)/middle)*((original_nut - maxv)/middle);
       }
   }
  return current_feasibility -original_partial_infeasibility_day + new_partial_infeasibility_day; 
}
void MPP::update_inc_day(vector<vector<double>> &nut_info, Neighbor &new_neighbor, vector<int> &current_solution)
{
   vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;
   vector<int> &v_constraint_day = MPP_problem->v_constraint_day;
   int opt = new_neighbor.variable;
   for(int b = 0; b < MPP_problem->opt_conf[opt].size(); b++)
   {	
       int a = MPP_problem->opt_conf[opt][b];
       for(int i = 0; i < v_constraint_day.size(); i++)
       {
          int j =  v_constraint_day[i];
         nut_info[a][j] = nut_info[a][j]-v_opt_dishes[opt][current_solution[new_neighbor.variable]].v_nutrient_value[j] + v_opt_dishes[opt][new_neighbor.newValue].v_nutrient_value[j];
       }
   }
   current_solution[opt] = new_neighbor.newValue;
}
void MPP::calculateVariability(vector<int> &sol, vector<double> &objs_var)
{
   vector<vector<infoDishes> > &v_opt_dishes = MPP_problem->v_opt_dishes;
   double variability_global = 0.0, variability_cat_day=0.0;
   int &max_description_id = MPP_problem->max_description_id;
   double v_global = 0, v_global_id = 0, v_global_cat = 0;
   badDaysVar.clear();
   vector<int> last_day_seen(N_TIMES*(max_description_id+1), -1), last_day_seen_cat(N_TIMES*(max_description_id+1), -1);
   vector<pair<int, int>> min_dcn(N_TIMES, make_pair(nDias+1, 0)), min_dcn_cat(N_TIMES, make_pair(nDias+1, 0));
   vector< vector<int> > min_time_id(N_TIMES, vector<int> (max_description_id+1, nDias+1));
   int min_dcn_day = nDias+1;
   for(int d = 0; d < nDias; d++)
   {
	for(int i = 0; i < MPP_problem->unique_opt_time.size(); i++)
	{
	   for(int j = 0; j < MPP_problem->unique_opt_time[i].size(); j++)
           {
	  	int opt = MPP_problem->unique_opt_time[i][j];
		struct infoDishes &dish = v_opt_dishes[opt][sol[d*N_OPT_DAY + opt]];
           	if(last_day_seen[i*N_TIMES+dish.description] != -1)
           	{
           	   int diff = d-last_day_seen[i*N_TIMES+dish.description];
	   	   diff = min((dish.favorite)?DAYS_FAVORITE:DAYS_NO_FAVORITE, diff);
	   	   update_dcn_pair(diff, min_dcn[i]);
		   if( min_dcn_day == diff)
		   {
		     badDaysVar.insert(d);
		     badDaysVar.insert(last_day_seen[i*N_TIMES + dish.description]);
		   }
		   else if(min_dcn_day >diff)
		   {
		     min_dcn_day = diff;
		     badDaysVar.clear();
		     badDaysVar.insert(d);
		     badDaysVar.insert(last_day_seen[i*N_TIMES + dish.description]);
		   }
           	}
           	last_day_seen[i*N_TIMES + dish.description] = d;
	   //	if(last_day_seen_cat[i*N_TIMES + dish.category] != -1)
           //	{
           //	   int diff = d-last_day_seen_cat[i*N_TIMES + dish.category];
	   //	   update_dcn_pair(diff, min_dcn_cat[i]);
           //	}
           //	last_day_seen_cat[i*N_TIMES + dish.category] = d;
	   }
	}
   }
  
  for(int i = 0; i < N_TIMES; i++)
  {
    objs_var[i+2] = f(min_dcn[i]);
	//cout << " i: "<<i+1<< " "<<min_dcn[i].first << " " <<min_dcn[i].second<<endl;
//    v_global_id += MPP_problem->weights[i]*f(min_dcn[i]);
//    v_global_cat += f(min_dcn_cat[i]);
  }
//  return W_VAR_GLOBAL*v_global_id + W_VAR_GLOBAL_CAT*v_global_cat;
}
