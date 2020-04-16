void MPP::full_search()
{
 vector<vector<infoDishes> > times = MPP_problem->v_opt_dishes;
 vector<constraint_nutrient> &v_constraints = (MPP_problem->v_constraints);
 vector< vector<int> > feasible_solutions;
 vector<pair<double, double > > fit_sol;
 long int cont = 0, max_perm =1;
 //information to get each permutation of the options-feasible space..
 fill(x_var.begin(), x_var.end(),0);
 vector<int> v_max_opt;
 for(int max_opt = 0; max_opt < MPP_problem->opt_conf.size(); max_opt++)
 {
   if( MPP_problem->opt_conf[max_opt].empty()) continue;
   int opt_s = (int)MPP_problem->v_opt_dishes[max_opt].size();
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
      for(unsigned int k = 0; k < MPP_problem->opt_conf.size(); k++)
      {
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
  }
   ofs.close();
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
