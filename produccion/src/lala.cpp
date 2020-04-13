double MPP::inc_eval_var_time(vector<int> &sol, Neighbor &new_neighbor, double current_var, vector< vector< vector < int > > > &time_id_day_table, vector< vector<int> > &time_diff)
{
  vector<vector<infoDishes> > &v_opt_dishes = (MPP_problem->v_opt_dishes);
  int day =  new_neighbor.variable/N_OPT_DAY;
  int opt = new_neighbor.variable%N_OPT_DAY;
  int time = MPP_problem->inv_unique_opt_time[opt];
  struct infoDishes &dish_out = v_opt_dishes[opt][sol[day*N_OPT_DAY + opt]];
  struct infoDishes &dish_in = v_opt_dishes[opt][new_neighbor.newValue];
  bool fav_in = dish_in.favorite;
  bool fav_out = dish_out.favorite;
  int id_day_in = dish_in.description;
  int id_day_out = dish_out.description;
  if( id_day_out == id_day_in) return current_var;
  vector<int> &old_id_days_in = time_id_day_table[time][id_day_in], &old_id_days_out = time_id_day_table[time][id_day_out]; 
  vector<int> new_id_days_in, new_id_days_out; 
   pair<int, int> min_after(nDias+1, 0), min_before(nDias+1, 0);
  for(int i = 0; i < old_id_days_out.size(); i++)
  {
	if( i > 0)
	{
	   int diff = min((fav_out)?DAYS_FAVORITE:DAYS_NO_FAVORITE, old_id_days_out[i]-old_id_days_out[i-1]);
	   update_dcn_pair(diff, min_before);
	}
	if(old_id_days_out[i] != day)
           new_id_days_out.push_back(old_id_days_out[i]);   
  }
  for(int i = 0; i < old_id_days_in.size(); i++)
  {
	if( i > 0)
	{
	   int diff = min((fav_in)?DAYS_FAVORITE:DAYS_NO_FAVORITE, old_id_days_in[i]-old_id_days_in[i-1]);
	   update_dcn_pair(diff, min_before);
	}
	new_id_days_in.push_back(old_id_days_in[i]);
	if( day > old_id_days_in[i])
	new_id_days_in.push_back(day), day =-1 ;
  }
  for(int i = 1; i < new_id_days_in.size(); i++)
  {
    int diff = min((fav_in)?DAYS_FAVORITE:DAYS_NO_FAVORITE, new_id_days_in[i]-new_id_days_in[i-1]);
    update_dcn_pair(diff, min_after);
  }
  for(int i = 1; i < new_id_days_out.size(); i++)
  {
    int diff = min((fav_out)?DAYS_FAVORITE:DAYS_NO_FAVORITE, new_id_days_out[i]-new_id_days_out[i-1]);
    update_dcn_pair(diff, min_after);
  }
  double newvar = current_var;
  if( old_id_days_in.size() >1 || old_id_days_out.size()>1)
    newvar -=W_VAR_GLOBAL*f(min_before);
  if( new_id_days_in.size() >1 || new_id_days_out.size()>1)
    newvar +=W_VAR_GLOBAL*f(min_after);
  //cout <<"--"<< newvar<<endl;
 return newvar;
}
void MPP::update_inc_var(vector<int> &sol, double &current_var, Neighbor &new_neighbor, vector< vector< vector < int > > > &time_id_day_table, vector< vector<int> > &time_diff)
{
  vector<vector<infoDishes> > &v_opt_dishes = (MPP_problem->v_opt_dishes);
  int day =  new_neighbor.variable/N_OPT_DAY;
  int opt = new_neighbor.variable%N_OPT_DAY;
  int time = MPP_problem->inv_unique_opt_time[opt];
  struct infoDishes &dish_out = v_opt_dishes[opt][sol[day*N_OPT_DAY + opt]];
  struct infoDishes &dish_in = v_opt_dishes[opt][new_neighbor.newValue];
  bool fav_in = dish_in.favorite;
  bool fav_out = dish_out.favorite;
  int id_day_in = dish_in.description;
  int id_day_out = dish_out.description;
  if( id_day_out == id_day_in) return ;
  vector<int> &old_id_days_in = time_id_day_table[time][id_day_in], &old_id_days_out = time_id_day_table[time][id_day_out]; 
  vector<int> new_id_days_in, new_id_days_out; 
   pair<int, int> min_after(nDias+1, 0), min_before(nDias+1, 0);
  for(int i = 0; i < old_id_days_out.size(); i++)
  {
	if( i > 0)
	{
	   int diff = min((fav_out)?DAYS_FAVORITE:DAYS_NO_FAVORITE, old_id_days_out[i]-old_id_days_out[i-1]);
	   update_dcn_pair(diff, min_before);
	}
	if(old_id_days_out[i] != day)
           new_id_days_out.push_back(old_id_days_out[i]);   
  }
  for(int i = 0; i < old_id_days_in.size(); i++)
  {
	if( i > 0)
	{
	   int diff = min((fav_in)?DAYS_FAVORITE:DAYS_NO_FAVORITE, old_id_days_in[i]-old_id_days_in[i-1]);
	   update_dcn_pair(diff, min_before);
	}
	new_id_days_in.push_back(old_id_days_in[i]);
	if( day > old_id_days_in[i])
	new_id_days_in.push_back(day), day =-1 ;
  }
  for(int i = 1; i < new_id_days_in.size(); i++)
  {
    int diff = min((fav_in)?DAYS_FAVORITE:DAYS_NO_FAVORITE, new_id_days_in[i]-new_id_days_in[i-1]);
    update_dcn_pair(diff, min_after);
  }
  for(int i = 1; i < new_id_days_out.size(); i++)
  {
    int diff = min((fav_out)?DAYS_FAVORITE:DAYS_NO_FAVORITE, new_id_days_out[i]-new_id_days_out[i-1]);
    update_dcn_pair(diff, min_after);
  }
  if( old_id_days_in.size() >1 || old_id_days_out.size()>1)
    current_var -=W_VAR_GLOBAL*f(min_before);
  if( new_id_days_in.size() >1 || new_id_days_out.size()>1)
    current_var +=W_VAR_GLOBAL*f(min_after);
  old_id_days_in = new_id_days_in;
  old_id_days_out = new_id_days_out;
}

