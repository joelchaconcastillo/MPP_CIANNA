void MPP::calculateFeasibilityDegree(vector<int> &sol, double &feas){
        feas = 0.0;
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
	   	   dayNutr += MPP_problem->v_opt_dishes[k][sol[j*N_OPT_DAY+k]].v_nutrient_value[i];
		 }
	   	globalNutr += dayNutr;
	   	if(v_constraints[i].type == DIARIA)
             	{
	               if(dayNutr  < minv) feas += ((minv - dayNutr)/middle)*((minv - dayNutr)/middle)*WEIGHT_DAY, badDays.insert(j);
	               else if (dayNutr > maxv) feas +=((dayNutr - maxv)/middle)*((dayNutr - maxv)/middle)*WEIGHT_DAY, badDays.insert(j);
           	}
	      }
	      if(v_constraints[i].type == GLOBAL)
              {
	        if(globalNutr < minv)
	        {
	          double v = ((minv - globalNutr)/middle)*((minv - globalNutr)/middle);
	          feas += v;
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
	          feas +=v;
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
