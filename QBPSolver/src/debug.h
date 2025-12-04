  void DEBUG_CHECK_PROPQ(int pick) {
    /*
    cerr << "\n -- propQ-warning:" << (propQ[0].v >> 1) << " " << pick << " " << level_finished[decisionLevel()-1] << " " << level_finished[decisionLevel()-2] << " " << decisionLevel() << endl;
    if (propQ.size() > 0 && assigns[propQ[0].v >> 1] == extbool_Undef) {
      //assert(propQ.size() == 1);
      for (int tt = 0; tt < propQ.size();tt++) {
	cerr << "++ inspect " << (propQ[tt].v >> 1) << "reason in level " << vardata[propQ[tt].v >> 1].level << endl;
	if (tt==0) assert(assigns[propQ[tt].v >> 1] == extbool_Undef);
	else assert(assigns[propQ[tt].v >> 1] != extbool_Undef);
      }
      //assert dass alle Member of propQ gesetzt sind, ausser propQ[0]
      cerr << "\n ++ propQ-warning:" << (propQ[0].v >> 1) << " " << pick << " " << level_finished[decisionLevel()-1] << " " << level_finished[decisionLevel()-2] << " " << decisionLevel() << endl;
    }
    */
  }

void  DEBUG_OUT_5SCENARIOS() {
    if (info_level == 6) {  
      cerr << endl;         
      for (int i = 0; i < 5/*scenario.size()*/; i++) 
	if (i < scenario.size()) cerr << (int)assigns[scenario[i]];
	else cerr << " ";
      for (int i = 1; i < trail_lim.size();i++)
	if (eas[trail[trail_lim[i]-1]]==EXIST) cerr << (int)assigns[trail[trail_lim[i]-1]];
	else cerr << (int)assigns[trail[trail_lim[i]-1]]+2;
      cerr << endl;
    }
}

void DEBUG_VARIABLE_MISSING() {
	      cerr << "Error: Variables missing: ";
  	      std::vector<int> ar;
	      for (int u=0;u<nVars();u++)
	        ar.push_back(0);
	      for (int u=0;u<trail.size();u++)
	        ar[trail[u]]=ar[trail[u]]+1;
	      for (int u=0;u<nVars();u++)
	        if(ar[u]==0 || ar[u] > 1) cerr << u << "," << (int)assigns[u] << "," << isFixed(u) << "," << ar[u] << " ";
	      cerr << " --> " << trail.size() << "," << nVars() << endl;
}
