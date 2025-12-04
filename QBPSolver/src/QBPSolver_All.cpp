/*
*
* Yasol: QBPSolver_All.cpp -- Copyright (c) 2012-2018 Ulf Lorenz, Michael Hartisch
*
* Permission is hereby granted, free of charge, to any person obtaining a
* copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
* LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
* OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
* WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
#include <iostream>
    using namespace std;

    #include "QBPSolver.h"
#include "yInterface.h"
    #include <cmath>
#include <queue>
#include <vector>
#include <set>

	#include "FeasibilityPump.h"
#define NO_DEBUG 1
	#define PRINT_PROP 0

void QBPSolver::WriteSolutionFile(coef_t value, double time, string status, string filename, bool removeDummy){ 
  if(!getWriteOutputFile()) return;
  double gap=0;
  if(value==n_infinity || value==-n_infinity) return;
    //status="INFEASIBLE";// (universal player can force violation of existential constraint system)";
   if(objInverted) gap = abs(100.0*(+global_dual_bound - incumbentBest) / (abs(incumbentBest)+1e-10) );
   else gap = abs(100.0*(-global_dual_bound + incumbentBest) / (abs(incumbentBest)+1e-10) );  
   //if(objInverted) gap = abs(100.0*(+global_dual_bound - value) / (abs(value)+1e-10) );
   //else gap = abs(100.0*(-global_dual_bound + value) / (abs(value)+1e-10) );

  if(objInverted) value*=-1;

  CommPrint Sol;
  string Final="";
  if (filename == "") Final = getInputFileName()+".sol";
  else Final = filename;
  const char *Start=Final.c_str();
  ifstream myfile (Start);
  int Counter=0;
  string TmpA=Final;
  while (myfile.is_open()){
    myfile.close();       
    Counter++;
    
    TmpA=Final+std::to_string(Counter);
    const char *tmp=TmpA.c_str();
    myfile.open(tmp);
    if(myfile.is_open()) continue;
  }
  Final=TmpA;
  if(getShowAnything()) cerr << "solution file written to " << Final << endl;
  //toWrite is a string that is appended to the solution file
  //this string is written in (pseudo)XML
  //the write stage is divided in stages for header+quality, linearConstraints and variables
  //values have to be added between \". example: \"value\"
  string toWrite = "";
  toWrite += "<?xml version = \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
  toWrite += "<YasolSolution version=\"1\">\n";
  toWrite += " <header\n";
  toWrite += "   ProblemName=\""+ getInputFileName()+ "\"\n";
  toWrite += "   SolutionName=\""+Final+"\"\n";
  if (status.compare("FEASIBLE")!=0)
    toWrite += "   ObjectiveValue=\""+std::to_string((double)value)+"\"\n";
  //toWrite += " ObjectiveOffset=\""+std::to_string(finalOffset)+"\"\n";
  //toWrite += " solutionStatusString=\""+Status+"\"\n";
  toWrite += "   Runtime=\""+(time==-1?"TIMEOUT":(std::to_string(time).substr(0,std::to_string(time).find(".") + 4)+" seconds"))+"\"\n";
  toWrite += "   DecisionNodes=\""+std::to_string((int)getNumberOfDecisions())+"\"\n";
  toWrite += "   PropagationSteps=\""+std::to_string((int)getNumberOfPropagationSteps())+"\"\n";
  toWrite += "   LearntConstraints=\""+std::to_string((int)getNumberOfLearntConstraints())+"\"/>\n";
  toWrite += " <quality>\n";
  toWrite += "   SolutionStatus=\""+status+"\"\n";
  if (status.compare("FEASIBLE")!=0){
    toWrite += "   Gap=\""+std::to_string(gap)+"\"\n";
  	toWrite += "   Dual=\""+std::to_string(objInverted?global_dual_bound:-global_dual_bound)+"\"\n";
  }
  toWrite += " </quality>\n";
  Sol.solprint(toWrite,Final.c_str());

  //reset toWrite because method solprint is in append mode
  toWrite = " <linearConstraints>\n";
  toWrite += " </linearConstraints>\n";
  Sol.solprint(toWrite,Final.c_str());
  toWrite = " <variables>\n";
  for(int i=0; i<nVars()-removeDummy; i++){
    if (((yInterface*)yIF)->integers[i].bitcnt > 1) {
      //get real name without underscore from binarized variables
      //create char array from name string
      const char *name = ((yInterface*)yIF)->integers[i].name.c_str();
      //positon of underscore in char array
      int bpoint = 0;
      //get position of underscore
      for(int j=0; j<=((yInterface*)yIF)->integers[i].name.length(); j++){
        if(name[j]=='_'){
          bpoint = j;
        }
      }
      
      //create new name for binarized variable with substring
      std::string noBinName = ((yInterface*)yIF)->integers[i].name.substr(0,bpoint);

      //change to real value
      if (((yInterface*)yIF)->integers[i].pt2leader == i) {
        int res = 0;
	for (int z = 0;z < ((yInterface*)yIF)->integers[i].bitcnt;z++) {
	  if(block[i] == 1){
	    res = 2*res + getFirstStageSolutionValue(i+z);
	  }else {
	    res = 2*res + PV[0][i+z];
	  }
	}
	//index of binarized variable
	std::string vIndex = std::to_string(i) + "-" + std::to_string(((yInterface*)yIF)->integers[i].bitcnt -1 + i);
	i+=((yInterface*)yIF)->integers[i].bitcnt;
	i--;
	//add for writing to string
	toWrite += "  <variable name=\"" + noBinName + "\" index=\"" + vIndex + "\" value=\"" + std::to_string(res) + "\" block=\"" + std::to_string(block[i]) + "\"/>\n";
      } 
    } else{
	//round value if i is of type binary
      if(type[i]==BINARY){
	if(block[i]==1)
          toWrite += "  <variable name=\"" + ((yInterface*)yIF)->integers[i].name + "\" index=\"" + std::to_string(i) + "\" value=\"" + std::to_string((int)floor(getFirstStageSolutionValue(i)+0.5)) + "\" block=\"" + std::to_string(block[i]) + "\"/>\n";
	else 
	  toWrite += "  <variable name=\"" + ((yInterface*)yIF)->integers[i].name + "\" index=\"" + std::to_string(i) + "\" value=\"" + std::to_string((int)floor(PV[0][i]+0.5)) + "\" block=\"" + std::to_string(block[i]) + "\"/>\n";
      }else{
        toWrite += "  <variable name=\"" + ((yInterface*)yIF)->integers[i].name + "\" index=\"" + std::to_string(i) + "\" value=\"" + std::to_string(PV[0][i]) + "\" block=\"" + std::to_string(block[i]) + "\"/>\n";
      }
    }
  //toWrite += " <variable name=\"" + noBinName + "\" index=\"" + std::to_string(i) + "\" value=\"" + std::to_string(PV[0][i]) + "\" block=\"" + std::to_string(block[i]) + "\"/>\n";
  }
  toWrite += " </variables>\n";
  toWrite += "</YasolSolution>";
  Sol.solprint(toWrite,Final.c_str());
}
  
/*void QBPSolver::WriteSolutionFile(coef_t value, double gap,std::string Status){ 
//cerr<<"Write Outputfile"<<endl;
//cerr <<"OFFSET IS" << finalOffset<<endl;
	if(abs(value)==abs(n_infinity)) Status="INFEASIBLE (universal player can force violation of existential constraint system)";
	else if (abs(value)==(AllInfeasible)) Status="FEASIBLE (existential player can force violation of universal constraint system)";
	if(objInverted) value*=-1;
	//creating CommPrint object to use solprint for solution file
	//setting path for solution file
	CommPrint Sol;
	
//	Sol.solFilename =  getInputFileName()+".sol";
	string Final(getInputFileName()+".sol");
        const char *Start=Final.c_str();
	ifstream myfile (Start);
	int Counter=0;
        while (myfile.is_open()){
	  myfile.close();	
	  Counter++;
	  string TmpA(getInputFileName()+".sol");
	  Final=TmpA+std::to_string(Counter);
	  const char *tmp=Final.c_str();
	  //ifstream 
	  myfile.open(tmp);
	  if(myfile.is_open()) continue;
	}	
	cerr << "solution file written to " << Final << endl;



	//toWrite is a string that is appended to the solution file
	//this string is written in (pseudo)XML
	//the write stage is divided in stages for header+quality, linearConstraints and variables
	//values have to be added between \". example: \"value\"
	string toWrite = "";
	toWrite += "<?xml version = \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
	toWrite += "<YasolSolution version=\"1\">\n";
	toWrite += " <header\n";
	toWrite += "   ProblemName=\""+ getInputFileName()+ "\"\n";
	toWrite += "   SolutionName=\""+Final+"\"\n";
	toWrite += "   ObjectiveValue=\""+std::to_string((double)value+finalOffset)+"\"\n";
	//toWrite += "   ObjectiveOffset=\""+std::to_string(finalOffset)+"\"\n";
	//toWrite += "   solutionStatusString=\""+Status+"\"\n";
	//toWrite += "   SolutionMethodString=\"qmip\"\n";
	toWrite += "   Runtime=\""+std::to_string((int)time(NULL)-ini_time)+" seconds\"\n";
	toWrite += "   DecisionNodes=\""+std::to_string((int)getNumberOfDecisions())+"\"\n";
	toWrite += "   PropagationSteps=\""+std::to_string((int)getNumberOfPropagationSteps())+"\"\n";
	toWrite += "   LearntConstraints=\""+std::to_string((int)getNumberOfLearntConstraints())+"\"/>\n";
	toWrite += " <quality\n";
	toWrite += "   SolutionStatus=\""+Status+"\"\n";
	toWrite += "   Gap=\""+std::to_string(gap)+"\"/>\n";
	Sol.solprint(toWrite,Final.c_str());

	//reset toWrite because method solprint is in append mode
	toWrite = " <linearConstraints>\n";
	toWrite += " </linearConstraints>\n";
	Sol.solprint(toWrite,Final.c_str());

	//reset toWrite because method solprint is in append mode
	toWrite = " <variables>\n";
	for(int i=0; i<nVars(); i++){




		//loop from QBPSolver_ctrl.cpp
		if (((yInterface*)yIF)->integers[i].bitcnt > 1) {

				//changed by MH and FB
				//get real name without underscore from binarized variables

				//create char array from name string
				const char *name = ((yInterface*)yIF)->integers[i].name.c_str();

				//positon of underscore in char array
				int bpoint = 0;

				//get position of underscore
				for(int j=0; j<=((yInterface*)yIF)->integers[i].name.length(); j++){
					if(name[j]=='_'){
						bpoint = j;
					}
				}

				//create new name for binarized variable with substring
				std::string noBinName = ((yInterface*)yIF)->integers[i].name.substr(0,bpoint);

				//change end

			//change to real value
			if (((yInterface*)yIF)->integers[i].pt2leader == i) {
				int res = 0;
				for (int z = 0;z < ((yInterface*)yIF)->integers[i].bitcnt;z++) {
					if(block[i] == 1){
						res = 2*res + getFirstStageSolutionValue(i+z);
					}else {
						res = 2*res + PV[0][i+z];
					}
					

				}


				//index of binarized variable
				std::string vIndex = std::to_string(i) + "-" + std::to_string(((yInterface*)yIF)->integers[i].bitcnt -1 + i);


				i+=((yInterface*)yIF)->integers[i].bitcnt;
				i--;

				
				//add for writing to string
				toWrite += "  <variable name=\"" + noBinName + "\" index=\"" + vIndex + "\" value=\"" + std::to_string(res) + "\" block=\"" + std::to_string(block[i]) + "\"/>\n";

			} 
		} else{

			//round value if i is of type binary
			if(type[i]==BINARY){
				toWrite += "  <variable name=\"" + ((yInterface*)yIF)->integers[i].name + "\" index=\"" + std::to_string(i) + "\" value=\"" + std::to_string((int)floor(PV[0][i]+0.5)) + "\" block=\"" + std::to_string(block[i]) + "\"/>\n";
			}else{
				toWrite += "  <variable name=\"" + ((yInterface*)yIF)->integers[i].name + "\" index=\"" + std::to_string(i) + "\" value=\"" + std::to_string(PV[0][i]) + "\" block=\"" + std::to_string(block[i]) + "\"/>\n";
			}

			
		}

		

		
		

		//toWrite += " <variable name=\"" + noBinName + "\" index=\"" + std::to_string(i) + "\" value=\"" + std::to_string(PV[0][i]) + "\" block=\"" + std::to_string(block[i]) + "\"/>\n";
	}
	toWrite += " </variables>\n";
        toWrite += "</YasolSolution>";
	Sol.solprint(toWrite,Final.c_str());
}*/

void QBPSolver::BuildLegalScenario(){	//Only for Polyhedral Uncertainty Set; NO INTERDEPENDENCE, since best case should no longer be consdidered for existential variables in universal constraints
	//cerr << "Call Build Scenario" <<endl;
	if(!UniversalConstraintsExist) return;
	for(int AllStage=0;AllStage<PermutationOfAllStages.size();AllStage++){
		  	for(int in=0;in<PermutationOfAllStages[AllStage].size();in++){
		  		SparseScenario[PermutationOfAllStages[AllStage][in]]=extbool_Undef;
			  	ScenarioProp[PermutationOfAllStages[AllStage][in]]=extbool_Undef;
			}
		}
	for (int j=0;j<ALLconstraints.size();j++){
		Constraint &c =ALLconstraintallocator[ALLconstraints[j]];
		c.header.localbest=c.header.btch1.best;
	}

	if(!getIsSimplyRestricted()){
	    for (int v=0;v<VarsPresentInAllConstraints.size();v++){
	    	int i=VarsPresentInAllConstraints[v];
		if(VarsInAllConstraints[i].size()==0) continue; //Should be deprecated due to new datastructure VarsPresentInAllConstraints
		SetBoundsInExternalIPSolver(AllSolver,i);
        	//if(decisionLevel()==1 && block[i]==1) continue;
        	/*if (assigns[i]==extbool_Undef || 
                	((decisionLevel()==1&& block[i]>1) || (block[i]> block[trail[getDecisionIndexInTrail(decisionLevel())]] ))){
                    AllSolver->setVarLB(i,0);
                     AllSolver->setVarUB(i,1);
                }
                else if (assigns[i]==1){ AllSolver->setVarLB(i,1);AllSolver->setVarUB(i,1);}
                else if (assigns[i]==0){ AllSolver->setVarLB(i,0);AllSolver->setVarUB(i,0);}
                */
            }
            AllSolver->solve();
            //assert(AllSolver->getSolutionStatus()==extSol::QpExternSolver::QpExtSolSolutionStatus::OPTIMAL);
            if(AllSolver->getSolutionStatus()!=extSol::QpExternSolver::QpExtSolSolutionStatus::OPTIMAL){
		if(getShowWarning()) cerr << "Warning: IP of universal constraint System is infeasible when trying to find a valid scenario" <<endl;
		return;
            }
            std::vector<data::QpNum> AllSol;
	    AllSolver->getValues(AllSol);
	    for (int i=0;i<nVars();i++)
    	        if(eas[i]==UNIV){ 
		    if(UniversalPolytope || (!UniversalPolytope && VarsInAllConstraints[i].size()==0))
		        SparseScenario[i]=(int)AllSol[i].asDouble();     
	        }
	}
	else if(/*!getIsSimplyRestricted() ||*/!UniversalPolytope){
		//i.e. existential variables are present within some universal constraints
		if(AllPropQStore.empty()) return;
		for (int bl=0;bl < AllPropQStore.size();bl++){
			for (int i=0;i<AllPropQStore[bl].size();i++){
				int va=AllPropQStore[bl][i].first;
				int val=AllPropQStore[bl][i].second;
				if(ScenarioProp[va]==extbool_Undef && assigns[va]==extbool_Undef){
	          		/*bool legal=false;
	          		bool Flipped =false;
	          		while(!legal){
	          			legal=true;
		          		for (int j=0;j<VarsInAllConstraints[va].size();j++){
		          			Constraint &c = ALLconstraintallocator[VarsInAllConstraints[va][j].cr];
		          		    int pos = VarsInAllConstraints[va][j].pos;
			          		if(c.saveInfeasForAllScenario(pos, val, va)){
				  				cerr << "Error: Propagated Universal Variable cannot be set this way"
				  				assert(0);
				  				//legal=false;
				  				//val=1-val;
				  				//assert(!Flipped);
				  				//Flipped=true;
				  				//break;
				  			}
				  		}
				  	}*/
				  	//cerr << "Try Setting x_" << va <<"="<<val << endl;
				  	SparseScenario[va]=val;
				  	PropagateForScenario(va,val);
			  	}
			}
		}
		for (int va=0;va<ScenarioProp.size();va++){
			if(assigns[va]!=extbool_Undef)  SparseScenario[va]=assigns[va];
			else if(ScenarioProp[va]!=extbool_Undef) {
		  		SparseScenario[va]=ScenarioProp[va];
		  		ScenarioProp[va]=extbool_Undef;
		  	}
		}
	}
	else{
		assert(UniversalPolytope);
		CreateAndUpdateShuffle();

		//Already initialized up top
		/*for (int j=0;j<ALLconstraints.size();j++){
			Constraint &c =ALLconstraintallocator[ALLconstraints[j]];
			c.header.localbest=c.header.btch1.best;
		}

		//More efficient here! Not needed????????
		for(int AllStage=0;AllStage<PermutationOfAllStages.size();AllStage++){
		  	for(int in=0;in<PermutationOfAllStages[AllStage].size();in++){
		  		SparseScenario[PermutationOfAllStages[AllStage][in]]=extbool_Undef;
			  	ScenarioProp[PermutationOfAllStages[AllStage][in]]=extbool_Undef;
			}
		}*/

		for(int AllStage=0;AllStage<PermutationOfAllStages.size();AllStage++){
		  	for(int in=0;in<PermutationOfAllStages[AllStage].size();in++){
		  		int va=PermutationOfAllStages[AllStage][in];
				bool PropOK=true;
		  		if(ScenarioProp[va]==extbool_Undef && assigns[va]==extbool_Undef){
		  			int val;
		  			if (killer[va] == 0 || killer[va]==1)
	            		val = killer[va];
	          		else val = (p_activity[va]<n_activity[va] ? 0 : 1);
	          		bool legal=false;
	          		bool Flipped =false;
				while(!legal){
	          			legal=true;
		          		for (int j=0;j<VarsInAllConstraints[va].size();j++){
		          			Constraint &c = ALLconstraintallocator[VarsInAllConstraints[va][j].cr];
		          		    int pos = VarsInAllConstraints[va][j].pos;
			          		if(c.saveInfeasForAllScenario(pos, val, va)){
				  				legal=false;
				  				val=1-val;
				  				assert(!Flipped);
				  				Flipped=true;
				  				break;
				  			}
				  		}
				  	}
				  	//cerr << "Try Setting x_" << va <<"="<<val << endl;
				  	SparseScenario[va]=val;
				  	PropOK=PropagateForScenario(va,val);
				}
			  	//else if(assigns[va]!=extbool_Undef)  SparseScenario[va]=assigns[va];
				else if(assigns[va]!=extbool_Undef){
			  		//cerr << "Is set x_" << va <<"="<<(int) assigns[va] << endl;
			  		SparseScenario[va]=assigns[va];
			  		ScenarioProp[va]=extbool_Undef;	
					PropOK=PropagateForScenario(va,assigns[va]);
			  	} 				  	
				else{
		                            //cerr << "Is propped x_" << va <<"="<<(int) ScenarioProp[va] << endl;

			  		//If already propageted in "PropagateForScenario"  
			  		SparseScenario[va]=ScenarioProp[va];
			  		ScenarioProp[va]=extbool_Undef;
					PropOK=PropagateForScenario(va,SparseScenario[va]);
			  	}
				if(!PropOK){
                                    	//cerr<<"PropagateForScenario failed" <<endl;
                                        for(int AllStage=0;AllStage<PermutationOfAllStages.size();AllStage++){
                                    	    for(int in=0;in<PermutationOfAllStages[AllStage].size();in++){
                                            	SparseScenario[PermutationOfAllStages[AllStage][in]]=extbool_Undef;
                                            	ScenarioProp[PermutationOfAllStages[AllStage][in]]=extbool_Undef;
                                    	    }
                                	}
//                                        cerr <<"Will return"<<endl;
                                        return;                                 
                                 }

		  	}
		}
	}
}

void QBPSolver::AdaptConstraint( ca_vec<CoeVar>& ps,bool ClearIt, bool Print){
if(ClearIt){ps.clear();
return;
}

//if(!UniversalConstraintsExist || UniversalPolytope) return;
    if(!UniversalConstraintsExist || !UniversalMultiBlockConstraints) return;
    bool Debug =false;//true;//Print;

    if(Debug){
	cerr << "Input Constraint: " << endl;
	int rhsV=1;
	double LHSV=0;
	double OptLHSV=0;
	for (int i=0;i<ps.size();i++){
	   cerr << (sign(ps[i])?" -":" +") << ps[i].coef <<"x_" << var(ps[i])<< "("<<(int)assigns[var(ps[i])] << "," <<getFixed(var(ps[i]))<<")";
	   rhsV-=sign(ps[i]);
	   if (/*assigns[var(ps[i])]!=extbool_Undef &&*/sign(ps[i])){
		LHSV-=ps[i].coef*assigns[var(ps[i])];
	       	if(USE_TRACKER) OptLHSV-=ps[i].coef*optSol[var(ps[i])];
	   }
	   else/* if (assigns[var(ps[i])]!=extbool_Undef)*/{
  		LHSV+=ps[i].coef*assigns[var(ps[i])];
		if(USE_TRACKER) OptLHSV+=ps[i].coef*optSol[var(ps[i])];
	    }
	}
	cerr << " >= " << rhsV << endl;
	cerr << "Current assigns: "<<LHSV << ">=?" << rhsV << endl;
	if(USE_TRACKER) cerr << "optsol: " << OptLHSV << ">=?" << rhsV << endl;
        if(USE_TRACKER && OptLHSV<rhsV){
		cerr << "ACHTUNG!!!!!" << endl;
	}
    }
    // Collect Data of Constraint: Max- & Min-Blocks of Existential and Universal Variables
    int MaxBlockE=0;
    int MinBlockE=maxBlock;
    int MaxBlockA=0;
    int MinBlockA=maxBlock;
    VarInCutSeen.clear();
    AllConstraintsSeen.clear();
    bool UniversalVariablesPresent=false;
    //if(Debug) cerr << "VarInCut.size() is " << VarInCut.size() << " and nVars() is " << nVars() << endl;
    assert(VarInCut.size()>=nVars());
    for (int i=0;i<ps.size();i++){
	//if(Debug)cerr << i <<" var(ps[i]) is " << var(ps[i]) << " and VarInCut[var(ps[i])] is "<< VarInCut[var(ps[i])] << endl;
    	if(VarInCut[var(ps[i])]>=0 && sign(ps[i])!=sign(ps[VarInCut[var(ps[i])]]) ) {if(getShowWarning()) cerr << "Warning: Variable x_" << var(ps[i]) << " already in Cut at position " << (int)VarInCut[var(ps[i])]<< " with different sign!" << endl;}
	else{ 	
	    VarInCutSeen.push(var(ps[i]));
	    VarInCut[var(ps[i])]=i;
	}
	if(eas[var(ps[i])]==EXIST /*&& !isFixed(var(ps[i]))*/){
	    if(block[var(ps[i])]>MaxBlockE) MaxBlockE=block[var(ps[i])];
	    if(block[var(ps[i])]<MinBlockE) MinBlockE=block[var(ps[i])];
	}else if(eas[var(ps[i])]==UNIV){
	    UniversalVariablesPresent=true;
	    if(block[var(ps[i])]>MaxBlockA) MaxBlockA=block[var(ps[i])];
      	    if(block[var(ps[i])]<MinBlockA) MinBlockA=block[var(ps[i])];
	}
    }
    if(UniversalVariablesPresent){
	ca_vec<CRef> ALLconstraint_seen_stack;
	ALLconstraint_seen_stack.clear();
	//if(/*?*/MinBlockA < MaxBlockE || MinBlockE<MaxBlockA){ //If existential variables with lower block than some universal variables are present
	if(1||MinBlockE<MaxBlockA){	
	    if(ClearIt){ps.clear();
		goto CleanUp;}
	    if(Debug){
		cerr << "Existential Variables from blocks " << MinBlockE << "-" << MaxBlockE << endl;
		cerr << "Universal Variables from blocks " << MinBlockA << "-" << MaxBlockA << endl;
	    }
	    for (int i=0;i<ps.size();i++){
// For each universal variable present in learntConstraint check whether there is a universal constraint that might affect its assignment
                if(eas[var(ps[i])]==UNIV){
		    if(Debug) cerr << "Universal variable x_" << var(ps[i]) <<" appears in " <<VarsInAllConstraints[var(ps[i])].size() << " universal constraints!"<<  endl;
		    for (int j=0;j<VarsInAllConstraints[var(ps[i])].size();j++){
			CRef cr = VarsInAllConstraints[var(ps[i])][j].cr;
			if(ALLconstraintallocator[cr].header.mark==0){
			    if(Debug)cerr << "Untersuche diese Constraint " << endl;
			    ALLconstraintallocator[cr].header.mark=1;
			    AllConstraintsSeen.push(cr);
        		    Constraint &c = ALLconstraintallocator[cr];
		            int sig = sign(c[VarsInAllConstraints[var(ps[i])][j].pos]);
			    if(sig!=assigns[var(ps[i])] && assigns[var(ps[i])]!=extbool_Undef) continue; //LOOKHERE
		//	    if(assigns[var(ps[i])]==extbool_Undef) continue;                             //LOOKHERE
			    // Add information of the universal constraint into the learntConstraint
			    for (int k=0;k<c.size();k++){
				if(eas[var(c[k])]!=UNIV||block[var(c[k])] <= MaxBlockE /*?*/ /* ||assigns[var(c[k])]!=extbool_Undef*/ ){
				    if(Debug)cerr<< "Of Variable  x_" << var(c[k]) << "AllpropVals are " <<   AllpropQlimiter[(var(c[k])<<1)|1] << " and " << AllpropQlimiter[(var(c[k])<<1)|0] << endl;
				    if( VarInCut[var(c[k])]>=0 && sign(ps[VarInCut[var(c[k])]]) == assigns[var(c[k])] )
				    {	
					if(Debug)cerr << "Variable x_" << var(c[k]) <<"="<<(int)assigns[var(c[k])] << "already set appropriately in Cut at position " << VarInCut[var(c[k])] << endl;
				        if(sign(ps[VarInCut[var(c[k])]])!=assigns[var(c[k])]){
					    //if(Debug || info_level>=-1)
					    if(getShowWarning()) cerr << "WARNING: Universal variable variable already present in cut; but set the other way.. " << endl;
					}
				    }
				    else if(VarInCut[var(c[k])]>=0 && sign(ps[VarInCut[var(c[k])]]) != assigns[var(c[k])]){
					if(getShowWarning()) cerr << "WARNING: Universal variable variable already present in cut; but set the other way.. " << endl;
				    }
				    //{
				    else{
					if(Debug)cerr << "Variable x_" << var(c[k]) <<"=" <<(int)assigns[var(c[k])]<< " not yet in Cut! " << endl;
					if((/*1||*/assigns[var(c[k])]!=sign(c[k])) && assigns[var(c[k])]!=extbool_Undef){
					    VarInCutSeen.push(var(c[k]));
					    VarInCut[var(c[k])]=ps.size(); 
					    CoeVar q = mkCoeVar(var(c[k]),1.0,assigns[var(c[k])]);
					    if (Debug) cerr << "Added x_"<<var(c[k]) <<"="<<(int)assigns[var(c[k])] << " to the cut " << endl;
					    // ? Helpful varBumpActivity(var(c[k]), 1-assigns[var(c[k])],0);
                            		    ps.push(q);
					}
					else if(0){
					    VarInCutSeen.push(var(c[k]));
                                            VarInCut[var(c[k])]=ps.size(); 
                                            CoeVar q = mkCoeVar(var(c[k]),1.0,sign(ps[VarInCut[var(c[k])]]));
					    if (Debug) cerr << "Added unsassigned  x_"<<var(c[k]) <<" to the cut " << endl;
					    ps.push(q);
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
    else{ 
	if(Debug) cerr << "No universal variables present; nothing to do" << endl;
    }
CleanUp:
    while(AllConstraintsSeen.size() > 0) {
	ALLconstraintallocator[AllConstraintsSeen.last()].header.mark = 0;
        AllConstraintsSeen.pop();
    }
    while(VarInCutSeen.size() > 0) {
	//if(Debug) cerr << "Deleted Variable x_" << (int) VarInCutSeen.last() << " in VariInCut!" << endl;
        VarInCut[VarInCutSeen.last()] = -1;
        VarInCutSeen.pop();
    }
if(Debug){
        cerr << "Output Constraint: " << endl;
        int rhsV=1;
        double LHSV=0;
        double OptLHSV=0;
        for (int i=0;i<ps.size();i++){
           cerr << (sign(ps[i])?" -":" +") << ps[i].coef <<"x_" << var(ps[i])<< "("<<(int)assigns[var(ps[i])] << "," <<getFixed(var(ps[i]));
           if(USE_TRACKER) cerr <<","<<optSol[var(ps[i])];
           cerr <<",B"<<block[var(ps[i])] <<  ")";
           rhsV-=sign(ps[i]);
           if (/*assigns[var(ps[i])]!=extbool_Undef &&*/sign(ps[i])){ 
                LHSV-=ps[i].coef*assigns[var(ps[i])];
                if(USE_TRACKER)OptLHSV-=ps[i].coef*optSol[var(ps[i])];
           }
           else/* if (assigns[var(ps[i])]!=extbool_Undef)*/{  
                LHSV+=ps[i].coef*assigns[var(ps[i])];
                if(USE_TRACKER) OptLHSV+=ps[i].coef*optSol[var(ps[i])];
           }
	}
        cerr << " >= " << rhsV << endl;
        cerr << "Current assigns: "<<LHSV << ">=?" << rhsV << endl;
	if(USE_TRACKER) cerr << "optsol: " << OptLHSV << ">=?" << rhsV << endl;
        if(USE_TRACKER && OptLHSV<rhsV) assert(0);// cerr<< "ERROR: CUTS OFF OPTIMAL" << endl;
    }
}


bool QBPSolver::PropagateForScenario(int va, extbool val){
	//cerr <<"Start PropagateForScenario x_" <<va << "="<<(int)val <<endl;
	if(!UniversalConstraintsExist) return false;
		//assert(UniversalPolytope);

	for (int i=0; i < VarsInAllConstraints[va].size();i++) {
        int pos = VarsInAllConstraints[va][i].pos;
        CRef cr = VarsInAllConstraints[va][i].cr;
        Constraint &c = ALLconstraintallocator[cr];
        coef_t rhs = c.header.rhs;
        assert(c.size()>0);
		if(0){
	        	cerr <<"Print Constraint:";
	        for (int j=0;j<c.size();j++){
	        	if(sign(c[j]))cerr<<-c[j].coef << "x_"<< var(c[j]) << "(";
	        	else cerr<<c[j].coef << "x_"<< var(c[j]) << "(";
	        	cerr << (int) assigns[var(c[j])] << "/"<<(int)SparseScenario[var(c[j])] <<"/" <<(int)ScenarioProp[var(c[j])]<<") ";
	        	
	        }
	        cerr << " >= " << c.header.rhs << endl;
	        cerr << "local best: " << c.header.localbest << endl; 
	        cerr << "Largest: " << c.header.largest << endl;
	    }
        //massert(ConstraintIsWellFormed(c));
        coef_t coef = c[pos].coef;
        massert(coef >= (coef_t)0);
        if (/*1||!c.header.isSat*/ assigns[va]==extbool_Undef) {
            if (sign(c[pos])) {
                if (val==1)  c.header.localbest -= coef;
            } else {
                if (val==0) c.header.localbest -= coef;
            }
            massert(c.header.localbest >= c.header.rhs);
            /*if (c.header.wtch2.worst >= c.header.rhs && cr != constraints[0]) {
                SwapOut(va,c);
                continue;
            }*/
            if (!( (sign(c[pos]) && val==1) || (!sign(c[pos]) && val==0)  )) continue;
            //if (c.header.btch1.best <= c.header.rhs-LP_EPS) continue;
            int l=c.header.largest;
            while (l < c.size() && c.header.localbest - c[l].coef < c.header.rhs) {
                    /*    cerr << "l"<<l<<endl;
                        cerr << "var(c[l])="<<var(c[l])<<endl;
                        cerr << "ScenarioProp[var(c[l])]="<<(int)ScenarioProp[var(c[l])]<<endl;
                        cerr << "1-(c[l].x&1)="<<(int)1-(c[l].x&1)<<endl;*/
                    if ( var(c[l]) != va &&assigns[ var(c[l]) ] == extbool_Undef && type[var(c[l])] == BINARY &&SparseScenario[var(c[l]) ]==extbool_Undef) {
             	    //cerr <<"Propagate x_" << var(c[l]) <<"="<<1-(c[l].x&1)<<endl;
                    //cerr <<(int)SparseScenario[var(c[l]) ] << " and "<<  (int)ScenarioProp[var(c[l]) ]<<endl;
                    if(!(ScenarioProp[var(c[l]) ]==extbool_Undef || ScenarioProp[var(c[l]) ]==1-(c[l].x&1))){/*cerr << "Ret1"<<endl;*/ return false;}
                    if(!(SparseScenario[var(c[l]) ]==extbool_Undef || SparseScenario[var(c[l]) ]==1-(c[l].x&1))){/*cerr <<"Ret2"<<endl;*/ return false;}
                    if (SparseScenario[var(c[l]) ]==extbool_Undef &&  ScenarioProp[var(c[l]) ]==extbool_Undef){
			ScenarioProp[var(c[l]) ]=1-(c[l].x&1);
                    	//cerr << "Added in SparseSce setting "<< var(c[l]) << " to " << (1-(c[l].x&1))<< endl;
                    	//if(!PropagateForScenario(var(c[l]), (1-(c[l].x&1)))){cerr << "deep Prop Failed: x_"<<var(c[l])<<"="<<(int)(1-(c[l].x&1) <<endl;  return false;}
		    }
		    else if ( ScenarioProp[var(c[l]) ] != extbool_Undef &&  ScenarioProp[var(c[l]) ]!=1-(c[l].x&1)){
		/*	cerr << "l"<<l<<endl;
			cerr << "var(c[l])="<<var(c[l])<<endl;
			cerr << "ScenarioProp[var(c[l])]="<<(int)ScenarioProp[var(c[l])]<<endl;
			cerr << "1-(c[l].x&1)="<<(int)1-(c[l].x&1)<<endl;*/
		    	if(getShowWarning()) cerr <<"WARNING: Propgation for finding a scenario failed; will return" <<endl;
			return false;
		    }
		}
				l++;
            }
        } 
    }
	return true;
}

    void QBPSolver::KeepAllCleanUnassign(int var, extbool val){
    	//return;
    	if(!UniversalConstraintsExist) return;
		for (int i=0; i < VarsInAllConstraints[var].size();i++) {
	        Constraint &c = ALLconstraintallocator[VarsInAllConstraints[var][i].cr];
	        massert(!c.header.deleted);
	        //massert(ConstraintIsWellFormed(c));
	        if (1||!c.header.isSat) {
	            c.header.btch1.best = VarsInAllConstraints[var][i].btch1.best;
	            c.header.wtch2.worst = VarsInAllConstraints[var][i].wtch2.worst;
	            if (VarsInAllConstraints[var][i].pos < c.header.largest) c.header.largest = VarsInAllConstraints[var][i].pos;
	            //massert(ConstraintIsWellFormed(c));
	        } 
	        /*else {
	            c.header.btch1.watch1 = VarsInConstraints[var][i].btch1.watch1;
	            c.header.wtch2.watch2 = VarsInConstraints[var][i].wtch2.watch2;
	            //if (VarsInConstraints[var][i].pos < c.header.largest) c.header.largest = VarsInConstraints[var][i].pos;
	            //massert(ConstraintIsWellFormed(c));
	        }*/
    	}

    	if(AllpropQlimiter[(var<<1)|(1-val)]<-nVars()){
    		// If this variable was propagated due to variable AllpropQlimiter[(var<<1)|(1-val)]+nVars()
    		if(PRINT_PROP) cerr <<"Reenter setting " << var << " to " << (int)val << " Current AllPropQLimiter: " <<AllpropQlimiter[(var<<1)|(1-val)] << " and future: " << AllpropQlimiter[(var<<1)|(1-val)]+nVars() << endl; 
    		AllPropQStore[block[var]-1].push_back( make_pair(var, val));
			AllpropQlimiter[(var<<1)|(1-val)]=AllpropQlimiter[(var<<1)|(1-val)]+nVars();
    	}
			if(AllVarsProppedDueTo.size()<nVars())     AllVarsProppedDueTo= vector<vector<int>>(nVars(), vector<int>(0));

    	while(AllVarsProppedDueTo[var].size()>0){
    		//Delete All those variables that have been propagated due to var itself
		 	if(PRINT_PROP) cerr << "Deleted, that x_" << (AllVarsProppedDueTo[var].back()>>1) << " is propagated due to " << var << endl;

			//AllpropQlimiter[AllVarsProppedDueTo[var].back()]=0;//-2*nVars()-1;
    		int Found =false;
    		for (int i=0;i< AllPropQStore[block[AllVarsProppedDueTo[var].back()>>1]-1].size();i++){
    			if(AllPropQStore[block[AllVarsProppedDueTo[var].back()>>1]-1][i].first== (AllVarsProppedDueTo[var].back()>>1)){
    				if(PRINT_PROP) cerr << "Found and deleted in AllPropQStore" << endl;
					AllPropQStore[block[AllVarsProppedDueTo[var].back()>>1]-1][i].second=AllPropQStore[block[AllVarsProppedDueTo[var].back()>>1]-1].back().second;
    				AllPropQStore[block[AllVarsProppedDueTo[var].back()>>1]-1][i].first=AllPropQStore[block[AllVarsProppedDueTo[var].back()>>1]-1].back().first;
    				AllPropQStore[block[AllVarsProppedDueTo[var].back()>>1]-1].pop_back();
				 	AllpropQlimiter[AllVarsProppedDueTo[var].back()]=0; //Here not the variable, but variable with orientation!
    				Found=true;
    				//break;
    			}
    		}
    		if(!Found && getShowWarning()) cerr << "WARNING: Variable not Found in AllPropQStore" << endl;
    		AllVarsProppedDueTo[var].pop_back();
    	}

    	//for (int i=0;i<AllPropQStore.size();i++){
		//	AllPropQStore[i].clear();
		//}
		EmptyAllPropQ();

    }
    void QBPSolver::KeepAllClean(int va, extbool val){
    	//cerr << "STart KeepClean" << endl;
		
		//EmptyAllPropQ();
		//return;
		if(!UniversalConstraintsExist) return;
    	for (int i=0; i < VarsInAllConstraints[va].size();i++) {
	        int pos = VarsInAllConstraints[va][i].pos;
	        CRef cr = VarsInAllConstraints[va][i].cr;
	        Constraint &c = ALLconstraintallocator[cr];
	        coef_t rhs = c.header.rhs;
	        
	        //massert(ConstraintIsWellFormed(c));
	        coef_t coef = c[pos].coef;
	        massert(coef >= (coef_t)0);
	        if (1||!c.header.isSat) {
	        	//Is Clique?
	            VarsInAllConstraints[va][i].btch1.best = c.header.btch1.best;
	            VarsInAllConstraints[va][i].wtch2.worst = c.header.wtch2.worst;
	            if (sign(c[pos])) {
	                if (val==0 ) c.header.wtch2.worst += coef;
	                if (val==1)  c.header.btch1.best -= coef;
	            } else {
	                if (val==0) c.header.btch1.best -= coef;
	                if (val==1) c.header.wtch2.worst += coef;
	            }
	            massert(c.header.btch1.best >= c.header.rhs);
	            /*if (c.header.wtch2.worst >= c.header.rhs && cr != constraints[0]) {
	                SwapOut(va,c);
	                continue;
	            }*/
	            //if (!( (sign(c[pos]) && val==1) || (!sign(c[pos]) && val==0)  )) continue;
	            //if (c.header.btch1.best <= c.header.rhs-LP_EPS) continue;
	            while ( c.header.largest < c.size() && assigns[ var(c[c.header.largest]) ] != extbool_Undef) c.header.largest++;
	            int l=c.header.largest;
	            while (l < c.size() && c.header.btch1.best - c[l].coef < c.header.rhs) {
	                if (eas[var(c[l])]==UNIV && assigns[ var(c[l]) ] == extbool_Undef && type[var(c[l])] == BINARY) {
	                	//cerr <<"Folgere x_"<<var(c[l])<<"="<< (1-(c[l].x&1))<<endl;
	                    massert(VarsInAllConstraints[va][i].cr != CRef_Undef);
	                    //constraintBumpActivity(constraintallocator[VarsInConstraints[va][i].cr]);
	                    if (AllpropQlimiter[c[l].x] == 0) { //==0!!! <0 is interpreted as variable which forced the progation
		 	 if(PRINT_PROP)cerr << "Current PropQSize= " << AllPropQ.size() << endl;          
	          	AllPropQ.push_back(std::make_pair(var(c[l]),1-(c[l].x&1)));
	                        AllpropQlimiter[c[l].x] = -va-1;//HERE ERROR: AllPropQ.size();		
       if(PRINT_PROP) cerr << "Added prop var " << var(c[l]) << " to " << 1-(c[l].x&1) << "Limiter: " <<AllpropQlimiter[c[l].x] <<endl;
		
	                        NumAllPropsPush++;
	                        //cerr << "Setze PropQ Eintrag " << (AllPropQ.size()-1) << " auf " << var(c[l]) << "/" <<(1-(c[l].x&1))<< endl;
	                        //cerr << "Setze Limiter Eintrag " << c[l].x << " auf " << AllPropQ.size() << endl;
	                        //cerr << "Dort steht " << AllPropQ[AllpropQlimiter[c[l].x]-1].first<< " " <<  AllPropQ[AllpropQlimiter[c[l].x]-1].second<<endl;
	                        if (AllpropQlimiter[c[l].x^1] > 0) {
	                        	//cerr << (c[l].x^1) << " Limiter Entry of Opposite: " << AllpropQlimiter[c[l].x^1] << endl;
	                        	//cerr << "AllPQ " << AllPropQ[AllpropQlimiter[c[l].x^1]-1].first<< " " <<  AllPropQ[AllpropQlimiter[c[l].x^1]-1].second<<endl;
	                        	//cerr <<"AND " <<var(c[l])<< " " << (1-(c[l].x&1)) << endl;
		                        /*Assertion Quatsch!*/ assert(AllPropQ[AllpropQlimiter[c[l].x^1]-1].first==var(c[l]));
	                        	if(PRINT_PROP) cerr<<"Warning: conflict in Universal Variables" <<endl;
	                        	//cerr << "check entry "<< var(c[l]) << " and " <<1-(c[l].x&1) <<"-"<<(c[l].x^1)<<endl;

	                       
	                            /*conflict = true;
	                            ValueConstraintPair tmp1=propQ[propQlimiter[c[l].x  ]-1];
	                            ValueConstraintPair tmp2=propQ[propQlimiter[c[l].x^1]-1];
	                            EmptyPropQ();
	                            PROPQ_PUSH(va,val,tmp1);
	                            propQlimiter[c[l].x  ] = 1;
	                            propQlimiter[c[l].x^1] = 2;
	                            PROPQ_PUSH(va,val,tmp2);
	                            ix1 = propQlimiter[c[l].x  ]-1;
	                            ix2 = propQlimiter[c[l].x^1]-1;
	                            //cerr <<"2";
	                            break;*/
	                        }
	                    //HIER: WAS MACHT DAS ELSE???
	                    } //else AllPropQ[AllpropQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l);
	                }
	                l++;
	            }
	        } 
	        /*else if (c.header.btch1.watch1 > -2) {
	            VarsInConstraints[va][i].btch1.watch1 = c.header.btch1.watch1;
	            VarsInConstraints[va][i].wtch2.watch2 = c.header.wtch2.watch2;
	            if ((sign(c[pos]) && val == 0) || (!sign(c[pos]) && val == 1)) {
	                // SAT Klausel ist erf�llt
	                SATswapOut(va,c);
	                c.header.btch1.watch1 = c.header.wtch2.watch2 = -2;
	                continue;
	            } else {
	                // SAT Klausel ist nicht erf�llt
	                // FindNewSatWatcher(va,c);
	                assert(va == var(c[c.header.btch1.watch1]) || va == var(c[c.header.wtch2.watch2]));
	                bool isSatisfied = false;
	                bool newWatcherFound = false;
	                int start_i = (c.header.wtch2.watch2 < c.header.btch1.watch1 ? c.header.btch1.watch1+1 : c.header.wtch2.watch2+1);
	                for (int ii = start_i; ii < c.size(); ii++) {
	                    if (assigns[var(c[ii])] == extbool_Undef) {
	                        SATAddWatcher(c, cr, va, ii); // watcher wird dort auch umgesetzt
	                        newWatcherFound = true;
	                        break;
	                    } else {
	                        if ( (assigns[var(c[ii])] == 1 && !sign(c[ii])) ||
	                            (assigns[var(c[ii])] == 0 &&  sign(c[ii]))	) {
	                            isSatisfied = true;
	                            SATswapOut(va,c);
	                            c.header.btch1.watch1 = c.header.wtch2.watch2 = -2;
	                            break;
	                        }
	                    }
	                }
	                if (conflict) continue;
	                if (isSatisfied) continue;
	                if (!newWatcherFound) {
	                    if (assigns[ var(c[c.header.btch1.watch1]) ] != extbool_Undef && assigns[ var(c[c.header.wtch2.watch2]) ] != extbool_Undef) {
	                        // kann passieren, wenn bei der Initialisierung nicht genau genug gearbeitet wurde
	                        if ((sign(c[c.header.btch1.watch1]) && assigns[var(c[c.header.btch1.watch1])] == 0) || (!sign(c[c.header.btch1.watch1]) && assigns[var(c[c.header.btch1.watch1])] == 1) ||
	                            (sign(c[c.header.wtch2.watch2]) && assigns[var(c[c.header.wtch2.watch2])] == 0) || (!sign(c[c.header.wtch2.watch2]) && assigns[var(c[c.header.wtch2.watch2])] == 1)  ) {
	                            // SAT Klausel ist erf�llt
	                            SATswapOut(va,c);
	                            c.header.btch1.watch1 = c.header.wtch2.watch2 = -2;
	                            continue;
	                        } else massert(0);
	                    }
	                    int l=0;
	                    if (assigns[ var(c[c.header.btch1.watch1]) ] == extbool_Undef) l = c.header.btch1.watch1;
	                    else if (assigns[ var(c[c.header.wtch2.watch2]) ] == extbool_Undef) l = c.header.wtch2.watch2;
	                    if (assigns[ var(c[l]) ] == extbool_Undef  && type[var(c[l])] == BINARY) {
	                        massert(VarsInConstraints[va][i].cr != CRef_Undef);
	                        constraintBumpActivity(constraintallocator[VarsInConstraints[va][i].cr]);
	                        if (propQlimiter[c[l].x] <= 0) {
	                            if (propQ.size() < PROPQ_LIMITER || feasPhase || !SUPPRESS_IMPLICATIONS) {
	                                PROPQ_PUSH(va,val,ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l));
	                                propQlimiter[c[l].x] = propQ.size();
	                                if (propQlimiter[c[l].x^1] > 0) {
	                                    conflict = true;
	                                    ValueConstraintPair tmp1(propQ[propQlimiter[c[l].x  ]-1].cr,propQ[propQlimiter[c[l].x  ]-1].v,propQ[propQlimiter[c[l].x  ]-1].pos);
	                                    ValueConstraintPair tmp2(propQ[propQlimiter[c[l].x^1]-1].cr,propQ[propQlimiter[c[l].x^1]-1].v,propQ[propQlimiter[c[l].x^1]-1].pos);
	                                    EmptyPropQ();
	                                    PROPQ_PUSH(va,val,tmp1);
	                                    propQlimiter[c[l].x  ] = 1;
	                                    propQlimiter[c[l].x^1] = 2;
	                                    PROPQ_PUSH(va,val,tmp2);
	                                    ix1 = propQlimiter[c[l].x  ]-1;
	                                    ix2 = propQlimiter[c[l].x^1]-1;
	                                }
	                            } // falsch! else propQ[propQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l);
	                        } else propQ[propQlimiter[c[l].x]-1] = ValueConstraintPair(VarsInConstraints[va][i].cr,c[l].x,l);
	                    }
	                }
	            }
	        }*/
	    }
    }


    void QBPSolver::CheckConstraint(ca_vec<CoeVar> & c, double rhs){
    	return;
    	double sum=0;
    	for (int i=0;i<c.size();i++){
			if(var(c[i])==4||
				var(c[i])==10||
				var(c[i])==16||
				var(c[i])==28||
				var(c[i])==22||
				var(c[i])==15||
				var(c[i])==25){
				if(sign(c[i])) sum-=c[i].coef;
				else sum += c[i].coef;
			}

		}   
		cerr <<"check: " << sum << "<="<< rhs << endl;
		if(sum<rhs){
			cerr <<"ERROR!"<<endl;
			cin.get();
		}
	}

	void QBPSolver::CreateAndUpdateShuffle(){
	if(UniversalConstraintsExist){	
		bool Print=false;
		if(PermutationOfAllStages.empty()){
			vector<int> NextStage;
			bool WriteThis=false;
			for (int i=0;i<nVars();i++){
				if(eas[i]==EXIST){ 
					if(WriteThis){
						if(!NextStage.empty()){
							PermutationOfAllStages.push_back(NextStage);
							NextStage.clear();
						}
					 	WriteThis=false;
					}
					continue;
				}
				if(eas[i]!=EXIST){
					WriteThis=true;
				 	NextStage.push_back(i);
				}
			}
		}
		for(int k=0;k<PermutationOfAllStages.size();k++){
			std::random_shuffle ( PermutationOfAllStages[k].begin(), PermutationOfAllStages[k].end() );
			if(Print){
				cerr << "Permutation of All Block " <<k<< endl;
				for (int i=0;i<PermutationOfAllStages[k].size();i++){
					cerr <<PermutationOfAllStages[k][i] << " ";
				}
				cerr << endl;
			}
		}
	}


	}
	void QBPSolver::PrintAssigns(){
		return;
		for (int i=0;i<nVars();i++){
			cerr << "x_" << i << "="; 
			if(assigns[i]!=extbool_Undef) cerr << (int)assigns[i];
			cerr<<endl;
		}	
	}

	void QBPSolver::setBranchingDecision(int &pick, int &left, int &right) {
		if(((yInterface*)yIF)->integers[pick].bitcnt==1 && abs(((yInterface*)yIF)->integers[pick].org_ub.asDouble()-((yInterface*)yIF)->integers[pick].org_lb.asDouble())<1e-6){
		   left=((yInterface*)yIF)->integers[pick].org_ub.asDouble();
		   right=left;
		   if(getShowWarning()) std::cerr << "WARNING: Fixed universal variable x_"<<pick <<" to " << left << " as given by its bounds"<<std::endl;
		}
		if(UniversalConstraintsExist){
			if(block[pick]-1<AllPropQStore.size()){
				while(!AllPropQStore[block[pick]-1].empty()){
					CoeVar XTO=mkCoeVar(AllPropQStore[block[pick]-1].back().first, 0, 1-AllPropQStore[block[pick]-1].back().second);
					assert(block[pick]==block[AllPropQStore[block[pick]-1].back().first]);

					if(AllpropQlimiter[XTO.x]==-2*nVars()-1 || AllpropQlimiter[XTO.x]<=-nVars()-1){
						if(PRINT_PROP) cerr << "Propagation of var " << AllPropQStore[block[pick]-1].back().first << " to " << ( 1-AllPropQStore[block[pick]-1].back().second) << " still in PropStore. "<< AllpropQlimiter[XTO.x]<< "  -> Deleted" << endl;
						AllpropQlimiter[XTO.x]=0;
						AllPropQStore[block[pick]-1].pop_back();
						continue;
					}
					else if (AllpropQlimiter[XTO.x]<0){
					   insertVarOrder(pick);
				           pick=AllPropQStore[block[pick]-1].back().first;
				           left=AllPropQStore[block[pick]-1].back().second;
				           if(PRINT_PROP) cerr << "Force Var " << pick << " to " << left << " due to " << -(AllpropQlimiter[XTO.x]+1) << endl;
				           right=left;//1-AllPropQStore[block[pick]-1].back().second;
				           AllPropQStore[block[pick]-1].pop_back();
   				           assert(assigns[pick]==extbool_Undef);
				           assert(AllpropQlimiter[XTO.x]<0);
				           AllpropQlimiter[XTO.x]=AllpropQlimiter[XTO.x]-nVars();
				           break;
				        }
				        else if (AllpropQlimiter[XTO.x]==0){
				    	    AllPropQStore[block[pick]-1].pop_back();
				        }

				    else if(getShowWarning()) cerr << "WARNING: SetBranchingDecision strange value: " <<  AllpropQlimiter[XTO.x] << " " << nVars() <<endl;
			    }
	      	}
		}
	}

    void QBPSolver::moveDown(int d, int decvar, int decpol, int pick,string where) {
    	return;
//    	cerr << "DOWN: x_" <<decvar <<" = " << decpol <<  " pick= " <<pick  << " from " << where<< endl;
	for(int i=0;i<trail.size();i++) cerr <<"x_"<<trail[i]<<"="<<(int)assigns[trail[i]] << " L:"<<vardata[trail[i]].level << " F:"<<isFixed(trail[i])<<" FL:"<<fixdata[trail[i]].level << endl;   
	//cerr <<"Stack Status: " <<search_stack.stack[search_stack.stack_pt].status << endl;
    	int8_t *val;
	stack_container &STACKz = search_stack.stack[d-2];
    	val = STACKz.val;//&stack_val[(d-1)<<1];
    	int8_t &val_ix = STACKz.val_ix;//stack_val_ix[(d-1)];

    	//cerr << "Val_ix: "<<(int)val_ix << "  val[0]: " << (int)val[0] << " val[1]: " <<(int)val[1] << endl;

    	//if (decvar==15) PrintAssigns();
    	return;
    	if (!UniversalConstraintsExist) return;

        //((yInterface*)yIF)->moveDown(d, decvar, decpol, pick);
    	/*std::cerr <<"MOVE DOWN " << d << " " << decvar << " " << decpol << " " << pick<< std::endl;
    	int8_t *val;
    	val = &stack_val[(d-1)<<1];
    	int8_t &val_ix = stack_val_ix[(d-1)];

    	cerr << (int)val_ix << " " << (int)val[0] << " " <<(int)val[1] << endl;
*/
    	assert(d>0);

    }

    void QBPSolver::moveUp(coef_t & v, coef_t b, int status) {
    	 return;
stack_container &STACK = search_stack.stack[search_stack.stack_pt-1];
    	//cerr << "DL:"<<getLastDecisionLevel()<<endl;
    	//cerr <<"Stack_Pick.size()"<<stack_pick.size() << endl;
    	int DecLeva = getLastDecisionLevel();
    	//cerr << "Stack pick " <<STACK.pick << endl;
    	int DecVara = search_stack.stack[DecLeva-1].pick;
    	//cerr << "Dec Var : " << DecVara << endl;
    	int DecPola=assigns[DecVara];

    	cerr << "Up: x_" <<DecVara <<" = " << DecPola <<  " v= " << v << " b=" << b << endl;

    	return;
    	bool MPrint=true;
    	if (!UniversalConstraintsExist) return;

    	//  Fix the already set variables in the relaxation
    	for (int hh = 0; hh < nVars();hh++) {
			if (type[hh] != BINARY) continue;
			if ((getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef)|| vardata[hh].level==getLastDecisionLevel()) {
				QlpStSolve->setVariableLB(hh,0,type.getData());
				QlpStSolve->setVariableUB(hh,1,type.getData());
			} else if (assigns[hh] != extbool_Undef) {
				QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
			} else {
				QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
			}

			updateStageSolver(maxLPStage,hh,hh);
		}


    	int DecLev = getLastDecisionLevel();
    	int DecVar = search_stack.stack[DecLev-1].pick;//
    	int DecPol=assigns[DecVar];

    	if(nVars()==trail.size()){    	//If at a leaf:
    		assert(order_heap.empty());
    		assert(v!=n_infinity);
    		if(v!=n_infinity) ExistLegalUntil=DecLev;
    		for (int v=0;v<VarsPresentInAllConstraints.size();v++){
    			int i=VarsPresentInAllConstraints[v];
    				if(VarsInAllConstraints[i].size()==0) continue;
				if (assigns[i]==extbool_Undef) {
					AllSolver->setVarLB(i,0);
					AllSolver->setVarUB(i,1);
				}
				else if (assigns[i]==1){ AllSolver->setVarLB(i,1);AllSolver->setVarUB(i,1);}
				else if (assigns[i]==0){ AllSolver->setVarLB(i,0);AllSolver->setVarUB(i,0);}
			}
			AllSolver->solve();
			if(AllSolver->getSolutionStatus()==extSol::QpExternSolver::QpExtSolSolutionStatus::OPTIMAL)
				AllLegalUntil=DecLev;
			else{
				cerr<<"All Infeasible at leaf "<< endl;
				v=AllInfeasible;
			}
		}

    	if (v!=n_infinity && v!=dont_know&& AllLegalUntil<DecLev && v!=AllInfeasible){// normal Exist-Strategy exists for this node
    		assert(ExistLegalUntil==DecLev);
    		//Check whether universal system is broken
    		for (int v=0;v<VarsPresentInAllConstraints.size();v++){
		    	int i=VarsPresentInAllConstraints[v];
		    		if(VarsInAllConstraints[i].size()==0) continue;
				if (assigns[i]==extbool_Undef) {
					AllSolver->setVarLB(i,0);
					AllSolver->setVarUB(i,1);
				}
				else if (assigns[i]==1){ AllSolver->setVarLB(i,1);AllSolver->setVarUB(i,1);}
				else if (assigns[i]==0){ AllSolver->setVarLB(i,0);AllSolver->setVarUB(i,0);}
			}
			AllSolver->solve();
			if(AllSolver->getSolutionStatus()==extSol::QpExternSolver::QpExtSolSolutionStatus::OPTIMAL){
				cerr << "Universal System OK" << endl;
				AllLegalUntil=DecLev;
			}
			else{
				cerr<<"All Infeasibility Found "<< endl;
				v=AllInfeasible;
			}
    	}



    	cerr << "All-System fine until " << AllLegalUntil << endl;
    	if((AllLegalUntil==DecLev && eas[DecVar]==UNIV )||(ExistLegalUntil==DecLev && eas[DecVar]==EXIST))
    		cerr << "Legality Ensured at " << DecLev << endl;

    	//getLastDecisionVariable();
    	//((yInterface*)yIF)->moveUp(v, b, status);
    	if (MPrint){
    		std::cerr <<"MOVE UP " << v << " " << b << " " << status << ": " << DecVar << " set to " << (int)assigns[DecVar] << " at Level "<<DecLev << std::endl;
		stack_container &STACKz = search_stack.stack[decisionLevel()-2];
        	int8_t *val;
		val = STACKz.val;//&stack_val[(decisionLevel()-1)<<1];
		int8_t &val_ix = STACKz.val_ix;//stack_val_ix[(decisionLevel()-1)];

        	cerr << (int)val_ix << " " << (int)val[0] << " " <<(int)val[1] << endl;
    	}
    	assert (getBranchingVariable(DecLev)==DecVar);

    	if (DecLev >0){
			assert(DecLev+1==search_stack.stack.size());
			assert(assigns[DecVar]!=extbool_Undef);

			if (DecVar != -1 &&  v==AllInfeasible && ExistLegalUntil<DecLev){
				if (MPrint) std::cerr << "Existential Decision Variable Resulting in All_Infeas -> Check legality of setting existential variable " << DecVar << " in Level " << DecLev<< std::endl;
				std::vector<data::QpNum> IntegerSolution;
				double IntegerScore;
				int CompTime=-1;
				algorithm::Algorithm::SolutionStatus LPstatus=algorithm::Algorithm::UNKNOWN;

				for (int hh = 0; hh < nVars();hh++) {
					if (type[hh] != BINARY) continue;
					if (getFixed(hh) == extbool_Undef && assigns[hh] == extbool_Undef) {
						QlpStSolve->setVariableLB(hh,0,type.getData());
						QlpStSolve->setVariableUB(hh,1,type.getData());
					} else if (assigns[hh] != extbool_Undef) {
						QlpStSolve->setVariableFixation(hh,(double)assigns[hh],type.getData());
					} else {
						QlpStSolve->setVariableFixation(hh,(double)getFixed(hh),type.getData());
					}

					updateStageSolver(maxLPStage,hh,hh);
				}
				QLPSTSOLVE_SOLVESTAGE(fmax((double)constraintallocator[constraints[0]].header.rhs,0),maxLPStage, LPstatus, lb, ub, IntegerSolution,algorithm::Algorithm::WORST_CASE,-1,/*-1*/feasPhase?-1:/*3*/3-4 /*simplex iterationen*/);
				if (LPstatus == algorithm::Algorithm::INFEASIBLE) cerr << "Found Infeasible " << endl;
				FeasibilityPump FP(IntegerSolution,QlpStSolve,((yInterface*)yIF)->qlp, maxLPStage, type.getData(), nVars());
				if(FP.FindIntegerSolution(IntegerSolution, IntegerScore, type.getData(),1,/*AR*/0,1,CompTime,this)){
					ExistLegalUntil=DecLev;
					if (MPrint){
						std::cerr << "Existential System Still Satisfiable" << std::endl;
						for (int i=0;i<nVars();i++)
							cerr << i << " " << IntegerSolution[i] << endl;
					}
				}
				else{

					//TODO: Find out for sure whether existential system is satisfiable

				}

			}
			else if(DecVar != -1  && v==n_infinity &&  AllLegalUntil<DecLev){
				for (int i=0;i<nVars();i++){
					if(VarsInAllConstraints[i].size()==0) continue;
					if (assigns[i]==extbool_Undef) {
						AllSolver->setVarLB(i,0);
						AllSolver->setVarUB(i,1);
					}
					else if (assigns[i]==1){ AllSolver->setVarLB(i,1);AllSolver->setVarUB(i,1);}
					else if (assigns[i]==0){ AllSolver->setVarLB(i,0);AllSolver->setVarUB(i,0);}

				}

					if (MPrint){
				std::vector<data::QpNum> LBvec;
				std::vector<data::QpNum> UBvec;
				AllSolver->getLB(LBvec);
				AllSolver->getUB(UBvec);
				cerr << "Print Variable Bounds in All-Solver:" << endl;
				for (int i=0;i<AllSolver->getVariableCount();i++)
					cerr << i << " " << LBvec[i].asDouble() << " " <<  UBvec[i].asDouble() << endl;
				}
				AllSolver->solve();
				if (MPrint) cerr << "Solved " << AllSolver->getSolutionStatus() << endl;
				if(AllSolver->getSolutionStatus()==extSol::QpExternSolver::QpExtSolSolutionStatus::INFEASIBLE){
					if (MPrint) cerr << "Universal System is Infeasible " << endl;
					v=AllInfeasible;
				}
				else{
					AllLegalUntil=DecLev;
					if (MPrint){
						std::vector<data::QpNum> AllSol;
						AllSolver->getValues(AllSol);
						for (int i=0;i<nVars();i++)
						cerr << i << " " << AllSol[i] << " ("<<(int)assigns[i] <<")" << endl;
						cerr <<"All Sol end" << endl;
					}
				}
			}
    	}
    	if ( ExistLegalUntil>=DecLev) ExistLegalUntil=DecLev-1;
    	if  (AllLegalUntil>=DecLev) AllLegalUntil=DecLev-1;
    }


bool QBPSolver::AllConstraintsStillSatisfiable(int va){
    //Checks satisfiability of Universal Contraints in which variable "va" is present
    //Problem must be simply restricted; normally knowing that a constraint is still satisfiable does not let you say anything about the system.
    for (int i=0; i < VarsInAllConstraints[va].size();i++) {
        Constraint &c = ALLconstraintallocator[VarsInAllConstraints[va][i].cr];
        if(c.header.btch1.best<c.header.rhs) return false;
    }
    return true;
}

bool QBPSolver::AllConstraintsStillSatisfiable(int va,double val){
    //Checks whether setting variable "va" to "val" immediately results in a violation of the constraints this variable is present in
    //Problem must be simply restricted; normally knowing that a constraint is still satisfiable does not let you say anything about the system.
    for (int i=0; i < VarsInAllConstraints[va].size();i++) {
        Constraint &c = ALLconstraintallocator[VarsInAllConstraints[va][i].cr];
        int pos = VarsInAllConstraints[va][i].pos;
        int s = sign(c[pos]);
        //if (c.header.isSat && c.header.btch1.watch1 == -2) continue;
        massert(!c.header.deleted);
        c.header.isSat=false;
        if (c.saveInfeas(pos, val, /*getFixed(va),*/ va, assigns, eas)) {
            /*cout <<"ASSIGN GIBT " << VarsInAllConstraints[va][i].cr << endl;
            cerr << c.header.btch1.best << " is best" << endl;
            c.print(c,assigns,false);
            cout << VarsInAllConstraints[va][i].cr << " var=" << va << " und value=" << val << endl;
            for (int j=0; j < c.size();j++) {
                cout << "-- " << vardata[var(c[j])].reason << "--" << vardata[var(c[j])].level << endl;
            }*/
            return false;
        }
    }
    return true;
}

void QBPSolver::SetBoundsInExternalIPSolver(extSol::QpExternSolver* Solver, int i){
    if(decisionLevel()==1){
        if(block[i]==1 || trail.size()==nVars()){
            if (assigns[i]==1){ Solver->setVarLB(i,1);Solver->setVarUB(i,1);}
            else if (assigns[i]==0){ Solver->setVarLB(i,0);Solver->setVarUB(i,0);}
            else if (isFixed(i)) { Solver->setVarLB(i,getFixed(i));Solver->setVarUB(i,getFixed(i));}
            else  {Solver->setVarLB(i,0); Solver->setVarUB(i,1);}
        }
        else {Solver->setVarLB(i,0); Solver->setVarUB(i,1);}
    }
    else{
        if ( ( assigns[i]==extbool_Undef && !isFixed(i)) || (trail.size()!=nVars() && block[i]> block[trail[getDecisionIndexInTrail(decisionLevel())]] ) ){
	    //Do not fix variables of later blocks in the external IP solver; they might disturb the feasibility. Unless the trail is already full
            Solver->setVarLB(i,0);
            Solver->setVarUB(i,1);
        }
        else if (assigns[i]==1){ Solver->setVarLB(i,1);Solver->setVarUB(i,1);}
        else if (assigns[i]==0){ Solver->setVarLB(i,0);Solver->setVarUB(i,0);}
        else if (isFixed(i)) { Solver->setVarLB(i,getFixed(i));Solver->setVarUB(i,getFixed(i));}
        else{
	    if(getShowWarning()) cerr << "Error: Trying to set bounds of simultaneously assigned and unassigned variable. Freed it." << endl;
            // Propagated variables of later blocks should not be considered here!
            Solver->setVarLB(i,0);
            Solver->setVarUB(i,1);
        }
    }
}
bool QBPSolver::ExistIPStillFeasible(){
    //Checks whether there still exists a solution for the existential constraint system
    for (int i=0;i<nVars();i++){
	SetBoundsInExternalIPSolver(ExistSolver,i);
        /*if (assigns[i]==extbool_Undef){
            ExistSolver->setVarLB(i,0);
            ExistSolver->setVarUB(i,1);
        }
        else if (assigns[i]==1){ ExistSolver->setVarLB(i,1);ExistSolver->setVarUB(i,1);}
        else if (assigns[i]==0){ ExistSolver->setVarLB(i,0);ExistSolver->setVarUB(i,0);}
	*/
    }
    ExistSolver->solve();
    if(ExistSolver->getSolutionStatus()==extSol::QpExternSolver::QpExtSolSolutionStatus::OPTIMAL){
        if(!NO_DEBUG){ 
            cerr<<"Exist-IP still feasible"<<endl;
            std::vector<data::QpNum> Sol;
            ExistSolver->getValues(Sol);
            for (int i=0;i<nVars();i++)
                cerr << i << " " << Sol[i].asDouble() << " ("<<(int)assigns[i] <<")" << endl;
        }
        return true;
    }
    else{
        if(!NO_DEBUG) cerr<<"Exist-IP infeasible"<<endl;
        return false;
    }
}

bool QBPSolver::AllIPStillFeasible(){

    if(!NO_DEBUG)  cerr << "Check after DV: x_"<<trail[getDecisionIndexInTrail(decisionLevel())] << " " <<trail.size()<<endl;
    for (int v=0;v<VarsPresentInAllConstraints.size();v++){
    	int i=VarsPresentInAllConstraints[v];
        if(VarsInAllConstraints[i].size()==0) continue; //Should be deprecated due to new datastructure VarsPresentInAllConstraints
	SetBoundsInExternalIPSolver(AllSolver,i);
	/*if(decisionLevel()==1){
	    if (assigns[i]==1){ AllSolver->setVarLB(i,1);AllSolver->setVarUB(i,1);}
            else if (assigns[i]==0){ AllSolver->setVarLB(i,0);AllSolver->setVarUB(i,0);}
	    else if (isFixed(i)) { AllSolver->setVarLB(i,getFixed(i));AllSolver->setVarUB(i,getFixed(i));}
	    else  {AllSolver->setVarLB(i,0); AllSolver->setVarUB(i,1);}
	}
	else{
            if (assigns[i]==extbool_Undef || 
            	    (trail.size()!=nVars() &&
        	 
        	    ((decisionLevel()==1&& block[i]>1) || (block[i]> block[trail[getDecisionIndexInTrail(decisionLevel())]])))){
            // Propagated variables of later blocks should not be considered here!
                AllSolver->setVarLB(i,0);
                AllSolver->setVarUB(i,1);
            }
            else if (assigns[i]==1){ AllSolver->setVarLB(i,1);AllSolver->setVarUB(i,1);}
            else if (assigns[i]==0){ AllSolver->setVarLB(i,0);AllSolver->setVarUB(i,0);}
	}*/
    }
    AllSolver->solve();
    if(AllSolver->getSolutionStatus()==extSol::QpExternSolver::QpExtSolSolutionStatus::OPTIMAL){
        if(!NO_DEBUG){
            cerr<<"All-IP still feasible"<<endl;
            std::vector<data::QpNum> Sol;
            AllSolver->getValues(Sol);
            for (int i=0;i<nVars();i++)
                cerr << i << " " << Sol[i].asDouble() << " ("<<(int)assigns[i] <<")" << endl;
        }
        return true;
    }
    else{
        if(!NO_DEBUG) cerr<<"All-IP infeasible"<<endl;
        return false;
    }
}

bool QBPSolver::AllIPStillFeasible(std::vector<data::QpNum>& sol){
    assert (sol.size()==nVars());
    for (int v=0;v<VarsPresentInAllConstraints.size();v++){
        int i = VarsPresentInAllConstraints[v];
        AllSolver->setVarLB(i,sol[i]);
        AllSolver->setVarUB(i,sol[i]);
	}
    AllSolver->solve();
    if(AllSolver->getSolutionStatus()==extSol::QpExternSolver::QpExtSolSolutionStatus::OPTIMAL){
        if(!NO_DEBUG){
            cerr<<"All-IP still feasible"<<endl;
            std::vector<data::QpNum> Sol;
            AllSolver->getValues(Sol);
            for (int i=0;i<nVars();i++)
                cerr << i << " " << Sol[i].asDouble() << " ("<<(int)assigns[i] <<")" << endl;
        }
        return true;
    }
    else{
        if(!NO_DEBUG) cerr<<"All-IP infeasible"<<endl;
        return false;
    }
}

bool QBPSolver::AllSystemSatisfied(){
    // Currently unused
    for (int j=0;j<ALLconstraints.size();j++){
        Constraint &c =ALLconstraintallocator[ALLconstraints[j]];
        if(c.saveInfeas()) return false;
    }
    return true;
}

bool QBPSolver::CheckAllFeasibility(int va, int val){
    if(!UniversalConstraintsExist) return true;
    if(!NO_DEBUG){
        cerr << "check x_"<<va<<"="<<val<< " " << VarsInAllConstraints[va].size()<<endl;
        for (int i=0;i<trail.size();i++)
            cerr << "x_"<<trail[i]<<"="<<(int)assigns[trail[i]]<<endl;
    }
    if(getIsSimplyRestricted()){
        if(!AllConstraintsStillSatisfiable(va,val)) return false;
        if(!isUniversalPolytope()){
   	        for (int v=trail.size()-1;v>=0 ;v--){
   	            if(eas[trail[v]]==UNIV) break;
   	            //Check existential variables of previous block; Maybe they alredy destroyed the universal constraint system?
                else if(VarsInAllConstraints[trail[v]].size()>0 && !AllConstraintsStillSatisfiable(trail[v])) return false;
   	        }
   	    }
        return true;
    }
    else{
        for (int v=0;v < VarsPresentInAllConstraints.size();v++){
            int i = VarsPresentInAllConstraints[v];
            if(VarsInAllConstraints[i].size()==0) continue;
            if (assigns[i]==extbool_Undef || block[i]>block[va]){
             	AllSolver->setVarLB(i,0);
             	AllSolver->setVarUB(i,1);
                if(!NO_DEBUG) cerr << "Set bounds of x_" << i<< " to {0,1}" << endl;
            }
            else if (assigns[i]==1){/*bdLB[i]=1;bdUB[i]=1;*/ AllSolver->setVarLB(i,1);AllSolver->setVarUB(i,1); if(!NO_DEBUG) cerr << "Set bounds of x_" << i<< " to 1" << endl;}
            else if (assigns[i]==0){/*bdLB[i]=0;bdUB[i]=0*/; AllSolver->setVarLB(i,0);AllSolver->setVarUB(i,0);if(!NO_DEBUG) cerr << "Set bounds of x_" << i<< " to 0" << endl;}
            if(i==va) {/*bdLB[i]=val;bdUB[i]=val;*/ AllSolver->setVarLB(i,val);AllSolver->setVarUB(i,val);if(!NO_DEBUG)cerr << "Set bounds of x_" << i<< " to "<<val << endl;}
		}
        AllSolver->solve();
        if(AllSolver->getSolutionStatus()==extSol::QpExternSolver::QpExtSolSolutionStatus::OPTIMAL){
            if(!NO_DEBUG) cerr <<"Check: All Feasible"<< endl;
	    	return true;
	}
        else{
            if(!NO_DEBUG) cerr<<"Check: All Infeasible"<< endl;
	    	return false;
        }
    }
}

	/*bool QBPSolver::CheckAllFeasibility(int va, int val){
    	// std::cerr << "x0 " << isFixed(0) <<std::endl;
    	//std::cerr <<"Start Check! " << va << " " <<val <<std::endl;
    	//for (int i=0;i<trail.size();i++) std::cerr << "x"<<trail[i]<<"="<<(int)assigns[trail[i]]<<std::endl;

    	//QUICK CHECK, WHEATER SOME CONSTRAINTS CANNOT BE SATISFIED EVEN IN FOR BEST CASE ASSIGNMENT
    	bool Print =false;
    	//int iFull=-1;
    	for (int i=0; i < VarsInAllConstraints[va].size();i++) {
    	    int pos = VarsInAllConstraints[va][i].pos;
    	    CRef cr = VarsInAllConstraints[va][i].cr;
    	    Constraint &c =ALLconstraintallocator[cr];
    	    
    	    	

    	    coef_t rhs = c.header.rhs;
    	    coef_t bc=0.0;

    	    for(int h=0;h<c.size();h++){
    	    	if(var(c[h])==va) bc+=(sign(c[h])?-c[h].coef*val:c[h].coef*val);
    	    	else if(assigns[var(c[h])]==extbool_Undef ){
    	    		if(!sign(c[h])) bc+=c[h].coef;
    	    	}
    	    	else{
    	    		bc+=(sign(c[h])?-c[h].coef*assigns[var(c[h])]:c[h].coef*assigns[var(c[h])]);
    	    	}
    	    }
    	    if(bc<rhs){
    	    	if(Print){
    	    	for(int h=0;h<c.size();h++){
    	    		if(sign(c[h])) cerr << " -"<<c[h].coef <<"x"<<var(c[h])<< "(" << (int)assigns[var(c[h])]<<")";
    	    		else  cerr << " +"<<c[h].coef <<"x"<<var(c[h]) << "(" << (int)assigns[var(c[h])]<<")";
    	    	}
    	    	cerr <<">= " << c.header.rhs<<endl;
    	    }
    	    	if (Print)
    	    		cerr <<"Illegal Universal Move setting x"<<va <<"="<<val << endl;
    	    	return false;
    	    } 

    	    //cerr << "bc: " << bc << "and rhs " << rhs << endl;
    	}
    	if (Print) cerr <<"LEGAL Move setting x"<<va <<"="<<val << endl;    	
    	return true;
    }*/

 void QBPSolver::EmptyAllPropQ(bool Single, bool PrintWarning,bool OnlyPop){
	//cerr <<"Clear PropQ for All " << AllPropQ.size()<<endl;
	if (Single){
		if(!OnlyPop)
				AllpropQlimiter[(AllPropQ[AllPropQ.size()-1].first<<1)|(1-AllPropQ[AllPropQ.size()-1].second)] = 0;
		//cerr << "Delete " << AllPropQ.back().first<<"/"<<AllPropQ.back().second<<endl;
		//cerr << "Set Limiter " << ((AllPropQ[AllPropQ.size()-1].first<<1)|(1-AllPropQ[AllPropQ.size()-1].second)) << " to 0" << endl; 
		    AllPropQ.pop_back ();
	}
	else{
		while (AllPropQ.size() > 0) {
			//cerr << "delete entry "<< AllPropQ[AllPropQ.size()-1].first << "-" <<(AllPropQ[AllPropQ.size()-1].first<<1) << " and " <<AllPropQ[AllPropQ.size()-1].second <<"-"<<((AllPropQ[AllPropQ.size()-1].first<<1)|(1-AllPropQ[AllPropQ.size()-1].second))<<endl;
			if (PrintWarning)
                std::cerr << "Warning: AllPropQ not empty." << std::endl;
			if(!OnlyPop)
				AllpropQlimiter[(AllPropQ[AllPropQ.size()-1].first<<1)|(1-AllPropQ[AllPropQ.size()-1].second)] = 0;
			//cerr << "Delete " << AllPropQ.back().first<<"/"<<AllPropQ.back().second<<endl;
			//cerr << "Set Limiter " << ((AllPropQ[AllPropQ.size()-1].first<<1)|(1-AllPropQ[AllPropQ.size()-1].second)) << " to 0" << endl; 
			AllPropQ.pop_back ();
		}
	}
}

int QBPSolver::PropagateUniversals(){
	//return 0, if existential system fails
	//reutrn 1, if both systems ok
	//return 2, if universal system fails
        float alpha=(float)global_score;

	if(!UniversalPolytope){
//		cerr << " No Universal Propagation" << endl;
//	 return true;
	}
  	//EmptyAllPropQ();
  	//return;
	if (decisionLevel()<=1){
		EmptyAllPropQ();
	 	return 1;
	}
	while (!AllPropQ.empty()){
		int var=AllPropQ[AllPropQ.size()-1].first;
		int val=AllPropQ[AllPropQ.size()-1].second;
		EmptyAllPropQ(true);
		if (assigns[var]!= extbool_Undef){
                    if(assigns[var]!=val){
			//for (int i=0;i<trail.size();i++) cerr << "x_"<<trail[i]<<"=" <<(int) assigns[trail[i]] << endl;
                	if(getShowError()) cerr <<"Error: Universal Variable already assigned THE WRONG WAY: x_" << var << " cannot be propagated to be " << val << endl;
			assert(0);
		    }
		    else{ 
			//cerr <<"Universal Variable already assigned" << endl;
                        continue;
                    }
                }
		if((!getIsSimplyRestricted() /*||UniversalMultiBlockConstraints*/)&&block[search_stack.stack[getLastDecisionLevel()-1].pick]!=block[var]){ 
			while(AllPropQStore.size()<max(maxBlock,block[var])){
				std::vector<std::pair<int,bool>> NewVec;
				AllPropQStore.push_back(NewVec);
			}
			if(AllpropQlimiter[(var<<1)|(1-val)]!=0){
				//cerr <<"Does this happen? A universal variable that is already marked as propagted!?" << endl;
				//Should not happen, as only variables without that property should be pushed into the AllPropQ
				continue;
			}
			AllPropQStore[block[var]-1].push_back( make_pair(var, val));
			if(PRINT_PROP){
				cerr << "Omitted Propagation of Later Variable"<<endl;
				cerr << block[search_stack.stack[getLastDecisionLevel()-1].pick] << " " << block[var] << endl;
				cerr << getLastDecisionLevel() << " " << var << " to " << val << endl;
			}
			//EmptyAllPropQ();
	 		//return true;
			if(AllVarsProppedDueTo.size()<nVars())     AllVarsProppedDueTo= vector<vector<int>>(nVars(), vector<int>(0));
			AllVarsProppedDueTo[search_stack.stack[getLastDecisionLevel()-1].pick].push_back((var<<1)|(1-val));
			if(PRINT_PROP) cerr << "Remembered, that x_" << var << "=" << val << " is propagated due to " << search_stack.stack[getLastDecisionLevel()-1].pick << endl;
			AllpropQlimiter[(var<<1)|(1-val)]=-search_stack.stack[getLastDecisionLevel()-1].pick-1;
			 // Will be forced, due to Variable as indicated negative in AllpropQlimiter
 			continue;
		}

		//if(assigns[var] != extbool_Undef && assigns[var] !=val){
		if(assigns[var]==1-val || getFixed(var) == 1-val){
			if(getShowError()) cerr << "Error: Propagated universal variable x_"<<var<<" cannot be set to " << val << endl;
			assert(0);
		} 
      	int ix1, ix2;
      	bool conflict=false;
      	int oob = assign(alpha, var,val,trail.size(),CRef_Undef, conflict, ix1, ix2, false);
      	//int oob = assign(var, val,trail.size(),CRef_Undef,true);
      	//increaseDecisionLevel();
      	if(PRINT_PROP)cerr << "Propagated Universal Variable x_" << var<<"="<<val << endl;
      if(oob==ASSIGN_UNIV_FAIL){
	    EmptyAllPropQ();
	   //cerr << "Propagated Universal Variable x_" << var<<"="<<val << endl;
	   // for (int i=0;i<nVars();i++) if (eas[i]==UNIV) cerr << "x_"<<i<<"="<<(int)assigns[i]<<endl;
	   // for (int i = 0; i< trail.size();i++) {
	   // cerr << i << ":" << trail[i] << "(" << vardata[trail[i]].level << ") ";
	   //}
	    if(getShowWarning()) cerr <<"WARNING: Universal propagation led to contradiction in universal constraint system!?" << endl;
	    return 2;
      }
          	    		NumAllProps++;
      if(oob!=ASSIGN_OK){
	//if(getShowWarning()) cerr <<"WARNING: Universal Propagation led to contradiction" << endl;
      	EmptyAllPropQ();
	//HIER GEÄNDERT
	EmptyPropQ();
      	return 0;
      }
      else fixdata[var].reason=0; //This way a universal variable can be identified as propagated

	}
	return 1;
}
bool QBPSolver::CheckValForScenario(int va, int & val, Scenario_t & sc){
	//cerr <<"CheckValForScenAll1"<<endl;
	if(!UniversalConstraintsExist) return true;
		assert(UniversalPolytope);

	for (int i=0; i < VarsInAllConstraints[va].size();i++) {
    	    CRef cr = VarsInAllConstraints[va][i].cr;
    	    Constraint &c =ALLconstraintallocator[cr];

    	    coef_t rhs = c.header.rhs;
    	    coef_t bc=0.0;

    	    for(int h=0;h<c.size();h++){
    	    	if(var(c[h])==va) bc+=(sign(c[h])?-c[h].coef*val:c[h].coef*val);
    	    	else{
    	    		bool Found=false;
    	    		for (int spot=0;spot < sc.scen_var.size();spot++){
						if(sc.scen_var[spot]==var(c[h])){
							bc+=(sign(c[h])?-c[h].coef*sc.scen_val[spot]:c[h].coef*sc.scen_val[spot]);
							Found=true;
							break;
						}
					}
					if(!Found){
	    	    		if(!sign(c[h])) 
	    	    			bc+=c[h].coef;
	    	    	}
    	    	} 

    	    }
    	    if (bc>=rhs) continue;
    	    else return false; 	   
    	}
    	//cerr <<"LEGAL Move setting x"<<va <<"="<<val << endl;

    	return true;
	
}

bool QBPSolver::CheckValForScenario(int va, int & val, std::vector<std::pair<int,bool>> & sc){
	//cerr <<"CheckValForScenAll"<<endl;
	if(!UniversalConstraintsExist) return true;
		assert(UniversalPolytope);


	for (int i=0; i < VarsInAllConstraints[va].size();i++) {
    	    CRef cr = VarsInAllConstraints[va][i].cr;
    	    Constraint &c =ALLconstraintallocator[cr];

    	    coef_t rhs = c.header.rhs;
    	    coef_t bc=0.0;

    	    for(int h=0;h<c.size();h++){
    	    	if(var(c[h])==va) bc+=(sign(c[h])?-c[h].coef*val:c[h].coef*val);
    	    	else if(assigns[var(c[h])]!=extbool_Undef)
    	    		bc+=(sign(c[h])?-c[h].coef*assigns[var(c[h])]:c[h].coef*assigns[var(c[h])]);
    	    	else{
    	    		bool Found=false;
    	    		for (int spot=0;spot < sc.size();spot++){
						if(sc[spot].first==var(c[h])){
							bc+=(sign(c[h])?-c[h].coef*sc[spot].second:c[h].coef*sc[spot].second);
							Found=true;
							break;
						}
					}
					if(!Found){
	    	    		if(!sign(c[h])) 
	    	    			bc+=c[h].coef;
	    	    	}
    	    	} 

    	    }
    	    if (bc>=rhs) continue;
    	    else return false; 	   
    	}
    	//cerr <<"LEGAL Move setting x"<<va <<"="<<val << endl;
    	//cerr << "ok"<<endl;
    	return true;
	
}
    bool QBPSolver::CheckScenario(Scenario_t & sc, int dl){
    		assert(UniversalPolytope);

	//cerr <<"CheckScenario"<<endl;

		for (int j=0;j<ALLconstraints.size();j++){
			double LhsV=0;
			Constraint &c =ALLconstraintallocator[ALLconstraints[j]];
			coef_t rhs = c.header.rhs;
			//cerr <<"NewConstraint"<<endl;
			for(int h=0;h<c.size();h++){
				for (int spot=0;spot < sc.scen_var.size();spot++){
					if(sc.scen_var[spot]==var(c[h])){
						//cerr <<"x_"<<var(c[h]) << "=" << sc.scen_val[spot] << " Coef="<<(sign(c[h])?-c[h].coef:c[h].coef) << endl;
						LhsV+=(sign(c[h])?-c[h].coef*sc.scen_val[spot]:c[h].coef*sc.scen_val[spot]);
						break;
					}
				}
    	    }
    	    if (LhsV< c.header.rhs) return false;
		}
		return true;

		/*std::vector<data::QpNum> LBvec;AllSolver->getLB(LBvec);
		std::vector<data::QpNum> UBvec;AllSolver->getUB(UBvec);
		cerr << "Print Variable Bounds in All-Solver:" << endl;
		 for (int i=0;i<AllSolver->getVariableCount();i++)
			cerr << i << " " << LBvec[i].asDouble() << " " <<  UBvec[i].asDouble() << endl;
*/



    	//cerr <<"LEGAL Move setting x"<<va <<"="<<val << endl;
		// Mit All-IP-Solver
		/*AllSolver->solve();
		if(AllSolver->getSolutionStatus()==extSol::QpExternSolver::QpExtSolSolutionStatus::OPTIMAL){
			cerr << "Scenario fine" << endl;
			return true;

		}
		else{
			cerr << "Scenario bad" << endl;
			return false;

		}*/
        }


int QBPSolver::GenerateReformulationCut( extSol::QpExternSolver& externSolver, vector< vector< data::IndexedElement > > &listOfCutsLhs,
		       vector< data::QpNum > &listOfCutsRhs, vector<int> &listOfCutsVars,
					  int treedepth, int currentBlock, bool &global_valid, std::vector<unsigned int> &candidates, int cuttype, int*types, int8_t* assigns, unsigned int initime, int* solu, int *fixs, int *blcks, int orgN){
		//return 0;
		//return listOfCutsLhs.size();
 		const double eps = 1e-7;
		bool PrintSum=true;
		bool Print=false;
		bool Strengthening=true;
		// Constraint Orientation of Constraints in constraintallocator: \sum ax >= b

		if (externSolver.getSolutionStatus() != extSol::QpExternSolver::OPTIMAL){
					if(getShowWarning()) cerr << "Warning: invalid solution status -> no reformulation cuts" << endl;
					//cin.get();
					return listOfCutsRhs.size();
				}
		const unsigned int n = externSolver.getVariableCount();
		 vector<data::QpNum> objVals(n);
		externSolver.getValues(objVals);

		//cerr << " x_493= " << 19-objVals[493].asDouble()*158 << "x_513=" << 54-objVals[513].asDouble()*139 <<  " x_533= " << 42-objVals[533].asDouble()*85 <<  " x_553= " << 43-objVals[553].asDouble()*43 << endl;
		//cerr << objVals[253].asDouble()+objVals[273].asDouble()+objVals[293].asDouble()+objVals[313].asDouble() << " " << objVals[493].asDouble()*158 +objVals[513].asDouble()*139 +objVals[533].asDouble()*85 +objVals[553].asDouble()*43 << endl;
		//cin.get();
		if(Print) cerr <<"_RC_ "<< n << " " << nVars() << endl;
		int numCuts=n;
		/*static*/ vector< pair< vector< pair<unsigned int, double> >, double > > cuts(numCuts, make_pair(vector< pair<unsigned int, double> >(), 0));
		cuts.clear();

		int K=15;
		/*Flip[i]=1, if Constraint i has positive rhs (in the Ax<=b sense), and -1 if otherwise
		 * Flipping is done if Ax<=-b resulting in Ax-s<=b*/
		vector<int> Flip;

		/*Slack[i]=s, from flipped Constraint i*/
		vector<double> Slack;


		/*VectorIn[i]= (R+)-Variable, i.e. inflow in Constraint i, connecting Constraint i-1 with i*/
		vector<int> VectorIn;

		/*VectorOut[i]= (R-)-Variable, i.e. outflow from Constraint i, connecting Constraint i with i+1*/
		vector<int> VectorOut;

		/*Vector Containing Outflow-Variables (N-) in some Constraint; possibly becoming (R-)-Variables*/
		vector<int> VectorPotentialRMinus;

		/*Vector Containing Outflow-Variables (R-)-Variables in any Constraint*/
		vector<int> VectorRMinus;
		vector<int> VectorPotentialRMinusFails;

		/*InAt[i]= ConstraintNumber in which Variable i is an Inflow-Variable (in (N+) or (R+))*/
		vector<int> InAt;

		/*OutAt[i]= ConstraintNumber in which Variable i is an Outflow-Variable (in (N-) or (R-))*/
		vector<int> OutAt;

		/*InCoef[i]= Coefficient of Inflow-Variable i in Constraint InAt[i]*/
		vector<coef_t> InCoef;

		/*OutCoef[i]= Coefficient of Outflow-Variable i in Constraint OutAt[i]*/
		vector<coef_t> OutCoef;

		/*VectorCoefOut[i]=Coefficient of (R-)-Variable VectorOut[i] of Constraint i*/
		vector<coef_t> VectorCoefOut;

		/*Multiplier in order to have equal in- & outflow coefficients for succeeding Constraints*/
		vector<double> VectorRescale;

		/*Vector Constaining the Constraint of the built network*/
		vector<CRef> VectorOfConstraints;
		for (int i=0;i<nVars();i++){
			InAt.push_back(-1);
			OutAt.push_back(-1);
			InCoef.push_back(0);
			OutCoef.push_back(0);
		}
		for (int CurrentCandidate=0; CurrentCandidate <nVars();CurrentCandidate++) {
			for (int i=0;i<nVars();i++){
				InAt[i]=-1;
				OutAt[i]=-1;
				InCoef[i]=0;
				OutCoef[i]=0;
			}
			if(Print)cerr << "Candidate " << CurrentCandidate << "/" << nVars() << endl;// " " << candidates[CurrentCandidate] << endl;
			Flip.clear();
			Slack.clear();
			VectorIn.clear();
			VectorOut.clear();
			VectorRescale.clear();
			VectorOfConstraints.clear();
			VectorCoefOut.clear();

			VectorPotentialRMinus.clear();
			VectorRMinus.clear();
			VectorPotentialRMinusFails.clear();
			double remSumB=0.0;
			double SumB=0;
			int Attempts=0;
			int CurFlip=0;
			while(Attempts<VarsInConstraints[CurrentCandidate].size()){// ceil(log(VarsInConstraints[CurrentCandidate].size()))){
				int ConIn=irand(random_seed,VarsInConstraints[CurrentCandidate].size());
				CRef StartingConstraint =VarsInConstraints[CurrentCandidate][ConIn].cr;
				if (StartingConstraint==constraints[0]){
					Attempts++;
					continue;
				}
				Constraint &c = constraintallocator[StartingConstraint];
				CurFlip=(c.header.rhs>0)?-1:1;
				double CompVal=-1;
				int Out=-1;
				coef_t CoefOut;
				for (int i=0; i<c.size();i++){
					//if(type[var(c[i])]==BINARY) continue;
					//if(var(c[i])>=nVars() && Print) cerr << "Slack-Variable?"<<endl;
					if(((CurFlip==1 &&!sign(c[i]))||(CurFlip==-1 &&sign(c[i]))) && c[i].coef*objVals[var(c[i])].asDouble()>CompVal){
						Out=var(c[i]);
						CoefOut=c[i].coef;
					}
				}
				if (Out==-1) Attempts++;	//No variable with negative coefficients (positive Coefficients in constraintallocator)
				else{
					 if(CurFlip==1)
                                                Slack.push_back(0);
                                        else{
                                                double s= -c.header.rhs;
                                                for (int i=0; i<c.size();i++){
                                                        if(sign(c[i])){
                                                                s-=c[i].coef*objVals[var(c[i])].asDouble();
                                                        }
                                                        else{
                                                                s+=c[i].coef*objVals[var(c[i])].asDouble();
                                                        }

                                                }
                                                if (Print) cerr << "s " << s << endl;
                                        	if(s<-eps){
							continue;
						}
                                                Slack.push_back(s);
                                        }

					for (int i=0; i<c.size();i++){
						if(((CurFlip==1 &&!sign(c[i]))||(CurFlip==-1 &&sign(c[i])))){// if(!sign(c[i]) ){
							VectorPotentialRMinus.push_back(var(c[i]));
							OutAt[var(c[i])]=VectorOfConstraints.size();
							OutCoef[var(c[i])]=c[i].coef;

						}
					}
					VectorOfConstraints.push_back(StartingConstraint);
					VectorOut.push_back(Out);
					VectorIn.push_back(-1);
					Flip.push_back(CurFlip);
					VectorRescale.push_back(1);
					VectorCoefOut.push_back(CoefOut);
					remSumB=abs(c.header.rhs);
					assert(remSumB -c.header.rhs==0 || Flip.back()==1);
					if (PrintSum && Print) cerr << "remSum_0 " << -c.header.rhs << " "<< remSumB << endl;
					break;
				}
			}
			if(VectorOut.empty()){
				if(Print) cerr << "No initial Constraint Found" <<endl;
				continue;
			}


			while(VectorOfConstraints.size()<K && VectorOut.back()!=-1){
				if(Print) cerr << "N"<< endl;
				if(Print) cerr << "val " <<VectorOut.back() << " "<< VarsInConstraints[VectorOut.back()].size() << " " << ceil(log(VarsInConstraints[VectorOut.back()].size())) << " " << VectorOfConstraints.size()<<" " << K   << endl;

				// Lower Bounds = 0 required?
				double Left=0.0;
				double Right=0.0;
				//Select Next Constraint
				Attempts=0;
				int NumTries =VarsInConstraints[VectorOut.back()].size();//ceil(log(VarsInConstraints[VectorOut.back()].size()));
				while(Attempts< NumTries){
					int ConIn=irand(random_seed,VarsInConstraints[VectorOut.back()].size());
					if(Print)cerr << "Rand " << ConIn << endl;
					CRef NextConstraint =VarsInConstraints[VectorOut.back()][ConIn].cr;
					if (NextConstraint==constraints[0]){
						Attempts++;
						continue;
					}
					Constraint &c = constraintallocator[NextConstraint];
					CurFlip=(c.header.rhs>0)?-1:1;
					if((CurFlip==1&&!sign(c[VarsInConstraints[VectorOut.back()][ConIn].pos]))||(CurFlip==-1&&sign(c[VarsInConstraints[VectorOut.back()][ConIn].pos]))){//if(!sign(c[VarsInConstraints[VectorOut.back()][ConIn].pos])){
						//Previous Out Variable (and thus now mandatory In-Variable) is not an In variable (i.e. has positive sign) in the new constraint
						if(Print)cerr << "Wrong Constraint 'Input'" << endl;
						Attempts++;
						continue;
					}
					if(std::find(VectorOfConstraints.begin(), VectorOfConstraints.end(), NextConstraint) != VectorOfConstraints.end()){
						if(Print)cerr <<"Already part of network " << endl;
						Attempts++;
						continue;
					}
					double CompVal=-1;
					int Out=-1;
					coef_t CoefOut;
					coef_t CoefIn=-1;
					for (int i=0; i<c.size();i++){
						if(((CurFlip==1 &&sign(c[i]))||(CurFlip==-1 &&!sign(c[i])))){// if(sign(c[i]) ){
							//This variable occurs Positive in the considered Constraint

							if(std::find(VectorRMinus.begin(), VectorRMinus.end(), var(c[i])) != VectorRMinus.end() ){
								//Already Inflow in other constraint -> not good / Not sure what happens
								Out=-2;
								break;
							}
						}
						else if(std::find(VectorPotentialRMinusFails.begin(), VectorPotentialRMinusFails.end(), var(c[i])) != VectorPotentialRMinusFails.end()){
							//Twice Marked as possible Outflow -> not good
							Out=-2;
							break;
						}
						if(((CurFlip==1 &&!sign(c[i]))||(CurFlip==-1 &&sign(c[i])))){// if(!sign(c[i]) ){
							for (int p=0;p<VarsInConstraints[var(c[i])].size();p++){
								for (int ch=0;ch<VectorOfConstraints.size();ch++){
									if (VectorOfConstraints[ch]==VarsInConstraints[var(c[i])][p].cr){
										Constraint &testc = constraintallocator[VarsInConstraints[var(c[i])][p].cr];
										int posi = VarsInConstraints[var(c[i])][p].pos;
										if(((Flip[ch]==1 &&sign(testc[posi]))||(Flip[ch]==-1 &&!sign(testc[posi])))){// if (sign(testc[posi])){
											Out=-2;
											if (Print) cerr <<"Break cause of  backward arc" <<endl;
										}
									}
								}

							}
						}
						if(Out==-2) break;
						if(((CurFlip==1 &&!sign(c[i]))||(CurFlip==-1 &&sign(c[i])))&& c[i].coef*objVals[var(c[i])].asDouble()>CompVal){
							Out=var(c[i]);
							CoefOut=c[i].coef;
						}
						if(var(c[i])==VectorOut.back()) CoefIn=c[i].coef;
					}
					/*if (Out==-1 && VectorOfConstraints.size()<K-1 && VectorOfConstraints.size()<3){
						//No variable with negative coefficients=> No potential Outflow => Network too small :(
						if (Print) cerr << "No Out Variables " << VectorOfConstraints.size()<< endl;
						Attempts++;
						//break;

					}
					else*/ if(Out==-2){
						if(Print) cerr << "Constraint contains already (R+)-marked variables " << endl;
						Attempts++;
					}
					/*else if (remSumB-c.header.rhs<0){
						Attempts++;
					}*/
					else{
						for (int i=0; i<c.size();i++){
							//cerr << (sign(c[i])?"-":"+")<<c[i].coef<<"x_"<<var(c[i])<<"("<<objVals[var(c[i])].asDouble()<<")";
							if(((CurFlip==1 &&sign(c[i]))||(CurFlip==-1 &&!sign(c[i])))){//if(sign(c[i]) ){
								//positive Variable
								if (std::find(VectorPotentialRMinus.begin(), VectorPotentialRMinus.end(), var(c[i])) != VectorPotentialRMinus.end()){
									//Variable is Outflow (-) in some other Constraint and is now an Inflow (+) => Variable moved from N_ to R+
									VectorRMinus.push_back(var(c[i]));
									InAt[var(c[i])]=VectorOfConstraints.size();	//Correct, as the actual constraint is pushed back later
									InCoef[var(c[i])]=(double)VectorCoefOut.back()/(double)CoefIn*c[i].coef;

									if(InCoef[var(c[i])]!=OutCoef[var(c[i])]){
										Attempts=NumTries;
										break;
									}
									if(Print) cerr << "Added RMinus Element " << endl;
								}
							}
							else{
								assert((CurFlip==1 &&!sign(c[i]))||(CurFlip==-1 &&sign(c[i])));
								//negative Variable
								if( std::find(VectorPotentialRMinus.begin(), VectorPotentialRMinus.end(), var(c[i])) != VectorPotentialRMinus.end()){
									// Variable was in previous constraint an Outflow and now again => possible trouble
									VectorPotentialRMinusFails.push_back(var(c[i]));
									if(std::find(VectorRMinus.begin(), VectorRMinus.end(), var(c[i])) == VectorRMinus.end()){
										assert(OutAt[var(c[i])]!=-1);
										OutAt[var(c[i])]=-1;

									}
								}
								else{
									VectorPotentialRMinus.push_back(var(c[i]));
									OutAt[var(c[i])]=VectorOfConstraints.size();
									OutCoef[var(c[i])]=(double)VectorCoefOut.back()/(double)CoefIn*c[i].coef;

								}
							}
						}
						if(CoefIn<0) continue;
						assert(Attempts==NumTries || std::find(VectorRMinus.begin(), VectorRMinus.end(), VectorOut.back()) != VectorRMinus.end());
						if(Print)cerr << "Added Constraint " << VectorOfConstraints.size()+1<< endl;
						VectorOfConstraints.push_back(NextConstraint);
						VectorIn.push_back(VectorOut.back());
						VectorOut.push_back(Out);

						VectorRescale.push_back((double)VectorCoefOut.back()/(double)CoefIn);
						Flip.push_back(CurFlip);
						if(CurFlip==1)
							Slack.push_back(0);
						else{
							double s= -VectorRescale.back()*c.header.rhs;
							for (int i=0; i<c.size();i++){
								if(sign(c[i])){
									s-=VectorRescale.back()*c[i].coef*objVals[var(c[i])].asDouble();
								}
								else{
									s+=VectorRescale.back()*c[i].coef*objVals[var(c[i])].asDouble();
								}

							}
							if(s<-eps && Print){
								cerr << "Negative Slack!? " << s << " " <<VectorRescale.back()<< endl;
								for (int i=0; i<c.size();i++){
									 cerr << (sign(c[i])?"-":"+")<<c[i].coef<<"x_"<<var(c[i])<<"("<<objVals[var(c[i])].asDouble()<<")";
								}
								cerr << ">= "<< c.header.rhs << endl;
							}
							if(s>=-eps) {
							} else {
								Attempts=NumTries;
								continue;
							}
							assert(s>=-eps);
							Slack.push_back(s);
						}
						if (VectorRescale.size() <= 0) {
						  if(getShowWarning()) cerr << "Warning: VectorRescale.size() <= 0" << endl;
						  continue;
						}
						if (VectorRescale.back() <= 0) {
						  if(getShowWarning()) cerr << "Warning: VectorRescale.back() <= 0" << endl;
						  return 0;
						}
						assert(VectorRescale.back() > 0);
						VectorCoefOut.push_back((double)VectorRescale.back()*(double)CoefOut);
						if(abs((double)VectorRescale.back()*c.header.rhs)>eps)
							remSumB+=abs(VectorRescale.back()*c.header.rhs);
						if (PrintSum&& Print )cerr << "remSum_" <<VectorIn.size() << "="<< remSumB << " " << -VectorRescale.back()*c.header.rhs << " "<< -(double)VectorRescale.back()*(double)c.header.rhs << " "<<VectorRescale.back() << " " << c.header.rhs <<endl;

						for (int r=0;r<VectorOfConstraints.size();r++){
							if(Print)cerr << "Constraint " << r << "/" << VectorOfConstraints.size() << "with Flip:" << Flip.back()<<endl;
							Constraint & c = constraintallocator[VectorOfConstraints[r]];
							for (int i=0; i<c.size();i++){
								if(Print) cerr << (sign(c[i])?"-":"+")<<VectorRescale[r]*c[i].coef<<"x_"<<var(c[i])<<"("<<objVals[var(c[i])].asDouble()<<")";
							}
							if(Print)cerr << ">= "<< VectorRescale[r]*c.header.rhs << endl;
							if(Print)cerr << "In " << VectorIn[r] << endl;
							if(Print)cerr << "Out " << VectorOut[r] << endl;
						}
						break;
					}
				}
				if(Print)cerr <<"B: " << Attempts <<" " << VectorOfConstraints.size()<<endl;
				if (Attempts==NumTries){
					if(Print)cerr << "Failed to find another constraint " << VarsInConstraints[VectorOut.back()].size() << endl;
					break;}
				if(Print)cerr <<"A: " << NumTries << endl;

				SumB=remSumB;
				vector<pair<int,int>> Cplus;
				Cplus.clear();
				if (PrintSum && Print)cerr  <<"Start " << SumB << endl;

				//Check weather violation can occur
				double 	SumY=0.0;
				double SumU=0.0;

				for (int k=0;k<VectorOfConstraints.size();k++){
					if(Flip[k]==-1) Right+=Slack[k];
					Constraint & c = constraintallocator[VectorOfConstraints[k]];
					if(Strengthening && VectorOut[k]!=-1){
						if (type[VectorOut[k]]==CONTINUOUS && UpperBoundVar[VectorOut[k]]!=-1)
							SumY+=objVals[UpperBoundVar[VectorOut[k]]].asDouble();
						else if(type[VectorOut[k]]==CONTINUOUS)
							SumY+=1;
						else{
							assert(type[VectorOut[k]]==BINARY);
							SumY+=objVals[VectorOut[k]].asDouble();
						}
					}
					for (int i=0; i<c.size();i++){
						if(abs(objVals[var(c[i])].asDouble())<eps) continue;

						if((k!=VectorOfConstraints.size()-1 && var(c[i])==VectorOut[k]) || var(c[i])==VectorIn[k] ||
								(std::find(VectorRMinus.begin(), VectorRMinus.end(), var(c[i])) != VectorRMinus.end() && ((((Flip[k]==1 &&!sign(c[i]))||(Flip[k]==-1 &&sign(c[i]))) && OutAt[var(c[i])]==k && k!=VectorOfConstraints.size()-1)||(((Flip[k]==1 &&sign(c[i]))||(Flip[k]==-1 &&!sign(c[i])))&& InAt[var(c[i])]==k) )) ) continue;
						if(((Flip[k]==1 &&!sign(c[i]))||(Flip[k]==-1 &&sign(c[i])))){
							//negative Variable
							if(Strengthening&& (double)c[i].coef*objVals[var(c[i])].asDouble()>SumY*(double)c[i].coef*(double)upperBounds[var(c[i])] ){ // Rescale on both sides!?

								assert(upperBounds[var(c[i])]>=0);
								if(Print)cerr << "Strengthening possible? " << objVals[var(c[i])].asDouble() << " > " <<SumY*(double)upperBounds[var(c[i])]  << " " << SumY << " "<<(double)upperBounds[var(c[i])] <<endl;
								SumU+=(double)VectorRescale[k]*upperBounds[var(c[i])]*(double)c[i].coef ;
								if(Print)cerr << "SUMU+ " << (double)upperBounds[var(c[i])]*(double)c[i].coef << " " << k << " " <<var(c[i]) << endl;
							}
							if(Print)cerr << "Add Right " << k << " " << var(c[i]) << " " <<c[i].coef << endl;
							Right+=VectorRescale[k]*c[i].coef*objVals[var(c[i])].asDouble();
							if(objVals[var(c[i])].asDouble()<0){
								if(Print)cerr <<"Kleiner als 0 " << objVals[var(c[i])].asDouble() << " Var x_" <<var(c[i]) <<endl;
							}
						}
						else if (type[var(c[i])]==BINARY){
							//positive (N+)
							double v=VectorRescale[k]*(double)c[i].coef*objVals[var(c[i])].asDouble()-SumB*objVals[var(c[i])].asDouble();
							if(v>eps){
								if(Print)cerr << "Add LeftB " << k << " " << var(c[i]) << " " <<c[i].coef << endl;
								Left+=v;
								Cplus.push_back(make_pair(k,i));
							}

						}
						else{
							assert(type[var(c[i])]==CONTINUOUS);
							double v;
							double backpart;
							if(UpperBoundVar[var(c[i])]!=-1)
								backpart=objVals[UpperBoundVar[var(c[i])]].asDouble()*SumB;
							else
								backpart=SumB;

							v=VectorRescale[k]*(double)c[i].coef*objVals[var(c[i])].asDouble()-backpart;

							if(v>eps){
								if(Print)cerr << "Add LeftC " << k << " " << c[i].coef << " x_" << var(c[i])  <<"="<<objVals[var(c[i])].asDouble()<< " " << v << " " << SumB  << " " << VectorRescale[k]<<endl;
								if ( UpperBoundVar[var(c[i])]!=-1 && Print) cerr << "Info: x_" << UpperBoundVar[var(c[i])]<<"="<<objVals[UpperBoundVar[var(c[i])]].asDouble() << " " << backpart << " " << v <<  endl;
								Left+=v;
								Cplus.push_back(make_pair(k,i));

							}
						}
					}
					if(abs((double)VectorRescale[k]*c.header.rhs)>eps)
						SumB-=abs(VectorRescale[k]*(double)c.header.rhs);
					if (PrintSum && Print ) cerr <<"Size "<<VectorRescale.size() << " " << VectorRescale.back() << endl;
					if (PrintSum && Print) cerr<< k <<" " << VectorRescale[k]*(double)c.header.rhs << " " <<  SumB << " " <<VectorRescale[k] << " " << c.header.rhs << " " << VectorOut[k] << " " <<VectorIn[k] << endl;
				}
				if (PrintSum && Print) cerr<<"Final " << SumB << endl;
				if (abs(SumB)>=eps){
                                  Attempts=NumTries;
				  if(Print) cerr << "Warning: epsilon assertion avoided. " << endl;
				  break;
				}

				assert(abs(SumB)<eps);
				if(Print) cerr << "LR " << Left << " " << Right << endl;
				double SaveSumY=SumY;
				if(Left>Right && Right>=0){
					assert(!Cplus.empty());
					assert(!VectorOfConstraints.empty());
					if(Print)cerr <<"Cut Found! "<< Left << " " << Right <<endl;
					int t=0;
					double final_rhs=0.0;
					SumB=remSumB;

					//PRINT START
					for (int r=0;r<VectorOfConstraints.size();r++){
						if(Print)cerr << "Constraint " << r << "/" << VectorOfConstraints.size() << " with Flip" << Flip[r] << endl;
						Constraint & c = constraintallocator[VectorOfConstraints[r]];
						for (int i=0; i<c.size();i++){
							if(Print) cerr << (sign(c[i])?"-":"+")<<VectorRescale[r]*c[i].coef<<"x_"<<var(c[i])<<"("<<objVals[var(c[i])].asDouble()<<"/"<<(type[var(c[i])]==BINARY?"B/":"C") << UpperBoundVar[var(c[i])]<<")";
						}
						if(Print)cerr << ">= "<< VectorRescale[r]*c.header.rhs << endl;
						if(Print)cerr << "In " << VectorIn[r] << endl;
						if(Print)cerr << "Out " << VectorOut[r] << endl;
					}
					//PRINT End

					vector<pair<unsigned int, double> > cutvec;
					SumY=0.0;
					double MinusSumU=0.0;
					for (int i=0;i<Cplus.size();i++){
						while(t!=Cplus[i].first){
							MinusSumU=0.0;
							if(Strengthening && VectorOut[t]!=-1){
								if (type[VectorOut[t]]==CONTINUOUS && UpperBoundVar[VectorOut[t]]!=-1)
									SumY+=objVals[UpperBoundVar[VectorOut[t]]].asDouble();
								else if(type[VectorOut[t]]==CONTINUOUS)
									SumY+=1;
								else{
									assert(type[VectorOut[t]]==BINARY);
									SumY+=objVals[VectorOut[t]].asDouble();
								}
							}
							Constraint & c = constraintallocator[VectorOfConstraints[t]];

							for(int j=0;j<c.size();j++){
								if(((Flip[t]==1 &&!sign(c[j]))||(Flip[t]==-1 &&sign(c[j]))) && (var(c[j])!=VectorOut[t]|| t==VectorOfConstraints.size()-1 )
										//&& (std::find(VectorRMinus.begin(), VectorRMinus.end(), var(c[j])) == VectorRMinus.end()|| (OutAt[var(c[j])]<t/* && t!=VectorOfConstraints.size()-1))*/
										){

									if(std::find(VectorRMinus.begin(), VectorRMinus.end(), var(c[j])) != VectorRMinus.end()) continue;
									//(N-)-Variable, in particular not (R-)

									if(Print) cerr << " -" << VectorRescale[t]*c[j].coef<< "x_"<<var(c[j])<<"("<<objVals[var(c[j])].asDouble()<<")";
									//Strengthening

									if(Strengthening&&(double)c[j].coef*objVals[var(c[j])].asDouble()>SumY*(double)c[j].coef*(double)upperBounds[var(c[j])]  && objVals[var(c[j])].asDouble()>eps ){
									  if (info_level > 2) cerr<<"-S";
										MinusSumU+=VectorRescale[t]*(double)upperBounds[var(c[j])]*(double)c[j].coef;
												if(Print){cerr << "MINUSSUMU+ " << (double)upperBounds[var(c[j])]*(double)c[j].coef << " " << t << " "<< var(c[j]) << endl;

										//cerr << "Strengthening possible? " << objVals[var(c[j])].asDouble() << " > " <<SumY*(double)upperBounds[var(c[j])]  << " " << SumY << " "<<(double)upperBounds[var(c[j])] <<endl;
												}
									}
									else{
							      cutvec.push_back( make_pair( var(c[j]), RoundIfClose(+VectorRescale[t]*c[j].coef,eps) ) );
									}
								}

							}
							if(abs((double)VectorRescale[t]*c.header.rhs)>eps)
								SumB-=abs(VectorRescale[t]*(double)c.header.rhs);
							t++;
							if(Strengthening)SumU-=MinusSumU;
						}
						//if(abs(SumB - floor(SumB+.5))<eps) SumB=floor(SumB+.5);
						assert(t==Cplus[i].first);
						Constraint & c = constraintallocator[VectorOfConstraints[t]];
						assert(sign(c[Cplus[i].second])||Flip[t]==-1);
						assert(var((c[Cplus[i].second]))!=VectorOut[t]&&var((c[Cplus[i].second]))!=VectorIn[t]);
						if(Print) cerr << " +" << VectorRescale[t]*c[Cplus[i].second].coef<< "x_"<<var(c[Cplus[i].second])<<"("<<objVals[var(c[Cplus[i].second])].asDouble()<<")";
					      cutvec.push_back( make_pair( var(c[Cplus[i].second]), RoundIfClose(-VectorRescale[t]*c[Cplus[i].second].coef,eps) ) );

						if(type[var(c[Cplus[i].second])]==CONTINUOUS && UpperBoundVar[var(c[Cplus[i].second])]==-1){
							final_rhs+=SumB+SumU;
						}
						else if(type[var(c[Cplus[i].second])]==CONTINUOUS ){
							if(Print) cerr << " -" << SumB+SumU<< "x_"<<UpperBoundVar[var(c[Cplus[i].second])]<<"("<<objVals[UpperBoundVar[var(c[Cplus[i].second])]].asDouble()<<")";
						      cutvec.push_back( make_pair(UpperBoundVar[var(c[Cplus[i].second])], RoundIfClose(SumB+SumU,eps)) );

						}
						else if (type[var(c[Cplus[i].second])]==BINARY){
							if(Print) cerr << " -" << SumB+SumU<< "x_"<<var(c[Cplus[i].second])<<"("<<objVals[var(c[Cplus[i].second])].asDouble()<<")";
						      cutvec.push_back( make_pair(var(c[Cplus[i].second]), RoundIfClose(SumB+SumU,eps)) );

						}
						else assert(0);
					}

					MinusSumU=0.0;
					while(t<VectorOfConstraints.size()){
						Constraint & c = constraintallocator[VectorOfConstraints[t]];
						if(Strengthening && VectorOut[t]!=-1){
							if (type[VectorOut[t]]==CONTINUOUS && UpperBoundVar[VectorOut[t]]!=-1)
								SumY+=objVals[UpperBoundVar[VectorOut[t]]].asDouble();
							else if(type[VectorOut[t]]==CONTINUOUS)
								SumY+=1;
							else{
								assert(type[VectorOut[t]]==BINARY);
								SumY+=objVals[VectorOut[t]].asDouble();
							}
						}
						for(int j=0;j<c.size();j++){
							if(((Flip[t]==1 &&!sign(c[j]))||(Flip[t]==-1 &&sign(c[j]))) && (var(c[j])!=VectorOut[t]|| t==VectorOfConstraints.size()-1 )

									//&& (std::find(VectorRMinus.begin(), VectorRMinus.end(), var(c[j])) == VectorRMinus.end()|| (OutAt[var(c[j])]<t/* && t!=VectorOfConstraints.size()-1))*/
									){
								if(std::find(VectorRMinus.begin(), VectorRMinus.end(), var(c[j])) != VectorRMinus.end()) continue;
								//(N-); not (R-)
								if(Print) cerr << " -" << VectorRescale[t]*c[j].coef<< "x_"<<var(c[j])<<"("<<objVals[var(c[j])].asDouble()<<")";
								if(Strengthening&&(double)c[j].coef*objVals[var(c[j])].asDouble()>SumY*(double)c[j].coef*(double)upperBounds[var(c[j])] && objVals[var(c[j])].asDouble()>eps ){
									MinusSumU+=VectorRescale[t]*(double)upperBounds[var(c[j])]*(double)c[j].coef;
									if(Print){cerr << "MINUSSUMU(!)+ " << (double)upperBounds[var(c[j])]*(double)c[j].coef << " " << t << " " <<var(c[j]) << endl;

									//cerr << "Strengthening possible!? " << objVals[var(c[j])].asDouble() << " > " <<SumY*(double)upperBounds[var(c[j])]  << " " << SumY << " "<<(double)upperBounds[var(c[j])] <<endl;
									}
								}
								else
									cutvec.push_back( make_pair( var(c[j]), RoundIfClose(VectorRescale[t]*c[j].coef,eps) ) );
							}

						}
						if(abs((double)VectorRescale[t]*c.header.rhs)>eps)
							SumB-=abs(VectorRescale[t]*(double)c.header.rhs);
						t++;
					}
					if(abs(SumY-SaveSumY)>eps) {
					  //if(getShowError()) cerr << "Error SumY " << SumY << " "<< SaveSumY << endl;
					  return 0;
					}
					assert (abs(SumY-SaveSumY)<=eps);
					if(abs(SumU-MinusSumU)>eps) {
					  //if(getShowError()) cerr << "Error SumU " << SumU << " "<< MinusSumU << endl;
					  return 0;
					}
					assert(abs(SumU-MinusSumU)<=eps);
					for (int k=0;k<VectorOfConstraints.size();k++){
						if(Flip[k]==-1){
							Constraint & c = constraintallocator[VectorOfConstraints[k]];
							final_rhs-=VectorRescale[k]*c.header.rhs;
							for (int i=0;i<c.size();i++){
								//Vorzeichen umdrehen! Daher -/+
								if(sign(c[i])){
									if(Print) cerr << " +" << VectorRescale[k]*c[i].coef<< "x_"<<var(c[i])<<"("<<objVals[var(c[i])].asDouble()<<";" << cutvec.size() <<")";
									cutvec.push_back( make_pair( var(c[i]), RoundIfClose(-VectorRescale[k]*c[i].coef,eps) ) );
								}
								else{
									if(Print) cerr << " -" << VectorRescale[k]*c[i].coef<< "x_"<<var(c[i])<<"("<<objVals[var(c[i])].asDouble()<<";" << cutvec.size() <<")";
									cutvec.push_back( make_pair( var(c[i]), RoundIfClose(VectorRescale[k]*c[i].coef,eps) ) );
								}

							}

						}
					}
					assert (t==VectorOfConstraints.size());
					/*for (int r=0;r<VectorOfConstraints.size();r++){
						Constraint & c = constraintallocator[VectorOfConstraints[r]];
						for (int i=0; i<c.size();i++){
							//RICHTIG=????????
							if(!sign(c[i]) &&
									(std::find(VectorRMinus.begin(), VectorRMinus.end(), var(c[i])) == VectorRMinus.end() || ((OutAt[var(c[i])]<r && r!=VectorOfConstraints.size()-1)))){
								cerr << " -" << "x_"<<var(c[i])<<"("<<objVals[var(c[i])].asDouble()<<")";

							}
						}

					}*/
					final_rhs+=+1e-12 +abs(final_rhs)*1e-12;
					if(Print) cerr << " <= " <<  final_rhs << "(" << cutvec.size() << ")" << endl ;
					cuts.push_back(make_pair(cutvec, -final_rhs));

					break;
				}
				else{
					if(Print) cerr <<"Try to Add a New Constraint" << endl;
				}
			}
		}

		/*listOfCutsLhs.reserve( cuts.size() );
			listOfCutsRhs.reserve( cuts.size() );

			for( unsigned int i = 0; i < cuts.size(); i++ ){
				vector< data::IndexedElement > addcut;
				addcut.clear();
				double rrhs=cuts[i].second;
				addcut.reserve( addcut.size()+cuts[i].first.size() );
							for( unsigned int jj = 0; jj < cuts[i].first.size(); ++jj ){
								addcut.push_back( data::IndexedElement( cuts[i].first[jj].first, data::QpNum( cuts[i].first[jj].second ) ) );
							}
			    if (addcut.size()>0) {
					listOfCutsLhs.push_back( addcut );
					listOfCutsRhs.push_back( rrhs );
				}

			}

		return listOfCutsRhs.size();*/

		listOfCutsLhs.reserve( cuts.size()+listOfCutsLhs.size() );
		listOfCutsRhs.reserve( cuts.size()+listOfCutsLhs.size() );

		if (info_level > 2) cerr << "CUTS OUT " << endl;
		for( unsigned int i = 0; i < cuts.size(); ++i ){
			vector< data::IndexedElement > addcut;
			addcut.clear();
			int j =0;
			double rrhs=0.0;

			for ( ; j < 1 && i+j < cuts.size();j++) {
				if(Print){
					cerr << "Round " <<i << " "<<j<<endl;
					for (int pp=0;pp<cuts.at( i+j ).first.size();pp++){
						cerr << cuts.at( i+j ).first[pp].second <<"x_"<<cuts.at( i+j ).first[pp].first << endl;
					}
					cerr << ">= " << cuts.at( i+j ).second << endl;
				}
				vector< pair<unsigned int, double> >& cut = cuts.at( i+j ).first;

				bool found=false;
				if (computeCutRatio1(cut) > 1000000.0) continue;

				rrhs = rrhs + cuts[i+j].second;
				addcut.reserve( addcut.size()+cut.size() );
				for( unsigned int jj = 0; jj < cut.size(); ++jj ){
					addcut.push_back( data::IndexedElement( cut.at( jj ).first, data::QpNum( cut.at( jj ).second ) ) );
				}
				std::sort(addcut.begin(), addcut.end(), [](data::IndexedElement e1, data::IndexedElement e2){return e1.index < e2.index;});

				// run through and add
				int ii=0, jj=0;
				for (ii = 0; ii < addcut.size();ii++) {
					jj = ii+1;
					while (jj < addcut.size() && addcut[ii].index == addcut[jj].index) {
						addcut[ii].value += addcut[jj].value;
						addcut[jj].value = 0.0;
						jj++;
					}
					ii = jj-1	;
				}
				// delete 0s
				ii=0, jj=0;
				for (ii = 0; ii < addcut.size();ii++) {
					if (abs((addcut[ii].value.asDouble()))>eps) continue;
					jj = addcut.size()-1;
					while (jj>=ii && abs((addcut[jj].value.asDouble()))<eps) {
						jj--;
						addcut.pop_back();
					}
					if (ii >= jj) break;
					else {
						addcut[ii] = addcut[jj];
						addcut[jj].value = 0.0;
						addcut.pop_back();
					}
				}
				if (!found) break;
			}
			i += j;
			j=0;
		    if (addcut.size()>0) {
				listOfCutsLhs.push_back( addcut );
				listOfCutsRhs.push_back( rrhs );
				if(Print){for (int pp=0;pp<addcut.size();pp++){
					cerr << addcut[pp].value.asDouble() <<"x_"<<addcut[pp].index << endl;
				}
				cerr << ">= " << rrhs << endl;
				}
			}

		}
		if (info_level > 2) cerr << cuts.size() << " " << listOfCutsRhs.size() << endl;
//cin.get();
		return listOfCutsRhs.size(); //number of cuts

	}

	double QBPSolver::computeCutRatio1(vector< pair<unsigned int, double> >& cut) {
		if (cut.size() < 1) return 10000000.0;
		double mx=abs(cut[0].second);
		double mn=mx;
	    for (int k=1; k < cut.size();k++) {
	    	if (abs(cut[k].second) > mx) mx = cut[k].second;
	    	if (abs(cut[k].second) < mn && abs(cut[k].second) > 0) mn = cut[k].second;
	    }
	    return mx / mn;
	}

	double QBPSolver::RoundIfClose(double Val, double EPS){
		//return Val;
		if(abs(floor(Val+0.5)-Val)<=EPS) return floor(Val+0.5);
		else return Val;
	}

	double QBPSolver::ConstraintScore(const CRef& Cons,const int& n,const double& slackOfThis, const double& RowNorm, const double& eps){
		if (Cons==constraints[0]) return -1;
		Constraint &c = constraintallocator[Cons];

		//Constraints that are Varialbe Lower- or Upper-Bound Constraints or Constraints containing integer-to-binary translation, e.g. 4x+2y+z<=6, are not of interest
		if(c.header.isIntBnd || c.header.isVarBnd )
			return -1;
		if(0&&c.header.learnt/* && !InputConstraintOk(c, eps*10, 1000, 30)*/){
			//cerr << "Not OK Learnt Constraint" << endl;
			return -1;
		}
		//return drand(random_seed);
		double density=c.size()/(double)n;
		double SlackPart=slackOfThis/(double)(max(RowNorm,0.1));
		double ret=0.0001*(1-density)+0.001*(1-SlackPart);

		// No dual information available -> If row is tight we assume a random dual value of this constraint
		if(abs(slackOfThis)<eps) ret+=max(drand(random_seed)/1000,0.0001);
		else ret+=0.0001;

		//cerr << "Ret: " << ret << endl;
		//assert(ret>=0);
		return pow(0.9,c.header.usedinAgg)*ret;
	}

	bool QBPSolver::ContCountOK(const Constraint& c,const int& minC, const int& maxC,const vector<data::QpNum>& Solution, const double& eps){
		//Count reals strictly between bounds
		int count=0;
		for (int i=0;i<c.size();i++){
		    if(var(c[i])>=nVars()){/*cerr << "Too Large Index"<<endl;*/ return false;}
				if(type[var(c[i])] == CONTINUOUS && abs(Solution[var(c[i])].asDouble()-lowerBounds[var(c[i])])>eps && abs(Solution[var(c[i])].asDouble()-upperBounds[var(c[i])])>eps )
					count++;
		}
		return (count<=maxC && count >=minC);
	}

	bool QBPSolver::ContCountOK(vector<data::IndexedElement> &Row,const int& minC, const int& maxC,const vector<data::QpNum>& Solution, const double& eps, int FirstSlack){
		//Count reals strictly between bounds
		int count=0;
		for (int i=0;i<Row.size();i++){
		  int index = Row[i].index;
		  //assert(index < nVars());
		  if (index>=FirstSlack) continue; //Ignore Slacks here
		  if (index < nVars());
		  else {
		    //return false;
		    //cerr << "i" << index << ", " << type.size() << ", " << nVars() << ", " << resizer.v_ids.size() << endl;
		    index = resizer.getShadowProjection(index);
		  }
		  if (index==nVars()) return false;
		  assert(index>=0 && index <nVars());
		  if(type[index] == CONTINUOUS && abs(Solution[index].asDouble() - lowerBounds[index])>eps &&abs(Solution[index].asDouble() - upperBounds[index])>eps)
		  //if(type[Row[i].index] == CONTINUOUS && abs(Solution[Row[i].index].asDouble() - lowerBounds[Row[i].index])>eps &&abs(Solution[Row[i].index].asDouble() - upperBounds[Row[i].index])>eps)
					count++;
		}
		return (count<=maxC && count >=minC);
}


	bool QBPSolver::ClearDoubleEntries(vector<data::IndexedElement>& Constraint, double eps){
		std::sort(Constraint.begin(), Constraint.end(), [](data::IndexedElement e1, data::IndexedElement e2){return e1.index < e2.index;});
		// run through and add
		int ii=0, jj=0;
		for (ii = 0; ii < Constraint.size();ii++) {
			jj = ii+1;
			while (jj < Constraint.size() && Constraint[ii].index == Constraint[jj].index) {

				Constraint[ii].value += Constraint[jj].value;
				Constraint[jj].value = 0.0;
				jj++;
			}
			ii = jj-1	;
		}
		// delete 0s
		ii=0, jj=0;
		for (ii = 0; ii < Constraint.size();ii++) {

			if (abs((Constraint[ii].value.asDouble()))>eps) continue;
			jj = Constraint.size()-1;
			while (jj>=ii && abs((Constraint[jj].value.asDouble()))<eps) {
				jj--;
				Constraint.pop_back();
			}
			if (ii >= jj) break;
			else {
				Constraint[ii] = Constraint[jj];
				Constraint[jj].value = 0.0;
				Constraint.pop_back();
			}
		}
		return true;
	}

	double QBPSolver::F(double alpha, double d){
		assert(alpha>=0 && alpha<1);
		return floor(d)+fmax(0,d-floor(d)-alpha)/(1-alpha);
	}
	double QBPSolver::FBar(double alpha, double d){
			assert(alpha>=0 && alpha<1);
			return fmin(0,d/(1-alpha));
		}


	bool QBPSolver::CollectIntegerInfos(vector<data::IndexedElement>& Row,vector<OrgInfo>& USet, vector<OrgInfo>& TSet, vector<data::QpNum>& sol, vector<double>& NStar, int FirstSlack, const double& eps){
		//WARUM? In ClearDoubleEntries beibehalten?
				std::sort(Row.begin(), Row.end(), [](data::IndexedElement e1, data::IndexedElement e2){return e1.index < e2.index;});
			double LargestCoef=0;
			for (int i=0;i<Row.size();i++){
			assert(i==Row.size()-1 || Row[i].index< Row[i+1].index);	//Assert that the vector is sorted and no double entries
			if(Row[i].index<FirstSlack&& type[Row[i].index]!=CONTINUOUS ){
				OrgInfo Info;
				if(Row[i].index!=((yInterface*)yIF)->integers[Row[i].index].pt2leader) return false;
				Info.StartIndex=Row[i].index;

				//Assert that upper bound is integer
				//cerr <<((yInterface*)yIF)->integers[Row[i].index].org_ub.asDouble()<< " BB " << ((yInterface*)yIF)->integers[Row[i].index].org_lb.asDouble() << endl;

				//assert(abs(((yInterface*)yIF)->integers[Row[i].index].org_ub.asDouble()-((yInterface*)yIF)->integers[Row[i].index].org_lb.asDouble())-floor((((yInterface*)yIF)->integers[Row[i].index].org_ub.asDouble()-((yInterface*)yIF)->integers[Row[i].index].org_lb.asDouble()+.5))<eps);
				Info.Bound=((yInterface*)yIF)->integers[Row[i].index].org_ub-((yInterface*)yIF)->integers[Row[i].index].org_lb;
				if(((yInterface*)yIF)->integers[Row[i].index].number!=-1 && Binarized[((yInterface*)yIF)->integers[Row[i].index].number].ub <Info.Bound.asDouble()){
					//cerr <<"FOUND BETTER BOUND!!!!!!" << endl;
					Info.Bound =  Binarized[((yInterface*)yIF)->integers[Row[i].index].number].ub;
				}
				Info.XStar=0;
				bool AsBool=false;
				for(int k=0;k<((yInterface*)yIF)->integers[Row[i].index].bitcnt;k++){
					// Assert that all binaries corresponding to the single integer variable are present
					// Because otherwise not clear yet what to do
					  if(i+k>=Row.size()||(Row[i].index+k>=nVars()) || (Row[i+k].index!=Row[i].index+k) ||
					    ((yInterface*)yIF)->integers[Row[i].index].pt2leader!=
					    ((yInterface*)yIF)->integers[Row[i+k].index].pt2leader){
					    //cerr << "Warning: Binarized Var of Integer appear not grouped" << endl;
					    //break;
					    AsBool=true;
					    break;
					    //return false;
					  }
					if(k>0 && abs(Row[i+k-1].value.asDouble()-2*Row[i+k].value.asDouble())>eps){
						//cerr <<"Ups. Found binarized integer variable with wacky coefficients." << endl;
						AsBool=true;
						break;
						//return false;
					}
					if(((yInterface*)yIF)->integers[Row[i].index].pt2leader!=((yInterface*)yIF)->integers[Row[i+k].index].pt2leader) return false;

					//Calculate XStar
					Info.XStar+=sol[Row[i].index+k].asDouble()*pow(2,((yInterface*)yIF)->integers[Row[i].index].bitcnt-1-k);
				}
				if(!AsBool){	//Interpret Ints as Int and Bool as Bool
					Info.SeeAsBool=false;
					if(((yInterface*)yIF)->integers[Row[i].index].bitcnt<1){
						if(getShowError()) cerr <<"Error: Wrong BitCount?? "<<i << " " << Row[i].index << " " << ((yInterface*)yIF)->integers[Row[i].index].bitcnt << endl;
						//cin.get();
					}
					i=i+((yInterface*)yIF)->integers[Row[i].index].bitcnt-1; //Goto binary with coef 0 in order to retrieve original Coefficient
					Info.Coef=Row[i].value;
					//NUMERICAL EPS?
					if(abs(Info.Coef.asDouble())>LargestCoef) LargestCoef=abs(Info.Coef.asDouble());


					if(Info.XStar>0 && Info.XStar<Info.Bound.asDouble() && Info.Coef!=0 && std::find(NStar.begin(), NStar.end(),Info.Coef.asDouble())==NStar.end()){
						NStar.push_back(abs(Info.Coef.asDouble()));

					}
					if (Info.XStar>=Info.Bound.asDouble()/(double)2) USet.push_back(Info);
					else TSet.push_back(Info);
				}
				else{	//Interpret this Int as Bool
				            for(int k=0;k<((yInterface*)yIF)->integers[Row[i].index].bitcnt;k++){
					        if(i+k>=Row.size()||(Row[i].index+k>=nVars()) || (Row[i+k].index!=Row[i].index+k) ||
						   ((yInterface*)yIF)->integers[Row[i].index].pt2leader!=
						   ((yInterface*)yIF)->integers[Row[i+k].index].pt2leader){
						  //cerr << "Warning: Binarized Var of Integer appear not grouped" << endl;                  
						  //break;                                                                                   
						  return false;
					        }
						//if(i+k>=Row.size() || ((yInterface*)yIF)->integers[Row[i].index].pt2leader!=((yInterface*)yIF)->integers[Row[i+k].index].pt2leader) break;
						if(1||Row[i+k].value.asDouble()>0){
							OrgInfo NewInfo;
							NewInfo.Bound = 1;
							NewInfo.Coef=Row[i+k].value.asDouble();
							NewInfo.SeeAsBool=true;
							NewInfo.StartIndex=i+k;
							NewInfo.XStar=sol[Row[i].index+k].asDouble();
							if(abs(NewInfo.Coef.asDouble())>LargestCoef) LargestCoef=abs(NewInfo.Coef.asDouble());
							if(NewInfo.XStar>0 && NewInfo.XStar<NewInfo.Bound.asDouble() && NewInfo.Coef!=0 && std::find(NStar.begin(), NStar.end(),NewInfo.Coef.asDouble())==NStar.end()){
								NStar.push_back(abs(NewInfo.Coef.asDouble()));

							}
							if (NewInfo.XStar>=NewInfo.Bound.asDouble()/(double)2) USet.push_back(NewInfo);
							else TSet.push_back(NewInfo);
						}
					}
					i=i+((yInterface*)yIF)->integers[Row[i].index].bitcnt-1; //Goto binary with coef 0 in order to retrieve original Coefficient
				}
			}
		}
			if(LargestCoef>eps)NStar.push_back(1+LargestCoef);
		return true;
	}


	void QBPSolver::CollectRowInfo(CRef Cons,  vector<data::QpNum>& Solution, double &Slack, double &Euclid, const double & eps, bool & OnlyNN){
		Constraint &c = constraintallocator[Cons];
		//ax>=b  =>  ax-s=b  => s=ax-b
		Slack= -c.header.rhs;
		Euclid=0;
		for (int i=0;i<c.size();i++){
			if(lowerBounds[var(c[i])]<0) OnlyNN=false;
			Euclid += c[i].coef*c[i].coef;
			if(sign(c[i]))
				Slack-=c[i].coef*Solution[var(c[i])].asDouble();
			else
				Slack+=c[i].coef*Solution[var(c[i])].asDouble();
		}

		if (abs(Slack)<eps) Slack=0;
		//else if(Slack<0){
		//	cerr <<"Weird Slack: " <<Slack<< endl;
		//}
		Euclid = sqrt(Euclid);
	}

	bool QBPSolver::FindConstraintForAggregation(vector<CRef>& cons,const int& index, double & BestAggrScore, CRef & BestCons, vector<data::QpNum>& Solution, double &Coef, vector<pair<bool,CRef>>& Equations,const double & eps){
		//Check For Variable x_index, whether another row exists with x_index that can be used for aggregation
		//If yes: return the one with best Score
		BestAggrScore=n_infinity;
		assert(cons.size()==Equations.size());
		for(int i=0;i<VarsInConstraints[index].size();i++){
			//Do not use Variable-Bound-Constraint or integer-to-binary-Constraint or the objective function
			if(constraintallocator[VarsInConstraints[index][i].cr].header.isIntBnd
					|| constraintallocator[VarsInConstraints[index][i].cr].header.isVarBnd
					|| VarsInConstraints[index][i].cr==constraints[0])
				continue;
			bool AlreadyUsed=false;
			for(int k=0;k<cons.size();k++){
				// If this row is already part of the aggregation or if the counterpart of initial equation is already part of aggregation
				if(cons[k]==VarsInConstraints[index][i].cr || (Equations[k].first && Equations[k].second==VarsInConstraints[index][i].cr) ){
					AlreadyUsed=true;
					break;
				}
			}
			if(!AlreadyUsed){
				double Sl=0;
				double Norm=0;

				//Would be nice to reuse older calculations ...
				//KRITISCH
				bool OnlyNonNegative=true;
				CollectRowInfo(VarsInConstraints[index][i].cr,Solution, Sl, Norm,eps,OnlyNonNegative);
				if (Sl<0||  !OnlyNonNegative) continue;

				//Scores change due to "UsedInAgg"; But main part is actuall fix, i.e. Sl, Norm and "dual"
				double Sc=ConstraintScore(constraints[i],Solution.size(),Sl,Norm, eps);
				if(Sc>BestAggrScore && Sc>-1){
					BestAggrScore=Sc;
					BestCons=VarsInConstraints[index][i].cr;
					if(sign(constraintallocator[VarsInConstraints[index][i].cr][VarsInConstraints[index][i].pos]))
					Coef=-constraintallocator[VarsInConstraints[index][i].cr][VarsInConstraints[index][i].pos].coef;
					else Coef=constraintallocator[VarsInConstraints[index][i].cr][VarsInConstraints[index][i].pos].coef;
				}
			}
		}
		if(BestAggrScore!=n_infinity)
			return true;
		else return false;
	}

	bool QBPSolver::CheckIfEquation(const CRef& cr, CRef &CP, const double & eps){
		Constraint &c= constraintallocator[cr];
		//KRITISCH
		//Select one of the first two variables: The one with fewer constraints

		int FirstVar;
		if(c.size()<=1)
			FirstVar=var(c[0]);
		else if (VarsInConstraints[var(c[0])].size()<= VarsInConstraints[var(c[1])].size())
			FirstVar=var(c[0]);
		else FirstVar=var(c[1]);
		for(int i=0;i<VarsInConstraints[FirstVar].size();i++){
			if(VarsInConstraints[FirstVar][i].cr!=cr){
				Constraint &comp =constraintallocator[VarsInConstraints[FirstVar][i].cr];
				if(comp.size() == c.size() &&abs((double)comp.header.rhs +(double)c.header.rhs)<eps && abs(comp[VarsInConstraints[FirstVar][i].pos].coef-c[0].coef)<eps &&  sign(comp[VarsInConstraints[FirstVar][i].pos])!= sign(c[0])){
					// Check if comp is counterpart as they previously might have been an equation
					bool IsCounterpart=true;
					for (int k=0;k<c.size();k++){
						if(abs(c[k].coef-comp[k].coef)>eps ||  sign(c[k]) ==sign(comp[k])){
							IsCounterpart=false;
							break;
						}
					}
					if(IsCounterpart){
						CP = VarsInConstraints[FirstVar][i].cr;
						return true;
					}
				}
			}
		}
		return false;
	}

	bool QBPSolver::CutVecOk(vector<std::pair<unsigned int,double>> vec, double eps, double ratio, int CutSize){
			if (vec.size()>CutSize){
				//cerr << "Cut too long" << endl;
				return false;
			}
			double AMax=abs(vec[0].second);
			double AMin=abs(vec[0].second);
			for(int i=0;i<vec.size();i++){
				//et Max/Min value of all but last entry
				if(abs(vec[i].second)>AMax) AMax=abs(vec[i].second);
				if(abs(vec[i].second)<AMin) AMin=abs(vec[i].second);
				if(abs(vec[i].second)<eps){
					//cerr << "Coefficient too small" << endl;
					return false;
				}
			}
			if(AMax/AMin > ratio){
				//cerr << "Ratio of Cut Coefficients too large" << endl;
				return false;
			}
			return true;
		}

		bool QBPSolver::InputConstraintOk(Constraint& vec, double eps, double ratio, int CutSize){
			if (vec.size()>CutSize){
				//cerr << "Constraint too long " << vec.size() << endl;
				return false;
			}
			double AMax=abs(vec[0].coef);
			double AMin=abs(vec[0].coef);
			for(int i=0;i<vec.size();i++){
				//et Max/Min value of all but last entry
				if(abs(vec[i].coef)>AMax) AMax=abs(vec[i].coef);
				if(abs(vec[i].coef)<AMin) AMin=abs(vec[i].coef);
				if(abs(vec[i].coef)<eps){
					//cerr << "Coefficient too small " << abs(vec[i].coef) << endl;
					return false;
				}
			}
			if(AMax/AMin > ratio){
				//cerr << "Ratio of Constraint Coefficients too large " << AMax/AMin <<  endl;
				return false;
			}
			//cerr << "Fine " <<  AMax/AMin << " " << vec.size()  << endl;
			return true;
		}

		bool QBPSolver::InputConstraintOk(vector<data::IndexedElement> vec, double eps, double ratio, int CutSize){
				if (vec.size()>CutSize){
					//cerr << "Constraint too long " << vec.size() << endl;
					return false;
				}
				double AMax=abs(vec[0].value.asDouble());
				double AMin=abs(vec[0].value.asDouble());
				for(int i=0;i<vec.size();i++){
					//et Max/Min value of all but last entry
					if(abs(vec[i].value.asDouble())>AMax) AMax=abs(vec[i].value.asDouble());
					if(abs(vec[i].value.asDouble())<AMin) AMin=abs(vec[i].value.asDouble());
					if(abs(vec[i].value.asDouble())<eps){
						//cerr << "Coefficient too small " << abs(vec[i].coef) << endl;
						return false;
					}
				}
				if(AMax/AMin > ratio){
					//cerr << "Ratio of Constraint Coefficients too large " << AMax/AMin <<  endl;
					return false;
				}
				//cerr << "Fine " <<  AMax/AMin << " " << vec.size()  << endl;
				return true;
			}

		bool QBPSolver::AggrValOk(vector<double> vec, double eps, double ratio){
			double AMax=abs(vec[0]);
			double AMin=abs(vec[0]);
			for(int i=0;i<vec.size();i++){
				//et Max/Min value of all but last entry
				if(abs(vec[i])>AMax  && i!=vec.size()-1) AMax=abs(vec[i]);
				if(abs(vec[i])<AMin && i!=vec.size()-1) AMin=abs(vec[i]);
				if(abs(vec[i])<eps){
					//cerr << "Aggregation Value too small" << endl;
					return false;
				}
			}
			if(AMax/abs(vec[vec.size()-1]) > ratio){
				//cerr << "Ration of Aggregation Values too large" << endl;
				return false;
			}
			if(abs(vec[vec.size()-1])/AMin > ratio){
				//cerr << "Ration of Aggregation Values too large" << endl;
				return false;
			}
			return true;
		}

	int QBPSolver::GenerateCMIRCut( extSol::QpExternSolver& externSolver, vector< vector< data::IndexedElement > > &listOfCutsLhs,
			       vector< data::QpNum > &listOfCutsRhs, vector<int> &listOfCutsVars,
						  int treedepth, int currentBlock, bool &global_valid, std::vector<unsigned int> &candidates, int cuttype, int*types, int8_t* assigns, unsigned int initime, int* solu, int *fixs, int *blcks, int orgN){
	  //return 0;
	  //return listOfCutsLhs.size();

			bool Print=false;

			if(Print) cerr <<"StartSCMIR " << constraints.size()<< endl;
	 		const double eps = 1e-7;
			// Constraint Orientation of Constraints in constraintallocator: \sum ax >= b

			if (externSolver.getSolutionStatus() != extSol::QpExternSolver::OPTIMAL){
				if(getShowWarning()) cerr << "Warning: invalid solution status -> no c-MIR-cuts" << endl;
				//cin.get();
				return listOfCutsRhs.size();
			}
			const unsigned int n = externSolver.getVariableCount();
			const unsigned int m = constraints.size();
			int FirstSlack = max(n,m);
			if (FirstSlack > 100000) return listOfCutsLhs.size();
			 vector<data::QpNum> objVals(n);
			externSolver.getValues(objVals);

			int MAXCUTS=100;
			int MAXCONTS=20;
			int MAXFAILS = 150;
			int MAXAGGR = 6;
			int Fails=0;
			int MAXTESTDELTA=10;
			double MINFRAC=0.05;
			double MAXFRAC=0.95;
			/*static*/ vector< pair< vector< pair<unsigned int, double> >, double > > cuts(MAXCUTS, make_pair(vector< pair<unsigned int, double> >(), 0));
			cuts.clear();


			/*Vector Containing the Aggregated Constraints (CRef)*/
			vector<CRef> VectorOfConstraints;

			/*Vector holding information on aggregated Constraints ax>=b. Are they actually constraints, i.e. ax=b? If yes: .second is Counterpart -ax>=-b*/
			vector<pair<bool,CRef>> ConstraintIsEquation;

			//Vector containing the initial Constraint-Scores used to select starting constraints
			vector<pair<CRef,double>> AggrScore;

			//Vector containing the Slack Value of each Row
			vector<pair<CRef,double>> SlacksOfRows;

			//Vector containing multiplier for each row used for aggregation (Needed for stability)
			vector<double> VectorAggregationValues;

			//Select rows eligible as starting row
			for (int i=0;i<constraints.size(); i++){
				double Sl=0;
				double Norm=0;
				bool OnlyNonNegative=true;

				// Get current slack (s=ax-b) and euclidean norm
				//KRITISCH! Wird für eine Zeile öfter ausgeführt
				CollectRowInfo(constraints[i],objVals, Sl, Norm,eps,OnlyNonNegative);
				if(Print && !OnlyNonNegative) cerr << "Nonnegative Variable-> CMIR-Danger" << endl;

				//If some slack is negative then something is weird: This contraints[i] was probably not part of the relaxation
				//Hence this row is not used for the Aggregation
				//Also: If some variable in this row has negative lower bound we won't use it
				if (Sl<0 || !OnlyNonNegative) continue;

				// Slack is stored for future use; Currently unused
				SlacksOfRows.push_back(std::make_pair(constraints[i],Sl));

				//Calculate Score of this row according to AGGRSCORE (see Wolter, page 23)
				double Sc=ConstraintScore(constraints[i],n,Sl,Norm,eps);
				if(Print) cerr << i << " " << Sc << " " << constraintallocator[constraints[i]].header.isIntBnd << " " <<  constraintallocator[constraints[i]].header.isVarBnd <<  endl;
				if(Sc>-1 && ContCountOK(constraintallocator[constraints[i]],0,MAXCONTS,objVals,eps))
				AggrScore.push_back(make_pair(constraints[i] ,Sc));
				//cerr<<"Slack " << i << " " <<SlacksOfRows[i].first << "  " <<SlacksOfRows[i].second;
			}
			//cin.get();
			if (Print) cerr << "Number of Allowed Constraints: " << AggrScore.size() << endl;
			//Sort Constraints; largest score first
			std::sort(AggrScore.begin(), AggrScore.end(), [](pair<CRef,double> e1,pair<CRef,double> e2){return e1.second > e2.second;});


			for (int StartingConstraintIndex=0; StartingConstraintIndex <AggrScore.size();StartingConstraintIndex++) {
				//cerr << "FAILS " << Fails <<endl;
				bool FoundOne=false;
				if(Print) cerr <<"CMIR " << StartingConstraintIndex<<"/"<<MAXCUTS << " Score " << AggrScore[StartingConstraintIndex].second << endl;
				VectorOfConstraints.clear();
				ConstraintIsEquation.clear();
				VectorAggregationValues.clear();
				//If too many CONSECUTIVE Fails!
				if(Fails>=MAXFAILS+max(MAXFAILS-2*StartingConstraintIndex,0)) break;
				if(cuts.size() >MAXCUTS) break;
				//Select and Check Starting row
				CRef StartingConstraint = AggrScore[StartingConstraintIndex].first;

				if (StartingConstraint==constraints[0]){
					MAXCUTS++;
					//Fails++; //Not really a fail, is it?
					continue;
				}

				vector< data::IndexedElement > AggregatedRow;
				AggregatedRow.clear();

				//Slacks of the aggregated rows
				vector<double> SlackOfConstraints;
				double AggregatedRHS=0.0;

				while(VectorOfConstraints.size()<MAXAGGR)
				{
					//Variable used for aggregation
					int FinalVar=-1;
					vector< data::IndexedElement > XMK;
					XMK.clear();
					CRef NextConstraint;
					double AggregationValue=1;

					if(VectorOfConstraints.size()==0){
						//Add Starting Constraint to XMK
						NextConstraint=StartingConstraint;
						AggregationValue=1;
						VectorAggregationValues.push_back(AggregationValue);
					}
					else
					{
						//Find Row to Aggregate

						// See Wolter, Algorithm 3.2., line 1-19
						if(AggregatedRow.size()==0){
							Fails++;
							break;
						}
						double GlobalBestAggrScore=n_infinity;
						double BestBoundDist=0;

						double CoefToAggr;
						ContiInfo FinalVarInfo;
						for (int i=0;i<AggregatedRow.size();i++){
							if (AggregatedRow[i].index<FirstSlack && type[AggregatedRow[i].index]==CONTINUOUS){
								// Check each available real variables for aggregation
								ContiInfo VarI;
								VarI.index =AggregatedRow[i].index;
								VarI.LocationInAggr=i;
								VarI.AggrCoef=AggregatedRow[i].value.asDouble();
								VarI.lbStar=lowerBounds[VarI.index];
								VarI.ubStar=upperBounds[VarI.index];

								double xstar=0;
								if(VariableBound[VarI.index].ActiveLB || VariableBound[VarI.index].ActiveUB){
									for (int k=VariableBound[VarI.index].VarIndexFirst;k<= VariableBound[VarI.index].VarIndexLast;k++){
										xstar+=objVals[k].asDouble()*pow(2,VariableBound[VarI.index].VarIndexLast-k);
									}
								}
								if(VariableBound[VarI.index].ActiveLB&& VariableBound[VarI.index].lb*xstar>VarI.lbStar)
									VarI.lbStar=VariableBound[VarI.index].lb*xstar;
								if(VariableBound[VarI.index].ActiveUB&& VariableBound[VarI.index].ub*xstar<VarI.ubStar)
									VarI.ubStar=VariableBound[VarI.index].ub*xstar;
								VarI.d=min(objVals[VarI.index].asDouble()-VarI.lbStar,VarI.ubStar-objVals[VarI.index].asDouble());
								if(abs(VarI.d)<eps) VarI.d=0;
								double BestScore;
								CRef BestCons;
								double Coef;
								bool OthersExist=FindConstraintForAggregation(VectorOfConstraints, VarI.index, BestScore, BestCons, objVals,Coef,ConstraintIsEquation,eps);
								if(VarI.d>eps && OthersExist && ((VarI.d>=BestBoundDist&&BestScore>GlobalBestAggrScore) ||  VarI.d>BestBoundDist)){//>eps
									//Unclear: Wolter, Algo 3.2, line 14 vs. middle page 22: select variable with greatest bound distance VarI.d AND THEN the constraint with greatest AGGRSCORE
									GlobalBestAggrScore=BestScore;
									BestBoundDist=VarI.d;
									FinalVar=VarI.index;
									NextConstraint=BestCons;
									FinalVarInfo=VarI;
									CoefToAggr=Coef;
								}
							}
						}
						if(FinalVar>=0 && FinalVar<nVars()){
							//cerr << "Final Var: " << FinalVar << "and info "<< FinalVarInfo.index << endl;
							assert(FinalVarInfo.index==FinalVar);
							AggregationValue=-FinalVarInfo.AggrCoef/CoefToAggr;
							VectorAggregationValues.push_back(AggregationValue);
						}
						else
						{
							//Aggregation Failed; No Further Constraint available
							Fails++;
							break;
						}
						//Find Aggregatable Row
							//NextConstraint=
					}
					if(!AggrValOk(VectorAggregationValues,eps,10000)){
						Fails++;
						break;
					}
					Constraint &con= constraintallocator[NextConstraint];
					CRef Counterpart = NextConstraint;
					bool IsEquation=CheckIfEquation(NextConstraint, Counterpart, eps);
					ConstraintIsEquation.push_back(std::make_pair(IsEquation,Counterpart));
					if (IsEquation ){
						assert(Counterpart!=NextConstraint);
						if (Print) cerr << "Constraint " << NextConstraint << " has Counterpart " << Counterpart << endl;
					}
					else if(Print) cerr << "Constraint " << NextConstraint << " has NO Counterpart " << endl;
					//--------------Aggregate Rows
					if(Print){
						cerr << "Before Aggregation: " << endl;
						for (int i=0; i<AggregatedRow.size();i++){
							if(AggregatedRow[i].index<FirstSlack)
							 cerr << " +" <<AggregatedRow[i].value.asDouble()<<"x_"<<AggregatedRow[i].index<<"("<<objVals[AggregatedRow[i].index].asDouble()<<","<<(type[AggregatedRow[i].index]==CONTINUOUS?"c,":"b,") <<" ["<<lowerBounds[AggregatedRow[i].index]<< ","<<upperBounds[AggregatedRow[i].index]<<"] - "<<VariableBound[AggregatedRow[i].index].VarIndexFirst << "[" << VariableBound[AggregatedRow[i].index].lb << "," << VariableBound[AggregatedRow[i].index].ub
									<<"])";
							else  cerr <<" +"<<AggregatedRow[i].value.asDouble()<<"x_"<<AggregatedRow[i].index<<"(Slack)";
						}
						cerr << "= "<< AggregatedRHS << endl;
					}
					if(Print) cerr << "Aggregate with Faktor " << AggregationValue << " due to variable " << FinalVar<< endl;
					/*if (FinalVar < 0 || FinalVar >= nVars()) {
						Fails++;
						break;
					}*/
					for (int i=0; i<con.size();i++){
						if(Print)cerr << " +" <<(sign(con[i])?data::QpNum(-con[i].coef):data::QpNum(con[i].coef)) << "x_" << var(con[i]) ;
						AggregatedRow.push_back(data::IndexedElement( var(con[i]), (sign(con[i])?data::QpNum(-AggregationValue*con[i].coef):data::QpNum(AggregationValue*con[i].coef)) ));
					}
					if(Print && !IsEquation)cerr << "-x_"<< FirstSlack+VectorOfConstraints.size() << " = " <<con.header.rhs << endl;

					//If it is equation no slack is needed
					if (!IsEquation) AggregatedRow.push_back(data::IndexedElement( FirstSlack+VectorOfConstraints.size(),data::QpNum(-AggregationValue*1)));

					//Sort AggregatedRow, aggregate multiple variable entries and delete zeros
					ClearDoubleEntries(AggregatedRow,eps);


					if(!ContCountOK(AggregatedRow,0,MAXCONTS,objVals,eps,FirstSlack)){
					        if (Print) cerr << "Too many reals after Aggregation" << endl;
						assert(VectorOfConstraints.size()>0);
						Fails++;
						break;
					}
					VectorOfConstraints.push_back(NextConstraint);
					AggregatedRHS+=AggregationValue*con.header.rhs;

					for(int i=0;i<AggregatedRow.size();i++)
						XMK.push_back(AggregatedRow[i]);
					double Slack;

					if(IsEquation){
						Slack=0;
					}
					else{
						Slack=-con.header.rhs;
						for (int i=0;i<con.size();i++){
							if(sign(con[i]))
								Slack-=objVals[var(con[i])].asDouble()*con[i].coef;
							else
								Slack+=objVals[var(con[i])].asDouble()*con[i].coef;
						}
					}
					if(Print) cerr <<"Slack: " << Slack << endl;
					if(Slack<-eps){
						Fails++;
						break;
					}
					con.header.usedinAgg++;
					SlackOfConstraints.push_back(Slack);
					// Constraints are of Type ax>=b
					// Introduce Slack s>=0 such that ax-s=b
					// Slack Variable Corresponding to Row VectorOfConstraints[i] has Index constraints.size()+i


					data::QpNum XMKRHS=AggregatedRHS;


					double s_star=0;
					vector<std::pair<int,string>> SubType;
					if(Print){
						cerr << "Aggregated: "<< endl;
						for (int i=0; i<AggregatedRow.size();i++){
							if(AggregatedRow[i].index<FirstSlack)
							 cerr << " +" <<AggregatedRow[i].value.asDouble()<<"x_"<<AggregatedRow[i].index<<"("<<objVals[AggregatedRow[i].index].asDouble()<<","<<(type[AggregatedRow[i].index]==CONTINUOUS?"c,":"b,") <<" ["<<lowerBounds[AggregatedRow[i].index]<< ","<<upperBounds[AggregatedRow[i].index]<<"] - "<<VariableBound[AggregatedRow[i].index].VarIndexFirst << "[" << VariableBound[AggregatedRow[i].index].lb << "," << VariableBound[AggregatedRow[i].index].ub
									<<"])";
							else  cerr <<" +"<<AggregatedRow[i].value.asDouble()<<"x_"<<AggregatedRow[i].index<<"(Slack)";
						}
						cerr << "= "<< AggregatedRHS << endl;
					}

					// BOUND SUBSTITUTION
					if (Print) cerr <<"Bound Substitution" << endl;

					for (int i =0;i<AggregatedRow.size();i++){
						int ind=AggregatedRow[i].index;
						if(ind>=FirstSlack){
							//Slack stay the way they are, or actually they are substituted with LB, i.e. nothing changes;

							if(AggregatedRow[i].value<0)
								s_star-= AggregatedRow[i].value.asDouble()*SlackOfConstraints[ind-FirstSlack];
							else{
								assert(XMK[i].value>=0);
								XMK[i].value=0;
							}
							continue;
						}

						if(type[ind]==CONTINUOUS){
							// See Wolters, Algo 3.3, line 4-7
							double lbjstar=0;
							double ubjstar=0;
							bool VarLB=false;
							bool VarUB=false;
							double xstar=0;
							if (VariableBound[ind].ActiveLB){
								//Calculate l*x*
								if(VariableBound[ind].VarIndexFirst!= VariableBound[ind].VarIndexLast){
									//If VarUBVar is Integer Variable
									for (int k=VariableBound[ind].VarIndexFirst;k<= VariableBound[ind].VarIndexLast;k++){
										xstar+=objVals[k].asDouble()*pow(2,VariableBound[ind].VarIndexLast-k);
									}
								}
								else{
									xstar=objVals[VariableBound[ind].VarIndexFirst].asDouble();
								}
								double lx=VariableBound[ind].lb*xstar;

								if(lowerBounds[ind]>lx){
									lbjstar=lowerBounds[ind];
									VarLB=false;
								}
								else{
									lbjstar=lx;
									VarLB=true;
								}
							}
							else{
								lbjstar=lowerBounds[ind];
								VarLB=false;
							}

							if (VariableBound[ind].ActiveUB){
								//Calculate u*x*
								if(!VariableBound[ind].ActiveLB){	//If xstar was not calculated previously
									if(VariableBound[ind].VarIndexFirst!= VariableBound[ind].VarIndexLast){
										//If VarUBVar is Integer Variable
										for (int k=VariableBound[ind].VarIndexFirst;k<= VariableBound[ind].VarIndexLast;k++){
											xstar+=objVals[k].asDouble()*pow(2,VariableBound[ind].VarIndexLast-k);
										}
									}
									else{
										xstar=objVals[VariableBound[ind].VarIndexFirst].asDouble();
									}
								}
								double ux=VariableBound[ind].ub*xstar;

								if(upperBounds[ind]<ux){
									ubjstar=upperBounds[ind];
									VarUB=false;
								}
								else{
									ubjstar=ux;
									VarUB=true;
								}
							}
							else{
								ubjstar=upperBounds[ind];
								VarUB=false;
							}

							if(objVals[ind].asDouble()-lbjstar<=ubjstar-objVals[ind].asDouble()){
								// See Wolters, Algo 3.3, line 8-13
								//substitute x_ind=lbj + y_j
								if(VarLB ){	//Substitute with variable Bound
									SubType.push_back(make_pair(ind,"VarLB"));
									for (int k=VariableBound[ind].VarIndexFirst;k<= VariableBound[ind].VarIndexLast;k++){
										XMK.push_back(data::IndexedElement(k,data::QpNum(pow(2,VariableBound[ind].VarIndexLast-k)*AggregatedRow[i].value.asDouble()*VariableBound[ind].lb)));
									}
								}
								else{	//Substitute with normal Bound
									SubType.push_back(make_pair(ind,"LB"));
									XMKRHS -= AggregatedRow[i].value*lowerBounds[ind];
								}
								if(XMK[i].value<0) s_star-= XMK[i].value.asDouble()*(objVals[ind].asDouble()-lbjstar);
							//Store Information on Substitution; Dont Forget Solution in objVals!!!!!!!!
							}
							else{
								// See Wolters, Algo 3.3, line 14-19
								if(VarUB){	//Substitute with variable Bound
									SubType.push_back(make_pair(ind,"VarUB"));
									XMK[i].value=-AggregatedRow[i].value.asDouble();
									for (int k=VariableBound[ind].VarIndexFirst;k<= VariableBound[ind].VarIndexLast;k++){
										XMK.push_back(data::IndexedElement(k,data::QpNum(pow(2,VariableBound[ind].VarIndexLast-k)*AggregatedRow[i].value.asDouble()*VariableBound[ind].ub)));
									}
								}
								else{	//Substitute with normal Bound
									SubType.push_back(make_pair(ind,"UB"));
									XMK[i].value= -AggregatedRow[i].value.asDouble();
									XMKRHS -= AggregatedRow[i].value.asDouble()*upperBounds[ind];
								}
								if(XMK[i].value<0) s_star-= XMK[i].value.asDouble()*(ubjstar-objVals[ind].asDouble());
							}
							if(XMK[i].value>=0) XMK[i].value=0;
						}
					}
					if(Print){
						cerr << "After Bound Substitution: "<< endl;
						int co=0;
						for (int i=0; i<XMK.size();i++){
							if(XMK[i].index<FirstSlack &&type[XMK[i].index]==CONTINUOUS){
								assert (SubType[co].first ==XMK[i].index);
								cerr << " "<<XMK[i].value.asDouble()<<"x_"<<XMK[i].index<<"("<<SubType[co].second<<")";
								co++;
							}
							else{
								cerr << " " <<XMK[i].value.asDouble()<<"x_"<<XMK[i].index;
							}
						}

						cerr << "= "<< XMKRHS << endl;
					}

					ClearDoubleEntries(XMK, eps);
					if (Print){
						cerr <<"Cleared XMK:" << endl;
						for (int i=0; i<XMK.size();i++){
							if(XMK[i].index<FirstSlack &&type[XMK[i].index]==CONTINUOUS){
								cerr << " "<<XMK[i].value.asDouble()<<"x_"<<XMK[i].index;
							}
							else{
								cerr << " " <<XMK[i].value.asDouble()<<"x_"<<XMK[i].index;
							}
						}
						cerr << "<= "<< XMKRHS << endl;
					}
					// XMK is Mixed Knapsack Set

					if(abs(XMKRHS.asDouble()-floor(XMKRHS.asDouble()+.5))<eps) XMKRHS=floor(XMKRHS.asDouble()+.5);
					if(abs(s_star-floor(s_star+.5))<eps) s_star=floor(s_star+.5);
					if (Print) cerr << "SStar: " << s_star << endl;
					//Now see Wolters Algo 3.4
					vector<OrgInfo> U,T;
					vector<double> NStar;
					NStar.clear();
					U.clear();
					T.clear();

					if(!CollectIntegerInfos(XMK,U,T,objVals,NStar,FirstSlack,eps)){
						Fails++;
						break;
					}
					if(Print){
						cerr <<"Collected Integer Infos" << endl;
						cerr <<"U:" << endl;
						for (int i=0;i<U.size();i++){
							cerr <<i<< "StartIndex: " << U[i].StartIndex << "; Coef: " << U[i].Coef.asDouble() << "; Bound: " << U[i].Bound.asDouble() << "; XStar: " << U[i].XStar<<endl;
						}
						cerr <<"T:" << endl;
						for (int i=0;i<T.size();i++){
							cerr << i<<" StartIndex: " << T[i].StartIndex << "; Coef: " << T[i].Coef.asDouble() << "; Bound: " << T[i].Bound.asDouble() << "; XStar: " << T[i].XStar<<endl;
						}
						cerr <<"NStar:" << endl;
						for (int i=0;i<NStar.size();i++){
							cerr << i<< " " << NStar[i] <<endl;
						}
					}
					int SizeT=T.size();
					bool DeltaFound=false;
					double VBest=n_infinity;
					double DeltaBest=1;

					int ell=0;

					//Attention: Sorted in nondecreasing order
					//hence go backwards through T, as pop_back is used
					std::sort(T.begin(), T.end(), [](OrgInfo e1, OrgInfo e2){return (e1.XStar-e1.Bound.asDouble()/(double)2) > (e2.XStar-e2.Bound.asDouble()/(double)2);});

					do{
						if(Print)cerr << "ELL: " << ell << endl;
						if (ell>0 && T.size()>0&& T[T.size()-ell].XStar>eps){
							U.push_back(T[T.size()-ell]);
							T[T.size()-ell]=T[T.size()-1];
							T.pop_back();
							SizeT--;
						}
						else if(ell>0){
							ell++; continue;
						}
						for (int del=0;del<NStar.size() && del <MAXTESTDELTA;del++){
							double delta=NStar[del];
							if(Print) cerr <<del << " Delta " << delta << endl;
							double beta=XMKRHS.asDouble();
							if(Print) cerr <<"Beta_Start: " << beta<<endl;

							for(int j=0;j<U.size();j++){
								beta-=U[j].Coef.asDouble()*U[j].Bound.asDouble();
								if(Print) cerr <<"Beta_Change: " << -U[j].Coef.asDouble()*U[j].Bound.asDouble()<<endl;
							}
							beta=(double) beta/delta;
							if(Print) cerr <<"Beta_Final: " << beta<<endl;

							double fbeta;
							if(abs(beta-floor(beta+.5))<eps) fbeta=0;
							else fbeta=beta-floor(beta);
							if(fbeta<MINFRAC  || fbeta>MAXFRAC) continue;
							if(Print) cerr << "f_beta " << fbeta << endl;
							if(Print) cerr << "Continue? " << abs(beta-floor(beta+.5)) << endl;
							if(abs(beta-floor(beta+.5))<eps){
								continue;
							}
							else{
								if(Print) cerr <<"Found Something!" << endl;
								DeltaFound=true;
								//+Eps?
								double v= -floor(beta+eps)-s_star/(double)(delta*(1-fbeta));
								if(Print) cerr << "v_Start: " << v << endl;
								for (int tind=0;tind<T.size();tind++){
									v+=F(fbeta,T[tind].Coef.asDouble()/(double)delta)*T[tind].XStar;
									if(Print) cerr << "v(T) changed by x_" << T[tind].StartIndex  << " : "<< F(fbeta,T[tind].Coef.asDouble()/(double)delta)*T[tind].XStar << endl;
								}
								for (int uind=0;uind<U.size();uind++){
									v+=F(fbeta,-U[uind].Coef.asDouble()/(double)delta)*(U[uind].Bound.asDouble()-U[uind].XStar);
									if(Print) cerr << "v(U) changed by x_" << U[uind].StartIndex  << " : "<< F(fbeta,-U[uind].Coef.asDouble()/(double)delta)*(U[uind].Bound.asDouble()-U[uind].XStar)<< endl;

								}
								if (Print) cerr << "Compare v=" << v << " vs. vbest=" << VBest << endl;
								if(v>VBest){
									VBest=v;
									DeltaBest=delta;
								}
							}

						}
						ell++;

					} while(DeltaFound!=true && ell<SizeT && T.size()>0);

					if(DeltaFound){
						//Improve Violation by modifying delta; Wolter Algo 3.4, line 18ff
						double Delta_Bar=DeltaBest;
						for(int div=2;div<=8;div=div*2){
							if (Print)cerr <<"div " << div << endl;
							double delta=Delta_Bar/(double)div;
							double beta=XMKRHS.asDouble();
							for(int j=0;j<U.size();j++){
								beta-=U[j].Coef.asDouble()*U[j].Bound.asDouble();
							}
							beta=(double) beta/delta;
							if(Print) cerr <<"New Beta " << beta << endl;
							double fbeta;
							if(abs(beta-floor(beta+.5))<eps) fbeta=0;
							else fbeta=beta-floor(beta);
							if(fbeta<MINFRAC || fbeta>MAXFRAC) continue;

							if(Print) cerr <<"New FBeta " << fbeta << endl;
							double v= -floor(beta+eps)-s_star/(double)(delta*(1-fbeta));
							if(Print) cerr << "v_Start: " << v << endl;

							for (int tind=0;tind<T.size();tind++){
								v+=F(fbeta,T[tind].Coef.asDouble()/(double)delta)*T[tind].XStar;
								if(Print) cerr << "v(T) changed by x_" << T[tind].StartIndex  << " : "<< F(fbeta,T[tind].Coef.asDouble()/(double)delta)*T[tind].XStar << endl;

							}
							for (int uind=0;uind<U.size();uind++){
								v+=F(fbeta,-U[uind].Coef.asDouble()/(double)delta)*(U[uind].Bound.asDouble()-U[uind].XStar);
								if(Print) cerr << "v(U) changed by x_" << U[uind].StartIndex  << " : "<< F(fbeta,-U[uind].Coef.asDouble()/(double)delta)*(U[uind].Bound.asDouble()-U[uind].XStar)<< endl;

							}
							if(Print) cerr << "Did v improve? " << v << " . Best was " << VBest << endl;
							if(v>VBest){
								VBest=v;
								DeltaBest=delta;
							}
						}



						if(VBest>eps){	//If Violated Constraint found?
							//Improve Violation by complementing additional variables
							//see Wolters, Algo 3.4., line 24-32

							if(Print) cerr <<"IMPROVE VIOLATION BY COMPLEMENTING" << endl;
							std::sort(T.begin(), T.end(), [](OrgInfo e1, OrgInfo e2){return (e1.XStar-e1.Bound.asDouble()/(double)2) < (e2.XStar-e2.Bound.asDouble()/(double)2);});
							for (int ell1=T.size()-1; ell1>=0;ell1--){
								assert (ell1<T.size() && ell1>=0);
								if(Print)cerr << "ELL1: " << ell1 << endl;
								if(T[ell1].XStar>0){
									// Implizit! U.push_back(T[ell1]);
									//T.pop_back();
									double beta=XMKRHS.asDouble();
									if(Print) cerr <<"Beta_Start: " << beta<<endl;

									for(int j=0;j<U.size();j++){
										beta-=U[j].Coef.asDouble()*U[j].Bound.asDouble();
										if(Print) cerr <<"Beta_Change: " << -U[j].Coef.asDouble()*U[j].Bound.asDouble()<<endl;
									}
									for(int j=T.size()-1;j>=ell1;j--){
										beta-=T[j].Coef.asDouble()*T[j].Bound.asDouble();
										if(Print) cerr <<"Beta_Change (prevT): " << -T[j].Coef.asDouble()*T[j].Bound.asDouble()<<endl;
									}
									beta=(double) beta/(double) DeltaBest;
									if(Print) cerr <<"Beta_Final: " << beta<<endl;

									double fbeta;
									if(abs(beta-floor(beta+.5))<eps) fbeta=0;
									else fbeta=beta-floor(beta);
									if(fbeta<MINFRAC  || fbeta>MAXFRAC) continue;
									if(Print) cerr << "f_beta " << fbeta << endl;
									if(Print) cerr << "Continue? " << abs(beta-floor(beta+.5)) << endl;
									if(abs(beta-floor(beta+.5))<eps){
										continue;
									}
									else{
										if(Print) cerr <<"Found Something!" << endl;
										DeltaFound=true;
										//+Eps?
										double v= -floor(beta+eps)-s_star/(double)(DeltaBest*(1-fbeta));
										if(Print) cerr << "v_Start: " << v << endl;
										for (int tind=0;tind<ell1;tind++){
											v+=F(fbeta,T[tind].Coef.asDouble()/(double)DeltaBest)*T[tind].XStar;
											if(Print) cerr << "v(T) changed by x_" << T[tind].StartIndex  << " : "<< F(fbeta,T[tind].Coef.asDouble()/(double)DeltaBest)*T[tind].XStar << endl;
										}
										for (int tind=ell1;tind<T.size();tind++){
											v+=F(fbeta,-T[tind].Coef.asDouble()/(double)DeltaBest)*(T[tind].Bound.asDouble()-T[tind].XStar);
											if(Print) cerr << "v(U) (prevT) changed by x_" << T[tind].StartIndex  << " : "<< F(fbeta,-T[tind].Coef.asDouble()/(double)DeltaBest)*(T[tind].Bound.asDouble()-T[tind].XStar)<< endl;
										}

										for (int uind=0;uind<U.size();uind++){
											v+=F(fbeta,-U[uind].Coef.asDouble()/(double)DeltaBest)*(U[uind].Bound.asDouble()-U[uind].XStar);
											if(Print) cerr << "v(U) changed by x_" << U[uind].StartIndex  << " : "<< F(fbeta,-U[uind].Coef.asDouble()/(double)DeltaBest)*(U[uind].Bound.asDouble()-U[uind].XStar)<< endl;

										}
										if (Print) cerr << "Compare v=" << v << " vs. vbest=" << VBest << endl;
										if(v>VBest){
											if (Print) cerr <<"IMPROVED!!! " << v << " vs " << VBest << endl;

											VBest=v;
											for(int move=T.size()-1;move>=ell1;move--){
												U.push_back(T[move]);
												T.pop_back();
											}

										}

									}
								}
							}


							if(Print) cerr << "Collect Cut" << endl;


							vector<pair<unsigned int, double> > cutvec;
							double beta=XMKRHS.asDouble();
							for(int j=0;j<U.size();j++){
								beta-=U[j].Coef.asDouble()*U[j].Bound.asDouble();
							}
							beta=(double) beta/DeltaBest;
							if(Print) cerr <<"Final Beta " << beta << endl;
							double fbeta;
							if(abs(beta-floor(beta+.5))<eps) fbeta=0;
							else fbeta=beta-floor(beta);
							assert(fbeta>=MINFRAC &&  fbeta<=MAXFRAC);

							if(Print) cerr <<"Final FBeta " << fbeta << endl;
							double final_rhs;
							if(abs(beta-floor(beta+.5))<eps) final_rhs=-floor(beta+.5);
							else final_rhs=-floor(beta);		//minus since größergleich
							if(Print) cerr << "Starting RHS: " <<final_rhs << endl;

							for (int tind=0;tind<T.size();tind++){
								double F_coef=F(fbeta,T[tind].Coef.asDouble()/(double)DeltaBest);
								if(abs(F_coef)>eps){
									if(T[tind].SeeAsBool){
										cutvec.push_back( make_pair( T[tind].StartIndex, -F_coef) );
										if(Print) cerr << "Added Entry(T)/aB: " <<cutvec.back().second <<"x_" << cutvec.back().first << endl;
										//cerr << "Really Used It " << endl;

									}
									else{
										for(int i=0;i<((yInterface*)yIF)->integers[T[tind].StartIndex].bitcnt;i++){
											//Minus, since (F(a/d)(x) and größergleich
											cutvec.push_back( make_pair( T[tind].StartIndex+i, -F_coef*pow(2,((yInterface*)yIF)->integers[T[tind].StartIndex].bitcnt-1-i) ) );
											if(Print) cerr << "Added Entry(T): " <<cutvec.back().second <<"x_" << cutvec.back().first << endl;

										}
									}
								}

							}
							for (int uind=0;uind<U.size();uind++){
								double F_coef=F(fbeta,-U[uind].Coef.asDouble()/(double)DeltaBest);
								if(abs(F_coef)>eps){
									if(U[uind].SeeAsBool){
										cutvec.push_back( make_pair( U[uind].StartIndex, F_coef) );
										if(Print) cerr << "Added Entry(U)/aB: " <<cutvec.back().second <<"x_" << cutvec.back().first << endl;
										//cerr << "Really Used It " << endl;
									}
									else{
										for(int i=0;i<((yInterface*)yIF)->integers[U[uind].StartIndex].bitcnt;i++){
										//Plus, since -F(-a/d)(x) and größergleich
											cutvec.push_back( make_pair( U[uind].StartIndex+i, F_coef*pow(2,((yInterface*)yIF)->integers[U[uind].StartIndex].bitcnt-1-i) ) );
											if(Print) cerr << "Added Entry(U): " <<cutvec.back().second <<"x_" << cutvec.back().first << endl;
										}
									}
									//Plus, since F(-a/d)(U) and brought to rhs and größergleich

									final_rhs +=  F_coef*U[uind].Bound.asDouble();
									if(Print) cerr << "Altered RHS: " << F_coef*U[uind].Bound.asDouble() << endl;


								}
							}


							if (Print)cerr << "Transform s" << endl;
							//turn s back into the actual continuous variable
							double ContiCoef=1/(double)(DeltaBest*(1-fbeta));
							if(abs(ContiCoef-floor(ContiCoef+.5))<eps) ContiCoef=floor(ContiCoef+.5);
							if(Print) cerr << "CoefOfS: " << ContiCoef << endl;
							for (int i=0;i<XMK.size();i++){
								int ind=XMK[i].index;
								if(ind>=FirstSlack){
									//If Slack
									assert (!ConstraintIsEquation[ind-FirstSlack].first);
									if(Print) cerr <<"Slack " << ind << endl;
									Constraint &c = constraintallocator[ VectorOfConstraints[ind-FirstSlack]];
									assert(XMK[i].value==-1 || ind>FirstSlack);
									if (Print)cerr << "SlackCoef " << XMK[i].value.asDouble();
									if(Print) cerr<< "Changed RHS Due to slack "<<  -XMK[i].value.asDouble()*ContiCoef*c.header.rhs<<endl;

									// Bleibt drüben, - von  s=-sum gamma*ybar - von größergleich
									if(abs(-XMK[i].value.asDouble()*ContiCoef*c.header.rhs-floor(-XMK[i].value.asDouble()*ContiCoef*c.header.rhs+.5))<eps)
										final_rhs += floor(-XMK[i].value.asDouble()*ContiCoef*c.header.rhs+.5); 	//Minus? s=-sum gamma*ybar; slack = -rhs+ax ; Mal -1 wegen größergleich; also 3*(-1)=-1
									else final_rhs += -XMK[i].value.asDouble()*ContiCoef*c.header.rhs;
									for (int ii=0;ii<c.size();ii++){
										if(sign(c[ii])){
											//Slack-=objVals[var(c[i])]*c[i].coef;
											//Plus: -1 von s=-sum gamma*ybar; -1 von rüber auf die linke Seite; -1 von sign; -1 von größergleich
											double UsedCoef;
											if(abs(-XMK[i].value.asDouble()*ContiCoef*c[ii].coef-floor(-XMK[i].value.asDouble()*ContiCoef*c[ii].coef+.5))<eps)
												UsedCoef=-floor(-XMK[i].value.asDouble()*ContiCoef*c[ii].coef+.5);
											else
												UsedCoef=+XMK[i].value.asDouble()*ContiCoef*c[ii].coef;
											if(Print) cerr<< "NewCutEntryFromSlack "<<  UsedCoef << "x_"<<var(c[ii])<<endl;
											cutvec.push_back( make_pair( var(c[ii]), UsedCoef ) );
										}
										else{
											//Slack+=objVals[var(c[i])]*c[i].coef;
											double UsedCoef;
										if(abs(-XMK[i].value.asDouble()*ContiCoef*c[ii].coef-floor(-XMK[i].value.asDouble()*ContiCoef*c[ii].coef+.5))<eps)
											UsedCoef=floor(-XMK[i].value.asDouble()*ContiCoef*c[ii].coef+.5);
										else
											UsedCoef=-XMK[i].value.asDouble()*ContiCoef*c[ii].coef;
											//Miuns: -1 von s=-sum gamma*ybar; -1 von rüber auf die linke Seite;  -1 von größergleich
											if(Print) cerr<< "NewCutEntryFromSlack "<<  UsedCoef  << "x_"<<var(c[ii])<<endl;

											cutvec.push_back( make_pair( var(c[ii]), UsedCoef ) );
										}
									}
								}
								else if(type[ind]==CONTINUOUS){
									if(Print) cerr << "ContVar: " << ind <<endl;
									assert (XMK[i].value<0);
									double UsedCoef;
									if(abs(ContiCoef*XMK[i].value.asDouble()-floor(ContiCoef*XMK[i].value.asDouble()+.5))<eps)
										UsedCoef=floor(ContiCoef*XMK[i].value.asDouble()+.5);
									else
										UsedCoef=ContiCoef*XMK[i].value.asDouble();
									std::vector<std::pair<int,string>>::iterator it = std::find_if(SubType.begin(), SubType.end(), [&](pair<int,string> e){
																																	return (e.first==ind);});
									if (it != SubType.end()){
										string Subbed=it->second;
										if(Subbed=="VarLB"){//y_bar=y-lx
												//for (int k=VariableBound[ind].VarIndexFirst;k<= VariableBound[ind].VarIndexLast;k++){
												//	XMK.push_back(data::IndexedElement(k,data::QpNum(pow(2,VariableBound[ind].VarIndexLast-k)*AggregatedRow[i].value*lowerBounds[ind])));
												//}

												//-lx
												for (int k=VariableBound[ind].VarIndexFirst;k<= VariableBound[ind].VarIndexLast;k++){
													// : -lx, move to lhs, größergleich
													cutvec.push_back(make_pair(k, UsedCoef*pow(2,VariableBound[ind].VarIndexLast-k)*VariableBound[ind].lb));
													if(Print) cerr << "NewCutEntry(VARLB): " <<cutvec.back().second <<"x_" << cutvec.back().first << endl;
												}
												cutvec.push_back( make_pair( ind, -UsedCoef ) );
												if(Print) cerr << "NewCutEntry(VARLB;conti): " <<cutvec.back().second <<"x_" << cutvec.back().first << endl;


											}
										else if(Subbed== "LB"){//y_bar=y-l
												//Minus: -1 von s=-sum gamma*ybar; -1 von größergleich;-1 von y_bar=y-l
												final_rhs -= UsedCoef*lowerBounds[ind];
												if(Print) cerr << "Changed rhs(LB): " <<-UsedCoef*lowerBounds[ind] << endl;

												//Minus: -1 von s=-sum gamma*ybar; -1 von größergleich, -1 von auf lhs bringen
												cutvec.push_back( make_pair( ind, -UsedCoef ) );
												if(Print) cerr << "NewCutEntry(LB;conti): " <<cutvec.back().second <<"x_" << cutvec.back().first << endl;

										}
										else if(Subbed== "VarUB"){//y_bar=ux-y
												for (int k=VariableBound[ind].VarIndexFirst;k<= VariableBound[ind].VarIndexLast;k++){
													// : move to lhs, größergleich
													cutvec.push_back(make_pair(k, -UsedCoef*pow(2,VariableBound[ind].VarIndexLast-k)*VariableBound[ind].ub));
													if(Print) cerr << "NewCutEntry(VARUB): " <<cutvec.back().second <<"x_" << cutvec.back().first << endl;

												}
												cutvec.push_back( make_pair( ind, UsedCoef ) );
												if(Print) cerr << "NewCutEntry(VARUB;conti): " <<cutvec.back().second <<"x_" << cutvec.back().first << endl;

										}
										else if(Subbed== "UB"){//y_bar=u-y
												//Plus: -1 von s=-sum gamma*ybar; -1 von größergleich
												final_rhs += UsedCoef*upperBounds[ind];
												if(Print) cerr << "Changed rhs(UB): " <<+UsedCoef*upperBounds[ind] << endl;

												//Plus: -1 von s=-sum gamma*ybar; -1 von größergleich, -1 von auf lhs bringen; -1 von  y_bar=u-y
												cutvec.push_back( make_pair( ind,UsedCoef ) );
												if(Print) cerr << "NewCutEntry(UB;conti): " <<cutvec.back().second <<"x_" << cutvec.back().first << endl;

										}
										else assert(0);
									}
									else assert(0);

								}
							}


							if(abs(final_rhs-floor(final_rhs+.5))<eps) final_rhs=floor(final_rhs+.5);

							if(Print) cerr << " Prev Final "  << final_rhs<< endl ;
							//final_rhs-=1e-10 ;
							final_rhs+=-1e-12 -abs(final_rhs)*1e-12;
							int MaxVars=100;
							if(sqrt(n)>MaxVars) MaxVars=sqrt(n);
							if(!CutVecOk(cutvec,eps,10000,MaxVars)){
								Fails++;
								break;
							}
													if(Print) cerr << " >= " <<  final_rhs << "(" << cutvec.size() << ")" << endl ;
													cuts.push_back(make_pair(cutvec, final_rhs));
							if(Print)cerr <<"CUT FOUND! Number of Aggregated Rows " << VectorOfConstraints.size()<< endl;
							FoundOne=true;
							Fails=0;
							//if(VectorOfConstraints.size()>2) cin.get();
							break;
							//cin.get();
						}
					}
				}
				if(VectorOfConstraints.size()==MAXAGGR && !FoundOne)
					Fails++;
				else if(FoundOne) Fails=0;

			}

	for( unsigned int i = 0; i < cuts.size(); ++i ){
	  vector< pair<unsigned int, double> >& cut = cuts.at( i ).first;
	  vector< data::IndexedElement > addcut;

	  if (1){
	    double max_c = 0.0;
	    for( unsigned int j = 0; j < cut.size(); ++j ){
	      if( fabs( cut.at( j ).second ) > max_c ){
		max_c = fabs( cut.at( j ).second );
	      }
	    }
	    if (cut.size() > 0 && max_c > 1e-9) {
	      //std::cerr << "cut before: ";
	      double scale = 2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0*2.0;
	      double unscale=0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5*0.5;
	      for( unsigned int j = 0; j < cut.size(); ++j ){
		cut[j].second = ceil(cut[j].second / max_c * scale) * unscale;
	      }
	      cuts[i].second = floor(cuts[i].second / max_c * scale) * unscale;
	    }
	  }
	  
	  /*for (int j=0;j<cut.size();j++)
	    addcut.push_back( data::IndexedElement( cut.at( j ).first, data::QpNum( cut.at( j ).second ) ) );
	  double rhs = cuts[i].second;
	  listOfCutsLhs.push_back(addcut);
	  listOfCutsRhs.push_back(rhs);
	  */
	}
	listOfCutsLhs.reserve( cuts.size()+listOfCutsLhs.size() );
	listOfCutsRhs.reserve( cuts.size()+listOfCutsLhs.size() );

			if (info_level > 2) cerr << "CUTS OUT " << endl;
			//cerr <<"Info: Created " << cuts.size() << " c-MIR cuts"<< endl;
			if(1)for( unsigned int i = 0; i < cuts.size(); ++i ){
				vector< data::IndexedElement > addcut;
				addcut.clear();
				int j =0;
				double rrhs=0.0;

				for ( ; j < 1 && i+j < cuts.size();j++) {
					if(Print){
						cerr << "Round " <<i << " "<<j<<endl;
						for (int pp=0;pp<cuts.at( i+j ).first.size();pp++){
							cerr << cuts.at( i+j ).first[pp].second <<"x_"<<cuts.at( i+j ).first[pp].first << endl;
						}
						cerr << ">= " << cuts.at( i+j ).second << endl;
					}
					vector< pair<unsigned int, double> >& cut = cuts.at( i+j ).first;

					bool found=false;
					if (computeCutRatio1(cut) > 1000000.0) continue;

					rrhs = rrhs + cuts[i+j].second;
					addcut.reserve( addcut.size()+cut.size() );
					for( unsigned int jj = 0; jj < cut.size(); ++jj ){
						addcut.push_back( data::IndexedElement( cut.at( jj ).first, data::QpNum( cut.at( jj ).second ) ) );
					}
					std::sort(addcut.begin(), addcut.end(), [](data::IndexedElement e1, data::IndexedElement e2){return e1.index < e2.index;});

					// run through and add
					int ii=0, jj=0;
					for (ii = 0; ii < addcut.size();ii++) {
						jj = ii+1;
						while (jj < addcut.size() && addcut[ii].index == addcut[jj].index) {
							addcut[ii].value += addcut[jj].value;
							addcut[jj].value = 0.0;
							jj++;
						}
						ii = jj-1	;
					}
					// delete 0s
					ii=0, jj=0;
					for (ii = 0; ii < addcut.size();ii++) {
						if (abs((addcut[ii].value.asDouble()))>eps) continue;
						jj = addcut.size()-1;
						while (jj>=ii && abs((addcut[jj].value.asDouble()))<eps) {
							jj--;
							addcut.pop_back();
						}
						if (ii >= jj) break;
						else {
							addcut[ii] = addcut[jj];
							addcut[jj].value = 0.0;
							addcut.pop_back();
						}
					}
					if (!found) break;
				}
				i += j;
				j=0;
			    if (addcut.size()>0) {
					listOfCutsLhs.push_back( addcut );
					listOfCutsRhs.push_back( rrhs );
					if(Print){for (int pp=0;pp<addcut.size();pp++){
						cerr << addcut[pp].value.asDouble() <<"x_"<<addcut[pp].index << endl;
					}
					cerr << ">= " << rrhs << endl;
					}
				}

			}
			if (info_level > 2) cerr << "CMIRS:"<< cuts.size() << " " << listOfCutsRhs.size() << endl;
	//cin.get();
			return listOfCutsRhs.size(); //number of cuts

		}

	void QBPSolver::ProcessConstraintWithIntKnowledge(){
				return;
			IniBinarizedVec();
				cerr << "Initialized Binarized Vektor " << Binarized.size() << endl;
				//Try Find Slacks
				double eps=+1e-7;
				bool Print=false;
				bool UseLifting=true;
				bool UseIntLifting=false;
				bool UseBoundStrengthening =true;
				coef_t RowMax;
				coef_t RowMin;
				//coef_t OldMax;
				//coef_t OldMin;
				vector<CRef> CheckAgain;
				CheckAgain.clear();
				int ActionCounter=0;
				for (int i = 0; i < constraints.size()  +CheckAgain.size() && i<3*constraints.size();i++) {

					//cerr << "Run number " << i << " of " << constraints.size() +CheckAgain.size() << endl;
					Constraint &c =(i<constraints.size()?constraintallocator[constraints[i]]:constraintallocator[CheckAgain[i-constraints.size()]]);
					if(i<constraints.size() && constraints[i]==constraints[0]) continue;
					if(i>constraints.size()&& CheckAgain[i-constraints.size()]==constraints[0]) continue;
					RowMax=0.0;
					RowMin=0.0;
					//OldMax=0.0;
					//OldMin=0.0;
					bool SomeFixed=false;
					bool Interrupted=false;
					if(UseLifting || UseBoundStrengthening||UseIntLifting){
						if (c.size()<=1 || c.header.isVarBnd || c.header.isIntBnd || i==0 || c.header.isSat /*|| c.header.learnt*/ ) continue;
						coef_t rhs = c.header.rhs;
						coef_t LargestAbsCoef = -1;
						int IndexOfLargest = -1;
						for(int h=0 ; h<c.size() ; h++){
							if(upperBounds[var(c[h])]>=999999){
								Interrupted=true;
								break;
							}
							/*if(abs(upperBounds[var(c[h])]-lowerBounds[var(c[h])])<eps){
								Interrupted=true;
								break;
							}*/
							if(Print) cerr <<" +" <<(sign(c[h])?-c[h].coef:c[h].coef) << "x_" << var(c[h]) << "["<<lowerBounds[var(c[h])]<<","<<upperBounds[var(c[h])]<< "]";
							if(assigns[var(c[h])]!=extbool_Undef || isFixed(var(c[h]))){
								SomeFixed=true;
								if(type[var(c[h])]!=BINARY){
									Interrupted=true;
									break;
								}
								if(assigns[var(c[h])]!=extbool_Undef){
									if(((yInterface*)yIF)->integers[var(c[h])].bitcnt==1){
										if(Print) cerr << "Use assigned var x_" << var(c[h]) <<"="<<(int)assigns[var(c[h])]<< endl;

										if(sign(c[h])){
											RowMax-=c[h].coef*assigns[var(c[h])];
											RowMin-=c[h].coef*assigns[var(c[h])];
										}
										else{
											RowMax+=c[h].coef*assigns[var(c[h])];
											RowMin+=c[h].coef*assigns[var(c[h])];
										}
									}
									else{
										int BinInd=((yInterface*)yIF)->integers[var(c[h])].number;
										if(Binarized[BinInd].ub>=999999){
											Interrupted=true;
											break;
										}
										if(Print) cerr << "Dealt with assigned INT" << " OldLB of Int: " << Binarized[BinInd].lb << " OldUB of INT: " << Binarized[BinInd].ub  << endl;
										int Nlb=0;
										int Nub=0;
										for(int bv=Binarized[BinInd].FirstIndex;bv<=Binarized[BinInd].LastIndex;bv++){
											int pot=pow(2,Binarized[BinInd].LastIndex-bv);
											if(assigns[bv]==1){
												Nlb+=pot;
												Nub+=pot;
											}
											else if(assigns[bv]==extbool_Undef){
												Nub+=pot;
											}
										}
										if(Nub<=Binarized[BinInd].ub)
											Binarized[BinInd].ub=Nub;
										if(Nlb>=Binarized[BinInd].lb)
											Binarized[BinInd].lb=Nlb;
										if(Print)cerr <<"NLB: " << Nlb << " NUB: " << Nub << endl;
										if(Print)cerr << "x_" << var(c[h]) << " assigned to " << (int)assigns[var(c[h])] << " LB of Int: " << Binarized[BinInd].lb << " UB of INT: " << Binarized[BinInd].ub  << endl;
										bool FoundRep=false;
										if(((yInterface*)yIF)->integers[var(c[h])].pt2leader+((yInterface*)yIF)->integers[var(c[h])].bitcnt-1==var(c[h])){
											FoundRep=true;
											if(Print)cerr <<"Int Var:  " << ((yInterface*)yIF)->integers[var(c[h])].pt2leader <<"-" <<Binarized[BinInd].LastIndex << "with "<< Binarized[BinInd].lb << ","<<Binarized[BinInd].ub <<endl;
											assert(Binarized[BinInd].LastIndex==var(c[h]));
											assert(Binarized[BinInd].ub>=Binarized[BinInd].lb);
											if(sign(c[h])){
												RowMax-=c[h].coef*Binarized[BinInd].lb;
												RowMin-=c[h].coef*Binarized[BinInd].ub;
											}
											else{
												RowMax+=c[h].coef*Binarized[BinInd].ub;
												RowMin+=c[h].coef*Binarized[BinInd].lb;
											}
										}
										if(!FoundRep){
											for(int co=0;co<c.size();co++){
												if(var(c[co])==((yInterface*)yIF)->integers[var(c[h])].pt2leader+((yInterface*)yIF)->integers[var(c[h])].bitcnt-1){
													FoundRep=true;
													break;
												}
											}
										}
										if(!FoundRep){
											//cerr <<"Representative Missing. Constraint skipped." << endl;
											break;
										}
									}
								}
								else{
								Interrupted=true;
								if (Print) cerr <<"(ASSIGNED: x_"<< var(c[h]) << " "<<(int)assigns[var(c[h])] << ")";
								break;
								}
							}
							else{
								/*if(sign(c[h])){
									OldMax-=c[h].coef*lowerBounds[var(c[h])];
									OldMin-=c[h].coef*upperBounds[var(c[h])];
								}
								else{
									OldMax+=c[h].coef*upperBounds[var(c[h])];
									OldMin+=c[h].coef*lowerBounds[var(c[h])];
								}*/
								if(((yInterface*)yIF)->integers[var(c[h])].bitcnt==1){
									if(upperBounds[var(c[h])]>=999999){
										Interrupted=true;
										break;
									}
									if(sign(c[h])){
										RowMax-=c[h].coef*lowerBounds[var(c[h])];
										RowMin-=c[h].coef*upperBounds[var(c[h])];
									}
									else{
										RowMax+=c[h].coef*upperBounds[var(c[h])];
										RowMin+=c[h].coef*lowerBounds[var(c[h])];
									}
								}
								else{

									if(((yInterface*)yIF)->integers[var(c[h])].pt2leader+((yInterface*)yIF)->integers[var(c[h])].bitcnt-1==var(c[h])){
										int BinInd=((yInterface*)yIF)->integers[var(c[h])].number;
										if(Binarized[BinInd].ub>=999999){
											Interrupted=true;
											break;
										}
										if(!PresentAsInt(c,((yInterface*)yIF)->integers[var(c[h])].pt2leader,Binarized[BinInd].LastIndex,sign(c[h]),c[h].coef)){
											Interrupted=true;
											break;
										}
										if (Print) cerr <<"Int Var:  " << ((yInterface*)yIF)->integers[var(c[h])].pt2leader <<"-" <<Binarized[BinInd].LastIndex << "with "<< Binarized[BinInd].lb << ","<<Binarized[BinInd].ub <<endl;

										assert(Binarized[BinInd].LastIndex==var(c[h]));
										assert(Binarized[BinInd].ub>=Binarized[BinInd].lb);
										if(sign(c[h])){
											RowMax-=c[h].coef*Binarized[BinInd].lb;
											RowMin-=c[h].coef*Binarized[BinInd].ub;
										}
										else{
											RowMax+=c[h].coef*Binarized[BinInd].ub;
											RowMin+=c[h].coef*Binarized[BinInd].lb;
										}
									}
									else{
										int BinInd=((yInterface*)yIF)->integers[var(c[h])].number;
										if(!PresentAsInt(c,((yInterface*)yIF)->integers[var(c[h])].pt2leader,Binarized[BinInd].LastIndex,sign(c[h]),c[h].coef/(double)(pow(2,Binarized[BinInd].LastIndex-var(c[h]))))){
											Interrupted=true;
											break;
										}
									}
								}

								if(UseLifting){
									if(c[h].coef>LargestAbsCoef && type[var(c[h])]==BINARY &&  ((yInterface*)yIF)->integers[var(c[h])].bitcnt==1){ // && !IsGeneral(var(c[h]))
										LargestAbsCoef=c[h].coef;
										IndexOfLargest = h;
									}

								}
							}
						}
						if(Print)cerr << ">= " << rhs << endl;
						if(Interrupted) continue;

						if(Print) cerr << "Row ana: Max: " << c.header.btch1.best << " >= "  << RowMax << endl;
						if(Print)cerr << "Row ana: Min: " << c.header.wtch2.worst << " <= "  << RowMin << endl;
						if(!c.header.isSat){
							if(c.header.btch1.best < RowMax) continue;
							if(c.header.wtch2.worst > RowMin) continue;
							assert(c.header.btch1.best >= RowMax );
							assert(c.header.wtch2.worst <= RowMin );
						}

						bool SomeBoundStrengthened=false;
						bool CoefStrengthened=false;
						if(UseIntLifting){
							for( int IndexStrengthen=0;IndexStrengthen<c.size(); IndexStrengthen++){
								if(assigns[var(c[IndexStrengthen])]!=extbool_Undef || isFixed(var(c[IndexStrengthen]))) continue;
								bool BoundStrengthened=false;
								if (((yInterface*)yIF)->integers[var(c[IndexStrengthen])].bitcnt>1){
									if(((yInterface*)yIF)->integers[var(c[IndexStrengthen])].pt2leader+((yInterface*)yIF)->integers[var(c[IndexStrengthen])].bitcnt-1==var(c[IndexStrengthen])){
										//cerr<< "IntIfno "<< ((yInterface*)yIF)->integers[var(c[IndexStrengthen])].bitcnt<< " "<<((yInterface*)yIF)->integers[var(c[IndexStrengthen])].org_ub << endl;
										int BinInd=((yInterface*)yIF)->integers[var(c[IndexStrengthen])].number;
										assert(Binarized[BinInd].LastIndex==var(c[IndexStrengthen]));
										assert(Binarized[BinInd].ub<999999);

										if(sign(c[IndexStrengthen])){
											RowMin+=c[IndexStrengthen].coef*Binarized[BinInd].ub;
										}
										else{
											RowMin-=c[IndexStrengthen].coef*Binarized[BinInd].lb;
										}

										if(sign(c[IndexStrengthen])){
											double d= -c.header.rhs + RowMin - c[IndexStrengthen].coef*(Binarized[BinInd].ub-1);
											if(d>eps && c[IndexStrengthen].coef>=d){
												bool NumFound=0;
												for(int ii=0;ii<c.size();ii++){
													if(ii!=IndexStrengthen &&((yInterface*)yIF)->integers[var(c[IndexStrengthen])].pt2leader==((yInterface*)yIF)->integers[var(c[ii])].pt2leader){
														NumFound=1;
														if(sign(c[ii])!=sign(c[IndexStrengthen]) ||
																c[ii].coef !=c[IndexStrengthen].coef*pow(2,((yInterface*)yIF)->integers[var(c[ii])].pt2leader+((yInterface*)yIF)->integers[var(c[ii])].bitcnt-1- var(c[ii]) )){
															if(Print)cerr << "Wacky Int Coefficients" << endl;
															NumFound=-1;
															break;
														}
														//c[ii].coef=(c[IndexStrengthen].coef-d)*pow(2,((yInterface*)yIF)->integers[var(c[ii])].pt2leader+((yInterface*)yIF)->integers[var(c[ii])].bitcnt-1- var(c[ii]));
														//if(Print)cerr << "New Coef for x_" <<  var(c[ii])  << " is " << c[ii].coef << endl;
													}

												}
												if(NumFound==((yInterface*)yIF)->integers[var(c[IndexStrengthen])].bitcnt-1){
													for(int ii=0;ii<c.size();ii++){
														if(ii!=IndexStrengthen &&((yInterface*)yIF)->integers[var(c[IndexStrengthen])].pt2leader==((yInterface*)yIF)->integers[var(c[ii])].pt2leader){
															c[ii].coef=(c[IndexStrengthen].coef-d)*pow(2,((yInterface*)yIF)->integers[var(c[ii])].pt2leader+((yInterface*)yIF)->integers[var(c[ii])].bitcnt-1- var(c[ii]));
															if(Print)cerr << "New Coef for x_" <<  var(c[ii])  << " is " << c[ii].coef << endl;
														}

													}
													if(Print)cerr <<"Coeff Strengthening(+)! d=" << d << endl;

													c[IndexStrengthen].coef=(c[IndexStrengthen].coef-d);
													if(Print)cerr << "New Coef for x_" <<  var(c[IndexStrengthen])  << " is " << c[IndexStrengthen].coef << endl;

													c.header.rhs=c.header.rhs-d*Binarized[BinInd].ub;
													if(Print)cerr << "NEW RHS: "<< c.header.rhs << endl;
													CoefStrengthened=true;
													ActionCounter++;

												}
											}
										}
										else{
											double d= -c.header.rhs + RowMin + c[IndexStrengthen].coef*(Binarized[BinInd].lb+1);
											if(d>eps && c[IndexStrengthen].coef>=d ){
												bool NumFound=0;
												for(int ii=0;ii<c.size();ii++){
													if(ii!=IndexStrengthen &&((yInterface*)yIF)->integers[var(c[IndexStrengthen])].pt2leader==((yInterface*)yIF)->integers[var(c[ii])].pt2leader){
														NumFound=1;
														if(sign(c[ii])!=sign(c[IndexStrengthen]) ||
																c[ii].coef !=c[IndexStrengthen].coef*pow(2,((yInterface*)yIF)->integers[var(c[ii])].pt2leader+((yInterface*)yIF)->integers[var(c[ii])].bitcnt-1- var(c[ii]) )){
															if(Print) cerr << "Wacky Int Coefficients" << endl;
															NumFound=-1;
															break;
														}
														//c[ii].coef=(c[IndexStrengthen].coef-d)*pow(2,((yInterface*)yIF)->integers[var(c[ii])].pt2leader+((yInterface*)yIF)->integers[var(c[ii])].bitcnt-1- var(c[ii]));
														if(Print) cerr << "New Coef for x_" <<  var(c[ii])  << " is " << c[ii].coef << endl;
													}
												}
												if(NumFound==((yInterface*)yIF)->integers[var(c[IndexStrengthen])].bitcnt-1){
													if(Print) cerr <<"Coeff Strengthening(-)! d=" << d << endl;
													for(int ii=0;ii<c.size();ii++){
														if(ii!=IndexStrengthen &&((yInterface*)yIF)->integers[var(c[IndexStrengthen])].pt2leader==((yInterface*)yIF)->integers[var(c[ii])].pt2leader){
															c[ii].coef=(c[IndexStrengthen].coef-d)*pow(2,((yInterface*)yIF)->integers[var(c[ii])].pt2leader+((yInterface*)yIF)->integers[var(c[ii])].bitcnt-1- var(c[ii]));
															if(Print) cerr << "New Coef for x_" <<  var(c[ii])  << " is " << c[ii].coef << endl;
														}

													}
													c[IndexStrengthen].coef=(c[IndexStrengthen].coef-d);
													if(Print)cerr << "New Coef for x_" <<  var(c[IndexStrengthen])  << " is " << c[IndexStrengthen].coef << endl;

													c.header.rhs=c.header.rhs-d*Binarized[BinInd].lb;
													if(Print)cerr << "NEW RHS: "<< c.header.rhs << endl;
													CoefStrengthened=true;
													ActionCounter++;
												}
											}

										}
										if(sign(c[IndexStrengthen])){
											RowMin-=c[IndexStrengthen].coef*Binarized[BinInd].ub;
										}
										else{
											RowMin+=c[IndexStrengthen].coef*Binarized[BinInd].lb;
										}
									}
									if(CoefStrengthened){
									if(i<constraints.size())
										CheckAgain.push_back(constraints[i]);
									else
										CheckAgain.push_back(CheckAgain[i-constraints.size()]);
									}
								}
							}
						}
						if(UseBoundStrengthening &&!CoefStrengthened){
							for( int IndexStrengthen=0;IndexStrengthen<c.size(); IndexStrengthen++){
								if(assigns[var(c[IndexStrengthen])]!=extbool_Undef || isFixed(var(c[IndexStrengthen]))) continue;
								bool BoundStrengthened=false;
								if (((yInterface*)yIF)->integers[var(c[IndexStrengthen])].bitcnt>1){
									if(((yInterface*)yIF)->integers[var(c[IndexStrengthen])].pt2leader+((yInterface*)yIF)->integers[var(c[IndexStrengthen])].bitcnt-1==var(c[IndexStrengthen])){
										//cerr<< "IntIfno "<< ((yInterface*)yIF)->integers[var(c[IndexStrengthen])].bitcnt<< " "<<((yInterface*)yIF)->integers[var(c[IndexStrengthen])].org_ub << endl;
										int BinInd=((yInterface*)yIF)->integers[var(c[IndexStrengthen])].number;
										assert(Binarized[BinInd].LastIndex==var(c[IndexStrengthen]));
										assert(Binarized[BinInd].ub<999999);
										if(!PresentAsInt(c,((yInterface*)yIF)->integers[var(c[IndexStrengthen])].pt2leader,Binarized[BinInd].LastIndex,sign(c[IndexStrengthen]),c[IndexStrengthen].coef)){
											cerr << var(c[IndexStrengthen]) << " " << ((yInterface*)yIF)->integers[var(c[IndexStrengthen])].pt2leader << " " << Binarized[BinInd].LastIndex << " " <<sign(c[IndexStrengthen]) << " " <<c[IndexStrengthen].coef;

											cin.get();
										}

										if(sign(c[IndexStrengthen])){
											RowMax+=c[IndexStrengthen].coef*Binarized[BinInd].lb;
										}
										else{
											RowMax-=c[IndexStrengthen].coef*Binarized[BinInd].ub;
										}
										//cerr << "Check Integer Bound " << endl;
										if(sign(c[IndexStrengthen]) &&Binarized[BinInd].ub> floor((RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps)){
											if(Print) cerr <<"FOUND INTEGER BOUND TO BE STRENGTHENED!!!!! x_" << var(c[IndexStrengthen]) << " OrgUB: " << Binarized[BinInd].ub <<
													" NewUB: " << floor((RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps) << endl;
											ActionCounter++;
											BoundStrengthened=true;
											Binarized[BinInd].ub=floor((RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps);
											if(  Binarized[BinInd].ubCon==0){
												if(Print) cerr << "before LS restriction: " << constraints.size() << endl;

												int remPrimalRestriction = constraints.size()-1;
												ca_vec<CoeVar> LHS;
												for(int va= Binarized[BinInd].FirstIndex ;va<=Binarized[BinInd].LastIndex;va++){
													LHS.push( mkCoeVar(va, pow(2,Binarized[BinInd].LastIndex-va), true));
												}
												addOrgConstraint(LHS,-Binarized[BinInd].ub-(1e-10)-abs(Binarized[BinInd].ub)*1e-10,0,false,-1,false);
												Binarized[BinInd].ubCon = constraints[constraints.size()-1];
												 assert(constraintallocator[Binarized[BinInd].ubCon].header.isIntBnd);
												 assert(!constraintallocator[Binarized[BinInd].ubCon].header.learnt);

												 if(Print) cerr << "after LS restriction: " << constraints.size() << endl;

											}
											else{
												constraintallocator[Binarized[BinInd].ubCon].header.rhs=-Binarized[BinInd].ub-(1e-10)-abs(Binarized[BinInd].ub)*1e-10;
												if(Print){
													cerr <<"Changed UBConstraint: " << endl;
												Constraint &t= constraintallocator[Binarized[BinInd].ubCon];
												for (int pp=0;pp<t.size();pp++){
													if(sign(t[pp])) cerr <<"-" << t[pp].coef<<"x_" <<var(t[pp]);
													else cerr <<"+" << t[pp].coef<<"x_" <<var(t[pp]);
												}
												cerr << ">=" << t.header.rhs << endl;
												}
											}


										}
										else if(!sign(c[IndexStrengthen]) && Binarized[BinInd].lb< ceil((c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps)){
											if(Print) cerr <<"FOUND INTEGER BOUND TO BE STRENGTHENED!!!!! x_" << var(c[IndexStrengthen]) << " OrgLB: " << Binarized[BinInd].lb <<
													" NewLB: " <<  ceil((c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps) << endl;
											ActionCounter++;
											BoundStrengthened=true;
											Binarized[BinInd].lb= ceil((c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps);
											if(  Binarized[BinInd].lbCon==0){
												if(Print)   cerr << "before LS restriction: " << constraints.size() << endl;

												 int remPrimalRestriction = constraints.size()-1;
												 ca_vec<CoeVar> LHS;
												 for(int va= Binarized[BinInd].FirstIndex ;va<=Binarized[BinInd].LastIndex;va++){
													 LHS.push( mkCoeVar(va, pow(2,Binarized[BinInd].LastIndex-va), false));
												 }
												 addOrgConstraint(LHS,Binarized[BinInd].lb-(1e-10)-abs(Binarized[BinInd].lb)*1e-10,0,false,-1,false);
												 Binarized[BinInd].lbCon = constraints[constraints.size()-1];
												 if(Print)  cerr << "after LS restriction: " << constraints.size() << endl;

												 assert(constraintallocator[Binarized[BinInd].lbCon].header.isIntBnd);
												 assert(!constraintallocator[Binarized[BinInd].lbCon].header.learnt);

											}
											else constraintallocator[Binarized[BinInd].lbCon].header.rhs=Binarized[BinInd].lb;
										}



										if(sign(c[IndexStrengthen])){
											RowMax-=c[IndexStrengthen].coef*Binarized[BinInd].lb;
										}
										else{
											RowMax+=c[IndexStrengthen].coef*Binarized[BinInd].ub;
										}
									}
								}
								else{
									assert(upperBounds[var(c[IndexStrengthen])]<999999);

									if(sign(c[IndexStrengthen])){
										RowMax+=c[IndexStrengthen].coef*lowerBounds[var(c[IndexStrengthen])];
									}
									else{
										RowMax-=c[IndexStrengthen].coef*upperBounds[var(c[IndexStrengthen])];
									}
									if(type[var(c[IndexStrengthen])]==CONTINUOUS){
										if(sign(c[IndexStrengthen]) &&upperBounds[var(c[IndexStrengthen])]> (RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps){
											if(abs(upperBounds[var(c[IndexStrengthen])]- ((RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps))>1000*eps && (RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps<1e8){
												if (Print) cerr <<" Found Bound To Be Strengthened (UB) in .h! For x_" << var(c[IndexStrengthen]) << ". Old Bound: " <<upperBounds[var(c[IndexStrengthen])] << ". NewBound: "<< (double)(RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef<< endl;
												ActionCounter++;
												setUpperBound(var(c[IndexStrengthen]),(RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps);
												BoundStrengthened=true;
												QlpStSolve->setVariableUB(var(c[IndexStrengthen]),(RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps,0);
												updateStageSolver(maxLPStage,var(c[IndexStrengthen]),var(c[IndexStrengthen]));
											}
										}
										else if(!sign(c[IndexStrengthen]) && lowerBounds[var(c[IndexStrengthen])]< (c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps){
											if(abs(lowerBounds[var(c[IndexStrengthen])]-( (c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps))>1000*eps  &&(c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps<1e8){
												if(Print) cerr <<" Found Bound To Be Strengthened (LB) in .h! For x_" << var(c[IndexStrengthen]) << ". Old Bound: " <<lowerBounds[var(c[IndexStrengthen])] << ". NewBound: "<< (c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef<< endl;
												ActionCounter++;
												setLowerBound(var(c[IndexStrengthen]),(c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps);
												BoundStrengthened=true;
												QlpStSolve->setVariableLB(var(c[IndexStrengthen]),(c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps,0);
												updateStageSolver(maxLPStage,var(c[IndexStrengthen]),var(c[IndexStrengthen]));
											}
										}
									}
									else{
										if(sign(c[IndexStrengthen]) &&upperBounds[var(c[IndexStrengthen])]> floor((RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps)){
											if (Print) cerr <<" Found Bound To Be Strengthened (UB) in .h! For x_" << var(c[IndexStrengthen]) << ". Old Bound: " <<upperBounds[var(c[IndexStrengthen])] << ". NewBound: "<< (double)(RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef<< endl;
											ActionCounter++;
											setUpperBound(var(c[IndexStrengthen]),floor((RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps));
											BoundStrengthened=true;
											QlpStSolve->setVariableUB(var(c[IndexStrengthen]),floor((RowMax-c.header.rhs)/(double)c[IndexStrengthen].coef+eps),0);
											updateStageSolver(maxLPStage,var(c[IndexStrengthen]),var(c[IndexStrengthen]));
											//if(type[A[IndexConti].index]==BINARY && abs(qbp->getUpperBound(A[IndexConti].index)-qbp->getLowerBound(A[IndexConti].index))<eps)
											//	qbp->setFixed(A[IndexConti].index,qbp->getUpperBound(A[IndexConti].index));
										}
										else if(!sign(c[IndexStrengthen]) && lowerBounds[var(c[IndexStrengthen])]< ceil((c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps)){
											if(Print) cerr <<" Found Bound To Be Strengthened (LB) in .h! For x_" << var(c[IndexStrengthen]) << ". Old Bound: " <<lowerBounds[var(c[IndexStrengthen])] << ". NewBound: "<< (c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef<< endl;
											ActionCounter++;
											setLowerBound(var(c[IndexStrengthen]),ceil((c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps));
											BoundStrengthened=true;
											QlpStSolve->setVariableLB(var(c[IndexStrengthen]),ceil((c.header.rhs-RowMax)/(double)c[IndexStrengthen].coef-eps),0);
											updateStageSolver(maxLPStage,var(c[IndexStrengthen]),var(c[IndexStrengthen]));
											//if(type[A[IndexConti].index]==BINARY && abs(qbp->getUpperBound(A[IndexConti].index)-qbp->getLowerBound(A[IndexConti].index))<eps)
											//	qbp->setFixed(A[IndexConti].index,qbp->getUpperBound(A[IndexConti].index));
										}
									}
									if(sign(c[IndexStrengthen])){
										RowMax-=c[IndexStrengthen].coef*lowerBounds[var(c[IndexStrengthen])];
									}
									else{
										RowMax+=c[IndexStrengthen].coef*upperBounds[var(c[IndexStrengthen])];
									}
								}
								if(BoundStrengthened){
									SomeBoundStrengthened=true;
									for (int k=0;k<VarsInConstraints[var(c[IndexStrengthen])].size();k++){
										CheckAgain.push_back(VarsInConstraints[var(c[IndexStrengthen])][k].cr);
									}
								}

							}
						}

						if(UseLifting && !SomeBoundStrengthened &&!CoefStrengthened && IndexOfLargest!=-1 &&!SomeFixed){
							bool Lifted=false;
							assert(  ((yInterface*)yIF)->integers[var(c[IndexOfLargest])].org_ub==1);
							assert(  ((yInterface*)yIF)->integers[var(c[IndexOfLargest])].number==-1);
							if(!sign(c[IndexOfLargest])){
								//Positive Variable; Check, wether setting it to 1 results in nonnegative slack

								if(RowMin+LargestAbsCoef>c.header.rhs+eps){
									if (Print)cerr << "Found something to lift(+)! Var " << IndexOfLargest << " is x_" << var(c[IndexOfLargest]) << ", RowMin is " <<RowMin << " and Slack is " << RowMin+LargestAbsCoef-c.header.rhs << " and new Coef is " <<  c[IndexOfLargest].coef-RowMin-LargestAbsCoef+c.header.rhs << endl;
									ActionCounter++;
									if(c[IndexOfLargest].coef-(RowMin+LargestAbsCoef-c.header.rhs)<0){
										cerr <<"What happened?" <<endl;
										continue;
									}
									c[IndexOfLargest].coef -= RowMin+LargestAbsCoef-c.header.rhs;
									assert(c[IndexOfLargest].coef>0);
									if(abs(c[IndexOfLargest].coef-floor(c[IndexOfLargest].coef+.5))<eps) c[IndexOfLargest].coef=floor(c[IndexOfLargest].coef+.5);
									Lifted=true;
								}
							}
							else if(sign(c[IndexOfLargest]) && RowMin-c.header.rhs>eps){
								coef_t Sla=RowMin-c.header.rhs;
								if(Print)cerr << "Found something to lift(-)! Var " << IndexOfLargest << " is x_" << var(c[IndexOfLargest]) << ", L is " <<RowMin << " and Slack is " << Sla<< endl;
								ActionCounter++;
								if(c[IndexOfLargest].coef-Sla<0){
									cerr <<"What happened?" <<endl;
									continue;
								}
								c[IndexOfLargest].coef-=Sla;
								assert(c[IndexOfLargest].coef>0);
								if(abs(c[IndexOfLargest].coef-floor(c[IndexOfLargest].coef+.5))<eps) c[IndexOfLargest].coef=floor(c[IndexOfLargest].coef+.5);
								c.header.rhs=RowMin-(1e-12)-abs(RowMin)*1e-12;
								Lifted=true;
							}


							if(Lifted){

							if(i<constraints.size())
								CheckAgain.push_back(constraints[i]);
							else
								CheckAgain.push_back(CheckAgain[i-constraints.size()]);
							}
							//CheckAgain.push_back(c.cr_index_varix);
						}
					}
				}
				//cerr <<"Actions performed: " << ActionCounter << endl;
			}
		void QBPSolver::IniBinarizedVec(){
		  //return;
			//assert(nVars()==((yInterface*)yIF)->integers.size());
			if(Binarized.size()==0){
				for(int i=0;i<nVars();i++){
					if(((yInterface*)yIF)->integers[i].bitcnt>1){
						IntInfo VarI;
						VarI.org_index = ((yInterface*)yIF)->integers[i].org_ind;
						VarI.ubCon=0;
						VarI.lbCon=0;
						VarI.ub=((yInterface*)yIF)->integers[i].org_ub.asDouble()-((yInterface*)yIF)->integers[i].org_lb.asDouble();
						VarI.lb =0;
						VarI.FirstIndex=i;
						while(((yInterface*)yIF)->integers[i].org_ind==VarI.org_index && i<nVars()){
							((yInterface*)yIF)->integers[i].number=Binarized.size();
							VarI.LastIndex=i;
							i++;
						}
						Binarized.push_back(VarI);

						i--;
					}
				}
			}
			else{
				for (int i=0;i<Binarized.size();i++){
					Binarized[i].lbCon=0;
					Binarized[i].ubCon=0;
				}
			}

			int count=0;
			for (int i = 0; i < constraints.size() ;i++) {
				Constraint &c =constraintallocator[constraints[i]];
				if(c.header.isIntBnd){
					count++;
					int index=((yInterface*)yIF)->integers[var(c[0])].number;
					assert(index!=-1);
					for(int k=0;k<c.size();k++)
						assert(((yInterface*)yIF)->integers[var(c[k])].number==index);
					if(sign(c[0])){
						assert(c.header.rhs<=0);
						Binarized[index].ubCon=constraints[i];
					}
					else{
						assert(c.header.rhs>=0);

						Binarized[index].lbCon=constraints[i];
					}
				}
				else{
					continue;
					//assert that it really is no Int Bound
					if(c.size()>1 && !c.header.isSat){
						int FirstVar=((yInterface*)yIF)->integers[var(c[0])].org_ind;
						bool FoundOther=false;
						for(int k=0;k<c.size();k++){
							if (FirstVar!=((yInterface*)yIF)->integers[var(c[k])].org_ind){
								FoundOther=true;
								break;
							}
						}
						if(!FoundOther){
							for(int k=0;k<c.size();k++){
								cerr << "+"<<(sign(c[k])?-c[k].coef:c[k].coef) <<"x_"<<var(c[k]);
							}
							cerr <<">= " << c.header.rhs << endl;
							cerr << "Found Not marked IntBound" << endl;
							cin.get();
						}
					}
				}
			}
			//cerr <<"Found " << count << " IntBoundConstraints" << endl;

			if(0){
				for(int i=0;i<Binarized.size();i++){
					if(Binarized[i].ubCon==0){
						/*cerr << i << " Org: x_" << Binarized[i].org_index << "[" <<Binarized[i].lb << "," << Binarized[i].ub << "] " << Binarized[i].FirstIndex << " " << Binarized[i].LastIndex <<endl;
						for(int k=Binarized[i].FirstIndex;k<=Binarized[i].LastIndex;k++){
							cerr << "  Now is: " << pow(2,Binarized[i].LastIndex-k) << "x_" << k << "[" << upperBounds[k] << ","<<lowerBounds[k] <<"]" <<endl;
							//cerr << ((yInterface*)yIF)->integers[k+1].number << " " << ((yInterface*)yIF)->integers[k+1].org_ind << endl;
						}
						*/
						//assert(Binarized[i].ub-Binarized[i].lb+1==pow(2,Binarized[i].LastIndex-Binarized[i].FirstIndex+1));
					}

				}
			}

		}

		bool QBPSolver::PresentAsInt(Constraint& c, int First, int Last, bool sgn ,coef_t coef){
			int NFound=0;
			for(int i=0;i<c.size() && NFound<=Last-First+1;i++){
				if(var(c[i])<=Last && var(c[i])>=First){
					if(sign(c[i])==sgn && c[i].coef == coef*pow(2,Last-var(c[i])))
						NFound++;
					else
					{
						//cerr << "Does not fit " << endl;
						NFound=-1;
						break;
					}

				}
			}
			if(NFound==Last-First+1)
				return true;
			else return false;
		}

		double QBPSolver::getObj(int var){
			double ObjCoeff=0;
			for (int i=0; i < VarsInConstraints[var].size();i++) {
			    Constraint &c = constraintallocator[VarsInConstraints[var][i].cr];
			    if (VarsInConstraints[var][i].cr != constraints[0]) {
			      	continue;
			    }
			    else{
			      	int pos = VarsInConstraints[var][i].pos;
			      	if(sign(c[pos])) ObjCoeff=-c[pos].coef;
			      	else ObjCoeff=c[pos].coef;
			    }
			}
			return ObjCoeff;
		}

		bool QBPSolver::PresentAsInt(vector<data::IndexedElement >& c, int First, int Last, coef_t coef){
			int NFound=0;
			for(int i=0;i<c.size() && NFound<=Last-First+1;i++){
				if(c[i].index<=Last && c[i].index>=First){
					if(c[i].value == coef*pow(2,Last-c[i].index))
						NFound++;
					else
					{
						//cerr << "Does not fit " << endl;
						NFound=-1;
						break;
					}

				}
			}
			if(NFound==Last-First+1)
				return true;
			else return false;
		}

