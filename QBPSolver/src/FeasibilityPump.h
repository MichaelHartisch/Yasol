/*
*
* Yasol: FeasibilityPump.h -- Copyright (c) 2017 Ulf Lorenz, Michael Hartisch
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

#ifndef FEASIBILITYPUMP_H_
#define FEASIBILITYPUMP_H_

#include "Utilities/QlpStageSolver.hpp" 


#ifdef USE_NBD_CPLEX_C
#include "ExternSolvers/QpExtSolCplexC.hpp"
#include "cmath"
#endif

static	vector<double> PrevIPSol;


class FeasibilityPump {

	//Maximum Number of LPs solved = Number of Pumping Rounds
	int MaxIter;
	int MaxTime;
	int MaxLong;
	int MaxDead;
	int MaxMultiTime;
	int TimeLong; //Time defined as "long run"

	bool SolutionFound;
	//Storing the old Coefficients and the solution of the LP relaxation
	std::vector<data::QpNum> SaveObjCoeffs;
	std::vector<data::QpNum> OriginalLPSolution;

	//include the original objective coefficients in the search; 
	//alpha specifies the importance of the original objective
	double alpha;
	//Decay of alpha per round
	double Decay;
	//maximum difference in alpha values for the dead ends to be compared
	double delta_alpha;

	double NormC;			//Norm of objective coefficient vector
	double NormDelta;		//sqrt(Number of binary variables)
	double Normalization;	//Normalization for including the original objective
	double MaxNormObj;		//max(abs(c_i))

	bool IsInteger;
	double RandSeed;
	double epsNull;
	unsigned int n;



	//For the StageSolver:
	int stage;
	utils::QlpStageSolver *StSolver;
	algorithm::Algorithm::SolutionStatus status;
	data::QpNum lb;
	data::QpNum ub;	
	
	// For each variable .first the .second value stores |x-[x]|, i.e. how far this variable is away from its rounded value
	// Vector for the absolute values |x_j-\tilde x_j|; needed when flipping in case of stalling
	std::vector< std::pair<int, double> > Fractionality; 

	//Vector for sorting variable order  
	std::vector< std::pair<int, double> > RandomOrder; 


	//If hard-restarting a variable i is flipped with propability: fractionality(i)+Restartp
	double RestartP;

	//If normal-perturbing flip x variables with T/2<=x<=3T/2, 
	int T;

	//Number of memorized previous runs
	int MaxStalls;

	//Queue storing the last MaxStalls Overall-Fractionalities
	std::deque<double> PreviousFrac;

		// Vectors to check, wheather we stopped in the same dead end again
	vector<vector<double>> DeadEnds;
	vector<double> DeadAlpha;
	vector<int> DeadCount;


public:
	FeasibilityPump(std::vector<data::QpNum> &LPSol,utils::QlpStageSolver *qss,data::Qlp &qlp, int stageIndex, int *types, int nVa, int mI=1000)  {
		const std::vector<data::QpNum>& SaveObjCoeffsTmp = qlp.getObjectiveFunctionValues();
		SaveObjCoeffs.clear();
		for (int i = 0; i < nVa;i++) SaveObjCoeffs.push_back(SaveObjCoeffsTmp[i]);
		RestartP=0.0;
		MaxIter=mI;
		MaxTime=1800;
		MaxLong=mI/10;
		MaxMultiTime=10;
		MaxDead=5;
		TimeLong=10;
		SolutionFound=false;
		RandSeed=time(NULL);
		MaxNormObj=0;
		MaxStalls=200;
		IsInteger=false;
		StSolver=qss;
		epsNull = 1e-7;
		stage=stageIndex;
		n=nVa;//StSolver->getExternSolver(stageIndex).getVariableCount();
		OriginalLPSolution=LPSol;
		T = 20;
		alpha=1;
		delta_alpha=0.005;
		Decay=1.0;//0.999;

		//Computing normalization term
		NormC=0;
		NormDelta=0;
		for (int i=0;i<n;i++){
			NormDelta+=(int)(types[i]==0);
			NormC+=SaveObjCoeffs[i].asDouble()*SaveObjCoeffs[i].asDouble();
			if(abs(SaveObjCoeffs[i].asDouble())>MaxNormObj)
				MaxNormObj=abs(SaveObjCoeffs[i].asDouble());
		}
		NormC=sqrt(NormC);
		NormDelta=sqrt(NormDelta);
		Normalization=NormDelta/NormC;

		for (int i = 0; i < n; i++){
			//ObjCoeff.push_back(0);
			Fractionality.push_back(make_pair(i, 0));
			RandomOrder.push_back(make_pair(i,drand(RandSeed)));
		}


	}

//Random Number Generator
 	static inline double drand(double& seed) {
 		assert(seed>0);
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;        
		return seed / 2147483647; }

// Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline int irand(double& seed, int size) {
        return (int)(drand(seed) * size); }

// AlternativeRounding Helper (see Achterberg) 
        double rho(double z, double ObjCoeff){
        	double first=0;
        	if(z<=0.5)
        		first= 2*z*(1-z);
        	else
        		first= 1-2*z*(1-z);
        	//	if(MaxNormObj>0)
        	//	first=first+ObjSense*1/4*alpha*ObjCoeff/MaxNormObj;
        	if(first<0) return 0;
        	else if(first>1) return 1;
        	else return first;
        	
        }

        // AlternativeRounding Helper (see Achterberg)
        double rho2(double z, double IPsol, double LPsol){
        	        double X = 1 - fabs(IPsol-LPsol);
        	        if (z < X) {
        	        	return IPsol;
        	        } else {
        	        	return 1.0-IPsol;
						return floor(LPsol+rho(drand(RandSeed),1.0));
        	        	//return round(LPsol);//1.0-IPsol;
        	        }
                }

        double SumValues(vector< pair<int, double> > X){
    	double R=0;
    	for (int i=0;i<X.size();i++) R+=X[i].second;
    	return R;
    }
   
	bool FindIntegerSolution(std::vector<data::QpNum> &Output, double &IntScore, int *types, bool UseObjectiveDecay, bool AlternativeRounding, bool UsePropagation, int &Time, QBPSolver *qbp){
		SolutionFound=false;
		DeadEnds.clear();
		DeadAlpha.clear();
		DeadCount.clear();
		Output.clear();

        if (qbp->decisionLevel()>=1 && qbp->getBranchingVariable(qbp->decisionLevel()) >=0) {
		   if(qbp->getInfoLevel()>1) cerr<<"Start Pump in level " << qbp->decisionLevel() << " and block " << qbp->getBlock( qbp->getBranchingVariable(qbp->decisionLevel()) ) << endl;
        } else {
            if(qbp->getInfoLevel()>1) cerr<<"Start Pump in level " << qbp->decisionLevel() << " and level " << qbp->decisionLevel() << endl;
        }
		PreviousFrac.clear();
		
		if (AlternativeRounding) {
			if (PrevIPSol.size() > 0) AlternativeRounding = true;
			else AlternativeRounding = false;
		}
		// Set influence factory of objective function
		if(UseObjectiveDecay) alpha=1.0;//0.9;
		else alpha=0;

		bool NewX=false;
		int CountIterations=0;

		// for variables rounded to 1 (\tilde x_j=1) it is |x_j-\tilde x_j|=1-x_j; CountOnes stores the constants, i.e. the Ones
		int CountOnes=0;		
		
		// Variables for the regression line
		double D_v=0;
		double Da_v=0;
		double sumx=0;
		double sumy=0;
		double growth=0;

		std::vector<data::QpNum> PreviousX(n);
		std::vector<data::QpNum> ObjCoeff(n);		
		std::vector<data::QpNum> CurrentLPSol(n);	//LP Solution vector
		std::vector<data::QpNum> PreviousLPSol(n);	//LP Solution vector

		int CountTooLong=0;	//Counter for rounds that took too long (only a limited number allowed)
		int CountRES=0;		//Count Restarts; limited number allowed until restart
		int CountPerturbs=0;

		unsigned int StartingTime = time(NULL);

		extSol::QpExternSolver::QpExtSolBase rembase;
		rembase.variables.clear();
		rembase.constraints.clear();
		StSolver->getExternSolver(stage).getBase(rembase);

		bool PerturbIntSolution=false;

		//Variables needed for Propagation
		CRef confl = CRef_Undef;
		CRef confl_partner = CRef_Undef;
		int confl_var = -1;

		while(CountIterations<MaxIter && (Time==-1||(time(NULL)-StartingTime)<=MaxMultiTime*Time) && (time(NULL)-StartingTime)<MaxTime&& CountTooLong<MaxLong){	//Hinzufügen: || Anstieg der Gerade durch die PreviousFrac Werte ist noch negativ genug (Hoffnung besteht)
			CountIterations++;

			int IterTime=time(NULL);
			NewX=false;
			bool HardRestart=false;

			//Indicator whether propagation is still possible (becomes impossible, if contradiction reached)
			bool PropagateOK=true;
			int CountProp=0;

			int StartingTrailSize=qbp->getTrailSize();
			int StartingDecisionLevel=qbp->decisionLevel();

			if(qbp->getInfoLevel()>2) cerr<<"P";

			//cerr << CountIterations << " " << SumValues(Fractionality)<< " " << alpha <<endl;
			if(1||CountIterations>1){
			  //#ifdef USE_NBD_CLP
				if(1||CountIterations!=1){ 
					if (rembase.variables.size() != StSolver->getExternSolver(stage).getVariableCount()) {
						rembase.variables.resize(StSolver->getExternSolver(stage).getVariableCount());
					}
					StSolver->getExternSolver(stage).setBase(rembase);
				}
				//#endif

				//Solve the new LP
#ifdef EIGEN_LAZY
                std::vector< int > LPsnapshot_row_indices;
                HTCutentry *HTCe;
                pair<coef_t, uint64_t> hash;
                do {
                    //LPsnapshot_row_indices.clear();
                    for (int j = 0; j < /*QlpStSolve.getExternSolver(maxLPstage).getLProws_snapshot()*/LPsnapshot_row_indices.size();j++) {
                        int i = LPsnapshot_row_indices[j];
                        
                        hash = qbp->HTC->computeHash((*StSolver->getExternSolver(stage).getRowLhs_snapshot(i)),
                                                (*StSolver->getExternSolver(stage).getRowRhs_snapshot())[i].getValue().asDouble(),
                                                (*StSolver->getExternSolver(stage).getRowRhs_snapshot())[i].getRatioSign());
                        if (1||!HTC->getEntry(&HTCe,hash.second, hash.first)) {
                            qbp->listOfEnteredCuts.push( std::make_pair(StSolver->addUserCut(stage,
                                                (*StSolver->getExternSolver(stage).getRowLhs_snapshot(i)),
                                                (*StSolver->getExternSolver(stage).getRowRhs_snapshot())[i].getRatioSign(),
                                                (*StSolver->getExternSolver(stage).getRowRhs_snapshot())[i].getValue()), i) );
                            qbp->listOfEnteredCutHashs.push(hash);
                            //HTC->setEntry(hash.first, hash.second);
                        }
                    }
                    if (LPsnapshot_row_indices.size() > 0) cerr << "ta" << LPsnapshot_row_indices.size() <<";";
                    LPsnapshot_row_indices.clear();
                    StSolver->solveStage(stage, status, lb, ub, CurrentLPSol, algorithm::Algorithm::WORST_CASE, -1, -1);
                } while (status == algorithm::Algorithm::FEASIBLE && StSolver->getExternSolver( stage ).getLazyRows( LPsnapshot_row_indices, CurrentLPSol, 1e-12 ) );
#endif
		int m = StSolver->getExternSolver( stage ).getRowCount();
		if (m < rembase.constraints.size()) rembase.constraints.resize(m);
		else if (m > rembase.constraints.size()) {
			  for (int i = rembase.constraints.size(); i < m; i++)
			    rembase.constraints.push_back(extSol::QpExternSolver::NotABasicStatus);
		} 
		if (rembase.variables.size()>0) {
		   StSolver->getExternSolver( stage ).setBase(rembase);
		}

                //StSolver->solveStage(stage, status, lb, ub, CurrentLPSol, algorithm::Algorithm::WORST_CASE, -1, -1);
                qbp->QLPSTSOLVE_SOLVESTAGE(-1e100, stage, status, lb, ub, CurrentLPSol, algorithm::Algorithm::WORST_CASE,1, -1, CountIterations > 1 || SolutionFound ? -10 : -1, false);

		if (qbp->getInfoLevel()>-8) cerr << "CountIterations=" << CountIterations << " alpha=" << alpha << endl;
		//#ifdef USE_NBD_CLP
				rembase.variables.clear();
				rembase.constraints.clear();
				StSolver->getExternSolver(stage).getBase(rembase);
				//#endif
			}
			else{
				for (int i = 0; i < n; i++)
					CurrentLPSol[i]=OriginalLPSolution[i];
			}			

			if(CountIterations==1 || status == algorithm::Algorithm::FEASIBLE){
				if(Time!=-1 && CountIterations==1){
					// if it is not the first call of the Pump; but the first round of an improvement attempt
					
					//Either use a hardrestart
					//HardRestart=true;

					//Or: Perturb n/4 of the variables; promising ones first (with regard to thei objective coeff)
					PerturbIntSolution=true;
					NewX=false;
					goto Perturb;
				}	
					
				IsInteger=true;
				//Find out the new order; especially needed when propagation is used
				//if(Time!=-1) //in second rounds
					std::sort(RandomOrder.begin(), RandomOrder.end(), [](std::pair<int, double> p1, std::pair<int, double> p2){ return p1.second > p2.second; });

				for (int k=0;k<n;k++){
					RandomOrder[k].second=drand(RandSeed);	//Set new importance
					int i=RandomOrder[k].first;				//Select next variable
					ObjCoeff[i]=0;
					Fractionality[i].first=i;				//needs to be set again; was sorted

					//If the variable values are neither close to one nor close to zero
					if(types[i]==0){
						
						//get fractionality of this variable
						if (round(CurrentLPSol[i].asDouble())==0) Fractionality[i].second=CurrentLPSol[i].asDouble();
						else Fractionality[i].second=1-CurrentLPSol[i].asDouble();

						//Check whether this variable is already integer
						if (IsInteger && !(abs(CurrentLPSol[i].asDouble()-1)<epsNull ||abs(CurrentLPSol[i].asDouble())<epsNull))
							IsInteger=false;

						//Not sure if this is useful..
						if(!NewX && UseObjectiveDecay&& CountIterations!=1 && types[i]==0 && abs(PreviousLPSol[i].asDouble()-CurrentLPSol[i].asDouble())>epsNull) 
							NewX=true;
						PreviousLPSol[i]=CurrentLPSol[i];

						//If we use propagation and this variable is now already set in assigns
						if(UsePropagation&&qbp->getAssignment(i) != extbool_Undef){
							CountProp++;
							if(!NewX&&PreviousX[i]!=qbp->getAssignment(i)){
									NewX=true;
								}
							if(qbp->getAssignment(i)==0){
								ObjCoeff[i]=1;
								PreviousX[i]=0;
							}
							else{
								ObjCoeff[i]=-1;
								PreviousX[i]=1;

							}
						} 
						else{	//i.e. if variable is not set yet
							//if(drand(RandSeed)>0.5 && PrevIPSol.size()>0){
							if(AlternativeRounding){//drand(RandSeed)	>0.25){//drand(RandSeed)	>0.5)	{
								//double NewVal= floor(CurrentLPSol[i].asDouble()+rho(drand(RandSeed),SaveObjCoeffs[i].asDouble()));
								double NewVal= rho2(drand(RandSeed),PrevIPSol[i],CurrentLPSol[i].asDouble());
								if(!NewX&&PreviousX[i]!=NewVal){
										NewX=true;
									}
								if(NewVal==0){
									ObjCoeff[i]=1;
								}
								else{
									ObjCoeff[i]=-1;
								}
								PreviousX[i]=NewVal;
							}
							else{
								if (round(CurrentLPSol[i].asDouble())<0.5) {
									if(!NewX&&PreviousX[i]==1)
										NewX=true;
									ObjCoeff[i]=1;
									PreviousX[i]=0;
								}
								else if (round(CurrentLPSol[i].asDouble())>0.5){
									if(!NewX&&PreviousX[i]==0){
										NewX=true;
									}
									ObjCoeff[i]=-1;
									PreviousX[i]=1;
								} 
								else{
									if(1||qbp->getInfoLevel()>2) cerr <<"ERROR WHILE ROUNDING "<< round(CurrentLPSol[i].asDouble())<<endl;
									return false;
								}
							}

							if(UsePropagation && PropagateOK){
								//cerr<<"A";
								if(qbp->isFixed(i) && qbp->getFixed(i)!=((PreviousX[i].asDouble() > 0.5 ? 1 : 0))){
									if (/*NewX&&*/PreviousX[i].asDouble()>0.5) {
										ObjCoeff[i]=1;
										PreviousX[i]=0;
									}
									else {
										ObjCoeff[i]=-1;
										PreviousX[i]=1;
									} 
								}
								if (qbp->assign(i, ((PreviousX[i].asDouble() > 0.5 ? 1 : 0)), qbp->getTrailSize(), CRef_Undef, true, true) == ASSIGN_OK){ //confl? Or CRef_Undef
									qbp->increaseDecisionLevel();
									//cerr<<"-PR-";
									if(!qbp->propagate(confl, confl_var, confl_partner, false, false, 1000000)) {
										//cerr <<"-STOP! " << (double)qbp->trail.size()/(double)n<<endl;
										PropagateOK=false;
									}
								}
								else PropagateOK=false;

							}
						}
						if(UseObjectiveDecay){
						  double help=(/*1-*/alpha)*ObjCoeff[i].asDouble()+(1-alpha)*Normalization*SaveObjCoeffs[i].asDouble();
								ObjCoeff[i]=help;
							}
					}
					else{
						Fractionality[i].second=0;
						ObjCoeff[i]=alpha*Normalization*SaveObjCoeffs[i].asDouble();
					}
					StSolver->changeObjFuncCoeff(stage, i,ObjCoeff[i]);	
				}
				if(UsePropagation){
					//unassign and trail clearing
					while (qbp->getTrailSize() > StartingTrailSize) {
						qbp->unassign(qbp->getTrailsLast()/*trail.last()*/);
					}
					while (qbp->decisionLevel() > StartingDecisionLevel) {
						qbp->decreaseDecisionLevel();
						qbp->resolveFixed(qbp->decisionLevel(),false);
					}

					/*while (qbp->decisionLevel() > StartingDecisionLevel){
						int l=qbp->getTrailSize() - 1;
						int dl=qbp->decisionLevel() - 1;
						if (qbp->getTrailSize() >= 1) {
							while (qbp->getAssignmentLevel(qbp->getTrailElement(l)) > dl && l > 0) {
								qbp->unassign(qbp->getTrailsLast());
								l--;
							}
							qbp->unassign(qbp->getTrailsLast());
						}
						qbp->decreaseDecisionLevel();

					}*/
					assert(StartingTrailSize==qbp->getTrailSize()&& StartingDecisionLevel==qbp->decisionLevel());
				}

				alpha=alpha*Decay;
				if(alpha<epsNull) alpha=0;

				if(IsInteger && status==algorithm::Algorithm::FEASIBLE){
				        double oldIS = IntScore;
					IntScore=0;
					Output.clear();
					//assert(Solution.size()==0);
					for (int i=0;i<n;i++){
						IntScore+=CurrentLPSol[i].asDouble()*SaveObjCoeffs[i].asDouble();
						if(types[i]==0)	assert(abs(CurrentLPSol[i].asDouble()-round(CurrentLPSol[i].asDouble()))<=epsNull);
						Output.push_back(CurrentLPSol[i]); //Round first ?
					}						
					if(Time==-1)
						Time=time(NULL)-StartingTime; 
					if(Time==0) 
						Time=1;	
					SolutionFound=true;
					if(oldIS > IntScore && qbp->getShowInfo()) cerr <<endl <<"info: heuristic integer solution with objective value (integer Test) "<< IntScore << " in " << (time(NULL)-StartingTime) <<" second(s) and " << CountIterations << " LPCalls and " << CountPerturbs << " perturbations. " << SolutionFound << " " << Output.size() << endl;
					if(oldIS > IntScore) CountIterations = 2;
					if (1||oldIS > IntScore) alpha = alpha * 0.8; 
					if (alpha < 0.1) break;

					//cin.get();
					//data::QpNum ObjVal=data::QpNum(IntScore);
					//StSolver->tightenObjFuncBound(stage, ObjVal);
					////break;
				}
				
				PreviousFrac.push_back(SumValues(Fractionality));

				if (PreviousFrac.size() > 1 && PreviousFrac[PreviousFrac.size()-1] < PreviousFrac[PreviousFrac.size()-2]) CountIterations--;
				
				//If we have 200 previous rounds we calculate the regression line (least sum of squares)
				//If the growth is negative we still seem to find better solutions
				if(PreviousFrac.size()>=MaxStalls){
					D_v=0;
					Da_v=0;
					sumx=0;
					sumy=0;
					for (int i=0;i<PreviousFrac.size();i++){
						D_v+=i*i;
						Da_v+=PreviousFrac[i]*i;
						sumx+=i;
						sumy+=PreviousFrac[i];
						}
					growth=(PreviousFrac.size()*Da_v-sumx*sumy)/(double)(PreviousFrac.size()*D_v-sumx*sumx);
					
						//cerr <<"Anstieg "<< (double)growth << " " << PreviousFrac.size()<< endl;
					if(PreviousFrac.size()>MaxStalls){
					PreviousFrac.pop_front();
					}
				}
				//if fractionality does not get any better: Continue with normal FP
				//if(PreviousFrac.size()==MaxStalls&&PreviousFrac.back()/PreviousFrac.front()>=0.9){

				if(PreviousFrac.size()==MaxStalls && growth>-epsNull/*PreviousFrac.size()==MaxStalls&&PreviousFrac.back()/PreviousFrac.front()>=0.9*/){
					//cerr <<"FRAC ZERO "<< PreviousFrac.back() << " " <<PreviousFrac.front()<< endl;
					if(1||qbp->getInfoLevel()>2) cerr<<"F";
					//CountIterations=0;
					
				 	if(UseObjectiveDecay) alpha=alpha-0.1+drand(RandSeed)*0.2;
				 	HardRestart=true;
				 	NewX=false;
				 	PreviousFrac.clear();
					if (SolutionFound) break;
				}
				
				Perturb:
				//cerr << NewX << " " << OldDistance<<" " <<SolValue<<endl;
				if(!NewX /*&&LastObjValues.size()==MaxStalls &&(abs(LastObjValues.front()-LastObjValues.back())<=0.1 || *//*&&abs(OldDistance-SolValue)<=epsNull*/ ){
					//cerr <<"P "<<DeadEnds.size()<<endl;
					CountPerturbs++;
					if(qbp->getInfoLevel()>-8) cerr<<"p";

					bool OldEnd=false;
					int NumDead=0;
					for(int i=0;i<DeadEnds.size() && !OldEnd;i++){
						if(DeadEnds[i][0]==SumValues(Fractionality) && DeadAlpha[i]-alpha<=delta_alpha){
							for(int j=0;j<n;j++){
								if(DeadEnds[i][j+1]!=PreviousX[j].asDouble()){
									break;
								}
								if(j==n-1){
									DeadCount[i]++;
									NumDead=DeadCount[i];
									OldEnd=true;
								} 
							}
						}
					}
					if(OldEnd&&NumDead>=MaxDead){
						if(qbp->getInfoLevel()>-8) cerr<<endl<<"Found Old "<< alpha << " " << CountIterations<<endl;
						HardRestart=true;
						if(alpha!=0) alpha=alpha-0.1+drand(RandSeed)*0.2;//alpha=drand(RandSeed);
						if (alpha < 0.0) alpha = 0.0;
						if (alpha > 1.0) alpha = 1.0;
					}
					else if(!OldEnd){
						vector<double> Save(n+1,0);
						Save[0]=SumValues(Fractionality);
						for (int i=1;i<=n;i++)
							Save[i]=PreviousX[i-1].asDouble();
						DeadAlpha.push_back(alpha);
						DeadEnds.push_back(Save);
						DeadCount.push_back(1);
					}
					CountRES++;
					//if (getInfoLevel() >= 2)
					int TT=irand(RandSeed,T/2);
					if(CountRES%200==0){
						//cerr <<"Perturb harder"<<endl;
						HardRestart=true;
					}
					if(PerturbIntSolution) TT = n/4-T/2;
					if(HardRestart){
						CountRES=0;
						if(qbp->getInfoLevel()>-8) cerr<<"h";

				 		PreviousFrac.clear();
				 		if(UseObjectiveDecay) alpha=alpha-0.1+drand(RandSeed)*0.2;//alpha=drand(RandSeed);
						if (alpha > 1.0) alpha=1.0; 
						if (alpha < 0.0) alpha=0.0; 
						TT=n-T/2;
					} 
									
					//cerr <<"p";
					vector< pair<int, double> > ImprovementFinder;
					//cerr <<"Perturb "<< (double)TT/(double)n<<endl;

					if(!HardRestart){
						assert(Fractionality.size()==n);
						//cerr <<"Sort "<< Fractionality.size() << endl;
						std::sort(Fractionality.begin(), Fractionality.end(), [](std::pair<int, double> p1, std::pair<int, double> p2){ return p1.second > p2.second; });
						/*for (int i=0;i<Fractionality.size();i++)
							if(Fractionality[i].second>0)
							cerr << Fractionality[i].first << " " << Fractionality[i].second << endl;
					*/
					}
					if(PerturbIntSolution && Output.size()==n){
						for (int i = 0; i < n; i++){
							if((SaveObjCoeffs[i]>=0 && Output[i]==1)||(SaveObjCoeffs[i]<=0 && Output[i]==0))
								ImprovementFinder.push_back(make_pair(i, SaveObjCoeffs[i].asDouble()));
						}
						//Output.clear();
						sort(ImprovementFinder.begin(), ImprovementFinder.end(), [](pair<int, double> p1, pair<int, double> p2){ return p1.second > p2.second; });
					}

					//ZF Coeff Setzen: StSolver.changeObjFuncCoeff(stage, obj_lhs[i].index, obj_lhs[i].value);
					//		cerr << TT << endl;
					int c=0;
					int i=0;
					while(c<T/2+TT && i!=n){
						//cerr << i << " " <<Fractionality[i].second <<endl;
						if(types[Fractionality[i].first]==0){
							
							if(HardRestart){
								if(drand(RandSeed)<Fractionality[i].second+RestartP){
									if (round(CurrentLPSol[Fractionality[i].first].asDouble())==0){
										StSolver->changeObjFuncCoeff(stage, Fractionality[i].first,-1);
									}
									else{
										StSolver->changeObjFuncCoeff(stage, Fractionality[i].first,1);
									}
									c++;
								}
							}
							/*if(HardRestart){
								if(drand(RandSeed)>0.5){
									StSolver->changeObjFuncCoeff(stage, Fractionality[i].first,-1);
									}
								else StSolver->changeObjFuncCoeff(stage, Fractionality[i].first,1);
								c++;
							}*/
							else if(PerturbIntSolution){
								if(i==ImprovementFinder.size()) break;
								if(ImprovementFinder[i].second>0){
									StSolver->changeObjFuncCoeff(stage, ImprovementFinder[i].first,1);
									}
								else StSolver->changeObjFuncCoeff(stage, ImprovementFinder[i].first,-1);
								c++;
							}
							else{
								if(1||Fractionality[i].second>0.2){
							//cerr << "f";
								//if(CurrentLPSol[xAbs[i].first].asDouble()>0.2){
									if (round(CurrentLPSol[Fractionality[i].first].asDouble())==0){
										StSolver->changeObjFuncCoeff(stage, Fractionality[i].first,-1);
									}
									else{
										StSolver->changeObjFuncCoeff(stage, Fractionality[i].first,1);
									}
									c++;
								}
								else break;
							}
								PreviousX[i]=-1;//-PreviousX[i].asDouble();
						}
						else { 
							StSolver->changeObjFuncCoeff(stage, Fractionality[i].first,alpha*Normalization*SaveObjCoeffs[Fractionality[i].first].asDouble());
							c++;
						}
						i++;
						//}
					}
					PerturbIntSolution=false;
					HardRestart=false;

				}
			}
			/*else if(status==algorithm::Algorithm::IT_LIMIT){
				cerr <<"LIMIT"<<endl;
				for (int i=0;i<n;i++){
					Fractionality[i].first=i;

					ObjCoeff[i]=0;
					//If the variable values are neither close to one nor close to zero
					if(types[i]==0){
						//Check whether this variable is already integer
						if(CurrentLPSol[i].asDouble()<0.5){
							//cerr<<i << " " <<CurrentLPSol[i].asDouble() << endl;
							Fractionality[i].second=abs(CurrentLPSol[i].asDouble());
							ObjCoeff[i]=1;

						}
						else{
							Fractionality[i].second=abs(1-CurrentLPSol[i].asDouble());

							CountOnes++;
							ObjCoeff[i]=-1;
						}
						if(UseObjectiveDecay){
							double help=(1-alpha)*ObjCoeff[i].asDouble()+alpha*Normalization*SaveObjCoeffs[i].asDouble();
							ObjCoeff[i]=help;
						}
					}
					else{
						ObjCoeff[i]=alpha*Normalization*SaveObjCoeffs[i].asDouble();
					}
					StSolver->changeObjFuncCoeff(stage, i,ObjCoeff[i]);
				}
				//cin.get();	
				alpha=alpha*Decay;
				if(alpha<epsNull) alpha=0;alpha=alpha*Decay;
				if(alpha<epsNull) alpha=0;
			}*/
			else break;

			IterTime=time(NULL)-IterTime;
			if(IterTime>TimeLong) CountTooLong++;

				
		}
		//cin.get();
			//restore objective function in Solver 
		for (int i=0;i<n;i++){
			if( n ==StSolver->getExternSolver(stage).getVariableCount() )
				StSolver->changeObjFuncCoeff(stage, i, SaveObjCoeffs[i].asDouble());
			else{
				StSolver->changeObjFuncCoeff(stage, i, 0);
			}
		}
 		//---------FALSCH!
		if (Output.size()>0&&SolutionFound){
		  if(qbp->getInfoLevel()>-8) cerr << "Solution ok pump " << Output.size() << " " << SolutionFound << endl;
		  PrevIPSol.clear();
		  for (int i = 0; i < Output.size();i++) {
		    PrevIPSol.push_back(Output[i].asDouble());
		  }
		  return true;
		}
		else{
		  if(qbp->getShowInfo()) cerr << "Warning: Solution not ok pump " << Output.size() << " " << SolutionFound << endl;
		  if(qbp->getInfoLevel()>2&& CountTooLong==MaxLong) cerr <<"Too many long lasting rounds" <<endl;
		  if(qbp->getInfoLevel()>2&& CountIterations==MaxIter) cerr <<"Iterations limit reached "<<endl;
		  if(qbp->getInfoLevel()>2&&(time(NULL)-StartingTime)>MaxMultiTime*Time&& Time!=-1) cerr<< "Time limit reached" <<endl;
		  return false;
		}
	}
};


#endif /* FEASIBILITYPUMP_H_ */
