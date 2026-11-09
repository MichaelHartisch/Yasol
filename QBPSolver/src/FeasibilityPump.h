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
#ifdef USE_NBD_HIGHS
#include "ExternSolvers/QpExtSolHighs.hpp"
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
		Decay=0.999;//1.0;//0.999;

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
		Normalization = 1.0/(MaxNormObj+1.0);
 
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
   
bool FindIntegerSolution(
        std::vector<data::QpNum> &Output,
        double &IntScore,
        int *types,
        bool UseObjectiveDecay,
        bool AlternativeRounding,
        bool UsePropagation,
        int &Time,
	data::Qlp &qlp,
        QBPSolver *qbp)
{
    SolutionFound = false;

    DeadEnds.clear();
    DeadAlpha.clear();
    DeadCount.clear();
    Output.clear();

    if (qbp->decisionLevel() >= 1 &&
        qbp->getBranchingVariable(qbp->decisionLevel()) >= 0) {

        if (qbp->getInfoLevel() > 1)
            cerr << "Start Pump in level "
                 << qbp->decisionLevel()
                 << " and block "
                 << qbp->getBlock(
                        qbp->getBranchingVariable(
                            qbp->decisionLevel()))
                 << endl;
    }
    else {
        if (qbp->getInfoLevel() > 1)
            cerr << "Start Pump in level "
                 << qbp->decisionLevel()
                 << " and level "
                 << qbp->decisionLevel()
                 << endl;
    }

    PreviousFrac.clear();


    // ------------------------------------------------------------
    // Alternative rounding only makes sense if an old IP solution
    // is available.
    // ------------------------------------------------------------

    if (AlternativeRounding) {
        if (PrevIPSol.size() == 0)
            AlternativeRounding = false;
    }


    // ------------------------------------------------------------
    // Objective influence
    // ------------------------------------------------------------

    if (UseObjectiveDecay)
        alpha = 1.0;
    else
        alpha = 0.0;


    bool NewX = false;

    // NOTE:
    // This is the OLD algorithmic counter. It is intentionally kept,
    // even though it is not monotone.
    int CountIterations = 0;

    int CountOnes = 0;


    // ------------------------------------------------------------
    // Variables for regression on fractionality
    // ------------------------------------------------------------

    double D_v = 0;
    double Da_v = 0;
    double sumx = 0;
    double sumy = 0;
    double growth = 0;


    std::vector<data::QpNum> PreviousX(n);
    std::vector<data::QpNum> ObjCoeff(n);
    std::vector<data::QpNum> CurrentLPSol(n);
    std::vector<data::QpNum> PreviousLPSol(n);


    int CountTooLong = 0;
    int CountRES = 0;
    int CountPerturbs = 0;
    int NumRoundingChanges = 0;

    // ============================================================
    // NEW: independent effort accounting
    // ============================================================

    long long NumPumpLPCalls = 0;
    long long PumpLPIterations = 0;
    double PumpLPTime = 0.0;

    const long long MaxPumpLPCalls = 100;
    const long long MaxPumpLPIterations = 50000;
    int RoundingChangesForNextLP = 0;
    bool EffortLimitReached = false;
    std::vector<double> LastBuiltPumpObj;
    
    double AlphaUsedForBuild = alpha;
    double LastBuiltAlpha = alpha;
    
    int BuiltObjChanges = 0;
    int BuiltSignChanges = 0;
    int BuiltMagnitudeChanges = 0;
    
    double BuiltObjL1Change = 0.0;
    double BuiltObjMaxChange = 0.0;
    
 // ------------------------------------------------------------
    // Direct access to HiGHS
    // ------------------------------------------------------------

    //Highs *M =
    //    (Highs*)StSolver->getExternSolver(stage).getSolverEnv();


    // ============================================================
    // NEW:
    // Keep our OWN image of the objective.
    //
    // This avoids comparing against M->getLp().col_cost_ after
    // QLPSTSOLVE_SOLVESTAGE, since the solver machinery may modify
    // / transform the objective internally.
    // ============================================================

    std::vector<double> CurrentPumpObj(n, 0.0);
    std::vector<double> PreviousSolvedPumpObj;

    const std::vector<data::QpNum>& SaveObjCoeffsTmp = qlp.getObjectiveFunctionValues();
    if (SaveObjCoeffsTmp.size() >= n) {
        for (int i = 0; i < n; ++i)
	  CurrentPumpObj[i] = SaveObjCoeffsTmp[i].asDouble();
    }
    //if (M != NULL && M->getNumCol() >= (int)n) {
    //    const std::vector<double> &cost = M->getLp().col_cost_;
    //    for (int i = 0; i < n; ++i)
    //        CurrentPumpObj[i] = cost[i];
    //}
    else {
        for (int i = 0; i < n; ++i)
            CurrentPumpObj[i] = SaveObjCoeffs[i].asDouble();
    }

    // ------------------------------------------------------------
    // IMPORTANT:
    // From now on, every objective modification inside this routine
    // should go through this lambda.
    // ------------------------------------------------------------

    auto SetPumpObjCoeff = [&](int col, double value)
        {
            StSolver->changeObjFuncCoeff(stage, col, value);
            if (col >= 0 && col < (int)CurrentPumpObj.size())
                CurrentPumpObj[col] = value;
        };

    unsigned int StartingTime = time(NULL);

    // ------------------------------------------------------------
    // Save basis
    // ------------------------------------------------------------

    extSol::QpExternSolver::QpExtSolBase rembase;

    rembase.variables.clear();
    rembase.constraints.clear();

    StSolver->getExternSolver(stage).getBase(rembase);

    bool PerturbIntSolution = false;

    // ------------------------------------------------------------
    // Variables needed for propagation
    // ------------------------------------------------------------

    CRef confl = CRef_Undef;
    CRef confl_partner = CRef_Undef;
    int confl_var = -1;

    // ============================================================
    // MAIN FP LOOP
    // ============================================================

    while (
        CountIterations < MaxIter &&
        (Time == -1 ||
         (time(NULL) - StartingTime) <= MaxMultiTime * Time) &&
        (time(NULL) - StartingTime) < MaxTime &&
        CountTooLong < MaxLong)
    {
        CountIterations++;
	NumRoundingChanges = 0;
	
        int IterTime = time(NULL);

        NewX = false;
        bool HardRestart = false;

        bool PropagateOK = true;
        int CountProp = 0;

        int StartingTrailSize = qbp->getTrailSize();

        int StartingDecisionLevel = qbp->decisionLevel();

        if (qbp->getInfoLevel() > 2) cerr << "P";

        // ========================================================
        // Restore previous basis
        // ========================================================

        if (rembase.variables.size() != StSolver->getExternSolver(stage).getVariableCount()) {
            rembase.variables.resize( StSolver->getExternSolver(stage).getVariableCount());
        }
        StSolver->getExternSolver(stage).setBase(rembase);

        // ========================================================
        // Adjust stored basis to current row count
        // ========================================================

        int m = StSolver->getExternSolver(stage).getRowCount();

        if (m < rembase.constraints.size()) {
            rembase.constraints.resize(m);
        } else if (m > rembase.constraints.size()) {
            for (int i = rembase.constraints.size(); i < m; i++) {
	      rembase.constraints.push_back( extSol::QpExternSolver::Basic/*NotABasicStatus*/);
            }
        }

        if (rembase.variables.size() > 0)
            StSolver->getExternSolver(stage).setBase(rembase);

        // ========================================================
        // NEW: Determine how much the objective changed since the
        // previous ACTUALLY SOLVED pump LP.
        // ========================================================

        int ChangedObj = 0;
        if (!PreviousSolvedPumpObj.empty()) {
            int nn = std::min( (int)CurrentPumpObj.size(), (int)PreviousSolvedPumpObj.size());
            for (int i = 0; i < nn; ++i) {
                if (abs( CurrentPumpObj[i] - PreviousSolvedPumpObj[i]) > 1e-9) {
                    ChangedObj++;
                }
            }
        }

        // This is exactly the objective about to be solved.
        PreviousSolvedPumpObj = CurrentPumpObj;

        // ========================================================
        // NEW: time the actual LP call
        // ========================================================

        auto lpStart = std::chrono::steady_clock::now();


        qbp->QLPSTSOLVE_SOLVESTAGE( -1e100, stage, status, lb, ub, CurrentLPSol, algorithm::Algorithm::WORST_CASE,
            1, -1, CountIterations > 1 || SolutionFound ? -10 : -1, false);

        auto lpEnd = std::chrono::steady_clock::now();

        double ThisLPTime = std::chrono::duration<double>( lpEnd - lpStart).count();

        int ThisLPIterations = 0;

        //if (M != NULL) {
        //    ThisLPIterations = M->getInfo().simplex_iteration_count;
        //}

        // ========================================================
        // NEW: monotone counters
        // ========================================================

        NumPumpLPCalls++;

        PumpLPIterations += (long long)ThisLPIterations;

        PumpLPTime += ThisLPTime;

	if (0) cerr
	  << "FP-LP #" << NumPumpLPCalls
	  << " roundingChanges=" << RoundingChangesForNextLP
	  << " changedObj=" << ChangedObj
	  << " simplexIters=" << ThisLPIterations
	  << " time=" << ThisLPTime
	  << " sec"
	  << endl;
 
        // ========================================================
        // NEW: hard effort guard.
        //
        // Do NOT break immediately: first inspect this LP solution.
        // ========================================================

        if (NumPumpLPCalls >= MaxPumpLPCalls || PumpLPIterations >= MaxPumpLPIterations) {
            EffortLimitReached = true;
        }

        if (qbp->getInfoLevel() > -8)
            cerr << "CountIterations=" << CountIterations << " alpha=" << alpha << endl;

        // ========================================================
        // Save new basis
        // ========================================================

        rembase.variables.clear();
        rembase.constraints.clear();

        StSolver->getExternSolver(stage).getBase(rembase);

        // ========================================================
        // Process LP solution
        // ========================================================

        if (CountIterations == 1 || status == algorithm::Algorithm::FEASIBLE)
        {
            if (Time != -1 && CountIterations == 1)
            {
                PerturbIntSolution = true;
                NewX = false;
                goto Perturb;
            }

            IsInteger = true;

            // ----------------------------------------------------
            // New random order
            // ----------------------------------------------------

            std::sort( RandomOrder.begin(), RandomOrder.end(), [](std::pair<int,double> p1, std::pair<int,double> p2)
                {
                    return p1.second > p2.second;
                });

            // ====================================================
            // Build rounded solution and new FP objective
            // ====================================================
	    // Dieses alpha wird beim Aufbau der JETZT erzeugten
	    // Zielfunktion tatsächlich benutzt.
	    AlphaUsedForBuild = alpha;
	    
            for (int k = 0; k < n; k++)
            {
                RandomOrder[k].second = drand(RandSeed);

                int i = RandomOrder[k].first;

                ObjCoeff[i] = 0;

                Fractionality[i].first = i;

                if (types[i] == 0)
                {
                    // --------------------------------------------
                    // Fractionality
                    // --------------------------------------------
                    if (round( CurrentLPSol[i].asDouble()) == 0) {
                        Fractionality[i].second = CurrentLPSol[i].asDouble();
                    } else {
                        Fractionality[i].second = 1 - CurrentLPSol[i].asDouble();
                    }

                    if (IsInteger && !(abs( CurrentLPSol[i].asDouble() - 1) < epsNull ||
                          abs( CurrentLPSol[i].asDouble()) < epsNull)) {
                        IsInteger = false;
                    }

                    if (!NewX && UseObjectiveDecay && CountIterations != 1 &&
                        abs(PreviousLPSol[i].asDouble() - CurrentLPSol[i].asDouble()) > epsNull) {
                        NewX = true;
                    }

                    PreviousLPSol[i] = CurrentLPSol[i];

                    // --------------------------------------------
                    // Already assigned by propagation
                    // --------------------------------------------

                    if (UsePropagation && qbp->getAssignment(i) != extbool_Undef) {
                        CountProp++;

                        if (!NewX && PreviousX[i] != qbp->getAssignment(i)) {
                            NewX = true;
                        }

                        if (qbp->getAssignment(i) == 0) {
                            ObjCoeff[i] = 1;
                            PreviousX[i] = 0;
                        } else {
                            ObjCoeff[i] = -1;
                            PreviousX[i] = 1;
                        }
                    }

                    // --------------------------------------------
                    // Variable not assigned
                    // --------------------------------------------

                    else {
                        if (AlternativeRounding) {
                            double NewVal = rho2(drand(RandSeed), PrevIPSol[i], CurrentLPSol[i].asDouble());

                            if (!NewX && PreviousX[i] != NewVal) {
                                NewX = true;
                            }


                            if (NewVal == 0)
                                ObjCoeff[i] = 1;
                            else
                                ObjCoeff[i] = -1;

                            PreviousX[i] = NewVal;
                        } else {

			  if (round(CurrentLPSol[i].asDouble()) < 0.5) {

			    if (PreviousX[i].asDouble() == 1.0)
			      NumRoundingChanges++;

			    if (!NewX && PreviousX[i] == 1)
			      NewX = true;

			    ObjCoeff[i] = 1;
			    PreviousX[i] = 0;
			  }
			  else if (round(CurrentLPSol[i].asDouble()) > 0.5) {

			    if (PreviousX[i].asDouble() == 0.0)
			      NumRoundingChanges++;

			    if (!NewX && PreviousX[i] == 0)
			      NewX = true;

			    ObjCoeff[i] = -1;
			    PreviousX[i] = 1;
			  } else {
                                cerr
                                    << "ERROR WHILE ROUNDING "
                                    << round(
                                        CurrentLPSol[i].
                                            asDouble())
                                    << endl;

                                return false;
                            }
                        }

                        // ----------------------------------------
                        // Propagation
                        // ----------------------------------------

                        if (UsePropagation && PropagateOK) {
                            if (qbp->isFixed(i) && qbp->getFixed(i) != (PreviousX[i].asDouble() > 0.5 ? 1 : 0)) {
                                if (PreviousX[i].asDouble() > 0.5)
                                {
                                    ObjCoeff[i] = 1;
                                    PreviousX[i] = 0;
                                }
                                else
                                {
                                    ObjCoeff[i] = -1;
                                    PreviousX[i] = 1;
                                }
                            }

                            if (qbp->assign(i, PreviousX[i].asDouble() > 0.5 ? 1 : 0,
                                    qbp->getTrailSize(), CRef_Undef, true, true) == ASSIGN_OK) {
                                qbp->increaseDecisionLevel();

                                if (!qbp->propagate( confl, confl_var, confl_partner, false, false, 1000000)) {
                                    PropagateOK = false;
                                }
                            }
                            else {
                                PropagateOK = false;
                            }
                        }
                    }

                    // --------------------------------------------
                    // Original-objective mixing
                    // ORIGINAL FORMULA LEFT UNCHANGED
                    // --------------------------------------------

                    if (UseObjectiveDecay) {
                        double help = alpha * ObjCoeff[i].asDouble() +
                            (1 - alpha) * Normalization * SaveObjCoeffs[i].asDouble();
                        ObjCoeff[i] = help;
                    }
                } else {
                    Fractionality[i].second = 0;

                    ObjCoeff[i] = (1.0-alpha) /*alpha*/ * Normalization * SaveObjCoeffs[i].asDouble();
                }

                // IMPORTANT:
                // all objective changes go through our tracker
                SetPumpObjCoeff( i, ObjCoeff[i].asDouble());
            }

	    RoundingChangesForNextLP = NumRoundingChanges;

	    BuiltObjChanges = 0;
	    BuiltSignChanges = 0;
	    BuiltMagnitudeChanges = 0;

	    BuiltObjL1Change = 0.0;
	    BuiltObjMaxChange = 0.0;

	    if (!LastBuiltPumpObj.empty()) {

	      for (int i = 0; i < n; ++i) {

		double oldc = LastBuiltPumpObj[i];
		double newc = CurrentPumpObj[i];

		double diff = newc - oldc;
		double adiff = abs(diff);

		if (adiff > 1e-9) {

		  BuiltObjChanges++;

		  BuiltObjL1Change += adiff;

		  if (adiff > BuiltObjMaxChange)
		    BuiltObjMaxChange = adiff;

		  // Hat das Vorzeichen wirklich gewechselt?
		  if ((oldc < -1e-9 && newc > 1e-9) ||
		      (oldc >  1e-9 && newc < -1e-9)) {
		    BuiltSignChanges++;
		  } else {
		    // Koeffizient hat nur seinen Betrag geändert
		    BuiltMagnitudeChanges++;
		  }
		}
	      }
	    }

	    if (0)cerr
		    << "FP-OBJ:"
		    << " nextLP=" << (NumPumpLPCalls + 1)
		    << " alpha=" << AlphaUsedForBuild
		    << " dAlpha="
		    << (LastBuiltPumpObj.empty()
			? 0.0
			: AlphaUsedForBuild - LastBuiltAlpha)
		    << " roundingChanges=" << NumRoundingChanges
		    << " objChanges=" << BuiltObjChanges
		    << " signChanges=" << BuiltSignChanges
		    << " magnitudeChanges=" << BuiltMagnitudeChanges
		    << " L1change=" << BuiltObjL1Change
		    << " maxChange=" << BuiltObjMaxChange
		    << endl;


	    LastBuiltPumpObj = CurrentPumpObj;
	    LastBuiltAlpha = AlphaUsedForBuild;
 
            // ====================================================
            // Undo temporary propagation assignments
            // ====================================================

            if (UsePropagation) {
                while ( qbp->getTrailSize() > StartingTrailSize) {
                    qbp->unassign(qbp->getTrailsLast());
                }

                while (qbp->decisionLevel() > StartingDecisionLevel) {
                    qbp->decreaseDecisionLevel();

                    qbp->resolveFixed( qbp->decisionLevel(), false);
                }

                assert(StartingTrailSize == qbp->getTrailSize() &&
                    StartingDecisionLevel == qbp->decisionLevel());
            }

            alpha = alpha * Decay;

            if (alpha < epsNull)
                alpha = 0;

            // ====================================================
            // Integer solution found
            // ====================================================

            if (IsInteger &&
                status == algorithm::Algorithm::FEASIBLE) {
                double oldIS = IntScore;

                IntScore = 0;
                Output.clear();

                for (int i = 0; i < n; i++)
                {
                    IntScore += CurrentLPSol[i].asDouble() * SaveObjCoeffs[i].asDouble();

                    if (types[i] == 0) {
                        assert( abs( CurrentLPSol[i].asDouble() - round( CurrentLPSol[i].asDouble())) <= epsNull);
                    }

                    Output.push_back(CurrentLPSol[i]);
                }


                if (Time == -1)
                    Time = time(NULL) - StartingTime;

                if (Time == 0)
                    Time = 1;


                SolutionFound = true;


                if (oldIS > IntScore && qbp->getShowInfo()) {
                    cerr
                        << endl
                        << "info: heuristic integer solution "
                        << "with objective value "
                        << IntScore
                        << " in "
                        << (time(NULL) -
                            StartingTime)
                        << " second(s), "
                        << NumPumpLPCalls
                        << " REAL LP calls, "
                        << PumpLPIterations
                        << " simplex iterations and "
                        << CountPerturbs
                        << " perturbations."
                        << endl;
                }


                // ORIGINAL algorithmic behaviour
                if (oldIS > IntScore)
                    CountIterations = 2;


                if (oldIS <= IntScore + fabs(IntScore * 0.05)) {
		  alpha = alpha * 0.8;
		  /*
		  if (alpha > 0.64 + 1e-9) {
		    alpha = alpha * 0.8;
		    
		    if (alpha < 0.64 + 1e-9)
		      alpha = 0.64;
		    alpha = 0.8; // LOOKHERE
		  }
		  */
                }

                if (alpha < 0.1)
                    break;
            }

            // ====================================================
            // Fractionality history
            // ====================================================

            PreviousFrac.push_back( SumValues(Fractionality));


            // IMPORTANT:
            // This is precisely why CountIterations is NOT a real
            // LP-call counter.
            if (PreviousFrac.size() > 1 && PreviousFrac[ PreviousFrac.size() - 1] < PreviousFrac[ PreviousFrac.size() - 2]) {
                CountIterations--;
            }


            // ====================================================
            // Regression over recent fractionalities
            // ====================================================

            if (PreviousFrac.size() >=  MaxStalls) {
                D_v = 0;
                Da_v = 0;
                sumx = 0;
                sumy = 0;

                for (int i = 0; i < PreviousFrac.size(); i++) {
                    D_v += i * i;
                    Da_v += PreviousFrac[i] * i;
                    sumx += i;
                    sumy += PreviousFrac[i];
                }

                growth = (PreviousFrac.size() * Da_v - sumx * sumy) /
                    (double)(PreviousFrac.size() * D_v - sumx * sumx);

                if (PreviousFrac.size() > MaxStalls) {
                    PreviousFrac.pop_front();
                }
            }

            if (PreviousFrac.size() == MaxStalls && growth > -epsNull) {
                if ( qbp->getInfoLevel() > 2)
                    cerr << "F";

                if (UseObjectiveDecay)
                {
                    alpha = alpha - 0.1 + drand(RandSeed) * 0.2;
                }

                HardRestart = true;
                NewX = false;

                PreviousFrac.clear();

                if (SolutionFound)
                    break;
            }


            // ====================================================
            // Perturbation
            // ====================================================

Perturb:

            if (!NewX)
            {
                CountPerturbs++;

                if (qbp->getInfoLevel() > -8)
                    cerr << "p";

                bool OldEnd = false;
                int NumDead = 0;

                // ------------------------------------------------
                // Search memorised dead ends
                // ------------------------------------------------

                for (int i = 0; i < DeadEnds.size() && !OldEnd; i++) {
                    if (DeadEnds[i][0] ==  SumValues( Fractionality) &&
                        DeadAlpha[i] - alpha <= delta_alpha) {
                        for (int j = 0; j < n; j++) {
                            if (DeadEnds[i][j + 1] != PreviousX[j].asDouble()) {
                                break;
                            }

                            if (j == n - 1) {
                                DeadCount[i]++;
                                NumDead = DeadCount[i];
                                OldEnd = true;
                            }
                        }
                    }
                }


                if (0&&OldEnd) {
                    cerr
                        << endl
                        << "FP: old dead-end found"
                        << ", DeadEnds.size()="
                        << DeadEnds.size()
                        << ", repeat count="
                        << NumDead
                        << ", CountPerturbs="
                        << CountPerturbs
                        << ", CountIterations="
                        << CountIterations
                        << endl;
                }


                if (OldEnd && NumDead >= MaxDead) {
                    if (qbp->getInfoLevel() > -8) {
                        cerr
                            << "Found Old "
                            << alpha
                            << " "
                            << CountIterations
                            << endl;
                    }

                    HardRestart = true;

                    if (alpha != 0) {
                        alpha = alpha - 0.1 + drand(RandSeed) * 0.2;
                    }

                    if (alpha < 0.0)
                        alpha = 0.0;

                    if (alpha > 1.0)
                        alpha = 1.0;
                }

                else if (!OldEnd) {
                    vector<double> Save( n + 1, 0);

                    Save[0] = SumValues( Fractionality);

                    for (int i = 1; i <= n; i++) {
                        Save[i] = PreviousX[i - 1].asDouble();
                    }

                    DeadAlpha.push_back(alpha);
                    DeadEnds.push_back(Save);
                    DeadCount.push_back(1);

                    if (0) cerr
                        << endl
                        << "FP: new dead-end"
                        << ", DeadEnds.size()="
                        << DeadEnds.size()
                        << ", CountPerturbs="
                        << CountPerturbs
                        << ", CountIterations="
                        << CountIterations
                        << endl;
                }

                CountRES++;

                int TT = irand( RandSeed, T / 2);

                if (CountRES % 200 == 0) {
                    HardRestart = true;
                }

                if (PerturbIntSolution) {
                    TT = n / 4 - T / 2;
                }

                if (HardRestart) {
                    CountRES = 0;
                    if (qbp->getInfoLevel() > -8) {
                        cerr << "h";
                    }

                    PreviousFrac.clear();

                    if (UseObjectiveDecay)
                    {
                        alpha = alpha - 0.1 + drand(RandSeed) * 0.2;
                    }

                    if (alpha > 1.0)
                        alpha = 1.0;

                    if (alpha < 0.0)
                        alpha = 0.0;

                    TT = n - T / 2;
                }

                vector<pair<int,double> > ImprovementFinder;

                if (!HardRestart) {
                    assert( Fractionality.size() == n);

                    std::sort( Fractionality.begin(), Fractionality.end(), [](std::pair<int,double> p1, std::pair<int,double> p2)
                        {
                            return p1.second > p2.second;
                        });
                }


                if (PerturbIntSolution && Output.size() == n)
                {
                    for (int i = 0; i < n; i++) {
                        if (
                            (SaveObjCoeffs[i] >= 0 && Output[i] == 1)
                            ||
                            (SaveObjCoeffs[i] <= 0 && Output[i] == 0))
                        {
                            ImprovementFinder.push_back( make_pair(i, SaveObjCoeffs[i].asDouble()));
                        }
                    }

                    sort(ImprovementFinder.begin(), ImprovementFinder.end(), [](pair<int,double> p1, pair<int,double> p2)
                        {
                            return p1.second > p2.second;
                        });
                }

                int c = 0;
                int i = 0;

                while ( c < T / 2 + TT && i != n) {
                    if (types[ Fractionality[i].first] == 0) {
                        if (HardRestart) {
                            if (drand(RandSeed) < Fractionality[i].second + RestartP) {
                                if (round( CurrentLPSol[ Fractionality[i].first].asDouble()) == 0) {
                                    SetPumpObjCoeff(Fractionality[i].first, -1);
                                } else {
                                    SetPumpObjCoeff( Fractionality[i].first, 1);
                                }

                                c++;
                            }
                        } else if (
                            PerturbIntSolution) {
                            if (ImprovementFinder.size() == 0)
                                break;

                            if (i == ImprovementFinder.size())
                                break;

                            if (ImprovementFinder[i].second > 0) {
                                SetPumpObjCoeff( ImprovementFinder[i].first, 1);
                            } else {
                                SetPumpObjCoeff(ImprovementFinder[i].first, -1);
                            }
                            c++;
                        } else {
                            if (1 || Fractionality[i].second > 0.2) {
                                if (round( CurrentLPSol[ Fractionality[i].first].asDouble()) == 0) {
                                    SetPumpObjCoeff( Fractionality[i].first, -1);
                                } else {
                                    SetPumpObjCoeff( Fractionality[i].first, 1);
                                }
                                c++;
                            } else {
                                break;
                            }
                        }

                        // ORIGINAL behaviour intentionally retained.
                        PreviousX[Fractionality[i].first] = -1;//PreviousX[i] = -1;
                    } else {
                        SetPumpObjCoeff( Fractionality[i].first, alpha * Normalization * SaveObjCoeffs[ Fractionality[i].first].asDouble());
                        c++;
                    }
                    i++;
                }

                PerturbIntSolution = false;
                HardRestart = false;
            }
        } else {
            break;
        }

        // ========================================================
        // OLD duration statistic
        // ========================================================

        IterTime = time(NULL) - IterTime;

        if (IterTime > TimeLong)
            CountTooLong++;

        // ========================================================
        // NEW:
        // Now that the current LP solution has been inspected,
        // obey the independent effort limit.
        // ========================================================

        if (EffortLimitReached)
        {
	  if (0)cerr
                << endl
                << "FP: effort limit reached"
                << " LPcalls="
                << NumPumpLPCalls
                << " simplexIters="
                << PumpLPIterations
                << " LPtime="
                << PumpLPTime
                << " sec"
                << endl;

            break;
        }
    }

    // ============================================================
    // NEW SUMMARY
    // ============================================================

    if (0)cerr
        << endl
        << "FP SUMMARY:"
        << " realLPcalls="
        << NumPumpLPCalls
        << " simplexIters="
        << PumpLPIterations
        << " avgIters="
        << (NumPumpLPCalls > 0
            ? (double)PumpLPIterations /
                (double)NumPumpLPCalls
            : 0.0)
        << " LPtime="
        << PumpLPTime
        << " avgLPtime="
        << (NumPumpLPCalls > 0
            ? PumpLPTime /
                (double)NumPumpLPCalls
            : 0.0)
        << " CountIterations="
        << CountIterations
        << " perturbations="
        << CountPerturbs
        << " deadEnds="
        << DeadEnds.size()
        << " solutionFound="
        << SolutionFound
        << endl;


    // ============================================================
    // Restore original objective
    // ============================================================

    for (int i = 0; i < n; i++) {
        if (n == StSolver->getExternSolver(stage).getVariableCount()) {
            StSolver->changeObjFuncCoeff( stage, i, SaveObjCoeffs[i].asDouble());
        } else {
            StSolver->changeObjFuncCoeff( stage, i, 0);
        }
    }

    // ============================================================
    // Return
    // ============================================================

    if (Output.size() > 0 && SolutionFound)
    {
        if (qbp->getInfoLevel() > -8)
        {
            cerr
                << "Solution ok pump "
                << Output.size()
                << " "
                << SolutionFound
                << endl;
        }

        PrevIPSol.clear();

        for (int i = 0; i < Output.size(); i++) {
            PrevIPSol.push_back( Output[i].asDouble());
        }
        return true;
    } else {
        if (qbp->getShowInfo())
        {
            cerr
                << "Warning: Solution not ok pump "
                << Output.size()
                << " "
                << SolutionFound
                << endl;
        }


        if (qbp->getShowInfo() && EffortLimitReached)
        {
            cerr
                << "FP stopped by independent effort limit."
                << endl;
        }

        else if (qbp->getShowInfo() && qbp->getInfoLevel() > 2 &&
                 CountTooLong == MaxLong)
        {
            cerr
                << "Too many long lasting rounds"
                << endl;
        }

        else if (qbp->getShowInfo() && qbp->getInfoLevel() > 2 &&
                 CountIterations == MaxIter)
        {
            cerr
                << "Iterations limit reached"
                << endl;
        }

        else if (qbp->getShowInfo() && 
            qbp->getInfoLevel() > 2 &&
            Time != -1 &&
            (time(NULL) - StartingTime) >
                MaxMultiTime * Time)
        {
            cerr
                << "Time limit reached"
                << endl;
        }


        return false;
    }
}
};

#endif /* FEASIBILITYPUMP_H_ */
