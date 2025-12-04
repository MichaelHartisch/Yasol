/*
 *
 * Yasol: Yasol.cpp -- Copyright (c) 2012-2017 Ulf Lorenz
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
#include <cstdlib>
using namespace std;
#include "yInterface.h"
#include "../../mpiClass/nomp.h"


void printUsage(const std::string& programName) {
    std::cout << "Usage: " << programName << " <instanceName> [OPTIONS]\n \n"
              << "instanceName     - first argument must be path to instance; must be in QLP file format \n"
              << "\n========For Users========================================================================\n"
              << " --help|-h              - show this help message\n"
              << " --timeLimit=i          - integer to indicate time limit t in seconds\n"
              << " --showSolution=b       - binary to turn off/on printing solution to console (default 0)\n"
              << " --outputFile=b         - binary to turn off/on the solver writing an output file to same\n"
              << "                          location as the instance file with ending '.solX' (default 0)\n"
              << " --useGMI=i             - integer i to - deactivate generation of GMI cuts (0) (default)\n"
              << "                                       - activate GMI at root (1)\n"
              << "                                       - activate GMI at inner nodes (2)\n"
              << "                                       - activate GMI at both root and inner nodes (3)\n"
              << " --useCover=b           - binary to turn off/on cover cuts (default 1)\n"
              << " --usePump=b            - binary to turn off/on feasibility pump (default 1)\n"
              << " --useMonotones=b       - binary to turn off/on the exploitation of monotes (default 1)\n"
              << " --isSimplyRestricted=b - binary to turn off/on the exploitation of a simple structure\n"
              << "                          of the universal constraints (the uncertainty set) (default 0)\n"
              << " --useCglRootCuts=b     - binary to turn off/on cut generation from CGL library (default 0)\n"
              << "\n========For Developers===================================================================\n"
              << " --showInfo=b               - binary to turn off/on info from solution process (default 0) \n"
              << " --showWarning=b            - binary to turn off/on warning from solution process (default 0)\n"
              << " --showError=b              - binary to turn off/on errors from solution process (default 0)\n"
              << " --maintainPv=b             - binary to turn off/on tracking of the value defining path (default 1)\n"
              << " --useLPCuts=b              - binary to turn off/on lift-and-project-cuts (default 0)\n"
              << " --useLazyLP=b              - binary to turn off/on lazy handling of lp constraints (default 1)\n"
              << " --useShadow=b              - binary to turn off/on exploitation of scenarios in  (default 1)\n"
              << " --useLimitedLP=b           - binary to turn off/on reduced LP calls at inner nodes (default 1)\n"
              << " --useAlphabeta=b           - binary to turn off/on alpha beta as main search algorithm (default 1)\n"
              << " --useStrongBranching=b     - binary to turn off/on strong branching (default 1)\n"
              << " --reduceStrongBranching=i  - integer to indicate reduction of usage of strong branching (default 2)\n"
              << " --useConflictGraph=i       - integer to indicate exploitation level of a conflict graph (default 2)\n";
}

int main(int argc, char** argv) {

  int myrank=0;
  MPI_Status status;

  if( !(argc-1) ){
    std::cout << "Usage error: provide an instance file or run with '-h' for help." << std::endl;
    return 0;
  }

  if (argc == 2 && (std::string(argv[1]) == "--h" || std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
      printUsage(argv[0]);
      return 0;
  }


  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  yInterface* yasolPt=new yInterface();
  yInterface& yasol = *yasolPt;
  coef_t result;
  time_t ini_time;
  int timeLimit=-1;

  while(argc > 0 && argv[argc-1][0]=='-' && argv[argc-1][1]=='-') {
    char buf[100];
    strcpy(buf,&(argv[argc-1][2]));
    if (strncmp(buf,"timeLimit=",10)==0) {
      timeLimit = std::stoi(std::string(buf + 10));
    }
    else if (strncmp(buf,"showSolution=",13)==0) {
      yasol.setShowSolution(std::stoi(std::string(buf + 13)));
    }
    else if (strncmp(buf,"showInfo=",9)==0) {
      yasol.yip.showInfo = buf[9] - '0';  
    }
    else if (strncmp(buf,"showWarning=",12)==0) {
      yasol.yip.showWarning = buf[12] - '0';      
    }
    else if (strncmp(buf,"showError=",10)==0) {
      yasol.yip.showError = buf[10] - '0';      
    }
    else if (strncmp(buf,"outputFile=",11)==0) {
      yasol.yip.writeOutputFile = buf[11] - '0';      
    }
    else if (strncmp(buf,"maintainPv=",11)==0) {
      yasol.yip.maintainPv = buf[11] - '0';            
    }
    else if (strncmp(buf,"learnDualCuts=",14)==0) {
      yasol.yip.learnDualCuts = buf[14] - '0';                  
    }
    else if (strncmp(buf,"useGMI=",7)==0) {
      yasol.yip.useGMI = buf[7] - '0';                        
    }
    else if (strncmp(buf,"useCover=",9)==0) {
      yasol.yip.useCover = buf[9] - '0';                              
    }
    else if (strncmp(buf,"useLPCuts=",10)==0) {
      yasol.yip.useLiftAndProjectCuts = buf[22] - '0';
    }
    else if (strncmp(buf,"useLazyLP=",10)==0) {
      yasol.yip.useLazyLP = buf[10] - '0';
    }
    else if (strncmp(buf,"useShadow=",10)==0) {
      yasol.yip.useShadow = buf[10] - '0';
    }
    else if (strncmp(buf,"useLimitedLP=",13)==0) {
      yasol.yip.useLimitedLP = buf[13] - '0';
    }
    else if (strncmp(buf,"useCglRootPreprocess=",21)==0) {
      yasol.yip.useCglRootPreprocess = buf[21] - '0';      
    }
    else if (strncmp(buf,"useCglRootCuts=",15)==0) {
      yasol.yip.useCglRootCuts = buf[15] - '0';      
    }
    else if (strncmp(buf,"useHighsHeuristic=",18)==0) {
      yasol.yip.useHighsHeuristic = buf[18] - '0';      
    }
    else if (strncmp(buf,"usePump=",8)==0) {
      yasol.yip.usePump = buf[8] - '0';      
    }
    else if (strncmp(buf,"useMonotones=",13)==0) {
      yasol.yip.useMonotones = buf[13] - '0';      
    }
    else if (strncmp(buf,"useAlphaCuts=",13)==0) {
      yasol.yip.useAlphaCuts = buf[13] - '0';      
    }
    else if (strncmp(buf,"useAlphabeta=",13)==0) {
      yasol.yip.useAlphabeta = buf[13] - '0';      
    }
    else if (strncmp(buf,"useMiniSearch=",14)==0) {
      yasol.yip.useMiniSearch = buf[14] - '0';      
    }
    else if (strncmp(buf,"useConflictGraph=",17)==0) {
      yasol.yip.useConflictGraph = buf[17] - '0';      
    }
    else if (strncmp(buf,"useStrongBranching=",19)==0) {
      yasol.yip.useStrongBranching = buf[19] - '0';      
    }
    else if (strncmp(buf,"reduceStrongBranching=",22)==0) {
      yasol.yip.reduceStrongBranching = buf[22] - '0';      
    }
    else if (strncmp(buf,"useFastFix=",11)==0) {
      yasol.yip.useFastFix = buf[11] - '0';      
    }
    else if (strncmp(buf,"isSimplyRestricted=",19)==0) {
      yasol.yip.isSimplyRestricted = buf[19] - '0';      
    }
    else cout << "Warning: Unknown option " << buf << endl;
    
    argc--;
  }
  if (argc > 2 && argv[2][0]!='-' && argv[2][0]!='R' ) {
    yasol.setInfoLevel(atoi(argv[2]));
  }

  if (argc > 3  && argv[3][0]!='-') {
    timeLimit = atoi(argv[3]);
  }
  
  std::string in = std::string(argv[1]);
  if (myrank == 0) cerr << "reading " <<in << endl;

  int DoneAfterSimplification=yasol.yInit(in);  
  if (DoneAfterSimplification){
    double ResInf;
    std::cout << "Solution process aborted; infeasible constraint system detected."<<std::endl;
    if(DoneAfterSimplification==1){
      std::cout << "The universal constraint system has no solution. The existential constraint system might be feasible."<< std::endl;
      ResInf=-((int64_t)1<<61)/1.5;
    }
    else{
      std::cout << "The existential constraint system has no solution. The universal constraint system might be feasible."<< std::endl;
      ResInf=((int64_t)1<<61);
    }
    if (yasol.objInverted) ResInf *= -1;
    std::cout << "RESULT: "<< ResInf  << std::endl;
    std::cout << "NODES: 0, #RELAXATIONS: 0";
    std::cout << ", LEARNED: 0" << std::endl;
    std::cout << "TIME: 0" << std::endl;
    delete yasolPt;
    return 0;
  }

  yasol.ySetProcessNo(myrank);

  if (myrank % 2 == 1) {
    ini_time = time(NULL);
    result = yasol.ySolve(ini_time);
    MPI_Finalize();
    delete yasolPt;
    return 0;
  }

  if(argc>2 &&(string)argv[2]=="Reduce"){
    yasol.Reduce(in);
    cerr << "created reduced QIP from QIPID " <<endl;
    delete yasolPt;
    return 0;
  }

  /*
  double pvsol[] = {0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,0.000000,0.000000,0.000000,1.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,0.000000,1.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,1.000000,0.000000,1.000000,0.000000,0.000000,0.000000,264.000000};
  int pvSolCnt=72;
  
  std::vector<data::QpNum> solu;
  //if(0)for (int i = 0; i < strlen(pvsol);i++)
  //    solu.push_back((double)((int)(pvsol[i]-'0')));
  if(0)for (int i = 0; i < pvSolCnt;i++)
    solu.push_back(pvsol[i]);
  //yasol.saveDepSolution(solu);
  */

  const int DEPfirst=0;
  if (DEPfirst) {
    cout << "Solving DEP..." << endl;
    algorithm::Algorithm::QlpSolution sol = yasol.solveDEP(true);

    cout << "DEP Solution: " << sol.getSolutionStatusString() << endl;
    cout << "DEP Objective Value: " << sol.getObjFunctionValue() << endl;
    cout << "DEP Vars: "
	 << data::QpNum::vecToString(sol.getSolutionVector()) << endl;
    std::vector<data::QpNum> solu = sol.getSolutionVector();
    yasol.saveDepSolution(solu);
  }

  ini_time = time(NULL);
  clock_t start, end;
  start = clock();
  bool useInteractiveMode=false;
  if (argc > 2 && argv[argc-1][0] == 'i') {
    useInteractiveMode=true;
  }
  if (timeLimit>=0) {
    yasol.setTimeout(time(NULL) + timeLimit);
  }
  if (argc > 3 && yasol.getTimeout() <= time(NULL)) {
    std::cout << "RESULT: " << "timeout" << " : " << in << std::endl;
    std::cout << "TIME: " << time(NULL)-ini_time << std::endl;
  } else if (yasol.getInfoLevel() > -10) {
    //cerr << "Objective Shift = " << yasol.offset_shift << endl;
    result = yasol.ySolve(ini_time);


    if (yasol.getInfoLevel() > 0) {
      if (yasol.getInfoLevel() >= 2) std::cout << "Und LOS #rows=" << yasol.getmRows() << " #cols=" << yasol.getnVars() << std::endl;
      if (yasol.getInfoLevel() >= 2) std::cout << "Objective Shift = " << yasol.offset_shift << endl;
      double time_taken = double(clock() - start) / double(CLOCKS_PER_SEC);
      if (yasol.getStatus() == YASOL_OPTIMAL){
        yasol.doShowSolution();
        yasol.qbp->WriteSolutionFile(+result+yasol.qbp->getFinalOffset(),time_taken, "OPTIMAL");
        std::cout << "Solution Status: OPTIMAL" << std::endl;
        if (yasol.objInverted) std::cout << "RESULT: " << -result-yasol.qbp->getFinalOffset() << " : " << in << std::endl;
        else std::cout << "RESULT: " << result+yasol.qbp->getFinalOffset() << " : " << in << std::endl;
      }
      else if (yasol.getStatus() == YASOL_FEASIBLE){
         yasol.doShowSolution();
        std::cout << "Solution Status: FEASIBLE" << std::endl;
        yasol.qbp->WriteSolutionFile(+result+yasol.qbp->getFinalOffset(),time_taken, "FEASIBLE");
      }
      else if (yasol.getStatus() == YASOL_INCUMBENT){
        yasol.doShowSolution();
        std::cout << "Solution Status: INCUMBENT" << std::endl;
        if (yasol.objInverted) std::cout << "RESULT: " << -result-yasol.qbp->getFinalOffset() << " : " << in << std::endl;
        else std::cout << "RESULT: " << result+yasol.qbp->getFinalOffset() << " : " << in << std::endl;
        std::cout << "GAP: " << yasol.getGap() << std::endl;

        yasol.qbp->WriteSolutionFile(+result+yasol.qbp->getFinalOffset(),time_taken, "INCUMBENT");
      }
      else if (yasol.getStatus() == YASOL_UNSAT){
        std::cout << "Solution Status: UNSAT" << std::endl;
      }
      else if (yasol.getStatus() == YASOL_UNKNOWN){
        std::cout << "terminated without incumbent — feasibility not established." << std::endl;
        std::cout << "Solution Status: UNKNOWN" << std::endl;
      }
      else if (yasol.getStatus() == YASOL_ERROR){
        std::cout << "Something went wrong..." << std::endl;
        std::cout << "Solution Status: ERROR" << std::endl;
      }
      else {
        std::cout << "Solution Status: I do not even know this status" << std::endl;
        std::cout << "Something went terribly wrong..." << std::endl;
      }
      std::cout << "FINAL NODES: " << yasol.getNumberOfDecisions()+yasol.getNumberOfPropagationSteps() << ", #RELAXATIONS:"<< yasol.getNumberOfQlpCalls();
      std::cout << ", LEARNED: " << yasol.getNumberOfLearntConstraints() << std::endl;
      std::cout << "TIME: " << time_taken << std::endl;
    } else {
      cerr << "result=" << result << " yasol.getNinf() = " << yasol.getNinf() << endl;
      if (yasol.getNinf() == result || yasol.getNinf() == -result) cout << "0" << endl;
      else {
	cout << "1" << endl << endl;
	cout << result << endl << endl;
	cout << "0" << endl;
      }
    }
    if (useInteractiveMode) {
      char input[1000];
      bool finished=false;
      while (!finished) {
	fgets(input,999,stdin);
	if (!strncmp(input,"quit",4)) {
	  cerr << "info quit" << endl;
	  finished = true;
	} else if (!strncmp(input,"get nodeinfo ",13)) {
	  int nodeID = atoi(&input[13]);
	  cerr << "info get nodeinfo " << nodeID << endl;
	  yasol.write_nodeinfo(nodeID);
	} else if (!strncmp(input,"get successors ",15)) {
	  int nodeID = atoi(&input[15]);
	  cerr << "info get successors of node " << nodeID << endl;
	  yasol.write_successors(nodeID);
	} else if (!strncmp(input,"get varname ",12)) {
	  int v = atoi(&input[12]);
	  yasol.writeRealName(v);
	} else if (!strncmp(input,"update node ",12)) {
	  int v = atoi(&input[12]);
	  yasol.updateNode(v);
	}
      }
    }
  }
    
  if (yasol.getInfoLevel() > 5 || yasol.getInfoLevel() < -5) {
    ini_time = time(NULL);
    cerr << "Solving DEP..." << endl;
    if (yasol.getInfoLevel() > -1000) {
      if (yasol.getInfoLevel() > -50) {
	algorithm::Algorithm::QlpSolution sol = yasol.solveDEP();
	cout << "DEP Solution: " << sol.getSolutionStatusString() << endl;
	cout << "DEP Objective Value: " << sol.getObjFunctionValue() << endl;
	cout << "DEP Vars: " << data::QpNum::vecToString(sol.getSolutionVector()) << endl;
	std::cout << "RESULT: " << sol.getObjFunctionValue().asDouble() << " : " << in << std::endl;
	std::cout << "TIME: " << time(NULL) - ini_time << std::endl;
        //yasol.qbp->WriteSolutionFile(0.0,0.0, "OPTIMAL");
      } else {
	algorithm::Algorithm::QlpSolution sol = yasol.solveDEP();
	if (sol.getSolutionStatusString().compare("Feasible")) {
	  cout << "1" << endl << endl;
	  cout << sol.getObjFunctionValue() << endl << endl;
	  cout << "0" << endl;
	} else {
	  cout << "0" << endl;
	}
      }
    } else {
      yasol.dummyDEP();
      std::cout << "RESULT: " << " no result" << " : " << in << std::endl;
      std::cout << "TIME: " << time(NULL) - ini_time << std::endl;
    }
  }
  MPI_Finalize();
  delete yasolPt;
  return 0;
}

