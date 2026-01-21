/*
*
* Yasol: yInterface.cc -- Copyright (c) 2012-2017 Ulf Lorenz
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

#include "yInterface.h"
#include "CutAdder.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void yInterface::yReadIniFile() {
	  string line;
	  ifstream myfile ("./Yasol.ini");
	  std::vector<std::string> tokens;
	  size_t p = string::npos;
	  if (myfile.is_open())
	  {
	    while ( getline (myfile,line) )
	    {
	      //cerr << line << '\n';
	      if (line.compare("END") == 0) break;
	      tokens.clear();
	      p = line.find_first_of("=", 0);
		  tokens.push_back(line.substr(0, p));
		  tokens.push_back(line.substr(p+1, 1000));
		  //cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;

		  if (tokens[0].compare("learnDualCuts") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setLearnDualCuts(true);
			  else qbp->setLearnDualCuts(false);
          } else if (tokens[0].compare("useShadow") == 0) {
              if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
              if (tokens[1].compare("1") == 0) qbp->setUseShadow(true);
              else qbp->setUseShadow(false);
          } else if (tokens[0].compare("useLazyLP") == 0) {
              if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
              if (tokens[1].compare("1") == 0) qbp->setUseLazyLP(true);
              else qbp->setUseLazyLP(false);
		  } else if (tokens[0].compare("usePump") == 0) {
                          if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
                          if (tokens[1].compare("1") == 0) qbp->setUsePump(true);
                          else qbp->setUsePump(false);
		  } else if (tokens[0].compare("useMiniSearch") == 0) {
                          if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
                          if (tokens[1].compare("1") == 0) qbp->setUseMiniSearch(1);
			  else if (tokens[1].compare("0") == 0) qbp->setUseMiniSearch(0);
			  else if (tokens[1].compare("2") == 0) qbp->setUseMiniSearch(2);
			  else if (tokens[1].compare("3") == 0) qbp->setUseMiniSearch(3);
			  else if (tokens[1].compare("4") == 0) qbp->setUseMiniSearch(4);
			  else if (tokens[1].compare("5") == 0) qbp->setUseMiniSearch(5);
			  else if (tokens[1].compare("6") == 0) qbp->setUseMiniSearch(6);
			  else if (tokens[1].compare("7") == 0) qbp->setUseMiniSearch(7);
			  else if (tokens[1].compare("8") == 0) qbp->setUseMiniSearch(8);
			  else if (tokens[1].compare("9") == 0) qbp->setUseMiniSearch(9);
			  else if (tokens[1].compare("10") == 0) qbp->setUseMiniSearch(10);
			  else if (tokens[1].compare("11") == 0) qbp->setUseMiniSearch(11);
			  else if (tokens[1].compare("12") == 0) qbp->setUseMiniSearch(12);
			  else if (tokens[1].compare("13") == 0) qbp->setUseMiniSearch(13);
			  else if (tokens[1].compare("14") == 0) qbp->setUseMiniSearch(14);
			  else if (tokens[1].compare("15") == 0) qbp->setUseMiniSearch(15);
                          else qbp->setUseMiniSearch(0);
			  //Meaning: - 0 : no minisearch
			  //         - 1 : minisearch yes, uviRELAX = false
			  //         - 2 : minisearch no, uviRELAX = yes
			  //         - 3 : minisearch yes, uviRELAX = yes
			  //         - >3: bit 1 and 2 analogously, but early finish when minisearch too long
		  } else if (tokens[0].compare("useConflictGraph") == 0) {
                          if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
                          if (tokens[1].compare("1") == 0) qbp->setUseConflictGraph(1);
			  else if (tokens[1].compare("2") == 0) qbp->setUseConflictGraph(2);
                          else qbp->setUseConflictGraph(0);
		  } else if (tokens[0].compare("maintainPv") == 0) {
                          if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
                          if (tokens[1].compare("1") == 0) qbp->setMaintainPv(true);
                          else qbp->setMaintainPv(false);
		  } else if (tokens[0].compare("useGMI") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("0") == 0) qbp->setUseGMI(0);
			  else if (tokens[1].compare("1") == 0) qbp->setUseGMI(1);
			  else if (tokens[1].compare("2") == 0) qbp->setUseGMI(2);
			  else if (tokens[1].compare("3") == 0) qbp->setUseGMI(3);			  
		  } else if (tokens[0].compare("useHighsHeuristic") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (qbp->getUseHighsH() >= 0) {
			    if (tokens[1].compare("1") == 0) qbp->setUseHighsH(true);
			    else qbp->setUseHighsH(false);
			  } else qbp->setUseHighsH(false);
		  } else if (tokens[0].compare("useAlphaCuts") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseAlphaCuts(true);
			  else qbp->setUseAlphaCuts(false);
		  } else if (tokens[0].compare("showInfo") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("0") == 0) qbp->setShowInfo(0);
			  else if (tokens[1].compare("1") == 0) qbp->setShowInfo(1);
			  else if (tokens[1].compare("2") == 0) qbp->setShowInfo(2);
			  else if (tokens[1].compare("3") == 0) qbp->setShowInfo(3);
			  else if (tokens[1].compare("4") == 0) qbp->setShowInfo(4);
			  else if (tokens[1].compare("5") == 0) qbp->setShowInfo(5);
			  else if (tokens[1].compare("6") == 0) qbp->setShowInfo(6);
			  else if (tokens[1].compare("7") == 0) qbp->setShowInfo(7);
			  else if (tokens[1].compare("8") == 0) qbp->setShowInfo(8);
			  else if (tokens[1].compare("9") == 0) qbp->setShowInfo(9);
			  
		  } else if (tokens[0].compare("showWarning") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setShowWarning(true);
			  else qbp->setShowWarning(false);
		  } else if (tokens[0].compare("showError") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setShowError(true);
			  else qbp->setShowError(false);
		  } else if (tokens[0].compare("useCover") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseCover(true);
			  else qbp->setUseCover(false);
		  } else if (tokens[0].compare("useMcts") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseMcts(true);
			  else qbp->setUseMcts(false);
		  } else if (tokens[0].compare("useCglRootCuts") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseCglRootCuts(true);
			  else qbp->setUseCglRootCuts(false);
		  } else if (tokens[0].compare("useCglRootPreprocess") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseCglRootPreprocess(true);
			  else qbp->setUseCglRootPreprocess(false);
		  } else if (tokens[0].compare("reduceStrongBranching") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setReduceStrongBranching(1);
			  else if (tokens[1].compare("0") == 0) qbp->setReduceStrongBranching(0);
			  else qbp->setReduceStrongBranching(2);
		  } else if (tokens[0].compare("useLiftAndProjectCuts") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseLaP(true);
			  else qbp->setUseLaP(false);
		  } else if (tokens[0].compare("useFarJumping") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseBendersBackJump(true);
			  else qbp->setUseBendersBackJump(false);
		  } else if (tokens[0].compare("useFastFix") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseFastFix(true);
			  else qbp->setUseFastFix(false);
		  } else if (tokens[0].compare("useImplications") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseImplications(true);
			  else qbp->setUseImplications(false);
		  } else if (tokens[0].compare("useLimitedLP") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseLimitedLP(true);
			  else qbp->setUseLimitedLP(false);
		  } else if (tokens[0].compare("useStrongBranching") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseStrongBranching(true);
			  else qbp->setUseStrongBranching(false);
		  } else if (tokens[0].compare("useEarlyBackjump") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseEarlyBackjump(true);
			  else qbp->setUseEarlyBackjump(false);
		  } else if (tokens[0].compare("useBestLevelExtraction") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseBestLevelExtraction(true);
			  else qbp->setUseBestLevelExtraction(false);
		  } else if (tokens[0].compare("useUniversalBackjump") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseUniversalBackjump(true);
			  else qbp->setUseUniversalBackjump(false);
		  } else if (tokens[0].compare("maxBaCLevel") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  int x;
			  if ( ! (istringstream(tokens[1]) >> x) ) x = 0;
			  qbp->setMaxBaCLevel(x);
		  } else if (tokens[0].compare("useScout") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseScout(true);
			  else qbp->setUseScout(false);
		  } else if (tokens[0].compare("useAlphabeta") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseAlphabeta(true);
			  else qbp->setUseAlphabeta(false);
		  } else if (tokens[0].compare("useMonotones") == 0) {
			  if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
			  if (tokens[1].compare("1") == 0) qbp->setUseMonotones(true);
			  else qbp->setUseMonotones(false);
		  } else if (tokens[0].compare("isSimplyRestricted") == 0) {
                          if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
                          if (tokens[1].compare("1") == 0) qbp->setIsSimplyRestricted(true);
                          else qbp->setIsSimplyRestricted(false);
                  } else if (tokens[0].compare("writeOutputFile") == 0) {
                          if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
                          if (tokens[1].compare("1") == 0) qbp->setWriteOutputFile(true);
                          else qbp->setWriteOutputFile(false);
                  } else if (tokens[0].compare("gzipExists") == 0) {
                          if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
                          if (tokens[1].compare("1") == 0) gzipExists = true;
                          else                             gzipExists = false;
                  } else if (tokens[0].compare("qdimacs2qipConverter") == 0) {
                          if (info_level >= 10) cerr << "Str1=" << tokens[0] << " und Str2=" << tokens[1] << endl;
                          setQdimacs2QIPPath(tokens[1]);
                  }                  
	    }
	    myfile.close();
	  } else if(supressAllOutput) cerr << "Info: Unable to find Yasol.ini. I use standard parameter setting instead."<<endl;
	  if (yip.showInfo >= 0)
	    qbp->setShowInfo(yip.showInfo);
	  if (yip.showWarning >= 0)
	    qbp->setShowWarning(yip.showWarning);	    
	  if (yip.showError >= 0)
	    qbp->setShowError(yip.showError);	    
	  if (yip.writeOutputFile >= 0)
	    qbp->setWriteOutputFile(yip.writeOutputFile);
	  if (yip.maintainPv >= 0)
	    qbp->setMaintainPv(yip.maintainPv);
	  if (yip.learnDualCuts >= 0)
	    qbp->setLearnDualCuts(yip.learnDualCuts);
	  if (yip.useGMI >= 0)
	    qbp->setUseGMI(yip.useGMI);
	  if (yip.useCover >= 0)
	    qbp->setUseCover(yip.useCover);
	  if (yip.useLiftAndProjectCuts >= 0)
	    qbp->setUseLaP(yip.useLiftAndProjectCuts);
	  if (yip.useLazyLP >= 0)
	    qbp->setUseLazyLP(yip.useLazyLP);
	  if (yip.useShadow >= 0)
	    qbp->setUseShadow(yip.useShadow);
	  if (yip.useLimitedLP >= 0)
	    qbp->setUseLimitedLP(yip.useLimitedLP);
	  if (yip.useCglRootPreprocess >= 0)
	    qbp->setUseCglRootPreprocess(yip.useCglRootPreprocess);
	  if (yip.useCglRootCuts >= 0)
	    qbp->setUseCglRootCuts(yip.useCglRootCuts);
	  if (yip.useHighsHeuristic >= 0)
	    qbp->setUseHighsH(yip.useHighsHeuristic);
	  if (yip.usePump >= 0)
	    qbp->setUsePump(yip.usePump);
	  if (yip.useMonotones >= 0)
	    qbp->setUseMonotones(yip.useMonotones);
	  if (yip.useAlphaCuts >= 0)
	    qbp->setUseAlphaCuts(yip.useAlphaCuts);
	  if (yip.useAlphabeta >= 0)
	    qbp->setUseAlphabeta(yip.useAlphabeta);
	  if (yip.useMiniSearch >= 0)
	    qbp->setUseMiniSearch(yip.useMiniSearch);
	  if (yip.useConflictGraph >= 0)
	    qbp->setUseConflictGraph(yip.useConflictGraph);
	  if (yip.useStrongBranching >= 0)
	    qbp->setUseStrongBranching(yip.useStrongBranching);
	  if (yip.reduceStrongBranching >= 0)
	    qbp->setReduceStrongBranching(yip.reduceStrongBranching);
	  if (yip.useFastFix >= 0)
	    qbp->setUseFastFix(yip.useFastFix);
	  if (yip.isSimplyRestricted >= 0)
	    qbp->setIsSimplyRestricted(yip.isSimplyRestricted);
}

void yInterface::yReadInput(int htable_size) {
	yReadIniFile();
	unsigned int i,j;
	int k,l,block,cblock;

	//get qlp-handles for constraints and variables
	int numVars = qlp.getVariableCount();
	int numConstraints = qlp.getConstraintCount();
	std::vector<const data::QpVar *> varVec = qlp.getVariableVectorConst();
	std::vector<const data::QpVar *> C_varVec = qlpRelax.getVariableVectorConst();
	std::vector<const data::QpRhs *> rhsVec = qlp.getRhsVecConst();
	std::vector<const data::Constraint *> conVec = qlp.getConstraintVecConst();
	//const data::QpRhs& rhs = conVec[0]->getRhs();
	qbp->initializeCliqueManager(varVec.size());
	qbp->setObjInverted(objInverted);

	if(qbp->getShowInfo()) cerr << "QLP ready." << endl;

    qbp->setInfoLevel(info_level);
    qbp->initializeHashTables(varVec.size(), htable_size,10000000);

	if(qbp->getShowInfo()) cerr << "Hashtable ready" << endl;

	for (i=0; i < varVec.size(); i++) {
		if (i == 0) {
		  block = 1;
		  cblock = 2;
		} else {
		  if (varVec[i]->getQuantifier() != varVec[i-1]->getQuantifier() ) block++;
		  if (C_varVec[i]->getQuantifier() != C_varVec[i-1]->getQuantifier() ) cblock+=2;
		}
		if (varVec[i]->getQuantifier() == data::QpVar::all) qbp->createVar(UNIV,block, cblock, varVec.size(),varVec[i]->getNumberSystem()==data::QpVar::real, varVec[i]->getLowerBound().asDouble(),varVec[i]->getUpperBound().asDouble());
		else if (varVec[i]->getQuantifier() == data::QpVar::exists) qbp->createVar(EXIST,block, cblock, varVec.size(),varVec[i]->getNumberSystem()==data::QpVar::real, varVec[i]->getLowerBound().asDouble(),varVec[i]->getUpperBound().asDouble());
		if (qbp->getType(i) == CONTINUOUS) qbp->setCblock(qbp->nVars()-1,qbp->getCblock(qbp->nVars()-1) + 1);
	}
	bool lastVarIsUniv=false;
	if(varVec[varVec.size()-1]->getQuantifier() == data::QpVar::all){
		//block++;
		//cblock+=2;
		//qbp->createVar(EXIST,block, cblock, varVec.size()+1,false, 0,1);
	  cerr << "Error!: Final Variable in ORDER has to be an existential variable.\nAdd a (binary) existential dummy variable to the end of the instance's ORDER!" << endl;
	  assert(0);
	  lastVarIsUniv=true;
	}
	qbp->setMaxBlock(block);
	qbp->setMaxLPBlock(cblock / 2);
	qbp->setMaxLPStage(cblock / 4);
	if(qbp->getShowInfo()){
	    cerr << "Variables ready" << endl;
	    cerr << "Blocks:" << qbp->getMaxBlock() << " ready" << endl;
	    cerr << "C-Blocks:" << qbp->getMaxLPBlock() << " Stages:" << qbp->getMaxLPStage() << " ready" << endl;
	}
	qbp->clearReverseImplicationQueue();

	//initialize objective
    CoeVar c;
    coef_t r;
	ca_vec<CoeVar> newConstraint;
	//if( qlp.getObjective() == data::QpObjFunc::max ) qlp.reverseObjFunc();
    data::QpNum objOffSet=qlp.getObjFuncOffset();
    std::vector<data::QpNum> objCoeffVec = qlp.getObjectiveFunctionValues(); //Dense Vector

    if (true) {
		newConstraint.clear();
		coef_t lb=(coef_t)0, ub=(coef_t)0;
		int cnt_negs=0;
		int cnt_elems=0;
		for( j = 0; j < objCoeffVec.size(); j++) {
			if (!objCoeffVec[j].isZero()) {
				cnt_elems++;
				c = mkCoeVar(j,
					            objCoeffVec[j].asDouble() >= 0 ? (coef_t)objCoeffVec[j].asDouble() : -(coef_t)objCoeffVec[j].asDouble(),
					            objCoeffVec[j].asDouble() >= 0 ? true : false); // negieren, weil wir maximieren wollen
				newConstraint.push(c);
				//std::cout << (sign(c) ? "-" : "") << c.coef << "x" << var(c) << endl;
				//std::cout << "O " <<  objCoeffVec[j].asDouble() << "x" << var(c) << endl;
				if (sign(c)) {
					cnt_negs++;
					lb -= c.coef*qbp->getUpperBound(var(c));
					ub -= c.coef*qbp->getLowerBound(var(c));
				} else {
					ub += c.coef*qbp->getUpperBound(var(c));
					lb += c.coef*qbp->getLowerBound(var(c));
				}
			}
		}
		qbp->setFeasibilityOnly(false);
		if (cnt_elems == 0) {
			int j = 0;
			for (int i = 0; i < qbp->nVars();i++)
				if (qbp->getMaxBlock() == qbp->getBlock(i)) {
					j = i;
					break;
				}
			cnt_elems++;
			c = mkCoeVar(j, 1.0, true);
			newConstraint.push(c);
			//std::cout << (sign(c) ? "-" : "") << c.coef << "x" << var(c) << " ";
			if (sign(c)) {
				cnt_negs++;
				lb -= c.coef;
			} else {
				ub += c.coef;
			}
			qbp->setFeasibilityOnly(true);
		}
		if (cnt_elems > 0) {
			qbp->setHasObjective(true);
			qbp->definePositiveInfinity(ub);
			qbp->defineNegativeInfinity((coef_t)(-((int64_t)1<<61)));
			qbp->defineAllInfeasible();
			qbp->defineDontKnowValue((coef_t)(-((int64_t)1<<60)));
			if (lb <= qbp->getDontKnowValue()) {
				lb = (coef_t)(-((int64_t)1<<59));
				cerr << "Note: lower bound of objective artificially limited by -2^59" << endl;
			}
			//std::cerr << " lb=" << lb << " ub=" << ub << " size of new constraint: " << newConstraint.size() << endl;;
			qbp->addObjective(newConstraint,lb); // non-sat!
			//cerr << "#Constraints=" << qbp->constraints.size() << endl;
			//cerr << "-inf=" << qbp->n_infinity << " +inf=" << qbp->p_infinity << " dont_know=" << qbp->dont_know << " lb=" << lb << " ub=" << ub << endl;
		} else {
			qbp->setHasObjective(false);
			qbp->definePositiveInfinity((coef_t)1);
			qbp->defineNegativeInfinity((coef_t)-2);
			qbp->defineAllInfeasible();
			qbp->defineDontKnowValue((coef_t)0);
		}
    }

    if(qbp->getShowInfo()) cerr << "Objective ready" << endl;

	//initialize constraints
    bool universalConstraintsExist=false;
    for (i = 0; i < numConstraints;i++) {
        bool isUniversal;
        bool isBndCon = false;
        int cntReals = 0;
        int rvar=-1;
        //cerr << "isUni:" << (int)rhsVec[i]->getResponsibility() << endl;
        if (rhsVec[i]->getResponsibility()==data::Constraint::UNIVERSAL) {
            isUniversal = true;
            if(qbp->getShowInfo())cerr << "Info: Found universal player responsibility." << endl;
	    universalConstraintsExist=true;
        } else if (rhsVec[i]->getResponsibility()==data::Constraint::UNDEF) {
            isUniversal = false;
            //cerr << "Error: Resonsibilty undefined." << endl;
        } else {
            isUniversal = false;
        }
        int countR=0;
        while(countR<2){

	    	std::vector<data::IndexedElement> lhs = conVec[i]->getElements();

	    	if(rhsVec[i]->getRatioSign()==data::QpRhs::smallerThanOrEqual) {
	        	newConstraint.clear();
	            cntReals = 0;
	            isBndCon = false;
	            int rvar=-1;
	        	for( j = 0; j < lhs.size(); j++) {
	        		c = mkCoeVar(lhs[j].index,
	        				            lhs[j].value.asDouble() >= 0 ? (coef_t)lhs[j].value.asDouble() : -(coef_t)lhs[j].value.asDouble(),
	        				            lhs[j].value.asDouble() >= 0 ? true : false);
	        		newConstraint.push(c);
	        		//std::cout << (sign(c) ? "-" : "") << c.coef << "x" << var(c) << " ";
	                if (varVec[ lhs[j].index ]->getNumberSystem()==data::QpVar::real) {
	                    cntReals++;
	                    rvar = lhs[j].index;
	                }
	        	}
	            if (cntReals == 1) {
	                isBndCon = true;
	                //cerr << "lbc1 ";
	            }
	            r = -(rhsVec[i]->getValue()).asDouble();
	        	qbp->addOrgConstraint(newConstraint,r,0,isUniversal,rvar,isBndCon);
	    	} else if(rhsVec[i]->getRatioSign()==data::QpRhs::equal) {
		  newConstraint.clear();
		  cntReals = 0;
		  isBndCon = false;
		  int rvar=-1;
		  int rvarCoef=0.0;
		  for( j = 0; j < lhs.size(); j++) {
		    c = mkCoeVar(lhs[j].index,
				 lhs[j].value.asDouble() >= 0 ? (coef_t)lhs[j].value.asDouble() : -(coef_t)lhs[j].value.asDouble(),
				 lhs[j].value.asDouble() >= 0 ? false : true);
		    newConstraint.push(c);
		    if (varVec[ lhs[j].index ]->getNumberSystem()==data::QpVar::real) {
		      cntReals++;
		      rvar = lhs[j].index;
		    }
		  }
		  if (cntReals == 1) {
		    isBndCon = true;
		    //cerr << "lbc2 ";
		  }
		  r = (rhsVec[i]->getValue()).asDouble();
		  qbp->addOrgConstraint(newConstraint,r,0,isUniversal,rvar,isBndCon);
		  newConstraint.clear();
		  cntReals = 0;
		  isBndCon = false;
		  /*int*/ rvar=-1;
		  for( j = 0; j < lhs.size(); j++) {
		    c = mkCoeVar(lhs[j].index,
				 lhs[j].value.asDouble() >= 0 ? (coef_t)lhs[j].value.asDouble() : -(coef_t)lhs[j].value.asDouble(),
				 lhs[j].value.asDouble() >= 0 ? true : false);
		    newConstraint.push(c);
		    if (varVec[ lhs[j].index ]->getNumberSystem()==data::QpVar::real) {
		      cntReals++;
		      rvar = lhs[j].index;
		      rvarCoef = lhs[j].value.asDouble();
		    }
		  }
		  if (cntReals == 1) {
		    isBndCon = true;
		    //cerr << "lbc3 ";
		  }
		  r = -(rhsVec[i]->getValue()).asDouble();
		  qbp->addOrgConstraint(newConstraint,r,0,isUniversal,rvar,isBndCon);
		  if(0&&qbp->getShowInfo()){
  		    if (cntReals>0) cerr << "info: c" << i << " reals:" << cntReals << "; ";
		    if (cntReals==2) {
		      for( int j = 0; j < lhs.size(); j++) {
		        if (varVec[ lhs[j].index ]->getNumberSystem()==data::QpVar::real) {
			  cerr << " x" << lhs[j].index;
		        }
		      }
		    }
  		    if (cntReals>0) cerr << endl;
		  }
		  if (cntReals > 0/*== 1*/) {
		    if(0&&qbp->getShowInfo()) cerr << "info: can eliminate a continuous variable" << endl;
		    int foundCand=0;
		    rvar = -1;
		    for( int j = 0; j < lhs.size(); j++) {
		      if (varVec[ lhs[j].index ]->getNumberSystem()==data::QpVar::real /*|| 
											 (qbp->getIsIntegerBit(lhs[j].index))*/ ) {
			if (varVec[ lhs[j].index ]->getNumberSystem()==data::QpVar::real) foundCand = 1;
			else if (qbp->getIsIntegerBit(lhs[j].index)) foundCand = 2;
			else foundCand = 3;
			if(0&&qbp->getShowInfo()) cerr << "info: can eliminate x" << lhs[j].index << " with mode=" << foundCand << endl;
			qbp->updateReplacement(lhs[j].value.asDouble(), lhs[j].index, lhs, rhsVec[i]->getValue().asDouble(),foundCand,cntReals);
			if (rvar == -1 || lhs[j].index < rvar) {
			  foundCand=1;
			  rvar = lhs[j].index;
			  rvarCoef = lhs[j].value.asDouble();
			}
		      }
		    }
		    if (rvar<0/*&& lhs.size() < 10*/) {
		      for( int j = 0; j < lhs.size(); j++) {
			if (fabs(fabs(lhs[j].value.asDouble())-1.0) < 1e-9 && qbp->getIsIntegerBit(lhs[j].index)) {
			  if (rvar == -1 || lhs[j].index < rvar) {
			    if (qbp->getBlock(lhs[j].index) == qbp->getMaxBlock()) {
			      foundCand=2;
			      rvar = lhs[j].index;
			      rvarCoef = lhs[j].value.asDouble();
			    }
			  }
			}
		      }
		    }
		    if (rvar<0 /*&& lhs.size() < 4*/) {
		      for( int j = 0; j < lhs.size(); j++) {
			if (fabs(fabs(lhs[j].value.asDouble())-1.0) < 1e-9) {
			  if (rvar == -1 || lhs[j].index < rvar) {
			    if (qbp->getBlock(lhs[j].index) == qbp->getMaxBlock()) {
			      foundCand=3;
			      rvar = lhs[j].index;
			      rvarCoef = lhs[j].value.asDouble();
			    }
			  }
			}
		      }
		    }
		    //if (rvar >= 0) qbp->updateReplacement(rvarCoef, rvar, lhs, rhsVec[i]->getValue().asDouble(),foundCand,cntReals);
		  }
	    	} else {
		  newConstraint.clear();
		  cntReals = 0;
		  isBndCon = false;
		  int rvar=-1;
		  for( j = 0; j < lhs.size(); j++) {
		    c = mkCoeVar(lhs[j].index,
				 lhs[j].value.asDouble() >= 0 ? (coef_t)lhs[j].value.asDouble() : -(coef_t)lhs[j].value.asDouble(),
				 lhs[j].value.asDouble() >= 0 ? false : true);
		    newConstraint.push(c);
		    if (varVec[ lhs[j].index ]->getNumberSystem()==data::QpVar::real) {
		      cntReals++;
		      rvar = lhs[j].index;
		    }
		  }
		  if (cntReals == 1) {
		    isBndCon = true;
		    //cerr << "lbc4 ";
		  }
		  r = (rhsVec[i]->getValue()).asDouble();
		  qbp->addOrgConstraint(newConstraint,r,0,isUniversal,rvar, isBndCon);
	    	}
	    	if(isUniversal){
	    		countR=2;
	    		//isUniversal=false;
	    	}
	    	else countR=2;
	    }
    }

    if (universalConstraintsExist && lastVarIsUniv) {
      if (lastVarIsUniv) {
	if(qbp->getShowError()) cerr << "Error: Final Variable in ORDER has to be an existential variable.\nAdd a (binary) existential dummy variable to the end of the instance's ORDER!" << endl;
	assert(0);
      }
    } else if (lastVarIsUniv) {
      if(qbp->getShowError()) cerr << "Error II: Final Variable in ORDER has to be an existential variable.\nAdd a (binary) existential dummy variable to the end of the instance's ORDER!" << endl;
      assert(0);
    }

    
    if(qbp->getShowInfo()) cerr << "Constraints ready" << endl;

    qbp->initializeMrows();
    //qbp->initializeConstraintWatcher();
    //yReadIniFile();
    //qbp->transferBoundsVars2Constraints();

    if(qbp->getShowInfo()) cerr << "info: prepared QBP object." << endl;
}


double yInterface::computeCutRatio(vector< pair<unsigned int, double> >& cut) {
	if (cut.size() < 1) return 10000000.0;
	double mx=fabs(cut[0].second);
	double mn=mx;
    for (int k=1; k < cut.size();k++) {
      if (fabs(cut[k].second) > mx) mx = fabs(cut[k].second);
      if (fabs(cut[k].second) < mn && fabs(cut[k].second) > 0) mn = fabs(cut[k].second);
    }
    return mx / mn;
}

#define isOne_(a) ((a) > 1.0-1e-12 && (a) < 1.0+1e-12)
#define isZero_(a) ((a) > 0.0-1e-12 && (a) < 0.0+1e-12)

// Cuts must be of type >=
#include "ExternSolvers/QpExtSolCLP.hpp"

int yInterface::GenerateCutAndBranchCuts( extSol::QpExternSolver& externSolver, vector< vector< data::IndexedElement > > &listOfCutsLhs,
		       vector< data::QpNum > &listOfCutsRhs, vector<int> &listOfCutsVars,
					  int treedepth, int currentBlock, bool &global_valid, std::vector<unsigned int> &candidates, int cuttype, int*types, int8_t* assigns, unsigned int initime, int* solu, int *fixs, int *blcks, int *eass, int orgN, double intLB) {
	static bool root_once = false;
	global_valid = false;
	root_once = true;
	const int CGL_KNAPSACK = 16;
	const int CGL_SIMROUND = 32;
	const int CGL_GMI      = 64;
	const int CGL_LIFTAPRO = 128;
	const int CGL_FLOWCOVR = 256;
	const int CGL_REDSPLIT = 512;
	const int CGL_MIR2     = 1024;
	const int CGL_TWOMIR   = 2048;
	const int CGL_PROB     = 4096;
	const int MIRsmart     = 8192;

	vector< pair< vector< pair<unsigned int, double> >, double > > *cuts_pt;
	vector< pair< vector< pair<unsigned int, double> >, double > > single;
	cuts_pt = &single;
	if (cuttype & CGL_LIB) {
	  //64 GMI
	  if (qbp->getInfoLevel()>-8) cerr << " first entry." << endl;
	  int useCglPreProcess=0;
	  if (qbp->getUseCglRootPreprocess()) 
	    useCglPreProcess = 2;
#ifndef NO_CGL
	  if (qbp->getInfoLevel()>-8) cerr << "Row Count is " << cutGen.solver->CutGenSolver->getRowCount() <<endl;
	  cuts_pt = cutGen.solver->CutGenSolver->CreateCuts(externSolver, types, assigns, treedepth, initime, solu, fixs, blcks, orgN,
					       (cuttype & 1) ? //1 means: break the loop
					       //CGL_KNAPSACK|CGL_SIMROUND|CGL_GMI|CGL_FLOWCOVR|CGL_REDSPLIT|CGL_MIR2|1:
					       CGL_KNAPSACK|CGL_SIMROUND|CGL_GMI|CGL_FLOWCOVR|CGL_REDSPLIT|CGL_MIR2|CGL_TWOMIR|CGL_PROB|useCglPreProcess:
					       //CGL_KNAPSACK|CGL_SIMROUND|CGL_GMI|CGL_FLOWCOVR|CGL_REDSPLIT|CGL_MIR2|1,true);
					       CGL_KNAPSACK|CGL_SIMROUND|CGL_GMI|CGL_FLOWCOVR|CGL_REDSPLIT|CGL_MIR2|CGL_TWOMIR|CGL_PROB|1|useCglPreProcess,true, -intLB);
#else
	  if(qbp->getShowInfo()) std::cerr << "Info: Cgl is not implemented." << std::endl;
#endif
	  if (qbp->getInfoLevel()>-8) cerr <<"Size of Cuts after: "<<cuts_pt->size() <<endl;
	} else if (cuttype & GMI) {
	  cuts_pt = CutAdder::getGMICuts( externSolver, candidates, types, assigns, treedepth, initime, listOfCutsVars, solu , fixs, blcks, eass, orgN);
        vector< pair< vector< pair<unsigned int, double> >, double > > &cuts = *cuts_pt;
	for( unsigned int i = 0; i < cuts.size(); ++i ){
	  vector< pair<unsigned int, double> >& cut = cuts.at( i ).first;
	  vector< data::IndexedElement > addcut;

	    if (0&&cut.size()<7 && fabs(cuts[i].second+0.25) < 0.01) {
	      cerr << "before gen GMI cut as comes in: ";
	      for (int v=0;v<cut.size();v++) {
		cerr << cut[v].second << "y" << cut[v].first << " + ";
	      }
	      cerr << " >= " << cuts[i].second <<endl;
	    }
	}
	} else if (cuttype & LaP) {
		cuts_pt = CutAdder::getLPCuts( externSolver, candidates, types, assigns, treedepth, initime, listOfCutsVars );
	} else if (cuttype & MirSmart) {
	  cuts_pt = CutAdder::getMIRsmartEQ( externSolver, types, assigns, treedepth, initime, solu, fixs, blcks, eass, orgN );
        } else if (cuttype & Cover) {
	  cuts_pt = CutAdder::getCoverCuts( externSolver, types, assigns, treedepth, initime, solu, fixs, blcks, eass, orgN );
	} else if (cuttype & UserCut) {
	  //cuts_pt = CutAdder::getAggCmirCuts( externSolver, types, assigns, treedepth, initime, solu, fixs, blcks,eass, orgN );
	}

	vector< pair< vector< pair<unsigned int, double> >, double > > &cuts = *cuts_pt;
	listOfCutsLhs.reserve( listOfCutsLhs.size() + cuts.size() );
	listOfCutsRhs.reserve( listOfCutsRhs.size() + cuts.size() );

	for( unsigned int i = 0; i < cuts.size(); ++i ){
		vector< data::IndexedElement > addcut;
		addcut.clear();
		int j =0;
		double rrhs=0.0;

		for ( ; j < 20 && i+j < cuts.size();j++) {
			vector< pair<unsigned int, double> >& cut = cuts.at( i+j ).first;
			bool found=false;
			if (computeCutRatio(cut) > 1000000.0) continue;
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
				ii = jj-1;
			}
			// delete 0s
			ii=0, jj=0;
			for (ii = 0; ii < addcut.size();ii++) {
				if (!isZero_(addcut[ii].value.asDouble())) continue;
				jj = addcut.size()-1;
				while (jj>=ii && isZero_(addcut[jj].value.asDouble())) {
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
			if (!(cuttype & GMI)) break;
			if (!found) break;
		}
		i += j;
		j=0;
	    if (addcut.size()>0) {
			listOfCutsLhs.push_back( addcut );
			listOfCutsRhs.push_back( rrhs );
		}

	}

	return listOfCutsRhs.size(); //number of cuts
}

double yInterface::getReducedCostDj(int ix, bool &atLower, bool &atUpper, bool &loUp) {
    atUpper = atLower = false;
    //cerr << "enter getDj with x" << ix << endl;
    if (basis.variables[ix] == extSol::QpExternSolver::AtLower) {
        atLower = true;
    	double d_j = reducedCost[ix].asDouble();
    	if (fabs(d_j) > LP_EPS) {
	  loUp = true;
	  //cerr << "return getDjlo with d_j" << d_j << endl;
	  return d_j;
    	} //else cerr << "getdj Lo d_j=" << d_j << endl;
    } else if (basis.variables[ix] == extSol::QpExternSolver::AtUpper) {
        atUpper = true;
    	double d_j = reducedCost[ix].asDouble();
    	if (fabs(d_j) > LP_EPS) {
	  loUp = true;
	  //cerr << "return getDjup with d_j" << d_j << endl;
	  return d_j;
    	} //else cerr << "getdj up  d_j=" << d_j << endl;
    } //else cerr << "neither lo nor up in getdj" << endl;
    loUp=false;
    return -1.0;
}
void yInterface::getRCandB(extSol::QpExternSolver& externSolver) {
	externSolver.getBase( basis );
	has_basis = true;
	if (basis.variables.size() > 0) {
	  //std::cerr << "+";
	  externSolver.getReducedCosts( reducedCost );
	} else {
	  std::cerr << "Warning: no basis although feasible.";
	  has_basis = false;
	  reducedCost.clear();
	}
}

int yInterface::isReducedCostFixed( double z, double lpval, int ix, double solx, double &rc) {
  //cerr << "enter isRCF with x" << ix << endl;
  if (!has_basis) return 3;
    if (basis.variables[ix] == extSol::QpExternSolver::AtLower) {
    	double d_j = fabs(reducedCost[ix].asDouble());
    	if (d_j > LP_EPS) {
        	double delta = floor(fabs(z-lpval) / d_j);
            if (delta < 1.0-LP_EPS) {
	      // z-lpval < d_j  <=> (z-lpval) / d_j < 1.0
            	//cerr << "V";
	      //cerr << "return isRCFlo with d_j" << d_j << endl;
            	return 0;
            } //else cerr << "pcf lo2 d_j=" << d_j << endl;
    	} //else cerr << "pcf Lo d_j=" << d_j << endl;
    } else if (basis.variables[ix] == extSol::QpExternSolver::AtUpper) {
    	double d_j = fabs(reducedCost[ix].asDouble());
    	if (d_j > LP_EPS) {
        	double delta = floor(fabs(z-lpval) / d_j);
            if (delta < 1.0-LP_EPS) {
            	//cerr << "W";
	      //cerr << "return isRCFup with d_j" << d_j << endl;
            	return 1;
            } //else cerr << "pcf Up2 d_j=" << d_j << endl;
    	} //else cerr << "pcf Up d_j=" << d_j << endl;
    } //else cerr << "neither lo nor up" << endl;
    return 2;
}

void yInterface::setBranchingDecision(int &pick, int &left, int &right) {

}

void yInterface::moveUp(coef_t &value, coef_t uBound, int status) {

}
void yInterface::moveDown(int toLevel, int pick, int val, int val_ix) {

}

std::vector<double> yInterface::getFirstStageSolution(){
	std::vector<double> firstStageSolution;
    qbp->getFirstStageSolution( firstStageSolution );
    return firstStageSolution;
}

#define min(a,b) ((a)<(b)?(a):(b))
int64_t yInterface::ggt(int64_t a, int64_t b) {
	if (a == 0) return b;
	else if (b==0) return a;
	if (a < 0) a = -a;
	if (b < 0) b = -b;
	if (a < b) {
		int64_t t = a; a = b; b = t;
	}
	return ggt(b,(a % b));
}


#define USE_LP_REDUCTION_OUT 0
void yInterface::simplifyCoefficients_01(std::vector<data::IndexedElement> &lhs,std::vector<data::IndexedElement> &org_lhs, data::QpRhs &rhs, data::QpRhs &org_rhs, int8_t *ba) {
	//return;
	struct Comp
	{
	   bool operator()(const std::pair<double,int>& s1, const std::pair<double,int>& s2)
	   {
	       return fabs(s1.first) < fabs(s2.first);
	   }
	};
	typedef std::priority_queue<std::pair<double,int>,std::vector<std::pair<double,int> >, Comp> mypq_type;
	if (rhs.getRatioSign() == data::QpRhs::equal) return;
	if (lhs.size() < 2) return;
	bool flipAll = false;
	if (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) flipAll = true;

	mypq_type myHeap;
	double alpha_bar = 0.0;
	if (USE_LP_REDUCTION_OUT) cerr << "Constraint: ";
	for (int ii = 0; ii < lhs.size();ii++) {
		if (flipAll) lhs[ii].set(lhs[ii].index, -lhs[ii].value.asDouble());
		if (USE_LP_REDUCTION_OUT) cerr << lhs[ii].value.asDouble() << "x" << lhs[ii].index << " ";
		myHeap.push(std::pair<double,int>(lhs[ii].value.asDouble(),ii));
		if (lhs[ii].value.asDouble() > 0) alpha_bar += lhs[ii].value.asDouble();
	}
	if (flipAll) rhs.set(data::QpRhs::smallerThanOrEqual,-rhs.getValue().asDouble());
	if (USE_LP_REDUCTION_OUT) cerr << " <= " << rhs.getValue().asDouble() << " ---- TOP:" << myHeap.top().first << endl;
    assert(rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual);
    bool changed;
    do {
    	changed = false;
    	std::pair<double,int> CVpair; // int index in der Constraint, nicht der Variablenname!
    	CVpair = myHeap.top();
    	if (USE_LP_REDUCTION_OUT) cerr << " alpha_bar= " << alpha_bar << " :: alpha_bar-topX= " << alpha_bar - fabs(CVpair.first) << endl;
    	if (alpha_bar > rhs.getValue().asDouble() && alpha_bar - fabs(CVpair.first) <= rhs.getValue().asDouble()) {
    		double d = alpha_bar - rhs.getValue().asDouble();
    		if (CVpair.first > 0) {
            	myHeap.pop();
    			if (fabs(rhs.getValue().asDouble() - (rhs.getValue().asDouble() - (CVpair.first - d) * 1.0)) > 1e-12) changed = true;
    			if (USE_LP_REDUCTION_OUT) cerr << "change rhs from " << rhs.getValue().asDouble() << " to " << rhs.getValue().asDouble() - (CVpair.first - d) * 1.0 << endl;
    			rhs.setValue(rhs.getValue().asDouble() - (CVpair.first - d) * 1.0);
    			lhs[CVpair.second].set(lhs[CVpair.second].index,d);
    			if (USE_LP_REDUCTION_OUT) cerr << "change x" << lhs[CVpair.second].index << " to " << d << endl;
        		myHeap.push(std::pair<double,int>(d,CVpair.second));
        		if (USE_LP_REDUCTION_OUT) cerr << "change alpha_bar from " << alpha_bar << " to " << alpha_bar-CVpair.first+d << endl;
    			alpha_bar -= CVpair.first;
    			alpha_bar += d;
    		} else if (CVpair.first < 0){
            	myHeap.pop();
            	double old_alpha_bar = alpha_bar;
            	alpha_bar = alpha_bar - lhs[CVpair.second].value.asDouble();
            	double beta_bar = -lhs[CVpair.second].value.asDouble() + rhs.getValue().asDouble();
    			if (fabs(beta_bar - (beta_bar - (-CVpair.first - d) * 1.0)) > 1e-12) changed = true;
    			if (USE_LP_REDUCTION_OUT) cerr << "2)change rhs from " << rhs.getValue().asDouble() << " to " << beta_bar - (-CVpair.first - d) * 1.0 -d << endl;
    			rhs.setValue(beta_bar - (-CVpair.first - d) * 1.0 -d);
    			lhs[CVpair.second].set(lhs[CVpair.second].index,-d);
    			if (USE_LP_REDUCTION_OUT) cerr << "2) change x" << lhs[CVpair.second].index << " to " << -d << endl;
    			myHeap.push(std::pair<double,int>(-d,CVpair.second));
    			alpha_bar = old_alpha_bar;
    			if (USE_LP_REDUCTION_OUT) cerr << "2) change alpha_bar from " << alpha_bar << " to " << alpha_bar/*+CVpair.first+d*/ << endl;
    			//alpha_bar += CVpair.first;
    			//alpha_bar += d;
    		}
    	}
    } while (changed);
}

void yInterface::analyzeAndChangeIntRow(std::vector<data::IndexedElement> &lhs,std::vector<data::IndexedElement> &org_lhs, data::QpRhs &rhs, data::QpRhs &org_rhs, int8_t *ba) {
  int64_t pow_max=(int64_t)1;
  int64_t pow_min=(int64_t)1;
  //int p1 = (int)((pow_max) >> 32);
  //int p2 = (int)((pow_max & 0xffffffffffffffff));
  //cerr << "p1=" << p1 << " p2=" << p2 << " sint=" << sizeof(int) << " sint64=" << sizeof(int64_t) << endl;
	double tmp_rhs = rhs.getValue().asDouble();
	if (0&&lhs.size() < 2) {
		return;
	}

    /*
    if (lhs.size()==2 && (lhs[0].index == 1157 || lhs[1].index == 1157)) cerr << "vor analyze: " << endl;
    if (lhs.size()==2 && lhs[0].index == 1157) for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " : " <= ") << rhs.getValue() << endl;
    if (org_lhs.size()==2 && (org_lhs[0].index == 1157 || org_lhs[1].index == 1157)) cerr << "vor analyze ORG: " << endl;
    if (org_lhs.size()==2 && org_lhs[0].index == 1157) for (int i = 0; i < org_lhs.size();i++) cerr << org_lhs[i].value.asDouble() << "x" << org_lhs[i].index << " ";
    if (org_lhs.size()==2 && org_lhs[0].index == 1157) cerr << (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " :
                                                                (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual ? " <= " : " = " )) << org_rhs.getValue(). asDouble() << " ba[1157]="<< (int)ba[1157] << " ba[1844]="<< (int)ba[1844]<< endl;
     */
    //if (rhs.getRatioSign() == data::QpRhs::equal) return;

	for (int i = 0; i < lhs.size();i++) {
		if (ba[lhs[i].index] == 1 || ba[lhs[i].index] == 0) {
			if (ba[lhs[i].index] == 1) {
				tmp_rhs = tmp_rhs - lhs[i].value.asDouble();
			}
            /*
            bool fits=false;
            for (int ii=0; ii < lhs.size() && ii < 4;ii++) {
                if (lhs[ii].index == 1157) fits=true;
            }
            if (lhs.size() <= 4 && fits) {
              cerr << "index " << i << "muss weg. ba=" << (int)ba[lhs[i].index] << endl;
              for (int ii = 0; ii < lhs.size();ii++) cerr << lhs[ii].value.asDouble() << "x" << lhs[ii].index << " ";
              cerr << (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " :
                     (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual ? " <= " : " = ") ) << rhs.getValue().asDouble() << endl;
            }
             */
			lhs[i] = lhs[lhs.size()-1];
			lhs.pop_back();
			org_lhs[i] = org_lhs[org_lhs.size()-1];
			org_lhs.pop_back();
			i--;
		}
	}
	//rhs.setValue(tmp_rhs);
	//org_rhs.setValue(tmp_rhs);
	//return;
	if(0)if (rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
		rhs.set(data::QpRhs::RatioSign::smallerThanOrEqual, -rhs.getValue().asDouble());
		org_rhs.set(data::QpRhs::RatioSign::smallerThanOrEqual, -org_rhs.getValue().asDouble());
		for (int i = 0; i < lhs.size();i++) {
			lhs[i].value = -lhs[i].value.asDouble();
			org_lhs[i].value = -org_lhs[i].value.asDouble();
		}
		tmp_rhs = rhs.getValue().asDouble();
	}

    /*
    if (lhs.size()==2 && (lhs[0].index == 1157 || lhs[1].index == 1157)) cerr << "vor step2: " << endl;
    if (lhs.size()==2 && lhs[0].index == 1157) for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " :
                                                        (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual ? " <= " : " = ")) << rhs.getValue() << endl;
    if (org_lhs.size()==2 && (org_lhs[0].index == 1157 || org_lhs[1].index == 1157)) cerr << "vor step2 ORG: " << endl;
    if (org_lhs.size()==2 && org_lhs[0].index == 1157) for (int i = 0; i < org_lhs.size();i++) cerr << org_lhs[i].value.asDouble() << "x" << org_lhs[i].index << " ";
    if (org_lhs.size()==2 && org_lhs[0].index == 1157) cerr << (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " :
                                                                (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual ? " <= " : " = " )) << org_rhs.getValue(). asDouble() << " ba[1157]="<< (int)ba[1157] << " ba[1844]="<< (int)ba[1844]<< endl;

    if (lhs.size()==2 && (lhs[0].index == 1157 || lhs[1].index == 1157)) {
        cerr << "tmp_rhs:" << tmp_rhs << " und signs:" << (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " :
                                                            (org_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual ? " <= " : " = " )) <<
                                                          (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " :
                                                           (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual ? " <= " : " = ")) << endl;
    }
     */
    
    rhs.setValue(tmp_rhs);
	org_rhs.setValue(tmp_rhs);

	if (lhs.size() < 2) {
		return;
	}

    /*
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << "vor simplify coeffs: " << endl;
    if (lhs.size()==2 && lhs[0].index == 1157) for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " : " <= ") << rhs.getValue() << endl;
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << "vor simplify coeffs ORG: " << endl;
    if (lhs.size()==2 && lhs[0].index == 1157) for (int i = 0; i < org_lhs.size();i++) cerr << org_lhs[i].value.asDouble() << "x" << org_lhs[i].index << " ";
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " : " <= ") << org_rhs.getValue() << endl;
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << (org_rhs.getRatioSign() == data::QpRhs::equal ? " NO ITS EQual!!" : "") << endl;
     */
	simplifyCoefficients_01(lhs,org_lhs, rhs, org_rhs, ba);
    /*
     if (lhs.size()==2 && lhs[0].index == 1157) cerr << "nach simplify coeffs: " << endl;
    if (lhs.size()==2 && lhs[0].index == 1157) for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " : " <= ") << rhs.getValue() << endl;
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << "nach simplify coeffs ORG: " << endl;
    if (lhs.size()==2 && lhs[0].index == 1157) for (int i = 0; i < org_lhs.size();i++) cerr << org_lhs[i].value.asDouble() << "x" << org_lhs[i].index << " ";
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << (org_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual ? " >= " : " <= ") << org_rhs.getValue() << endl;
    if (lhs.size()==2 && lhs[0].index == 1157) cerr << (org_rhs.getRatioSign() == data::QpRhs::equal ? " NO ITS EQual!!" : "") << endl;
     */

    for (int i = 0; i < lhs.size();i++) {
      //cerr << "coeff=" << lhs[i].value.asDouble() << endl;
      int64_t pow = (int64_t)1;
      //int p1 = (int)((pow) >> 32);
      //int p2 = (int)((pow & 0xffffffffffffffff));
      //cerr << "P1=" << p1 << " P2=" << p2 << endl;
		for (int j=0; j < 32;j++) {
			double coe = floor(fabs(lhs[i].value.asDouble()) * ( (double)pow ) + 1e-12);
			double frac = fabs(lhs[i].value.asDouble()) * ( (double)pow ) - coe;
			if (frac / (double)pow > 1e-12) {
			  //cerr << "cu pmi=" << pow_min << " pma=" << pow_max << " pow=" << pow << " coe=" << coe << " frac=" << frac  << endl;
				pow = pow * 10;
			} else {
				if (pow > pow_max) pow_max = pow;
				else if (pow < pow_min) pow_min = pow;
				//cerr << " pmi=" << pow_min << " pma=" << pow_max << " pow=" << pow << endl;
				break;
			}
			if (j > 25) {
			  //cerr << "could not find correct power. Abort => Warning:" << lhs[i].value.asDouble() << endl;
				return;
			}
		}
	}

	std::vector<int64_t> ilhs(lhs.size());

	for (int i = 0; i < lhs.size(); i++) {
		if (lhs[i].value.asDouble() >= 0)
  		    ilhs[i] = (int64_t)(fabs(lhs[i].value.asDouble()) * ( (double)pow_max ) );
		else
  		    ilhs[i] = -(int64_t)(fabs(lhs[i].value.asDouble()) * ( (double)pow_max ) );
		//cerr << "made " << (double)ilhs[i] << endl;
	}
	tmp_rhs = rhs.getValue().asDouble() * ((double)pow_max);

	int64_t gt = ggt(ilhs[0], ilhs[1]);
	for (int i = 2; i < ilhs.size();i++) gt = ggt(ilhs[i],gt);

	if (gt <= 1) return;
	if (USE_LP_REDUCTION_OUT) cerr << "Preprocess ggt = " << gt << " rhs=" << tmp_rhs / ((double)gt) << endl;

	for (int i = 0; i < ilhs.size();i++) ilhs[i] = ilhs[i] / gt;
	tmp_rhs = tmp_rhs / ((double)gt);

	//double tmp_rhs_frac = fabs(tmp_rhs) - floor(fabs(tmp_rhs));

	if (USE_LP_REDUCTION_OUT) for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
	if (USE_LP_REDUCTION_OUT) cerr << " # " << rhs.getValue() << endl;

    if (rhs.getRatioSign() == data::QpRhs::equal) {
    	//return;

		if (!(fabs((tmp_rhs)-floor((tmp_rhs)+0.5)) <= (1e-9 > 1e-15 * fabs(tmp_rhs) ? 1e-9 : 1e-15 * fabs(tmp_rhs)))/*tmp_rhs_frac > 1e-10*/) {
			cerr << "Preprocess: should be infeasible." << endl;
		}
		return;
		rhs.setValue(/*floor*/(tmp_rhs));
		for (int i = 0; i < ilhs.size();i++)
			lhs[i].set(lhs[i].index,((double)ilhs[i]));
		if (USE_LP_REDUCTION_OUT) for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
		if (USE_LP_REDUCTION_OUT) cerr << " = " << rhs.getValue() << endl;
		return;
	}
	if (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
		rhs.setValue(ceil(tmp_rhs)-max(1e-15,fabs(tmp_rhs)*1e-18));
		//rhs.setValue(tmp_rhs/*-1e-18*/);
		for (int i = 0; i < ilhs.size();i++)
			lhs[i].set(lhs[i].index,((double)ilhs[i]));
		if (USE_LP_REDUCTION_OUT) for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
		if (USE_LP_REDUCTION_OUT) cerr << " >= " << rhs.getValue() << endl;
		return;
	}

	if (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
		rhs.setValue(floor(tmp_rhs)+max(1e-15,fabs(tmp_rhs)*1e-18));
		//rhs.setValue(tmp_rhs/*+1e-18*/);
		for (int i = 0; i < ilhs.size();i++)
			lhs[i].set(lhs[i].index,((double)ilhs[i]));
		if (USE_LP_REDUCTION_OUT) for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
		if (USE_LP_REDUCTION_OUT) cerr << " <= " << rhs.getValue() << endl;
		return;
	}
}

bool yInterface::exactAvail(std::vector<data::IndexedElement> &table_lhs, std::vector<data::IndexedElement> &lhs, data::QpRhs table_rhs, data::QpRhs rhs) {
	bool isExisting1 = true;
	bool isExisting2 = true;
    if (table_lhs.size() != lhs.size()) return false;
    if (lhs.size() == 2 && (lhs[0].index == 1157 || lhs[1].index == 1157)) {
        cerr << "VOR SORT" << endl;
        if (lhs[0].index == 1157 && lhs[1].index == 1844) {
            cerr << "in GENAU VERGLEICH : " << lhs.size() << endl;
            for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
            if (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << rhs.getValue().asDouble() << endl;
        }
        if (table_lhs[0].index == 1157 && table_lhs[1].index == 1844) {
            cerr << "in GENAU VERGLEICH2 : " << table_lhs.size() << endl;
            for (int i = 0; i < table_lhs.size();i++) cerr << table_lhs[i].value.asDouble() << "x" << table_lhs[i].index << " ";
            if (table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << table_rhs.getValue().asDouble() << endl;
        }
        if (lhs[1].index == 1157 && lhs[0].index == 1844) {
            cerr << "in GENAU VERGLEICH : " << lhs.size() << endl;
            for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
            if (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << rhs.getValue().asDouble() << endl;
        }
        if (table_lhs[1].index == 1157 && table_lhs[0].index == 1844) {
            cerr << "in GENAU VERGLEICH2 : " << table_lhs.size() << endl;
            for (int i = 0; i < table_lhs.size();i++) cerr << table_lhs[i].value.asDouble() << "x" << table_lhs[i].index << " ";
            if (table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << table_rhs.getValue().asDouble() << endl;
        }
    }
    std::sort(table_lhs.begin(),table_lhs.end());
    std::sort(lhs.begin(),lhs.end());

    for (int i=0; i < table_lhs.size();i++) {
        if (table_lhs[i].index != lhs[i].index) return false;
    }
    
    if (lhs.size() == 2 && (lhs[0].index == 1157 || lhs[1].index == 1157)) {
        cerr << "NACH SORT" << endl;
        if (lhs[0].index == 1157 && lhs[1].index == 1844) {
            cerr << "in GENAU VERGLEICH : " << lhs.size() << endl;
            for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
            if (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << rhs.getValue().asDouble() << endl;
        }
        if (table_lhs[0].index == 1157 && table_lhs[1].index == 1844) {
            cerr << "in GENAU VERGLEICH2 : " << table_lhs.size() << endl;
            for (int i = 0; i < table_lhs.size();i++) cerr << table_lhs[i].value.asDouble() << "x" << table_lhs[i].index << " ";
            if (table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << table_rhs.getValue().asDouble() << endl;
        }
        if (lhs[1].index == 1157 && lhs[0].index == 1844) {
            cerr << "in GENAU VERGLEICH : " << lhs.size() << endl;
            for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
            if (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << rhs.getValue().asDouble() << endl;
        }
        if (table_lhs[1].index == 1157 && table_lhs[0].index == 1844) {
            cerr << "in GENAU VERGLEICH2 : " << table_lhs.size() << endl;
            for (int i = 0; i < table_lhs.size();i++) cerr << table_lhs[i].value.asDouble() << "x" << table_lhs[i].index << " ";
            if (table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << table_rhs.getValue().asDouble() << endl;
        }
    }

    double factor1 = fabs(1000 / lhs[0].value.asDouble());
    double factor2 = fabs(1000 / table_lhs[0].value.asDouble());
    if (fabs(lhs[0].value.asDouble() - 1000.0) < 1e-12) {
		for (int i=0; i < table_lhs.size();i++) {
			lhs[i].value = lhs[i].value.asDouble() * factor1;
		}
        rhs.setValue(rhs.getValue().asDouble() * factor1);
    }
	if (fabs(table_lhs[0].value.asDouble() - 1000.0) < 1e-12) {
		for (int i=0; i < table_lhs.size();i++) {
			table_lhs[i].value = table_lhs[i].value.asDouble() * factor2;
		}
		table_rhs.setValue(table_rhs.getValue().asDouble() * factor2);
	}

    if (lhs.size() == 2 && (lhs[0].index == 1157 || lhs[1].index == 1157)) {
        cerr << "NACH NORMIERUNG" << endl;
        if (lhs[0].index == 1157 && lhs[1].index == 1844) {
            cerr << "in GENAU VERGLEICH : " << lhs.size() << endl;
            for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
            if (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << rhs.getValue().asDouble() << endl;
        }
        if (table_lhs[0].index == 1157 && table_lhs[1].index == 1844) {
            cerr << "in GENAU VERGLEICH2 : " << table_lhs.size() << endl;
            for (int i = 0; i < table_lhs.size();i++) cerr << table_lhs[i].value.asDouble() << "x" << table_lhs[i].index << " ";
            if (table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << table_rhs.getValue().asDouble() << endl;
        }
        if (lhs[1].index == 1157 && lhs[0].index == 1844) {
            cerr << "in GENAU VERGLEICH : " << lhs.size() << endl;
            for (int i = 0; i < lhs.size();i++) cerr << lhs[i].value.asDouble() << "x" << lhs[i].index << " ";
            if (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< rhs.getValue().asDouble() << endl;
            if (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << rhs.getValue().asDouble() << endl;
        }
        if (table_lhs[1].index == 1157 && table_lhs[0].index == 1844) {
            cerr << "in GENAU VERGLEICH2 : " << table_lhs.size() << endl;
            for (int i = 0; i < table_lhs.size();i++) cerr << table_lhs[i].value.asDouble() << "x" << table_lhs[i].index << " ";
            if (table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< table_rhs.getValue().asDouble() << endl;
            if (table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << table_rhs.getValue().asDouble() << endl;
        }
    }
    

    
    
    for (int i=0; i < table_lhs.size();i++) {
    	if (fabs(table_lhs[i].value.asDouble() - lhs[i].value.asDouble()) > 1e-9) {
    		isExisting1 = false;
    		break;
    	}
    }
    if (isExisting1) {
    	if ((table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual || table_rhs.getRatioSign() == data::QpRhs::equal) && rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual &&
    			table_rhs.getValue().asDouble() <= rhs.getValue().asDouble() ) ;
    	else if ((table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual|| table_rhs.getRatioSign() == data::QpRhs::equal) && rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual &&
    			table_rhs.getValue().asDouble() >= rhs.getValue().asDouble() ) ;
    	else if (fabs(table_rhs.getValue().asDouble() - rhs.getValue().asDouble()) > 1e-9) isExisting1 = false;
        else if (table_rhs.getRatioSign() != rhs.getRatioSign()) isExisting1 = false;
    }
	//if (fabs(table_rhs.getValue().asDouble() - rhs.getValue().asDouble()) > 1e-9) isExisting1 = false;

    for (int i=0; i < table_lhs.size();i++) {
    	if (fabs(-table_lhs[i].value.asDouble() - lhs[i].value.asDouble()) > 1e-9) {
    		isExisting2 = false;
    		break;
    	}
    }
    
    if (isExisting2) { // lhs ist -1 mal table_lhs
        data::QpRhs tmp_rhs;
        // bild artificial rhs
        tmp_rhs.setValue( -rhs.getValue().asDouble()  );
        if (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) tmp_rhs.setRatioSign(data::QpRhs::greaterThanOrEqual);
        else if (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) tmp_rhs.setRatioSign(data::QpRhs::smallerThanOrEqual);
        else tmp_rhs.setRatioSign(data::QpRhs::equal);
    
    if (lhs.size() == 2 && (lhs[0].index == 1157 || lhs[1].index == 1157)) {
        cerr << "CHECK EXISTING: " << isExisting1 << " " << isExisting2 << endl;
    }

        if ((table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual || table_rhs.getRatioSign() == data::QpRhs::equal) && tmp_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual &&
            table_rhs.getValue().asDouble() <= tmp_rhs.getValue().asDouble() ) ;
        else if ((table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual|| table_rhs.getRatioSign() == data::QpRhs::equal) && tmp_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual &&
                 table_rhs.getValue().asDouble() >= tmp_rhs.getValue().asDouble() ) ;
        else if (fabs(table_rhs.getValue().asDouble() - tmp_rhs.getValue().asDouble()) > 1e-9) isExisting2 = false;
        else if (table_rhs.getRatioSign() != tmp_rhs.getRatioSign()) isExisting2 = false;
        /*
    	if ((table_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual || table_rhs.getRatioSign() == data::QpRhs::equal) && rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual &&
    			table_rhs.getValue().asDouble() >= rhs.getValue().asDouble() ) ;
    	else if ((table_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual|| table_rhs.getRatioSign() == data::QpRhs::equal) && rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual &&
    			table_rhs.getValue().asDouble() <= rhs.getValue().asDouble() ) ;
    	else if (fabs(table_rhs.getValue().asDouble() - rhs.getValue().asDouble()) > 1e-9) isExisting2 = false;
         */
    }
	//if (fabs(-table_rhs.getValue().asDouble() - rhs.getValue().asDouble()) > 1e-9) isExisting2 = false;
    if (USE_LP_REDUCTION_OUT) cerr << "year, have found a constraint 1:1 !" << endl;
    if (isExisting1 || isExisting2) {
    	return true;
    	if (table_rhs.getRatioSign() == data::QpRhs::equal && rhs.getRatioSign() == data::QpRhs::equal) return true;
    }
    if (table_rhs.getRatioSign() == data::QpRhs::equal || rhs.getRatioSign() == data::QpRhs::equal) return false;
    if (isExisting1 && table_rhs.getRatioSign() == rhs.getRatioSign()) return true;
    else if (isExisting2 && table_rhs.getRatioSign() != rhs.getRatioSign()) return true;
    return false;
}

static bool doSos=true;
static int SHOW_DETAILS = 0;
//2

bool yInterface::preprocessConstraint(std::vector<int> &virtual_ids, std::vector<data::IndexedElement> &LHS_chg, data::QpRhs &RHS_chg,
		data::Qlp &qmip, utils::QlpStageSolver& QlpStSolve, HCTable * hct, std::pair<coef_t,uint64_t> &hp, int8_t *binary_assignments_, coef_t rhsbnd, int *types,
		int maxLPstage, std::vector<std::pair<int,double> > &cpropQ, coef_t *lowerBounds, coef_t *upperBounds,
		bool feasPhase, std::vector< std::pair<int,int> > &clist, int *block, int *eas) {

    data::QpRhs new_rhs = RHS_chg;
	data::QpRhs org_rhs = RHS_chg;
	std::vector<data::IndexedElement> org_lhs;// = conVec[i]->getElements(); //?
	std::vector<data::IndexedElement> int_lhs;// = conVec[i]->getElements(); //?
	for (int i = 0; i < LHS_chg.size();i++) {
		org_lhs.push_back(LHS_chg[i]);
		int_lhs.push_back(LHS_chg[i]);
	}
	LHS_chg.clear();
	this->QlpStSolve = &QlpStSolve;
	std::vector<int> mon_pos(qmip.getVariableVectorConst().size(),0);
    std::vector<int> mon_neg(qmip.getVariableVectorConst().size(),0);
	std::vector<const data::QpVar *> varVec = qmip.getVariableVectorConst();
    std::vector<int8_t> binary_assignments(varVec.size());
    int cntVars = qmip.getVariableCount(); //qmip.getVariableVectorConst().size();
    std::cerr << "cntVar=" << cntVars << " VIRT.SIZE=" << virtual_ids.size() << std::endl;
    for (int z=0; z < cntVars;z++) {
      assert(binary_assignments_[virtual_ids[z]] == 0 || binary_assignments_[virtual_ids[z]] == 1 || binary_assignments_[virtual_ids[z]] == 2);
    	binary_assignments[z] = (binary_assignments_[virtual_ids[z]]);
    }
    int cntSOS=0;
    SOSvars.clear();
    SOSconstraints.clear();

	if (SHOW_DETAILS > 1) {
		std::cerr << "original constraint:"<< std::endl;
		for (int j = 0; j < org_lhs.size();j++) {
			std::cerr << org_lhs[j].value.asDouble() << "x" << virtual_ids[org_lhs[j].index] << " + ";
		}
		std::cerr << (org_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
					  (org_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << org_rhs.getValue() << std::endl;
	}

	double rhs_offset=0.0;
	for (int j = 0; j < org_lhs.size();j++) {
        if (virtual_ids[org_lhs[j].index] <0 || virtual_ids[org_lhs[j].index] >= cntVars) virtual_ids[org_lhs[j].index] = org_lhs[j].index;
		if (fabs(qmip.getVariableByIndex(virtual_ids[org_lhs[j].index]).getUpperBound().asDouble()-qmip.getVariableByIndex(virtual_ids[org_lhs[j].index]).getLowerBound().asDouble()) < 0.0000001) {
			rhs_offset = rhs_offset + org_lhs[j].value.asDouble() * ((qmip.getVariableByIndex(virtual_ids[org_lhs[j].index]).getUpperBound().asDouble()+qmip.getVariableByIndex(virtual_ids[org_lhs[j].index]).getLowerBound().asDouble())/2);
			org_lhs[j].value = 0.0;
			int_lhs[j].value = 0.0;
		}
	}
	org_rhs.setValue(org_rhs.getValue().asDouble()-rhs_offset);
	new_rhs.setValue(new_rhs.getValue().asDouble()-rhs_offset);
	if (SHOW_DETAILS > 1) {
		std::cerr << "original constraint 2:"<< std::endl;
		for (int j = 0; j < org_lhs.size();j++) {
			std::cerr << org_lhs[j].value.asDouble() << "x" << virtual_ids[org_lhs[j].index] << " + ";
		}
		std::cerr << (org_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
					  (org_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << org_rhs.getValue() << std::endl;
	}

	double lhs_min = 0.0, lhs_max = 0.0;
	for (int j = 0; j < org_lhs.size();j++) {
		if (types[virtual_ids[org_lhs[j].index]] == 0/*BINARY*/ && binary_assignments[virtual_ids[org_lhs[j].index]] == 0) {
			;
		} else if (types[virtual_ids[org_lhs[j].index]] == 0/*BINARY*/ && binary_assignments[virtual_ids[org_lhs[j].index]] == 1) {
			lhs_min += org_lhs[j].value.asDouble();
			lhs_max += org_lhs[j].value.asDouble();
		} else if (types[virtual_ids[org_lhs[j].index]] == 0/*BINARY*/ && org_lhs[j].value.asDouble() < 0.0) {
        	lhs_min += org_lhs[j].value.asDouble();
        } else if (types[virtual_ids[org_lhs[j].index]] == 0/*BINARY*/ && org_lhs[j].value.asDouble() >= 0.0) {
        	lhs_max += org_lhs[j].value.asDouble();
        } else if (org_lhs[j].value.asDouble() < 0.0) {
        	lhs_min += org_lhs[j].value.asDouble()*qmip.getVariableByIndex(virtual_ids[org_lhs[j].index]).getUpperBound().asDouble();
        	lhs_max += org_lhs[j].value.asDouble()*qmip.getVariableByIndex(virtual_ids[org_lhs[j].index]).getLowerBound().asDouble();
        } else if (org_lhs[j].value.asDouble() >= 0.0) {
        	lhs_min += org_lhs[j].value.asDouble()*qmip.getVariableByIndex(virtual_ids[org_lhs[j].index]).getLowerBound().asDouble();
        	lhs_max += org_lhs[j].value.asDouble()*qmip.getVariableByIndex(virtual_ids[org_lhs[j].index]).getUpperBound().asDouble();
        }
	}
	if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
		if (lhs_min > new_rhs.getValue().asDouble()) return false;
	} else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
		if (lhs_max < new_rhs.getValue().asDouble()) return false;
	}
	// does the constraint contain only integer variables?
	bool contReal = false;
	bool isSOS=true;
	isSOS = false;
	int numnegs=0;
	for (int j = 0; j < org_lhs.size();j++) {
		int old_index = virtual_ids[org_lhs[j].index];
		if (types[old_index] == 5000) {
			contReal = true;
			break;
		}
	}
	// if yes, make coeffs. integer and use gcd to strengthen.
	if (contReal == false) {
		if (isSOS==true && new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal
				&& new_rhs.getValue().asDouble() <= 1.0 - numnegs + LP_EPS && new_rhs.getValue().asDouble() >= 1.0 - numnegs - LP_EPS) {
			if (doSos) isSOS = extractSOS(int_lhs, new_rhs, cntVars);
			else isSOS = false;
		    //std::cerr << "EX!" << std::endl;
			cntSOS++;
		} else isSOS = false;//else  std::cerr << "no EX!" << std::endl;

		if (USE_LP_REDUCTION_OUT) cerr << "old size=" << org_lhs.size();
		if (isSOS==false) analyzeAndChangeIntRow(int_lhs, org_lhs, new_rhs, org_rhs, binary_assignments.data());
		if (USE_LP_REDUCTION_OUT) cerr << "  new size=" << org_lhs.size() << endl;
		hp = hct->computeHash(int_lhs, new_rhs.getValue().asDouble(), new_rhs.getRatioSign());
		if (USE_LP_REDUCTION_OUT) cerr << "hash=" << hp.second << ", rhs:" << hp.first << endl;
		if (USE_LP_REDUCTION_OUT && org_lhs.size() <= 0) {
			cerr << "size==0?:" << (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) << "," << new_rhs.getValue().asDouble() << ";"
								<< (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) << "," << new_rhs.getValue().asDouble() << ";"
								<< (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal) << "," <<new_rhs.getValue().asDouble() << endl;
		}
		HTCutentry *htce;
		if (hct->getEntry(&htce, hp.second, hp.first, new_rhs.getRatioSign()) == true /*&& htce->index < LHSs.size() && exactAvail(LHSs[htce->index], int_lhs, RHSs[htce->index], new_rhs)*/) {
			if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint" << endl;
            cerr << "E";
			return false;
		} else {
			int cnt_triv = 0;
			coef_t lhssum=0.0;
			if (int_lhs.size() == 1) {
				if (fabs(int_lhs[0].value.asDouble()) < 1e-12 ) int_lhs.pop_back();
				else {
					if (int_lhs[0].value.asDouble() < 0.0) {
						if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual)
							new_rhs.setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
						else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual)
							new_rhs.setRatioSign(data::QpRhs::RatioSign::greaterThanOrEqual);
					}
					new_rhs.setValue(new_rhs.getValue() / int_lhs[0].value.asDouble());
					int_lhs[0].value = 1.0;
				}
			}

			if (int_lhs.size() == 1) {
				if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && fabs(new_rhs.getValue().asDouble()) <= 1e-12) {
					binary_assignments[virtual_ids[int_lhs[0].index]] = 0;
					cpropQ.push_back(std::pair<int,double>(virtual_ids[int_lhs[0].index],0.0));
					cnt_triv++;
				} else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual && fabs(1.0-new_rhs.getValue().asDouble()) <= 1e-12) {
					binary_assignments[virtual_ids[int_lhs[0].index]] = 1;
					cpropQ.push_back(std::pair<int,double>(virtual_ids[int_lhs[0].index],1.0));
					cnt_triv++;
				} else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal && (fabs(1.0-new_rhs.getValue().asDouble()) <= 1e-12 || fabs(new_rhs.getValue().asDouble()) <= 1e-12)) {
					binary_assignments[virtual_ids[int_lhs[0].index]] = (int)(0.5+new_rhs.getValue().asDouble());
					cpropQ.push_back(std::pair<int,double>(virtual_ids[int_lhs[0].index],new_rhs.getValue().asDouble()));
					cnt_triv++;
				}
			}
			//if (cnt_triv > 0) continue;

			if (int_lhs.size() <= 0 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual && new_rhs.getValue().asDouble() <= 0.0) {
				if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint II" << endl;
				return false;
			} else if (0&&int_lhs.size() <= 1 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && new_rhs.getValue().asDouble() >= 1.0) {
				if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint III" << endl;
				return false;
			} else if (int_lhs.size() <= 0 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && new_rhs.getValue().asDouble() >= 0.0) {
				if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint III" << endl;
				return false;
			} else if (int_lhs.size() <= 0 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal && fabs(new_rhs.getValue().asDouble()) <= 1e-12) {
				if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint IV" << endl;
					return false;
			} else if (int_lhs.size() <= 1 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && new_rhs.getValue().asDouble() >= 1.0) {
				if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint V" << endl;
					return false;
			} else if (int_lhs.size() <= 1 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual && new_rhs.getValue().asDouble() <= 0.0) {
				if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint VI" << endl;
					return false;
			} //else hct->setEntry(hp.first, hp.second, LHSs.size());
			//if (cnt_triv > 0) cerr << "T";
		}
		for (int j = 0; j < int_lhs.size();j++) {
			double ce = fabs(int_lhs[j].value.asDouble());
			if (ce < 1.0-LP_EPS || ce > 1.0+LP_EPS) isSOS=false;
			if (int_lhs[j].value.asDouble() < 0.0) numnegs++;
		}
		if (int_lhs.size() == 0 || new_rhs.getRatioSign() != data::QpRhs::RatioSign::equal ||
				new_rhs.getValue().asDouble() > 1.0 - numnegs + LP_EPS || new_rhs.getValue().asDouble() < 1.0 - numnegs - LP_EPS) isSOS = false;
        //cerr << "Kennzahlen:" << new_rhs.getValue().asDouble() << "," << 1.0-numnegs << "," << (int)(new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal) << endl;

		//assert(int_lhs.size() > 0);
		for (int zz=0; zz < int_lhs.size();zz++) {
			  //if (int_lhs[zz].index == 272) is272 = true;
			  if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && int_lhs[zz].value.asDouble() >= 0) mon_neg[virtual_ids[int_lhs[zz].index]]++;
			  else if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && int_lhs[zz].value.asDouble() < 0) mon_pos[virtual_ids[int_lhs[zz].index]]++;
			  else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && int_lhs[zz].value.asDouble() >= 0) mon_pos[virtual_ids[int_lhs[zz].index]]++;
			  else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && int_lhs[zz].value.asDouble() < 0) mon_neg[virtual_ids[int_lhs[zz].index]]++;
			  else if (new_rhs.getRatioSign()==data::QpRhs::equal) { mon_pos[virtual_ids[int_lhs[zz].index]]++;mon_neg[virtual_ids[int_lhs[zz].index]]++;}
		}
		/*if (is272) {
 	        for (int f=0;f<int_lhs.size();f++) {
 	        	cerr << int_lhs[f].value.asDouble() << "x" << int_lhs[f].index << "(" << (int)binary_assignments[int_lhs[f].index] << ")" << " + ";
 	        }
 	        if (new_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr << " <= ";
 	        if (new_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= ";
 	        if (new_rhs.getRatioSign() == data::QpRhs::equal) cerr << " == ";
 	        cerr << new_rhs.getValue().asDouble() << endl;
		}*/

		bool isTrivFullfilled = false;
		if (new_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
			double row_bound = 0.0;
			for (int zz=0; zz < int_lhs.size();zz++) {
				if (int_lhs[zz].value.asDouble() > 0) row_bound = row_bound + int_lhs[zz].value.asDouble();
			}
			//sleep(40);
			//cerr << row_bound << endl;
			//sleep(4);
			//cerr << new_rhs.getValue().asDouble() << endl;
			//sleep(4);
			if (row_bound <= new_rhs.getValue().asDouble()) isTrivFullfilled = true;
		} else if (new_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
			double row_bound = 0.0;
			for (int zz=0; zz < int_lhs.size();zz++) {
				if (int_lhs[zz].value.asDouble() < 0) row_bound = row_bound + int_lhs[zz].value.asDouble();
			}
			if (row_bound >= new_rhs.getValue().asDouble()) isTrivFullfilled = true;
		}

		if (!isSOS && !isTrivFullfilled && int_lhs.size() > 0) {
			//LHSs.push_back(/*org_lhs*/int_lhs);
			//RHSs.push_back(/*org_rhs*/new_rhs);
			//assert(LHSs[LHSs.size()-1].size() <= cntVars);
			assert(int_lhs.size() <= cntVars);
			LHS_chg.clear();
			for (int zz = 0; zz < int_lhs.size();zz++) {
				LHS_chg.push_back(int_lhs[zz]);
			}
			RHS_chg = new_rhs;
        } else return false;
		if (SHOW_DETAILS > 1) {
			std::cerr << "preprocessed1 constraint:"<< std::endl;
			for (int j = 0; j < int_lhs.size();j++) {
				std::cerr << int_lhs[j].value.asDouble() << "x" << int_lhs[j].index << " + ";
			}
			std::cerr << (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << new_rhs.getValue() << std::endl;
		}
	} else {
		isSOS = false;
		for (int zz=0; zz < org_lhs.size();zz++) {
			  //if (org_lhs[zz].index == 272) is272 = true;
			  if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && org_lhs[zz].value.asDouble() >= 0) mon_neg[virtual_ids[org_lhs[zz].index]]++;
			  else if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && org_lhs[zz].value.asDouble() < 0) mon_pos[virtual_ids[org_lhs[zz].index]]++;
			  else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && org_lhs[zz].value.asDouble() >= 0) mon_pos[virtual_ids[org_lhs[zz].index]]++;
			  else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && org_lhs[zz].value.asDouble() < 0) mon_neg[virtual_ids[org_lhs[zz].index]]++;
			  else if (new_rhs.getRatioSign()==data::QpRhs::equal) { mon_pos[virtual_ids[org_lhs[zz].index]]++;mon_neg[virtual_ids[org_lhs[zz].index]]++;}
		}
		if (!isSOS && org_lhs.size() > 0) {
			//LHSs.push_back(org_lhs);
			//RHSs.push_back(new_rhs);
			assert(org_lhs.size() <= cntVars);
			LHS_chg.clear();
			for (int zz = 0; zz < int_lhs.size();zz++) {
				LHS_chg.push_back(int_lhs[zz]);
			}
			RHS_chg = new_rhs;
        } else return false;
		if (SHOW_DETAILS > 1) {
			std::cerr << "preprocessed2 constraint:"<< std::endl;
			for (int j = 0; j < org_lhs.size();j++) {
				std::cerr << org_lhs[j].value.asDouble() << "x" << virtual_ids[org_lhs[j].index] << " + ";
			}
			std::cerr << (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << new_rhs.getValue() << std::endl;
		}
		assert(LHS_chg.size() <= cntVars);
	}
    if (LHS_chg.size() > 10*log2(cntVars)) return false;
	return true;
}

void yInterface::updateConstraints(data::Qlp &qmip, utils::QlpStageSolver& QlpStSolve, int8_t *binary_assignments_, coef_t rhsbnd, int *types,
		int maxLPstage, std::vector<std::pair<int,double> > &cpropQ, coef_t *lowerBounds, coef_t *upperBounds, bool feasPhase, std::vector< std::pair<int,int> > &clist, int *block, int *eas) {
#ifndef FIND_BUG
	//return;
#endif

	static bool oo=true;
	if (0&&!oo) return;
	else oo=false;
	this->QlpStSolve = &QlpStSolve;
	std::vector<int> mon_pos(qmip.getVariableVectorConst().size(),0);
    std::vector<int> mon_neg(qmip.getVariableVectorConst().size(),0);
	std::vector<const data::QpVar *> varVec = qmip.getVariableVectorConst();
    std::vector<int8_t> binary_assignments(varVec.size());
    //std::memcpy((void*)binary_assignments.data(),(void*)binary_assignments_,varVec.size());
    int cntVars = qmip.getVariableCount(); //qmip.getVariableVectorConst().size();
    for (int z=0; z < cntVars;z++)
    	binary_assignments[z] = (binary_assignments_[z]);
    int cntSOS=0;
    SOSvars.clear();
    SOSconstraints.clear();
    //sortRows(maxLPstage);

	for (int hhh=0; hhh < 1; hhh++) {
	    cntSOS=0;
	    SOStabus.clear();
		//if (7*QlpStSolve.getExternSolver(maxLPstage).getRowCount() < varVec.size()) return;
		HCTable * hct = new HCTable(qmip.getVariableVectorConst().size(), 5*qmip.getConstraintCount());
		int numConstraints = qmip.getConstraintCount();
		std::vector<const data::QpRhs *> rhsVec = qmip.getRhsVecConst();
		std::vector<const data::Constraint *> conVec = qmip.getConstraintVecConst();
		std::vector<data::QpRhs> RHSs;
		std::vector< std::vector<data::IndexedElement> > LHSs;
		if (info_level >= 2) cerr << "BEGINN get Constraint count = " << numConstraints << endl;
		for (int i = 0; i < numConstraints;i++) {
			data::QpRhs new_rhs = *rhsVec[i];
			data::QpRhs org_rhs = *rhsVec[i];
			std::vector<data::IndexedElement> org_lhs = conVec[i]->getElements(); //?
			std::vector<data::IndexedElement> int_lhs = conVec[i]->getElements(); //?
            
            if (org_rhs.getResponsibility() == data::QpRhs::UNIVERSAL) continue;

			if (SHOW_DETAILS > 1) {
				std::cerr << "original constraint:"<< std::endl;
				for (int j = 0; j < org_lhs.size();j++) {
					std::cerr << org_lhs[j].value.asDouble() << "x" << org_lhs[j].index << " + ";
				}
				std::cerr << (org_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
							  (org_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << org_rhs.getValue() << std::endl;
			}

			double rhs_offset=0.0;
			for (int j = 0; j < org_lhs.size();j++) {
				if (fabs(qmip.getVariableByIndex(org_lhs[j].index).getUpperBound().asDouble()-qmip.getVariableByIndex(org_lhs[j].index).getLowerBound().asDouble()) < 0.0000001) {
					rhs_offset = rhs_offset + org_lhs[j].value.asDouble() * ((qmip.getVariableByIndex(org_lhs[j].index).getUpperBound().asDouble()+qmip.getVariableByIndex(org_lhs[j].index).getLowerBound().asDouble())/2);
					org_lhs[j].value = 0.0;
					int_lhs[j].value = 0.0;
				}
			}
			org_rhs.setValue(org_rhs.getValue().asDouble()-rhs_offset);
			new_rhs.setValue(new_rhs.getValue().asDouble()-rhs_offset);
			if (SHOW_DETAILS > 1) {
				std::cerr << "original constraint 2:"<< std::endl;
				for (int j = 0; j < org_lhs.size();j++) {
					std::cerr << org_lhs[j].value.asDouble() << "x" << org_lhs[j].index << " + ";
				}
				std::cerr << (org_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
							  (org_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << org_rhs.getValue() << std::endl;
			}

			double lhs_min = 0.0, lhs_max = 0.0;
			for (int j = 0; j < org_lhs.size();j++) {
				if (types[org_lhs[j].index] == 0/*BINARY*/ && binary_assignments[org_lhs[j].index] == 0) {
					;
				} else if (types[org_lhs[j].index] == 0/*BINARY*/ && binary_assignments[org_lhs[j].index] == 1) {
					lhs_min += org_lhs[j].value.asDouble();
					lhs_max += org_lhs[j].value.asDouble();
				} else if (types[org_lhs[j].index] == 0/*BINARY*/ && org_lhs[j].value.asDouble() < 0.0) {
                	lhs_min += org_lhs[j].value.asDouble();
                } else if (types[org_lhs[j].index] == 0/*BINARY*/ && org_lhs[j].value.asDouble() >= 0.0) {
                	lhs_max += org_lhs[j].value.asDouble();
                } else if (org_lhs[j].value.asDouble() < 0.0) {
                	lhs_min += org_lhs[j].value.asDouble()*qmip.getVariableByIndex(org_lhs[j].index).getUpperBound().asDouble();
                	lhs_max += org_lhs[j].value.asDouble()*qmip.getVariableByIndex(org_lhs[j].index).getLowerBound().asDouble();
                } else if (org_lhs[j].value.asDouble() >= 0.0) {
                	lhs_min += org_lhs[j].value.asDouble()*qmip.getVariableByIndex(org_lhs[j].index).getLowerBound().asDouble();
                	lhs_max += org_lhs[j].value.asDouble()*qmip.getVariableByIndex(org_lhs[j].index).getUpperBound().asDouble();
                }
			}
			if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) {
				if (lhs_min > new_rhs.getValue().asDouble()) continue;
			} else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) {
				if (lhs_max < new_rhs.getValue().asDouble()) continue;
			}
			// does the constraint contain only integer variables?
			bool contReal = false;
			bool isSOS=true;
			isSOS = false;
			int numnegs=0;
			for (int j = 0; j < org_lhs.size();j++) {
				int old_index = org_lhs[j].index;
				if (types[old_index] == 5000) {
					contReal = true;
					break;
				}
			}
			// if yes, make coeffs. integer and use gcd to strengthen.
			if (contReal == false) {
				if (isSOS==true && new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal
						&& new_rhs.getValue().asDouble() <= 1.0 - numnegs + LP_EPS && new_rhs.getValue().asDouble() >= 1.0 - numnegs - LP_EPS) {
					if (doSos) isSOS = extractSOS(int_lhs, new_rhs, cntVars);
					else isSOS = false;
				    //std::cerr << "EX!" << std::endl;
					cntSOS++;
				} else isSOS = false;//else  std::cerr << "no EX!" << std::endl;

				if (USE_LP_REDUCTION_OUT) cerr << "old size=" << org_lhs.size();
                if (0&&org_lhs.size()==2) {
                    if (int_lhs.size()!=2) cerr << "size unklar." << endl;
                    else {
                      if (int_lhs[0].index == 1157 && int_lhs[1].index == 1844) {
                        if (int_lhs[0].value.asDouble() < 0) cerr << "VOR ANALYSE int_lhs: 1157 < 0" << endl;
                      }
                      if (int_lhs[1].index == 1157 && int_lhs[0].index == 1844) {
                        if (int_lhs[1].value.asDouble() < 0) cerr << "VOR ANALYSE int_lhs2: 1157 < 0" << endl;
                      }
                      if (org_lhs[0].index == 1157 && org_lhs[1].index == 1844) {
                        if (org_lhs[0].value.asDouble() < 0) cerr << "VOR ANALYSE org_lhs: 1157 < 0" << endl;
                      }
                      if (org_lhs[1].index == 1157 && org_lhs[0].index == 1844) {
                        if (org_lhs[1].value.asDouble() < 0) cerr << "VOR ANALYSE org_lhs2: 1157 < 0" << endl;
                      }
                    }
                }

				if (isSOS==false) analyzeAndChangeIntRow(int_lhs, org_lhs, new_rhs, org_rhs, binary_assignments.data());
                if (0&&int_lhs.size()==2) {
                    if (int_lhs[0].index == 1157 && int_lhs[1].index == 1844) {
                        if (int_lhs[0].value.asDouble() < 0) cerr << "NACH ANALYSE int_lhs: 1157 < 0" << endl;
                    }
                    if (int_lhs[1].index == 1157 && int_lhs[0].index == 1844) {
                        if (int_lhs[1].value.asDouble() < 0) cerr << "NACH ANALYSE int_lhs2: 1157 < 0" << endl;
                    }
                    if (org_lhs[0].index == 1157 && org_lhs[1].index == 1844) {
                        if (org_lhs[0].value.asDouble() < 0) cerr << "NACH ANALYSE org_lhs: 1157 < 0" << endl;
                    }
                    if (org_lhs[1].index == 1157 && org_lhs[0].index == 1844) {
                        if (org_lhs[1].value.asDouble() < 0) cerr << "NACH ANALYSE org_lhs2: 1157 < 0" << endl;
                    }
                    //cerr << "sizes:" << int_lhs.size() << " " << org_lhs.size() << endl;
                }

				if (USE_LP_REDUCTION_OUT) cerr << "  new size=" << org_lhs.size() << endl;
				std::pair<coef_t,uint64_t> hp = hct->computeHash(int_lhs, new_rhs.getValue().asDouble(), new_rhs.getRatioSign());
				if (USE_LP_REDUCTION_OUT) cerr << "hash=" << hp.second << ", rhs:" << hp.first << endl;
				if (USE_LP_REDUCTION_OUT && org_lhs.size() <= 0) {
					cerr << "size==0?:" << (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual) << "," << new_rhs.getValue().asDouble() << ";"
										<< (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual) << "," << new_rhs.getValue().asDouble() << ";"
										<< (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal) << "," <<new_rhs.getValue().asDouble() << endl;
				}
				HTCutentry *htce;
				if (hct->getEntry(&htce, hp.second, hp.first, new_rhs.getRatioSign()) == true && htce->index < LHSs.size() && exactAvail(LHSs[htce->index], int_lhs, RHSs[htce->index], new_rhs)) {
					if (/*USE_LP_REDUCTION_OUT*/int_lhs.size()==2 && int_lhs[0].index == 1157 && int_lhs[1].index == 1844) cerr << "excluded constraint" << endl;
                    int_lhs.clear();
                    new_rhs.setValue(0.0);
					continue;
				} else {
					int cnt_triv = 0;
					coef_t lhssum=0.0;
					if (int_lhs.size() == 1) {
						if (fabs(int_lhs[0].value.asDouble()) < 1e-12 ) int_lhs.pop_back();
						else {
							if (int_lhs[0].value.asDouble() < 0.0) {
								if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual)
									new_rhs.setRatioSign(data::QpRhs::RatioSign::smallerThanOrEqual);
								else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual)
									new_rhs.setRatioSign(data::QpRhs::RatioSign::greaterThanOrEqual);
							}
							new_rhs.setValue(new_rhs.getValue() / int_lhs[0].value.asDouble());
							int_lhs[0].value = 1.0;
						}
					}

					if (int_lhs.size() == 1) {
						if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && fabs(new_rhs.getValue().asDouble()) <= 1e-12) {
							binary_assignments[int_lhs[0].index] = 0;
							cpropQ.push_back(std::pair<int,double>(int_lhs[0].index,0.0));
							cnt_triv++;
						} else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual && fabs(1.0-new_rhs.getValue().asDouble()) <= 1e-12) {
							binary_assignments[int_lhs[0].index] = 1;
							cpropQ.push_back(std::pair<int,double>(int_lhs[0].index,1.0));
							cnt_triv++;
						} else if (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal && (fabs(1.0-new_rhs.getValue().asDouble()) <= 1e-12 || fabs(new_rhs.getValue().asDouble()) <= 1e-12)) {
							binary_assignments[int_lhs[0].index] = (int)(0.5+new_rhs.getValue().asDouble());
							cpropQ.push_back(std::pair<int,double>(int_lhs[0].index,new_rhs.getValue().asDouble()));
							cnt_triv++;
						}
					}
					//if (cnt_triv > 0) continue;

					if (int_lhs.size() <= 0 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual && new_rhs.getValue().asDouble() <= 0.0) {
						if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint II" << endl;
						continue;
					} else if (0&&int_lhs.size() <= 1 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && new_rhs.getValue().asDouble() >= 1.0) {
						if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint III" << endl;
						continue;
					} else if (int_lhs.size() <= 0 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && new_rhs.getValue().asDouble() >= 0.0) {
						if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint III" << endl;
						continue;
					} else if (int_lhs.size() <= 0 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal && fabs(new_rhs.getValue().asDouble()) <= 1e-12) {
						if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint IV" << endl;
							continue;
					} else if (int_lhs.size() <= 1 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && new_rhs.getValue().asDouble() >= 1.0) {
						if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint V" << endl;
							continue;
					} else if (int_lhs.size() <= 1 && new_rhs.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual && new_rhs.getValue().asDouble() <= 0.0) {
						if (USE_LP_REDUCTION_OUT) cerr << "excluded constraint VI" << endl;
							continue;
					} else hct->setEntry(hp.first, hp.second, LHSs.size());
					//if (cnt_triv > 0) cerr << "T";
				}
				for (int j = 0; j < int_lhs.size();j++) {
					double ce = fabs(int_lhs[j].value.asDouble());
					if (ce < 1.0-LP_EPS || ce > 1.0+LP_EPS) isSOS=false;
					if (int_lhs[j].value.asDouble() < 0.0) numnegs++;
				}
				if (int_lhs.size() == 0 || new_rhs.getRatioSign() != data::QpRhs::RatioSign::equal ||
						new_rhs.getValue().asDouble() > 1.0 - numnegs + LP_EPS || new_rhs.getValue().asDouble() < 1.0 - numnegs - LP_EPS) isSOS = false;
                //cerr << "Kennzahlen:" << new_rhs.getValue().asDouble() << "," << 1.0-numnegs << "," << (int)(new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal) << endl;

				//assert(int_lhs.size() > 0);
				for (int zz=0; zz < int_lhs.size();zz++) {
					  //if (int_lhs[zz].index == 272) is272 = true;
					  if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && int_lhs[zz].value.asDouble() >= 0) mon_neg[int_lhs[zz].index]++;
					  else if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && int_lhs[zz].value.asDouble() < 0) mon_pos[int_lhs[zz].index]++;
					  else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && int_lhs[zz].value.asDouble() >= 0) mon_pos[int_lhs[zz].index]++;
					  else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && int_lhs[zz].value.asDouble() < 0) mon_neg[int_lhs[zz].index]++;
					  else if (new_rhs.getRatioSign()==data::QpRhs::equal) { mon_pos[int_lhs[zz].index]++;mon_neg[int_lhs[zz].index]++;}
				}
				/*if (is272) {
	     	        for (int f=0;f<int_lhs.size();f++) {
	     	        	cerr << int_lhs[f].value.asDouble() << "x" << int_lhs[f].index << "(" << (int)binary_assignments[int_lhs[f].index] << ")" << " + ";
	     	        }
	     	        if (new_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr << " <= ";
	     	        if (new_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= ";
	     	        if (new_rhs.getRatioSign() == data::QpRhs::equal) cerr << " == ";
	     	        cerr << new_rhs.getValue().asDouble() << endl;
				}*/

				bool isTrivFullfilled = false;
				if (new_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
					double row_bound = 0.0;
					for (int zz=0; zz < int_lhs.size();zz++) {
						if (int_lhs[zz].value.asDouble() > 0) row_bound = row_bound + int_lhs[zz].value.asDouble();
					}
					//sleep(40);
					//cerr << row_bound << endl;
					//sleep(4);
					//cerr << new_rhs.getValue().asDouble() << endl;
					//sleep(4);
					if (row_bound <= new_rhs.getValue().asDouble()) isTrivFullfilled = true;
				} else if (new_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
					double row_bound = 0.0;
					for (int zz=0; zz < int_lhs.size();zz++) {
						if (int_lhs[zz].value.asDouble() < 0) row_bound = row_bound + int_lhs[zz].value.asDouble();
					}
					if (row_bound >= new_rhs.getValue().asDouble()) isTrivFullfilled = true;
				}

				if (!isSOS && !isTrivFullfilled && int_lhs.size() > 0) {
                    //cerr << "%";
					LHSs.push_back(/*org_lhs*/int_lhs);
					RHSs.push_back(/*org_rhs*/new_rhs);
					assert(LHSs[LHSs.size()-1].size() <= cntVars);
                    /*
                    if (int_lhs.size()==2) {
                        if (int_lhs[0].index == 1157 && int_lhs[1].index == 1844) {
                            if (int_lhs[0].value.asDouble() < 0) cerr << "int_lhs: 1157 < 0" << endl;
                        }
                        if (int_lhs[1].index == 1157 && int_lhs[0].index == 1844) {
                            if (int_lhs[1].value.asDouble() < 0) cerr << "int_lhs2: 1157 < 0" << endl;
                        }
                        if (org_lhs[0].index == 1157 && org_lhs[1].index == 1844) {
                            if (org_lhs[0].value.asDouble() < 0) cerr << "org_lhs: 1157 < 0" << endl;
                            cerr << "vor analyze lhs mit : " << int_lhs.size() << endl;
                            for (int i = 0; i < int_lhs.size();i++) cerr << int_lhs[i].value.asDouble() << "x" << int_lhs[i].index << " ";
                            if (new_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) cerr << " >= " << new_rhs.getValue() << endl;
                            if (new_rhs.getRatioSign() == data::QpRhs::equal) cerr << " = "<< new_rhs.getValue() << endl;
                            if (new_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) cerr <<  " <= " << new_rhs.getValue() << endl;
                            cerr << "index ist " << RHSs.size()-1 << endl;
                        }
                        if (org_lhs[1].index == 1157 && org_lhs[0].index == 1844) {
                            if (org_lhs[1].value.asDouble() < 0) cerr << "org_lhs2: 1157 < 0" << endl;
                        }

                    }
                     */
				} //else cerr << "?";
				if (SHOW_DETAILS > 1) {
					std::cerr << "preprocessed1 constraint:"<< std::endl;
					for (int j = 0; j < int_lhs.size();j++) {
						std::cerr << int_lhs[j].value.asDouble() << "x" << int_lhs[j].index << " + ";
					}
					std::cerr << (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
								  (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << new_rhs.getValue() << std::endl;
				}
			} else {
                //assert(0);
				isSOS = false;
				for (int zz=0; zz < org_lhs.size();zz++) {
					  //if (org_lhs[zz].index == 272) is272 = true;
					  if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && org_lhs[zz].value.asDouble() >= 0) mon_neg[org_lhs[zz].index]++;
					  else if (new_rhs.getRatioSign()==data::QpRhs::smallerThanOrEqual && org_lhs[zz].value.asDouble() < 0) mon_pos[org_lhs[zz].index]++;
					  else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && org_lhs[zz].value.asDouble() >= 0) mon_pos[org_lhs[zz].index]++;
					  else if (new_rhs.getRatioSign()==data::QpRhs::greaterThanOrEqual && org_lhs[zz].value.asDouble() < 0) mon_neg[org_lhs[zz].index]++;
					  else if (new_rhs.getRatioSign()==data::QpRhs::equal) { mon_pos[org_lhs[zz].index]++;mon_neg[org_lhs[zz].index]++;}
				}
				if (!isSOS && org_lhs.size() > 0) {
					LHSs.push_back(org_lhs);
					RHSs.push_back(new_rhs);
				}
				if (SHOW_DETAILS > 1) {
					std::cerr << "preprocessed2 constraint:"<< std::endl;
					for (int j = 0; j < org_lhs.size();j++) {
						std::cerr << org_lhs[j].value.asDouble() << "x" << org_lhs[j].index << " + ";
					}
					std::cerr << (new_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
								  (new_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << new_rhs.getValue() << std::endl;
				}
				assert(LHSs[LHSs.size()-1].size() <= cntVars);
			}
		}
        //cerr << "END LOOP" << endl;
#define DO_CHS
#ifdef DO_CHS
		//if (doSos) replaceSOSvars(cntVars);
		//replaceSOSvars(LHSs,RHSs,cntVars);
		if (info_level >= 2) cerr << "DELETE ALL CONSTRAINTS:" << QlpStSolve.getExternSolver(maxLPstage).getRowCount() << endl;
        if (QlpStSolve.getExternSolver(maxLPstage).getRowCount() > 0) {
		   QlpStSolve.removeUserCutsFromCut(maxLPstage);
		   QlpStSolve.getExternSolver(maxLPstage).clearLP_snapshot();
        }
		if (info_level >= 2) cerr << "DELETED ALL CONSTRAINTS:" << QlpStSolve.getExternSolver(maxLPstage).getRowCount() << endl;
		int targetConstraintsSize = 0, j=0;
		//std::vector< std::pair<int, std::pair<double,double>  > > replacer(qmip.getVariableVectorConst().size());
		//for (int zz = 0; zz < qmip.getVariableVectorConst().size();zz++) replacer[zz].first = -1;

		//totalRecallQlp = qlp;
		//totalRecallQlp.deleteAllRows();

		for (int i = 0; i < LHSs.size();i++)
			if (LHSs[i].size() > 0) { //durch sos koennen welche verkuerzt worden sein
			  QlpStSolve.getExternSolver(maxLPstage).addLProw_snapshot(LHSs[i], RHSs[i]);
			}

		const std::vector<data::QpNum>& tmpObjVec = qmip.getObjectiveFunctionValues();
		//data::Constraint& c = qlp.createRhsConstraint(rhs);
		std::vector<data::IndexedElement> obj_lhs;
		for (unsigned int i = 0; i < tmpObjVec.size(); i++) {
			if (!tmpObjVec[i].isZero())
				//c.createConstraintElement(i, tmpObjVec[i]);
				obj_lhs.push_back(data::IndexedElement(i, tmpObjVec[i]));
		}
		if (obj_lhs.size() > 0) {
			data::QpRhs obj_offset;
			obj_offset.setValue(0.0);
#ifdef FIND_BUG
			replaceSOSvarsInObj( obj_lhs, obj_offset , cntVars);
#endif
			for (int zz=0; zz < obj_lhs.size();zz++) {
				if (binary_assignments[obj_lhs[zz].index] == extbool_Undef) {
				    if (obj_lhs[zz].value.asDouble() >= 0) mon_neg[obj_lhs[zz].index]++;
				    else if (obj_lhs[zz].value.asDouble() < 0) mon_pos[obj_lhs[zz].index]++;
				}
			}
			data::QpRhs obj_rhs(data::QpRhs::smallerThanOrEqual,rhsbnd);
			QlpStSolve.addUserCut(maxLPstage, obj_lhs, data::QpRhs::smallerThanOrEqual, rhsbnd);
			QlpStSolve.getExternSolver(maxLPstage).addLPobj_snapshot(obj_lhs, obj_rhs);
			QlpStSolve.setObjIndex(targetConstraintsSize);
			for (int i = 0; i < cntVars;i++) {
				data::QpRhs zero;
				zero.setValue(0.0);
				QlpStSolve.changeObjFuncCoeff(maxLPstage, i, 0.0);
			}
			if(SHOW_DETAILS > 0) cerr << "objective: ";
			for (int i = 0; i < obj_lhs.size();i++) {
				if( cntVars ==QlpStSolve.getExternSolver(maxLPstage).getVariableCount() )				
					QlpStSolve.changeObjFuncCoeff(maxLPstage, obj_lhs[i].index, obj_lhs[i].value);
				if(SHOW_DETAILS > 0) cerr << obj_lhs[i].value << "x" << obj_lhs[i].index << " + ";
			}
			LPoffset = -obj_offset.getValue().asDouble();
			if(SHOW_DETAILS > 0) cerr << " + 0 <= " << LPoffset<< endl;
		}

        for (int i = 0; i < QlpStSolve.getExternSolver(maxLPstage).getLProws_snapshot();i++) {
	    if ((*QlpStSolve.getExternSolver(maxLPstage).getRowLhs_snapshot(i)).size() == 0) continue;
            if(1)QlpStSolve.addUserCut(maxLPstage,
                                       (*QlpStSolve.getExternSolver(maxLPstage).getRowLhs_snapshot(i)),
                                       (*QlpStSolve.getExternSolver(maxLPstage).getRowRhs_snapshot())[i].getRatioSign(),
                                       (*QlpStSolve.getExternSolver(maxLPstage).getRowRhs_snapshot())[i].getValue());
            QlpStSolve.getExternSolver(maxLPstage).setLazyStatus(i,false);
        }
        
        QlpStSolve.getExternSolver(maxLPstage).reinitLPcols_snapshot();
		sortCols(maxLPstage);

		if (info_level >= 2) cerr << "NEW CONSTRAINTS in extern Solver:" << QlpStSolve.getExternSolver(maxLPstage).getRowCount() << endl;
		if (info_level >= 2) cerr << "Size of relaxation:" << QlpStSolve.getExternSolver(maxLPstage).getLProws_snapshot() << "x" << cntVars << endl;
		if (info_level >= 2) cerr << "SOS CONSTRAINTS:" << cntSOS << " SOSvars:" << SOSvars.size() << endl;
		int qmipsz = qmip.getVariableVectorConst().size();
		//for (int zz = 0; zz < qmipsz;zz++) {
		//	cerr << "[" << zz << "," << SOSvars.count(zz) << "]";
		//}
		if (info_level >= 2) cerr << endl;
		int minBlock = qmipsz + 10;
		for (int zz = 0; zz < qmipsz;zz++)
			if (block[zz] < minBlock) minBlock = block[zz];
		if(0)for (int zz = 0; zz < qmipsz;zz++) {
			if (mon_neg[zz] > 0 && mon_pos[zz] == 0) {
				if (SOSvars.count(std::pair<int,int> (zz,0)) > 0) continue;
				//cerr << "MonNeg=" << mon_neg[zz] << " x" << zz << endl;
				if (types[zz] == 0 /*BINARY*/ && eas[zz] == EXIST && block[zz]==minBlock) {
					binary_assignments[zz] = 0;
					cpropQ.push_back(std::pair<int,double>(-zz-1, 0.0));
				}
			} else if (mon_pos[zz] > 0 && mon_neg[zz] == 0) {
				if (SOSvars.count(std::pair<int,int> (zz,0)) > 0) continue;
				//cerr << "MonPos=" << mon_pos[zz] << " x" << zz << endl;
				if (types[zz] == 0 /*BINARY*/ && eas[zz] == EXIST && block[zz]==minBlock) {
					binary_assignments[zz] = 1;
					cpropQ.push_back(std::pair<int,double>(-zz-1, 1.0));

				}
			}
		}
		updateContBounds(qmip, QlpStSolve, binary_assignments.data(), types, maxLPstage, LHSs, RHSs , cntVars, lowerBounds, upperBounds, cpropQ);
#endif

		//data::Qlp *stsoQlp = QlpStSolve.getQlp();
		//qmip = *stsoQlp;
	//qmip = tmpQlp;
		if(0)for (int z = 0;z < qmipsz;z++) {
			if (types[z] == 0) {
				//if (binary_assignments[z] == 0) totalRecallQlp.setVariableValue(z,0.0);
				//else if (binary_assignments[z] == 1) totalRecallQlp.setVariableValue(z,1.0);
			}
		}
		//utils::ToolBox::writeToFile("./test.qlp", totalRecallQlp.toQlpFileString(false));
	    delete hct;
	}

	qbp->saveNumberOfLinesOfLP(QlpStSolve.getExternSolver(maxLPstage).getRowCount());

    //delete all constraints
    //add saved constraints
    //delete hct;
	if(0)for (int i = 0; i < cntVars;i++) {
		if (QlpStSolve.getExternSolver(maxLPstage).getObj_snapshot(i).value.asDouble() != 0.0)
			cerr << " " << QlpStSolve.getExternSolver(maxLPstage).getObj_snapshot(i).value.asDouble() << "x" << QlpStSolve.getExternSolver(maxLPstage).getObj_snapshot(i).index;
	}
	//cerr << endl;
	if(0)for (int j = 0; j < QlpStSolve.getExternSolver(maxLPstage).getLProws_snapshot();j++) {
		std::vector<data::IndexedElement> * lhs = QlpStSolve.getExternSolver(maxLPstage).getRowLhs_snapshot(j);
		for (int i = 0; i < lhs->size();i++) {
			cerr << (*lhs)[i].value.asDouble() << "x" << (*lhs)[i].index << " + ";
		}
        cerr << endl;
	}
	int sz = cntVars;
	if(0)for (int i = 0; i < sz;i++) {
    	std::vector<data::IndexedElement> * col = QlpStSolve.getExternSolver(maxLPstage).getCol_snapshot(sortcols[i]);
    	cerr << sortcols[i] << ": " << QlpStSolve.getExternSolver(maxLPstage).getObj_snapshot(sortcols[i]).value.asDouble() << " | ";
	    for (int jj = 0; jj < col->size();jj++) cerr << " " << (*col)[jj].value.asDouble() << "u" << (*col)[jj].index;
	    cerr << endl;
	}
    //findSymmetries(QlpStSolve, maxLPstage, true, clist, types, block, binary_assignments.data());
}

#ifdef FIND_BUG
bool yInterface::extractSOS(std::vector<data::IndexedElement> &org_lhs, data::QpRhs &org_rhs, int cntVars) {
	//std::cerr << "BEGIN SOS EXTRACTION" << std::endl;
	// erstmal die schon gefundenen SOS Variablen ersetzen

	if(0)for (int i = 0; i < org_lhs.size(); i++) {
		if (( SOSvars.find(std::pair<int,int>(org_lhs[i].index,0))) != SOSvars.end()) {
			//replaceX_iInC1throughC2(org_lhs[i].index,org_lhs,org_rhs,SOSconstraints[index]);
			//if (!isSOS(org_lhs,org_rhs) return false;
		}
	}

	//  suche die Variable x_i mit gr��tem Index, die noch nicht in SOSvars drin ist
	int x_i = cntVars+2;
	int x_i_index = cntVars+2;
	if (org_lhs.size() == 0)  return false;
	if (org_rhs.getValue().asDouble() <= qbp->getDontKnowValue()) return false;
	assert(org_lhs.size() > 0);
	//std::cerr << " new line ";
	for (int i = 0; i < org_lhs.size(); i++) {
		//std::cerr << x_i << "<?" << org_lhs[i].index << "," << (int)(isOne_(org_lhs[i].value.asDouble()))<< "," << (int)((signed int)org_lhs[i].index > x_i)<< "  ";
		if ((x_i == cntVars+2 || (signed int)org_lhs[i].index > x_i) &&
				isOne_(org_lhs[i].value.asDouble()) &&
				SOSvars.count(std::pair<int,int> (org_lhs[i].index,0)) <= 0 &&
				SOStabus.count(org_lhs[i].index) <= 0) {
			//std::cerr << " yep " << org_lhs[i].index;
			x_i = org_lhs[i].index;
			x_i_index = i;
		}
	}
	//std::cerr << std::endl;
	//  tu x_i in SOSvars
	if (x_i == cntVars + 2) return false;
	assert(x_i < cntVars+2);
	bool isAR=true;
	if (SOSvars.count(std::pair<int,int>(x_i,0)) > 0) {
		//std::cerr << "x" << x_i << " already in SOSvars." << std::endl;
		return false;
	}
	SOSvars.insert(std::pair<int,int>(x_i,SOSconstraints.size()));
	//  f�ge x_i und seinen Term in SOSconstraints hinzu
	std::pair< std::vector< std::pair<int,double> >,int > SOSc;
	SOSc.second = x_i;
	std::vector< std::pair<int,double> > SOSc_lhs(org_lhs.size());
	for (int i = 0; i < org_lhs.size();i++) {
		if(org_lhs[i].value.asDouble() != org_lhs[i].value.asDouble()) isAR=false;
		std::pair<int, double> SOSc_lhs_elem;
		if (i != x_i_index) {
			SOSc_lhs_elem.first = org_lhs[i].index;
			SOSc_lhs_elem.second = org_lhs[i].value.asDouble();
			SOSc_lhs[i] = SOSc_lhs_elem;
		} else {
			SOSc_lhs_elem.first = cntVars+2;
			SOSc_lhs_elem.second = -org_rhs.getValue().asDouble();
			SOSc_lhs[i] = SOSc_lhs_elem;
		}
	}
	if (!isAR) SHOW_DETAILS=2;
	SOSc.first = SOSc_lhs;
	SOSconstraints.push_back(SOSc);
	for(int k = 0; k < SOSc_lhs.size(); k++) {
		if (org_lhs[k].index <= cntVars) SOStabus.insert(org_lhs[k].index);
	}
	if (SHOW_DETAILS > 1)
		std::cerr << "extrahierte SOS: (" << org_lhs.size()<< ","<< SOSc.first.size() << ")" << std::endl;
	//org_lhs[x_i_index] = org_lhs[org_lhs.size()-1];
	//org_lhs.pop_back();
	org_rhs.set(data::QpRhs::equal/*data::QpRhs::smallerThanOrEqual*/,org_rhs.getValue().asDouble());
	if (SHOW_DETAILS > 1) {
		for (int k=0;k<SOSconstraints.size();k++) {
			for (int kk=0;kk<SOSconstraints[k].first.size();kk++) {
				std::cerr << SOSconstraints[k].first[kk].second << "x" << SOSconstraints[k].first[kk].first << " + ";
			}
			if (SHOW_DETAILS > 1)
				std::cerr << std::endl << "extrahierte Variable:" << SOSconstraints[k].second << std::endl;
		}
		std::cerr << std::endl << "Extraktion beendet" << std::endl;
	}
	assert(isAR);
	//  Ungleichung zur�ckgeben
	//std::cerr << "FINISH SOS EXTRACTION" << std::endl;
	return true;
}

void yInterface::addNegativeOfSosC2ToSosC1(std::vector< std::pair<int,double> > &c1, std::vector< std::pair<int,double> > &c2, int x_i) {
	//std::cerr << "pivot "<< x_i <<  std::endl;
    // is x_i in c1? Wenn ja, hinten dran haengen
	bool found = false;
	bool isAR=true;
	int oldC1Size = c1.size();
	for (int ii = 0;ii < c1.size();ii++)
		if (c1[ii].second != c1[ii].second) isAR = false;
	assert(isAR);
	for (int ii = 0;ii < c2.size();ii++)
		if (c2[ii].second != c2[ii].second) isAR = false;
	assert(isAR);

	for (int i = 0; i < oldC1Size; i++) {
		if (x_i == c1[i].first) {
			//std::cerr << "pivot found "<<  std::endl;
			if (isZero_(c1[i].second)) {
				;
			} else if (c1[i].second < 0) {
				for (int j = 0; j < c2.size();j++) {
					assert(fabs(-c1[i].second*c2[j].second)<100000.0);
					c1.push_back(c2[j]);
					assert(c1[c1.size()-1].second == c2[j].second);
					c1[c1.size()-1].second = -c1[i].second*c1[c1.size()-1].second;
				}
				c1[i].second = 0.0;
				std::sort(c1.begin(), c1.end(), [](std::pair<int,double> p1, std::pair<int,double> p2){return p1.first < p2.first;});
			} else {
				for (int j = 0; j < c2.size();j++) {
					assert(fabs(-c1[i].second*c2[j].second)<100000.0);
					c1.push_back(c2[j]);
					c1[c1.size()-1].second = -c1[i].second*c1[c1.size()-1].second;
				}
				c1[i].second = 0.0;
				std::sort(c1.begin(), c1.end(), [](std::pair<int,double> p1, std::pair<int,double> p2){return p1.first < p2.first;});
			}
			found  = true;
			break;
		}
	}
	for (int ii = 0;ii < c1.size();ii++)
		if (c1[ii].second != c1[ii].second) isAR = false;
	assert(isAR);
	for (int ii = 0;ii < c2.size();ii++)
		if (c2[ii].second != c2[ii].second) isAR = false;
	assert(isAR);
	if (found) {
		if (SHOW_DETAILS > 1) {
			std::cerr << "raw sum"<< std::endl;
			for (int j = 0; j < c1.size();j++) {
				std::cerr << c1[j].second << "x" << c1[j].first << " + ";
			}
			std::cerr << std::endl;
		}
		// durchlaufen und addieren
		int i=0, j=0;
		for (i = 0; i < c1.size();i++) {
			j = i+1;
			while (j < c1.size() && c1[i].first == c1[j].first) {
				c1[i].second += c1[j].second;
				c1[j].second = 0.0;
				j++;
			}
			i = j-1;
		}
		// 0-en wegl�schen
		i=0, j=0;
		for (i = 0; i < c1.size();i++) {
			if (!isZero_(c1[i].second)) continue;
			j = c1.size()-1;
			while (j>=i && isZero_(c1[j].second)) {
				j--;
				c1.pop_back();
			}
			if (i >= j) break;
			else {
				c1[i] = c1[j];
				c1[j].second = 0.0;
				c1.pop_back();
			}
		}
	}
    if (c1.size()==0) return;
	//---------------------------
	//simplifyCoefficients_01(lhs,org_lhs, rhs, org_rhs, ba);

	int pow_max = 1;
	int pow_min = 1;
    for (int i = 0; i < c1.size();i++) {
		int64_t pow = 1;
		for (int j=0; j < 32;j++) {
			double coe = floor(fabs(c1[i].second) * ( (double)pow ));
			double frac = fabs(c1[i].second) * ( (double)pow ) - coe;
			if (frac > 1e-20) {
				pow = pow * 10;
			} else {
				if (pow > pow_max) pow_max = pow;
				else if (pow < pow_min) pow_min = pow;
				break;
			}
			if (j > 25) {
				cerr << "could not find correct power. Abort => Warning:" << c1[i].second << endl;
				break;
				//return;
			}
		}
	}

    assert(pow_max < 2);

	std::vector<int64_t> ilhs(c1.size());

	for (int i = 0; i < c1.size(); i++) {
		if (c1[i].second >= 0)
  		    ilhs[i] = (int64_t)(fabs(c1[i].second) * ( (double)pow_max ) );
		else
  		    ilhs[i] = -(int64_t)(fabs(c1[i].second) * ( (double)pow_max ) );
	}

	int64_t gt = ggt(ilhs[0], ilhs[1]);
	for (int i = 2; i < ilhs.size();i++) gt = ggt(ilhs[i],gt);

	assert(gt<=1);
	if (gt <= 1) {
		//for (int i = 0; i < /*ilhs*/c1.size();i++) cerr << "gt" << c1[i].second<< ";"; //ilhs[i] = ilhs[i] / gt;
	} else {
		cerr << "GT" << gt << ";";
		for (int i = 0; i < /*ilhs*/c1.size();i++) c1[i].second = c1[i].second / (double)gt; //ilhs[i] = ilhs[i] / gt;
	}


	//---------------------------

	for (int ii = 0;ii < c1.size();ii++)
		if (c1[ii].second != c1[ii].second) isAR = false;
	if (isAR==false) {
		for (int kk=0;kk<c1.size();kk++) {
			std::cerr << c1[kk].second << "x" << c1[kk].first << " + ";
		}
		std::cerr << std::endl;
		for (int kk=0;kk<c2.size();kk++) {
			std::cerr << c2[kk].second << "x" << c2[kk].first << " + ";
		}
		std::cerr << std::endl;
	}
	assert(isAR);
	for (int ii = 0;ii < c2.size();ii++)
		if (c2[ii].second != c2[ii].second) isAR = false;
	assert(isAR);
}

bool yInterface::checkStabilityOfMult(std::vector< std::pair<int,double> > &c1, std::vector< std::pair<int,double> > &c2, int x_i) {
	for (int i = 0; i < c1.size(); i++) {
		if (x_i == c1[i].first) {
			double hiCo=0.0;
			if (isZero_(c1[i].second)) {
				;
			} else if (c1[i].second < 0) {
				for (int j = 0; j < c2.size();j++) { hiCo = fabs(-c1[i].second*c2[j].second); if (hiCo > 10000.0) return false;   }
			} else {
				for (int j = 0; j < c2.size();j++) { hiCo = fabs(-c1[i].second*c2[j].second); if (hiCo > 10000.0) return false;   }
			}
		}
	}
	return true;
}

void yInterface::replaceSOSvars( int cntVars ) {
	//doSos=false;
	bool isAR=true;
	for (int i = 0; i < SOSconstraints.size();i++) {
		for (int ii = 0;ii < SOSconstraints[i].first.size();ii++)
			if (SOSconstraints[i].first[ii].second != SOSconstraints[i].first[ii].second) isAR = false;
	}
	assert(isAR);
	if (SHOW_DETAILS > 1) {
		std::cerr << "BEGIN SOS STRUCTURING," << SOSconstraints.size() << std::endl;
	}
    // sort SOSconstraint nach Index der zu ersetzenden Variable, damit binary search m�glich wird
	std::sort( SOSconstraints.begin(), SOSconstraints.end(), []( std::pair< std::vector< std::pair<int,double> >,int > p1,
			std::pair< std::vector< std::pair<int,double> >,int > p2 ){
		return p1.second < p2.second;
	} );

	for (int i = 0; i < SOSconstraints.size();i++) {
		for (int ii = 0;ii < SOSconstraints[i].first.size();ii++)
			if (SOSconstraints[i].first[ii].second != SOSconstraints[i].first[ii].second) isAR = false;
	}
	assert(isAR);
    // for all SOS constraints
	int replaced=0;
	for (int i = 0; i < SOSconstraints.size();i++) {
		for (int ii = 0;ii < SOSconstraints[i].first.size();ii++)
			if (SOSconstraints[i].first[ii].second != SOSconstraints[i].first[ii].second) isAR = false;
        //    bestimme x_i und replace x_i in den anderen SOS constraints
		int x_i = SOSconstraints[i].second;
		if (x_i < 0) continue;
		//std::cerr << "S:" << x_i << std::endl;
		for (int j = /*i+1*/0; j < SOSconstraints.size();j++) {
			//std::cerr << "T:" << SOSconstraints[j].second << std::endl;
			if (i!=j) {
				for (int ii = 0;ii < SOSconstraints[i].first.size();ii++)
					if (SOSconstraints[i].first[ii].second != SOSconstraints[i].first[ii].second) isAR = false;
				assert(isAR);
				for (int ii = 0;ii < SOSconstraints[j].first.size();ii++)
					if (SOSconstraints[j].first[ii].second != SOSconstraints[j].first[ii].second) isAR = false;
				assert(isAR);
				//-----------------------

				//-----------------------
				bool stability = checkStabilityOfMult(SOSconstraints[j].first, SOSconstraints[i].first, x_i);
				if (stability) {
					replaced++;
					addNegativeOfSosC2ToSosC1(SOSconstraints[j].first, SOSconstraints[i].first, x_i);
				}
				else SOSconstraints[i].second = -fabs(x_i)-1;
				for (int ii = 0;ii < SOSconstraints[i].first.size();ii++)
					if (SOSconstraints[i].first[ii].second != SOSconstraints[i].first[ii].second) isAR = false;
				assert(isAR);
				for (int ii = 0;ii < SOSconstraints[j].first.size();ii++)
					if (SOSconstraints[j].first[ii].second != SOSconstraints[j].first[ii].second) isAR = false;
				assert(isAR);
			}
		}
	}
	cerr << "IN SOS STRUCTURING: made " << replaced << " replacement operations." << endl;
	replaced = 0;
	for (int i = 0; i < SOSconstraints.size();i++) {
		if (SOSconstraints[i].second >= 0) {
			if (i < 30) cerr << "to replace: x" << SOSconstraints[i].second << endl;
			replaced++;
		}
	}
	cerr << "IN SOS STRUCTURING: replaced " << replaced << " SOS constraints." << endl;

    //// --------------------
	if (isAR==false) SHOW_DETAILS = 2;
	if (SHOW_DETAILS > 0) {
		std::cerr << std::endl << "Final sos. size=" << SOSconstraints.size() << std::endl;
		for (int k=0;k<SOSconstraints.size();k++) {
			for (int kk=0;kk<SOSconstraints[k].first.size();kk++) {
				std::cerr << SOSconstraints[k].first[kk].second << "x" << SOSconstraints[k].first[kk].first << " + ";
			}
			std::cerr << std::endl;
		}
		std::cerr << "FINISH SOS STRUCTURING" << std::endl;
	}
	assert(isAR);
}

void yInterface::replaceSOSvars( std::vector< std::vector<data::IndexedElement> > &LHSs, std::vector<data::QpRhs> &RHSs , int cntVars) {
	if (SHOW_DETAILS > 1)
		std::cerr << "BEGIN SOS REPLACEMENT," << RHSs.size() << ","<< LHSs.size() << ","<< SOSconstraints.size() << std::endl;

	for (int z=0;z<SOSconstraints.size();z++)
		cerr << SOSconstraints[z].second << " ";
	cerr << endl;
    // for all constraints
	for (int i = 0; i < RHSs.size();i++) {
		//cerr << "|";
		//for all elements
		for (int j = LHSs[i].size()-1; j >= 0; j--) {
			//cerr << LHSs[i][j].value.asDouble() << "X" << LHSs[i][j].index << " ";
		    // binary search
			int v = LHSs[i][j].index;
			int jj=SOSconstraints.size()-1,ii=0;
			int match=-1;
			while (jj >= ii) {
				int m = (ii + jj)>>1;
				//if (v==24) cerr << ii <<"," << jj <<  "," << m << "," << SOSconstraints[m].second<< endl;
				int w = fabs(SOSconstraints[m].second);
				if (w==v && SOSconstraints[m].second >= 0) {
					match = m;
					break;
				} else {
					if (ii==jj) break;
					if (v>w) {
						ii=m+1;
					} else if (v<w) {
						jj=m-1;
					}
				}
			}
			if (match!=-1) {
				//assert(0);
				    //cerr << "match=" << match << ";"<< LHSs[i][j].value.asDouble() << "X" << LHSs[i][j].index << " ";
					if (isZero_(LHSs[i][j].value.asDouble())) {
						//cerr << "nothing" << endl;
						continue;
					} else if (LHSs[i][j].value.asDouble() < 0) {
						for (int jj = 0; jj < SOSconstraints[match].first.size();jj++) {
							data::IndexedElement IE(SOSconstraints[match].first[jj].first,-SOSconstraints[match].first[jj].second*LHSs[i][j].value.asDouble());  //(index,oefficient)
							LHSs[i].push_back(IE);
						}
						//std::sort(c1.begin(), c1.end(), [](std::pair<int,double> p1, std::pair<int,double> p2){return p1.first < p2.first;});
					} else {
						for (int jj = 0; jj < SOSconstraints[match].first.size();jj++) {
							data::IndexedElement IE(SOSconstraints[match].first[jj].first,-SOSconstraints[match].first[jj].second*LHSs[i][j].value.asDouble());
							LHSs[i].push_back(IE);
						}
						//std::sort(c1.begin(), c1.end(), [](std::pair<int,double> p1, std::pair<int,double> p2){return p1.first < p2.first;});
					}
					//cerr << "set it zero:" << i << "," << j << endl;
					LHSs[i][j].value = 0.0;
			} //else cerr << "." << endl;
		}

		std::sort(LHSs[i].begin(), LHSs[i].end(), [](data::IndexedElement e1, data::IndexedElement e2) { return e1.index < e2.index;} );
        if (SHOW_DETAILS > 2) {
			std::cerr << "after final sorting"<< std::endl;
			for (int j = 0; j < LHSs[i].size();j++) {
				std::cerr << LHSs[i][j].value.asDouble() << "x" << LHSs[i][j].index << " + ";
			}
			std::cerr << (RHSs[i].getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (RHSs[i].getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << RHSs[i].getValue() << std::endl;
        }
		// durchlaufen und addieren
		int ii=0, jj=0;
		for (ii = 0; ii < LHSs[i].size();ii++) {
			jj = ii+1;
			//cerr << "LHSs[ii]=" <<  LHSs[i][ii].index << "LHS[jj]=" << LHSs[i][jj].index << endl;
			while (jj < LHSs[i].size() && LHSs[i][ii].index == LHSs[i][jj].index) {
				LHSs[i][ii].value = LHSs[i][ii].value.asDouble() + LHSs[i][jj].value.asDouble();
				LHSs[i][jj].value = 0.0;
				jj++;
			}
			ii = jj-1; // wird sofort wieder erhoeht
		}
		for (ii = 0; ii < LHSs[i].size();ii++) {
			if (ii < LHSs[i].size()-1) assert(LHSs[i][ii].index != LHSs[i][ii+1].index || LHSs[i][ii+1].value.asDouble() == 0.0);
		}
		if (SHOW_DETAILS > 2) {
			std::cerr << "after addition"<< std::endl;
			for (int j = 0; j < LHSs[i].size();j++) {
				std::cerr << LHSs[i][j].value.asDouble() << "x" << LHSs[i][j].index << " + ";
			}
			std::cerr << (RHSs[i].getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (RHSs[i].getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << RHSs[i].getValue() << std::endl;
        }
		// 0-en wegl�schen
		ii=0, jj=0;
		for (ii = 0; ii < LHSs[i].size();ii++) {
			if (!isZero_(LHSs[i][ii].value.asDouble())) continue;
			jj = LHSs[i].size()-1;
			while (jj>=ii && isZero_(LHSs[i][jj].value.asDouble())) {
				jj--;
				LHSs[i].pop_back();
			}
			if (ii >= jj) break;
			else {
				LHSs[i][ii] = LHSs[i][jj];
				LHSs[i][jj].value = 0.0;
				LHSs[i].pop_back();
			}
		}
        if (SHOW_DETAILS > 2) {
			std::cerr << "after deleting zeros"<< std::endl;
			for (int j = 0; j < LHSs[i].size();j++) {
				std::cerr << LHSs[i][j].value.asDouble() << "x" << LHSs[i][j].index << " + ";
			}
			std::cerr << "<=>" << RHSs[i].getValue().asDouble() << std::endl;
        }
		// dann die "rechten-Seiten" (es sollte <= 1 davon geben) nach rechts bringen
		for (ii = 0; ii < LHSs[i].size()-1;ii++) {
			if (LHSs[i][ii].index > cntVars) {
				data::IndexedElement tmp = LHSs[i][ii];
				LHSs[i][ii] = LHSs[i][LHSs[i].size()-1];
				LHSs[i][LHSs[i].size()-1] = tmp;
				break;
			}
		}
        if (SHOW_DETAILS > 2) {
			std::cerr << "after reordering"<< std::endl;
			for (int j = 0; j < LHSs[i].size();j++) {
				std::cerr << LHSs[i][j].value.asDouble() << "x" << LHSs[i][j].index << " + ";
			}
			std::cerr << (RHSs[i].getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (RHSs[i].getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << RHSs[i].getValue() << std::endl;
        }
		// zum schluss noch rechte seite korrigieren
		for (ii = LHSs[i].size()-1;ii>=0;ii--) {
			if (isZero_(LHSs[i][ii].value.asDouble())) {
				LHSs[i].pop_back();
				continue;
			}
			if (LHSs[i][ii].index < cntVars + 2) break;
			//RHSs[i].set(data::QpRhs::RatioSign::smallerThanOrEqual,RHSs[i].getValue() - LHSs[i][ii].value);
			RHSs[i].setValue(RHSs[i].getValue() - LHSs[i][ii].value);
			LHSs[i].pop_back();
		}
	}
	// end: put (relaxed) SOS constraints into LHSs and RHSs
	for (int i = 0; i < SOSconstraints.size();i++) {
		std::pair< std::vector< std::pair<int,double> >,int > &SOSc = SOSconstraints[i];
		int x_i = SOSc.second;
		std::vector<data::IndexedElement> sos_lhs;
		data::QpRhs sos_rhs;
		if (x_i < 0) { // could not be replaced -> take original equation sum x_i = 1;
			x_i = -x_i - 1;
			if (x_i > cntVars+2) {
				cerr << cntVars << " " << x_i << endl;
				assert(0);
			}
			for (int ii = 0; ii < SOSc.first.size();ii++) {
				if (ii < SOSc.first.size()-1) assert(SOSc.first[ii].first != SOSc.first[ii+1].first);
				if (SOSc.first[ii].first >= cntVars+2) { // rhs
					sos_rhs.set(data::QpRhs::RatioSign::equal, -SOSc.first[ii].second);
				} else {
					data::IndexedElement sos_lhs_elem;
					sos_lhs_elem.value = SOSc.first[ii].second;
					sos_lhs_elem.index = SOSc.first[ii].first;
					sos_lhs.push_back(sos_lhs_elem);
				}
			}
			data::IndexedElement sos_lhs_elem;
			sos_lhs_elem.value = 1.0;
			sos_lhs_elem.index = x_i;
			sos_lhs.push_back(sos_lhs_elem);
			LHSs.push_back(sos_lhs);
			RHSs.push_back(sos_rhs);
			/*
			std::cerr << std::endl << "RECONSTRUCTION:" << std::endl;
			for (int k=0;k<sos_lhs.size();k++) {
					std::cerr << sos_lhs[k].value.asDouble() << "x" << sos_lhs[k].index << " + ";
			}
			std::cerr << (sos_rhs.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (sos_rhs.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << sos_rhs.getValue() << std::endl;
			std::cerr << std::endl << "RECONSTRUCTION beendet" << std::endl;
			*/
		} else {
			//assert(0);
			for (int ii = 0; ii < SOSc.first.size();ii++) {
				if (SOSc.first[ii].first >= cntVars+1) { // rhs
					sos_rhs.set(data::QpRhs::RatioSign::smallerThanOrEqual, -SOSc.first[ii].second);
				} else {
					data::IndexedElement sos_lhs_elem;
					sos_lhs_elem.value = SOSc.first[ii].second;
					sos_lhs_elem.index = SOSc.first[ii].first;
					sos_lhs.push_back(sos_lhs_elem);
				}
			}
			LHSs.push_back(sos_lhs);
			RHSs.push_back(sos_rhs);
		}
	}
	if (SHOW_DETAILS > 1) {
		std::cerr << "FINAL constraints" << std::endl;
		int i = 0;
		for (; i < RHSs.size();i++) {
			//for all elements
			for (int j = 0; j < LHSs[i].size();j++) {
				std::cerr << LHSs[i][j].value.asDouble() << "x" << LHSs[i][j].index << " + ";
			}
			std::cerr << (RHSs[i].getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (RHSs[i].getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << RHSs[i].getValue() << std::endl;
		}
		std::cerr << "after " << i << " FINISH SOS REPLACEMENT" << std::endl;
		//char a; cin >> a;
	}
}

void yInterface::replaceSOSvarsInObj( std::vector<data::IndexedElement> &OBJ, data::QpRhs &OBJBND , int cntVars) {

	for (int z=0;z<SOSconstraints.size();z++)
		cerr << SOSconstraints[z].second << " ";
	cerr << endl;
    // for objective
	{
		//cerr << "|";
		//for all elements
		for (int j = OBJ.size()-1; j >= 0; j--) {
			//cerr << LHSs[i][j].value.asDouble() << "X" << LHSs[i][j].index << " ";
		    // binary search
			int v = OBJ[j].index;
			int jj=SOSconstraints.size()-1,ii=0;
			int match=-1;
			while (jj >= ii) {
				int m = (ii + jj)>>1;
				//if (v==24) cerr << ii <<"," << jj <<  "," << m << "," << SOSconstraints[m].second<< endl;
				int w = fabs(SOSconstraints[m].second);
				if (w==v && SOSconstraints[m].second >= 0) {
					match = m;
					break;
				} else {
					if (ii==jj) break;
					if (v>w) {
						ii=m+1;
					} else if (v<w) {
						jj=m-1;
					}
				}
			}
			if (match!=-1) {
				//assert(0);
				    //cerr << "match=" << match << ";"<< LHSs[i][j].value.asDouble() << "X" << LHSs[i][j].index << " ";
					if (isZero_(OBJ[j].value.asDouble())) {
						//cerr << "nothing" << endl;
						continue;
					} else if (OBJ[j].value.asDouble() < 0) {
						for (int jj = 0; jj < SOSconstraints[match].first.size();jj++) {
							data::IndexedElement IE(SOSconstraints[match].first[jj].first,-SOSconstraints[match].first[jj].second*OBJ[j].value.asDouble());  //(index,oefficient)
							OBJ.push_back(IE);
						}
						//std::sort(c1.begin(), c1.end(), [](std::pair<int,double> p1, std::pair<int,double> p2){return p1.first < p2.first;});
					} else {
						for (int jj = 0; jj < SOSconstraints[match].first.size();jj++) {
							data::IndexedElement IE(SOSconstraints[match].first[jj].first,-SOSconstraints[match].first[jj].second*OBJ[j].value.asDouble());
							OBJ.push_back(IE);
						}
						//std::sort(c1.begin(), c1.end(), [](std::pair<int,double> p1, std::pair<int,double> p2){return p1.first < p2.first;});
					}
					//cerr << "set it zero:" << i << "," << j << endl;
					OBJ[j].value = 0.0;
			} //else cerr << "." << endl;
		}

		std::sort(OBJ.begin(), OBJ.end(), [](data::IndexedElement e1, data::IndexedElement e2) { return e1.index < e2.index;} );
        if (SHOW_DETAILS > 2) {
			std::cerr << "after final sorting"<< std::endl;
			for (int j = 0; j < OBJ.size();j++) {
				std::cerr << OBJ[j].value.asDouble() << "x" << OBJ[j].index << " + ";
			}
			std::cerr << (OBJBND.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (OBJBND.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << OBJBND.getValue() << std::endl;
        }
		// durchlaufen und addieren
		int ii=0, jj=0;
		for (ii = 0; ii < OBJ.size();ii++) {
			jj = ii+1;
			//cerr << "OBJ[ii]=" <<  OBJ[ii].index << "OBJ[jj]=" << OBJ[jj].index << endl;
			while (jj < OBJ.size() && OBJ[ii].index == OBJ[jj].index) {
				OBJ[ii].value = OBJ[ii].value.asDouble() + OBJ[jj].value.asDouble();
				OBJ[jj].value = 0.0;
				jj++;
			}
			ii = jj-1; // wird sofort wieder erhoeht
		}
		for (ii = 0; ii < OBJ.size();ii++) {
			if (ii < OBJ.size()-1) assert(OBJ[ii].index != OBJ[ii+1].index || OBJ[ii+1].value.asDouble() == 0.0);
		}
		if (SHOW_DETAILS > 2) {
			std::cerr << "after addition"<< std::endl;
			for (int j = 0; j < OBJ.size();j++) {
				std::cerr << OBJ[j].value.asDouble() << "x" << OBJ[j].index << " + ";
			}
			std::cerr << (OBJBND.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (OBJBND.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << OBJBND.getValue() << std::endl;
        }
		// 0-en wegl�schen
		ii=0, jj=0;
		for (ii = 0; ii < OBJ.size();ii++) {
			if (!isZero_(OBJ[ii].value.asDouble())) continue;
			jj = OBJ.size()-1;
			while (jj>=ii && isZero_(OBJ[jj].value.asDouble())) {
				jj--;
				OBJ.pop_back();
			}
			if (ii >= jj) break;
			else {
				OBJ[ii] = OBJ[jj];
				OBJ[jj].value = 0.0;
				OBJ.pop_back();
			}
		}
        if (SHOW_DETAILS > 2) {
			std::cerr << "after deleting zeros"<< std::endl;
			for (int j = 0; j < OBJ.size();j++) {
				std::cerr << OBJ[j].value.asDouble() << "x" << OBJ[j].index << " + ";
			}
			std::cerr << "<=>" << OBJBND.getValue().asDouble() << std::endl;
        }
		// dann die "rechten-Seiten" (es sollte <= 1 davon geben) nach rechts bringen
		for (ii = 0; ii < OBJ.size()-1;ii++) {
			if (OBJ[ii].index > cntVars) {
				data::IndexedElement tmp = OBJ[ii];
				OBJ[ii] = OBJ[OBJ.size()-1];
				OBJ[OBJ.size()-1] = tmp;
				break;
			}
		}
        if (SHOW_DETAILS > 2) {
			std::cerr << "after reordering"<< std::endl;
			for (int j = 0; j < OBJ.size();j++) {
				std::cerr << OBJ[j].value.asDouble() << "x" << OBJ[j].index << " + ";
			}
			std::cerr << (OBJBND.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (OBJBND.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << OBJBND.getValue() << std::endl;
        }
		// zum schluss noch rechte seite korrigieren
		for (ii = OBJ.size()-1;ii>=0;ii--) {
			if (isZero_(OBJ[ii].value.asDouble())) {
				OBJ.pop_back();
				continue;
			}
			if (OBJ[ii].index < cntVars + 2) break;
			//OBJBND.set(data::QpRhs::RatioSign::smallerThanOrEqual,OBJBND.getValue() - OBJ[ii].value);
			OBJBND.setValue(OBJBND.getValue() - OBJ[ii].value);
			OBJ.pop_back();
		}
	}
	if (SHOW_DETAILS > 1) {
		std::cerr << "FINAL objective" << std::endl;
		int i = 0;
		{
			//for all elements
			for (int j = 0; j < OBJ.size();j++) {
				std::cerr << OBJ[j].value.asDouble() << "x" << OBJ[j].index << " + ";
			}
			std::cerr << (OBJBND.getRatioSign() == data::QpRhs::RatioSign::equal ? " = " :
						  (OBJBND.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual ? " <= " : " >= ") ) << OBJBND.getValue() << std::endl;
		}
		std::cerr << "after FINISH SOS REPLACEMENT and new BOUND=" << OBJBND.getValue() << std::endl;
		//char a; cin >> a;
	}
}
#endif // FIND_BUG

static int SHOW_UCB_DETAILS = 0;
void yInterface::updateContBounds(data::Qlp &qmip, utils::QlpStageSolver& QlpStSolve, int8_t *binary_assignments,
		int *types, int maxLPstage, std::vector< std::vector<data::IndexedElement> > &LHSs, std::vector<data::QpRhs> &RHSs , int cntVars,
		coef_t *lowerBounds, coef_t *upperBounds, std::vector<std::pair<int,double> > &cpropQ) {
    // update bounds of Reals and Ints
    //return;
    std::vector<data::QpVar *> varVec = qmip.getVariableVector();
    //for (int z=0;z<varVec.size();z++) {
    //	cerr << varVec[z]->getName() << "(" << (int)varVec[z]->getNumberSystem() << "," << varVec[z]->getUpperBound() << ")";
    //}
    //cerr << endl;
    int numConstraints = qmip.getConstraintCount();
    std::vector<const data::QpRhs *> rhsVec = qmip.getRhsVecConst();
    std::vector<const data::Constraint *> conVec = qmip.getConstraintVecConst();
    
    // try to find sharper bounds for int and real
    for (int i = 0; i < numConstraints;i++) {
        data::QpRhs rhs = *rhsVec[i];
	if(rhs.getResponsibility()==data::Constraint::UNIVERSAL) continue;
        std::vector<data::IndexedElement> lhs = conVec[i]->getElements();
        int cntReals=0;
        int cntInts=0;
        int jOfInterest=-1;
        data::QpNum upper=0.0, lower=0.0;
        for (int j = 0; j < lhs.size();j++) {
            int index = lhs[j].index;
            double coeff = lhs[j].value.asDouble();
            
            if (varVec[index]->getNumberSystem()==data::QpVar::binaries) {
                if (binary_assignments[index] == 2) {
                   if (coeff > 0) lower = lower - coeff;     // "-" because lateron to rhs
                   else if (coeff < 0) upper = upper - coeff;
                } else if (binary_assignments[index] == 1) {
                    upper = upper - coeff;
                    lower = lower - coeff;
                }
            } else if (varVec[index]->getNumberSystem()==data::QpVar::real) {
                jOfInterest = j;
                cntReals++;
            } else {
                jOfInterest = j;
                cntInts++;
            }
        }
        if (cntInts + cntReals == 1 && binary_assignments[lhs[jOfInterest].index] == extbool_Undef) {
            lower = rhs.getValue() + lower;
            upper = rhs.getValue() + upper;
            int index = lhs[jOfInterest].index;
            double coeff = lhs[jOfInterest].value.asDouble();
            if (coeff >= 0) {
                lower = lower / coeff;
                upper = upper / coeff;
            } else {
                data::QpNum u = upper;
                upper = lower / coeff;
                lower = u / coeff;
            }
            lower = lower.asDouble() - 1e-9 -LP_EPS;//-1.0 - 1e-5;
            upper = upper.asDouble() + 1e-9 +LP_EPS;//+1.0 + 1e-5;
            if (/*varVec[index]->getNumberSystem()!=data::QpVar::binaries &&*/ varVec[index]->getNumberSystem()!=data::QpVar::real) {
	      lower = ceil(lower.asDouble()) -LP_EPS;// - 1e-9;
	      upper = floor(upper.asDouble()) + LP_EPS;//1e-9;
            }
            if (rhs.getRatioSign() == data::QpRhs::equal) {
                if (lower > varVec[index]->getLowerBound() ) {
		  if (info_level > 2) cerr << "eq increase lower bound of x" << index << "from " << varVec[index]->getLowerBound().asDouble() << " to " << lower.asDouble() << endl;
                    varVec[index]->setLowerBound(lower);
                }
                if (upper < varVec[index]->getUpperBound() ) {
		  if (info_level > 2) cerr << "eq decrease upper bound of x" << index << "from " << varVec[index]->getUpperBound().asDouble() << " to " << upper.asDouble() << endl;
                    varVec[index]->setUpperBound(upper);
                }
                if (lower > lowerBounds[index] ) {
		  if (info_level > 2) cerr << "eq2 increase lower bound of x" << index << "from " << lowerBounds[index] << " to " << lower.asDouble() << endl;
                    lowerBounds[index] = lower.asDouble();
                }
                if (upper < upperBounds[index] ) {
		  if (info_level > 2) cerr << "eq2 decrease upper bound of x" << index << "from " << upperBounds[index] << " to " << upper.asDouble() << endl;
                    upperBounds[index] = upper.asDouble();
                }
            } else if ( ((rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual && coeff >= 0) || (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual && coeff < 0)) ) {
                if (upper < varVec[index]->getUpperBound() ) {
		  if (info_level > 2) cerr << "decrease upper bound of x" << index << "from " << varVec[index]->getUpperBound().asDouble() << " to " << upper.asDouble() << endl;
                    varVec[index]->setUpperBound(upper);
                }
                if (upper < upperBounds[index] ) {
		  if (info_level > 2) cerr << "decrease2 upper bound of x" << index << "from " << upperBounds[index] << " to " << upper.asDouble() << endl;
                    upperBounds[index] = upper.asDouble();
                }
            } else if ( ((rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual && coeff >= 0) || (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual && coeff < 0)) ){
                if (lower > varVec[index]->getLowerBound() ) {
		  if (info_level > 2) cerr << "increase lower bound of x" << index << "from " << varVec[index]->getLowerBound().asDouble() << " to " << lower.asDouble() << endl;
                    varVec[index]->setLowerBound(lower);
                }
                if (lower > lowerBounds[index] ) {
		  if (info_level > 2) cerr << "increase2 lower bound of x" << index << "from " << lowerBounds[index] << " to " << lower.asDouble() << endl;
                    lowerBounds[index] = lower.asDouble();
                }
            }
        }
    }
    // end: try to find sharper bounds for int and real
    
    return;
}

#define OUT_COMPONENTS 0
void yInterface::findComponents(data::Qlp &qmip, int8_t *assigns, int *components, std::vector< std::vector<int> > &varsOfComponents) {
	std::vector<const data::QpVar *> varVec = qmip.getVariableVectorConst();
	std::vector<const data::QpRhs *> rhsVec = qmip.getRhsVecConst();
	int n = varVec.size();
	int present_component = 1;
	//std::vector<int> components(n);
	std::vector<int> stack;
	std::vector<data::IndexedElement> lhs;
	std::vector< std::vector<int> > conVec(n);
	std::vector<const data::Constraint *> rows = qmip.getConstraintVecConst();
	int rs = rows.size();
	std::vector<bool> rowInd(rs);
	for (int i = 0; i < rs;i++) {
		lhs.clear();
		lhs = rows[i]->getElements();
		for (int jj = 0; jj < lhs.size();jj++) {
			conVec[lhs[jj].index].push_back(i);
		}
		rowInd[i] = false;
	}
	for (int i = 0; i < n; i++) components[i] = 0;
	for (int i = 0; i < n; i++) {
		while ((components[i] > 0 || assigns[i] != 2) && i < n) i++;
		if (i >= n) break;
		stack.push_back(i);
		components[i] = present_component;
		if (OUT_COMPONENTS) cerr << " set var " << i << " to component " << present_component << ", stacksize = " << stack.size() << endl;
		while (stack.size() > 0) {
			int var = stack.back();
			stack.pop_back();
			//std::vector<data::Constraint*> conVec = qmip.getVariableConstraints(var);
			if (OUT_COMPONENTS) cerr << "Var " << var << " occurs in " << conVec[var].size() << " consttraints" << endl;
			for (int j = 0; j < conVec[var].size();j++) {
				lhs.clear();
				if (rowInd[conVec[var][j]] == true) continue;
				lhs = rows[conVec[var][j]]->getElements();//conVec[j]->getElements();
				rowInd[conVec[var][j]] = true;
				if (OUT_COMPONENTS) cerr << "check constraint: " << j << " : " << present_component << ":"<< lhs.size() << " : ";
				for (int jj = 0; jj < lhs.size();jj++) {
					if (assigns[lhs[jj].index] != 2) continue;
					if (OUT_COMPONENTS) cerr << lhs[jj].value.asDouble() << "*" << lhs[jj].index << " ";
					if (components[lhs[jj].index] == 0) {
						stack.push_back(lhs[jj].index);
						components[lhs[jj].index] = present_component;
					}
				}
				if (OUT_COMPONENTS) cerr << endl;
			}
		}
		present_component++;
	}
	for (int i = 0; i< varsOfComponents.size();i++) varsOfComponents[i].clear();
	varsOfComponents.resize(n+1);
	for (int i = 0; i < n;i++)
		if (components[i] > 0) varsOfComponents[components[i]].push_back(i);
	cerr << "The Instance has " << present_component-1 << " components." << endl;

	if(0)for (int i = 1; i< present_component;i++) { // for each component find the biconnected sub-components
		Component C(n,rs,varsOfComponents,conVec,rows, assigns);
		C.isBiConnected(i);
	}

}

#define OUT_STRONG_COMPONENTS 0
// Tarjan's algorithm, following
// https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm

void yInterface::strongConnect(CliqueManager &CM, std::vector<node> &V, int8_t *assigns, int *components, int *index, node v, int *c_identifier, int *types) {
    v.index = *index;
    v.lowlink = *index;
    *index = *index + 1;
    stack_S.push_back(v);
    v.onStack = true;
    int j;

    // Consider successors of v

    //cerr << "+" << v.i / 2 << "+";
        j = CM.FirstAdjacentInConflictGraph(v.i);
        if (j >= 0) {
                //cerr << "j=" << j << endl;
                int kk = CM.NextAdjacentInConflictGraph(j);
                while (kk >= 0) {
                          //cerr << "j=" << j << " kk=" << kk << " a(j)=" << CM.getAdjacent(j) << " a(kk)=" << CM.getAdjacent(kk) << endl;
                          node w = V[CM.getAdjacent(kk)^1];
                          assert(types[w.i/2] != /*CONTINUOUS*/5000);
                      if (w.index == -1) {
                        // Successor w has not yet been visited; recurse on it
                        strongConnect(CM, V, assigns, components, index, w, c_identifier, types);
                        v.lowlink = fmin(v.lowlink, w.lowlink);
                        //cerr << ".";
                      } else if (w.onStack) {
                        // Successor w is in stack S and hence in the current SCC
                        v.lowlink  = fmin(v.lowlink, w.index);
                        //cerr << ":";
                      } else {
                             //   cerr << "-";
                      }

                          kk = CM.NextAdjacentInConflictGraph(kk);
                }
        }
        // If v is a root node, pop the stack and generate an SCC
        if (v.lowlink == v.index) {
          node w;
          //start a new strongly connected component
          do {
            w = stack_S[stack_S.size()-1];
            stack_S.pop_back();
            w.onStack = false;
            //add w to current strongly connected component
            components[w.i] = *c_identifier;
          } while (w.i != v.i);
          //output the current strongly connected component
          //cerr << "component " << *c_identifier << ": ";
          //for (int z=0;z<V.size();z++)
            //  if (components[V[z].i]==*c_identifier) cerr << V[z].i << " ";
          //cerr << endl;
          (*c_identifier)++;
        }
    }


void yInterface::findStrongComponents(CliqueManager &CM, int8_t *assigns, int *components, int *c_identifier, int *types) {
    cerr << "N=" << getnVars() << endl;
    int cntVars2 = 0;
	for (int vi = 0; vi < 2*getnVars();vi++) {
    	int j = CM.FirstAdjacentInConflictGraph(vi);
        int cnt=0;
    	if (j >= 0) {
    		//cerr << "j=" << j << endl;
    		int kk = CM.NextAdjacentInConflictGraph(j);
    		while (kk >= 0) {
      			  //cerr << "j=" << j << " kk=" << kk << " a(j)=" << CM.getAdjacent(j) << " a(kk)=" << CM.getAdjacent(kk) << endl;
    			  cnt++;
    			  kk = CM.NextAdjacentInConflictGraph(kk);
    		}
    	}
    	if (cnt > 0) cntVars2++;
    	//cerr << vi << "-te Variable hat " << cnt << endl;
	}

    cerr << "there are " << cntVars2 << " Coevars with successors" << endl;

	std::vector<node> V(2*getnVars());
	*c_identifier=0;
	for (int i = 0; i < 2*getnVars(); i++) {
		V[i].index = -1;
		V[i].i = i;
		V[i].onStack = false;
		V[i].lowlink = -1;
		components[i] = -1;
	}
	int index = 0;
	stack_S.clear();
return;
	for (int i = 0; i < 2*getnVars();i++) {
		if (V[i].index == -1 && assigns[V[i].i >> 1] == 2) {
			//cerr << "starte mit coevar " << i << endl;
			strongConnect(CM, V, assigns, components, &index, V[i], c_identifier, types);
		}
	}

	std::vector<int> invComp(*c_identifier);
	for (int i=0; i < *c_identifier;i++) invComp[i]=0;
	for (int i=0; i < getnVars();i++) {
		invComp[components[2*i]]++;
		invComp[components[2*i+1]]++;
	}

	for (int i=0; i < *c_identifier;i++)
		if (invComp[i]>1) cerr << "Komponente " << i << " hat " << invComp[i] << " viele Teilnehmer." << endl;


	for (int i = 0; i < getnVars();i++) {
		//cerr << i << " " << (int)assigns[i] << " " << components[2*i] << " " << components[2*i+1] << endl;
		//assert(assigns[i] != 2 || components[2*i] > -1);
		//assert(assigns[i] != 2 || components[2*i+1] > -1);
		if (components[2*i] == components[2*i+1] && components[2*i+1] >= 0)
			assert(0);
	}

}

#define MAXABSINT 1e9 
#ifndef NEW_MKBIN
int yInterface::makeBinary(data::Qlp &qmip, data::Qlp &qbp) {

	// fill the qbp with the variables
	std::vector<data::QpVar *> varVec = qmip.getVariableVector();
	for (int i=0; i < varVec.size(); i++) {
		if (varVec[i]->getNumberSystem()!=data::QpVar::real && varVec[i]->getNumberSystem()!=data::QpVar::binaries) {
		  // is general
		  if (isZero_(varVec[i]->getLowerBound()) && isOne_(varVec[i]->getUpperBound())) {
		    varVec[i]->setNumberType(data::QpVar::binaries);
		  }
		}
	}

	std::vector<double> orgLBs;
	std::vector<double> orgUBs;
	for (int i = 0; i < varVec.size();i++) {
	  orgLBs.push_back(varVec[i]->getLowerBound().asDouble());
	  orgUBs.push_back(varVec[i]->getUpperBound().asDouble());
	}
	//for (int z=0;z<varVec.size();z++) {
	//	cerr << varVec[z]->getName() << "(" << (int)varVec[z]->getNumberSystem() << "," << varVec[z]->getUpperBound() << ")";
	//}
	//cerr << endl;
    int numConstraints = qmip.getConstraintCount();
    std::vector<const data::QpRhs *> rhsVec = qmip.getRhsVecConst();
    std::vector<const data::Constraint *> conVec = qmip.getConstraintVecConst();
    std::vector<bool> useless;
    
    // try to find sharper bounds for int and real
    for (int i = 0; i < numConstraints;i++) {
        data::QpRhs rhs = *rhsVec[i];
	data::QpRhs new_rhs = *rhsVec[i];
        std::vector<data::IndexedElement> lhs = conVec[i]->getElements();
        int cntReals=0;
        int cntInts=0;
        int jOfInterest=-1;
        data::QpNum upper=0.0, lower=0.0;
	double rhsVal = rhs.getValue().asDouble();
	useless.push_back(false);
	cerr.precision(17);
    LdejaVu:;
	if(1){
	  if (lhs.size() == 0) {
	     
	      if (new_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual && rhsVal <= 0.0){
	        if(new_rhs.getValue().asDouble()> 0.0) std::cerr << "ERROR: NOT FULFILLED CONSTRAINT 0 >= "<<new_rhs.getValue().asDouble() <<std::endl;
	        useless[useless.size()-1] = true;
	        }
	      if (new_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual && rhsVal >= 0.0){
	 	if(new_rhs.getValue().asDouble()< 0.0) std::cerr << "ERROR: NOT FULFILLED CONSTRAINT 0 <= " << new_rhs.getValue().asDouble() <<std::endl;
		useless[useless.size()-1] = true;
	      }
	      if (new_rhs.getRatioSign() == data::QpRhs::equal && fabs(rhsVal) <= 1e-9)
		useless[useless.size()-1] = true;
	  }
	  if (lhs.size() == 1 /*&& varVec[ lhs[0].index ]->getQuantifier() == data::QpVar::exists*/) {
	   if(info_level>=2){
	     std::cerr<< "Warning while reading: found constraint that simlified to " << lhs[0].value.asDouble()<< varVec[ lhs[0].index ]->getName();
	     if(new_rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) std::cerr<< ">=";
	     else if(new_rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual)std::cerr<< "<=";
	     else std::cerr<< "==";
	     std::cerr<<rhsVal<<std::endl;
  	    }
	    double rhsVal = new_rhs.getValue().asDouble();
	    double co0 = lhs[0].value.asDouble();
	    int old_index = lhs[0].index;
	    double locEps;
	    data::QpRhs::RatioSign rs = new_rhs.getRatioSign();
	    if(old_index == 1)  cerr << "V:" << varVec[ lhs[0].index ]->getName() << " " << co0 << "x" << old_index << (rs == data::QpRhs::smallerThanOrEqual ? " <= " : (rs == data::QpRhs::greaterThanOrEqual ? " >= " : " == " )) << rhsVal << endl;
            if (co0 < 0.0) {
	      //cerr << "VORHER:" << co0 << "x" << old_index << (rs == data::QpRhs::smallerThanOrEqual ? " <= " : (rs == data::QpRhs::greaterThanOrEqual ? " >= " : " == " )) << rhsVal << endl;
	      if (rs == data::QpRhs::greaterThanOrEqual) rs = data::QpRhs::smallerThanOrEqual; 
	      else if (rs == data::QpRhs::smallerThanOrEqual) rs = data::QpRhs::greaterThanOrEqual; 
	      rhsVal = -rhsVal;
	      co0 = -co0;
	      //cerr << "NACHHER:" << co0 << "x" << old_index << (rs == data::QpRhs::smallerThanOrEqual ? " <= " : (rs == data::QpRhs::greaterThanOrEqual ? " >= " : " == " )) << rhsVal << endl;
	    } 
	    if (varVec[old_index]->getNumberSystem()==data::QpVar::binaries) locEps = 0.0;
	    else locEps = 0.0;//1e-8;
	    if (fabs(co0) < 1e-10) {
	      if (rs == data::QpRhs::greaterThanOrEqual && rhsVal <= 0.0)
	        useless[useless.size()-1] = true;
	      if (rs == data::QpRhs::smallerThanOrEqual && rhsVal >= 0.0)
		useless[useless.size()-1] = true;
	      if (rs == data::QpRhs::equal && fabs(rhsVal) <= 1e-9)
		useless[useless.size()-1] = true;
	      else{
                if(info_level>=2){
  	           std::cerr<< "Warning while reading: Found potentially unsatisfiable constraint that simlified to " << co0<< varVec[ lhs[0].index ]->getName() << std::endl ;
	           if(rs == data::QpRhs::greaterThanOrEqual) std::cerr << " >= ";
	           else if(rs == data::QpRhs::smallerThanOrEqual) std::cerr << " <= ";
	           else std::cerr << " = ";
	           std::cerr <<  rhsVal << std::endl;
		}
	      } 
	    }

	    if (rs == data::QpRhs::greaterThanOrEqual || rs == data::QpRhs::equal) {
	      if (rhsVal/fabs(co0) > orgLBs[old_index] + locEps) {
	        if(rhsVal/fabs(co0) > orgUBs[old_index]+ locEps){
	          if(info_level>=2){
	            std::cerr << "Info while reading: simplified constraint of size 1 led to contradiction with variable bounds:" << std::endl;
	            std::cerr <<co0<< varVec[ lhs[0].index ]->getName()<<" >= " << rhsVal<< "; but "<< orgLBs[old_index] << " <= " << varVec[lhs[0].index ]->getName() << " <= " << orgUBs[old_index] <<std::endl;
		  }
		  if(new_rhs.getResponsibility()==data::QpRhs::UNIVERSAL) return 1;
		  else return 2;
	        }
	        else if( (new_rhs.getResponsibility()==data::QpRhs::EXISTENTIAL && varVec[ lhs[0].index ]->getQuantifier() == data::QpVar::exists)
	        	|| (new_rhs.getResponsibility()==data::QpRhs::UNIVERSAL && varVec[ lhs[0].index ]->getQuantifier() == data::QpVar::all))
	          //Only a universal constraint can force the fixation of a universally quantified constraint; the same holds for existentially quantified variables
	        {
		  //std::cerr << "Change x" << old_index << " lb From " << varVec[old_index]->getLowerBound().asDouble() << " to " << rhsVal/fabs(co0)-1e-10 << " ub=" << varVec[old_index]->getUpperBound().asDouble() << " " << (varVec[ lhs[0].index ]->getQuantifier() == data::QpVar::exists) <<std::endl;
		  orgLBs[old_index] = rhsVal/fabs(co0) - locEps; 
	  	  //varVec[old_index]->setLowerBound(rhsVal/fabs(co0) - locEps);
		  if (0&&new_rhs.getRatioSign() == data::QpRhs::equal) {
		    if (rhsVal/fabs(co0) < varVec[old_index]->getUpperBound().asDouble() ) {
		      varVec[old_index]->setUpperBound(rhsVal/fabs(co0));
		    }   
		  }  
		}
	      }
	      //useless[useless.size()-1] = true;
	    }
	    if (rs == data::QpRhs::smallerThanOrEqual || rs == data::QpRhs::equal) {
	      if (rhsVal/fabs(co0) < orgUBs[old_index] - locEps) {
		if(rhsVal/fabs(co0) < orgLBs[old_index]+ locEps){
                  if(info_level>=2){
                    std::cerr << "Info: simplified constraint of size 1 led to contradiction with variable bounds:" << std::endl;
		    std::cerr <<co0<< varVec[ lhs[0].index ]->getName()<<" <= " << rhsVal<< "; but "<< orgLBs[old_index] << " <= " << varVec[lhs[0].index ]->getName() << " <= " << orgUBs[old_index] <<std::endl;
		  }
	 	  if(new_rhs.getResponsibility()==data::QpRhs::UNIVERSAL) return 1;
		  else return 2;
		}
		else if( (new_rhs.getResponsibility()==data::QpRhs::EXISTENTIAL && varVec[ lhs[0].index ]->getQuantifier() == data::QpVar::exists)
		                        || (new_rhs.getResponsibility()==data::QpRhs::UNIVERSAL && varVec[ lhs[0].index ]->getQuantifier() == data::QpVar::all))
		  //Only a universal constraint can force the fixation of a universally quantified constraint; the same holds for existentially quantified variables
		  {
 		  //std::cerr << "Change x" << old_index << " ub From " << varVec[old_index]->getUpperBound().asDouble() << " to " << rhsVal/fabs(co0)+locEps << " lb=" << varVec[old_index]->getLowerBound().asDouble() << " " << (varVec[ lhs[0].index ]->getQuantifier() == data::QpVar::exists) <<std::endl;
		  orgUBs[old_index] = rhsVal/fabs(co0) + locEps; 
		  //varVec[old_index]->setUpperBound(rhsVal/fabs(co0) + locEps );
		  if (0&&new_rhs.getRatioSign() == data::QpRhs::equal) {
		    if (rhsVal/fabs(co0) > varVec[old_index]->getLowerBound().asDouble() ) {
		      varVec[old_index]->setLowerBound(rhsVal/fabs(co0) );
		    }
		  } 
		}
	      } 
	      //useless[useless.size()-1] = true;
	    } 
	  }

	  if (useless[useless.size()-1]) continue;
	  bool hasBeenShortened = false;
	  for (int j = 0; j < lhs.size();j++) {
            int old_index = lhs[j].index;
            double coeff = lhs[j].value.asDouble();
            if (varVec[old_index]->getNumberSystem()==data::QpVar::binaries) {
	      if (fabs( orgUBs[old_index] - orgLBs[old_index] ) < 0.99 ) {
		if (isOne_(orgUBs[old_index])) {
		  new_rhs.setValue( new_rhs.getValue().asDouble() - coeff  );
		}
		if (isZero_(orgLBs[old_index]) || isOne_(orgUBs[old_index]) ) {
		  if (j < lhs.size() -1) lhs[j] = lhs[lhs.size()-1];
		  lhs.pop_back(); 
		  hasBeenShortened = true;
		}
	      }
	    } else {
	      if (fabs( orgUBs[old_index] - orgLBs[old_index] ) < 1e-9 ) {
		new_rhs.setValue( new_rhs.getValue().asDouble() - coeff*(orgUBs[old_index] + orgLBs[old_index]) * 0.5  );
		if (j < lhs.size() -1) lhs[j] = lhs[lhs.size()-1];
		lhs.pop_back(); 
		hasBeenShortened = true;
	      }
	    }
	  }
	  if (hasBeenShortened) goto LdejaVu;
	}

        for (int j = 0; j < lhs.size();j++) {
            int index = lhs[j].index;
            double coeff = lhs[j].value.asDouble();

            if (varVec[index]->getNumberSystem()==data::QpVar::binaries) {
                if (coeff > 0) lower = lower - coeff;     // "-" because lateron to rhs
                else if (coeff < 0) upper = upper - coeff;
            } else if (varVec[index]->getNumberSystem()==data::QpVar::real) {
                jOfInterest = j;
                cntReals++;
            } else if (varVec[index]->getQuantifier() == data::QpVar::exists){
                jOfInterest = j;
                cntInts++;
            }
        }
        if (cntInts + cntReals == 1) {
            lower = rhs.getValue() + lower;
            upper = rhs.getValue() + upper;
            int index = lhs[jOfInterest].index;
            double coeff = lhs[jOfInterest].value.asDouble();
            if (coeff >= 0) {
                lower = lower / coeff;
                upper = upper / coeff;
            } else {
                data::QpNum u = upper;
                upper = lower / coeff;
                lower = u / coeff;
            }
            lower = lower - LP_EPS;//0001;// - 1e-5;
            upper = upper + LP_EPS;//0001;// + 1e-5;
            if (/*varVec[index]->getNumberSystem()!=data::QpVar::binaries &&*/ varVec[index]->getNumberSystem()!=data::QpVar::real) {
                lower = ceil(lower.asDouble());
                upper = floor(upper.asDouble())+LP_EPS;
            }
            if (rhs.getRatioSign() == data::QpRhs::equal) {
	      if (lower.asDouble() > /*varVec[index]->getLowerBound()*/orgLBs[index] ) {
		  if (info_level > 2) cerr << "eq increase lower bound of x" << index << "from " << varVec[index]->getLowerBound().asDouble() << " to " << lower.asDouble() << endl;
		  //varVec[index]->setLowerBound(lower);
		  orgLBs[index] = lower.asDouble();
                }
                if (upper.asDouble() < /*varVec[index]->getUpperBound()*/orgUBs[index] ) {
		  if (info_level > 2) cerr << "eq decrease upper bound of x" << index << "from " << varVec[index]->getUpperBound().asDouble() << " to " << upper.asDouble() << endl;
		  //varVec[index]->setUpperBound(upper);
		  orgUBs[index] = upper.asDouble();
                }
            } else if (((rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual && coeff >= 0) || (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual && coeff < 0)) ) {
	      if (upper.asDouble() < /*varVec[index]->getUpperBound()*/orgUBs[index] ) {
		  if (info_level > 2) cerr << "decrease upper bound of x" << index << "from " << varVec[index]->getUpperBound().asDouble() << " to " << upper.asDouble() << endl;
		  //varVec[index]->setUpperBound(upper);
		  orgUBs[index] = upper.asDouble();
                }
            } else if ( ((rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual && coeff >= 0) || (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual && coeff < 0)) ){
	      if (lower.asDouble() > /*varVec[index]->getLowerBound()*/orgLBs[index] ) {
		  if (info_level > 2) cerr << "increase lower bound of x" << index << "from " << varVec[index]->getLowerBound().asDouble() << " to " << lower.asDouble() << endl;
		  //varVec[index]->setLowerBound(lower);
		  orgLBs[index] = lower.asDouble();
                }
            }
        }
    }
    // end: try to find sharper bounds for int and real
    
	int index = 0;
	for (int i=0; i < varVec.size(); i++) {
		int leader_index = index;
		integers_lim.push_back(index);
		assert(integers_lim[integers_lim.size()-1] == index);
		assert(integers.size() == index);
		//if (varVec[i]->getNumberSystem()==data::QpVar::real) cerr << "x" << i << " is real" << endl;
		//else if (varVec[i]->getNumberSystem()==data::QpVar::binaries) cerr << "x" << i << " is binary" << endl;
		//else cerr << "x" << i << " is general" << endl;
		//if (varVec[i]->getNumberSystem()!=data::QpVar::real && varVec[i]->getNumberSystem()!=data::QpVar::binaries) {
		//  // is general
		//  if (isZero_(varVec[i]->getLowerBound()) && isOne_(varVec[i]->getUpperBound())) {
		//    //cerr << "change NumberType" << endl;
		//    varVec[i]->setNumberType(data::QpVar::binaries);
		//  }
		//}
		//if (varVec[i]->getNumberSystem()==data::QpVar::real) cerr << "y" << i << " is real" << endl;
		//else if (varVec[i]->getNumberSystem()==data::QpVar::binaries) cerr << "y" << i << " is binary" << endl;
		//else cerr << "y" << i << " is general" << endl;
		if (varVec[i]->getNumberSystem()==data::QpVar::real || varVec[i]->getNumberSystem()==data::QpVar::binaries) {
			data::QpVar var = data::QpVar(varVec[i]->getName(), index,
					varVec[i]->getLowerBound(), varVec[i]->getUpperBound(),
					varVec[i]->getNumberSystem(), varVec[i]->getQuantifier());
			qbp.createVariable(var);
			IntInfo varInfo;
			varInfo.index = index;
			varInfo.pt2leader = leader_index;
			varInfo.org_lb = varVec[i]->getLowerBound();
			varInfo.org_ub = varVec[i]->getUpperBound();
			varInfo.org_ind = i;
			varInfo.bitcnt = 1;
			varInfo.number=-1;
			varInfo.name = var.getName();
			integers.push_back(varInfo);
			index++;
		} else {
		  int num_bits;
		  if ( ceil(varVec[i]->getUpperBound().asDouble())-floor(varVec[i]->getLowerBound().asDouble()) > 0 ) 
		    num_bits = (int)( abs(log2(ceil(varVec[i]->getUpperBound().asDouble())-floor(varVec[i]->getLowerBound().asDouble()))) ) + 1;
		  else num_bits = 1;
		  //int num_bits = (int)( abs(log2(round(varVec[i]->getUpperBound().asDouble())-round(varVec[i]->getLowerBound().asDouble()))) ) + 1;
            if (num_bits > 40) num_bits = 40;
        	std::vector<data::IndexedElement> artif_lhs;
        	for (int j = num_bits-1; j >= 0; j--) {
    			data::QpVar var = data::QpVar(varVec[i]->getName() + "_b" + std::to_string(j), index,
    					0.0, 1.0,
    					data::QpVar::binaries, varVec[i]->getQuantifier());
    			qbp.createVariable(var);
			//cerr << "var " << index << " bitix=" << j << " BITCNT=" << num_bits << " name:" << var.getName() << endl;
    			IntInfo varInfo;
    			varInfo.index = index;
    			varInfo.pt2leader = leader_index;
    			varInfo.org_lb = varVec[i]->getLowerBound();
    			varInfo.org_ub = varVec[i]->getUpperBound();
    			varInfo.org_ind = i;
    			varInfo.name = var.getName();
    			varInfo.bitcnt = num_bits;
    			varInfo.number=-1;
			integers.push_back(varInfo);
            	data::IndexedElement artif_lhs_elem;
            	artif_lhs_elem.index = index;
                artif_lhs_elem.value = ipow2(j);
    			artif_lhs.push_back(artif_lhs_elem);
    			index++;
    			//if (info_level > 0) cerr << artif_lhs_elem.value << "x" << artif_lhs_elem.index << " + ";
            }
            // add artificial constraints sum a*2^i*x <= old_ub - old_lb
        	//if (info_level > 0) cerr << " <= " << varVec[i]->getUpperBound()-varVec[i]->getLowerBound() << ": vorher: lb=" << varVec[i]->getLowerBound() << " ub=" << varVec[i]->getUpperBound() << endl;
            if(varVec[i]->getQuantifier() == data::QpVar::exists){
	  	data::Constraint& d = qbp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, min(varVec[i]->getUpperBound()-varVec[i]->getLowerBound(),MAXABSINT), data::QpRhs::EXISTENTIAL);
                d.setElements(artif_lhs);
 		data::Constraint& d2 = qbp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, min(varVec[i]->getUpperBound()-varVec[i]->getLowerBound(),MAXABSINT), data::QpRhs::UNIVERSAL);
                d2.setElements(artif_lhs);
	    }
	    else if( fabs(min(varVec[i]->getUpperBound().asDouble()-varVec[i]->getLowerBound().asDouble(),MAXABSINT)-(pow(2,num_bits)-1) ) > 1e-6 ){
                //cerr << "Info: Added universal upper Bound Constraint with rhs "<<  min(varVec[i]->getUpperBound()-varVec[i]->getLowerBound(),MAXABSINT) << endl;
		data::Constraint& d = qbp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, min(varVec[i]->getUpperBound()-varVec[i]->getLowerBound(),MAXABSINT), data::QpRhs::UNIVERSAL);
                d.setElements(artif_lhs);
		data::Constraint& d2 = qbp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, min(varVec[i]->getUpperBound()-varVec[i]->getLowerBound(),MAXABSINT), data::QpRhs::EXISTENTIAL);
                d2.setElements(artif_lhs);
            }
            // end: add artificial constraints

		}
	}
	// end: fill the qbp with the variables


  int ifl = getInfoLevel();
  if(ifl >=2) std::cerr << "InfoLevel is = " << getInfoLevel() << std::endl;
  if (ifl <0) ifl = -ifl;
  if (ifl % 2 == 1 && getInfoLevel() <= -20) {
    std::cerr << "no binarization, although reached makeBinary" << std::endl;
    qbp = qmip;
    return 0;
  } else if(ifl >=2)     std::cerr << "DO makeBinary" << std::endl;


	std::vector<const data::QpVar *> new_varVec = qbp.getVariableVectorConst();
	//for (int z=0;z<new_varVec.size();z++)
	//	cerr << new_varVec[z]->getName() << " ";
	//cerr << endl;

	// convert and add constraints
    for (int i = 0; i < numConstraints;i++) {
    	data::QpRhs new_rhs = *rhsVec[i];
    	std::vector<data::IndexedElement> old_lhs = conVec[i]->getElements();
    	std::vector<data::IndexedElement> new_lhs;
	if (useless[i]) continue;
        for (int j = 0; j < old_lhs.size();j++) {
        	int old_index = old_lhs[j].index;
        	int new_index = integers[integers_lim[old_index]].index;
        	data::IndexedElement new_lhs_elem;
        	new_lhs_elem.index = new_index;
        	new_lhs_elem.value = old_lhs[j].value;
		if (fabs(integers[integers_lim[old_index]].org_lb.asDouble() - integers[integers_lim[old_index]].org_ub.asDouble()) < 1e-9 ) {
		  //new_rhs.setValue(new_rhs.getValue().asDouble() - (integers[integers_lim[old_index]].org_ub.asDouble() + integers[integers_lim[old_index]].org_lb.asDouble()) * 0.5 );
		  new_rhs.setValue(new_rhs.getValue().asDouble() - old_lhs[j].value.asDouble() *(integers[integers_lim[old_index]].org_ub.asDouble() + integers[integers_lim[old_index]].org_lb.asDouble()) * 0.5 );

		  continue; 
		}
        	if (integers[integers_lim[old_index]].index == integers[integers_lim[old_index]].pt2leader &&
        		integers[integers_lim[old_index]].bitcnt == 1) {
		  if (varVec[ old_index ]->getNumberSystem()!=data::QpVar::real && varVec[ old_index ]->getNumberSystem()!=data::QpVar::binaries) {
		    //cerr << "Var x" << old_index << " is not binary and should not be real. lb=" << integers[integers_lim[old_index]].org_lb.asDouble() << " and  ub=" << integers[integers_lim[old_index]].org_ub.asDouble() << endl;
		    //if (varVec[ old_index ]->getNumberSystem()==data::QpVar::binaries) cerr << " is binary." << endl;
		    //if (varVec[ old_index ]->getNumberSystem()==data::QpVar::real) cerr << " is real" << endl;
		    //else cerr << " is Integer" << endl;
		    new_rhs.setValue(new_rhs.getValue() - old_lhs[j].value.asDouble() * integers[integers_lim[old_index]].org_lb.asDouble());
		  }
		  new_lhs.push_back(new_lhs_elem);
        	} else {
		  double coeff = old_lhs[j].value.asDouble();
		  for (int k = integers[integers_lim[old_index]].bitcnt-1, k2=0; k >= 0; k--, k2++) {
		    new_lhs_elem.index = integers[integers[integers_lim[old_index]].pt2leader+k2].index;
		    assert(integers[integers[integers_lim[old_index]].pt2leader+k2].index == integers[integers_lim[old_index]].pt2leader+k2);
		    new_lhs_elem.value = coeff * ipow2(k);
		    new_lhs.push_back(new_lhs_elem);
		  }
		  new_rhs.setValue(new_rhs.getValue() - coeff * integers[integers_lim[old_index]].org_lb.asDouble());
        	}
        }
        // add constraint (new_lhs,new_rhs) to qbp
		data::Constraint& c = qbp.createRhsConstraint(new_rhs);
        c.setElements(new_lhs);
    }
	// end: convert and add constraints

	// convert and add objective
	//if( qmip.getObjective() == data::QpObjFunc::max ) qmip.reverseObjFunc();
    std::vector<data::IndexedElement> old_obj = qmip.getObjectiveFunction().getObjectiveElementsSparse();
    data::QpNum objOffSet=qmip.getObjFuncOffset();
	std::vector<data::IndexedElement> new_obj;
	double offset_movement = 0.0;
	for (int j = 0; j < old_obj.size();j++) {
		int old_index = old_obj[j].  index;
		int new_index = integers[integers_lim[old_index]].index;
		data::IndexedElement new_obj_elem;
		new_obj_elem.index = new_index;
		new_obj_elem.value = old_obj[j].value;
		if (integers[integers_lim[old_index]].index == integers[integers_lim[old_index]].pt2leader &&
			integers[integers_lim[old_index]].bitcnt == 1) {
		        if (varVec[ old_index ]->getNumberSystem()!=data::QpVar::real && varVec[ old_index ]->getNumberSystem()!=data::QpVar::binaries)
			  offset_movement = offset_movement + old_obj[j].value.asDouble() * integers[integers_lim[old_index]].org_lb.asDouble();
			new_obj.push_back(new_obj_elem);
		} else {
			double coeff = old_obj[j].value.asDouble();
			for (int k = integers[integers_lim[old_index]].bitcnt-1, k2=0; k >= 0; k--, k2++) {
				new_obj_elem.index = integers[integers[integers_lim[old_index]].pt2leader+k2].index;
				assert(integers[integers[integers_lim[old_index]].pt2leader+k2].index == integers[integers_lim[old_index]].pt2leader+k2);
				new_obj_elem.value = coeff * ipow2(k);
				new_obj.push_back(new_obj_elem);
			}
			offset_movement = offset_movement + coeff * integers[integers_lim[old_index]].org_lb.asDouble();
		}
	}

	qbp.setObjective(qmip.getObjective());

	data::QpObjFunc& nn_obj = qbp.getObjectiveFunction();
	nn_obj.setSize(new_varVec.size());
	for (int j=0; j < new_varVec.size();j++) {
	  qbp.setObjectiveFunctionElement(j, 0.0);
	}
	for (int j=0; j < new_obj.size();j++) {
	  int index = new_obj[j].index;
	  double value = new_obj[j].value.asDouble();
	  qbp.setObjectiveFunctionElement(index, value);
	}
	// end: convert and add objective
	offset_shift = offset_movement;
	if (info_level >= 2) std::cerr << "OFFSET_SHIFT IN OBJ:" << offset_shift << " " << objOffSet << std::endl;
	//this->qbp->setFinalOffset( objOffSet.asDouble());
	//finalize bounds
	std::vector<data::QpVar *> varVecBin = qbp.getVariableVector();
	for (int i = 0; i < varVec.size();i++) {
	  if (varVec[i]->getNumberSystem()==data::QpVar::binaries) {
	    //cerr << "change bnds for var " << i << endl;
	    int index = integers_lim[i];
	    //cerr << "change bounds for variable " << i << " with index " << index << " and name " << varVec[i]->getName() << endl;
	    bool sth_chg = false;
	    if (orgLBs[i] > varVecBin[index]->getLowerBound().asDouble()) {
	      varVecBin[index]->setLowerBound(orgLBs[i]);
	      sth_chg = true;
	    } 
	    if (orgUBs[i] < varVecBin[index]->getUpperBound().asDouble()) {
	      varVecBin[index]->setUpperBound(orgUBs[i]);
	      sth_chg = true;
	    } 
	    if (sth_chg && fabs(varVecBin[index]->getUpperBound().asDouble() - varVecBin[index]->getLowerBound().asDouble()) < 1e-6) {
	      varVecBin[index]->setUpperBound( min(1.0, varVecBin[index]->getUpperBound().asDouble() +1e-6) );
	      varVecBin[index]->setLowerBound( max(0.0, varVecBin[index]->getLowerBound().asDouble() -1e-6) );
	    }
	  } else if (0) {
	    int index = integers_lim[i];
	    if (fabs(varVecBin[index]->getUpperBound().asDouble() - varVecBin[index]->getLowerBound().asDouble()) < 1e-6) {
	      varVecBin[index]->setUpperBound( varVecBin[index]->getUpperBound().asDouble() +1e-7 );
	      varVecBin[index]->setLowerBound( varVecBin[index]->getLowerBound().asDouble() -1e-7 );
	    }
	  } 
	}
  return 0;
}
#else
void yInterface::makeBinary(data::Qlp &qmip, data::Qlp &qbp) {
	// fill the qbp with the variables
	std::vector<data::QpVar *> varVec = qmip.getVariableVector();
	//for (int z=0;z<varVec.size();z++) {
	//	cerr << varVec[z]->getName() << "(" << (int)varVec[z]->getNumberSystem() << "," << varVec[z]->getUpperBound() << ")";
	//}
	//cerr << endl;
    int numConstraints = qmip.getConstraintCount();
    std::vector<const data::QpRhs *> rhsVec = qmip.getRhsVecConst();
    std::vector<const data::Constraint *> conVec = qmip.getConstraintVecConst();
    
    // try to find sharper bounds for int and real
    for (int i = 0; i < numConstraints;i++) {
        data::QpRhs rhs = *rhsVec[i];
        std::vector<data::IndexedElement> lhs = conVec[i]->getElements();
        int cntReals=0;
        int cntInts=0;
        int jOfInterest=-1;
        data::QpNum upper=0.0, lower=0.0;
        for (int j = 0; j < lhs.size();j++) {
            int index = lhs[j].index;
            double coeff = lhs[j].value.asDouble();

            if (varVec[index]->getNumberSystem()==data::QpVar::binaries) {
                if (coeff > 0) lower = lower - coeff;     // "-" because lateron to rhs
                else if (coeff < 0) upper = upper - coeff;
            } else if (varVec[index]->getNumberSystem()==data::QpVar::real) {
                jOfInterest = j;
                cntReals++;
            } else {
                jOfInterest = j;
                cntInts++;
            }
        }
        if (cntInts + cntReals == 1) {
            lower = rhs.getValue() + lower;
            upper = rhs.getValue() + upper;
            int index = lhs[jOfInterest].index;
            double coeff = lhs[jOfInterest].value.asDouble();
            if (coeff >= 0) {
                lower = lower / coeff;
                upper = upper / coeff;
            } else {
                data::QpNum u = upper;
                upper = lower / coeff;
                lower = u / coeff;
            }
            lower = lower - LP_EPS - 1.0;//0001;// - 1e-5;
            upper = upper + LP_EPS + 1.0;//0001;// + 1e-5;
            if (/*varVec[index]->getNumberSystem()!=data::QpVar::binaries &&*/ varVec[index]->getNumberSystem()!=data::QpVar::real) {
                lower = ceil(lower.asDouble())-LP_EPS - 1e-9;
                upper = floor(upper.asDouble())+LP_EPS +1e-9;
            }
            if (rhs.getRatioSign() == data::QpRhs::equal) {
                if (lower > varVec[index]->getLowerBound() ) {
		  if (info_level > 2) cerr << "eq increase lower bound of x" << index << "from " << varVec[index]->getLowerBound().asDouble() << " to " << lower.asDouble() << endl;
                    varVec[index]->setLowerBound(lower);
                }
                if (upper < varVec[index]->getUpperBound() ) {
		  if (info_level > 2) cerr << "eq decrease upper bound of x" << index << "from " << varVec[index]->getUpperBound().asDouble() << " to " << upper.asDouble() << endl;
                    varVec[index]->setUpperBound(upper);
                }
            } else if ( ((rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual && coeff >= 0) || (rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual && coeff < 0)) ) {
                if (upper < varVec[index]->getUpperBound() ) {
		  if (info_level > 2) cerr << "decrease upper bound of x" << index << "from " << varVec[index]->getUpperBound().asDouble() << " to " << upper.asDouble() << endl;
                    varVec[index]->setUpperBound(upper);
                }
            } else if ( ((rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual && coeff >= 0) || (rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual && coeff < 0)) ){
                if (lower > varVec[index]->getLowerBound() ) {
		  if (info_level > 2) cerr << "increase lower bound of x" << index << "from " << varVec[index]->getLowerBound().asDouble() << " to " << lower.asDouble() << endl;
                    varVec[index]->setLowerBound(lower);
                }
            }
        }
    }
    // end: try to find sharper bounds for int and real
    
	int index = 0;
	for (int i=0; i < varVec.size(); i++) {
		int leader_index = index;
		integers_lim.push_back(index);
		assert(integers_lim[integers_lim.size()-1] == index);
		assert(integers.size() == index);
		if (varVec[i]->getNumberSystem()==data::QpVar::real || varVec[i]->getNumberSystem()==data::QpVar::binaries) {
			data::QpVar var = data::QpVar(varVec[i]->getName(), index,
					varVec[i]->getLowerBound(), varVec[i]->getUpperBound(),
					varVec[i]->getNumberSystem(), varVec[i]->getQuantifier());
			qbp.createVariable(var);
			IntInfo varInfo;
			varInfo.index = index;
			varInfo.pt2leader = leader_index;
			varInfo.org_lb = varVec[i]->getLowerBound();
			varInfo.org_ub = varVec[i]->getUpperBound();
			varInfo.org_ind = i;
			varInfo.bitcnt = 1;
			varInfo.name = var.getName();
			integers.push_back(varInfo);
			index++;
		} else {
            int num_bits = (int)( abs(log2(round(varVec[i]->getUpperBound().asDouble())-round(varVec[i]->getLowerBound().asDouble()))) ) + 1;
            if (num_bits > 40) num_bits = 40;
        	std::vector<data::IndexedElement> artif_lhs;
        	for (int j = num_bits-1; j >= 0; j--) {
    			data::QpVar var = data::QpVar(varVec[i]->getName() + "_b" + std::to_string(j), index,
    					0.0, 1.0,
    					data::QpVar::binaries, varVec[i]->getQuantifier());
    			qbp.createVariable(var);
    			IntInfo varInfo;
    			varInfo.index = index;
    			varInfo.pt2leader = leader_index;
    			varInfo.org_lb = varVec[i]->getLowerBound();
    			varInfo.org_ub = varVec[i]->getUpperBound();
    			varInfo.org_ind = i;
    			varInfo.name = var.getName();
    			varInfo.bitcnt = num_bits;
    			integers.push_back(varInfo);
            	data::IndexedElement artif_lhs_elem;
            	artif_lhs_elem.index = index;
                artif_lhs_elem.value = ipow2(j);
    			artif_lhs.push_back(artif_lhs_elem);
    			index++;
    			//if (info_level > 0) cerr << artif_lhs_elem.value << "x" << artif_lhs_elem.index << " + ";
            }
            // add artificial constraints sum a*2^i*x <= old_ub - old_lb
        	//if (info_level > 0) cerr << " <= " << varVec[i]->getUpperBound()-varVec[i]->getLowerBound() << ": vorher: lb=" << varVec[i]->getLowerBound() << " ub=" << varVec[i]->getUpperBound() << endl;
            data::Constraint& d = qbp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, min(varVec[i]->getUpperBound()-varVec[i]->getLowerBound(),MAXABSINT), data::QpRhs::EXISTENTIAL);
            d.setElements(artif_lhs);
            // end: add artificial constraints

		}
	}
	// end: fill the qbp with the variables
	std::vector<const data::QpVar *> new_varVec = qbp.getVariableVectorConst();
	//for (int z=0;z<new_varVec.size();z++)
	//	cerr << new_varVec[z]->getName() << " ";
	//cerr << endl;

	// convert and add constraints
    for (int i = 0; i < numConstraints;i++) {
    	data::QpRhs new_rhs = *rhsVec[i];
    	std::vector<data::IndexedElement> old_lhs = conVec[i]->getElements();
    	std::vector<data::IndexedElement> new_lhs;
        for (int j = 0; j < old_lhs.size();j++) {
        	int old_index = old_lhs[j].index;
        	int new_index = integers[integers_lim[old_index]].index;
        	data::IndexedElement new_lhs_elem;
        	new_lhs_elem.index = new_index;
        	new_lhs_elem.value = old_lhs[j].value;
        	if (integers[integers_lim[old_index]].index == integers[integers_lim[old_index]].pt2leader &&
        		integers[integers_lim[old_index]].bitcnt == 1) {
        		new_lhs.push_back(new_lhs_elem);
        	} else {
        		double coeff = old_lhs[j].value.asDouble();
                for (int k = integers[integers_lim[old_index]].bitcnt-1, k2=0; k >= 0; k--, k2++) {
                	new_lhs_elem.index = integers[integers[integers_lim[old_index]].pt2leader+k2].index;
                	assert(integers[integers[integers_lim[old_index]].pt2leader+k2].index == integers[integers_lim[old_index]].pt2leader+k2);
                	new_lhs_elem.value = coeff * ipow2(k);
                	new_lhs.push_back(new_lhs_elem);
                }
            	new_rhs.setValue(new_rhs.getValue() + coeff * integers[integers_lim[old_index]].org_lb.asDouble());
        	}
        }
        // add constraint (new_lhs,new_rhs) to qbp
		data::Constraint& c = qbp.createRhsConstraint(new_rhs);
        c.setElements(new_lhs);
    }
	// end: convert and add constraints

	// convert and add objective
	//if( qmip.getObjective() == data::QpObjFunc::max ) qmip.reverseObjFunc();
    std::vector<data::IndexedElement> old_obj = qmip.getObjectiveFunction().getObjectiveElementsSparse();
    data::QpNum objOffSet=qmip.getObjFuncOffset();
	std::vector<data::IndexedElement> new_obj;
	double offset_movement = 0.0;
	for (int j = 0; j < old_obj.size();j++) {
		int old_index = old_obj[j].  index;
		int new_index = integers[integers_lim[old_index]].index;
		data::IndexedElement new_obj_elem;
		new_obj_elem.index = new_index;
		new_obj_elem.value = old_obj[j].value;
		if (integers[integers_lim[old_index]].index == integers[integers_lim[old_index]].pt2leader &&
			integers[integers_lim[old_index]].bitcnt == 1) {
			new_obj.push_back(new_obj_elem);
		} else {
			double coeff = old_obj[j].value.asDouble();
			for (int k = integers[integers_lim[old_index]].bitcnt-1, k2=0; k >= 0; k--, k2++) {
				new_obj_elem.index = integers[integers[integers_lim[old_index]].pt2leader+k2].index;
				assert(integers[integers[integers_lim[old_index]].pt2leader+k2].index == integers[integers_lim[old_index]].pt2leader+k2);
				new_obj_elem.value = coeff * ipow2(k);
				new_obj.push_back(new_obj_elem);
			}
			offset_movement = offset_movement + coeff * integers[integers_lim[old_index]].org_lb.asDouble();
		}
	}
	qbp.setObjective(qmip.getObjective());
    for (int j=0; j < new_obj.size();j++) {
    	int index = new_obj[j].index;
    	double value = new_obj[j].value.asDouble();
    	qbp.setObjectiveFunctionElement(index, value);
    }
	// end: convert and add objective
	offset_shift=offset_movement;
    //offset_shift = objOffSet;//offset_movement;
    //this->qbp->setFinalOffset( objOffSet.asDouble());
}

#endif

