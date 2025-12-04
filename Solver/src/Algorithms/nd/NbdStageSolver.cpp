/*
*
* Solver: NbdStageSolver.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Algorithms/nd/NbdStageSolver.hpp"
#include "Algorithms/nd/NbdMaster.hpp"
#include "ExternSolvers/QpExternSolvers.hpp"

namespace algorithm {

std::string NbdStageSolver::LOG_TAG = "NbdStage ";
NbdStageSolver::NbdStageSolver(NbdMaster& nd, unsigned int depth) :
		nd(nd), masterDepth(depth), variables(0), origConstraints(0), cuts(0), thisStageVariableRange(nd.qlpSplitter.existMatrixIndexAtDepth[depth]), thisStageConstraintsRange(nd.qlpSplitter.conIndAtDepth[depth]), masterObjective(), masterVariables(), proposalMasterRhs(), currentMasterRhs(), recourseMasterColumns(), masterColumns(), recourseCutColumns(), masterCutColumns(), tmpQpNumVec(), ucidToConIndexMap(), uCutVec(),/* uCutColVec(),*/lastNode(
		NULL), solver(NULL), exactSolver(NULL) {

	solver = extSol::initNbdExternSolver();

#if defined(EXACT_ARITHMETIC) && defined(NESTED_BENDERS_EXACT)
	exactSolver = extSol::initExactExternSolver();
#endif
	this->initialize();
}

NbdStageSolver::~NbdStageSolver() {
	if (this->solver)
		delete this->solver;
	if (this->exactSolver)
		delete exactSolver;
}

bool NbdStageSolver::solveMaster(NbdTreeNode& node) {

	lastNode = &node;
	lastNode->solved++;
	nd.subProbSolved++;

	extSol::QpExternSolver::QpExtSolSolutionStatus status;

	if (AUTO_SOLVE_FIRST) {
		this->solver->adaptToSolverMode(extSol::QpExternSolver::DEFAULT);
		this->solver->solve(1.0e+6, 1.0e+6);
		if (SOLVE_NODE_LP_TWICE) {
			this->solver->solve(1.0e+6, 1.0e+6);
		}
		if ((status = this->solver->getSolutionStatus()) == extSol::QpExternSolver::INFEASIBLE) {
			this->solver->adaptToSolverMode(extSol::QpExternSolver::NBD);
			this->solver->solve(1.0e+6, 1.0e+6);
		}
	} else {
		this->solver->solve(nd.simplexIterationLimit, nd.timeLimit);
		if (SOLVE_NODE_LP_TWICE) {
			this->solver->solve(1.0e+6, 1.0e+6);
		}
	}

	if ((status = this->solver->getSolutionStatus()) != extSol::QpExternSolver::INFEASIBLE && status != extSol::QpExternSolver::OPTIMAL) {
//		std::string s("Trying to repair solution at depth: ");
//		s+=utils::ToolBox::convertToString(node.depth);
//		s+=", NodeNr.: ";
//		s+=utils::ToolBox::convertToString(node.nodeNumber);
//		s+=", Status: ";
//		s+=extSol::QpExternSolver::solutionStatusToString(status);
//		utils::Logger::globalLog(utils::LOG_ERROR,LOG_TAG,s);
//		solver->writeToFile("/tmp/","error.sav");
//		solver->readFromFile("/tmp/error.sav");
//		this->resetBase();
//		this->solver->solve(1.0e+6, 1.0e+6);
//		status = this->solver->getSolutionStatus();
//		utils::Logger::globalLog(utils::LOG_ERROR,LOG_TAG,"Status after resolving   : "+extSol::QpExternSolver::solutionStatusToString(status));
	}

	//if we use exact arithmetic, verify the solution using an ecaxt solver
	if (this->exactSolver) {
		extSol::QpExternSolver::QpExtSolBase base;
		this->solver->getBase(base);
		this->exactSolver->setBase(base);
		if (this->solver->getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE) {
			std::vector<data::QpNum> tmp;
			this->solver->getDualFarkas(tmp);
			this->exactSolver->setRayGuess(tmp);
		}
		this->exactSolver->solve(1.0e+6, 1.0e+6);
		status = this->exactSolver->getSolutionStatus();
	}

	//some debug output
	if (LOG_ND_STAGE) {
		std::string s;
		for (unsigned int i = 0; i < node.depth; i++)
			s += "\t";
		s += "< ";
		s += utils::ToolBox::convertToString(node.depth);
		s += " , ";
		//	s += utils::ToolBox::convertToString(node.nodeNumber);
		s += " , ";
		s += utils::ToolBox::convertToString(node.solved);
		s += " , ";
		s += extSol::QpExternSolver::solutionStatusToString(this->solver->getSolutionStatus());
		if (this->solver->getSolutionStatus() != extSol::QpExternSolver::INFEASIBLE) {
			s += " , ";
			s += this->solver->getObjValue().toString();
		}
		s += " >";
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, s);
	}

	if (status == extSol::QpExternSolver::OPTIMAL || status == extSol::QpExternSolver::INFEASIBLE) {
		return true;
	} else if (status == extSol::QpExternSolver::OPTIMAL_INFEAS || status == extSol::QpExternSolver::NUM_BEST) {
		if (DISPLAY_EXTSOL_WARNINGS)
			std::cout << "Numerical Difficulites ( " + extSol::QpExternSolver::solutionStatusToString(status) << " ) at node: " << node.toString() << std::endl;
		return true;
	} else if (status == extSol::QpExternSolver::ABORT_IT_LIM || status == extSol::QpExternSolver::ABORT_TIME_LIM || status == extSol::QpExternSolver::ABORT_OBJ_LIM) {
		if (DISPLAY_EXTSOL_WARNINGS)
			std::cout << "Breaking due to: " + extSol::QpExternSolver::solutionStatusToString(status) << " at node: " << node.toString() << std::endl;
		return true;
	} else {
		if (DISPLAY_EXTSOL_WARNINGS)
			std::cout << "Unsupported solution status ( " + extSol::QpExternSolver::solutionStatusToString(status) << " ) at node: " << node.toString() << std::endl;
		return false;
	}

}


bool NbdStageSolver::changeObjFuncCoeff(unsigned int index, const data::QpNum& coeff){
	return this->solver->changeObjFuncCoeff(index,coeff);
}


extSol::QpExternSolver::QpExtSolSolutionStatus NbdStageSolver::getSolutionStatus() {
	extSol::QpExternSolver* s;
	exactSolver ? s = exactSolver : s = solver;
	return s->getSolutionStatus();
}

data::QpNum NbdStageSolver::getObjValue() {
	extSol::QpExternSolver* s;
	exactSolver ? s = exactSolver : s = solver;
	if (s->getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE)
		return data::QpNum(false);
	return s->getObjValue() + this->masterObjective.getOffset();
}

bool NbdStageSolver::getValues(std::vector<data::QpNum>& values) {
	extSol::QpExternSolver* s;
	exactSolver ? s = exactSolver : s = solver;
	if (s->getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE)
		return false;
	s->getValues(values);
	if(EXT_SOL_EPSILON)
	data::QpNum::removeSmallValuesFromVec(values, EXT_SOL_EPSILON);
	return true;
}

bool NbdStageSolver::getOptCut(NbdTreeNode& n, BendersCut& cut) {
	extSol::QpExternSolver* s;
	exactSolver ? s = exactSolver : s = solver;
	s->getExtendedDuals(this->tmpQpNumVec);
	if(EXT_SOL_EPSILON)
	data::QpNum::removeSmallValuesFromVec(tmpQpNumVec, EXT_SOL_EPSILON);
	this->createCut(n, cut, BendersCut::OPTIMALITY, tmpQpNumVec);
	return true;
}

bool NbdStageSolver::getFeasCut(NbdTreeNode& n, BendersCut& cut) {
	extSol::QpExternSolver* s;
	exactSolver ? s = exactSolver : s = solver;
	if (s->getSolutionStatus() != extSol::QpExternSolver::INFEASIBLE)
		return false;
	s->getExtendedDualFarkas(tmpQpNumVec, this->masterColumns, this->masterCutColumns);
	if(EXT_SOL_EPSILON)
	data::QpNum::removeSmallValuesFromVec(tmpQpNumVec, EXT_SOL_EPSILON);
	this->createCut(n, cut, BendersCut::FEASIBILITY, tmpQpNumVec);
	return true;
}

void NbdStageSolver::changeMasterRhs() {

	/*if (thisStageConstraintsRange.first != -1) {

	 unsigned int j = 0;
	 std::vector<unsigned int> ind;
	 std::vector<data::QpNum> rhsVals;
	 ind.reserve(this->origConstraints+this->cuts);
	 rhsVals.reserve(this->origConstraints+this->cuts);
	 for (int i = this->thisStageConstraintsRange.first; i <= thisStageConstraintsRange.second; i++, j++) {
	 if (proposalMasterRhs[j] != nd.currScenRhs[i]) {
	 proposalMasterRhs[j] = nd.currScenRhs[i];
	 if (proposalMasterRhs[j].isMaxInf() || proposalMasterRhs[j].isMinInf())
	 throw utils::AlgorithmException("changeMasterRhs -->MAX/MIN");
	 ind.push_back(j);
	 rhsVals.push_back(proposalMasterRhs[j]);
	 }
	 }
	 this->solver->updateModel();


	 this->solver->changeRhsElements(ind, rhsVals);

	 if (this->exactSolver){
	 this->exactSolver->updateModel();
	 this->exactSolver->changeRhsElements(ind, rhsVals);
	 }
	 }*/

	if (thisStageConstraintsRange.first != -1) {

		unsigned int j = 0;
		std::vector<int> ind;
		std::vector<double> rhsVals;
		ind.reserve(this->origConstraints + this->cuts);
		rhsVals.reserve(this->origConstraints + this->cuts);
		for (int i = this->thisStageConstraintsRange.first; i <= thisStageConstraintsRange.second; i++, j++) {
			if (proposalMasterRhs[j] != nd.currScenRhs[i]) {
				proposalMasterRhs[j] = nd.currScenRhs[i];
				//if (proposalMasterRhs[j].isMaxInf() || proposalMasterRhs[j].isMinInf())
				//	throw utils::AlgorithmException("changeMasterRhs -->MAX/MIN");
				ind.push_back(j);
				rhsVals.push_back(proposalMasterRhs[j].asDouble());
			}
		}
		this->solver->updateModel();

		this->solver->changeRhsElements(ind, rhsVals);

		if (this->exactSolver) {
			this->exactSolver->updateModel();
			this->exactSolver->changeRhsElements(ind, rhsVals);
		}
	}

}

void NbdStageSolver::reloadMasterRhs(unsigned int index) {
	if (thisStageConstraintsRange.first != -1) {
		if (index >= this->thisStageConstraintsRange.first && index <= thisStageConstraintsRange.second) {
			unsigned int ind = index - this->thisStageConstraintsRange.first;
			data::QpNum val(nd.qlpSplitter.originalRhs[index]);
			this->currentMasterRhs[ind] = val;
			this->solver->changeRhsElement(ind, val);
			if (this->exactSolver)
				this->exactSolver->changeRhsElement(ind, val);
		}
	}
}

void NbdStageSolver::resetBase() {
	extSol::QpExternSolver::QpExtSolBase base;
	base.constraints.resize(solver->getRowCount(), extSol::QpExternSolver::Basic);
	base.variables.resize(solver->getVariableCount(), extSol::QpExternSolver::FreeOrSuperbasic);
	solver->setBase(base);
}

void NbdStageSolver::changeObjectiveOffSet(const data::QpNum& off) {
	//TODO besser direkt in extern solver (kann aber nicht jeder)
	this->masterObjective.setOffset(off);
}

void NbdStageSolver::setVariableLB(unsigned int i, const data::QpNum& lb) {
	this->masterVariables[i].setLowerBound(lb);
	this->solver->setVarLB(i, lb);
	if (this->exactSolver)
		this->exactSolver->setVarLB(i, lb);
}

void NbdStageSolver::setVariableUB(unsigned int i, const data::QpNum& ub) {
	this->masterVariables[i].setUpperBound(ub);
	this->solver->setVarUB(i, ub);
	if (this->exactSolver)
		this->exactSolver->setVarUB(i, ub);
}

bool NbdStageSolver::addStageCut(unsigned int ucid, BendersCut& c) {
	this->uCutVec.push_back(c);
	this->ucidToConIndexMap.insert(std::make_pair(ucid, this->uCutVec.size() - 1));

	/*
	 unsigned int itmp = this->uCutColVec.size();
	 UserCutColIndices uci;
	 uCutColVec.push_back(uci);

	 uCutColVec[itmp].rrInd.resize(c.recourseRow.size());
	 uCutColVec[itmp].mrInd.resize(c.masterRow.size());

	 for(unsigned int i = 0; i < c.masterRow.size(); i++){
	 uCutColVec[itmp].mrInd[i].first=c.masterRow[i].index;
	 uCutColVec[itmp].mrInd[i].second=this->masterCutColumns[c.masterRow[i].index].size();
	 }

	 for(unsigned int i = 0; i < c.recourseRow.size(); i++){
	 uCutColVec[itmp].rrInd[i].first=c.recourseRow[i].index;
	 uCutColVec[itmp].rrInd[i].second=this->recourseCutColumns[c.recourseRow[i].index].size();
	 }*/

	if (this->nd.mode == NbdMaster::SINGLE_STAGE) {
		this->addCut(c, nd.propVec);
	}
	return true;
}

void NbdStageSolver::removeStageCut(unsigned int ucid) {

	throw utils::AlgorithmException("NbdStageSolver::removeStageCut(unsigned int ucid) --> depreceated");

//	std::map<unsigned int, unsigned int>::iterator it, end;
//	it = this->ucidToConIndexMap.find(ucid);
//	end = this->ucidToConIndexMap.end();
//	if (it == end)
//		throw utils::DataStructureException("NbdExtStageSolCplex::removeUserCut --> Cut to remove doies not exist: " + utils::ToolBox::convertToString(ucid));
//	unsigned int index = it->second;
//	while (it != end) {
//		it->second--;
//		it++;
//	}
//	this->ucidToConIndexMap.erase(ucid);
//	this->userCuts.erase(this->userCuts.begin() + index);
//	if (nd.mode == NbdMaster::SINGLE_STAGE){
//		this->proposalMasterRhs.erase(this->proposalMasterRhs.begin()+this->origConstraints+index);
//		this->solver->removeCut(this->origConstraints+index);
//		this->cuts--;
//	}
}

void NbdStageSolver::removeStageCutsFromCut(unsigned int ucid) {

	std::map<unsigned int, unsigned int>::iterator it, end;
	it = this->ucidToConIndexMap.find(ucid);
	end = this->ucidToConIndexMap.end();
	if (it == end) {
#ifdef SHOW_EXTERN_SOLVER_WARNINGS
	  std::cerr << "Error: NbdExtStageSolCplex::removeUserCutsFromCut --> Cut to remove does not exist" << std::endl;
#endif
	  return;
	  //throw utils::DataStructureException("NbdExtStageSolCplex::removeUserCutsFromCut --> Cut to remove doies not exist: " + utils::ToolBox::convertToString(ucid));
	}
	unsigned int index = it->second;

	this->ucidToConIndexMap.erase(it, end);

	if (nd.mode == NbdMaster::SINGLE_STAGE) {
		this->proposalMasterRhs.erase(this->proposalMasterRhs.begin() + this->origConstraints + index, this->proposalMasterRhs.end());
		this->currentMasterRhs.erase(this->currentMasterRhs.begin() + this->origConstraints + index, this->currentMasterRhs.end());

		this->solver->removeCutsFromCut(this->origConstraints + index);

		for (unsigned int i = 0; i < masterCutColumns.size(); i++) {

			int eraser = 0;
			//std::cout << "Column before: "<< data::indexedElementVecToString(masterCutColumns[i]) << std::endl;
			for (int j = masterCutColumns[i].size() - 1; j >= 0; j--) {
				if (masterCutColumns[i][j].index == this->origConstraints + index) {
					eraser = j;
					break;
				}
				if (masterCutColumns[i][j].index < this->origConstraints + index) {
					eraser = j + 1;
					break;
				}
			}
			//std::cout << "Erasing with: " << eraser << " for index " << index << " "<< this->origConstraints << std::endl;
			this->masterCutColumns[i].erase(this->masterCutColumns[i].begin() + eraser, this->masterCutColumns[i].end());
			//std::cout << "Column after : "<< data::indexedElementVecToString(masterCutColumns[i]) << std::endl;
		}

		/*for(unsigned int i = 0; i < uCutColVec[index].rrInd.size(); i++){
		 this->recourseCutColumns[uCutColVec[index].rrInd[i].first].erase(recourseCutColumns[uCutColVec[index].rrInd[i].first].begin() + uCutColVec[index].rrInd[i].second, recourseCutColumns[uCutColVec[index].rrInd[i].first].end());
		 }

		 for(unsigned int i = 0; i < uCutColVec[index].mrInd.size(); i++){
		 this->masterCutColumns[uCutColVec[index].mrInd[i].first].erase(this->masterCutColumns[uCutColVec[index].mrInd[i].first].begin()+uCutColVec[index].mrInd[i].second,this->masterCutColumns[uCutColVec[index].mrInd[i].first].end());
		 }*/

	}

	this->cuts -= uCutVec.size();
	this->uCutVec.erase(this->uCutVec.begin() + index, this->uCutVec.end());
	//TODO this->uCutColVec.erase(this->uCutColVec.begin() + index, this->uCutColVec.end());
	this->cuts += uCutVec.size();

}

void NbdStageSolver::removeStageAllConstraints(/*unsigned int ucid*/) {

	/*
	 *
	std::map<unsigned int, unsigned int>::iterator it, end;
	it = this->ucidToConIndexMap.find(ucid);
	end = this->ucidToConIndexMap.end();
	if (it == end)
		throw utils::DataStructureException("NbdExtStageSolCplex::removeUserCutsFromCut --> Cut to remove doies not exist: " + utils::ToolBox::convertToString(ucid));
	unsigned int index = it->second;
	*/

	this->ucidToConIndexMap.clear();

	if (nd.mode == NbdMaster::SINGLE_STAGE) {
		this->proposalMasterRhs.clear();//erase(this->proposalMasterRhs.begin() + this->origConstraints + index, this->proposalMasterRhs.end());
		this->currentMasterRhs.clear();//erase(this->currentMasterRhs.begin() + this->origConstraints + index, this->currentMasterRhs.end());

		this->solver->removeCutsFromCut(0/*this->origConstraints + index*/);

		for (unsigned int i = 0; i < masterCutColumns.size(); i++) {

			masterCutColumns[i].clear();
			this->masterColumns[i].clear();
			//int eraser = 0;
			//std::cout << "Column before: "<< data::indexedElementVecToString(masterCutColumns[i]) << std::endl;
			/*for (int j = masterCutColumns[i].size() - 1; j >= 0; j--) {
				if (masterCutColumns[i][j].index == this->origConstraints + index) {
					eraser = j;
					break;
				}
				if (masterCutColumns[i][j].index < this->origConstraints + index) {
					eraser = j + 1;
					break;
				}
			}
			*/
			//std::cout << "Erasing with: " << eraser << " for index " << index << " "<< this->origConstraints << std::endl;
			//this->masterCutColumns[i].erase(this->masterCutColumns[i].begin() + eraser, this->masterCutColumns[i].end());
			//std::cout << "Column after : "<< data::indexedElementVecToString(masterCutColumns[i]) << std::endl;
		}

		/*for(unsigned int i = 0; i < uCutColVec[index].rrInd.size(); i++){
		 this->recourseCutColumns[uCutColVec[index].rrInd[i].first].erase(recourseCutColumns[uCutColVec[index].rrInd[i].first].begin() + uCutColVec[index].rrInd[i].second, recourseCutColumns[uCutColVec[index].rrInd[i].first].end());
		 }

		 for(unsigned int i = 0; i < uCutColVec[index].mrInd.size(); i++){
		 this->masterCutColumns[uCutColVec[index].mrInd[i].first].erase(this->masterCutColumns[uCutColVec[index].mrInd[i].first].begin()+uCutColVec[index].mrInd[i].second,this->masterCutColumns[uCutColVec[index].mrInd[i].first].end());
		 }*/

	}
    this->origConstraints = 0;
	this->cuts = 0; //-= uCutVec.size();
	this->uCutVec.clear();//erase(this->uCutVec.begin() + index, this->uCutVec.end());
	//TODO this->uCutColVec.erase(this->uCutColVec.begin() + index, this->uCutColVec.end());
	////this->cuts += uCutVec.size();

}

void NbdStageSolver::changeStageCutRhs(unsigned int ucid, const data::QpNum& newRhs) {
	std::map<unsigned int, unsigned int>::iterator it, end;
	it = this->ucidToConIndexMap.find(ucid);
	end = this->ucidToConIndexMap.end();
	if (it == end)
		throw utils::DataStructureException("NbdExtStageSolCplex::removeUserCut --> Cut to remove does not exist: " + utils::ToolBox::convertToString(ucid));
	unsigned int index = it->second;
	this->uCutVec[index].rhs.setValue(newRhs);
	if (nd.mode == NbdMaster::SINGLE_STAGE) {
		this->proposalMasterRhs[this->origConstraints + index] = newRhs;
		this->currentMasterRhs[this->origConstraints + index] = newRhs;
		this->solver->changeRhsElement(this->origConstraints + index, newRhs);
	}
}

void NbdStageSolver::removeStageCuts() {
	//this->cuts -= this->uCutVec.size();
	this->uCutVec.clear();
	//this->uCutColVec.clear();
	this->ucidToConIndexMap.clear();
	if (nd.mode == NbdMaster::SINGLE_STAGE) {
		this->proposalMasterRhs.resize(origConstraints);
		this->currentMasterRhs.resize(origConstraints);
		this->recourseCutColumns.clear();
		this->recourseMasterColumns.clear();
		this->solver->removeCuts();
	}
}

void NbdStageSolver::loadNode(NbdTreeNode& n, std::vector<data::QpNum>& currProp) {

	if (!n.depth && !n.childVec) {
		if (!this->cuts && !this->uCutVec.size())
			return;
	}

	if (SAVE_NODE_CUTS) {
		//Save last base of previous node
		if (SAVE_BASE) {
			if (lastNode) {
				if (!lastNode->base) {
					lastNode->base = new extSol::QpExternSolver::QpExtSolBase();
				}
				if (this->exactSolver) {
					exactSolver->getBase(*lastNode->base);
				} else {
					solver->getBase(*lastNode->base);
				}
			}
		}
		//Clear cuts of previous node
		this->clearCuts();
		//Load user cuts
		this->reloadCuts(this->uCutVec, currProp);
		//Load cuts of current node
		if (n.scenarioCuts)
			this->reloadCuts(*n.scenarioCuts, currProp);
		//Load last base of current node
		solver->updateModel();
		if (SAVE_BASE && n.base) {
			if (n.base->constraints.size() < solver->getRowCount()) {
				n.base->constraints.resize(solver->getRowCount(), 1);
			}
			solver->setBase(*n.base);
		}
	} else {
		//Clear cuts form previous iteration in StageSolver
		this->clearCuts();
		//Clear cuts form previous iteration in previous node
		if (lastNode && lastNode->scenarioCuts)
			lastNode->scenarioCuts->clear();
	}
}

void NbdStageSolver::loadSingleStageNode(NbdTreeNode& n, std::vector<data::QpNum>& currProp) {
	//Clear cuts of previous node
	this->clearCuts();
	//Load user cuts
	this->reloadCuts(this->uCutVec, currProp);
	solver->updateModel();
}

void NbdStageSolver::reloadSingleStageNode(std::vector<data::QpNum>& currProp) {

//	this->lastNode = NULL;
////	if (!this->cuts)
////		return;
//	this->currentMasterRhs.resize(this->origConstraints);
//	this->proposalMasterRhs.resize(this->origConstraints);
//	for (unsigned int i = 0, size = this->recourseCutColumns.size(); i < size; i++) {
//		this->recourseCutColumns[i].clear();
//	}
//	for (unsigned int i = 0, size = this->masterCutColumns.size(); i < size; i++) {
//		this->masterCutColumns[i].clear();
//	}
//	this->cuts = 0;
//
//
//	for (unsigned int i = 0; i < this->userCuts.size(); i++){
//
//		for (unsigned int j = 0, size = this->userCuts[i].recourseRow.size(); j < size; j++) {
//			this->recourseCutColumns[this->userCuts[i].recourseRow[j].index].push_back(
//					data::IndexedElement(this->origConstraints + cuts, this->userCuts[i].recourseRow[j].value));
//		}
//
//		data::QpNum frac(this->userCuts[i].rhs.getValue());
//
//		this->currentMasterRhs.push_back(frac);
//		for (unsigned int j = 0, size = this->userCuts[i].masterRow.size(); j < size; j++) {
//			this->masterCutColumns[this->userCuts[i].masterRow[j].index - thisStageVariableRange.first].
//			push_back(data::IndexedElement(this->origConstraints + cuts, this->userCuts[i].masterRow[j].value));
//		}
//		for (unsigned int j = 0, size = this->userCuts[i].recourseRow.size(); j < size; j++) {
//			frac -= this->userCuts[i].recourseRow[j].value * currProp[this->userCuts[i].recourseRow[j].index];
//		}
//		proposalMasterRhs.push_back(frac);
//
//		std::vector<data::IndexedElement> tmpVec = this->userCuts[i].masterRow;
//		for (unsigned int j = 0; j < this->userCuts[i].masterRow.size(); j++)
//			tmpVec[j].index -= thisStageVariableRange.first;
//
//		this->solver->changeRhsElement(this->origConstraints+i,frac);
//
//		cuts++;
//
//	}

}

void NbdStageSolver::reloadCuts(const std::vector<BendersCut>& cuts, std::vector<data::QpNum>& currProp) {
	for (unsigned int i = 0; i < cuts.size(); i++)
		this->addCut(cuts[i], currProp);
}

void NbdStageSolver::addCut(const BendersCut& cut, std::vector<data::QpNum>& currentProposal) {

	for (unsigned int i = 0, size = cut.recourseRow.size(); i < size; i++) {
		this->recourseCutColumns[cut.recourseRow[i].index].push_back(data::IndexedElement(this->origConstraints + cuts, cut.recourseRow[i].value));
	}

	data::QpNum frac(cut.rhs.getValue());

	this->currentMasterRhs.push_back(frac);
	for (unsigned int i = 0, size = cut.masterRow.size(); i < size; i++) {
		this->masterCutColumns[cut.masterRow[i].index - thisStageVariableRange.first].push_back(data::IndexedElement(this->origConstraints + cuts, cut.masterRow[i].value));
	}
	for (unsigned int i = 0, size = cut.recourseRow.size(); i < size; i++) {
		frac -= cut.recourseRow[i].value * currentProposal[cut.recourseRow[i].index];
	}
	proposalMasterRhs.push_back(frac);

	std::vector<data::IndexedElement> tmpVec = cut.masterRow;
	for (unsigned int i = 0; i < cut.masterRow.size(); i++)
		tmpVec[i].index -= thisStageVariableRange.first;
	this->solver->addCut(tmpVec, cut.rhs.getRatioSign(), frac);
	if (this->exactSolver)
		this->exactSolver->addCut(tmpVec, cut.rhs.getRatioSign(), frac);
	cuts++;
}

bool NbdStageSolver::addBendersCut(NbdTreeNode& node, BendersCut& c, std::vector<data::QpNum>& currentProposal, bool checkRedundant, bool adaptRedundant) {

//Check if cut already exists
	if (c.cType == BendersCut::EMPTY) {
		return false;
	}

	if (checkRedundant) {

		BendersCut *cut;
		std::vector<BendersCut> &cutVec = *node.scenarioCuts;

		for (unsigned int i = 0; i < cutVec.size(); i++) {
			if ((cut = &cutVec[i])->rhs.getRatioSign() == c.rhs.getRatioSign() && cut->equalsLHS(c)) {
				if (cut->rhs.getRatioSign() == c.rhs.getRatioSign()) {
					if (cut->rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
						if ((cut->rhs.getValue() >= c.rhs.getValue())) {
							return false;
						} else if (adaptRedundant) {
							this->modifyExistingBendersCut(node, c, i, false);
							return true;
						}
					} else if (cut->rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
						if ((cut->rhs.getValue() <= c.rhs.getValue())) {
							return false;
						} else if (adaptRedundant) {
							this->modifyExistingBendersCut(node, c, i, false);
							return true;
						}
					} else {
						if ((cut->rhs.getValue() == c.rhs.getValue())) {
							return false;
						}
					}
				}
			}
		}
	}

	node.scenarioCuts->push_back(c);
	this->addCut(c, currentProposal);
	return true;
}

void NbdStageSolver::clearStageSolver() {
	this->lastNode = NULL;
	this->clearCuts();
}

void NbdStageSolver::clearCuts() {

	if (!this->cuts)
		return;
	this->currentMasterRhs.resize(this->origConstraints);
	this->proposalMasterRhs.resize(this->origConstraints);
	for (unsigned int i = 0, size = this->recourseCutColumns.size(); i < size; i++) {
		this->recourseCutColumns[i].clear();
	}
	for (unsigned int i = 0, size = this->masterCutColumns.size(); i < size; i++) {
		this->masterCutColumns[i].clear();
	}
	this->solver->removeCuts();
	if (this->exactSolver)
		this->exactSolver->removeCuts();
	this->cuts = 0;

}

void NbdStageSolver::printCutStatus() {

	if (!cuts)
		return;

	std::vector<unsigned int> tmp;
	for (unsigned int i = 0, size = this->recourseCutColumns.size(); i < size; i++) {
		tmp.push_back(this->recourseCutColumns[i].size());
		;
	}
	std::cout << "recourseCutColumns: " << utils::ToolBox::vecToString(tmp) << std::endl;
	tmp.clear();
	for (unsigned int i = 0, size = this->masterCutColumns.size(); i < size; i++) {
		tmp.push_back(this->masterCutColumns[i].size());
		;
	}
	std::cout << "masterCutColumns: " << utils::ToolBox::vecToString(tmp) << std::endl;
}

//---------------------------------- Initialization Routines --------------------------------->
void NbdStageSolver::initialize() {

	const data::QpSparseMatrix& ieM = this->nd.qlpSplitter.fastExistMatrix;
	data::QpSparseMatrix mM, rM;
	std::vector<data::QpRhs> origMasterRhs;								//can be removed
	unsigned int subProbApproxVariables = (this->masterDepth != nd.qpTree.size() - 1) ? 1 : 0;
	this->masterObjective.setObjective(nd.qlpWork.getObjective());
	this->masterObjective.setSize(1 + thisStageVariableRange.second - thisStageVariableRange.first + subProbApproxVariables);
	if (masterDepth == 0)
		this->masterObjective.setOffset(nd.qlpWork.getObjectiveFunction().getOffset());

	if (thisStageConstraintsRange.first == -1 && thisStageConstraintsRange.second == -1) {
	} else {
		data::getSubMatrix(ieM, mM, thisStageConstraintsRange.first, thisStageConstraintsRange.second, thisStageVariableRange.first, thisStageVariableRange.second);
		if (masterDepth) {
			data::getSubMatrix(ieM, rM, thisStageConstraintsRange.first, thisStageConstraintsRange.second, 0, thisStageVariableRange.first == 0 ? 0 : (thisStageVariableRange.first - 1));
		}
	}

	int index = 0;
	data::QpNum varValue;
	for (unsigned int varIndex = thisStageVariableRange.first; index < (1 + thisStageVariableRange.second - thisStageVariableRange.first); index++, varIndex++) {
		masterVariables.push_back(nd.qlpSplitter.existVars[varIndex]);
		if (!(varValue = nd.qlpSplitter.existObjVec[varIndex]).isZero()) {
			masterObjective.setObjElement(index, varValue);
		}
	}

	data::QpNum lb(APPROX_LOWER_BOUND /*/ pow(10.0, (int) masterDepth)*/), ub(APPROX_UPPER_BOUND /*/ pow(10.0, (int) masterDepth)*/);

	if (COMPUTE_APPROX_BOUNDS && this->masterDepth != nd.qpTree.size() - 1) {
		data::QpNum lowerApproxBound, upperApproxBound;
		unsigned int startIndex = thisStageVariableRange.second + 1;
		const std::vector<data::QpNum> & objElemVec = this->nd.qlpSplitter.qlp.getObjectiveFunction().getObjectiveElementsDense();
		std::vector<data::QpVar*>& objElemVars = this->nd.qlpSplitter.vars;
		for (unsigned int i = startIndex; i < objElemVec.size(); i++) {
			if (!objElemVec[i].isZero()) {
				lowerApproxBound += (objElemVec[i] * (objElemVec[i] < 0 ? objElemVars[i]->getUpperBound() : objElemVars[i]->getLowerBound()));
				upperApproxBound += (objElemVec[i] * (objElemVec[i] < 0 ? objElemVars[i]->getLowerBound() : objElemVars[i]->getUpperBound()));
			}
		}
		if (lowerApproxBound.asDouble() != (long int) lowerApproxBound.asDouble())
			lowerApproxBound -= 1.0;
		if (lowerApproxBound > lb)
			lb = ((long int) lowerApproxBound.asDouble());
		if (upperApproxBound.asDouble() != (long int) upperApproxBound.asDouble())
			upperApproxBound += 1.0;
		if (upperApproxBound < ub) {
			if (upperApproxBound > 1000) {
				ub = ((long int) upperApproxBound.asDouble());
			} else {
				ub = 1000;
			}
		}
		//std::cout << "\t Depth: " << masterDepth << "\t LB: ( " << lowerApproxBound.toString() << " \t UB: " << upperApproxBound.toString() <<" , "<< ub.toString() <<" ) " << std::endl;
	}

	for (unsigned int i = 0; i < subProbApproxVariables; i++, index++) {
		data::QpVar var(std::string("z") += utils::ToolBox::convertToString(i), index, lb, ub);
		masterVariables.push_back(var);
		masterObjective.setObjElement(index, 1.0 / subProbApproxVariables);
	}

//--------------------- Create the rhs vectors --------------------------------------->
	std::vector<data::QpRhs> origRhs = nd.qlpWork.getRhsVec();
	this->proposalMasterRhs.reserve(ELEMS_TO_RESERVE);
	this->currentMasterRhs.reserve(ELEMS_TO_RESERVE);
	if (thisStageConstraintsRange.second != -1) {
		int con = thisStageConstraintsRange.first;
		while (con <= thisStageConstraintsRange.second) {
			proposalMasterRhs.push_back(origRhs[con].getValue());
			currentMasterRhs.push_back(origRhs[con].getValue());
			origMasterRhs.push_back(origRhs[con]);
			con++;
		}
	}

//------------------------------------------------------------------------------------>
	this->origConstraints = origMasterRhs.size();
	this->variables = masterVariables.size();

	unsigned int masterVariables = (1 + thisStageVariableRange.second - thisStageVariableRange.first);
	if (subProbApproxVariables)
		masterVariables += 1;
	unsigned int recourseVariables = thisStageVariableRange.first;
	for (unsigned int i = 0; i < masterVariables; i++) {
		masterColumns.push_back(std::vector<data::IndexedElement>());
		masterColumns[i].reserve(ELEMS_TO_RESERVE);
		masterCutColumns.push_back(std::vector<data::IndexedElement>());
		masterCutColumns[i].reserve(ELEMS_TO_RESERVE);
	}
	for (unsigned int i = 0; i < mM.size(); i++) {
		for (unsigned int j = 0; j < mM[i].size(); j++) {
			masterColumns[mM[i][j].index].push_back(data::IndexedElement(i, mM[i][j].value));
		}
	}

	for (unsigned int i = 0; i < recourseVariables; i++) {
		recourseMasterColumns.push_back(std::vector<data::IndexedElement>());
		recourseMasterColumns[i].reserve(ELEMS_TO_RESERVE);
		recourseCutColumns.push_back(std::vector<data::IndexedElement>());
		recourseCutColumns[i].reserve(ELEMS_TO_RESERVE);
	}
	for (unsigned int i = 0; i < rM.size(); i++) {
		for (unsigned int j = 0; j < rM[i].size(); j++) {
			recourseMasterColumns[rM[i][j].index].push_back(data::IndexedElement(i, rM[i][j].value));
		}
	}

	solver->init(this->masterObjective, this->masterVariables, mM, origMasterRhs);

	if (!this->masterDepth) {
		solver->adaptToSolverMode(extSol::QpExternSolver::NBD);
	} else {
		solver->adaptToSolverMode(extSol::QpExternSolver::NBD);
	}
	if (exactSolver) {
		exactSolver->init(this->masterObjective, this->masterVariables, mM, origMasterRhs);
		exactSolver->adaptToSolverMode(extSol::QpExternSolver::NBD);
	}

	if (ELEMS_TO_RESERVE)
		tmpQpNumVec.reserve(origMasterRhs.size() + 2 * this->masterVariables.size() + ELEMS_TO_RESERVE);

}

void NbdStageSolver::createCut(const NbdTreeNode& node, BendersCut& cut, BendersCut::CutType t, std::vector<data::QpNum>& tmpVec) {

	if (!node.depth)
		throw utils::AlgorithmException("NbdStageSolver::createCut(...) --> Stage 0");

	data::QpNum r, tmp;
	cut.clear();
	if ((cut.cType = t) == BendersCut::OPTIMALITY) {
		(r = this->masterObjective.getOffset());
	}

	std::vector<data::QpNum> currentCutBase(currentMasterRhs.begin(), currentMasterRhs.end());
	if (thisStageConstraintsRange.first >= 0) {
		for (int i = thisStageConstraintsRange.first, j = 0; i <= thisStageConstraintsRange.second; i++, j++)
			currentCutBase[j] -= node.data ? node.data->scenarioRhs[j] : nd.currEventRhs[i];
	}

	for (unsigned int i = 0, size = currentCutBase.size(); i < size; i++) {
		if (!tmpVec[i].isZero()) {
			(tmp = tmpVec[i]) *= currentCutBase[i];
			r += tmp;
		}
	}

	unsigned int variables = this->masterVariables.size();
	for (unsigned int i = currentMasterRhs.size(), index = 0; i < tmpVec.size(); i++, index++) {
		if (tmpVec[i].isZero())
			continue;
		if (index >= variables) {
			(tmp = tmpVec[i]) *= this->masterVariables[index % variables].getUpperBound();
			r += tmp;
		} else {
			(tmp = tmpVec[i]) *= this->masterVariables[index].getLowerBound();
			r += tmp;
		}
	}

	//if (!r.isZero() && std::fabs(r.asDouble()) < ABS_BREAK_VAL)
	//	r.setZero();

	cut.rhs.set(data::QpRhs::greaterThanOrEqual, r);
	data::IndexedElement ie, ie_val;
	std::vector<data::IndexedElement> *tmpVec1;
	int firstIndex = nd.qlpSplitter.existMatrixIndexAtDepth[masterDepth - 1].first;
	for (unsigned int i = 0, size = this->recourseMasterColumns.size(); i < size; i++) {
		ie_val.value.setZero();
		tmpVec1 = &recourseMasterColumns[i];
		for (unsigned j = 0, size_inner = tmpVec1->size(); j < size_inner; j++) {
			if ((ie = tmpVec1->operator [](j)).value.isZero())
				continue;
			(tmp = ie.value) *= tmpVec[ie.index];
			ie_val.value += tmp;
		}

		tmpVec1 = &recourseCutColumns[i];
		for (unsigned int j = 0, size_inner = tmpVec1->size(); j < size_inner; j++) {
			if ((ie = tmpVec1->operator [](j)).value.isZero())
				continue;
			( tmp = ie.value) *= tmpVec[ie.index];
			ie_val.value += tmp;
		}

		if (!ie_val.value.isZero()) {
			//if (std::fabs(ie_val.value.asDouble()) < ABS_BREAK_VAL)
			//	continue;
			ie_val.index = i;
			if (((int) i < firstIndex)) {
				cut.recourseRow.push_back(ie_val);
			} else {
				cut.masterRow.push_back(ie_val);
			}
		}
	}

	if (cut.cType == BendersCut::FEASIBILITY) {
		if (!cut.recourseRow.size() && !cut.masterRow.size()) {
			cut.cType = BendersCut::EMPTY;
		}
	} else if (t == BendersCut::OPTIMALITY) {
		cut.masterRow.push_back(data::IndexedElement(thisStageVariableRange.first, 1.0));
	}

	if (NORMALIZE_NBD_CUTS && CUT_COEF_EPSILON) {
		if (std::fabs(cut.rhs.getValue().asDouble()) < CUT_COEF_EPSILON)
			cut.rhs.setValue(0);
		for (unsigned int i = 0; i < cut.recourseRow.size(); i++) {
			if (std::fabs(cut.recourseRow[i].value.asDouble()) < CUT_COEF_EPSILON)
				cut.recourseRow[i].value.setZero();
		}
		for (unsigned int i = 0; i < cut.masterRow.size(); i++) {
			if (std::fabs(cut.masterRow[i].value.asDouble()) < CUT_COEF_EPSILON)
				cut.masterRow[i].value.setZero();
		}
		if (cut.masterRow.size() + cut.recourseRow.size() == 1) {
			if (cut.masterRow.size()) {
				cut.rhs.operator /=(cut.masterRow[0].value);
				cut.masterRow[0].value = 1;
			} else {
				cut.rhs.operator /=(cut.recourseRow[0].value);
				cut.recourseRow[0].value = 1;
			}
		}
	}

}

void NbdStageSolver::modifyExistingBendersCut(NbdTreeNode& node, BendersCut& cut, unsigned int i, bool equal) {
	data::QpNum f(cut.rhs.getValue() * (equal ? -1 : 1));
	unsigned int index = this->origConstraints + i;
	node.scenarioCuts->operator [](i).rhs.setValue(f);
	if (masterDepth) {
		data::QpNum diff(this->currentMasterRhs[index] - this->proposalMasterRhs[index]);
		(proposalMasterRhs[index] = f - diff);
	} else {
		this->proposalMasterRhs[index] = f;
	}
	this->currentMasterRhs[index] = f;
	this->solver->changeRhsElement(index, this->proposalMasterRhs[index]);
	if (this->exactSolver)
		this->exactSolver->changeRhsElement(index, this->proposalMasterRhs[index]);
}

void NbdStageSolver::modifyExistingStageCut(BendersCut& cut, unsigned int i) {
	throw utils::AlgorithmException("modifyExistingUserCut(..) --> not yet implemented");
}

void * NbdStageSolver::getLpSolverObj() {
	if (this->solver) {
		return this->solver->getSolverModel();
	} else {
		return NULL;
	}
}

void * NbdStageSolver::getLpSolverEnv() {
	if (this->solver) {
		return this->solver->getSolverEnv();
	} else {
		return NULL;
	}
}

}

