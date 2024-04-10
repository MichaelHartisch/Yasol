/*
*
* Solver: NbdStageSolver.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef NDMASTER_HPP_
#define NDMASTER_HPP_
#include "Algorithms/Algorithm.hpp"
#include "Datastructures/Datastructures.hpp"
#include "Utilities/QlpSplitter.hpp"

namespace algorithm {
class NbdMaster;
class NbdStageSolver {

	friend class NbdMaster;

public:

	NbdStageSolver(NbdMaster& nd, unsigned int depth);
	~NbdStageSolver();

	//Solve the node
	bool solveMaster(NbdTreeNode& node);
	//
	extSol::QpExternSolver::QpExtSolSolutionStatus getSolutionStatus();
	//
	data::QpNum getObjValue();
	//
	bool getValues(std::vector<data::QpNum>&);
	//
	bool getOptCut(NbdTreeNode& node, BendersCut&);
	//
	bool getFeasCut(NbdTreeNode& node, BendersCut&);

	//Set rhs with respect to the current scenario and proposals from higher stages
	void changeMasterRhs();
	//Change objective function offset (scenario specific if universal variable in objective function)
	void changeObjectiveOffSet(const data::QpNum&);
	//Changes lower bound of a variable
	void setVariableLB(unsigned int i, const data::QpNum& lb);
	//Changes upper bound of a vairbale
	void setVariableUB(unsigned int i, const data::QpNum& ub);
	//Loads a node into the stage solver (cut,base,...)
	void loadNode(NbdTreeNode& n, std::vector<data::QpNum>&);
	//Loads a node into the stage solver (cut,base,...)
	void loadSingleStageNode(NbdTreeNode& n, std::vector<data::QpNum>&);

	void reloadSingleStageNode(std::vector<data::QpNum>&);
	//Add a benders cut
	bool addBendersCut(NbdTreeNode& node, BendersCut&, std::vector<data::QpNum>&, bool checkRedundant, bool adaptRedundant);
	//Add a user cut
	bool addStageCut(unsigned int ucid, BendersCut&);
	//Remove user cut with given ucid
	void removeStageCut(unsigned int ucid);
	//Remove user cut with given ucid
	void removeStageCutsFromCut(unsigned int ucid);
	void removeStageAllConstraints();
	//Change rhs value of a user cut with given cutId
	void changeStageCutRhs(unsigned int ucid,const data::QpNum&);
	//Remove all user cuts
	void removeStageCuts();
	//Remove all cuts (benders cuts + user)
	void clearCuts();
	//Remove all cuts (benders cuts + user)
	void clearSingleStageCuts();
	//
	void clearStageSolver();
	//
	void reloadMasterRhs(unsigned int index);

	void resetBase();

	bool changeObjFuncCoeff(unsigned int index, const data::QpNum& coeff);

protected:

	void initialize();

	void reloadCuts(const std::vector<BendersCut>&,std::vector<data::QpNum>&);
	void addCut(const BendersCut&,std::vector<data::QpNum>&);

	void createCut(const NbdTreeNode&, BendersCut&, BendersCut::CutType, std::vector<data::QpNum>&);

	void modifyExistingBendersCut(NbdTreeNode& node, BendersCut&, unsigned int, bool equality);
	void modifyExistingStageCut(BendersCut&, unsigned int);

	void * getLpSolverObj();
	void * getLpSolverEnv();

	extSol::QpExternSolver& getExternSolver(){
		return *this->solver;
	}

	data::QpSparseMatrix& getMasterColumns(){
		return this->masterColumns;
	}

	data::QpSparseMatrix& getMasterCutColumns(){
		return this->masterCutColumns;
	}

	std::vector<data::QpNum>& getMasterRhs(){
		return this->currentMasterRhs;
	}

	std::vector<data::QpVar>& getMasterVariables(){
		return this->masterVariables;
	}

	//Hidden to avoid auto code from compiler
	NbdStageSolver(const NbdStageSolver&); //only declared to avoid auto creation by compiler
	NbdStageSolver& operator=(const NbdStageSolver&); //only declared to avoid auto creation by compiler

	void printCutStatus();

private:

	// Log string for debug output
	static std::string LOG_TAG;

	//Reference to ND Algorithm
	NbdMaster& nd;
	// The current depth in the QpTree of this NDMAster
	unsigned int masterDepth;
	//Number of variables at this stage
	unsigned int variables;
	//The number of original constraints at this stage
	unsigned int origConstraints;
	//Number of cuts currently added to the master
	unsigned int cuts;

	// Contains the index of the first and last variable from the exist matrix, this master if responsible to
	std::pair<int, int> thisStageVariableRange;
	// Contains the index of the first and last variable from the exist matrix, this master if responsible to
	std::pair<int, int> thisStageConstraintsRange;

	// Objective function parts of this stage
	data::QpObjFunc masterObjective;
	// Variables of this stage
	std::vector<data::QpVar> masterVariables;
	// Constains the Rhs when the proposal was added
	std::vector<data::QpNum> proposalMasterRhs;
	// Conatins the Rhs the original masterRhs of this stage
	std::vector<data::QpNum> currentMasterRhs;

	// Recourse part of matrix at this stage in column-wise format
	std::vector<std::vector<data::IndexedElement> > recourseMasterColumns;
	// Part of active matrix at this stage in column-wise format
	std::vector<std::vector<data::IndexedElement> > masterColumns;
	// Recourse part of cuts at this stage in column-wise format (changes for each node at this stage)
	std::vector<std::vector<data::IndexedElement> > recourseCutColumns;
	// Active parts of cuts at this stage in column-wise format (changes for each node at this stage)
	std::vector<std::vector<data::IndexedElement> > masterCutColumns;

	std::vector<data::QpNum> tmpQpNumVec;

	//Mapping from user-cut id to current constraint index (needed for partial removing of user cuts in other order as they were added)
	std::map<unsigned int, unsigned int> ucidToConIndexMap;

	//UserCuts that are currently added to this stage
	std::vector<BendersCut> uCutVec;

	//std::vector<UserCutColIndices > uCutColVec;

	//Pointer to the last node solved by this master
	NbdTreeNode* lastNode;

	//Pointer to the main LP solver
	extSol::QpExternSolver* solver;
	//Pointer to an exact LP solver, used to verify the inexact LP solution
	extSol::QpExternSolver* exactSolver;

};
}

#endif /*NDMASTER_HPP_*/
