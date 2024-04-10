/*
*
* Solver: NbdMaster.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef NESTEDDECOMPOSITION_HPP_
#define NESTEDDECOMPOSITION_HPP_
#include "Algorithms/Algorithm.hpp"
#include "Utilities/Timer.hpp"
#include "NbdStageSolver.hpp"
#include "NbdStructs.hpp"
#include "Utilities/QlpStageRelaxer.hpp"
#undef _ALPHA
namespace algorithm {
class NbdMaster: public Algorithm {

	friend class SAA;
	friend class NbdStageSolver;

public:

	typedef enum {
		LP, SINGLE_STAGE, TWO_STAGE, MULTI_STAGE
	} Mode;

	typedef enum {
		FAST_FORWARD, FAST_BACK, FAST_FORWARD_FAST_BACK, FAST_FORWARD_FAST_FEASIBLE_BACK
	} SequencingProtocol;

	NbdMaster(const data::Qlp&, SequencingProtocol = (SEQUENCING_PROTOCOL == 0 ? FAST_FORWARD : (SEQUENCING_PROTOCOL == 1 ? FAST_BACK : (SEQUENCING_PROTOCOL == 2 ? FAST_FORWARD_FAST_BACK : FAST_FORWARD_FAST_FEASIBLE_BACK))));

	virtual ~NbdMaster();

	QlpSolution solveQlp(SolutionCase);

	void setAdvancedStartInformation(const data::QpNum&, const data::QpNum&, const std::vector<data::QpNum>&);

	//------------------------ Some Methods used by QlpStageSolver for Yasol (Ulf's Solver) ------------------------------------------------>
	unsigned int addUserCuts(const std::vector<data::IndexedElement>&, const data::QpRhs::RatioSign&, const data::QpNum&);

	void removeUserCut(unsigned int cutId);

	void removeUserCutsFromCut(unsigned int cutId);

	void removeAllConstraints();

	void removeUserCuts();

	void clearCuts();

	void changeUserCutRhs(unsigned int cutId, const data::QpNum&);

	void changeRhsVal(unsigned int cIndex, const data::QpNum&);

	void setVariableBounds(unsigned int, const data::QpNum&, const data::QpNum&);

	void setVariableLB(unsigned int, const data::QpNum&);

	void setVariableUB(unsigned int, const data::QpNum&);

	const std::vector<BendersCut>& getFirstStageBendersCuts() const;

	const std::vector<BendersCut>& getFirstStageUserCuts() const;

	bool changeObjFuncCoeff(unsigned int index, const data::QpNum& coeff);


	//-------------------------------------------------------------------------------------------------------------------------------------->

	unsigned int getSubProbSolved() {
		return subProbSolved;
	}

	unsigned int getIterations() {
		return totalIterations;
	}

	const SolutionStatus& getSolutionStatus() const {
		return this->solutionStatus;
	}

	void setLpLimit(int limit) {
		this->lpLimit = limit < 0 ? 1.0e+9 : limit;
	}

	void setIterationLimit(int limit) {
		this->simplexIterationLimit = limit < 0 ? 1.0e+9 : limit;
	}

	void setTimeLimit(double limit) {
		this->timeLimit = limit < 0 ? (NBD_TIMER ? NBD_TIME_LIMIT : 1.0e+9) : limit;
	}

	std::pair<data::QpNum, data::QpNum> getBounds() const {
		return std::make_pair(this->boundVec[0].first, this->boundVec[0].second);
	}

	void * getRootLpSolverObj() {
		return this->masters[0]->getLpSolverObj();
	}


	void * getRootLpSolverEnv() {
		return this->masters[0]->getLpSolverEnv();
	}


	extSol::QpExternSolver& getExternSolver() {
			return this->masters[0]->getExternSolver();
	}

	extSol::QpExternSolver& getFirstStageExternSolver() {
		return this->masters[0]->getExternSolver();
	}

	const data::QpSparseMatrix& getFirstStageMasterColumns() {
		return this->masters[0]->getMasterColumns();
	}

	data::QpSparseMatrix& getFirstStageMasterCutColumns() {
		return this->masters[0]->getMasterCutColumns();
	}

	std::vector<data::QpNum>& getFirstStageMasterRhs() {
		return this->masters[0]->getMasterRhs();
	}

	std::vector<data::QpVar>& getFirstStageMasterVariables() {
		return this->masters[0]->getMasterVariables();
	}

	void setTwoStageParameters(bool stopInf, bool addInfOpt, bool addAll, bool addWc, bool mo) {
		_TWOSTAGE_STOP_INF_CUT = stopInf;
		_TWOSTAGE_ADD_INF_OPT_CUT = addInfOpt;
		_TWOSTAGE_ADD_ALL_CUTS = addAll;
		_TWOSTAGE_ADD_WC_OPT_CUT = addWc;
		_TWOSTAGE_MOVE_ORDERING = mo;
	}

	void setMultiStageParameters(bool stopInf, bool addInfOpt, bool addAll, bool addWc, bool mo, bool a, bool b, bool ub) {
		_MULTISTAGE_STOP_INF_CUT = stopInf;
		_MULTISTAGE_ADD_INF_OPT_CUT = addInfOpt;
		_MULTISTAGE_ADD_ALL_CUTS = addAll;
		_MULTISTAGE_ADD_WC_OPT_CUT = addWc;
		_MULTISTAGE_MOVE_ORDERING = mo;
		_ALPHA = a;
		_BETA = b;
		_COMP_REL_UB = ub;
	}


protected:

	//------------------------- Benders Decomposition main methods --------------------------->

	/**Single Stage Mode (LP, EA-Qlp or AE-Qlp)*/
	BendersSolution SingleStage(NbdTreeNode& node);

	/**The TwoStage Benders algorithm*/
	BendersSolution TwoStageUniversal(NbdTreeNode& node);
	BendersSolution TwoStageRandom(NbdTreeNode& node);

	/**The Nested Benders algorithm*/
	bool MultiStage(BendersSolution&, NbdTreeNode& node, const data::QpNum& alpha, const data::QpNum& beta, bool nonRedCut);
	bool MultiStageTmp(BendersSolution&, NbdTreeNode& node);
	BendersSolution NestedBendersAvg(NbdTreeNode& node, const data::QpNum& parentA, const data::QpNum& parentB);

	void updateProposalVector(unsigned int, const std::vector<data::QpNum>&);
	void updateProposalVector(unsigned int, const BendersSolution&);
	void updateScenarioRhs(const NbdTreeNode&);
	void updateProposalRhs(const NbdTreeNode&);

	bool isRedundant(std::vector<BendersCut>&, BendersCut&);
	bool checkUserBreak();

	int getExistVariableStage(unsigned int) const;
	int getExistVariableStageIndex(unsigned int) const;

	/**Initializes all data structures of this algorithm*/
	void checkInput() {
	}
	void restart(bool = true);

	// Print several status informations
	void printSolutionInformation() const;

	//Hidden, only declared not implemented
	NbdMaster(const NbdMaster&);
	NbdMaster& operator=(const NbdMaster&);

private:

	// Log string for debug output
	static std::string LOG_TAG;

	bool sampleAverageApproximation;

	//-------------------------------------------- Some informations about algorithm mode------------------------------------------------------------>
	// OneStage, TwoStage or MultiStage problem
	Mode mode;
	// The benders solutions status after termination of the algorithm
	SolutionStatus solutionStatus;
	// Solution mode (feasibility, worst case, average case)
	SolutionCase solutionCase;
	// Sequencing protocol ("Fast Forward" or "Fast Forward Fast Back")
	SequencingProtocol protocol;
	// The type of the input qlp
	data::Qlp::QlpType qlpType;

	//-------------------------------------------- The tree that is traversed - --------------------------------------------------------------------->
	// Tree of nodes that hold information needed during computation (e.g. scenarioRhs, cuts, base , ...)
	data::QpMatrix<NbdTreeNode>::Type qpTree;
	//Holds single parts of QLP and several meta-information abount the qlp
	utils::QlpSplitter qlpSplitter;

	//-------------------------------------------- Several Vectors that change during the solution process ------------------------------------------>
	/** This vector contains the values of the current rhs at a specific node in the scenario QpTree  */
	std::vector<data::QpNum> currEventRhs; //currEventRhs[i] += row[j].value * this->qlpSplitter.currScenVarString[row[j].index];
	/** The current rhs after proposal rhs and scenario rhs have been substracted from the original rhs*/
	std::vector<data::QpNum> currScenRhs;
	/** This vector contains the values that are substituted and subtracted form rhs due to the current scenario (univ matrix) */
	std::vector<data::QpNum> currPropRhs; //currPropRhs[i] += row[j].value * this->propVec[row[j].index];
	/** Contains current proposals of upper stages at active path*/
	std::vector<data::QpNum> propVec;

	//-------------------------------------------- Stage Solver that solves the underlying LPs ------------------------------------------------------>
	/** Each level of the scenario QpTree has its own NbdExtStageSolCplex */
	std::vector<algorithm::NbdStageSolver*> masters;
	utils::QlpStageRelaxer* relaxer;
	/** Current lower and upper bound at each stage */
	std::vector<std::pair<data::QpNum, data::QpNum> > boundVec;

	//-------------------------------------------- Advanced Start Information ----------------------------------------------------------------------->
	//lower and upper bound (if the last exists)
	data::QpNum lbVal, ubVal;
	// lower or upper bound variable allocation for wartm start
	std::vector<data::QpNum> varAlloc;
	//scaling factor objective function scaling
	data::QpNum scaleFactor; //(not yet used)

	//-------------------------------------------- Status Information about the solution process----------------------------------------------------->
	unsigned int totalIterations, subProbSolved, relaxProbSolved;

	unsigned int alphaCuts, betaCuts, ipcuts, deepAlphaCuts, deepBetaCuts;

	long int lpLimit;

	long int simplexIterationLimit;

	double timeLimit;

	unsigned int userCutId;

	unsigned int userCuts;

	std::map<unsigned int, unsigned int> userCutIdToStageMap;

	BendersSolution nbdSol;

	utils::Timer timer;

	unsigned int currentDepth;

	bool _TWOSTAGE_STOP_INF_CUT;
	bool _TWOSTAGE_ADD_ALL_CUTS;
	bool _TWOSTAGE_ADD_WC_OPT_CUT;
	bool _TWOSTAGE_ADD_INF_OPT_CUT;
	bool _TWOSTAGE_MOVE_ORDERING;

	bool _MULTISTAGE_STOP_INF_CUT;
	bool _MULTISTAGE_ADD_ALL_CUTS;
	bool _MULTISTAGE_ADD_WC_OPT_CUT;
	bool _MULTISTAGE_ADD_INF_OPT_CUT;
	bool _MULTISTAGE_MOVE_ORDERING;
	bool _ALPHA;
	bool _BETA;
	bool _COMP_REL_UB;

};
}

#endif /*NESTEDDECOMPOSITION_HPP_*/
