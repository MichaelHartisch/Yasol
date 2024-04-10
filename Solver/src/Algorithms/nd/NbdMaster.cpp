/*
*
* Solver: NbdMaster.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Algorithms/nd/NbdMaster.hpp"
#include "Utilities/QlpConverter.hpp"
#include "Utilities/QlpSplitter.hpp"

namespace algorithm {

std::string NbdMaster::LOG_TAG = "NbdMAster";
NbdMaster::NbdMaster(const data::Qlp& qp, SequencingProtocol s) :
		Algorithm(qp, this->nested_benders), sampleAverageApproximation(false), mode(), solutionStatus(UNKNOWN), solutionCase(), protocol(s), qlpType(data::Qlp::QLP_TYPE_ERROR), qpTree(), qlpSplitter(this->qlpWork), currEventRhs(), currScenRhs(), currPropRhs(), propVec(), masters(),
		relaxer(), boundVec(), lbVal(true), ubVal(false), varAlloc(), scaleFactor(1), totalIterations(0), subProbSolved(0), relaxProbSolved(0), alphaCuts(0), betaCuts(0), ipcuts(0), deepAlphaCuts(0), deepBetaCuts(0), lpLimit(1.0e+9), simplexIterationLimit(
				1.0e+9), timeLimit(NBD_TIMER ? NBD_TIME_LIMIT : 1.0e+9), userCutId(1), userCuts(0), userCutIdToStageMap(), nbdSol(), timer(), currentDepth(1), _TWOSTAGE_STOP_INF_CUT(TS_STOP_NON_RED_INF_CUT), _TWOSTAGE_ADD_ALL_CUTS(TS_ADD_ALL_CUTS), _TWOSTAGE_ADD_WC_OPT_CUT(
				TS_ADD_WC_OPT_CUT), _TWOSTAGE_ADD_INF_OPT_CUT(TS_ADD_INF_OPT_CUT), _TWOSTAGE_MOVE_ORDERING(TS_MOVE_ORDERING), _MULTISTAGE_STOP_INF_CUT(MS_STOP_INF_SUB_PROB), _MULTISTAGE_ADD_ALL_CUTS(MS_ADD_ALL_CUTS), _MULTISTAGE_ADD_WC_OPT_CUT(
				MS_ADD_WC_OPT_CUT), _MULTISTAGE_ADD_INF_OPT_CUT(MS_ADD_INF_OPT_CUT), _MULTISTAGE_MOVE_ORDERING(MS_MOVE_ORDERING_OPT), _BETA(MS_BETA), _COMP_REL_UB(MS_COMP_REL_UB) {
	_ALPHA =MS_ALPHA;

	//--------------------------------------------- CHECKING INPUT ----------------------------------------------->
	unsigned int univVars = qlpWork.getQuantifierCount(data::QpVar::all);
	unsigned int randVars = qlpWork.getQuantifierCount(data::QpVar::random);
	if (univVars && randVars) {
		throw utils::AlgorithmException("NbdMaster::NbdMaster(...) --> the number of universal variables exceeds the maximum in Settings.hpp");
	}

	if (univVars > MAX_UNIVERSAL_VARS || randVars > MAX_RANDOM_VARS) {
		throw utils::AlgorithmException("NbdMaster::NbdMaster(...) --> the number of universal or random variables exceeds the maximum in Settings.hpp");
	}

	unsigned int stages = this->qlpWork.getStageCount();
	if (stages == 1) {
		this->mode = SINGLE_STAGE;
	} else if (stages == 2) {
		this->mode = TWO_STAGE;
	} else {
		this->mode = MULTI_STAGE;
	}

	this->qlpType = this->qlpWork.getQlpType();

	if (qlpType == data::Qlp::LP || qlpType == data::Qlp::IP || qlpType == data::Qlp::MIP) {
		if (LOG_ND)
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Input not a quantified problem: " + data::Qlp::qlpTypeToString(qlpType));
	}

	if (qlpType != data::Qlp::QLP) {
		//TODO check if after first stage integral variables do exist
		//TODO check mixed universal quantification
	}

	//--------------------------------------------- INITIALIZATION ----------------------------------------------->
	if (this->qlpWork.getObjective() == data::QpObjFunc::max) {
		if (LOG_ND)
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Reverting objective function... (max -> min)");
		this->qlpWork.getObjectiveFunction().reverseObjectiveFuntion();
	}

	if (this->sampleAverageApproximation) {
		throw utils::AlgorithmException("NestedBendersDecomposition --> SAA not longer supported");
	} else {
		qlpSplitter.initSplitter(this->solutionCase);
		qlpSplitter.initializeNbdTree(qpTree);
		if (randVars)
			qlpSplitter.initializeNbdTreeScenarioData(qpTree);

	}

	currEventRhs = currScenRhs = currPropRhs = std::vector<data::QpNum>(qlpWork.getConstraintCount(), 0.0);
	propVec = std::vector<data::QpNum>(qlpWork.getQuantifierCount(data::QpVar::exists), 0.0);

	for (unsigned int i = 0; i < qlpSplitter.stages; i++) {
		boundVec.push_back(std::make_pair(data::QpNum(false), data::QpNum(true)));
	}

	for (unsigned int i = 0; i < qlpSplitter.stages; i++) {
		this->masters.push_back(new NbdStageSolver(*this, i));
	}

	if (MS_COMP_REL_UB && this->mode == MULTI_STAGE)
		this->relaxer = new utils::QlpStageRelaxer(qlpWork);

}

NbdMaster::~NbdMaster() {
	for (unsigned int i = 0; i < this->masters.size(); i++) {
		delete masters[i];
	}
	delete relaxer;
}

bool NbdMaster::changeObjFuncCoeff(unsigned int index, const data::QpNum& coeff){
	return this->masters[0]->changeObjFuncCoeff(index,coeff);
}

Algorithm::QlpSolution NbdMaster::solveQlp(SolutionCase s) {

	Algorithm::QlpSolution sol((this->solutionCase = s), qlpWork.getObjective());

	nbdSol.clear();

	if (this->mode == SINGLE_STAGE && !this->qlpWork.isQuantifiedProblem()) {
		nbdSol = this->SingleStage(qpTree[0][0]);
		sol.solution.ofVal = nbdSol.objFuncVal;
		sol.solution.status = nbdSol.status;
		sol.solution.varAlloc = nbdSol.varAlloc;
	} else if (this->mode == SINGLE_STAGE && this->qlpWork.isQuantifiedProblem()) {
		utils::QlpRelaxer rel(qlpWork);
		if (qlpWork.getVariableByIndex(0).getQuantifier() != data::QpVar::exists) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving A...AE...E-Problem...");
			if (rel.getLowerBound().isMinInf()) {
				this->solutionStatus = INFEASIBLE;
				sol.solution.status = extSol::QpExternSolver::INFEASIBLE;
			} else {
				this->solutionStatus = FEASIBLE;
				sol.solution.status = extSol::QpExternSolver::OPTIMAL;
				sol.solution.ofVal = rel.getLowerBound();
				sol.solution.varAlloc = rel.getLowerBoundVarAlloc();
			}
		} else {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving E...EA...A-Problem...");
			if (rel.getUpperBound().isMaxInf()) {
				this->solutionStatus = INFEASIBLE;
				sol.solution.status = extSol::QpExternSolver::INFEASIBLE;
			} else {
				this->solutionStatus = FEASIBLE;
				sol.solution.status = extSol::QpExternSolver::OPTIMAL;
				sol.solution.ofVal = rel.getUpperBound();
				sol.solution.varAlloc = rel.getUpperBoundVarAlloc();
			}
		}
		return sol;
	} else {

		if (this->totalIterations || this->subProbSolved) {
			this->restart();
		}

		if (this->mode == TWO_STAGE) {
			if (this->solutionCase == algorithm::Algorithm::AVERAGE_CASE) {
				nbdSol = this->TwoStageRandom(qpTree[0][0]);
			} else {
				nbdSol = this->TwoStageUniversal(qpTree[0][0]);
			}
		} else {
			if (this->solutionCase == algorithm::Algorithm::AVERAGE_CASE) {
				nbdSol = this->NestedBendersAvg(qpTree[0][0], data::QpNum(true), data::QpNum(false));
			} else {
				if (ITERATIVE_DEEPENING) {
					this->currentDepth = BOUNCING_STAGE;
					while (this->currentDepth < qpTree.size()) {
						this->MultiStage(nbdSol, qpTree[0][0], data::QpNum(true), data::QpNum(false), true);
						if ((currentDepth + 1) < qpTree.size()) {
							if ((currentDepth + BOUNCING_FACTOR) < qpTree.size()) {
								std::cout << "Increasing from " << currentDepth << " to " << currentDepth + BOUNCING_FACTOR << std::endl;
								currentDepth += BOUNCING_FACTOR;
							} else {
								std::cout << "Increasing from " << currentDepth << " to " << qpTree.size() - 1 << std::endl;
								currentDepth = qpTree.size() - 1;
							}
						} else {
							std::cout << "Increasing from " << currentDepth << " to " << currentDepth + 1 << std::endl;
							currentDepth++;
						}

					}
				} else {
					currentDepth = 100;
					this->MultiStage(nbdSol, qpTree[0][0], data::QpNum(true), data::QpNum(false), true);
					//this->MultiStageTmp(nbdSol, qpTree[0][0]);
				}

			}
		}

		if ((sol.solution.status = nbdSol.status) != extSol::QpExternSolver::INFEASIBLE) {
			if (this->solutionStatus == UNKNOWN) {
				this->solutionStatus = FEASIBLE;
			}
			sol.solution.ofVal = nbdSol.objFuncVal;

			for (unsigned int i = 0; i < qlpWork.getVariableCount(); i++) {
				if (qlpWork.getVariableByIndex(i).getQuantifier() != data::QpVar::exists)
					break;
				sol.solution.varAlloc.push_back(nbdSol.varAlloc[i]);
			}
		} else {
			if (this->solutionStatus == UNKNOWN) {
				this->solutionStatus = INFEASIBLE;
			} else if (this->solutionStatus == TIME_LIMIT) {
				sol.solution.status = extSol::QpExternSolver::ABORT_TIME_LIM;
			} else if (this->solutionStatus == LP_LIMIT) {
				sol.solution.status = extSol::QpExternSolver::ABORT_IT_LIM;
			}
		}
	}
	timer.stop();
	if (false) {
		unsigned int depth = this->qpTree.size() - 1;
		unsigned int width = this->qpTree[depth].size();
		unsigned int visited = 0;
		for (unsigned int i = 0; i < width; i++) {
			if (qpTree[depth][i].solved)
				visited++;
		}
		std::cout << "Visited: " << visited << " From: " << width << std::endl;
	}
	if (false)
		this->printSolutionInformation();
	return sol;
}

//----------------------------------------------------------- Nested Benders Decomposition algorithm ------------------------------------------------>
BendersSolution NbdMaster::SingleStage(NbdTreeNode& node) {

	//-------------------------- some preparations ------------------------------------------------------------->
	BendersSolution solution;
	this->boundVec[node.depth].first.setMinInf();
	this->boundVec[node.depth].second.setMaxInf();

	this->restart(false); //false = kein clearCuts() Aufruf

	/*
	for (unsigned int i = 0, size = this->masters[node.depth]->recourseCutColumns.size(); i < size; i++) {
		this->masters[node.depth]->recourseCutColumns[i].clear();
	}
	for (unsigned int i = 0, size = this->masters[node.depth]->masterCutColumns.size(); i < size; i++) {
		this->masters[node.depth]->masterCutColumns[i].clear();
	}

	this->masters[node.depth]->cuts = 0;

	for (unsigned int i = 0; i < this->masters[node.depth]->uCutVec.size(); i++) {

		const BendersCut& cut = this->masters[node.depth]->uCutVec[i];

		for (unsigned int j = 0, size = cut.recourseRow.size(); j < size; j++) {
			this->masters[node.depth]->recourseCutColumns[cut.recourseRow[j].index].push_back(data::IndexedElement(this->masters[node.depth]->origConstraints + this->masters[node.depth]->cuts, cut.recourseRow[j].value));
		}

		for (unsigned int j = 0, size = cut.masterRow.size(); j < size; j++) {
			this->masters[node.depth]->masterCutColumns[cut.masterRow[j].index - this->masters[node.depth]->thisStageVariableRange.first].push_back(
					data::IndexedElement(this->masters[node.depth]->origConstraints + this->masters[node.depth]->cuts, cut.masterRow[j].value));
		}

		this->masters[node.depth]->cuts++;

	}*/

	if (this->masters[node.depth]->solveMaster(node)) {
		if ((solution.status = this->masters[node.depth]->getSolutionStatus()) == extSol::QpExternSolver::INFEASIBLE) {
			this->solutionStatus = INFEASIBLE;
			solution.objFuncVal.setZero();
			return solution;
		} else {
			this->solutionStatus = FEASIBLE;
			this->boundVec[node.depth].first = this->boundVec[node.depth].second = solution.objFuncVal = this->masters[node.depth]->getObjValue();
			this->masters[node.depth]->getValues(solution.varAlloc);
		}

		if (solution.status == extSol::QpExternSolver::ABORT_IT_LIM) {
			this->solutionStatus = IT_LIMIT;
		}
		if (solution.status == extSol::QpExternSolver::ABORT_TIME_LIM) {
			this->solutionStatus = TIME_LIMIT;
		}
		if (solution.status == extSol::QpExternSolver::ABORT_OBJ_LIM) {
			this->solutionStatus = OBJ_LIMIT;
		}

	} else {
		//throw utils::AlgorithmException("Exception solving initial master problem");
		solution.status = extSol::QpExternSolver::QpExtSolSolutionStatus::ERROR;
	}

	return solution;
}

//----------------------------------------------------------- Nested Benders Decomposition algorithm ------------------------------------------------>
BendersSolution NbdMaster::TwoStageUniversal(NbdTreeNode & node) {

	//-------------------------- some preparations ------------------------------------------------------------->
	BendersSolution solution;
	extSol::QpExternSolver::QpExtSolSolutionStatus msStat, ssStat;
	std::vector<data::QpNum> proposal;
	data::QpNum msObj, ssObj, wcObj, z;
	BendersCut cut;

	bool infeasibleSubProblems = false, newCuts = false;
	bool advStart = this->varAlloc.size() ? true : false;
	int moveIndex = -1;
	std::vector<BendersCut> optCuts, feasCuts, infeasOptCuts;

	this->boundVec[node.depth].first.setMinInf();
	this->boundVec[node.depth].second.setMaxInf();

	//-------------------------- initially solve master problem ----------------------------------------------->
	this->masters[node.depth]->loadSingleStageNode(node, this->propVec);	//(!!! needed for QlpStageSolver user cuts)
	if (this->masters[node.depth]->solveMaster(node)) {
		if ((msStat = this->masters[node.depth]->getSolutionStatus()) == extSol::QpExternSolver::INFEASIBLE) {
			solution.status = extSol::QpExternSolver::INFEASIBLE;
			solution.objFuncVal.setZero();
			return solution;
		} else if (msStat == extSol::QpExternSolver::OPTIMAL) {
			msObj = this->masters[node.depth]->getObjValue();
			this->masters[node.depth]->getValues(proposal);
		}
	} else {
		throw utils::AlgorithmException("Exception solving initial master problem");
	}

	unsigned int sApprInd = proposal.size() - 1;

	//-------------------------- set an advanced start solution ----------------------------------------------->
	if (advStart) {
		msObj = 0;
		for (unsigned int i = 0, size = proposal.size() - 1; i < size; i++) {
			proposal[i] = this->varAlloc[i];
			msObj += (proposal[i] * this->qlpWork.getObjectiveFunctionElement(i));
		}
		proposal[sApprInd] = (lbVal - msObj);
	}

	//--------------- start loop and run until LB==UB || no new cuts || master becomes infeasible ------------->
	do {
		//----------------------------------------- Next iteration -------------------------------------------->

		this->totalIterations++;

		(z = msObj) -= (this->boundVec[node.depth].first = proposal[sApprInd]);

		if (checkUserBreak())
			break;

		this->updateProposalVector(node.depth, proposal);

		this->updateProposalRhs(*node.childVec->operator[](0));

		if (LOG_ND) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Iteration: " + utils::ToolBox::convertToString(this->totalIterations));
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "LB       : " + (this->boundVec[node.depth].first + z).toString(10));
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "UB       : " + (this->boundVec[node.depth].second + z).toString(10));
		}

		infeasibleSubProblems = newCuts = false;
		feasCuts.clear();
		infeasOptCuts.clear();
		optCuts.clear();
		wcObj.setMinInf();

		//------------------------------------------ Solve all Subproblems ------------------------------------>
		if (LOG_ND)
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving Subproblems ...");

		for (unsigned int i = 0, successors = node.moveOrder->size(); i < successors; i++) {

			if (checkUserBreak())
				break;

			if ((moveIndex = node.moveOrder->operator [](i)) < 0)
				continue;

			this->qlpSplitter.setScenarioBitVector(*node.childVec->operator[](moveIndex), moveIndex);
			this->qlpSplitter.setScenarioVariableVectorByBitVector(*node.childVec->operator[](moveIndex));

			if (qlpSplitter.univObjVec.size() || qlpSplitter.randObjVec.size()) {
				this->masters[node.depth + 1]->changeObjectiveOffSet(this->qlpSplitter.precomputeSzenarioOffSet(node.depth, qlpSplitter.currScenVarString));
			}

			this->updateScenarioRhs(*node.childVec->operator[](moveIndex));

			this->masters[node.depth + 1]->changeMasterRhs();

			masters[node.depth + 1]->loadNode(*node.childVec->operator[](moveIndex), this->propVec);	//needed for usercuts

			if (LOG_ND) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving subproblem: " + utils::ToolBox::convertToString(i));
			}

			if (!this->masters[node.depth + 1]->solveMaster(*(node.childVec->operator[](moveIndex)))) {
				throw utils::AlgorithmException("Error Solving Subproblem: " + utils::ToolBox::convertToString(moveIndex));
			}
			if ((ssStat = this->masters[node.depth + 1]->getSolutionStatus()) != extSol::QpExternSolver::INFEASIBLE) {

				if (!_TWOSTAGE_ADD_ALL_CUTS && infeasibleSubProblems)
					continue;

				if (wcObj < (ssObj = this->masters[node.depth + 1]->getObjValue())) {

					if (LOG_ND) {
						utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "New worst-case solution...");
					}

					wcObj = ssObj;

					if (_TWOSTAGE_ADD_WC_OPT_CUT) {
						this->masters[node.depth + 1]->getOptCut(*(node.childVec->operator[](moveIndex)), cut);
						if (!this->isRedundant(optCuts, cut) && !node.isRedundant(cut)) {
							newCuts = true;
							if (optCuts.size()) {
								optCuts[0] = cut;
							} else {
								optCuts.push_back(cut);
							}
							if (_TWOSTAGE_MOVE_ORDERING) {
								node.pushToFront(i);
							}
						} else {
							if (LOG_ND) {
								utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Optimality Cut redundant: " + cut.toString());
							}
						}
					}
				}

				if (!_TWOSTAGE_ADD_WC_OPT_CUT) {
					if (proposal[sApprInd] > ssObj)
						continue;
					this->masters[node.depth + 1]->getOptCut(*(node.childVec->operator[](moveIndex)), cut);

					if (!this->isRedundant(optCuts, cut) && !node.isRedundant(cut)) {
						newCuts = true;
						optCuts.push_back(cut);
						if (_TWOSTAGE_MOVE_ORDERING) {
							node.pushToFront(i);
						}
					}
				}
			} else {

				this->masters[node.depth + 1]->getFeasCut(*(node.childVec->operator[](moveIndex)), cut);

				if (!this->isRedundant(feasCuts, cut) && !node.isRedundant(cut)) {
					newCuts = true;
					infeasibleSubProblems = true;
					feasCuts.push_back(cut);
					if (_TWOSTAGE_ADD_INF_OPT_CUT) {
						this->masters[node.depth + 1]->getOptCut(*(node.childVec->operator[](moveIndex)), cut);
						infeasOptCuts.push_back(cut);
					}
					if (_TWOSTAGE_MOVE_ORDERING) {
						node.pushToFront(i);
					}
					if (_TWOSTAGE_STOP_INF_CUT /*&& (!advStart || node.solved != 1)*/) {
						break;
					}
				} else {
					if (LOG_ND) {
						utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Feasibility Cut redundant: " + cut.toString());
					}
				}
			}
		}

		if (checkUserBreak())
			break;

		//------------- Check if we can stop due to matching bounds ------------------------------------------->
		if (infeasibleSubProblems) {
			//this->boundVec[node.depth].second.setMaxInf();
		} else {
			this->boundVec[node.depth].second = wcObj;
			if (this->solutionCase == FEASIBILITY)
				break;
			if (TS_CHECK_BOUNDS) {
				data::QpNum diff(this->boundVec[node.depth].second - this->boundVec[node.depth].first);
				if (diff < TS_BREAK_EPSILON)
					break;
				if ((this->boundVec[node.depth].second - this->boundVec[node.depth].first) < 0.5) {
					diff /= this->boundVec[node.depth].second;
					if (diff < 0)
						diff *= -1.0;
					if (diff < 1.0e-6) {
						break;
					}
				}
			}
		}

		//------------- If there are no new cuts we can stop -------------------------------------------------->
		if (!newCuts) {
			if (LOG_ND) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Break: No New Cuts...");
			}
			if (TS_BREAK_NO_NEW_CUTS) {
				break;
			}
		}

		//--------------------------- Adding cuts --------------------------------------------------------------->
		if (LOG_ND) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Adding Cuts...");
		}
		if (infeasibleSubProblems) {
			for (unsigned int i = 0; i < feasCuts.size(); i++) {
				if (LOG_ND)
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Adding cut: " + feasCuts[i].toString());
				this->masters[node.depth]->addBendersCut(node, feasCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS);
			}
			for (unsigned int i = 0; i < infeasOptCuts.size(); i++) {
				if (LOG_ND)
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Adding cut: " + infeasOptCuts[i].toString());
				this->masters[node.depth]->addBendersCut(node, infeasOptCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS);
			}
		}
		if (_TWOSTAGE_ADD_ALL_CUTS || !infeasibleSubProblems) {
			for (unsigned int i = 0; i < optCuts.size(); i++) {
				if (LOG_ND)
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Adding cut: " + optCuts[i].toString());
				if (this->masters[node.depth]->addBendersCut(node, optCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS)) {
					if (_TWOSTAGE_ADD_WC_OPT_CUT)
						break;
				}
			}
		}

		//------------- Resolve the master to obtain new candidate solution ----------------------------------->
		if (LOG_ND)
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Resolving Master...");
		if (!this->masters[node.depth]->solveMaster(node)) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Error solving Master LP.");
		} else if ((msStat = this->masters[node.depth]->getSolutionStatus()) == extSol::QpExternSolver::INFEASIBLE) {
			if (LOG_ND) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Break: Master infeasible");
			}
			break;
		} else {
			msObj = this->masters[node.depth]->getObjValue();
			this->masters[node.depth]->getValues(proposal);
		}

	} while (true);

	this->boundVec[node.depth].first = msObj;
	if (!this->boundVec[node.depth].second.isMaxInf())
		this->boundVec[node.depth].second += msObj - proposal[sApprInd];

	//Create Solution
	if (msStat != extSol::QpExternSolver::INFEASIBLE) {
		if (infeasibleSubProblems && this->subProbSolved < this->lpLimit) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, +"MasterSolution feasible but infeasible subproblems at depth 0.Value: " + this->masters[node.depth]->getObjValue().toString());
			solution.status = extSol::QpExternSolver::OPTIMAL;
			solution.objFuncVal = this->masters[node.depth]->getObjValue();
			this->masters[node.depth]->getValues(solution.varAlloc);
		} else {
			solution.status = extSol::QpExternSolver::OPTIMAL;
			solution.objFuncVal = this->masters[node.depth]->getObjValue();
			this->masters[node.depth]->getValues(solution.varAlloc);
		}
	} else {
		solution.status = extSol::QpExternSolver::INFEASIBLE;
		solution.objFuncVal.setZero();
		solution.varAlloc.clear();
	}

	return solution;
}

//----------------------------------------------------------- Nested Benders Decomposition algorithm ------------------------------------------------>
bool NbdMaster::MultiStage(BendersSolution& msSol, NbdTreeNode& node, const data::QpNum& parentAlpha, const data::QpNum& parentBeta, bool nonRedParentCut) {

//Start the next iteration
	this->totalIterations++;

//The solution of the master at this stage
	msSol.clear();

//Get the CplexWrapper that is responsible for this stage
	NbdStageSolver* master = this->masters[node.depth];

//Reinitialize bound vector
	this->boundVec[node.depth].first.setMinInf();
	this->boundVec[node.depth].second.setMaxInf();

//Inner Node, check for Alpha-Relaxtion-Cut
	if (_COMP_REL_UB && node.depth && node.childVec != NULL && !parentAlpha.isMinInf()) {

		//Setting exist part in relaxer
		for (unsigned int i = 0, size = this->qlpSplitter.existMatrixIndexAtDepth[node.depth - 1].second; i <= size; i++) {
			relaxer->setVariableFixation(node.depth, this->qlpSplitter.existVars[i].getIndex(), this->propVec[i]);
		}
		//Setting forall part in relaxer
		for (unsigned int i = 0, size = this->qlpSplitter.scenarioVarIndexAtDepth[node.depth]; i < size; i++) {	//TODO warum <
			relaxer->setVariableFixation(node.depth, this->qlpSplitter.univVars[i].getIndex(), this->qlpSplitter.currScenVarString[i]);
		}
		//Solve relaxation and check for deep cut
		data::QpNum val;
		this->relaxProbSolved++;
		data::QpNum offset = this->qlpSplitter.precomputeSzenarioOffSet(node.parent->depth, qlpSplitter.currScenVarString);
		this->subProbSolved++;

		bool feas = relaxer->solveStage(node.depth, val);

		for (unsigned int i = 0, size = this->qlpSplitter.existMatrixIndexAtDepth[node.depth - 1].second; i <= size; i++) {
			relaxer->detachFixation(node.depth, this->qlpSplitter.existVars[i].getIndex());
		}
		//Setting forall part in relaxer
		for (unsigned int i = 0, size = this->qlpSplitter.scenarioVarIndexAtDepth[node.depth]; i < size; i++) {
			relaxer->detachFixation(node.depth, this->qlpSplitter.univVars[i].getIndex());
		}

		if (feas && ((parentAlpha >= (val += offset)) || this->solutionCase == FEASIBILITY)) {
			deepAlphaCuts++;
			this->qlpSplitter.revertScenarioVariableVector(node);
			this->qlpSplitter.revertScenarioBitVector(node);
			msSol.status = extSol::QpExternSolver::OPTIMAL;
			msSol.objFuncVal = val;
			return true;
		}

	}

//Load the current node to the NbdExtStageSolCplex, if we are not at the root and not at a leaf
	master->loadNode(node, this->propVec);

//Load the current node int the NbdExtStageSolCplex
	if (node.depth) {
		//Check if the offset must be changed (can only happen if univ var coeff are in the objective function)
		if (qlpSplitter.univObjVec.size() || qlpSplitter.randObjVec.size()) {
			master->changeObjectiveOffSet(this->qlpSplitter.precomputeSzenarioOffSet(node.parent->depth, qlpSplitter.currScenVarString));
		}
		//Check if this master holds a submatrix of the original input
		if (this->qlpSplitter.conIndAtDepth[node.depth].second) {
			master->changeMasterRhs();
		}
	}

//Solve this nodes master problem
	if (!master->solveMaster(node) && node.depth) {
		utils::Logger::globalLog(utils::LOG_ERROR, LOG_TAG, +"Exception Solving Master at Node: " + node.toString());
		msSol.status = extSol::QpExternSolver::INFEASIBLE;
		msSol.objFuncVal.setZero();
		return false;
	}

//Check whether the master solution is initially feasible, if not directly go back
	if ((msSol.status = master->getSolutionStatus()) != extSol::QpExternSolver::INFEASIBLE) {

		//If this is a leaf node we are done and pass back the solution
		if (node.childVec != NULL && node.depth < this->currentDepth) {

			//Get master objective function and proposal
			msSol.objFuncVal = master->getObjValue();
			master->getValues(msSol.varAlloc);
			//Warm-Start
			if (node.depth) {
				if (node.parent->solved == 1 && node.solved == 1 && this->varAlloc.size()) {
					msSol.objFuncVal = this->ubVal;
					for (unsigned int i = 0, size = msSol.varAlloc.size(), offset = this->qlpSplitter.proposalIndexAtDepth[node.depth]; i < size; i++) {
						msSol.varAlloc[i] = this->varAlloc[offset + i];
					}
				}
			} else {
				if (node.solved == 1 && this->varAlloc.size()) {
					msSol.objFuncVal = this->ubVal;
					for (unsigned int i = 0, size = msSol.varAlloc.size(), offset = this->qlpSplitter.proposalIndexAtDepth[node.depth]; i < size; i++) {
						msSol.varAlloc[i] = this->varAlloc[offset + i];
					}
				}
			}

			//Some datastructures we need
			BendersSolution wcSubSol, subSol;
			data::QpNum msObj, z, phi, alpha, beta, offset;
			std::vector<BendersCut> optCuts, feasCuts, infeasOptCuts;
			bool infeasibleSubProblems = false, newCut = false;
			unsigned int subApproxIndex = msSol.varAlloc.size() - 1;

			do {

				if (this->protocol == FAST_BACK && node.depth) {
					master->getOptCut(node, msSol.opt);
					if (nonRedParentCut || !node.parent->isRedundant(msSol.opt))
						break;
				}

				if (this->checkUserBreak())
					break;

				//get master objective function value and split it into in master and subproblem part
				msObj = msSol.objFuncVal;
				phi = msSol.varAlloc[subApproxIndex];
				z = msObj - phi;

				this->boundVec[node.depth].first = msObj;

				this->updateProposalVector(node.depth, msSol.varAlloc);

				this->updateProposalRhs(*node.childVec->operator[](0));

				//-----------------AlphaBeta-Heuristic------------------------------------------------->
				alpha.setMinInf();
				beta.setMaxInf();
				if (!this->boundVec[node.depth].second.isMaxInf()) {
					(beta = this->boundVec[node.depth].second) -= z;
				}

				infeasibleSubProblems = newCut = false;

				int index = -1;

				wcSubSol.objFuncVal.setMinInf();

				//------------- Solve all subproblems ------------------------------------------------------------->
				for (unsigned int i = 0, successors = node.moveOrder->size(); i < successors; i++) {

					if (checkUserBreak())
						break;

					if ((index = node.moveOrder->operator [](i)) < 0) {
						continue;
					}

					this->qlpSplitter.setScenarioBitVector(*node.childVec->operator[](index), index);
					this->qlpSplitter.setScenarioVariableVectorByBitVector(*node.childVec->operator[](index));

					this->updateScenarioRhs(*node.childVec->operator[](index));

					//Get SubSolution
					if (!MultiStage(subSol, *node.childVec->operator[](index), alpha, beta, newCut)) {
						continue;
						throw utils::AlgorithmException("NestedBendersException");
					}

					if (subSol.status != extSol::QpExternSolver::INFEASIBLE) {

						if (alpha < subSol.objFuncVal) {
							alpha = subSol.objFuncVal;
						}

						if (!_MULTISTAGE_ADD_ALL_CUTS && infeasibleSubProblems)
							continue;

						if (wcSubSol.objFuncVal < subSol.objFuncVal) {

							wcSubSol = subSol;

							if (_MULTISTAGE_ADD_WC_OPT_CUT) {

								if (!this->isRedundant(optCuts, subSol.opt) && !node.isRedundant(subSol.opt)) {
									newCut = true;
									if (optCuts.size()) {
										optCuts[0] = subSol.opt;
									} else {
										optCuts.push_back(subSol.opt);
									}
									if (_MULTISTAGE_MOVE_ORDERING) {
										node.pushToFront(i);
									}
								}
							}

						}

						if (!_MULTISTAGE_ADD_WC_OPT_CUT) {
							if (phi > subSol.objFuncVal)
								continue;
							if (!this->isRedundant(optCuts, subSol.opt) && !node.isRedundant(subSol.opt)) {
								newCut = true;
								optCuts.push_back(subSol.opt);
								if (_MULTISTAGE_MOVE_ORDERING) {
									node.pushToFront(i);
								}
							}
						}

					} else {

						if (!this->isRedundant(feasCuts, subSol.feas) && !node.isRedundant(subSol.feas)) {
							newCut = true;
							infeasibleSubProblems = true;
							feasCuts.push_back(subSol.feas);
							if (_MULTISTAGE_ADD_INF_OPT_CUT)
								infeasOptCuts.push_back(subSol.opt);
							if (_MULTISTAGE_MOVE_ORDERING) {
								node.pushToFront(i);
							}
							if (_MULTISTAGE_STOP_INF_CUT) {
								break;
								if (node.depth) {
									if (node.parent->solved != 1 || node.solved != 1 || !this->varAlloc.size())
										break;
								} else {
									if (node.solved != 1 || !this->varAlloc.size())
										break;
								}
							}
						} else {
							if (LOG_ND) {
								utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "Redundant worst-case cut detected ...");
							}
						}
					}
				}

				//------------- Check if we can stop due to matching bounds or AlphaBeta Break -------------------->
				if (infeasibleSubProblems) {
					this->boundVec[node.depth].second.setMaxInf();
				} else {

					if (this->solutionCase == FEASIBILITY) {
						break;
					}

					(this->boundVec[node.depth].second = z) += wcSubSol.objFuncVal;

					if ((this->boundVec[node.depth].second < this->boundVec[node.depth].first)) {
						this->boundVec[node.depth].second = this->boundVec[node.depth].first;
					}

					if (MS_CHECK_BOUNDS) {

						if (((this->boundVec[node.depth].second - this->boundVec[node.depth].first) < MS_BREAK_EPSILON)) {
							//|| (((this->boundVec[node.depth].second - this->boundVec[node.depth].first) / (fabs(this->boundVec[node.depth].first.asDouble()) + 1e-10)) < MULTISTAGE_BREAK_EPSILON)) {
							if (node.depth) {
								master->getOptCut(node, msSol.opt);
								if (nonRedParentCut || !node.parent->isRedundant(msSol.opt))
									break;
							} else {
								break;
							}
						}
					}

					if (_BETA && node.depth && parentBeta <= this->boundVec[node.depth].first) {
						master->getOptCut(node, msSol.opt);
						if (nonRedParentCut || !node.parent->isRedundant(msSol.opt)) {
							this->betaCuts++;
							break;
						}
					}

					if (_ALPHA && node.depth && this->boundVec[node.depth].second <= parentAlpha) {
						master->getOptCut(node, msSol.opt);
						if (nonRedParentCut || !node.parent->isRedundant(msSol.opt)) {
							this->alphaCuts++;
							break;
						}
					}
				}

				//------------- If there are no new cuts we can stop -------------------------------------------------->
				if (MS_BREAK_NO_NEW_CUTS && !newCut) {
					if (LOG_ND)
						utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Break: No New Cuts");
					break;
				}

				//--------------------------- Adding optimality cuts ----------------------------------------------------->
				if (infeasibleSubProblems) {
					for (unsigned int i = 0; i < feasCuts.size(); i++) {
						this->masters[node.depth]->addBendersCut(node, feasCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS);
					}
					for (unsigned int i = 0; i < infeasOptCuts.size(); i++) {
						this->masters[node.depth]->addBendersCut(node, infeasOptCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS);
					}
				}

				if (_MULTISTAGE_ADD_ALL_CUTS || !infeasibleSubProblems) {
					for (unsigned int i = 0; i < optCuts.size(); i++) {
						if (this->masters[node.depth]->addBendersCut(node, optCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS)) {
							if (_MULTISTAGE_ADD_WC_OPT_CUT)
								break;
						}
					}
				}

				feasCuts.clear();
				infeasOptCuts.clear();
				optCuts.clear();

				msSol.clear();
				subSol.clear();
				wcSubSol.clear();

				//------------- Resolve the master to obtain new candidate solution ----------------------------------->
				if (this->masters[node.depth]->solveMaster(node)) {
					if ((msSol.status = master->getSolutionStatus()) == extSol::QpExternSolver::INFEASIBLE) {
						break;
					}
				} else {
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Error solving LP at depth: " + utils::ToolBox::convertToString(node.depth));
					msSol.status = extSol::QpExternSolver::INFEASIBLE;
					msSol.objFuncVal.setZero();
					return false;
				}

				//------------- Do protocol dependent checks ---------------------------------------------------------->
				if (this->protocol == FAST_FORWARD) {
					//nothing to do, just continue
				} else if (this->protocol == FAST_FORWARD_FAST_BACK) {
					if (node.depth) {
						master->getOptCut(node, msSol.opt);
						if (nonRedParentCut || !node.parent->isRedundant(msSol.opt)) {
							break;
						}
					}
				} else if (this->protocol == FAST_FORWARD_FAST_FEASIBLE_BACK) {
					if (!infeasibleSubProblems && node.depth) {
						master->getOptCut(node, msSol.opt);
						if (nonRedParentCut || !node.parent->isRedundant(msSol.opt)) {
							break;
						}
					}
				} else if (this->protocol == FAST_BACK) {
					if (node.depth) {
						master->getOptCut(node, msSol.opt);
						if (nonRedParentCut || !node.parent->isRedundant(msSol.opt)) {
							break;
						}
					}
				} else {
					throw utils::AlgorithmException("Protocol not yet implemented...");
				}

				msSol.objFuncVal = master->getObjValue();
				master->getValues(msSol.varAlloc);

			} while (true);

			//----------------------------------- Passing back solution (depends on sequencing protocol)----------------------------------------->
			//------------- Do protocol dependent checks ---------------------------------------------------------->
			if (this->protocol == FAST_FORWARD) {
				if (msSol.status != extSol::QpExternSolver::INFEASIBLE) {
					//msSol.objFuncVal = master->getObjValue();
					if (!this->boundVec[node.depth].second.isMaxInf()) {
						if (node.depth) {
							msSol.objFuncVal = this->boundVec[node.depth].second;
						} else {
							msSol.objFuncVal = master->getObjValue();
						}
					} else {
						if (node.depth) {
							msSol.objFuncVal = 1000000000;
						} else {
							msSol.objFuncVal = master->getObjValue();
						}
					}
					if (node.depth) {
						if (msSol.opt.cType == BendersCut::EMPTY)
							master->getOptCut(node, msSol.opt);
					} else {
						master->getValues(msSol.varAlloc);
					}
				} else {
					if (node.depth)
						master->getFeasCut(node, msSol.feas);
				}
			} else if (this->protocol == FAST_FORWARD_FAST_BACK) {
				if (msSol.status != extSol::QpExternSolver::INFEASIBLE) {
					if (!this->boundVec[node.depth].second.isMaxInf()) {
						if (node.depth) {
							msSol.objFuncVal = this->boundVec[node.depth].second;
						} else {
							msSol.objFuncVal = master->getObjValue();
						}

					} else {
						if (node.depth) {
							msSol.objFuncVal = 1000000000;
						} else {
							msSol.objFuncVal = master->getObjValue();
						}
					}
					if (node.depth) {
						if (msSol.opt.cType == BendersCut::EMPTY)
							master->getOptCut(node, msSol.opt);
					} else {
						master->getValues(msSol.varAlloc);
					}
				} else {
					if (node.depth)
						master->getFeasCut(node, msSol.feas);
				}
			} else if (this->protocol == FAST_FORWARD_FAST_FEASIBLE_BACK) {
				if (msSol.status != extSol::QpExternSolver::INFEASIBLE) {
					if (!this->boundVec[node.depth].second.isMaxInf()) {
						//msSol.objFuncVal = this->boundVec[node.depth].second;
						if (node.depth) {
							msSol.objFuncVal = this->boundVec[node.depth].second;
						} else {
							msSol.objFuncVal = master->getObjValue();
						}
					} else {
						if (node.depth) {
							msSol.objFuncVal = 1000000000;
						} else {
							msSol.objFuncVal = master->getObjValue();
						}
					}
					if (node.depth) {
						if (msSol.opt.cType == BendersCut::EMPTY)
							master->getOptCut(node, msSol.opt);
					} else {
						master->getValues(msSol.varAlloc);
					}
				} else {
					if (node.depth)
						master->getFeasCut(node, msSol.feas);
				}
			} else if (this->protocol == FAST_BACK) {
				if (msSol.status != extSol::QpExternSolver::INFEASIBLE) {
					if (!this->boundVec[node.depth].second.isMaxInf()) {
						msSol.objFuncVal = this->boundVec[node.depth].second;
					} else {
						if (node.depth) {
							msSol.objFuncVal = 1000000000;
						} else {
							msSol.objFuncVal = master->getObjValue();
						}
					}
					if (node.depth) {
						if (msSol.opt.cType == BendersCut::EMPTY)
							master->getOptCut(node, msSol.opt);
					} else {
						master->getValues(msSol.varAlloc);
					}
				} else {
					if (node.depth)
						master->getFeasCut(node, msSol.feas);
				}
			} else {
				throw utils::AlgorithmException("Protocol not yet implemented...");
			}

			//revert all changes that have been made by this node
			for (unsigned int i = 0; i < this->qlpSplitter.proposalsAtDepth[node.depth]; i++)
				this->propVec[this->qlpSplitter.proposalIndexAtDepth[node.depth] + i].setZero();
			//revert the scenario variable vector
			this->qlpSplitter.revertScenarioVariableVector(node);
			//revert the scenario bit string
			this->qlpSplitter.revertScenarioBitVector(node);
		} else {
			if (LOG_ND) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, +"Leaf Node, nothing more to do. Depth: " + utils::ToolBox::convertToString(node.depth));
			}
			if (msSol.status != extSol::QpExternSolver::INFEASIBLE) {
				msSol.objFuncVal = master->getObjValue();
				if (node.depth)
					master->getOptCut(node, msSol.opt);
			} else {
				if (node.depth)
					master->getFeasCut(node, msSol.feas);
			}
		}
	} else {
		if (LOG_ND) {
			utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, +"The Master is initially infeasible. Returning Cut.");
		}
		if (node.depth)
			master->getFeasCut(node, msSol.feas);
	}
	if (LOG_ND) {
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, +"Passing back solution: " + msSol.toString());
	}
	return true;
}

bool NbdMaster::checkUserBreak() {
	if (this->timer.getCurrentSeconds() > timeLimit) {
		if (LOG_ND) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Breaking due to Timeout. Timer: " + utils::ToolBox::convertToString(this->timer.getCurrentSeconds()));
		}
		this->solutionStatus = TIME_LIMIT;
		return true;
	}
	if (this->subProbSolved >= this->lpLimit) {
		if (LOG_ND) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Number of LPs solved reached Limit. Subproblems solved: " + utils::ToolBox::convertToString(this->subProbSolved));
		}
		this->solutionStatus = LP_LIMIT;
		return true;
	}
	return false;
}

bool NbdMaster::isRedundant(std::vector<BendersCut>& vec, BendersCut& c) {
	if (c.cType == BendersCut::EMPTY)
		return true;
	for (unsigned int i = 0; i < vec.size(); i++) {
		if (vec[i].equalsLHS(c)) {
			if (vec[i].rhs.getRatioSign() == c.rhs.getRatioSign()) {
				if (vec[i].rhs.getValue() == c.rhs.getValue())
					return true;
				if (vec[i].rhs.getRatioSign() == data::QpRhs::greaterThanOrEqual) {
					if (vec[i].rhs.getValue() < c.rhs.getValue()) {
						vec[i].rhs.setValue(c.rhs.getValue());
					}
					return true;
				} else if (vec[i].rhs.getRatioSign() == data::QpRhs::smallerThanOrEqual) {
					if (vec[i].rhs.getValue() > c.rhs.getValue()) {
						vec[i].rhs.setValue(c.rhs.getValue());
					}
					return true;
				} else {
					utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "vec[i]: " + vec[i].toString());
					utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "cut   : " + c.toString());
					throw utils::AlgorithmException("NbdMaster::isRedundant --> ==");
				}
			}
		}
	}
	return false;
}

void NbdMaster::updateScenarioRhs(const NbdTreeNode& node) {

	if (qlpSplitter.conIndAtDepth[node.depth].second == -1)
		return;

	unsigned int startCons = qlpSplitter.conIndAtDepth[node.depth].first;
	unsigned int endCons = qlpSplitter.conIndAtDepth[node.depth].second;

	std::vector<std::vector<data::IndexedElement> > &m = (this->sampleAverageApproximation) ? this->qlpSplitter.fastRandMatrix : this->qlpSplitter.fastUnivMatrix;
	data::QpNum tmp;
	if (!node.data) {
	for (unsigned int i = startCons, size = endCons + 1; i < size; i++) {
			std::vector<data::IndexedElement>& row = m[i];
			currEventRhs[i].setZero();
			for (unsigned int j = 0, elems = row.size(); j < elems; j++) {
				(tmp = row[j].value) *= this->qlpSplitter.currScenVarString[row[j].index];
				currEventRhs[i] += tmp;
			}
		}
	}

	for (unsigned int j = 0, i = startCons, size = endCons + 1; i < size; i++, j++) {
		//((currScenRhs[i] = qlpSplitter.originalRhs[i]) -= currPropRhs[i]);
		//currScenRhs[i] -= node.data ? node.data->scenarioRhs[j] : currEventRhs[i];
		((tmp=qlpSplitter.originalRhs[i]) -= currPropRhs[i])-= node.data ? node.data->scenarioRhs[j] : currEventRhs[i];
		currScenRhs[i]=tmp;
	}
}

void NbdMaster::updateProposalRhs(const NbdTreeNode& node) {
	if (qlpSplitter.conIndAtDepth[node.depth].second == -1)
		return;
	for (unsigned int i = qlpSplitter.conIndAtDepth[node.depth].first, size = qlpSplitter.conIndAtDepth[node.depth].second + 1; i < size; i++) {
		currPropRhs[i].setZero();
		std::vector<data::IndexedElement>& row = this->qlpSplitter.fastExistMatrix[i];
		for (unsigned j = 0, elems = row.size(); j < elems; j++) {
			currPropRhs[i] += (row[j].value * this->propVec[row[j].index]);
		}
	}
}

void NbdMaster::updateProposalVector(unsigned int depth, const std::vector<data::QpNum>& proposal) {
	unsigned int offset = this->qlpSplitter.proposalIndexAtDepth[depth];
	for (unsigned int i = 0, size = this->qlpSplitter.proposalsAtDepth[depth]; i < size; i++) {
		this->propVec[offset + i] = proposal[i];
	}
}

void NbdMaster::updateProposalVector(unsigned int depth, const BendersSolution& sol) {
	unsigned int offset = this->qlpSplitter.proposalIndexAtDepth[depth];
	for (unsigned int i = 0, size = this->qlpSplitter.proposalsAtDepth[depth]; i < size; i++) {
		this->propVec[offset + i] = sol.varAlloc[i];
	}
}

unsigned int NbdMaster::addUserCuts(const std::vector<data::IndexedElement>& lhs, const data::QpRhs::RatioSign& sign, const data::QpNum& rhs) {

	int depth = -1;
	unsigned int lastIndex = lhs[lhs.size() - 1].index;
	for (unsigned int i = 0; i < this->qlpSplitter.existVarIndAtDepth.size(); i++) {
		if (lastIndex <= this->qlpSplitter.existVarIndAtDepth[i].second) {
			depth = i;
			break;
		}
	}
	if (depth < 0)
		throw utils::AlgorithmException("NbdMaster::addUserCuts(...) -> could not determine cut depth for " + data::indexedElementVecToString(lhs));
	if (!depth) {
		//utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Adding UserCut at Root Node");
	} else if (depth == this->qpTree.size() - 1) {
		//utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Adding UserCut at Leaf Node");
	} else {
		throw utils::AlgorithmException("NbdMaster::addUserCuts(...) -> Adding inner UserCut not yet implemented. Depth: " + utils::ToolBox::convertToString(depth));
	}

	BendersCut cut(BendersCut::USER, sign, rhs);

	unsigned int breakIndex = qlpSplitter.existVarIndAtDepth[depth].first;

	for (unsigned int i = 0; i < lhs.size(); i++) {
		if (lhs[i].index < breakIndex) {
			cut.recourseRow.push_back(lhs[i]);
		} else {
			cut.masterRow.push_back(lhs[i]);
		}
	}
	this->masters[depth]->addStageCut(userCutId, cut);
	this->userCutIdToStageMap.insert(std::make_pair(userCutId, depth));
	userCuts++;
	return userCutId++;

}

void NbdMaster::removeUserCut(unsigned int ucid) {
	unsigned int stage = this->userCutIdToStageMap.find(ucid)->second;
	this->masters[stage]->removeStageCut(ucid);
	this->userCutIdToStageMap.erase(ucid);
	userCuts--;
}


void NbdMaster::removeAllConstraints() {
	unsigned int stage = 0;//this->userCutIdToStageMap.find(ucid)->second;
	this->masters[stage]->removeStageAllConstraints();
	this->userCutIdToStageMap.clear();
	currEventRhs.clear();
	/** The current rhs after proposal rhs and scenario rhs have been substracted from the original rhs*/
	currScenRhs.clear();
	/** This vector contains the values that are substituted and subtracted form rhs due to the current scenario (univ matrix) */
	currPropRhs.clear();
	userCuts = 0;
}

void NbdMaster::removeUserCutsFromCut(unsigned int ucid) {
	unsigned int stage = this->userCutIdToStageMap.find(ucid)->second;
	this->masters[stage]->removeStageCutsFromCut(ucid);
	std::map<unsigned int, unsigned int>::iterator it, end;
	it = this->userCutIdToStageMap.find(ucid);
	end = this->userCutIdToStageMap.end();
	unsigned int size = 0;
	while (it != end) {
		size++;
		it++;
	}
	this->userCutIdToStageMap.erase(it, end);
	userCuts -= size;
}

void NbdMaster::changeUserCutRhs(unsigned int ucid, const data::QpNum& newRhs) {
	unsigned int stage = this->userCutIdToStageMap.find(ucid)->second;
	this->masters[stage]->changeStageCutRhs(ucid, newRhs);
}

void NbdMaster::changeRhsVal(unsigned int cIndex, const data::QpNum& val) {
	this->qlpSplitter.originalRhs[cIndex] = val;
	for (unsigned int i = 0; i < this->masters.size(); i++)
		this->masters[i]->reloadMasterRhs(cIndex);
}

void NbdMaster::removeUserCuts() {
	for (unsigned int i = 0; i < this->masters.size(); i++)
		this->masters[i]->removeStageCuts();
	this->userCutIdToStageMap.clear();
	userCuts = 0;
}

void NbdMaster::setVariableBounds(unsigned int i, const data::QpNum& lb, const data::QpNum& ub) {
	this->setVariableLB(i, lb);
	this->setVariableUB(i, ub);
}

void NbdMaster::setVariableLB(unsigned int i, const data::QpNum& lb) {
	data::QpVar& v = this->qlpWork.getVariableByIndex(i);
	if (v.getQuantifier() != data::QpVar::exists) {
		return;
		throw utils::AlgorithmException("NbdMaster::setVariableLB(unsigned int i, const data::QpNum& lb) --> NOT existential variable");
	}
	v.setLowerBound(lb);
	this->masters[getExistVariableStage(i)]->setVariableLB(getExistVariableStageIndex(i), lb);
}

void NbdMaster::setVariableUB(unsigned int i, const data::QpNum& ub) {
	data::QpVar& v = this->qlpWork.getVariableByIndex(i);
	if (v.getQuantifier() != data::QpVar::exists) {
		return;
		throw utils::AlgorithmException("NbdMaster::setVariableLB(unsigned int i, const data::QpNum& lb) --> NOT existential variable");
	}
	v.setLowerBound(ub);
	this->masters[getExistVariableStage(i)]->setVariableUB(getExistVariableStageIndex(i), ub);
}

int NbdMaster::getExistVariableStage(unsigned int index) const {
	int stage = -1;
	for (unsigned int i = 0; i < this->qlpSplitter.existVarIndAtDepth.size(); i++) {
		if (this->qlpSplitter.existVarIndAtDepth[i].second >= index) {
			stage = i;
			break;
		}
	}
	return stage;
}

int NbdMaster::getExistVariableStageIndex(unsigned int index) const {
	int existIndex = -1;
	for (unsigned int i = 0; i < this->qlpSplitter.existVarIndAtDepth.size(); i++) {
		if (this->qlpSplitter.existVarIndAtDepth[i].second >= index) {
			(existIndex = index) -= this->qlpSplitter.existVarIndAtDepth[i].first;
			break;
		}
	}
	return existIndex;
}

void NbdMaster::clearCuts() {
	if (this->mode == SINGLE_STAGE)
		return;
	//removes all BendersCuts, clears the bases and restores the original move order
	for (unsigned int i = 0; i < this->qpTree.size(); i++) {
		for (unsigned int j = 0; j < this->qpTree[i].size(); j++) {
			this->qpTree[i][j].clear();
		}
	}
}

void NbdMaster::restart(bool clearStageSolvers) {
	this->totalIterations = this->subProbSolved = this->relaxProbSolved = 0;
	this->alphaCuts = this->betaCuts = this->ipcuts = this->deepAlphaCuts = this->deepBetaCuts = 0;
	this->solutionStatus = UNKNOWN;
	for (unsigned int i = 0; i < currEventRhs.size(); i++) {
		currEventRhs[i].setZero();
		currScenRhs[i].setZero();
		currPropRhs[i].setZero();
	}
	for (unsigned int i = 0; i < propVec.size(); i++) {
		propVec[i].setZero();
	}

//prepare tree nodes for restart
	for (unsigned int i = 0; i < this->qpTree.size(); i++) {
		for (unsigned int j = 0; j < this->qpTree[i].size(); j++) {
			this->qpTree[i][j].solved = 0;
		}
	}

//prepare stage solvers for restart
	if (clearStageSolvers) {
		for (unsigned int i = 0; i < masters.size(); i++) {
			masters[i]->clearStageSolver();
		}
	}

	this->timer.restart();
}

const std::vector<BendersCut>& NbdMaster::getFirstStageBendersCuts() const {
	return *(this->qpTree[0][0].scenarioCuts);
}

const std::vector<BendersCut>& NbdMaster::getFirstStageUserCuts() const {
	return this->masters[0]->uCutVec;
}

void NbdMaster::setAdvancedStartInformation(const data::QpNum& lb, const data::QpNum& ub, const std::vector<data::QpNum>& varAlloc) {
	this->lbVal = lb;
	this->ubVal = ub;
	this->varAlloc = varAlloc;
}

void NbdMaster::printSolutionInformation() const {

	std::string st("Finished  : " + utils::ToolBox::convertToString(timer.getSeconds(), 2));
	st += ", ";
	st += extSol::QpExternSolver::solutionStatusToString(this->nbdSol.status);
	st += ", " + algorithm::Algorithm::solutionCaseToString(this->solutionCase);
	st += ", ";
	(this->protocol == FAST_FORWARD) ? st += "FF" : ((this->protocol == FAST_BACK) ? st += "FB" : (this->protocol == FAST_FORWARD_FAST_BACK) ? st += "FF_FB" : st += "FF_FFB");
	st += ", ";
	(this->mode == SINGLE_STAGE) ? st += "SingleStage" : (this->mode == TWO_STAGE) ? st += "TwoStage" : st += "MultiStage";
	st += ", Scenarios:  " + utils::ToolBox::convertToString(this->qpTree[qpTree.size() - 1].size());
	st += ", Iterations: " + utils::ToolBox::convertToString(totalIterations);
	st += ", LPs Solved: " + utils::ToolBox::convertToString(subProbSolved);

	if (this->mode == MULTI_STAGE) {
		if (MS_ALPHA)
			st += ", Alpha Cuts: " + utils::ToolBox::convertToString(this->alphaCuts);
		if (MS_BETA)
			st += ", Beta Cuts:  " + utils::ToolBox::convertToString(this->betaCuts);
		if (MS_COMP_REL_UB) {
			st += ", Deep-Alpha Cuts: " + utils::ToolBox::convertToString(deepAlphaCuts);
			st += ", Deep-Beta Cuts: " + utils::ToolBox::convertToString(deepBetaCuts);
			st += ", Relaxations solved: " + utils::ToolBox::convertToString(relaxProbSolved);
		}
	}
	if (this->qlpType == data::Qlp::QIP_BIN)
		st += ", IP-Cuts:  " + utils::ToolBox::convertToString(ipcuts);

	st += ", Diff: " + utils::ToolBox::convertToString((this->boundVec[0].second.asDouble() - this->boundVec[0].first.asDouble()), 5);
	st += ", Lb: " + this->boundVec[0].first.toString();
	st += ", Ub: " + this->boundVec[0].second.toString();
	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, st);
}

////---------------------------------------------------- BENDERS DECOMPOSITION FOR R-QLPs --------------------------------------------------------------------------------------------------->
BendersSolution NbdMaster::TwoStageRandom(NbdTreeNode & node) {

//-------------------------- some preparations ------------------------------------------------------------->
	BendersSolution solution;
	extSol::QpExternSolver::QpExtSolSolutionStatus msStat, ssStat;
	std::vector<data::QpNum> proposal, objVec;
	data::QpNum msObj, z, ssObj, cutRhsAvgVal, ub;
	BendersCut cut;

	bool infeasibleSubProblems = false, newCuts = false;
	int moveIndex = -1;
	std::vector<BendersCut> optCuts, feasCuts, infeasOptCuts;

	this->boundVec[node.depth].first.setMinInf();
	this->boundVec[node.depth].second.setMaxInf();

//-------------------------- initially solve master problem ----------------------------------------------->
	this->masters[node.depth]->loadNode(node, this->propVec);
	if (this->masters[node.depth]->solveMaster(node)) {
		if ((msStat = this->masters[node.depth]->getSolutionStatus()) == extSol::QpExternSolver::INFEASIBLE) {
			solution.status = extSol::QpExternSolver::INFEASIBLE;
			solution.objFuncVal.setZero();
			return solution;
		} else {
			msObj = this->masters[node.depth]->getObjValue();
			this->masters[node.depth]->getValues(proposal);
		}
	} else {
		throw utils::AlgorithmException("Exception solving initial master problem");
	}

	unsigned int sApprInd = proposal.size() - 1;
	std::vector<data::QpNum> aggrCut(proposal.size(), 0);

	if (this->mode == TWO_STAGE) {

		//--------------- start loop and run until LB==UB || no new cuts || master becomes infeasible ------------->
		do {
			//----------------------------------------- Next iteration -------------------------------------------->
			if (LOG_ND)
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Iteration: " + utils::ToolBox::convertToString(totalIterations));

			this->totalIterations++;

			this->boundVec[node.depth].first = msObj;
			(z = msObj) -= proposal[sApprInd];

			if (checkUserBreak())
				break;

			this->updateProposalVector(node.depth, proposal);

			this->updateProposalRhs(*node.childVec->operator[](0));

			//------------------------------------------ Solve all Subproblems ------------------------------------>
			if (LOG_ND)
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving subproblems ...");

			infeasibleSubProblems = newCuts = false;
			feasCuts.clear();
			optCuts.clear();
			objVec.clear();

			for (unsigned int i = 0, successors = node.moveOrder->size(); i < successors; i++) {

				if (checkUserBreak())
					break;

				if ((moveIndex = node.moveOrder->operator [](i)) < 0)
					continue;

				if (qlpSplitter.univObjVec.size() || qlpSplitter.randObjVec.size()) {
					this->masters[node.depth + 1]->changeObjectiveOffSet(node.childVec->operator[](moveIndex)->data->scenarioOffSet);
				}

				this->updateScenarioRhs(*node.childVec->operator[](moveIndex));

				masters[node.depth + 1]->loadNode(*node.childVec->operator[](moveIndex), this->propVec);

				this->masters[node.depth + 1]->changeMasterRhs();

				if (LOG_ND) {
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving subproblem: " + utils::ToolBox::convertToString(i));
				}

				if (!this->masters[node.depth + 1]->solveMaster(*(node.childVec->operator[](moveIndex)))) {
					throw utils::AlgorithmException("Error Solving Subproblem: " + utils::ToolBox::convertToString(moveIndex));
				}
				if ((ssStat = this->masters[node.depth + 1]->getSolutionStatus()) != extSol::QpExternSolver::INFEASIBLE) {
					if (infeasibleSubProblems)
						continue;
					this->masters[node.depth + 1]->getOptCut(*(node.childVec->operator[](moveIndex)), cut);
					objVec.push_back(this->masters[node.depth + 1]->getObjValue());
					optCuts.push_back(cut);
				} else {
					infeasibleSubProblems = true;
					this->masters[node.depth + 1]->getFeasCut(*(node.childVec->operator[](moveIndex)), cut);
					if (!this->isRedundant(feasCuts, cut) && !node.isRedundant(cut)) {
						newCuts = true;
						feasCuts.push_back(cut);
						if (TS_STOP_NON_RED_INF_CUT) {
							break;
						}
					} else {
						if (LOG_ND) {
							utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Feasibility Cut redundant: " + cut.toString());
						}
					}
				}
			}

			//--------------------------- Adding cuts --------------------------------------------------------------->
			if (LOG_ND) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Adding Cuts...");
			}

			if (infeasibleSubProblems) {
				for (unsigned int i = 0; i < feasCuts.size(); i++) {
					this->masters[node.depth]->addBendersCut(node, feasCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS);
				}
			} else {

				ub.setZero();
				cutRhsAvgVal.setZero();
				std::fill(aggrCut.begin(), aggrCut.end(), 0);
				double p = 0;

				for (unsigned int i = 0; i < optCuts.size(); i++) {
					p = node.childVec->operator[](i)->data->probability.asDouble(); //TODO might be a problem
					ub += (objVec[i] * p);
					cutRhsAvgVal += (optCuts[i].rhs.getValue() * p);
					for (unsigned int j = 0; j < optCuts[i].masterRow.size(); j++) {
						aggrCut[optCuts[i].masterRow[j].index] += (optCuts[i].masterRow[j].value * p);
					}
				}

				BendersCut c(BendersCut::OPTIMALITY, data::QpRhs::greaterThanOrEqual, cutRhsAvgVal);
				for (unsigned int i = 0; i < aggrCut.size(); i++) {
					if (!aggrCut[i].isZero())
						c.masterRow.push_back(data::IndexedElement(i, aggrCut[i]));
				}

				if (this->masters[node.depth]->addBendersCut(node, c, this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS)) {
					newCuts = true;
				}

			}
			//------------- Check if we can stop due to matching bounds ------------------------------------------->
			if (infeasibleSubProblems) {
				this->boundVec[node.depth].second.setMaxInf();
			} else {
				(this->boundVec[node.depth].second = ub) += z;
				if (this->solutionCase == FEASIBILITY)
					break;
				if (TS_CHECK_BOUNDS && (this->boundVec[node.depth].second - this->boundVec[node.depth].first) < TS_BREAK_EPSILON)
					break;
			}

			//------------- If there are no new cuts we can stop -------------------------------------------------->
			if (!newCuts) {
				if (LOG_ND)
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "No New Cuts");
				if (TS_BREAK_NO_NEW_CUTS) {
					break;
				}
			}

			//------------- Check for break due to timeout or max number of LPs to solve is reached --------------->
			if (checkUserBreak())
				break;
			//------------- Resolve the master to obtain new candidate solution ----------------------------------->
			if (LOG_ND)
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Resolving Master...");
			if (!this->masters[node.depth]->solveMaster(node)) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Error solving Master LP.");
			} else if ((msStat = this->masters[node.depth]->getSolutionStatus()) == extSol::QpExternSolver::INFEASIBLE) {
				if (LOG_ND) {
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Break: Master infeasible");
				}
				break;
			} else {
				msObj = this->masters[node.depth]->getObjValue();
				this->masters[node.depth]->getValues(proposal);
			}
			//----------------------------------------------------------------------------------------------------->
		} while (true);
	} else {

	}

//Create Solution
	if (msStat != extSol::QpExternSolver::INFEASIBLE) {
		if (infeasibleSubProblems && this->subProbSolved < this->lpLimit) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, +"MasterSolution feasible but infeasible subproblems at depth 0");
			solution.status = extSol::QpExternSolver::INFEASIBLE;
			solution.objFuncVal.setZero();
			solution.varAlloc.clear();
		} else {
			solution.status = extSol::QpExternSolver::OPTIMAL;
			solution.objFuncVal = this->masters[node.depth]->getObjValue();
			this->masters[node.depth]->getValues(solution.varAlloc);
		}
	} else {
		solution.status = extSol::QpExternSolver::INFEASIBLE;
		solution.objFuncVal.setZero();
		solution.varAlloc.clear();
	}
	return solution;
}

BendersSolution NbdMaster::NestedBendersAvg(NbdTreeNode& node, const data::QpNum& parentA, const data::QpNum& parentB) {

	this->totalIterations++;

//The solution of the master at this stage
	BendersSolution masterSol;

//Get the CplexWrapper that is responsible for this stage
	NbdStageSolver* master = this->masters[node.depth];

//Offset from current scenario due to univ vars in the objective function
	data::QpNum scenarioOffSet;

//Load the current node int the NbdExtStageSolCplex
	if (node.depth) {

		//Load the current node to the NbdExtStageSolCplex, if we are not at the leaf
		if (node.childVec != NULL)
			master->loadNode(node, this->propVec);

		//Check if the offst must be changed (can only happen if univ var coeff are in the objective function)
		if (qlpSplitter.univObjVec.size() || qlpSplitter.randObjVec.size()) {
			master->changeObjectiveOffSet(node.data->scenarioOffSet);
		}

		//Check if this master holds a submatrix of the original input
		if (this->qlpSplitter.conIndAtDepth[node.depth].second)
			master->changeMasterRhs();
	}

//Solve the Master
	if (this->masters[node.depth]->solveMaster(node)) {
		if ((masterSol.status = this->masters[node.depth]->getSolutionStatus()) == extSol::QpExternSolver::INFEASIBLE) {
			//Infeasible
			masterSol.objFuncVal.setZero();
			if (node.depth)
				this->masters[node.depth]->getFeasCut(node, masterSol.feas);
			return masterSol;
		} else {
			//Feasible
			masterSol.objFuncVal = this->masters[node.depth]->getObjValue();
			this->masters[node.depth]->getValues(masterSol.varAlloc);
			if (node.depth)
				this->masters[node.depth]->getOptCut(node, masterSol.opt);
		}
	} else {
		throw utils::AlgorithmException("Exception solving initial master problem");
	}

	if (LOG_ND) {
		std::string s("");
		s += "\tNode: " + utils::ToolBox::convertToString(node.nodeNumber);
		s += "\tDepth: " + utils::ToolBox::convertToString(node.depth);
		s += "\nMaster Solution: " + masterSol.toString();
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, s);
	}

	if (node.childVec != NULL) {

		//Check the solution of the master
		if (masterSol.status != extSol::QpExternSolver::INFEASIBLE) {

			this->boundVec[node.depth].first.setMinInf();
			this->boundVec[node.depth].second.setMaxInf();

			BendersSolution subSol;

			std::vector<BendersCut> feas, opt;
			std::vector<data::QpNum> objFuncVals;

			data::QpNum masterOffsetValue, ubOffset, avgSubSol, cutRhsAvgVal;

			bool infeasibleSubProblems = false, newCuts = false;

			do {
				//------------------------------------ Solve all Subproblems ----------------------------------------->

				this->boundVec[node.depth].first = masterSol.objFuncVal;

				this->updateProposalVector(node.depth, masterSol);

				this->updateProposalRhs(*node.childVec->operator[](0));

				newCuts = infeasibleSubProblems = false;

				int moveIndex = -1;

				feas.clear(), opt.clear();
				objFuncVals.clear();

				for (unsigned int i = 0, successors = node.moveOrder->size(); i < successors; i++) {

					if ((moveIndex = node.moveOrder->operator [](i)) < 0) {
						continue;
					}

					this->updateScenarioRhs(*node.childVec->operator[](moveIndex));

					subSol = NestedBendersAvg(*node.childVec->operator[](moveIndex), 0, 0);

					if (subSol.status != extSol::QpExternSolver::INFEASIBLE) {
						opt.push_back(subSol.opt);
						objFuncVals.push_back(subSol.objFuncVal);
					} else {
						infeasibleSubProblems = true;
						if (this->masters[node.depth]->addBendersCut(node, subSol.feas, this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS)) {
							newCuts = true;
							//if (MULTISTAGE_STOP_INF_SUB_PROB) {
							//	break;
							//}
						}
					}
				}
				//---------------------------------------------------------------------------------------------------->

				//------------------------------------ Create weighted cut   ----------------------------------------->
				if (!infeasibleSubProblems) {

					if (this->solutionCase == FEASIBILITY)
						break;

					avgSubSol.setZero();
					cutRhsAvgVal.setZero();

					std::vector<data::QpNum> aggrCutRecourse;

					if (node.depth && qlpSplitter.existVarIndAtDepth[node.depth].first > 0)
						aggrCutRecourse = std::vector<data::QpNum>(qlpSplitter.existMatrixIndexAtDepth[node.depth].first, 0);

					std::vector<data::QpNum> aggrCutMaster(2 + qlpSplitter.existVarIndAtDepth[node.depth].second - qlpSplitter.existVarIndAtDepth[node.depth].first, 0);

					double p = 0;
					for (unsigned int i = 0; i < opt.size(); i++) {

						p = node.childVec->operator[](i)->data->probability.asDouble();

						avgSubSol += (objFuncVals[i] * p);
						cutRhsAvgVal += (opt[i].rhs.getValue() * p);

						for (unsigned int j = 0; j < opt[i].recourseRow.size(); j++) {
							aggrCutRecourse[opt[i].recourseRow[j].index] += (opt[i].recourseRow[j].value * p);
						}
						for (unsigned int j = 0; j < opt[i].masterRow.size(); j++) {
							aggrCutMaster[opt[i].masterRow[j].index - qlpSplitter.existMatrixIndexAtDepth[node.depth].first] += (opt[i].masterRow[j].value * p);
						}
					}

					BendersCut c(BendersCut::OPTIMALITY, opt[0].rhs.getRatioSign(), cutRhsAvgVal);
					for (unsigned int i = 0; i < aggrCutRecourse.size(); i++) {
						if (!aggrCutRecourse[i].isZero()) {
							c.recourseRow.push_back(data::IndexedElement(i, aggrCutRecourse[i]));
						}
					}
					for (unsigned int i = 0; i < aggrCutMaster.size(); i++) {
						if (!aggrCutMaster[i].isZero()) {
							c.masterRow.push_back(data::IndexedElement(qlpSplitter.existMatrixIndexAtDepth[node.depth].first + i, aggrCutMaster[i]));
						}
					}

					if (this->masters[node.depth]->addBendersCut(node, c, this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS)) {
						newCuts = true;
					}
				}
				//---------------------------------------------------------------------------------------------------->

				//------------------------------------ Resolve updated master ---------------------------------------->
				if (!infeasibleSubProblems) {

					if (this->solutionCase == FEASIBILITY)
						break;
					data::QpNum tmp(masterSol.objFuncVal);
					tmp -= masterSol.varAlloc[masterSol.varAlloc.size() - 1];
					tmp += avgSubSol;
					this->boundVec[node.depth].second = tmp;
					if (MS_CHECK_BOUNDS && (this->boundVec[node.depth].second - this->boundVec[node.depth].first) < MS_BREAK_EPSILON)
						break;
				}

				if (!newCuts)
					break;

				if (!this->masters[node.depth]->solveMaster(node)) {
					utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Error solving LP at depth: " + utils::ToolBox::convertToString(node.depth));
				} else if ((masterSol.status = this->masters[node.depth]->getSolutionStatus()) == extSol::QpExternSolver::INFEASIBLE) {
					masterSol.clear();
					masterSol.status = extSol::QpExternSolver::INFEASIBLE;
					if (node.depth)
						this->masters[node.depth]->getFeasCut(node, masterSol.feas);
					break;
				} else {
					masterSol.objFuncVal = this->masters[node.depth]->getObjValue();
					this->masters[node.depth]->getValues(masterSol.varAlloc);
					if (node.depth)
						this->masters[node.depth]->getOptCut(node, masterSol.opt);
				}

				if (this->protocol == FAST_FORWARD) {
					//nothing to do
				} else if (this->protocol == FAST_FORWARD_FAST_BACK) {
					if (infeasibleSubProblems)
						continue;
					if (node.depth && !node.parent->isRedundant(masterSol.opt))
						break;
				} else {
					throw utils::AlgorithmException("Protocol not yet implemented...");
				}
				//---------------------------------------------------------------------------------------------------->
			} while (true);
			//revert all changes that have been made by this node, but keep level 0 proposals since we are interested in them
			if (node.depth) {
				for (unsigned int i = 0; i < this->qlpSplitter.proposalsAtDepth[node.depth]; i++) {
					this->propVec[this->qlpSplitter.proposalIndexAtDepth[node.depth] + i].setZero();
				}
			}
		} else {
			if (LOG_ND)
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, +"The Master is initially infeasible. Returning Cut.");
		}
	} else {
		if (LOG_ND)
			utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, +"\nLeaf Node, nothing more to do.");
	}
	if (LOG_ND) {
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, +"Passing back solution: " + masterSol.toString());
	}
	return masterSol;
}

//----------------------------------------------------------- Nested Benders Decomposition algorithm ------------------------------------------------>
bool NbdMaster::MultiStageTmp(BendersSolution& msSol, NbdTreeNode& node) {

	//Start the next iteration
	this->totalIterations++;

	//The solution of the master at this stage
	msSol.clear();

	//Get the CplexWrapper that is responsible for this stage
	NbdStageSolver* master = this->masters[node.depth];

	//Reinitialize bound vector
	this->boundVec[node.depth].first.setMinInf();
	this->boundVec[node.depth].second.setMaxInf();

	//Load the current node to the NbdExtStageSolCplex, if we are not at the root and not at a leaf
	master->loadNode(node, this->propVec);

	//Load the current node int the NbdExtStageSolCplex
	if (node.depth) {
		//Check if the offset must be changed (can only happen if univ var coeff are in the objective function)
		if (qlpSplitter.univObjVec.size() || qlpSplitter.randObjVec.size()) {
			master->changeObjectiveOffSet(this->qlpSplitter.precomputeSzenarioOffSet(node.parent->depth, qlpSplitter.currScenVarString));
		}
		//Check if this master holds a submatrix of the original input
		if (this->qlpSplitter.conIndAtDepth[node.depth].second) {
			master->changeMasterRhs();
		}
	}

	//Solve this nodes master problem
	if (!master->solveMaster(node) && node.depth) {
		utils::Logger::globalLog(utils::LOG_ERROR, LOG_TAG, +"Exception Solving Master at Node: " + node.toString());
		msSol.status = extSol::QpExternSolver::INFEASIBLE;
		msSol.objFuncVal.setZero();
		return false;
	}

	//Check whether the master solution is initially feasible, if not directly go back
	if ((msSol.status = master->getSolutionStatus()) != extSol::QpExternSolver::INFEASIBLE) {

		//If this is a leaf node we are done and pass back the solution
		if (node.childVec != NULL && node.depth < this->currentDepth) {

			//Get master objective function and proposal
			msSol.objFuncVal = master->getObjValue();
			master->getValues(msSol.varAlloc);

			//Some datastructures we need
			BendersSolution wcSubSol, subSol;
			data::QpNum msObj, ssObj, z, phi, alpha, beta, offset;
			std::vector<BendersCut> optCuts, feasCuts, infeasOptCuts;
			bool infeasibleSubProblems = false, newCut = false;
			unsigned int subApproxIndex = msSol.varAlloc.size() - 1;

			do {

				if (this->checkUserBreak())
					break;

				//get master objective function value and split it into in master and subproblem part
				msObj = msSol.objFuncVal;
				phi = msSol.varAlloc[subApproxIndex];
				z = msObj - phi;

				this->boundVec[node.depth].first = msObj;

				this->updateProposalVector(node.depth, msSol.varAlloc);

				this->updateProposalRhs(*node.childVec->operator[](0));

				infeasibleSubProblems = newCut = false;

				int index = -1;

				wcSubSol.objFuncVal.setMinInf();

				//------------- Solve all subproblems ------------------------------------------------------------->
				for (unsigned int i = 0, successors = node.moveOrder->size(); i < successors; i++) {

					index = i;

					this->qlpSplitter.setScenarioBitVector(*node.childVec->operator[](index), index);
					this->qlpSplitter.setScenarioVariableVectorByBitVector(*node.childVec->operator[](index));

					this->updateScenarioRhs(*node.childVec->operator[](index));

					//Get SubSolution
					if (!MultiStage(subSol, *node.childVec->operator[](index), alpha, beta, newCut)) {
						continue;
						throw utils::AlgorithmException("NestedBendersException");
					}

					if (subSol.status != extSol::QpExternSolver::INFEASIBLE) {

						if (wcSubSol.objFuncVal < subSol.objFuncVal) {

							wcSubSol = subSol;

						}
						optCuts.push_back(subSol.opt);
					} else {
						infeasibleSubProblems = true;
						feasCuts.push_back(subSol.feas);
						if (_MULTISTAGE_ADD_INF_OPT_CUT)
							infeasOptCuts.push_back(subSol.opt);
					}
				}

				//------------- Check if we can stop due to matching bounds or AlphaBeta Break -------------------->
				if (infeasibleSubProblems) {
					//this->boundVec[node.depth].second.setMaxInf();
				} else {

					if (this->solutionCase == FEASIBILITY) {
						break;
					}

					(this->boundVec[node.depth].second = z) += wcSubSol.objFuncVal;

					if (MS_CHECK_BOUNDS) {

						if ((this->boundVec[node.depth].second < this->boundVec[node.depth].first)) {
							this->boundVec[node.depth].second = this->boundVec[node.depth].first;
						}

						if (((this->boundVec[node.depth].second - this->boundVec[node.depth].first) < MS_BREAK_EPSILON)) {
							break;
						}
					}
				}

//				//------------- If there are no new cuts we can stop -------------------------------------------------->
//				if (MULTISTAGE_BREAK_NO_NEW_CUTS && !newCut) {
//					if (LOG_ND)
//						utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Break: No New Cuts");
//					break;
//				}

				//--------------------------- Adding optimality cuts ----------------------------------------------------->
				for (unsigned int i = 0; i < feasCuts.size(); i++) {
					this->masters[node.depth]->addBendersCut(node, feasCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS);
				}
				for (unsigned int i = 0; i < infeasOptCuts.size(); i++) {
					this->masters[node.depth]->addBendersCut(node, infeasOptCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS);
				}

				for (unsigned int i = 0; i < optCuts.size(); i++) {
					this->masters[node.depth]->addBendersCut(node, optCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS);
				}

				feasCuts.clear();
				infeasOptCuts.clear();
				optCuts.clear();

				msSol.clear();
				subSol.clear();
				wcSubSol.clear();

				//------------- Resolve the master to obtain new candidate solution ----------------------------------->
				if (this->masters[node.depth]->solveMaster(node)) {
					if ((msSol.status = master->getSolutionStatus()) == extSol::QpExternSolver::INFEASIBLE) {
						break;
					}
				} else {
					throw utils::AlgorithmException("Error solving LP at depth: " + utils::ToolBox::convertToString(node.depth));
				}

				//------------- Do protocol dependent checks ---------------------------------------------------------->
				if (this->protocol == FAST_FORWARD) {
					//nothing to do, just continue
				} else if (this->protocol == FAST_FORWARD_FAST_BACK) {
					if (node.depth) {
						master->getOptCut(node, msSol.opt);
						//if (nonRedParentCut || !node.parent->isRedundant(msSol.opt)) {
						break;
						//}
					}
				} else if (this->protocol == FAST_FORWARD_FAST_FEASIBLE_BACK) {
					if (!infeasibleSubProblems && node.depth) {
						master->getOptCut(node, msSol.opt);
						//if (nonRedParentCut || !node.parent->isRedundant(msSol.opt)) {
						break;
						//}
					}
				} else if (this->protocol == FAST_BACK) {
					if (node.depth) {
						master->getOptCut(node, msSol.opt);
						//if (nonRedParentCut || !node.parent->isRedundant(msSol.opt)) {
						break;
						//}
					}
				} else {
					throw utils::AlgorithmException("Protocol not yet implemented...");
				}

				msSol.objFuncVal = master->getObjValue();
				master->getValues(msSol.varAlloc);

			} while (true);

			//----------------------------------- Passing back solution (depends on sequencing protocol)----------------------------------------->
			//------------- Do protocol dependent checks ---------------------------------------------------------->
			if (this->protocol == FAST_FORWARD) {
				if (msSol.status != extSol::QpExternSolver::INFEASIBLE) {
					msSol.objFuncVal = master->getObjValue();
//					if (!this->boundVec[node.depth].second.isMaxInf()) {
					msSol.objFuncVal = this->boundVec[node.depth].second;
//					} else {
//						if (node.depth) {
//							msSol.objFuncVal = 1000000000;
//						} else {
//							msSol.objFuncVal = master->getObjValue();
//						}
//					}
					if (node.depth) {
						//if (msSol.opt.cType == BendersCut::EMPTY)
						master->getOptCut(node, msSol.opt);
					} else {
						master->getValues(msSol.varAlloc);
					}
				} else {
					if (node.depth)
						master->getFeasCut(node, msSol.feas);
				}
			} else if (this->protocol == FAST_FORWARD_FAST_BACK) {
				if (msSol.status != extSol::QpExternSolver::INFEASIBLE) {
					if (!this->boundVec[node.depth].second.isMaxInf()) {
						msSol.objFuncVal = this->boundVec[node.depth].second;
					} else {
						if (node.depth) {
							msSol.objFuncVal = this->boundVec[node.depth].second;			//1000000000;
						} else {
							msSol.objFuncVal = master->getObjValue();
						}
					}
					if (node.depth) {
						if (msSol.opt.cType == BendersCut::EMPTY)
							master->getOptCut(node, msSol.opt);
					} else {
						master->getValues(msSol.varAlloc);
					}
				} else {
					if (node.depth)
						master->getFeasCut(node, msSol.feas);
				}
			} else if (this->protocol == FAST_FORWARD_FAST_FEASIBLE_BACK) {
				if (msSol.status != extSol::QpExternSolver::INFEASIBLE) {
					if (!this->boundVec[node.depth].second.isMaxInf()) {
						msSol.objFuncVal = this->boundVec[node.depth].second;
					} else {
						if (node.depth) {
							msSol.objFuncVal = 1000000000;
						} else {
							msSol.objFuncVal = master->getObjValue();
						}
					}
					if (node.depth) {
						if (msSol.opt.cType == BendersCut::EMPTY)
							master->getOptCut(node, msSol.opt);
					} else {
						master->getValues(msSol.varAlloc);
					}
				} else {
					if (node.depth)
						master->getFeasCut(node, msSol.feas);
				}
			} else if (this->protocol == FAST_BACK) {

			} else {
				throw utils::AlgorithmException("Protocol not yet implemented...");
			}

			//revert all changes that have been made by this node
			for (unsigned int i = 0; i < this->qlpSplitter.proposalsAtDepth[node.depth]; i++)
				this->propVec[this->qlpSplitter.proposalIndexAtDepth[node.depth] + i].setZero();
			//revert the scenario variable vector
			this->qlpSplitter.revertScenarioVariableVector(node);
			//revert the scenario bit string
			this->qlpSplitter.revertScenarioBitVector(node);
		} else {
			if (LOG_ND) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, +"Leaf Node, nothing more to do. Depth: " + utils::ToolBox::convertToString(node.depth));
			}
			if (msSol.status != extSol::QpExternSolver::INFEASIBLE) {
				msSol.objFuncVal = master->getObjValue();
				if (node.depth)
					master->getOptCut(node, msSol.opt);
			} else {
				if (node.depth)
					master->getFeasCut(node, msSol.feas);
			}
		}
	} else {
		if (LOG_ND) {
			utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, +"The Master is initially infeasible. Returning Cut.");
		}
		if (node.depth)
			master->getFeasCut(node, msSol.feas);
	}
	if (LOG_ND) {
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, +"Passing back solution: " + msSol.toString());
	}
	return true;
}

}

/*bool NbdMaster::MultiStage(BendersSolution& masterSol, NbdTreeNode& node,const data::QpNum& parentAlpha, const data::QpNum& parentBeta, bool nonRedParentCut) {

 //Start the next iteration
 this->totalIterations++;

 //The solution of the master at this stage
 masterSol.clear();

 //Get the CplexWrapper that is responsible for this stage
 NbdStageSolver* master = this->masters[node.depth];

 //Reinitialize bound vector
 this->boundVec[node.depth].first.setMinInf();
 this->boundVec[node.depth].second.setMaxInf();

 //Inner Node, check for Alpha-Relaxtion-Cut
 if (COMP_REL_UB && node.depth && node.childVec != NULL && !parentAlpha.isMinInf()) {
 //Setting exist part in relaxer
 for (unsigned int i = 0, size = this->qlpSplitter.existMatrixIndexAtDepth[node.depth - 1].second; i <= size; i++) {
 relaxer->setVariableFixation(node.depth, this->qlpSplitter.existVars[i].getIndex(), this->propVec[i]);
 }
 //Setting forall part in relaxer
 for (unsigned int i = 0, size = this->qlpSplitter.scenarioVarIndexAtDepth[node.depth]; i < size; i++) {
 relaxer->setVariableFixation(node.depth, this->qlpSplitter.univVars[i].getIndex(), this->qlpSplitter.currScenVarString[i]);
 }
 //Solve relaxation and check for deep cut
 data::QpNum val;
 this->relaxProbSolved++;
 if (relaxer->solveStage(node.depth, val) && (parentAlpha > val)) {
 deepAlphaCuts++;
 this->qlpSplitter.revertScenarioVariableVector(node);
 this->qlpSplitter.revertScenarioBitVector(node);
 masterSol.status = extSol::QpExternSolver::OPTIMAL;
 masterSol.objFuncVal=val;
 return true;
 }
 }

 //Load the current node to the NbdExtStageSolCplex, if we are not at the root and not at a leaf
 master->loadNode(node, this->propVec);

 //Load the current node int the NbdExtStageSolCplex
 if (node.depth) {
 //Check if the offset must be changed (can only happen if univ var coeff are in the objective function)
 if (qlpSplitter.univObjVec.size() || qlpSplitter.randObjVec.size()) {
 master->changeObjectiveOffSet(SAVE_NODE_SCENARIO_DATA ? node.scenarioOffSet : this->qlpSplitter.computeSzenarioOffSet(*node.parent, qlpSplitter.currScenVarString));
 }
 //Check if this master holds a submatrix of the original input
 if (this->qlpSplitter.conIndAtDepth[node.depth].second) {
 master->changeMasterRhs();
 }
 }

 //Solve this nodes master problem
 if (!master->solveMaster(masterSol,node)) {
 utils::Logger::globalLog(utils::LOG_ERROR, LOG_TAG, +"Exception Solving Master at Node: " + node.toString());
 masterSol.status = extSol::QpExternSolver::INFEASIBLE;
 masterSol.objFuncVal.setZero();
 }

 //If this is a leaf node we are done and pass back the solution
 if (node.childVec != NULL) {
 //Check whether the master solution is initially feasible, if not directly go back
 if (masterSol.status != extSol::QpExternSolver::INFEASIBLE) {
 //Warm-Start
 if (node.depth) {
 if (node.parent->solved == 1 && node.solved == 1 && this->varAlloc.size()) {
 masterSol.objFuncVal = this->ubVal;
 for (unsigned int i = 0, size = masterSol.varAlloc.size(), offset = this->qlpSplitter.proposalIndexAtDepth[node.depth]; i < size; i++) {
 masterSol.varAlloc[i] = this->varAlloc[offset + i];
 }
 }
 } else {
 if (node.solved == 1 && this->varAlloc.size()) {
 masterSol.objFuncVal = this->ubVal;
 for (unsigned int i = 0, size = masterSol.varAlloc.size(), offset = this->qlpSplitter.proposalIndexAtDepth[node.depth]; i < size; i++) {
 masterSol.varAlloc[i] = this->varAlloc[offset + i];
 }
 }
 }

 //Some datastructures we need
 BendersSolution wcSubSol, subSolution;
 data::QpNum msObj, ssObj, z, phi, alpha(true), beta(false);
 std::vector<BendersCut> optCuts, feasCuts, infeasOptCuts;
 bool infeasibleSubProblems = false, newCut = false , breakIt = false;
 unsigned int subApproxIndex = masterSol.varAlloc.size() - 1;

 do {

 if (LOG_ND && !node.depth) {
 this->log(node, masterSol);
 }

 if (this->checkUserBreak())
 break;

 if (this->protocol == FAST_BACK && node.depth && (nonRedParentCut || !node.parent->isRedundant(masterSol.opt))) {
 break;
 }

 //get master objective function value and split it into in master and subproblem part
 msObj = masterSol.objFuncVal;
 phi = masterSol.varAlloc[subApproxIndex];
 z = msObj - phi;

 this->boundVec[node.depth].first = masterSol.objFuncVal;

 this->updateProposalVector(node.depth, masterSol);

 this->updateProposalRhs(*node.childVec->operator[](0));

 //-----------------AlphaBeta-Heuristic------------------------------------------------->
 alpha.setMinInf();
 beta.setMaxInf();
 if (!this->boundVec[node.depth].second.isMaxInf()) {
 (beta = this->boundVec[node.depth].second) -= z;
 }

 infeasibleSubProblems = newCut = breakIt = false;

 int index = -1;

 wcSubSol.objFuncVal.setMinInf();

 //------------- Solve all subproblems ------------------------------------------------------------->
 for (unsigned int i = 0, successors = node.moveOrder.size(); i < successors; i++) {

 if (checkUserBreak())
 break;

 if ((index = node.moveOrder[i]) < 0) {
 continue;
 }

 if (!SAVE_NODE_SCENARIO_DATA || COMP_REL_UB) {
 this->qlpSplitter.setScenarioBitVector(*node.childVec->operator[](index), index);
 this->qlpSplitter.setScenarioVariableVectorByBitVector(*node.childVec->operator[](index));
 }

 this->updateScenarioRhs(*node.childVec->operator[](index));

 if(!MultiStage(subSolution, *node.childVec->operator[](index), wcSubSol.objFuncVal, beta, newCut))
 throw utils::AlgorithmException("NestedBendersException");

 if (subSolution.status != extSol::QpExternSolver::INFEASIBLE) {

 if (ALPHA && (alpha < subSolution.objFuncVal)) {
 alpha = subSolution.objFuncVal;
 }

 if (!MULTISTAGE_ADD_ALL_CUTS && infeasibleSubProblems)
 continue;

 if ( wcSubSol.objFuncVal < subSolution.objFuncVal ) {

 wcSubSol = subSolution;

 if (subSolution.opt.cType == algorithm::BendersCut::EMPTY) {
 <<"Optimality Cut is EMPTY. SolutionStatus: "<< extSol::QpExternSolver::solutionStatusToString(subSolution.status) <<std::endl;
 continue;
 }

 if (MULTISTAGE_ADD_WC_OPT_CUT) {
 if (!this->isRedundant(optCuts, subSolution.opt) && !node.isRedundant(subSolution.opt)) {
 newCut = true;
 if (optCuts.size()) {
 optCuts[0] = subSolution.opt;
 } else {
 optCuts.push_back(subSolution.opt);
 }
 if (MULTISTAGE_MOVE_ORDERING) {
 node.pushToFront(i);
 }
 }
 }

 }

 if (!MULTISTAGE_ADD_WC_OPT_CUT) {
 if (!this->isRedundant(optCuts, subSolution.opt) && !node.isRedundant(subSolution.opt)) {
 newCut = true;
 optCuts.push_back(subSolution.opt);
 }
 }

 } else {
 infeasibleSubProblems = true;
 if (!this->isRedundant(feasCuts, subSolution.feas) && !node.isRedundant(subSolution.feas)) {
 newCut = true;
 feasCuts.push_back(subSolution.feas);
 if(MULTISTAGE_ADD_INF_OPT_CUT)
 infeasOptCuts.push_back(subSolution.opt);
 if (MULTISTAGE_MOVE_ORDERING) {
 node.pushToFront(i);
 }
 if (MULTISTAGE_STOP_INF_SUB_PROB) {
 if (node.depth) {
 if (node.parent->solved != 1 || node.solved != 1 || !this->varAlloc.size())
 break;
 } else {
 if (node.solved != 1 || !this->varAlloc.size())
 break;
 }
 }
 } else {
 if (LOG_ND) {
 utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "Redundant worst-case cut detected ...");
 }
 }
 }
 }

 //------------- Check if we can stop due to matching bounds or AlphaBeta Break -------------------->
 if (infeasibleSubProblems) {
 this->boundVec[node.depth].second.setMaxInf();
 } else {

 if (this->solutionCase == FEASIBILITY) {
 break;
 }

 (this->boundVec[node.depth].second = z) += wcSubSol.objFuncVal;

 if (MULTISTAGE_CHECK_BOUNDS) {
 if ((this->boundVec[node.depth].second - this->boundVec[node.depth].first) < BENDERS_BREAK_EPSILON) {
 break;
 }
 }

 if (BETA && node.depth && parentBeta < this->boundVec[node.depth].first) {
 //if (nonRedParentCut || !node.parent->isRedundant(masterSol.opt)) {
 this->betaCuts++;
 break;
 //}
 }

 if (ALPHA && node.depth && this->boundVec[node.depth].second < parentAlpha) {
 //if (nonRedParentCut || !node.parent->isRedundant(masterSol.opt)) {
 this->alphaCuts++;
 break;
 //}
 }
 }

 //------------- If there are no new cuts we can stop -------------------------------------------------->
 if (MULTISTAGE_BREAK_NO_NEW_CUTS && !newCut) {
 if (LOG_ND)
 utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Break: No New Cuts");
 break;
 }

 //--------------------------- Adding optimality cuts ----------------------------------------------------->
 if (infeasibleSubProblems) {
 for (unsigned int i = 0; i < feasCuts.size(); i++) {
 this->masters[node.depth]->addBendersCut(node, feasCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS);
 }
 for (unsigned int i = 0; i < infeasOptCuts.size(); i++) {
 this->masters[node.depth]->addBendersCut(node, infeasOptCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS);
 }
 }

 if (MULTISTAGE_ADD_ALL_CUTS || !infeasibleSubProblems) {
 for (unsigned int i = 0; i < optCuts.size(); i++) {
 if (this->masters[node.depth]->addBendersCut(node, optCuts[i], this->propVec, CHECK_REDUNDANT_NBD_CUTS, ADAPT_REDUNDANT_NBD_CUTS)) {
 if (MULTISTAGE_ADD_WC_OPT_CUT){
 break;
 }
 }
 }
 }

 feasCuts.clear();
 infeasOptCuts.clear();
 optCuts.clear();
 wcSubSol.clear();

 //------------- Resolve the master to obtain new candidate solution ----------------------------------->
 this->masters[node.depth]->updateModel();
 if (!this->masters[node.depth]->solveMaster(masterSol, node)) {
 utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Error solving LP at depth: " + utils::ToolBox::convertToString(node.depth));
 //TODO
 masterSol.clear();
 masterSol.objFuncVal.setMinInf();
 break;
 }

 if (masterSol.status == extSol::QpExternSolver::INFEASIBLE) {
 break;
 }

 //------------- Do protocol dependent checks ---------------------------------------------------------->
 if (this->protocol == FAST_FORWARD) {
 if (infeasibleSubProblems) {
 continue;
 }
 } else if (this->protocol == FAST_FORWARD_FAST_BACK) {
 if (node.depth && (nonRedParentCut || !node.parent->isRedundant(masterSol.opt))) {
 break;
 }
 } else if (this->protocol == FAST_FORWARD_FAST_FEASIBLE_BACK) {
 if (infeasibleSubProblems)
 continue;
 if (node.depth && (nonRedParentCut || !node.parent->isRedundant(masterSol.opt))) {
 break;
 }
 } else if (this->protocol == FAST_BACK) {
 if (node.depth && (nonRedParentCut || !node.parent->isRedundant(masterSol.opt))) {
 break;
 }
 } else {
 throw utils::AlgorithmException("Protocol not yet implemented...");
 }

 } while (true);

 //----------------------------------- Passing back solution ----------------------------------------->
 if (node.depth) {
 if (masterSol.status != extSol::QpExternSolver::INFEASIBLE) {
 if (infeasibleSubProblems) {
 masterSol.status = extSol::QpExternSolver::INFEASIBLE;
 masterSol.objFuncVal = this->boundVec[node.depth].second;
 } else {
 masterSol.objFuncVal = this->boundVec[node.depth].second;
 }
 }
 //revert all changes that have been made by this node, but keep level 0 proposals since we are interested in them
 for (unsigned int i = 0; i < this->qlpSplitter.proposalsAtDepth[node.depth]; i++) {
 this->propVec[this->qlpSplitter.proposalIndexAtDepth[node.depth] + i].setZero();
 }
 } else {
 if (masterSol.status != extSol::QpExternSolver::INFEASIBLE) {
 if (infeasibleSubProblems) {
 utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, +"MasterSolution feasible but infeasible subproblems at depth 0");
 masterSol.status = extSol::QpExternSolver::INFEASIBLE;
 masterSol.objFuncVal.setZero();
 masterSol.varAlloc.clear();
 } else {
 if (this->subProbSolved < this->maxSubProbSolve)
 masterSol.objFuncVal = this->boundVec[node.depth].second;
 }
 }
 }
 //revert the scenario bit string at all positions that were changed by this node
 if (!SAVE_NODE_SCENARIO_DATA || COMP_REL_UB) {
 this->qlpSplitter.revertScenarioVariableVector(node);
 }
 this->qlpSplitter.revertScenarioBitVector(node);
 } else {
 if (LOG_ND) {
 utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, +"The Master is initially infeasible. Returning Cut.");
 }
 }
 } else {
 if (LOG_ND) {
 utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, +"Leaf Node, nothing more to do.");
 }
 }
 if (LOG_ND) {
 utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, +"Passing back solution: " + masterSol.toString());
 }
 return true;
 }*/

