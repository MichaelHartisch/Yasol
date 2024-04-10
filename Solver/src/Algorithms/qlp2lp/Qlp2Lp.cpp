/*
*
* Solver: Qlp2Lp.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Algorithms/qlp2lp/Qlp2Lp.hpp"
#include "Utilities/QlpRelaxer.hpp"
//#include <stdio.h>
#include "ExternSolvers/QpExternSolvers.hpp"

namespace algorithm {
std::string Qlp2Lp::LOG_TAG = "Qlp2LP";
Qlp2Lp::Qlp2Lp(const data::Qlp& qp, utils::QlpConverter::DepType t) :
		Algorithm(qp, this->dep), timer(), extSol(NULL), externSolverExact(NULL), depFile(), type(t), tCreateDep(0), tLoadDep(0), tSolveDep(0) {
}

Qlp2Lp::Qlp2Lp(const data::Qlp& qp, const std::string& fName, utils::QlpConverter::DepType t) :
		Algorithm(qp, this->dep), timer(), extSol(NULL), externSolverExact(NULL), depFile(fName), type(t), tCreateDep(0), tLoadDep(0), tSolveDep(0) {
}

Qlp2Lp::~Qlp2Lp() {
	delete extSol;
}

Algorithm::QlpSolution Qlp2Lp::solveQlp(SolutionCase s) {

	if (!qlpWork.getVariableVectorByQuantifier(data::QpVar::all).size() && !qlpWork.getVariableVectorByQuantifier(data::QpVar::random).size()) {
		if (LOG_QLP2LP){ utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving LP/MIP..."); }
	} else if (qlpWork.getStageCount() == 1) {
		if (qlpWork.getVariableByIndex(0).getQuantifier() != data::QpVar::exists) {
			if (LOG_QLP2LP){ utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving A...AE...E-Problem..."); }
		} else {
			if (LOG_QLP2LP){ utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving E...EA...A-Problem..."); }
		}
	}

	timer.restart();

	Algorithm::QlpSolution sol(s, qlpWork.getObjective());

	if (this->depFile.empty()) {

		if (LOG_QLP2LP) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating DEP...");
		}

		//DEP parts
		data::QpObjFunc obj;
		std::vector<data::QpVar> vars;
		data::QpSparseMatrix matrix;
		std::vector<data::QpRhs> rhs;

		timer.restart();
		utils::QlpConverter::convertToLP(qlpWork, obj, vars, matrix, rhs, type, (s != AVERAGE_CASE) ? utils::QlpConverter::WORST : utils::QlpConverter::AVG);
		if (LOG_QLP2LP) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Loading DEP...");
		}
		timer.stop();
		tCreateDep = timer.getSeconds();
		timer.restart();

		if (LOG_QLP2LP) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Loading into inexact solver...");
		}

		this->extSol = extSol::initExternSolver();
		extSol->init(obj, vars, matrix, rhs);

		if (LOG_QLP2LP) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Loading into exact solver...");
		}

		#if defined(EXACT_ARITHMETIC) && defined(DEP_EXACT)
		externSolverExact = extSol::initExactExternSolver();
		externSolverExact->init(obj, vars, matrix, rhs);
		#endif

		timer.stop();
		tLoadDep = timer.getSeconds();

		if (LOG_QLP2LP) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Clearing Tmp Data...");
		}

	} else {

		if (LOG_QLP2LP) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Reading DEP-File: "+this->depFile);
		}

		timer.restart();

		this->extSol = extSol::initExternSolver();
		extSol->readFromFile(this->depFile);

		#if defined(EXACT_ARITHMETIC) && defined(DEP_EXACT)
		externSolverExact = extSol::initExactExternSolver();
		externSolverExact->readFromFile(this->depFile);
		#endif


		timer.stop();
		tLoadDep = timer.getSeconds();

	}

	if (LOG_QLP2LP) {
		std::string lpInfo("(Columns: ");
		lpInfo += utils::ToolBox::convertToString(extSol->getVariableCount());
		lpInfo += ", Rows: " + utils::ToolBox::convertToString(extSol->getRowCount());
		lpInfo += ")";
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving LP " + lpInfo);
	}

	timer.restart();
	if(qlpWork.getQlpType()==data::Qlp::QLP || qlpWork.getQlpType()==data::Qlp::LP){
		extSol->adaptToSolverMode(extSol::QpExternSolver::DEFAULT);
	}else{
		extSol->adaptToSolverMode(extSol::QpExternSolver::DEFAULT);
	}


	if (LOG_QLP2LP) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving with inexact solver...");
	}

	sol.solution.status = extSol->solve(2000000000, DEP_TIMER ? DEP_TIME_LIMIT : 2000000000);


	if (LOG_QLP2LP) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving with exact solver...");
	}

	//if we use exact arithmetic, verify the solution using an ecaxt solver
	if (this->externSolverExact) {
		extSol::QpExternSolver::QpExtSolBase base;
		this->extSol->getBase(base);
		this->externSolverExact->setBase(base);
		if (this->extSol->getSolutionStatus() == extSol::QpExternSolver::INFEASIBLE) {
			std::vector<data::QpNum> tmp;
			this->extSol->getDualFarkas(tmp);
			this->externSolverExact->setRayGuess(tmp);
		}
		this->externSolverExact->solve(1.0e+6, 1.0e+6);
		sol.solution.status = this->externSolverExact->getSolutionStatus();
	}


	timer.stop();
	tSolveDep = timer.getSeconds();

	if (sol.solution.status == extSol::QpExternSolver::INFEASIBLE) {
		sol.solution.ofVal.setZero();
		sol.solution.varAlloc.clear();
	} else if (sol.solution.status != extSol::QpExternSolver::INFEASIBLE) {
		sol.solution.ofVal = externSolverExact ? externSolverExact->getObjValue() : extSol->getObjValue();
		std::vector<data::QpNum> tmpVec;
		if(externSolverExact){externSolverExact->getValues(tmpVec);}else{extSol->getValues(tmpVec);}

		//If worst-case DEP, then the variable at position 0 is the objective function approximation variable and the start index is 1
		unsigned int offset = (s != AVERAGE_CASE) ? 1 : 0;

		//If input was an LP, offset must also be zero because objective function is not pushed to matrix
		if (!qlpWork.getVariableVectorByQuantifier(data::QpVar::all).size()
				&& !qlpWork.getVariableVectorByQuantifier(data::QpVar::random).size())
			offset = 0;
		for (unsigned int i = 0; i < tmpVec.size(); i++) {
			if (qlpWork.getVariableByIndex(i).getQuantifier() != data::QpVar::exists)
				break;
			if((i+offset)<tmpVec.size())
			sol.solution.varAlloc.push_back(tmpVec[i+offset]);
		}
	} else if (sol.solution.status == extSol::QpExternSolver::ABORT_TIME_LIM) {

	} else {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "ERROR. Unknown DEP SolutionStatus: " + extSol::QpExternSolver::solutionStatusToString(sol.solution.status));
	}
	return sol;
}
}
