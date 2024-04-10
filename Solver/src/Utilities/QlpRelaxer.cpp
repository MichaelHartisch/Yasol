/*
*
* Solver: QlpRelaxer.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Utilities/QlpRelaxer.hpp"
#include "Algorithms/Algorithm.hpp"
#include "Utilities/QlpConverter.hpp"
//#include "ilcplex/ilocplex.h"

namespace utils {
std::string QlpRelaxer::LOG_TAG = "QlpRelaxer";
QlpRelaxer::QlpRelaxer(const data::Qlp& source, bool LP) :
		qlpRelComp(false), lbComp(false), ubComp(false), qlp(source), qlpSplitter(qlp), qlpRelSol(), lbSol(), ubSol(), extSol(NULL), extSolLp(NULL), scenario(0), infeasible(0), wcFlag(false) {

	//qlp.sortQlp();//Needed here for some reason

	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Splitting equalities...");
	utils::QlpConverter::splitEqualities(qlp);

	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing Splitter...");
	qlpSplitter.initSplitter(algorithm::Algorithm::WORST_CASE);

	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing ExternSolver...");
	this->extSol = extSol::initExternSolver();
	this->extSol->init(data::QpObjFunc(this->qlpSplitter.objective, this->qlpSplitter.existObjVec, this->qlpSplitter.existOffset), this->qlpSplitter.existVars, this->qlpSplitter.fastExistMatrix, this->qlpSplitter.qlp.getRhsVec());
	this->extSol->adaptToSolverMode(extSol::QpExternSolver::RELAXER);

	this->extSolLp = extSol::initNbdExternSolver();
	this->extSolLp->init(qlp);
	this->extSolLp->adaptToSolverMode(extSol::QpExternSolver::RELAXER);

	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Done.");
}

QlpRelaxer::~QlpRelaxer() {
	delete extSol;
	delete extSolLp;
}

std::pair<data::QpNum, data::QpNum> QlpRelaxer::computeBounds() {
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Getting bounds...");
	data::QpNum ub = getUpperBound();
	data::QpNum lb = getLowerBound();
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished.");
	return std::make_pair(lb, ub);
}

data::QpNum QlpRelaxer::getLpRelaxationBound(){
	data::QpNum v;
	this->extSolLp->solve(1000000, 100000);
	if (this->extSolLp->getSolutionStatus() != extSol::QpExternSolver::INFEASIBLE) {
		v = this->extSolLp->getObjValue();
	}else{
		v.setMaxInf();
	}
	return v;
}

data::QpNum QlpRelaxer::getQlpRelaxationBound() {
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Getting Qlp relaxation bound...");
	if (!qlpRelComp) {
		this->computeQlpRelaxationBound();
	}
	return this->qlpRelSol.ofVal;
}

std::vector<data::QpNum> QlpRelaxer::getLpRelaxationVarAlloc(){
	std::vector<data::QpNum> v;
	this->extSolLp->solve(1000000, 100000);
	if (this->extSolLp->getSolutionStatus() != extSol::QpExternSolver::INFEASIBLE) {
			this->extSolLp->getValues(v);
	}
	return v;
}

const std::vector<data::QpNum>& QlpRelaxer::getQlpRelaxationVarAlloc() {
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Getting Qlp relaxation variable allocation...");
	if (!qlpRelComp)
		this->computeLowerBound();
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished.");
	return this->qlpRelSol.varAlloc;
}

data::QpNum QlpRelaxer::getLowerBound() {
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Getting lower bound...");
	if (!lbComp) {
		this->computeLowerBound();
	}
	return this->lbSol.ofVal;
}

data::QpNum QlpRelaxer::getUpperBound() {
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Getting upper bound...");
	if (!ubComp) {
		this->computeUpperBound();
	}
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished.");
	return this->ubSol.ofVal;
}

const std::vector<data::QpNum>& QlpRelaxer::getLowerBoundVarAlloc() {
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Getting lower bound variable allocation...");
	if (!lbComp)
		this->computeLowerBound();
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished.");
	return this->lbSol.varAlloc;
}

const std::vector<data::QpNum>& QlpRelaxer::getUpperBoundVarAlloc() {
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Getting upper bound variable allocation...");
	if (!ubComp)
		this->computeUpperBound();
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished.");
	return this->ubSol.varAlloc;
}

void QlpRelaxer::computeLowerBound() {
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Computing lower bound...");

	//Collect all universal variable ranges
	data::QpNumMatrix currUnivVarRanges;

	unsigned int genExistVars = this->qlp.getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::generals) + this->qlp.getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::binaries);

	for (unsigned int i = 0; i < this->qlpSplitter.univVarCount; i++) {

		if (genExistVars && this->qlpSplitter.univVars[i].getNumberSystem() == data::QpVar::generals) {
			std::vector<data::QpNum> tmpVec(1 + this->qlpSplitter.univVars[i].getUpperBound().asDouble() - this->qlpSplitter.univVars[i].getLowerBound().asDouble(), 0);
			int val = (int) this->qlpSplitter.univVars[i].getLowerBound().asDouble();
			for (unsigned int j = 0; j < tmpVec.size(); j++) {
				tmpVec[j] = val;
				val += 1;
			}
			currUnivVarRanges.push_back(tmpVec);
		} else {
			currUnivVarRanges.push_back(this->qlpSplitter.univVars[i].getVariableRange());
		}

	}
	std::vector<data::QpNum> v;
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Computing Scenarios...");
	this->solveScenarios(v, currUnivVarRanges.begin(), currUnivVarRanges.end());
	lbComp = true;
}

bool QlpRelaxer::solveScenarios(std::vector<data::QpNum>& currScen, std::vector<std::vector<data::QpNum> >::const_iterator currIt, std::vector<std::vector<data::QpNum> >::const_iterator end) {

	if (currIt == end) {
		if (LOG_QLPRELAXER)
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Computing ScenarioSoution: " + utils::ToolBox::convertToString(scenario));
		return this->computeScenarioSolution(currScen);
	}

	const std::vector<data::QpNum>& tmpVec = *currIt;
	for (std::vector<data::QpNum>::const_iterator it = tmpVec.begin(); it != tmpVec.end(); it++) {
		currScen.push_back(*it);
		if (!solveScenarios(currScen, currIt + 1, end))
			return false;
		currScen.pop_back();
	}
	return true;
}

bool QlpRelaxer::computeScenarioSolution(const std::vector<data::QpNum>& input) {
	extSol::QpExternSolver::QpSolution scenSol;
	std::vector<data::QpNum> scenRhs = this->qlpSplitter.originalRhs;

	data::QpNum value, offSet;

	if (this->qlpSplitter.univObjVec.size()) {
		for (unsigned int i = 0, size = this->qlpSplitter.univObjVec.size(); i < size; i++) {
			offSet += this->qlpSplitter.univObjVec[i] * input[i];
		}
	}

	for (unsigned int i = 0, size = this->qlpSplitter.fastUnivMatrix.size(); i < size; i++) {
if(this->qlpSplitter.fastUnivMatrix[i].size())
		for (unsigned int j = 0, elems = this->qlpSplitter.fastUnivMatrix[i].size(); j < elems; j++) {
			data::IndexedElement& ie = this->qlpSplitter.fastUnivMatrix[i][j];
			scenRhs[i] -= ie.value * input[ie.index];
		}
	}

	for (unsigned int i = 0; i < scenRhs.size(); i++) {
		this->extSol->changeRhsElement(i, scenRhs[i]);
	}

	if (LOG_QLPRELAXER) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Solving scenario: " + data::QpNum::vecToString(input));

	}
	if ((scenSol.status = this->extSol->solve()) == extSol::QpExternSolver::INFEASIBLE) {
		if (LOG_QLPRELAXER){
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Infeasible. Scenario: "+utils::ToolBox::convertToString(scenario));
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Infeasible. Input   : "+data::QpNum::vecToString(input));
		}
		lbSol.ofVal.setMinInf();
		infeasible++;
		scenario++;
		return false;
	} else {
		if (LOG_QLPRELAXER)
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Feasible: "+this->extSol->getObjValue().toString());
		scenSol.ofVal = this->extSol->getObjValue();
		scenSol.ofVal += offSet;
		if (scenario == 0) {
			this->lbSol = scenSol;
		} else if (this->qlpSplitter.objective == data::QpObjFunc::min) {
			if (scenSol.ofVal > lbSol.ofVal) {
				this->lbSol.ofVal = scenSol.ofVal;
				this->extSol->getValues(lbSol.varAlloc);

			}
		} else {
			if (scenSol.ofVal < lbSol.ofVal) {
				this->lbSol.ofVal = scenSol.ofVal;
				this->extSol->getValues(lbSol.varAlloc);
			}
		}
		scenario++;
		return true;
	}
}

void QlpRelaxer::computeUpperBound() {
	if (LOG_QLPRELAXER)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Computing upper bound...");

	data::QpNum offSet;

	if (!wcFlag) {
		std::vector<data::QpNum> wcRhs = this->qlpSplitter.qlp.getRhsValVec();
		const std::vector<data::QpRhs>& rhs = this->qlpSplitter.qlp.getRhsVec();
		std::pair<data::QpNum, data::QpNum> bounds;
		data::QpNum coeff;
		for (unsigned int i = 0, rows = this->qlpSplitter.fastUnivMatrix.size(); i < rows; i++) {
			for (unsigned int j = 0, cols = this->qlpSplitter.fastUnivMatrix[i].size(); j < cols; j++) {
				bounds = this->qlpSplitter.univVars[this->qlpSplitter.fastUnivMatrix[i][j].index].getBounds();
				coeff = qlpSplitter.fastUnivMatrix[i][j].value;
				if (rhs[i].getRatioSign() == data::QpRhs::smallerThanOrEqual) {
					wcRhs[i] -= ((coeff < 0) ? bounds.first : bounds.second) * coeff;
				} else {
					wcRhs[i] -= ((coeff < 0) ? bounds.second : bounds.first) * coeff;
				}
			}
		}

		for (unsigned int i = 0; i < wcRhs.size(); i++) {
			this->extSol->changeRhsElement(i, wcRhs[i]);
		}
		wcFlag = true;
	}

	if (this->qlpSplitter.univObjVec.size()) {
		for (unsigned int i = 0, size = this->qlpSplitter.univObjVec.size(); i < size; i++) {
			//if (!this->qlpSplitter.univObjVec[i].isZero()){
				data::QpNum tmp(this->qlpSplitter.univObjVec[i]);
				tmp*=((this->qlpSplitter.univObjVec[i] < 0) ? this->qlpSplitter.univVars[i].getLowerBound() : this->qlpSplitter.univVars[i].getUpperBound());
				offSet += tmp;
			//}
		}
	}

	if ((ubSol.status = this->extSol->solve()) == extSol::QpExternSolver::INFEASIBLE) {
		ubSol.ofVal.setMaxInf();
	} else {
		ubSol.ofVal = this->extSol->getObjValue();
		ubSol.ofVal += offSet;
		this->extSol->getValues(ubSol.varAlloc);
	}
	ubComp = true;
}

void QlpRelaxer::computeQlpRelaxationBound() {
	const std::vector<data::QpNum>& wcRhs = this->qlpSplitter.qlp.getRhsValVec();
	for (unsigned int i = 0; i < wcRhs.size(); i++)
		this->extSol->changeRhsElement(i, wcRhs[i]);
	if ((qlpRelSol.status = this->extSol->solve()) == extSol::QpExternSolver::INFEASIBLE) {
		qlpRelSol.ofVal.setMaxInf();
	} else {
		qlpRelSol.ofVal = this->extSol->getObjValue();
		this->extSol->getValues(qlpRelSol.varAlloc);
	}
	qlpRelComp = true;
}

}
