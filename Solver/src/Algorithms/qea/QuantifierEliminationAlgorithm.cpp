/*
*
* Solver: QuantifierEliminationAlgorithm -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Algorithms/qea/QuantifierEliminationAlgorithm.hpp"
#include "Utilities/QlpConverter.hpp"

namespace algorithm {

std::string QuantifierEliminationAlgorithm::LOG_TAG = "QEA";

QuantifierEliminationAlgorithm::QuantifierEliminationAlgorithm(data::Qlp& qp) :
		Algorithm(qp, this->qae), variable(), constraints(), cLessThanOrEqual(), cGreaterThanOrEqual(), cEquals(), cRest() {
	//TODO check the input
}

QuantifierEliminationAlgorithm::~QuantifierEliminationAlgorithm() {
	if (LOG_QEA)
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "Destroyed QuantifierEliminationAlgorithm");
}

Algorithm::QlpSolution QuantifierEliminationAlgorithm::solveQlp(SolutionCase s) {

	Algorithm::QlpSolution sol;
	sol.sc = s;
	sol.obj = qlpWork.getObjective();

	switch (s) {
	case FEASIBILITY: {
		bool feasible = this->checkFeasibility();
		sol.solution.ofVal = 0;
		sol.solution.status = feasible ? extSol::QpExternSolver::OPTIMAL : extSol::QpExternSolver::INFEASIBLE;
		break;
	}
	case WORST_CASE: {

		data::Qlp tmpQlp(this->qlpWork);

		utils::QlpConverter::pushObjectiveFunctionToMatrixFront(tmpQlp, this->qlpWork);

		bool feasible = this->checkFeasibility();

		sol.solution.ofVal = feasible ? this->qlpWork.getVariableByIndex(0).getLowerBound().toString() : 0;
		sol.solution.status = feasible ? extSol::QpExternSolver::OPTIMAL : extSol::QpExternSolver::INFEASIBLE;
		std::cout << "Status: " << extSol::QpExternSolver::solutionStatusToString(sol.solution.status) << " \nObj: " << sol.solution.ofVal.toString() << std::endl;
		break;
	}
	case AVERAGE_CASE: {
		throw utils::AlgorithmException("QuantifierEliminationAlgorithm::solveQlp(SolutionCase s) --> SolutionCase not supported");
	}
	default: {
		throw utils::AlgorithmException("QuantifierEliminationAlgorithm::solveQlp(SolutionCase s) --> unknown SolutionCase");
	}
	}
	return sol;
}

bool QuantifierEliminationAlgorithm::checkFeasibility() {

	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "---------------------------------------------------------------------------------------------------------------------------->\n");

	//Normalize each Constraint, ... + 1x_n holds for each last constraint variable coefficient
	qlpWork.normalizeConstraints();

	data::Constraint* tmp = NULL;
	data::Constraint* outer = NULL;
	data::QpNum lb;
	data::QpNum ub;

	unsigned int count = qlpWork.getVariableCount() - 1;

	for (unsigned int i = count, iteration = 0; i > 0; --i, iteration++) {

		unsigned int varIndex = i;

		variable = &qlpWork.getVariableByIndex(i);
		lb = variable->getLowerBound();
		ub = variable->getUpperBound();

		if (LOG_QEA) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "---------------------------------------> ITERATION: " + utils::ToolBox::convertToString(iteration));
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, " ");
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Variable to eliminate                    : " + variable->toString());
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Current number of constraints:Size before: " + utils::ToolBox::convertToString(qlpWork.getConstraintCount()));
		}

		outer = qlpWork.getFirstConstraint();
		while (outer) {
			if (outer->getElementCount() <= 1) {
				tmp = outer;
				outer = outer->getNextConstraintPtr();
				if (!this->checkConstraint(*tmp))
					return false;
			} else {
				outer = outer->getNextConstraintPtr();
			}
		}
		outer = NULL;

		if (variable->getQuantifier() == data::QpVar::exists) {
			qlpWork.createRhsConstraint(data::QpRhs::smallerThanOrEqual, ub).createConstraintElement(varIndex, 1);
			qlpWork.createRhsConstraint(data::QpRhs::greaterThanOrEqual, lb).createConstraintElement(varIndex, 1);
		}

		qlpWork.getSortedVariableConstraints(varIndex, constraints, cLessThanOrEqual, cGreaterThanOrEqual, cEquals, cRest);

		if (LOG_QEA) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "List sizes: ");
			this->printListSizes();
		}

		if(lb==ub){
			for(unsigned i = 0; i < cEquals.size();i++){
				cEquals[i]->eliminateConstraintElement(varIndex,lb);
			}
			for(unsigned i = 0; i < cLessThanOrEqual.size();i++){
				cLessThanOrEqual[i]->eliminateConstraintElement(varIndex,lb);
						}
			for(unsigned i = 0; i < cGreaterThanOrEqual.size();i++){
				cGreaterThanOrEqual[i]->eliminateConstraintElement(varIndex,lb);
			}
		}else if (cEquals.size()) {
			unsigned int size = constraints.size();
			data::Constraint* replacement = cEquals[0];
			data::Constraint* constraint;
			for (unsigned int i = 0; i < size; ++i) {
				if ((constraint = constraints[i]) != replacement) {
					constraint->eliminateConstraintElement(varIndex, *replacement);
					if (!this->checkConstraintConsistency(constraint))
						return false;
				}
			}
			cEquals[0]->removeConstraint();
		} else {
			if (variable->getQuantifier() == data::QpVar::exists) {
				if (!eliminateExistentialVariable(*variable))
					return false;
			} else if (variable->getQuantifier() == data::QpVar::all) {

				data::Constraint* tmp = NULL;
				data::QpNum lb = variable->getLowerBound();
				data::QpNum ub = variable->getUpperBound();

				unsigned int size = cLessThanOrEqual.size();
				for (unsigned int i = 0; i < size; ++i) {
					tmp = cLessThanOrEqual[i];
					tmp->eliminateConstraintElement(varIndex, ub);
					tmp->normalizeConstraint();
				}

				size = cGreaterThanOrEqual.size();
				for (unsigned int i = 0; i < size; ++i) {
					tmp = cGreaterThanOrEqual[i];
					tmp->eliminateConstraintElement(varIndex, lb);
					tmp->normalizeConstraint();
				}

			} else {
				throw utils::AlgorithmException("QuantifierEliminationAlgorithm::checkFeasibility() --> unsupported quantifier in elimination step");
			}
		}

		constraints.clear();
		cLessThanOrEqual.clear();
		cGreaterThanOrEqual.clear();
		cEquals.clear();
		cRest.clear();

		if (!qlpWork.getConstraintCount())
			break;

	}

	return this->checkResultingSystem();

}

bool QuantifierEliminationAlgorithm::eliminateExistentialVariable(data::QpVar& v) {

	unsigned int newSize = (cLessThanOrEqual.size() * cGreaterThanOrEqual.size());

	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "\tElements to create (with NULL): " + utils::ToolBox::convertToString(newSize));

	std::vector<std::vector<data::IndexedElement> > lessRows;

	std::vector<std::vector<data::IndexedElement> > greaterRows;

	std::vector<data::QpRhs*> newConRhs;

	std::vector<std::vector<data::IndexedElement>*> newConVec;

	unsigned int i, j, k, size1, size2, size3;
	for (i = 0, size1 = cLessThanOrEqual.size(); i < size1; ++i) {
		if (cLessThanOrEqual[i] != NULL)
			lessRows.push_back(cLessThanOrEqual[i]->getElements());
	}

	for (i = 0, size1 = cGreaterThanOrEqual.size(); i < size1; ++i) {
		if (cGreaterThanOrEqual[i] != NULL)
			greaterRows.push_back(cGreaterThanOrEqual[i]->getElements());
	}

	data::QpRhs tmpRhs;
	std::vector<data::QpNum> tmpCon(v.getIndex() + 1, 0.0);
	data::IndexedElement ie;

	data::QpRhs* newRhsPointer;
	std::vector<data::IndexedElement> *ieVecPointer;

	data::QpNum value, rhs;

	data::QpVar* var;

	unsigned int iteration = 0;

	unsigned int newElemIndex = 0, newElems = 0;
	for (i = 0, size1 = cLessThanOrEqual.size(); i < size1; ++i) {

		if (cLessThanOrEqual[i] == NULL)
			continue;

		for (j = 0, size2 = cGreaterThanOrEqual.size(); j < size2; ++j) {

			if (iteration % 1000000 == 0) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Iteration: " + utils::ToolBox::convertToString(iteration));
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "New Constraints: " + utils::ToolBox::convertToString(newElems));
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Redundant: " + utils::ToolBox::convertToString((1.0 - double(newElems) / iteration) * 100.0));
			}

			if (cGreaterThanOrEqual[j] == NULL)
				continue;

			iteration++;

			for (unsigned int tmpConPos = 0, tmpConSize = v.getIndex() + 1; tmpConPos < tmpConSize; tmpConPos++)
				tmpCon[tmpConPos].setZero();

			newRhsPointer = new data::QpRhs(cLessThanOrEqual[i]->getRhsValue() - cGreaterThanOrEqual[j]->getRhsValue(), data::QpRhs::smallerThanOrEqual);
			for (unsigned int k = 0, size3 = lessRows[i].size(); k < size3; k++)
				tmpCon[lessRows[i][k].index] = lessRows[i][k].value;

			for (unsigned int k = 0, size3 = greaterRows[j].size(); k < size3; k++)
				tmpCon[greaterRows[j][k].index] -= greaterRows[j][k].value;

			ieVecPointer = new std::vector<data::IndexedElement>;

			for (unsigned int k = 0, size3 = tmpCon.size(); k < size3; k++) {
				if (!((ie.value = tmpCon[k]).isZero())) {
					ie.index = k;
					ieVecPointer->push_back(ie);
				}
			}

			if (ieVecPointer->size() == 0) {
				if (newRhsPointer->getRatioSign() == data::QpRhs::smallerThanOrEqual) {
					if (newRhsPointer->getValue() < 0.0)
						return false;
				} else if (newRhsPointer->getRatioSign() == data::QpRhs::greaterThanOrEqual) {
					if (newRhsPointer->getValue() > 0.0)
						return false;
				} else {
					if (newRhsPointer->getValue() != 0.0)
						return false;
				}

				delete newRhsPointer;
				delete ieVecPointer;
				continue;
			} else {
				value = ieVecPointer->operator [](ieVecPointer->size() - 1).value;
				newRhsPointer->operator /=(value);
				for (k = 0, size3 = ieVecPointer->size(); k < size3; k++)
					ieVecPointer->operator [](k).value /= value;

			}

			if (ieVecPointer->size() == 1) {
				var = &qlpWork.getVariableByIndex(ieVecPointer->operator [](0).index);
				rhs = newRhsPointer->getValue();

				if (newRhsPointer->getRatioSign() == data::QpRhs::smallerThanOrEqual) {
					if (var->getUpperBound() > rhs) {
						if (var->getQuantifier() != data::QpVar::exists) {
							return false;
						}
						var->setUpperBound(rhs);
						if (var->getLowerBound() == rhs)
							qlpWork.eliminateVariable(var->getIndex(), rhs);
					}
				} else if (newRhsPointer->getRatioSign() == data::QpRhs::greaterThanOrEqual) {
					if (var->getLowerBound() < rhs) {
						if (var->getQuantifier() != data::QpVar::exists) {
							return false;
						}
						var->setLowerBound(rhs);
						if (rhs == var->getUpperBound())
							qlpWork.eliminateVariable(var->getIndex(), rhs);
					}
				} else {
					if (rhs != var->getLowerBound() || rhs != var->getUpperBound()) {
						if (var->getQuantifier() != data::QpVar::exists)
							return false;
						var->setLowerBound(rhs);
						var->setUpperBound(rhs);
						qlpWork.eliminateVariable(var->getIndex(), rhs);
					}
				}
				if (var->getLowerBound() > var->getUpperBound())
					return false;

				delete newRhsPointer;
				delete ieVecPointer;
				continue;
			}
			newConRhs.push_back(newRhsPointer);
			newConVec.push_back(ieVecPointer);
			newElems++;
			newElemIndex++;
		}
	}

	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Removing old ...");
	for (unsigned int i = 0; i < cLessThanOrEqual.size(); ++i) {
		if (cLessThanOrEqual[i])
			cLessThanOrEqual[i]->removeConstraint();
	}
	for (unsigned int i = 0; i < cGreaterThanOrEqual.size(); ++i) {
		if (cGreaterThanOrEqual[i])
			cGreaterThanOrEqual[i]->removeConstraint();
	}
	for (unsigned int i = 0; i < cEquals.size(); ++i) {
		if (cEquals[i])
			cEquals[i]->removeConstraint();
	}

	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Adding new...");
	data::Constraint* newConstraint;
	//add the resulting new constraints
	unsigned int added = 0;
	if ((size1 = newConVec.size())) {
		for (i = 0; i < size1 - 1; i++) {
			if (i % 1000000 == 0)
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Added: " + utils::ToolBox::convertToString(i));
			if (newConVec[i] == NULL)
				continue;
			added++;
			newConstraint = &qlpWork.createRhsConstraint(newConRhs[i]->getRatioSign(), newConRhs[i]->getValue());
			for (j = 0, size2 = newConVec[i]->size(); j < size2; j++) {
				newConstraint->createConstraintElement(newConVec[i]->operator [](j).index, newConVec[i]->operator [](j).value);
			}
		}
		delete newConVec[i];
		delete newConRhs[i];
	}

	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished. Added: " + utils::ToolBox::convertToString(added));
	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Redundant: " + utils::ToolBox::convertToString(iteration - newElems));
	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "New Size: " + utils::ToolBox::convertToString(qlpWork.getConstraintCount()));

	return true;
}

bool QuantifierEliminationAlgorithm::eliminateUniversalVariable(data::QpVar& v) {

	data::Constraint* tmp = NULL;

	unsigned int varIndex = v.getIndex();
	data::QpNum lb = variable->getLowerBound();
	data::QpNum ub = variable->getUpperBound();

	unsigned int size = cLessThanOrEqual.size();
	for (unsigned int i = 0; i < size; ++i) {
		tmp = cLessThanOrEqual[i];
		tmp->eliminateConstraintElement(varIndex, ub);
		tmp->normalizeConstraint();
	}

	size = cGreaterThanOrEqual.size();
	for (unsigned int i = 0; i < size; ++i) {
		tmp = cGreaterThanOrEqual[i];
		tmp->eliminateConstraintElement(varIndex, lb);
		tmp->normalizeConstraint();
	}
	return true;
}

bool QuantifierEliminationAlgorithm::changeBounds(data::Constraint& c) {

	if (c.getElementCount() != 1)
		throw utils::AlgorithmException("Preprocessor::changeBounds(data::Constraint& c) --> Constraint has not a single element");

	data::QpVar& var = qlpWork.getVariableByIndex(c.getFirstConstraintElementIndex());

	c.normalizeConstraint();
	data::QpNum rhs = c.getRhsValue();

	if (c.getRhsRatioSign() == data::QpRhs::smallerThanOrEqual) {
		if (var.getUpperBound() > rhs) {
			if (var.getQuantifier() != data::QpVar::exists) {
				return false;
			}
			var.setUpperBound(rhs);
			if (var.getLowerBound() == rhs)
				this->qlpWork.eliminateVariable(var.getIndex(), rhs);
		}
	} else if (c.getRhsRatioSign() == data::QpRhs::greaterThanOrEqual) {
		if (var.getLowerBound() < rhs) {
			if (var.getQuantifier() != data::QpVar::exists) {
				return false;
			}
			var.setLowerBound(rhs);
			if (rhs == var.getUpperBound())
				this->qlpWork.eliminateVariable(var.getIndex(), rhs);
		}
	} else {
		if (rhs != var.getLowerBound() || rhs != var.getUpperBound()) {
			if (var.getQuantifier() != data::QpVar::exists)
				return false;
			var.setLowerBound(rhs);
			var.setUpperBound(rhs);
			this->qlpWork.eliminateVariable(var.getIndex(), rhs);
		}
	}
	return (var.getLowerBound() <= var.getUpperBound());
}

bool QuantifierEliminationAlgorithm::checkConstraintConsistency(data::Constraint* c) {
	if (c->getElementCount())
		return true;
	if (c->getRhsRatioSign() == data::QpRhs::smallerThanOrEqual) {
		return (c->getRhsValue() >= 0.0);
	} else if (c->getRhsRatioSign() == data::QpRhs::greaterThanOrEqual) {
		return (c->getRhsValue() <= 0.0);
	} else {
		return (c->getRhsValue() == 0.0);
	}
}

void QuantifierEliminationAlgorithm::pruneConstraints(std::vector<data::Constraint*>& cVec) {
	if (LOG_QEA)
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "QuantifierEliminationAlgorithm::pruneConstraints(std::vector<data::Constraint*>& cVec)");
	data::Constraint* tmp;
	std::vector<data::Constraint*> tmpVec;

	unsigned int outerSize = cVec.size() - 1;
	unsigned int innerSize = outerSize + 1;

	for (unsigned int i = 0; i < outerSize; ++i) {
		if ((tmp = cVec[i]) == NULL)
			// the outer constraint has already been eliminated
			continue;
		for (unsigned int j = i + 1; j < innerSize; ++j) {
			if (cVec[j] == NULL)
				// the inner constraint has already been eliminated
				continue;
			if (tmp->operator ==(*cVec[j]) && (tmp != cVec[j])) {
				// prune the constraints to reduce the resulting amount of new constraints after
				// the elimination  step
				cVec[j]->removeConstraint();
				cVec[j] = NULL;
				continue;
			}
		}
		tmpVec.push_back(tmp);
	}
	if (cVec[cVec.size() - 1] != NULL)
		tmpVec.push_back(cVec[cVec.size() - 1]);
	cVec.assign(tmpVec.begin(), tmpVec.end());
}

//-----------------------------------------Constraint Validation--------------------------------------------->
bool QuantifierEliminationAlgorithm::checkConstraint(data::Constraint& c) {
	if (LOG_QEA)
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "QuantifierEliminationAlgorithm::checkConstraint(data::Constraint& c)");
	if (c.getElementCount() == 0) {
		if (!this->checkConstraintConsistency(&c)) {
			if (LOG_QEA)
				utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "\t\t checkConstraintConsistency(&c) == false : " + c.toString());
			return false;
		} else {
			if (LOG_QEA)
				utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "\t\t Removed empty constraint: " + c.toString());
			c.removeConstraint();
		}
	} else if (c.getElementCount() == 1) {
		if (!this->changeBounds(c)) {
			if (LOG_QEA)
				utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "pre.changeBounds(c) == false : " + c.toString());
			return false;
		} else {
			if (LOG_QEA)
				utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "\t\t Removed implicit bound: " + c.toString());
			c.removeConstraint();
		}
	}
	return true;
}

bool QuantifierEliminationAlgorithm::checkResultingSystem() {
	if (LOG_QEA)
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "QuantifierEliminationAlgorithm::checkResultingSystem(): " + qlpWork.toString());

	data::Constraint* c = qlpWork.getFirstConstraint();
	while (c != NULL) {
		this->checkConstraint(*c);
		c = c->getNextConstraintPtr();
	}
	data::QpVar & variable = qlpWork.getVariableByIndex(0);
	data::QpNum maxUb = variable.getUpperBound();
	data::QpNum minLb = variable.getLowerBound();
	if (LOG_QEA) {
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "Maximum: " + maxUb.toString());
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "Minimum: " + minLb.toString());
	}
	return maxUb >= minLb;
}

bool QuantifierEliminationAlgorithm::checkDuplicateConstraint(data::Constraint& c) {
	data::Constraint* tmp;
	for (unsigned int i = 0; i < cRest.size(); i++) {
		if ((tmp = cRest[i]) == NULL)
			continue;
		if (&c != tmp && c.operator ==(*tmp)) {
			if (LOG_QEA) {
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "\t\t Duplicate. New: " + c.toString());
				utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "\t\t Duplicate. Old: " + tmp->toString());
			}
			tmp->removeConstraint();
			cRest[i] = NULL;
			return true;
		}
	}
	return false;
}

void QuantifierEliminationAlgorithm::printListSizes() const {
	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "\tcEquals: " + utils::ToolBox::convertToString(cEquals.size()) + "\tcLessThanOrEqual: " + utils::ToolBox::convertToString(cLessThanOrEqual.size()) + "\tcGreaterThanOrEqual: " + utils::ToolBox::convertToString(cGreaterThanOrEqual.size()) + "\tcRest: " + utils::ToolBox::convertToString(cRest.size()));
}

void QuantifierEliminationAlgorithm::printListContent() const {
	utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "cEquals[]:");
	for (unsigned int i = 0; i < cEquals.size(); i++)
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "\t\tConstraint: " + (cEquals[i] ? cEquals[i]->toString() : "NULL"));
	utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "cLessThanOrEqual[]:");
	for (unsigned int i = 0; i < cLessThanOrEqual.size(); i++)
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "\t\tConstraint: " + (cLessThanOrEqual[i] ? cLessThanOrEqual[i]->toString() : "NULL"));
	utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "cGreaterThanOrEqual[]:");
	for (unsigned int i = 0; i < cGreaterThanOrEqual.size(); i++)
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "\t\tConstraint: " + (cGreaterThanOrEqual[i] ? cGreaterThanOrEqual[i]->toString() : "NULL"));
	for (unsigned int i = 0; i < cRest.size(); i++)
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "\t\tConstraint: " + (cRest[i] ? cRest[i]->toString() : "NULL"));
}

}
