/*
*
* Solver: QlpConverter.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Utilities/QlpConverter.hpp"

//#include <unordered_map>

#include "Datastructures/components/QpVar.hpp"
#include "Datastructures/qlp/Constraint.hpp"

#include "Algorithms/Algorithms.hpp"
#include "Utilities/ExternCode/QIP2QBP.h"

namespace utils {
std::string QlpConverter::LOG_TAG = "QlpConverter";

void QlpConverter::generateQlpFromLp(data::Qlp& in, data::Qlp& out) {

	utils::Rng randGen(time(NULL));

	unsigned int variables, blocks, variablesPerBlock;

	variables = UNIV_VARS;
	if ((blocks = UNIV_VAR_BLOCKS/*UNIV_VAR_PERCENTAGE_PER_BLOCK*/) <= variables) {
		variablesPerBlock = variables / blocks;
	} else {
		blocks = variablesPerBlock = variables;
	}

	std::vector<unsigned int> blockVarVec;
	unsigned int rest = variables;
	for (unsigned int i = 0; i < blocks; i++) {
		blockVarVec.push_back(variablesPerBlock);
		rest -= variablesPerBlock;
	}

	for (unsigned int i = 0, j = 0; i < rest; i++, j++)
		blockVarVec[(int) randGen.uniform(blockVarVec.size() - 1)]++;

	//utils::Logger::globalLog(utils::LOG_INFO, "QLPGEN", "blockVarVec: " + utils::ToolBox::vecToString(blockVarVec));

	std::vector<unsigned int> univBlockPos(in.getVariableCount(), 0);
	if(UNIV_VAR_BLOCKS > UNIV_VARS)
		throw utils::QlpSolverException("UNIV_VAR_BLOCKS > UNIV_VARS");

	for(unsigned int i = 1; i <= UNIV_VAR_BLOCKS; i++){
		unsigned int blockStartIndex = in.getVariableCount() / (UNIV_VAR_BLOCKS +1);
		univBlockPos[blockStartIndex*i] = 1;

	}



//	if (UNIV_VAR_PERCENTAGE_PER_BLOCK == 100) {
//		unsigned int blockStartIndex = in.getVariableCount() / 2;
//		univBlockPos[blockStartIndex] = 1;
//	} else if (UNIV_VAR_PERCENTAGE_PER_BLOCK == 50) {
//		unsigned int candidates = in.getVariableCount() / 3;
//		unsigned int blockStartIndex = candidates;
//		univBlockPos[blockStartIndex] = 1;
//		univBlockPos[blockStartIndex * 2] = 1;
//	} else if (UNIV_VAR_PERCENTAGE_PER_BLOCK == 20) {
//		unsigned int candidates = in.getVariableCount() / 6;
//		unsigned int blockStartIndex = candidates;
//		univBlockPos[blockStartIndex] = 1;
//		univBlockPos[blockStartIndex * 2] = 1;
//		univBlockPos[blockStartIndex * 3] = 1;
//		univBlockPos[blockStartIndex * 4] = 1;
//		univBlockPos[blockStartIndex * 5] = 1;
//	} else if (UNIV_VAR_PERCENTAGE_PER_BLOCK == 10) {
//		unsigned int candidates = in.getVariableCount() / 11;
//		unsigned int blockStartIndex = candidates;
//		univBlockPos[blockStartIndex] = 1;
//		univBlockPos[blockStartIndex * 2] = 1;
//		univBlockPos[blockStartIndex * 3] = 1;
//		univBlockPos[blockStartIndex * 4] = 1;
//		univBlockPos[blockStartIndex * 5] = 1;
//		univBlockPos[blockStartIndex * 6] = 1;
//		univBlockPos[blockStartIndex * 7] = 1;
//		univBlockPos[blockStartIndex * 8] = 1;
//		univBlockPos[blockStartIndex * 9] = 1;
//		univBlockPos[blockStartIndex * 10] = 1;
//	} else {
//		throw utils::AlgorithmException("UNIV_VAR_PERCENTAGE_PER_BLOCK Not supported");
//	}

	data::QpVar::NumberSystem numberSystem;
	if (UNIV_VAR_NUMBER_SYSTEM == "REAL")
		numberSystem = data::QpVar::real;
	else if (UNIV_VAR_NUMBER_SYSTEM == "GENERAL")
		numberSystem = data::QpVar::generals;
	else
		numberSystem = data::QpVar::binaries;

	// Compute universal bounds
	std::vector<data::QpNum> lowerBounds(variables, 0), upperBounds(variables, 0);
	for (unsigned int i = 0; i < variables; i++) {
		if (numberSystem == data::QpVar::binaries) {
			lowerBounds[i] = 0;
			upperBounds[i] = 1;
		} else {
			lowerBounds[i] = MIN_UNIV_VAR_BOUND;
			upperBounds[i] = MAX_UNIV_VAR_BOUND;
		}
	}

	//Create the existential and universal variables in the determined order
	unsigned int uIndex = 0;
	unsigned int existVariables = in.getVariableCount();
	std::string y;
	unsigned int blockVarIndex = 0, index = 0, univIndex = 0;
	for (unsigned int i = 0; i < existVariables; i++) {
		out.createVariable(in.getVariableNameByIndex(i), index++, data::QpVar::exists, in.getVariableByIndex(i).getNumberSystem(), in.getVariableByIndex(i).getLowerBound(), in.getVariableByIndex(i).getUpperBound());
		if (univBlockPos[i]) {
			for (unsigned int j = 0; j < blockVarVec[blockVarIndex]; j++) {
				y = "y" + utils::ToolBox::convertToString(uIndex++);
				out.createVariable(y, index++, data::QpVar::all, numberSystem, lowerBounds[univIndex], upperBounds[univIndex]);
				univIndex++;
			}
			blockVarIndex++;
		}
	}

	//Copy restrictions from input to output qlp
	data::Constraint *newCon, *oldCon, *tmpCon;
	std::vector<data::Constraint*> constraintVec = in.getConstraintVec();
	std::vector<data::IndexedElement> rowElements;
	for (unsigned int i = 0, size = constraintVec.size(); i < size; i++) {
		oldCon = constraintVec[i];
		rowElements = oldCon->getElements();
		if (rowElements.size()) {
			newCon = &out.createRhsConstraint(oldCon->getRhsRatioSign(), oldCon->getRhsValue(), oldCon->getResponsibility());
			for (unsigned int j = 0, rowSize = rowElements.size(); j < rowSize; j++) {
				if (!rowElements[j].value.isZero())
					newCon->createConstraintElement(out.getVariableIndexByName(in.getVariableNameByIndex(rowElements[j].index)), rowElements[j].value);
			}
		}
	}

	std::vector<data::QpVar*> univVarVec = out.getVariableVectorByQuantifier(data::QpVar::all);
	constraintVec = out.getConstraintVec();
	unsigned int constraints = constraintVec.size(), colElems, rowIndex, univVarIndex;
	data::QpNum value;
	unsigned int nonZeroUnivCoeffs = 0;
	for (unsigned int i = 0, size = univVarVec.size(); i < size; i++) {
		colElems = randGen.uniform((constraints * MIN_UNIV_VAR_COL_DENSITY_PERCENTAGE) / 100.0, (constraints * MAX_UNIV_VAR_COL_DENSITY_PERCENTAGE) / 100.0);
		if (colElems == 0)
			colElems = 1;

		if (LOG_QLPGEN)
			utils::Logger::globalLog(utils::LOG_INFO, "QlpGenerator", "\t ColsElems: " + utils::ToolBox::convertToString(colElems));

		for (unsigned int j = 0; j < colElems; j++) {

			do {
				rowIndex = randGen.uniform(constraints - 1);
			} while (rowIndex >= constraints);
			tmpCon = constraintVec[rowIndex];
			if (tmpCon->getRhsRatioSign() != data::QpRhs::equal || tmpCon->getLastConstraintElementIndex() > univVarVec[i]->getIndex()) {
				univVarIndex = univVarVec[i]->getIndex();
				if (tmpCon->getConstraintElementCoeff(univVarIndex).isMinInf()) {
					do {
						value = utils::ToolBox::round(randGen.uniform(100 * MIN_UNIV_VAR_MATRIX_COEFF, 100 * MAX_UNIV_VAR_MATRIX_COEFF) / 100.0, UNIV_MATRIX_COEFF_TAIL);
					} while (value == 0);
					tmpCon->createConstraintElement(univVarIndex, value);
					nonZeroUnivCoeffs++;
				} else {
					j--;
				}
			} else {
				j--;
			}
		}
	}

	out.getObjectiveFunction().setSize(out.getVariableCount());
	std::vector<data::QpVar*> varVec = out.getVariableVector();
	for (unsigned int i = 0, size = varVec.size(); i < size; i++) {
		if (varVec[i]->getQuantifier() == data::QpVar::exists) {
			index = in.getVariableIndexByName(out.getVariableNameByIndex(varVec[i]->getIndex()));
			out.setObjectiveFunctionElement(i, in.getObjectiveFunctionElement(index));
		} else {
			if (UNIV_VAR_OBJ_COEFFS) {
				do {
					value = utils::ToolBox::round(randGen.uniform(100 * MIN_UNIV_VAR_MATRIX_COEFF, 100 * MAX_UNIV_VAR_MATRIX_COEFF) / 100.0, UNIV_MATRIX_COEFF_TAIL);
				} while (value == 0);
				out.setObjectiveFunctionElement(i, value);
			}
		}
	}
}

void QlpConverter::addArtificialVariable(const data::Qlp& input, data::Qlp& out, const data::QpNum& penalty, bool pushToBack) {

	data::Qlp tmp(input);
	utils::QlpConverter::splitEqualities(tmp);
	std::vector<data::IndexedElement> elems;

	std::vector<data::Constraint*> cVec = tmp.getConstraintVec();
	for (unsigned int i = 0; i < cVec.size(); i++) {
		if (cVec[i]->getRhsRatioSign() == data::QpRhs::equal) {
			cVec[i]->getElementsSparse(elems);
			for (unsigned j = 0; j < elems.size(); j++) {
				if (tmp.getVariableByIndexConst(elems[j].index).getQuantifier() == data::QpVar::all) {
					data::QpNum rhsVal(cVec[i]->getRhsValue());
					cVec[i]->setRhsRatioSign(data::QpRhs::greaterThanOrEqual);
					data::Constraint* tmpCon = &tmp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, rhsVal, cVec[i]->getResponsibility());
					int index = cVec[i]->getFirstConstraintElementIndex();
					while (index >= 0) {
						tmpCon->createConstraintElement(index, cVec[i]->getConstraintElementCoeff(index));
						index = cVec[i]->getNextConstraintElementIndex(index);
					}
					cVec[i]->setRhsValue(rhsVal);
					tmpCon->setRhsValue(rhsVal);
					break;
				}
			}
		}
	}

	unsigned int zIndex = 0;
	if (pushToBack) {
		out = tmp;
		zIndex = tmp.getVariableCount();
		data::QpVar v("Z", zIndex, data::QpNum(0), data::QpNum(1000));
		out.createVariable(v);
		out.setObjectiveFunctionElement(zIndex, penalty);
	} else {
		out.clear();
		data::QpVar v("Z", zIndex, data::QpNum(0), data::QpNum(1000));
		out.createVariable(v);
		out.setObjectiveFunctionElement(zIndex, penalty);
		std::vector<data::QpVar*> vec = tmp.getVariableVector();
		for (unsigned int i = 0; i < vec.size(); i++) {
			vec[i]->setIndex(i + 1);
			out.createVariable(*vec[i]);
			out.setObjectiveFunctionElement(i + 1, tmp.getObjectiveFunctionElement(i));
		}

		std::vector<data::Constraint*> cVec = tmp.getConstraintVec();
		for (unsigned int i = 0; i < cVec.size(); i++) {
			data::Constraint& c = out.createRhsConstraint(cVec[i]->getRhs());
			std::vector<data::IndexedElement> elems = cVec[i]->getElements();
			for (unsigned j = 0; j < elems.size(); j++) {
				c.createConstraintElement(elems[j].index + 1, elems[j].value);
			}
		}
	}

	std::vector<unsigned int> univCons(tmp.getConstraintCount(), 1);
	cVec = tmp.getConstraintVec();
	for (unsigned int i = 0; i < cVec.size(); i++) {
		cVec[i]->getElementsSparse(elems);
		for (unsigned j = 0; j < elems.size(); j++) {
			if (tmp.getVariableByIndexConst(elems[j].index).getQuantifier() == data::QpVar::all) {
				univCons[i] = 1;
				break;
			}
		}
	}

	std::vector<data::Constraint*> vec = out.getConstraintVec();
	for (unsigned int i = 0; i < vec.size(); i++) {
		if (!univCons[i])
			continue;
		if (vec[i]->getRhsRatioSign() == data::QpRhs::smallerThanOrEqual) {
			vec[i]->createConstraintElement(zIndex, -1);
		} else if (vec[i]->getRhsRatioSign() == data::QpRhs::greaterThanOrEqual) {
			vec[i]->createConstraintElement(zIndex, 1);
		} else {
			vec[i]->createConstraintElement(zIndex, 1);
		}
	}
}

void QlpConverter::addRecourseVariable(const data::Qlp& in, data::Qlp& out, const data::QpNum& penalty) {

	out = in;

	//utils::QlpConverter::splitEqualities(out);

	std::vector<data::IndexedElement> elems;
	unsigned int recourseStartIndex = out.getVariableCount();
	std::vector<data::Constraint*> cVec = out.getConstraintVec();

	unsigned int zIndex = 0;

	unsigned int recourseVariableCount = 0;
	std::vector<unsigned int> univCons(out.getConstraintCount(), 0);
	cVec = out.getConstraintVec();
	for (unsigned int i = 0; i < cVec.size(); i++) {
		cVec[i]->getElementsSparse(elems);
		for (unsigned j = 0; j < elems.size(); j++) {
			if (out.getVariableByIndexConst(elems[j].index).getQuantifier() == data::QpVar::all) {
				univCons[i] = 1;
				recourseVariableCount++;
				if(cVec[i]->getRhsRatioSign() == data::QpRhs::equal){
					recourseVariableCount++;
				}

				break;
			}
		}
	}


	for(unsigned int i = 0; i < recourseVariableCount;i++){
	zIndex = out.getVariableCount();
	data::QpVar v("Z_"+utils::ToolBox::convertToString(i), zIndex, data::QpNum(0), data::QpNum(1000));
	out.createVariable(v);
	out.setObjectiveFunctionElement(zIndex, penalty);
	}

	//std::cout<< "Done."<<std::endl;

//	for (unsigned int i = 0; i < cVec.size(); i++) {
//		if (cVec[i]->getRhsRatioSign() == data::QpRhs::equal) {
//			cVec[i]->getElementsSparse(elems);
//			for (unsigned j = 0; j < elems.size(); j++) {
//				if (tmp.getVariableByIndexConst(elems[j].index).getQuantifier() == data::QpVar::all) {
//					data::QpNum rhsVal(cVec[i]->getRhsValue());
//					cVec[i]->setRhsRatioSign(data::QpRhs::greaterThanOrEqual);
//					data::Constraint* tmpCon = &tmp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, rhsVal);
//					int index = cVec[i]->getFirstConstraintElementIndex();
//					while (index >= 0) {
//						tmpCon->createConstraintElement(index, cVec[i]->getConstraintElementCoeff(index));
//						index = cVec[i]->getNextConstraintElementIndex(index);
//					}
//					cVec[i]->setRhsValue(rhsVal);
//					tmpCon->setRhsValue(rhsVal);
//					break;
//				}
//			}
//		}
//	}

//	unsigned int zIndex = 0;
//	if (pushToBack) {
//		out = tmp;
//		zIndex = tmp.getVariableCount();
//		data::QpVar v("Z", zIndex, data::QpNum(0), data::QpNum(1000));
//		out.createVariable(v);
//		out.setObjectiveFunctionElement(zIndex, penalty);
//	} else {
//		out.clear();
//		data::QpVar v("Z", zIndex, data::QpNum(0), data::QpNum(1000));
//		out.createVariable(v);
//		out.setObjectiveFunctionElement(zIndex, penalty);
//		std::vector<data::QpVar*> vec = tmp.getVariableVector();
//		for (unsigned int i = 0; i < vec.size(); i++) {
//			vec[i]->setIndex(i + 1);
//			out.createVariable(*vec[i]);
//			out.setObjectiveFunctionElement(i + 1, tmp.getObjectiveFunctionElement(i));
//		}
//
//		std::vector<data::Constraint*> cVec = tmp.getConstraintVec();
//		for (unsigned int i = 0; i < cVec.size(); i++) {
//			data::Constraint& c = out.createRhsConstraint(cVec[i]->getRhs());
//			std::vector<data::IndexedElement> elems = cVec[i]->getElements();
//			for (unsigned j = 0; j < elems.size(); j++) {
//				c.createConstraintElement(elems[j].index + 1, elems[j].value);
//			}
//		}
//	}



	std::vector<data::Constraint*> vec = out.getConstraintVec();
	unsigned int rsI = recourseStartIndex;

	for (unsigned int i = 0; i < vec.size(); i++) {
		if (!univCons[i])
			continue;
		if (vec[i]->getRhsRatioSign() == data::QpRhs::smallerThanOrEqual) {
			vec[i]->createConstraintElement(rsI, -1);
		} else if (vec[i]->getRhsRatioSign() == data::QpRhs::greaterThanOrEqual) {
			vec[i]->createConstraintElement(rsI, 1);
		} else {
			vec[i]->createConstraintElement(rsI, 1);
			rsI++;
			vec[i]->createConstraintElement(rsI, -1);
		}
		rsI++;
	}
}


void QlpConverter::relaxQlpQuantifiers(data::Qlp& qlp) {
	std::vector<data::QpVar*> vec = qlp.getVariableVectorByQuantifier(data::QpVar::all);
	for (unsigned int i = 0; i < vec.size(); i++)
		vec[i]->setQuantifier(data::QpVar::exists);
	vec = qlp.getVariableVectorByQuantifier(data::QpVar::random);
	for (unsigned int i = 0; i < vec.size(); i++)
		vec[i]->setQuantifier(data::QpVar::exists);
}

void QlpConverter::normalizeRatioSign(data::Qlp& qlp, data::QpRhs::RatioSign s) {
	std::vector<data::Constraint*> cVec = qlp.getConstraintVec();
	for (unsigned int i = 0; i < cVec.size(); i++) {
		if (cVec[i]->getRhsRatioSign() == data::QpRhs::equal)
			throw utils::AlgorithmException("Cannot normalize =-sign to >= or <=, call QlpConverter::splitEqualities(...) first");
		if (cVec[i]->getRhsRatioSign() != s)
			cVec[i]->multiplyConstraint(-1.0);
	}
}

void QlpConverter::setBounds(data::Qlp& qlp, data::QpNum lb, data::QpNum ub) {
	unsigned int size = qlp.getVariableCount();
	std::vector<data::QpVar*> vec = qlp.getVariableVector();
	for (unsigned int i = 0; i < size; i++) {
		if (vec[i]->getQuantifier() == data::QpVar::exists) {
			if (vec[i]->getUpperBound() > ub) {
				vec[i]->setUpperBound(ub);
			}
			if (vec[i]->getLowerBound() < lb) {
				vec[i]->setLowerBound(lb);
			}
		}
	}
}

void QlpConverter::setUpperBounds(data::Qlp& qlp, data::QpNum ub) {
	unsigned int size = qlp.getVariableCount();
	std::vector<data::QpVar*> vec = qlp.getVariableVector();
	for (unsigned int i = 0; i < size; i++) {
		if (vec[i]->getQuantifier() == data::QpVar::exists && vec[i]->getUpperBound() > ub)
			vec[i]->setUpperBound(ub);
	}
}

void QlpConverter::setLowerBounds(data::Qlp& qlp, data::QpNum lb) {
	unsigned int size = qlp.getVariableCount();
	std::vector<data::QpVar*> vec = qlp.getVariableVector();
	for (unsigned int i = 0; i < size; i++) {
		if (vec[i]->getQuantifier() == data::QpVar::exists && vec[i]->getLowerBound() < lb)
			vec[i]->setUpperBound(lb);
	}
}

void QlpConverter::normalizeExistentialBounds(data::Qlp& q) {

	std::vector<data::QpVar*> varVec = q.getVariableVectorByQuantifier(data::QpVar::exists);
	data::QpNum lb, coeff;
	data::Constraint* con;
	std::vector<data::Constraint*> cVec;
	data::QpNum offSet = 0;

	for (unsigned int i = 0; i < varVec.size(); i++) {

		if ((lb = varVec[i]->getLowerBound()) < 0 && !lb.isMinInf()) {
			cVec = q.getVariableConstraints(varVec[i]->getIndex());
			for (unsigned int j = 0; j < cVec.size(); j++) {
				con = cVec[j];
				if (!((coeff = con->getConstraintElementCoeff(varVec[i]->getIndex())).isMinInf())) {
					con->setRhsValue(con->getRhsValue() - (lb * coeff));
				}
			}

			offSet += (lb * q.getObjectiveFunctionElement(varVec[i]->getIndex()));
			varVec[i]->setLowerBound(0);
			varVec[i]->setUpperBound(varVec[i]->getUpperBound() - lb);
		}
	}
	if (!offSet.isZero())
		q.getObjectiveFunction().setOffset(q.getObjectiveFunction().getOffset() + offSet);
}

void QlpConverter::normalizeUniversalBounds(data::Qlp& q) {

	std::vector<data::QpVar*> varVec = q.getVariableVectorByQuantifier(data::QpVar::all);
	data::QpNum lb, offSet, coeff;
	data::Constraint* con;
	std::vector<data::Constraint*> cVec;
	data::QpNum offSetObj = 0;

	for (unsigned int i = 0; i < varVec.size(); i++) {
		if (!(lb = varVec[i]->getLowerBound()).isZero() || varVec[i]->getUpperBound() != 1.0) {
			offSet = (varVec[i]->getUpperBound() - lb);
			cVec = q.getVariableConstraints(varVec[i]->getIndex());
			for (unsigned int j = 0; j < cVec.size(); j++) {
				con = cVec[j];
				coeff = con->getConstraintElementCoeff(varVec[i]->getIndex());
				con->setRhsValue(con->getRhsValue() - (lb * coeff));
				con->deleteConstraintElement(varVec[i]->getIndex());
				con->createConstraintElement(varVec[i]->getIndex(), coeff * offSet);
			}
			offSetObj += (lb * q.getObjectiveFunctionElement(varVec[i]->getIndex()));
			varVec[i]->setLowerBound(0);
			varVec[i]->setUpperBound(1);
		}
	}
	if (!offSetObj.isZero())
		q.getObjectiveFunction().setOffset(q.getObjectiveFunction().getOffset() + offSetObj);
}

unsigned int QlpConverter::fixEqualityBounds(data::Qlp& q) {
	unsigned int fixed = 0;
	std::vector<data::QpVar*> varVec = q.getVariableVectorByQuantifier(data::QpVar::exists);
	for (unsigned int i = 0; i < varVec.size(); i++) {
		if (varVec[i]->getLowerBound() == varVec[i]->getUpperBound()) {
			fixed++;
			q.setVariableValue(varVec[i]->getIndex(), varVec[i]->getLowerBound());
		}
	}
	return fixed;
}

void QlpConverter::splitEqualities(data::Qlp& qlp) {
	std::vector<data::Constraint*> cVec = qlp.getConstraintVec();
	for (unsigned int i = 0; i < cVec.size(); i++) {
		if (cVec[i]->getRhsRatioSign() == data::QpRhs::equal) {
			data::QpNum rhsVal(cVec[i]->getRhsValue());
			cVec[i]->setRhsRatioSign(data::QpRhs::greaterThanOrEqual);
			data::Constraint* tmp = &qlp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, rhsVal, cVec[i]->getResponsibility());
			int index = cVec[i]->getFirstConstraintElementIndex();
			while (index >= 0) {
				tmp->createConstraintElement(index, cVec[i]->getConstraintElementCoeff(index));
				index = cVec[i]->getNextConstraintElementIndex(index);
			}
			cVec[i]->setRhsValue(rhsVal);
			tmp->setRhsValue(rhsVal);
		}
	}
}

void QlpConverter::splitEqualityBounds(data::Qlp& qlp) {

	std::vector<data::QpVar*> v =qlp.getVariableVector();
	for(unsigned int i = 0; i < v.size();i++){
		if(v[i]->getQuantifier()!=data::QpVar::exists)
			continue;
		if(v[i]->getUpperBound()==v[i]->getLowerBound()){
			data::QpNum rhsVal(v[i]->getUpperBound());
            qlp.createRhsConstraint(data::QpRhs::smallerThanOrEqual,rhsVal, (int)data::QpRhs::EXISTENTIAL).addConstraintElement(v[i]->getIndex(),1.0);
			qlp.createRhsConstraint(data::QpRhs::greaterThanOrEqual,rhsVal, (int)data::QpRhs::EXISTENTIAL).addConstraintElement(v[i]->getIndex(),1.0);
			v[i]->setLowerBound(rhsVal>=0 ? 0 : LB_QLP);
			v[i]->setUpperBound(UB_QLP);
		}
	}

}

void QlpConverter::removeEmptyConstraints(data::Qlp& qlp) {
	std::vector<data::Constraint*> cVec = qlp.getConstraintVec();
	for (unsigned int i = 0; i < cVec.size(); i++) {
		if (cVec[i]->getElementCount() == 0) {
			if (LOG_QLPCONV)
				utils::Logger::globalLog(utils::LOG_INFO, "QlpConverter::removeEmptyConstraints(data::Qlp& qlp)", "Removing: " + cVec[i]->toString());
			cVec[i]->removeConstraint();
		}
	}
}

void QlpConverter::relaxQlpNumberSystem(data::Qlp& qlp, bool all) {
	std::vector<data::QpVar*> vars = qlp.getVariableVector();
	unsigned int startIndex = 0;
	if (!all) {
		std::vector<const data::QpVar *> v = qlp.getVariableVectorByQuantifierConst(data::QpVar::all);
		if (v.size())
			startIndex = v[0]->getIndex();
	}
	for (unsigned int i = startIndex; i < vars.size(); i++)
		vars[i]->setNumberType(data::QpVar::real);
}

void QlpConverter::getMatrixPartByQuantifier(const data::Qlp& qlp, std::vector<std::vector<data::IndexedElement> >& m, data::QpVar::Quantifier q) {

	if (qlp.getConstraintCount() != m.size())
		throw utils::QlpSolverException("QlpConverter: can't create existential matrix since dimensions do not coincide");

	std::vector<const data::QpVar *> vars = qlp.getVariableVectorByQuantifierConst(q);
	std::vector<data::QpVar const *> vVec = qlp.getVariableVectorConst();
	std::vector<data::Constraint const *> cons = qlp.getConstraintVecConst();

	std::unordered_map<std::string, unsigned int> nameToVariableIndex;
	for (unsigned int i = 0; i < vars.size(); i++) {
		nameToVariableIndex.insert(std::pair<std::string, unsigned int>(vars[i]->getName(), i));
	}

	std::vector<data::QpNum> tmpCoeffs(qlp.getVariableCount());
	std::vector<data::IndexedElement> tmpIeCoeffs;

	data::QpNum val;
	unsigned int index = 0;
	for (unsigned int i = 0; i < cons.size(); i++) {
		m[i].clear();
		tmpIeCoeffs.clear();
		cons[i]->getElementsSparse(tmpIeCoeffs);
		for (unsigned int j = 0; j < tmpIeCoeffs.size(); j++) {
			if (vVec[(index = tmpIeCoeffs[j].index)]->getQuantifier() == q) {
				unsigned int col = nameToVariableIndex.find(vVec[index]->getName())->second;
				m[i].push_back(data::IndexedElement(col, tmpIeCoeffs[j].value));
			}
		}
	}
}

//---------------------------------------------------------------------------------- AFTER CLEANUP -------------------------------------------------------------------------------->
void QlpConverter::convertToLP(const data::Qlp& source, data::Qlp& target, DepType s, ProblemType type) {

	//Check if input is not a quantified problem
	if (!(source.getQuantifierCount(data::QpVar::random) + source.getQuantifierCount(data::QpVar::all))) {
		target = source;
		return;
	}

	//------------------------------ Initialize TMP QLP ---------------------------------------->
	data::Qlp tmp;
	QlpConverter::initDepConversion(source, tmp, type);

	//-------------------------------- Initialize Splitter ------------------------------------->
	data::QpMatrix<algorithm::DepTreeNode>::Type qpTree;
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Splitter...");
	utils::QlpSplitter qlpSplitter(tmp);
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing Splitter...");
	if (type == WORST) {
		qlpSplitter.initSplitter(algorithm::Algorithm::WORST_CASE);
	} else {
		qlpSplitter.initSplitter(algorithm::Algorithm::AVERAGE_CASE);
	}
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing ScenarioTree...");
	qlpSplitter.initializeDepTree(qpTree);

	//--------------------------------- Create DEP ---------------------------------------------->
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating DEP...");
	if (s == COMPACT_VIEW) {
		convertToCompactViewDEP(tmp, target, type, qlpSplitter, qpTree);
	} else {
		convertToSplitVariableDEP(tmp, target, type, qlpSplitter, qpTree);
	}
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished.");
}

void QlpConverter::convertToLP(const data::Qlp& source, data::QpObjFunc& obj, std::vector<data::QpVar>& vars, data::QpSparseMatrix& matrix, std::vector<data::QpRhs>& rhs, DepType s, ProblemType type) {

	//Check if input is not a quantified problem
	if (!(source.getQuantifierCount(data::QpVar::random) + source.getQuantifierCount(data::QpVar::all))) {
		source.getQlpParts(obj, vars, matrix, rhs);
		return;
	}

	//------------------------------ Initialize TMP QLP ---------------------------------------->
	data::Qlp tmp;
	QlpConverter::initDepConversion(source, tmp, type);

	//-------------------------------- Initialize Splitter ------------------------------------->
	data::QpMatrix<algorithm::DepTreeNode>::Type qpTree;
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Splitter...");
	utils::QlpSplitter qlpSplitter(tmp);
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing Splitter...");
	if (type == WORST) {
		qlpSplitter.initSplitter(algorithm::Algorithm::WORST_CASE);
	} else {
		qlpSplitter.initSplitter(algorithm::Algorithm::AVERAGE_CASE);
	}
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Initializing ScenarioTree...");
	qlpSplitter.initializeDepTree(qpTree);

	//--------------------------------- Create DEP --------------------------------------------->
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating DEP...");
	if (s == COMPACT_VIEW) {
		convertToCompactViewDEP(tmp, obj, vars, matrix, rhs, type, qlpSplitter, qpTree);
	} else {
		convertToSplitVariableDEP(tmp, obj, vars, matrix, rhs, type, qlpSplitter, qpTree);
	}
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Finished.");
}

void QlpConverter::initDepConversion(const data::Qlp& source, data::Qlp& tmp, ProblemType type) {

	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Pushing objective function to QLP...");

	if (type == utils::QlpConverter::WORST) {
		if (source.getVariableIndexByName("QpDepObjApprox") == -1)
			QlpConverter::pushObjectiveFunctionToMatrixFront(source, tmp);
	} else {
		if (source.getVariableIndexByName("QpDepObjApprox") == -1)
			QlpConverter::pushObjectiveFunctionToMatrixBack(source, tmp);
	}

	if (!source.getObjectiveFunction().getObjectiveElementsDense().size()) {
		tmp.getVariableByName("QpDepObjApprox").setBounds(0, 0);
	}

	if (tmp.getVariableVectorConst()[0]->getQuantifier() != data::QpVar::exists) {
		QlpConverter::addDummyVariable(tmp, true);
	}

	if (tmp.getVariableVectorConst()[tmp.getVariableVectorConst().size() - 1]->getQuantifier() != data::QpVar::exists) {
		QlpConverter::addDummyVariable(tmp, false);
	}

}

void QlpConverter::pushObjectiveFunctionToMatrixFront(const data::Qlp& s, data::Qlp& target) {

	data::Qlp tmp(s);

	//Create Objective Function approximation variable
	int index = 0;
	target.clear();
	target.createVariable("QpDepObjApprox", index, data::QpVar::exists).setBounds(APPROX_LOWER_BOUND, APPROX_UPPER_BOUND);
	target.setObjective(tmp.getObjective());
	target.getObjectiveFunction().setSize(tmp.getVariableCount());
	target.getObjectiveFunction().setOffset(tmp.getObjectiveFunction().getOffset());
	//Copy all variables with indices being shifted by +1
	data::QpVar const * var;
	std::vector<data::QpVar const *> varVec = tmp.getVariableVectorConst();
	data::QpNum coeff;
	for (unsigned int i = 0, index = 1; i < varVec.size(); i++, index++) {
		var = varVec[i];
		if (var->getQuantifier() != data::QpVar::random) {
			target.createVariable(tmp.getVariableNameByIndex(var->getIndex()), index, var->getQuantifier(), var->getNumberSystem(), var->getLowerBound(), var->getUpperBound());
		} else {
			target.createVariable(tmp.getVariableNameByIndex(var->getIndex()), index, var->getQuantifier(), var->getNumberSystem(), var->getVariableRange(), var->getVariableDistribution());
		}
		if (!(coeff = tmp.getObjectiveFunctionElement(i)).isZero())
			target.setObjectiveFunctionElement(index, coeff);
	}
	//Copy constraints with all indices being shifted by +1
	data::Constraint *created;
	std::vector<data::Constraint const *> cVec = tmp.getConstraintVecConst();
	std::vector<data::IndexedElement> coeffVec;
	for (unsigned int i = 0; i < cVec.size(); i++) {
		cVec[i]->getElementsSparse(coeffVec);
		created = &target.createRhsConstraint(cVec[i]->getRhsRatioSign(), cVec[i]->getRhsValue(), cVec[i]->getResponsibility());
		for (unsigned int j = 0; j < coeffVec.size(); j++) {
			created->createConstraintElement(coeffVec[j].index + 1, coeffVec[j].value);
		}
	}
	//Create Objective Function Constraints and add corresponding k Variable index to Objective Function
	data::Constraint& newCon = target.createRhsConstraint((target.getObjective() == data::QpObjFunc::min) ? data::QpRhs::smallerThanOrEqual : data::QpRhs::greaterThanOrEqual, target.getObjectiveFunction().getOffset() * (-1.0), (int)data::QpRhs::EXISTENTIAL);
	std::vector<data::QpNum> vec = target.getObjectiveFunctionValues();
	for (unsigned int i = 0; i < vec.size(); i++) {
		if (!vec[i].isZero())
			newCon.createConstraintElement(i, vec[i]);
	}
	newCon.createConstraintElement(target.getVariableIndexByName("QpDepObjApprox"), -1.0);
	target.getObjectiveFunction().clearElements();
	target.setObjectiveFunctionElement(0, 1.0);

}

void QlpConverter::pushObjectiveFunctionToMatrixBack(const data::Qlp& qlp, data::Qlp& target) {
	target = qlp;
	target.createVariable("QpDepObjApprox", qlp.getVariableCount(), data::QpVar::exists).setBounds(APPROX_LOWER_BOUND, APPROX_UPPER_BOUND);
	data::Constraint& newCon = target.createRhsConstraint((qlp.getObjective() == data::QpObjFunc::min) ? data::QpRhs::smallerThanOrEqual : data::QpRhs::greaterThanOrEqual, qlp.getObjectiveFunction().getOffset() * -1.0, (int)data::QpRhs::EXISTENTIAL);
	std::vector<data::QpNum> objVec = qlp.getObjectiveFunctionValues();
	for (unsigned j = 0; j < objVec.size(); j++) {
		if (!objVec[j].isZero()) {
			newCon.createConstraintElement(j, objVec[j]);
		}
	}
	newCon.createConstraintElement(target.getVariableIndexByName("QpDepObjApprox"), -1.0);
	target.getObjectiveFunction().clearElements();
	target.getObjectiveFunction().setOffset(0.0);
	target.getObjectiveFunction().setSize(target.getVariableCount());
	target.setObjectiveFunctionElement(target.getVariableCount() - 1, 1.0);
}

void QlpConverter::addDummyVariable(data::Qlp& qlp, bool front) {
	if (front) {

		if (qlp.getVariableIndexByName("DepDummyFront") != -1)
			throw utils::AlgorithmException("QlpConverter::addDummyVariable(data::Qlp& qlp, bool front = true) --> Variable DepDummyFront already exists");

		data::Qlp tmp;

		//Create Objective Function approximation variable
		int index = 0;

		tmp.createVariable("DepDummyFront", index, data::QpVar::exists).setBounds(APPROX_LOWER_BOUND, APPROX_UPPER_BOUND);
		tmp.setObjective(qlp.getObjective());
		tmp.getObjectiveFunction().setOffset(qlp.getObjectiveFunction().getOffset());

		//Copy all variables with indices being shifted by +1
		data::QpVar const * var;
		std::vector<data::QpVar const *> varVec = qlp.getVariableVectorConst();
		data::QpNum coeff;
		for (unsigned int i = 0, index = 1; i < varVec.size(); i++, index++) {
			var = varVec[i];
			if (var->getQuantifier() != data::QpVar::random) {
				tmp.createVariable(qlp.getVariableNameByIndex(var->getIndex()), index, var->getQuantifier(), var->getNumberSystem(), var->getLowerBound(), var->getUpperBound());
			} else {
				tmp.createVariable(qlp.getVariableNameByIndex(var->getIndex()), index, var->getQuantifier(), var->getNumberSystem(), var->getVariableRange(), var->getVariableDistribution());
			}
			if (!(coeff = qlp.getObjectiveFunctionElement(i)).isZero())
				tmp.setObjectiveFunctionElement(index, coeff);
		}

		//Copy constraints with all indices being shifted by +1
		data::Constraint *created;
		std::vector<data::Constraint const *> cVec = qlp.getConstraintVecConst();
		std::vector<data::IndexedElement> coeffVec;
		for (unsigned int i = 0; i < cVec.size(); i++) {
			cVec[i]->getElementsSparse(coeffVec);
			created = &tmp.createRhsConstraint(cVec[i]->getRhsRatioSign(), cVec[i]->getRhsValue(), cVec[i]->getResponsibility());
			for (unsigned int j = 0; j < coeffVec.size(); j++) {
				created->createConstraintElement(coeffVec[j].index + 1, coeffVec[j].value);
			}
		}

	} else {
		if (qlp.getVariableIndexByName("DepDummyFront") != -1)
			throw utils::AlgorithmException("QlpConverter::addDummyVariable(data::Qlp& qlp, bool front = false) --> Variable DepDummyBack already exists");
		qlp.createVariable("DepDummyBack", qlp.getVariableCount(), data::QpVar::exists).setBounds(0, 0);
	}
}

void QlpConverter::convertToCompactViewDEP(const data::Qlp& source, data::Qlp& target, ProblemType type, utils::QlpSplitter& qlpSplitter, data::QpMatrix<algorithm::DepTreeNode>::Type& qpTree) {

	std::string s, si;
	target.clear();

	//--------------- Create Variables and set Objective Function ------------------------------>
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Variables...");
	unsigned int index = 0;
	for (unsigned int i = 0; i < qlpSplitter.stages; i++) {
		for (unsigned int j = 0; j < qpTree[i].size(); j++) {
			for (unsigned int from = qlpSplitter.existVarIndAtDepth[i].first, to = qlpSplitter.existVarIndAtDepth[i].second; from <= to; from++) {
				const data::QpVar& v = source.getVariableByIndexConst(from);
				s.assign(v.getName());
				s += "_" + utils::ToolBox::convertToString(j);
				target.createVariable(s, index++, v.getQuantifier(), v.getNumberSystem(), v.getLowerBound(), v.getUpperBound());
			}
		}
	}
	target.setObjective(source.getObjective());
	target.setObjFuncOffset(source.getObjFuncOffset());
	target.getObjectiveFunction().setSize(target.getVariableCount());
	if (type == WORST) {
		target.setObjectiveFunctionElement(0, 1);
	} else {
		data::QpRationalMatrix distr, cartProdDistr;
		std::vector<data::QpRational> probVec;
		for (unsigned int i = 0; i < qlpSplitter.randVars.size(); i++) {
			distr.push_back(qlpSplitter.randVars[i].getVariableDistribution());
		}
		utils::QlpConverter::createCartesianProduct(distr, cartProdDistr);
		for (unsigned int i = 0, size = cartProdDistr.size(); i < size; i++)
			probVec.push_back(data::QpNum::prodQpRationalVec(cartProdDistr[i]));
		for (unsigned int i = 0; i < qpTree[qlpSplitter.stages - 1].size(); i++) {
			s = "QpDepObjApprox_";
			s += utils::ToolBox::convertToString(i);
			target.setObjectiveFunctionElement(target.getVariableIndexByName(s), probVec[i].asDouble());
		}
	}

	//-------------------------------------- Creating constraints ---------------------------------------------------->
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Constraints...");
	std::vector<data::Constraint const *> vec = source.getConstraintVecConst();
	std::vector<data::IndexedElement> ieVec;
	for (unsigned int i = 0; i < qlpSplitter.stages; i++) {
		if (qlpSplitter.conIndAtDepth[i].first < 0)
			continue; //no single constraints at this stage
		for (unsigned int j = 0; j < qpTree[i].size(); j++) {
			std::vector<unsigned int> currExitVarScenIndex(qlpSplitter.stages, 0);
			utils::QlpSplitter::createScenIndexVector(qpTree[i][j], qpTree, currExitVarScenIndex);
			for (int k = qlpSplitter.conIndAtDepth[i].first, m = 0; k <= qlpSplitter.conIndAtDepth[i].second; k++, m++) {
				data::Constraint& tmpCon = target.createRhsConstraint(vec[k]->getRhsRatioSign(), vec[k]->getRhsValue() - qpTree[i][j].scenarioRhs->operator [](m), vec[k]->getResponsibility());
				vec[k]->getElementsSparse(ieVec);
				for (unsigned int l = 0; l < ieVec.size(); l++) {
					const data::QpVar& v = source.getVariableByIndexConst(ieVec[l].index);
					if (v.getQuantifier() != data::QpVar::exists)
						continue;
					s.assign(v.getName()/*source.getVariableNameByIndex(ieVec[l].index)*/);
					for (int m = 0; m < qlpSplitter.existVarIndAtDepth.size(); m++) {
						if (ieVec[l].index <= qlpSplitter.existVarIndAtDepth[m].second) {
							si.assign(s);
							si += "_" + utils::ToolBox::convertToString(currExitVarScenIndex[m]);
							tmpCon.createConstraintElement(target.getVariableIndexByName(si), ieVec[l].value);
							break;
						}
					}
				}
			}
		}
	}
}

void QlpConverter::convertToCompactViewDEP(const data::Qlp& source, data::QpObjFunc& obj, std::vector<data::QpVar>& vars, data::QpSparseMatrix& matrix, std::vector<data::QpRhs>& rhs, ProblemType type, utils::QlpSplitter& qlpSplitter, data::QpMatrix<algorithm::DepTreeNode>::Type& qpTree) {

	obj.clear();
	vars.clear();
	matrix.clear();
	rhs.clear();

	//--------------- Create Variables and set Objective Function ------------------------------>
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Variables ...");
	std::unordered_map<std::string, unsigned int> nameToVariableIndex;
	std::vector<data::QpVar const *> varVec = source.getVariableVectorConst();
	std::string s;
	unsigned int index = 0;
	for (unsigned int i = 0; i < qlpSplitter.stages; i++) {
		for (unsigned int j = 0; j < qpTree[i].size(); j++) {
			for (unsigned int from = qlpSplitter.existVarIndAtDepth[i].first, to = qlpSplitter.existVarIndAtDepth[i].second; from <= to; from++) {
				data::QpVar const * v = varVec[from];
				s.assign(v->getName());
				s += "_" + utils::ToolBox::convertToString(j);
				nameToVariableIndex.insert(std::pair<std::string, int>(s, index));
				vars.push_back(data::QpVar(s, index, v->getLowerBound(), v->getUpperBound(), v->getNumberSystem(), v->getQuantifier()));
				index++;
			}
		}
	}

	obj.setObjective(source.getObjective());
	obj.setOffset(source.getObjFuncOffset());
	obj.setSize(vars.size());
	if (type == WORST) {
		obj.setObjElement(0, 1);
	} else {
		data::QpRationalMatrix distr, cartProdDistr;
		std::vector<data::QpRational> probVec;
		for (unsigned int i = 0; i < qlpSplitter.randVars.size(); i++) {
			distr.push_back(qlpSplitter.randVars[i].getVariableDistribution());
		}
		utils::QlpConverter::createCartesianProduct(distr, cartProdDistr);
		for (unsigned int i = 0, size = cartProdDistr.size(); i < size; i++)
			probVec.push_back(data::QpNum::prodQpRationalVec(cartProdDistr[i]));

		for (unsigned int i = 0; i < qpTree[qlpSplitter.stages - 1].size(); i++) {
			s.assign("QpDepObjApprox_" + utils::ToolBox::convertToString(i));
			obj.setObjElement(nameToVariableIndex.find(s)->second, probVec[i].asDouble());
		}
	}

	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Constraints...");
	std::vector<data::Constraint const *> vec = source.getConstraintVecConst();
	std::vector<data::IndexedElement> ieVec;
	data::IndexedElement ie;
	std::string s1;
	unsigned int rhsIndex = 0;
	for (unsigned int i = 0; i < qlpSplitter.stages; i++) {
		if (qlpSplitter.conIndAtDepth[i].first < 0)
			continue;
		for (unsigned int j = 0; j < qpTree[i].size(); j++) {
			std::vector<unsigned int> currExitVarScenIndex(qlpSplitter.stages, 0);
			utils::QlpSplitter::createScenIndexVector(qpTree[i][j], qpTree, currExitVarScenIndex);
			for (int k = qlpSplitter.conIndAtDepth[i].first, m = 0; k <= qlpSplitter.conIndAtDepth[i].second; k++, m++) {
				rhs.push_back(data::QpRhs(vec[k]->getRhsValue() - qpTree[i][j].scenarioRhs->operator [](m), vec[k]->getRhsRatioSign()));
				vec[k]->getElementsSparse(ieVec);
				matrix.push_back(std::vector<data::IndexedElement>());
				for (int l = 0; l < ieVec.size(); l++) {
					if (varVec[(ie = ieVec[l]).index]->getQuantifier() != data::QpVar::exists)
						continue;
					s.assign(varVec[ie.index]->getName()/*source.getVariableNameByIndex(ie.index)*/);
					for (int m = 0; m < qlpSplitter.existVarIndAtDepth.size(); m++) {
						if (ie.index <= qlpSplitter.existVarIndAtDepth[m].second) {
							s1.assign(s);
							s1 += "_" + utils::ToolBox::convertToString(currExitVarScenIndex[m]);
							matrix[rhsIndex].push_back(data::IndexedElement(nameToVariableIndex.find(s1)->second, ie.value));
							break;
						}
					}
				}
				rhsIndex++;
			}
		}
	}
}

void QlpConverter::convertToSplitVariableDEP(const data::Qlp& source, data::Qlp& target, ProblemType type, utils::QlpSplitter& qlpSplitter, data::QpMatrix<algorithm::DepTreeNode>::Type& qpTree) {

	target.clear();
	std::string s, s0, si;

	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Variables ...");
	unsigned int index = 0;
	for (unsigned int i = 0; i < qpTree[qlpSplitter.stages - 1].size(); i++) {
		for (unsigned int j = 0, size = qlpSplitter.existVars.size(); j < size; j++) {
			data::QpVar& v = qlpSplitter.existVars[j];
			s.assign(v.getName());
			s += "_" + utils::ToolBox::convertToString(i);
			target.createVariable(s, index++, v.getQuantifier(), v.getNumberSystem(), v.getLowerBound(), v.getUpperBound());
		}
	}

	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Obj ...");
	target.setObjective(source.getObjective());
	target.setObjFuncOffset(source.getObjFuncOffset());
	target.getObjectiveFunction().setSize(target.getVariableCount());
	if (type == WORST) {
		target.setObjectiveFunctionElement(0, 1);
	} else {
		data::QpRationalMatrix distr, cartProdDistr;
		std::vector<data::QpRational> probVec;
		for (unsigned int i = 0; i < qlpSplitter.randVars.size(); i++) {
			distr.push_back(qlpSplitter.randVars[i].getVariableDistribution());
		}
		utils::QlpConverter::createCartesianProduct(distr, cartProdDistr);
		for (unsigned int i = 0, size = cartProdDistr.size(); i < size; i++)
			probVec.push_back(data::QpNum::prodQpRationalVec(cartProdDistr[i]));
		for (unsigned int i = 0; i < qpTree[qlpSplitter.stages - 1].size(); i++) {
			std::string s("QpDepObjApprox_" + utils::ToolBox::convertToString(i));
			target.setObjectiveFunctionElement(target.getVariableIndexByName(s), probVec[i].asDouble());
		}
	}

	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Constraints ...");

	std::vector<data::Constraint const *> vec = source.getConstraintVecConst();
	std::vector<data::QpNum> scenRhs(vec.size(), 0);
	std::vector<data::IndexedElement> ieVec;
	for (unsigned int i = 0; i < qpTree[qlpSplitter.stages - 1].size(); i++) {

		computeCurrentScenRhs(scenRhs, qpTree[qlpSplitter.stages - 1][i]);

		for (unsigned int j = 0; j < vec.size(); j++) {

			data::Constraint & tmpCon = target.createRhsConstraint(vec[j]->getRhsRatioSign(), vec[j]->getRhsValue() - scenRhs[j], vec[j]->getResponsibility());

			ieVec = vec[j]->getElements();

			for (unsigned int l = 0; l < ieVec.size(); l++) {
				if (source.getVariableByIndexConst(ieVec[l].index).getQuantifier() != data::QpVar::exists)
					continue; //only add existential variables
				s.assign(source.getVariableNameByIndex(ieVec[l].index));
				s += "_" + utils::ToolBox::convertToString(i);
				tmpCon.createConstraintElement(target.getVariableIndexByName(s), ieVec[l].value);
			}
		}
	}

	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating NonAnticipativity ...");
	for (unsigned int i = 0; i < qlpSplitter.stages - 1; i++) {
		unsigned int mult = 0;

		for (unsigned j = 0; j < qpTree[i].size(); j++) {
			mult = qpTree[qpTree.size() - 1].size() / qpTree[i].size();
			for (unsigned int from = (i == 0 && type == AVG) ? qlpSplitter.existMatrixIndexAtDepth[i].first + 1 : qlpSplitter.existMatrixIndexAtDepth[i].first, to = qlpSplitter.existMatrixIndexAtDepth[i].second; from <= to; from++) {
				s0.assign(qlpSplitter.existVars[from].getName());
				s0 += "_" + utils::ToolBox::convertToString(j * mult);
				for (unsigned int k = 1; k < mult; k++) {
					si.assign(qlpSplitter.existVars[from].getName());
					si += "_" + utils::ToolBox::convertToString(j * mult + k);
					data::Constraint& tmpCon = target.createRhsConstraint(data::QpRhs::equal, 0, (int)data::QpRhs::EXISTENTIAL);
					tmpCon.createConstraintElement(target.getVariableIndexByName(s0), 1);
					tmpCon.createConstraintElement(target.getVariableIndexByName(si), -1);
				}
			}
		}
	}
}

void QlpConverter::convertToSplitVariableDEP(const data::Qlp&source, data::QpObjFunc& obj, std::vector<data::QpVar>& vars, data::QpSparseMatrix& matrix, std::vector<data::QpRhs>& rhs, ProblemType type, utils::QlpSplitter& qlpSplitter, data::QpMatrix<algorithm::DepTreeNode>::Type& qpTree) {

	obj.clear();
	vars.clear();
	matrix.clear();
	rhs.clear();

	//--------------- Create Variables and set Objective Function ------------------------------>
	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Variables ...");
	std::unordered_map<std::string, unsigned int> nameToVariableIndex;
	std::string s, s0, si;
	unsigned int index = 0;
	for (unsigned int i = 0; i < qpTree[qlpSplitter.stages - 1].size(); i++) {
		for (unsigned int j = 0, size = qlpSplitter.existVars.size(); j < size; j++) {
			data::QpVar& v = qlpSplitter.existVars[j];
			s.assign(v.getName());
			s += "_" + utils::ToolBox::convertToString(i);
			nameToVariableIndex.insert(std::pair<std::string, int>(s, index));
			vars.push_back(data::QpVar(s, index, v.getLowerBound(), v.getUpperBound(), v.getNumberSystem(), v.getQuantifier()));
			index++;
		}
	}

	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Obj ...");
	obj.setObjective(source.getObjective());
	obj.setOffset(source.getObjFuncOffset());
	obj.setSize(vars.size());
	if (type == WORST) {
		obj.setObjElement(0, 1);
	} else {
		data::QpRationalMatrix distr, cartProdDistr;
		std::vector<data::QpRational> probVec;
		for (unsigned int i = 0; i < qlpSplitter.randVars.size(); i++) {
			distr.push_back(qlpSplitter.randVars[i].getVariableDistribution());
		}
		utils::QlpConverter::createCartesianProduct(distr, cartProdDistr);
		for (unsigned int i = 0, size = cartProdDistr.size(); i < size; i++)
			probVec.push_back(data::QpNum::prodQpRationalVec(cartProdDistr[i]));
		for (unsigned int i = 0; i < qpTree[qlpSplitter.stages - 1].size(); i++) {
			std::string s("QpDepObjApprox_" + utils::ToolBox::convertToString(i));
			obj.setObjElement(nameToVariableIndex.find(s)->second, probVec[i].asDouble());
		}
	}

	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating Constraints ...");

	std::vector<data::Constraint const *> vec = source.getConstraintVecConst();
	std::vector<data::QpNum> scenRhs(vec.size(), 0);
	std::vector<data::IndexedElement> ieVec;

	for (unsigned int i = 0; i < qpTree[qlpSplitter.stages - 1].size(); i++) {
		computeCurrentScenRhs(scenRhs, qpTree[qlpSplitter.stages - 1][i]);
		for (unsigned int j = 0; j < vec.size(); j++) {
			rhs.push_back(data::QpRhs(vec[j]->getRhsValue() - scenRhs[j], vec[j]->getRhsRatioSign()));
			ieVec = vec[j]->getElements();
			matrix.push_back(std::vector<data::IndexedElement>());
			for (unsigned int l = 0; l < ieVec.size(); l++) {
				if (source.getVariableByIndexConst(ieVec[l].index).getQuantifier() != data::QpVar::exists)
					continue;
				s.assign(source.getVariableNameByIndex(ieVec[l].index));
				s += "_" + utils::ToolBox::convertToString(i);
				matrix[matrix.size() - 1].push_back(data::IndexedElement(nameToVariableIndex.find(s)->second, ieVec[l].value));
			}
		}
	}

	if (LOG_QLPCONV)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Creating NonAnticipativity ...");
	for (unsigned int i = 0; i < qlpSplitter.stages - 1; i++) {
		unsigned int mult = 0;
		for (unsigned j = 0; j < qpTree[i].size(); j++) {
			mult = qpTree[qpTree.size() - 1].size() / qpTree[i].size();
			for (unsigned int from = (i == 0 && type == AVG) ? qlpSplitter.existMatrixIndexAtDepth[i].first + 1 : qlpSplitter.existMatrixIndexAtDepth[i].first, to = qlpSplitter.existMatrixIndexAtDepth[i].second; from <= to; from++) {
				s0.assign(qlpSplitter.existVars[from].getName());
				s0 += "_" + utils::ToolBox::convertToString(j * mult);
				for (unsigned int k = 1; k < mult; k++) {
					si.assign(qlpSplitter.existVars[from].getName());
					si += "_" + utils::ToolBox::convertToString(j * mult + k);
					rhs.push_back(data::QpRhs(0, data::QpRhs::equal));
					matrix.push_back(std::vector<data::IndexedElement>());
					matrix[matrix.size() - 1].push_back(data::IndexedElement(nameToVariableIndex.find(s0)->second, 1));
					matrix[matrix.size() - 1].push_back(data::IndexedElement(nameToVariableIndex.find(si)->second, -1));
				}
			}
		}
	}
}

}
