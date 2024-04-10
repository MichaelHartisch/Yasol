/*
*
* Solver: Qlp.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QUANTIFIEDPROBLEM_HPP_
#define QUANTIFIEDPROBLEM_HPP_

//#include "Settings.hpp"
#include "Constraint.hpp"
#include "Datastructures/Datastructures.hpp"
#include "Datastructures/numbers/Numbers.hpp"
#include <map>

namespace data {

class Qlp {

	friend class QlpMatCoeff;
	friend class Constraint;

public:

	//Type of problem
	typedef enum {
		LP, IP, MIP, QLP, QIP, QIP_BIN, QMIP, QMIP_HDE, QMIP_HDA, RQLP, RQIP, RQIP_BIN,
		RQMIP, RQMIP_HDE, RQMIP_HDA, QLP_TYPE_ERROR
	} QlpType;

	//----------------- Constructors and Destructor -------------------------------->
	Qlp();

	Qlp(const data::QpObjFunc& objective,
			const std::vector<data::QpVar>& varVec,
			const data::QpSparseMatrix&,
			const std::vector<data::QpRhs>& rhsVec);
	~Qlp();
	//----------------- Copy Constructor and Assignment Opperator------------------->
	Qlp(const Qlp& rhs);
	Qlp& operator=(const Qlp&);
	Qlp& initWithQlpParts(const data::QpObjFunc&,
			const std::vector<data::QpVar>&,
			const data::QpSparseMatrix&,
			const std::vector<data::QpRhs>&);
	//----------------- Comparison Operators --------------------------------------->
	bool operator==(const Qlp&) const;
	bool operator!=(const Qlp&) const;

	//Clear the entire QLP
	void clear();
	//Get several information about the qlp
	unsigned int getNumberSystemByQuantifier(data::QpVar::Quantifier,
			data::QpVar::NumberSystem) const;
	unsigned int getQuantifierCount(data::QpVar::Quantifier) const;
	unsigned int getVariableCount() const;
	unsigned int getConstraintCount() const;
	unsigned int getMatrixElementCount() const;
	unsigned int getStageCount() const;
	unsigned int getScenarioCount() const;

	//Access single components of the Qlp
	//TODO remove this
	data::QpObjFunc& getObjectiveFunction();

	//---------------Methods to work with the QpObjFunc of the QLP --------->
	const data::QpObjFunc& getObjectiveFunction() const;

	const data::QpNum& getObjFuncOffset() const;

	void setObjFuncOffset(const data::QpNum&);

	void reverseObjFunc();

	QpObjFunc::Objective getObjective() const;

	void setObjective(QpObjFunc::Objective);

	const data::QpNum& getObjectiveFunctionElement(unsigned int) const;

	void setObjectiveFunctionElement(unsigned int, const data::QpNum&);

	const std::vector<QpNum>& getObjectiveFunctionValues() const;

	data::QpNum getObjectiveFunctionValue(const std::vector<data::QpNum>&) const;

	//---------------Methods to work with the QuantifiedVariables of the QLP --------->

	QpVar& getVariableByIndex(unsigned int);
	const QpVar& getVariableByIndexConst(unsigned int) const;

	QpVar& getVariableByName(const std::string&);
	const QpVar& getVariableByNameConst(const std::string&) const;

	std::vector<data::QpVar> getQuantifiedVariables() const;

	QpVar& createVariable(const data::QpVar&);

	QpVar& createVariable(const std::string&, const int, QpVar::Quantifier =
			QpVar::exists);

	QpVar& createVariable(const std::string&, const int, QpVar::Quantifier,
			QpVar::NumberSystem, const data::QpNum&, const data::QpNum&);

	QpVar& createVariable(const std::string&, const int, QpVar::Quantifier,
			QpVar::NumberSystem, const std::vector<data::QpNum>&,
			const std::vector<data::QpRational>&);

	std::vector<data::QpVar *> getVariableVector();
	std::vector<const data::QpVar *> getVariableVectorConst() const;

	std::vector<data::QpVar *> getVariableVectorByQuantifier(
			data::QpVar::Quantifier);
	std::vector<const data::QpVar *> getVariableVectorByQuantifierConst(
			data::QpVar::Quantifier) const;

	std::vector<Constraint*> getVariableConstraints(unsigned int index);
	std::vector<Constraint const *> getVariableConstraintsConst(
			unsigned int index) const;

	//---------------Methods to work with the Constraints of the QLP --------->
	Constraint& createRhsConstraint(QpRhs::RatioSign, const data::QpNum&, int resp);
	Constraint& createRhsConstraint(const data::QpRhs& rhs);
	Constraint& createConstraint(const data::Constraint& c);

	std::list<Constraint*> getConstraints() const;
	std::vector<Constraint*> getConstraintVec();
	std::vector<const Constraint *> getConstraintVecConst() const;

	Constraint* getFirstConstraint() const;
	Constraint* getLastConstraint() const;

	std::vector<data::QpNum> getRhsValVec() const;
	std::vector<data::QpRhs> getRhsVec() const;
	std::vector<const data::QpNum *> getRhsValVecConst() const;
	std::vector<const data::QpRhs *> getRhsVecConst() const;

	//---------------------- Various other methods ---------------------->
	void getCoeffMatrix(std::vector<std::vector<data::IndexedElement> >&m) const;


	// NEW FOR ALL SYSTEM
	void getCoeffMatrixByResp(std::vector<std::vector<data::IndexedElement> >&m, data::QpRhs::Responsibility resp) const;
	std::vector<data::QpRhs> getRhsVecByResp(data::QpRhs::Responsibility resp) const;


	bool isQuantifiedProblem() const;
	bool containsNumberSystem(data::QpVar::NumberSystem) const;
	bool isVariableUnused(unsigned int) const;
	unsigned int getVariableConstraintCount(unsigned int) const;
	unsigned int getNextMoveVariable() const;

	void eliminateVariables(const std::vector<data::IndexedElement>&);
	void eliminateVariable(unsigned int, const data::QpNum&);
	void setVariableValue(unsigned int, const data::QpNum&);

	void normalizeConstraints();
	void sortQlp();

	const std::string& getVariableNameByIndex(unsigned int) const;
	int getVariableIndexByName(const std::string&) const;

	std::string toQlpFileString(bool lp) const;
	std::string toString() const;
	std::string getQlpInfo() const;
	std::string getVariableStringByIndex(unsigned int index) const;

	void printVariables() const;

	void printColumnRowPointers();

	void getSortedVariableConstraints(unsigned int index, std::vector<
			data::Constraint*>&, std::vector<data::Constraint*>&, std::vector<
			data::Constraint*>&, std::vector<data::Constraint*>&, std::vector<
			data::Constraint*>&);

	QlpType getQlpType()const;
	static std::string qlpTypeToString(const QlpType&);

	//TODO Method to get single parts
	void getQlpParts(data::QpObjFunc& obj, std::vector<data::QpVar>& vars, data::QpSparseMatrix& matrix, std::vector<data::QpRhs>& rhs)const;

	void deleteAllRows();

protected:

	inline void setFirstRowPtr(Constraint* fr) {
		firstRow = fr;
	}

	inline Constraint* getFirstRowPtr() {
		return firstRow;
	}

	inline void setLastRowPtr(Constraint* lr) {
		lastRow = lr;
	}

	inline Constraint* getLastRowPtr() {
		return lastRow;
	}


	inline void setFirstColumnsRowPtr(QlpMatCoeff& p) {
			firstColumnsRowPointer[p.columnIndex] = &p;
	}


//	inline void setFirstColumnsRowPtr(unsigned int index, Constraint* row) {
//		firstColumnsRowPointer[index] = row;
//	}

	inline QlpMatCoeff* getFirstColumnsPointPtr(unsigned int index) const {
			if (index > firstColumnsRowPointer.size())
				return NULL;
			return firstColumnsRowPointer[index];
	}

	inline Constraint* getFirstColumnsRowPtr(unsigned int index) const {
		if (index > firstColumnsRowPointer.size())
			return NULL;
		if(!firstColumnsRowPointer[index])
		return NULL;

		return firstColumnsRowPointer[index]->pRow;
	}



//	inline Constraint* getFirstColumnsRowPtr(unsigned int index) const {
//		if (index > firstColumnsRowPointer.size())
//			return NULL;
//		return firstColumnsRowPointer[index];
//	}

	inline void setLastColumnsRowPtr(QlpMatCoeff& p) {
		if (firstColumnsRowPointer[p.columnIndex] == NULL)
			firstColumnsRowPointer[p.columnIndex] = &p;
		lastColumnsRowPointer[p.columnIndex] = &p;
	}

//	inline void setLastColumnsRowPtr(unsigned int index, Constraint* row) {
//		if (firstColumnsRowPointer[index] == NULL)
//			firstColumnsRowPointer[index] = row;
//		lastColumnsRowPointer[index] = row;
//	}

	inline QlpMatCoeff* getLastColumnsPointPtr(unsigned int index) const {
		if (index > lastColumnsRowPointer.size())
			return NULL;
		return lastColumnsRowPointer[index];
	}

	inline Constraint* getLastColumnsRowPtr(unsigned int index) const {
		if (index > lastColumnsRowPointer.size())
			return NULL;
		if(!lastColumnsRowPointer[index])
			return NULL;
		return lastColumnsRowPointer[index]->pRow;
	}

//	inline Constraint* getLastColumnsRowPtr(unsigned int index) const {
//		if (index > lastColumnsRowPointer.size())
//			return NULL;
//		return lastColumnsRowPointer[index];
//	}

	void deleteAllColumns();
	void deleteColumnByIndex(unsigned int);

	void copyQlpContent(const data::Qlp&);

	void copyQlpContent(const data::QpObjFunc&,
			const std::vector<data::QpVar>&, const data::QpSparseMatrix&,
			const std::vector<data::QpRhs>&);

	inline data::QpVar& getVar(unsigned int index)const{
		if ((unsigned) index >= variables.size()) {
			throw utils::QlpSolverException(
					"QuantifiedProblem::getVar(unsigned int) --> Index Exception i: "
							+ utils::ToolBox::convertToString(index));
		}
		return *this->variables[index];
	}

	inline const data::QpVar& getVarConst(unsigned int index) const {
		return getVar(index);
	}

private:

	// Log string for debug output
	static std::string LOG_TAG;

	// The objective function of the QLP
	data::QpObjFunc objFunc;
	// The list of variables
	std::vector<QpVar*> variables;
	// Each Column has a Pointer to the first Row that contains a QlpMatCoeff with the columns index
	// and whoose value is not equal to zero, if such a point does not exist, the pointer is null
	//std::vector<Constraint*> firstColumnsRowPointer;
	std::vector<QlpMatCoeff*> firstColumnsRowPointer;

	// Each Column has a Pointer to the last Row that contains a QlpMatCoeff with the columns index
	// and whoose value is not equal to zero, if such a point does not exist, the pointer is null
	//std::vector<Constraint*> lastColumnsRowPointer;
	std::vector<QlpMatCoeff*> lastColumnsRowPointer;


	// Pointer to the first Row Element of the Qlp
	Constraint* firstRow;
	// Pointer to the last Row Element of the Qlp
	Constraint* lastRow;
	std::map<std::string, unsigned int> nameToVariableIndex;
	std::map<unsigned int, std::string> variableIndexToName;

	inline std::string mapsToString() const {
		std::map<std::string, unsigned int>::const_iterator iter1;
		std::map<unsigned int, std::string>::const_iterator iter2;
		std::string strToReturn("\nnameToVariableIndex:\n");
		for (iter1 = nameToVariableIndex.begin(); iter1
				!= nameToVariableIndex.end(); iter1++) {
			strToReturn += "\t" + iter1->first;
			strToReturn += "=" + utils::ToolBox::convertToString(iter1->second);
			strToReturn += "\n";
		}
		strToReturn += "variableIndexToName:\n";
		for (iter2 = variableIndexToName.begin(); iter2
				!= variableIndexToName.end(); iter2++) {
			strToReturn += "\t" + utils::ToolBox::convertToString(iter2->first);
			strToReturn += "=" + iter2->second;
			strToReturn += "\n";
		}
		return strToReturn;
	}
};
}

#endif /* QUANTIFIEDPROBLEM_HPP_ */
