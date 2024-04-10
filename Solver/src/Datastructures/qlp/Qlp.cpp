/*
*
* Solver: Qlp.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Datastructures/qlp/Qlp.hpp"
#include "Datastructures/components/QpObjFunc.hpp"

namespace data {

std::string Qlp::LOG_TAG = "Qlp";

Qlp::Qlp() :
		objFunc(), variables(), firstColumnsRowPointer(), lastColumnsRowPointer(), firstRow(NULL), lastRow(NULL), nameToVariableIndex(), variableIndexToName() {
}

Qlp::Qlp(const data::QpObjFunc& o, const std::vector<data::QpVar>& v, const data::QpSparseMatrix& m, const std::vector<data::QpRhs>& r) :
		objFunc(), variables(), firstColumnsRowPointer(), lastColumnsRowPointer(), firstRow(), lastRow(), nameToVariableIndex(), variableIndexToName() {
	this->copyQlpContent(o, v, m, r);
}

Qlp::~Qlp() {
	this->clear();
}

Qlp& Qlp::operator =(const Qlp& rhs) {
	if (this == &rhs)
		return *this;
	this->clear();
	this->copyQlpContent(rhs);
	return *this;
}

Qlp& Qlp::initWithQlpParts(const data::QpObjFunc& o, const std::vector<data::QpVar>& v, const data::QpSparseMatrix& m, const std::vector<data::QpRhs>& r) {
	this->clear();
	this->copyQlpContent(o, v, m, r);
	return *this;
}

void Qlp::copyQlpContent(const data::Qlp& rhs) {
	this->objFunc = rhs.getObjectiveFunction();
	std::vector<data::QpVar const *> vV = rhs.getVariableVectorConst();
	this->variables = std::vector<data::QpVar*>(vV.size(), NULL);
	for (unsigned int i = 0; i < vV.size(); i++)
		this->variables[i] = new data::QpVar(*vV[i]);

	this->variableIndexToName = rhs.variableIndexToName;
	this->nameToVariableIndex = rhs.nameToVariableIndex;
	this->firstColumnsRowPointer = std::vector<QlpMatCoeff*>(vV.size(), NULL);
	this->lastColumnsRowPointer = std::vector<QlpMatCoeff*>(vV.size(), NULL);

	std::vector<Constraint const *> cV = rhs.getConstraintVecConst();
	for (unsigned int i = 0; i < cV.size(); i++) {
		this->createConstraint(*cV[i]);
	}
}

void Qlp::clear() {
	objFunc.clear();
	deleteAllRows();
	deleteAllColumns(); //TODO
	this->firstColumnsRowPointer.clear();
	this->lastColumnsRowPointer.clear();
	this->nameToVariableIndex.clear();
	this->variableIndexToName.clear();
}

void Qlp::copyQlpContent(const data::QpObjFunc& o, const std::vector<data::QpVar>& v, const data::QpSparseMatrix& m, const std::vector<data::QpRhs>& r) {

	if (m.size() != r.size())
		throw utils::DataStructureException("Qlp::Qlp(...) --> (m.size()!=rhsVec.size())");

	for (unsigned int i = 0, size = v.size(); i < size; i++) {
		this->createVariable(v[i]);
	}
	this->objFunc = o;
	for (unsigned int i = 0; i < r.size(); i++) {
		data::Constraint& c = this->createRhsConstraint(r[i]);
		c.createConstraintElements(m[i]);
	}
}

Qlp::QlpType Qlp::getQlpType() const {

	unsigned int exist, univ, rand;

	std::vector<const data::QpVar*> eV, uV, rV;
	exist = (eV = this->getVariableVectorByQuantifierConst(data::QpVar::exists)).size();
	univ = (uV = this->getVariableVectorByQuantifierConst(data::QpVar::all)).size();
	rand = (rV = this->getVariableVectorByQuantifierConst(data::QpVar::random)).size();

	if (exist && univ && rand) {
		utils::Logger::globalLog(utils::LOG_ERROR, LOG_TAG, "Mixed-Quantfied (randdom/all)");
		return QLP_TYPE_ERROR;
	}

	unsigned int eR = 0, eG = 0, eB = 0;
	for (unsigned int i = 0; i < eV.size(); i++) {
		data::QpVar::NumberSystem ns = eV[i]->getNumberSystem();
		ns == data::QpVar::real ? eR++ : (ns == data::QpVar::generals ? eG++ : eB++);
	}

	if (!(univ + rand)) {
		if (eR && !(eG || eB)) {
			return LP;
		}
		if (!eR && (eG || eB)) {
			return IP;
		}
		if (eR && (eG || eB)) {
			return MIP;
		}
	}

	unsigned int oR = 0, oG = 0, oB = 0;
	std::vector<const data::QpVar*>& oV = univ ? uV : rV;
	for (unsigned int i = 0; i < oV.size(); i++) {
		data::QpVar::NumberSystem ns = oV[i]->getNumberSystem();
		ns == data::QpVar::real ? oR++ : (ns == data::QpVar::generals ? oG++ : oB++);
	}

	if (eR && oR && !(eG || eB || oG || oB))
		return univ ? QLP : RQLP;
	if ((!eR && !oR)) {
		if (eB && oB && !(eG || oG)) {
			return univ ? QIP_BIN : RQIP_BIN;
		} else {
			return univ ? QIP : RQIP;
		}
	}

	if ((eR || oR) && (eG || eB || oG || oB)) {
		if (!eR)
			return QMIP_HDE;
		if (!oR)
			return QMIP_HDA;
		if ((eR && (eB || eG)) || (oR && (oB || oG)))
			return QMIP;
	}
	return QLP_TYPE_ERROR;

}

std::string Qlp::qlpTypeToString(const QlpType& t) {
	switch (t) {
	case LP:
		return "LP";
		break;
	case IP:
		return "IP";
		break;
	case MIP:
		return "MIP";
		break;
	case QLP:
		return "QLP";
		break;
	case QIP:
		return "QIP";
		break;
	case QIP_BIN:
		return "QIP (binary)";
		break;
	case QMIP:
		return "QMIP";
		break;
	case QMIP_HDE:
		return "QMIP (existentially half discrete)";
		break;
	case QMIP_HDA:
		return "QMIP (universally half discrete)";
		break;
	case RQLP:
		return "RQLP";
		break;
	case RQIP:
		return "RQIP";
		break;
	case RQIP_BIN:
		return "RQIP (binary)";
		break;
	case RQMIP:
		return "RQMIP";
		break;
	case RQMIP_HDE:
		return "RQMIP (existentially half discrete)";
		break;
	case RQMIP_HDA:
		return "RQMIP (universally half discrete)";
		break;
	case QLP_TYPE_ERROR:
		return "ERROR";
		break;
	default:
		return "Type not Not Supported";
	}
}

Qlp::Qlp(const Qlp& rhs) :
		objFunc(), variables(), firstColumnsRowPointer(), lastColumnsRowPointer(), firstRow(NULL), lastRow(NULL), nameToVariableIndex(), variableIndexToName() {
	this->copyQlpContent(rhs);
//	if(LOG_DEBUG)
//		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Copying Maps...");
//	utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Size: "+utils::ToolBox::convertToString(this->nameToVariableIndex.size()));
//	utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Size: "+utils::ToolBox::convertToString(this->variableIndexToName.size()));
//	this->nameToVariableIndex = rhs.nameToVariableIndex;
//	this->variableIndexToName = rhs.variableIndexToName;
//	if(LOG_DEBUG)
//		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Done.");
}

bool Qlp::operator==(const Qlp& rhs) const {
	if (this->objFunc.operator !=(rhs.objFunc)) {
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "this->objFunc != rhs.objFunc");
		return false;
	}

	for (unsigned int i = 0, size = this->getVariableCount(); i < size; i++) {
		if (this->getVarConst(i) != rhs.getVarConst(i)) {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "this->vars != rhs.vars");
			return false;
		}
	}

	const Constraint *lhsCon, *rhsCon;
	lhsCon = this->firstRow;
	rhsCon = rhs.firstRow;
	unsigned int i = 0;
	while (lhsCon && rhsCon) {
		if (lhsCon->operator ==(*rhsCon)) {
			lhsCon = lhsCon->getNextConstraintPtr();
			rhsCon = rhsCon->getNextConstraintPtr();
		} else {
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "index: " + utils::ToolBox::convertToString(i));
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "this->con != rhs.con");
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "this->con: " + lhsCon->toString());
			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "rhs->con: " + rhsCon->toString());
			return false;
		}
		i++;
	}
	return true;
}

bool Qlp::operator!=(const Qlp& rhs) const {
	return !this->operator ==(rhs);
}

////------------------------------------------------------------------------------>

unsigned int Qlp::getNumberSystemByQuantifier(data::QpVar::Quantifier q, data::QpVar::NumberSystem ns) const {
	unsigned int count = 0;
	std::vector<const data::QpVar *> vars = this->getVariableVectorByQuantifierConst(q);
	for (unsigned int i = 0; i < vars.size(); i++)
		if (vars[i]->getNumberSystem() == ns)
			count++;
	return count;
}

unsigned int Qlp::getQuantifierCount(data::QpVar::Quantifier q) const {
	unsigned int count = 0;
	for (unsigned int i = 0; i < variables.size(); i++) {
		if (variables[i]->getQuantifier() == q)
			count++;
	}
	return count;
}

unsigned int Qlp::getVariableCount() const {
	return this->variables.size();
}

unsigned int Qlp::getConstraintCount() const {
	unsigned int size = 0;
	if (!firstRow)
		return size;
	Constraint* tmp = firstRow;
	while (tmp) {
		size++;
		tmp = tmp->nextConstraint;
	}
	return size;
}

unsigned int Qlp::getMatrixElementCount() const {
	unsigned int size = 0;
	if (!firstRow)
		return size;
	Constraint* tmp = firstRow;
	while (tmp) {
		size += tmp->getElementCount();
		tmp = tmp->nextConstraint;
	}
	return size;
}

unsigned int Qlp::getStageCount() const {
	unsigned int stages = 1;
	std::vector<const data::QpVar *> vec = this->getVariableVectorByQuantifierConst(data::QpVar::exists);
	for (unsigned int i = 0; i < vec.size() - 1; i++)
		if (vec[i + 1]->getIndex() - vec[i]->getIndex() > 1)
			stages++;
	return stages;
}

unsigned int Qlp::getScenarioCount() const {
	return pow(2.0, (int) (this->variables.size() - this->getQuantifierCount(data::QpVar::exists)));
}

//-------------------- Methods to work with Objective Function ----------->
const data::QpNum& Qlp::getObjFuncOffset() const {
	return this->objFunc.getOffset();
}

void Qlp::setObjFuncOffset(const data::QpNum& val) {
	this->objFunc.setOffset(val);
}

void Qlp::reverseObjFunc() {
	this->objFunc.reverseObjectiveFuntion();
}

QpObjFunc::Objective Qlp::getObjective() const {
	return this->objFunc.getObjective();
}

void Qlp::setObjective(QpObjFunc::Objective obj) {
	this->objFunc.setObjective(obj);
}

const data::QpNum& Qlp::getObjectiveFunctionElement(unsigned int index) const {
	return this->objFunc[index];
}

void Qlp::setObjectiveFunctionElement(unsigned int index, const data::QpNum& coeff) {
	this->objFunc.setObjElement(index, coeff);
}

const std::vector<QpNum>& Qlp::getObjectiveFunctionValues() const {
	return this->objFunc.getObjectiveElementsDense();
}

data::QpNum Qlp::getObjectiveFunctionValue(const std::vector<data::QpNum>& allocation) const {
	data::QpNum value;
	if (allocation.size() == this->objFunc.getSize()) {
		for (unsigned int i = 0; i < allocation.size(); ++i) {
			value += (allocation[i] * objFunc[i]);
		}
	} else {
		throw utils::QlpSolverException("Allocation vector has not the same size as the objective funtion");
	}
	return value;
}

data::QpObjFunc& Qlp::getObjectiveFunction() {
	return this->objFunc;
}

const data::QpObjFunc& Qlp::getObjectiveFunction() const {
	return this->objFunc;
}

//-------------------- Methods to work with quantified variables ----------->
QpVar& Qlp::getVariableByIndex(unsigned int index) {
	return this->getVar(index);
}
const QpVar& Qlp::getVariableByIndexConst(unsigned int index) const {
	return this->getVarConst(index);
}

QpVar& Qlp::getVariableByName(const std::string& descr) {
	return this->getVar(nameToVariableIndex.find(descr)->second);
}
const QpVar& Qlp::getVariableByNameConst(const std::string& descr) const {
	return this->getVarConst(nameToVariableIndex.find(descr)->second);
}

std::vector<data::QpVar> Qlp::getQuantifiedVariables() const {
	std::vector<data::QpVar> varVec;
	for (unsigned int i = 0; i < variables.size(); i++) {
		varVec.push_back(*variables[i]);

	}
	return varVec;
}

std::vector<data::QpVar*> Qlp::getVariableVector() {
	return this->variables;
}

std::vector<data::QpVar const *> Qlp::getVariableVectorConst() const {
	std::vector<data::QpVar const *> varVec;
	for (unsigned int i = 0; i < variables.size(); i++) {
		varVec.push_back(variables[i]);

	}
	return varVec;
}

std::vector<data::QpVar*> Qlp::getVariableVectorByQuantifier(data::QpVar::Quantifier q) {
	std::vector<data::QpVar*> varVec;
	for (unsigned int i = 0; i < variables.size(); i++) {
		if (variables[i]->getQuantifier() == q)
			varVec.push_back(variables[i]);

	}
	return varVec;
}

std::vector<const data::QpVar *> Qlp::getVariableVectorByQuantifierConst(data::QpVar::Quantifier q) const {
	std::vector<data::QpVar const *> varVec;
	for (unsigned int i = 0; i < variables.size(); i++) {
		if (variables[i]->getQuantifier() == q)
			varVec.push_back(variables[i]);

	}
	return varVec;
}

QpVar& Qlp::createVariable(const data::QpVar& var) {
	QpVar& c = createVariable(var.getName(), var.getIndex(), var.getQuantifier());
	c.setNumberType(var.getNumberSystem());
	if (var.getQuantifier() != data::QpVar::random) {
		c.setLowerBound(var.getLowerBound());
		c.setUpperBound(var.getUpperBound());
	} else {
		c.setVariableRange(var.getVariableRange(), var.getVariableDistribution());
	}
	return c;
}

QpVar& Qlp::createVariable(const std::string& descr, int index, QpVar::Quantifier q) {

	if (nameToVariableIndex.count(descr)) {
		utils::Logger::globalLog(utils::LOG_INSANE, LOG_TAG, "Qlp: " + this->toString());
		throw utils::DataStructureException("Qlp::createVariable with name already exists: " + descr);
	}

	nameToVariableIndex.insert(std::pair<std::string, int>(descr, index));
	variableIndexToName.insert(std::pair<int, std::string>(index, descr));
	data::QpVar* v = new QpVar(descr, index, data::QpNum(true), data::QpNum(false));
	v->setQuantifier(q);
	// add the column to the vector of columns
	variables.push_back(v);
	// create space for the first and last row element for this column
	firstColumnsRowPointer.push_back(NULL);
	lastColumnsRowPointer.push_back(NULL);
	this->objFunc.setSize(this->objFunc.getSize() + 1); //TODO
	if (LOG_QLP)
		utils::Logger::globalLog(utils::LOG_INSANE, LOG_TAG, "Column " + descr + "added with index " + utils::ToolBox::convertToString(index));
	return *v;
}

QpVar& Qlp::createVariable(const std::string& descr, int index, QpVar::Quantifier q, QpVar::NumberSystem n, const data::QpNum& lb, const data::QpNum& ub) {
	QpVar& c = createVariable(descr, index, q);
	c.setNumberType(n);
	c.setLowerBound(lb);
	c.setUpperBound(ub);
	return c;
}

QpVar& Qlp::createVariable(const std::string& descr, int index, QpVar::Quantifier q, QpVar::NumberSystem n, const std::vector<data::QpNum>& vRange, const std::vector<data::QpRational>& vDistr) {
	QpVar& c = createVariable(descr, index, q);
	c.setNumberType(n);
	c.setVariableRange(vRange, vDistr);
	return c;
}

std::vector<Constraint*> Qlp::getVariableConstraints(unsigned int index) {
	std::vector<Constraint*> cVector;
	if (isVariableUnused(index))
		return cVector;
	QlpMatCoeff* firstColPoint = getFirstColumnsRowPtr(index)->getPointPtrByIndex(index);
	cVector.push_back(firstColPoint->pRow);
	while (!firstColPoint->isLastColElement()) {
		firstColPoint = firstColPoint->pNextColPoint;
		cVector.push_back(firstColPoint->pRow);
	}
	return cVector;
}

std::vector<Constraint const *> Qlp::getVariableConstraintsConst(unsigned int index) const {
	std::vector<const Constraint*> cVector;
	if (isVariableUnused(index))
		return cVector;
	QlpMatCoeff* firstColPoint = getFirstColumnsRowPtr(index)->getPointPtrByIndex(index);
	cVector.push_back(firstColPoint->pRow);
	while (!firstColPoint->isLastColElement()) {
		firstColPoint = firstColPoint->pNextColPoint;
		cVector.push_back(firstColPoint->pRow);
	}
	return cVector;
}

Constraint& Qlp::createConstraint(const data::Constraint& con) {
	Constraint& c = this->createRhsConstraint(con.getRhsRatioSign(), con.getRhsValue(), con.getResponsibility());
	c.createConstraintElements(con);
	return c;
}

Constraint& Qlp::createRhsConstraint(const data::QpRhs& rhs) {
	return this->createRhsConstraint(rhs.getRatioSign(), rhs.getValue(), rhs.getResponsibility());
}

Constraint& Qlp::createRhsConstraint(QpRhs::RatioSign r, const data::QpNum& f, int resp) {
	data::QpRhs rhs(f, r, (data::QpRhs::Responsibility)resp);
	Constraint* p = new Constraint(*this, rhs, lastRow);
	if (firstRow == NULL) {
		firstRow = p;
	}
	lastRow = p;
	if (LOG_QLP)
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "Row added. New RowCount: " + utils::ToolBox::convertToString(getConstraintCount()));
	//this->constraintCount++;
	return *p;
}

Constraint* Qlp::getFirstConstraint() const {
	return firstRow;
}

Constraint* Qlp::getLastConstraint() const {
	return lastRow;
}

std::list<Constraint*> Qlp::getConstraints() const {
	std::list<Constraint*> cList;
	if (this->getConstraintCount()) {
		Constraint* tmpRow = firstRow;
		if (tmpRow != NULL && !tmpRow->isLastConstraint()) {
			do {
				cList.push_back(tmpRow);
				tmpRow = (Constraint*) tmpRow->nextConstraint;
			} while (!tmpRow->isLastConstraint());
			cList.push_back(tmpRow);
		} else if (tmpRow->isLastConstraint()) {
			cList.push_back(tmpRow);
		}
	}
	return cList;
}

std::vector<Constraint*> Qlp::getConstraintVec() {
	std::vector<Constraint*> cVec;
	data::Constraint* c = this->getFirstConstraint();
	while (c) {
		cVec.push_back(c);
		c = c->getNextConstraintPtr();
	}
	return cVec;
}

std::vector<Constraint const *> Qlp::getConstraintVecConst() const {
	std::vector<Constraint const *> cVec;
	data::Constraint* c = this->getFirstConstraint();
	while (c) {
		cVec.push_back(c);
		c = c->getNextConstraintPtr();
	}
	return cVec;

}

std::vector<QpNum> Qlp::getRhsValVec() const {
	std::vector<QpNum> vec;
	data::Constraint* c = this->getFirstConstraint();
	while (c != NULL) {
		vec.push_back(c->getRhsValue());
		c = c->getNextConstraintPtr();
	}
	return vec;
}

std::vector<data::QpRhs> Qlp::getRhsVec() const {
	std::vector<data::QpRhs> vec;
	data::Constraint* c = this->getFirstConstraint();
	while (c != NULL) {
		data::QpRhs r(c->getRhsValue(), c->getRhsRatioSign());
		vec.push_back(r);
		c = c->getNextConstraintPtr();
	}
	return vec;
}

std::vector<const data::QpNum *> Qlp::getRhsValVecConst() const {
	std::vector<const data::QpNum *> vec;
	data::Constraint* c = this->getFirstConstraint();
	while (c != NULL) {
		vec.push_back(&(c->getRhs().getValue()));
		c = c->getNextConstraintPtr();
	}
	return vec;
}

std::vector<const data::QpRhs *> Qlp::getRhsVecConst() const {
	std::vector<const data::QpRhs *> vec;
	data::Constraint* c = this->getFirstConstraint();
	while (c != NULL) {
		vec.push_back(&(c->getRhs()));
		c = c->getNextConstraintPtr();
	}
	return vec;
}

void Qlp::getCoeffMatrix(std::vector<std::vector<data::IndexedElement> >&m) const {
	m.clear();
	data::Constraint* c = this->getFirstConstraint();
	while (c != NULL) {
		m.push_back(c->getElements());
		c = c->getNextConstraintPtr();
	}
}

void Qlp::getCoeffMatrixByResp(std::vector<std::vector<data::IndexedElement> >&m, data::QpRhs::Responsibility resp) const {
	m.clear();
	data::Constraint* c = this->getFirstConstraint();
	while (c != NULL) {
		if (c->getResponsibility()==resp){
			m.push_back(c->getElements());
		}
		c = c->getNextConstraintPtr();

	}
}

std::vector<data::QpRhs> Qlp::getRhsVecByResp(data::QpRhs::Responsibility resp) const{
	std::vector<data::QpRhs> vec;
	data::Constraint* c = this->getFirstConstraint();
	while (c != NULL) {
		if (c->getResponsibility()==resp){
			data::QpRhs r(c->getRhsValue(), c->getRhsRatioSign());
			vec.push_back(r);
		}
		c = c->getNextConstraintPtr();

	}
	return vec;
}


//data::QpSparseMatrix Qlp::getSparseCoeffMatrix() const {
//
//	data::QpSparseMatrix
//			m(this->getConstraintCount(), this->getVariableCount());
//
//	data::Constraint* c = this->getFirstConstraint();
//	unsigned int i = 0;
//	while (c != NULL) {
//		m.insert(i, c->getElements());
//		c = c->getNextConstraintPtr();
//		i++;
//	}
//
//	return m;
//}
bool Qlp::isQuantifiedProblem() const {
	return (this->getQuantifierCount(data::QpVar::all) || this->getQuantifierCount(data::QpVar::random));
}
bool Qlp::containsNumberSystem(data::QpVar::NumberSystem n) const {
	for (unsigned int i = 0; i < variables.size(); i++) {
		if (variables[i]->getNumberSystem() == n)
			return true;
	}
	return false;
//	std::list<data::QpVar*>::const_iterator it = variables.begin();
//	std::list<data::QpVar*>::const_iterator end = variables.end();
//	while (it != end) {
//		if ((*it)->getNumberSystem() == num)
//			return true;
//		it++;
//	}
//	return false;
}

bool Qlp::isVariableUnused(unsigned int index) const {
	return !(getFirstColumnsRowPtr(index) && getLastColumnsRowPtr(index));
}

unsigned int Qlp::getVariableConstraintCount(unsigned int index) const {
	if (isVariableUnused(index))
		return 0;
	QlpMatCoeff* firstColPoint = getFirstColumnsRowPtr(index)->getPointPtrByIndex(index);
	int size = 1;
	while (!firstColPoint->isLastColElement()) {
		if ((firstColPoint = firstColPoint->pNextColPoint) == NULL)
			throw utils::QlpSolverException("NullPointer. Calling method -> Col::getColSIze()");
		size++;
	}
	return size;
}

void Qlp::eliminateVariables(const std::vector<data::IndexedElement>& allocations) {
	for (unsigned int i = 0; i < allocations.size(); i++)
		eliminateVariable(allocations[i].index, allocations[i].value);
}

void Qlp::eliminateVariable(unsigned int index, const data::QpNum& value) {
	data::QpNum coeff;

	if ((coeff = this->getObjectiveFunctionElement(index)) != 0) {
		this->getObjectiveFunction().setOffset(this->getObjectiveFunction().getOffset() + coeff * value);
		this->setObjectiveFunctionElement(index, 0.0);
	}

	setVariableValue(index, value);
	//this->getVar(index).setEliminated(); //TODO
	firstColumnsRowPointer[index] = lastColumnsRowPointer[index] = 0;
}

void Qlp::setVariableValue(unsigned int index, const data::QpNum& value) {
	QpVar& c = this->getVar(index);
	std::vector<Constraint*> constraints = this->getVariableConstraints(c.getIndex());
	for (unsigned int i = 0; i < constraints.size(); i++) {
		constraints[i]->eliminateConstraintElement(index, value);
	}
}

void Qlp::normalizeConstraints() {
	if (getConstraintCount()) {
		Constraint* tmp = firstRow;
		while (tmp) {
			tmp->normalizeConstraint();
			tmp = (Constraint*) tmp->getNextConstraintPtr();
		}
	}
}

void Qlp::getQlpParts(data::QpObjFunc& obj, std::vector<data::QpVar>& vars, data::QpSparseMatrix& matrix, std::vector<data::QpRhs>& rhs) const {

	obj.setObjective(this->objFunc.getObjective());
	obj.setOffset(this->objFunc.getOffset());
	obj.setObjElements(this->objFunc.getObjectiveElementsDense());

	for (unsigned int i = 0; i < this->variables.size(); i++)
		vars.push_back(*this->variables[i]);

	std::vector<const data::QpRhs *> rhsSource = this->getRhsVecConst();
	for (unsigned int i = 0; i < rhsSource.size(); i++)
		rhs.push_back(*rhsSource[i]);

	this->getCoeffMatrix(matrix);

}

void Qlp::sortQlp() {
	Constraint* tmp = this->firstRow;
	if (!tmp)
		return;
	bool sort = false;
	tmp = this->firstRow->getNextConstraintPtr();
	while (tmp) {
		if (tmp->getLastConstraintElementIndex() < tmp->getPreviousConstraintPtr()->getLastConstraintElementIndex()) {
			sort = true;
			break;
		}
		tmp = tmp->nextConstraint;
	}

	if (!sort) {
		return;
	}

	std::vector<data::QpRhs> rhsVec;
	std::vector<std::vector<data::IndexedElement> > matrix;
	unsigned int index = 0;

	std::multimap<int, unsigned int> conSizeMap;
	std::multimap<int, unsigned int>::iterator it = conSizeMap.begin();
	tmp = this->firstRow;
	while (tmp) {
		it = conSizeMap.insert(it, std::pair<int, unsigned int>(tmp->getLastConstraintElementIndex(), index));
		rhsVec.push_back(data::QpRhs(tmp->getRhs()));
		matrix.push_back(std::vector<data::IndexedElement>());
		tmp->getElementsSparse(matrix[index]);
		tmp = tmp->nextConstraint;
		index++;
	}
	this->deleteAllRows();
	for (it = conSizeMap.begin(); it != conSizeMap.end(); it++) {
		data::Constraint& c = this->createRhsConstraint(rhsVec[(*it).second]);
		std::vector<data::IndexedElement> row = matrix[(*it).second];
		for (unsigned int i = 0; i < row.size(); i++) {
			c.createConstraintElement(row[i].index, row[i].value);
		}
	}
}

const std::string& Qlp::getVariableNameByIndex(unsigned int index) const {
	std::map<unsigned int, std::string>::const_iterator iter = variableIndexToName.begin();
	if ((iter = variableIndexToName.find(index)) != variableIndexToName.end()) {
		return iter->second;
	} else {
		throw utils::DataStructureException("std::string& Qlp::getVariableNameByIndex(unsigned int index) --> not in map: " + utils::ToolBox::convertToString(index));
	}
}

int Qlp::getVariableIndexByName(const std::string& varName) const {
	std::map<std::string, unsigned int>::const_iterator iter = nameToVariableIndex.begin();
	if ((iter = nameToVariableIndex.find(varName)) != nameToVariableIndex.end()) {
		return iter->second;
	} else {
		return -1;
	}
}

std::string Qlp::toQlpFileString(bool lp) const {

	if (LOG_QLP)
		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Qlp::toQlpFileString(bool lp)");

	std::string begin(data::QpObjFunc::objectiveString(this->getObjective()));
	begin += "\n";

	std::vector<data::IndexedElement> ieVec = this->objFunc.getObjectiveElementsSparse();
	for (unsigned int i = 0; i < ieVec.size(); i++) {
		if (i != 0)
			begin += "+ ";
		begin += ieVec[i].value.toString();
		begin += getVariableNameByIndex(ieVec[i].index);
		begin += " ";
	}

	begin += "\nSUBJECT TO";
	if (getConstraintCount()) {
		Constraint* tmp = firstRow;
		while (tmp) {
			begin += "\n " + tmp->getDoubleString();
			tmp = tmp->nextConstraint;
		}
	}

	begin += "\nBOUNDS";
	std::vector<data::QpVar*>::const_iterator it = this->variables.begin();
	std::vector<data::QpVar*>::const_iterator end = this->variables.end();

	unsigned int i = 0;
	bool integral = false, binaries = false;
	while (it != end) {
		if ((*it)->getNumberSystem() != data::QpVar::real) {
			if ((*it)->getNumberSystem() == data::QpVar::generals)
				integral = true;
			if ((*it)->getNumberSystem() == data::QpVar::binaries)
				binaries = true;
		}
		begin += "\n";
		begin += (*it)->getLowerBound().toString();
		begin += "<=";
		begin += this->getVariableNameByIndex(i);
		begin += "<=";
		begin += (*it)->getUpperBound().toString();
		it++;
		i++;
	}
	if (integral) {
		begin += "\nGENERAL\n";
		i = 0;
		it = this->variables.begin();
		while (it != end) {
			if ((*it)->getNumberSystem() == data::QpVar::generals) {
				begin += getVariableNameByIndex(i) + "\n ";
			}
			it++;
			i++;
		}
	}

	if (binaries) {
		begin += "\nBINARY\n";
		i = 0;
		it = this->variables.begin();
		while (it != end) {
			if ((*it)->getNumberSystem() == data::QpVar::binaries) {
				begin += getVariableNameByIndex(i) + "\n ";
			}
			it++;
			i++;
		}
	}

	if (lp) {
		begin += "\nEND";
		return begin;
	}

	begin += "\nEXISTS\n";
	i = 0;
	it = this->variables.begin();
	while (it != end) {
		if ((*it)->getQuantifier() == data::QpVar::exists) {
			begin += getVariableNameByIndex(i) + " ";
		}
		it++;
		i++;
	}

	if (this->getQuantifierCount(data::QpVar::all)) {
		begin += "\nALL\n";
		i = 0;
		it = this->variables.begin();
		while (it != end) {
			if ((*it)->getQuantifier() == data::QpVar::all) {
				begin += getVariableNameByIndex(i) + " ";
			}
			it++;
			i++;
		}
	}

	if (this->getQuantifierCount(data::QpVar::random)) {
		begin += "\nRANDOM\n";
		i = 0;
		it = this->variables.begin();
		while (it != end) {
			if ((*it)->getQuantifier() == data::QpVar::random) {
				begin += getVariableNameByIndex(i) + " ";
			}
			it++;
			i++;
		}
	}

	begin += "\nORDER\n";
	i = 0;
	it = this->variables.begin();
	while (it != end) {
		begin += getVariableNameByIndex(i) + " ";
		it++;
		i++;
	}

	begin += "\nEND";
	return begin;
}

std::string Qlp::getQlpInfo() const {
	std::string s("Qlp [ ");
	s += "REAL: " + utils::ToolBox::convertToString(this->getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::real));
	s += ", INT: " + utils::ToolBox::convertToString(this->getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::generals));
	s += ", BIN: " + utils::ToolBox::convertToString(this->getNumberSystemByQuantifier(data::QpVar::exists, data::QpVar::binaries));
	s += ", Constraints: " + utils::ToolBox::convertToString(this->getConstraintCount());
	s += "]";
	return s;
}

std::string Qlp::toString() const {
	if (LOG_QLP)
		utils::Logger::globalLog(utils::LOG_DEBUG, LOG_TAG, "Trying to print Qlp");

	const data::QpVar* v;

	std::string begin("");
	begin += "\n\tObjective Function: \n\t\t";
	if (variables.size() > 1000000) {
		begin += data::QpObjFunc::objectiveString(this->objFunc.getObjective());
		begin += " to many variables to display(";
		begin += utils::ToolBox::convertToString(variables.size());
		begin += ")";
	} else {
		this->objFunc.getObjective() == data::QpObjFunc::max ? begin += "max " : begin += "min ";
		if (!objFunc.getOffset().isZero())
			begin += objFunc.getOffset().toString();
		std::vector<data::IndexedElement> iList = this->objFunc.getObjectiveElementsSparse();
		bool first = true;
		for (unsigned int i = 0; i < iList.size(); i++) {
			if (first && !objFunc.getOffset().isZero())
				begin += " +";
			if (!first && iList[i].value > 0.0)
				begin += " +";
			begin += iList[i].value.toString();
			begin += this->getVariableNameByIndex(iList[i].index);
			first = false;
		}

	}
	begin += "\n\tQuantifiers: \n";
	begin += "\t\tExists: " + utils::ToolBox::convertToString(this->getQuantifierCount(data::QpVar::exists)) + "\n";
	begin += "\t\tAll: " + utils::ToolBox::convertToString(this->getQuantifierCount(data::QpVar::all)) + "\n";
	begin += "\t\tRandom: " + utils::ToolBox::convertToString(this->getQuantifierCount(data::QpVar::random)) + "\n";

	begin += "\n\tList sizes: \n";
	begin += "\t\tColumns: " + utils::ToolBox::convertToString(variables.size()) + "\n";
	begin += "\t\tFirstColumnsRowPointer: " + utils::ToolBox::convertToString(firstColumnsRowPointer.size()) + "\n";
	begin += "\t\tLastColumnsRowPointer: " + utils::ToolBox::convertToString(lastColumnsRowPointer.size()) + "\n";
	begin += "\t\tRows: " + utils::ToolBox::convertToString(getConstraintCount()) + "\n";

	if (variables.size() == 0) {
		begin += "\n\tVariables: No columns assigned";
	} else {
		begin += "\n\tVariables: ";
	}
	if (variables.size() > 1000000) {
		begin += "to many variables to display(";
		begin += utils::ToolBox::convertToString(variables.size());
		begin += ")";
	} else {
		for (unsigned int i = 0; i < variables.size(); i++) {
			v = &this->getVarConst(i);
			begin += "\n\t\t ";
			begin += this->getVariableNameByIndex(i);
			begin += " (" + utils::ToolBox::convertToString(i);
			begin += ")[";
			begin += data::QpNum::vecToString(v->getVariableRange());
			begin += data::QpNum::vecToString(v->getVariableDistribution());

			begin += "] ";	//
			if (v->getQuantifier() == data::QpVar::exists) {
				begin += "[exist,";
			} else if (v->getQuantifier() == data::QpVar::all) {
				begin += "[all,";
			} else {
				begin += "[random,";
			}
			if (v->getNumberSystem() == data::QpVar::real) {
				begin += "real]";
			} else if (v->getNumberSystem() == data::QpVar::generals) {
				begin += "general]";
			} else {
				begin += "binary]";
			}
		}
	}

	if (getConstraintCount()) {
		begin += "\n\tRows:";
		if (getConstraintCount() > 5000000) {
			begin += "to many rows to display";
			begin += utils::ToolBox::convertToString(getConstraintCount());
			begin += ")";
			return begin;
		}
		begin += "\n";
		int i = 1;
		Constraint* tmp = firstRow;
		while (tmp != lastRow) {
			begin += "\t\tc" + utils::ToolBox::convertToString(i);
			begin += " [" + utils::ToolBox::convertToString(tmp->elements);
			begin += "]";
			begin += ":\t";
			begin += " " + tmp->toString();
			begin += "\n";
			tmp = (Constraint*) tmp->nextConstraint;
			i++;
		}
		begin += "\t\tc" + utils::ToolBox::convertToString(i);
		begin += " [" + utils::ToolBox::convertToString(lastRow->elements);
		begin += "]";
		begin += ":\t";
		begin += " " + lastRow->toString();
		begin += "\n";
	} else {
		begin += "\n\tRows: No Rows assigned \n";
	}
	return begin;
}

//--------------------------------- PROTECTED METHOD IMPLEMENTATIONS ----------------------------------------->
void Qlp::deleteAllRows() {
	if (firstRow) {
		Constraint* tmp = firstRow;
		while (!tmp->isLastConstraint()) {
			tmp = (Constraint*) tmp->nextConstraint;
			tmp->previousConstraint->removeConstraint();
		}
		tmp->removeConstraint();
	}

//
//	if(firstRow){
//		data::Constraint* tmp = firstRow,*prev;
//		while(tmp){
//			prev=tmp;
//
//			QlpMatCoeff* tmpP = prev->pFirstRowPoint,*tmpPp;
//			while(tmpP){
//				tmpPp=tmpP;
//				tmpP=tmpP->pNextRowPoint;
//				delete tmpPp;
//			}
//
//			tmp=tmp+1;
//			delete prev;
//		}
//	}
	for (unsigned int i = 0; i < this->firstColumnsRowPointer.size(); i++) {
		this->firstColumnsRowPointer[i] = NULL;
		this->lastColumnsRowPointer[i] = NULL;
	}
	this->firstRow = this->lastRow = NULL;
}

void Qlp::deleteAllColumns() {
	for (unsigned int i = 0; i < variables.size(); i++) {
		delete variables[i];
	}
	variables.clear();
	this->nameToVariableIndex.clear();
	this->variableIndexToName.clear();
	this->firstColumnsRowPointer.clear();
	this->lastColumnsRowPointer.clear();

}

void Qlp::deleteColumnByIndex(unsigned int index) {

	if (index >= variables.size())
		throw utils::DataStructureException("Qlp::deleteColumnByIndex(int index)");

	data::QpVar& v = getVar(index);
	Constraint* tmpFirst = firstColumnsRowPointer[index]->pRow;

	QlpMatCoeff *tmpPoint, *tmpPointToDelete;

	if (tmpFirst && (tmpPoint = tmpFirst->getPointPtrByIndex(index))) {
		while (!tmpPoint->isLastColElement()) {
			if (!(tmpPoint = tmpPoint->pNextColPoint))
				throw utils::QlpSolverException("Calling method -> deleteColumnByIndex()");
			if (!(tmpPointToDelete = tmpPoint->pNextColPoint))
				throw utils::QlpSolverException("Calling method -> deleteColumnByIndex()[Point to delete]");
			tmpPointToDelete->deleteFromRowAndCol();
		}
		tmpPoint->deleteFromRowAndCol();
	}
	delete &v;
}

//void Qlp::printColumnRowPointers() {
//
//	utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "PrintColumnRowPointers");
//
//	if (firstColumnsRowPointer.size() != this->variables.size())
//		throw utils::DataStructureException(
//				"firstColumnsRowPointer.size()!=this->variables.size()");
//	if (firstColumnsRowPointer.size() != this->variables.size())
//		throw utils::DataStructureException(
//				"firstColumnsRowPointer.size()!=this->variables.size()");
//
//	for (unsigned int i = 0; i < variables.size(); i++) {
//
//		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG, "Variable: "
//				+ getVar(i).toString());
//		if (firstColumnsRowPointer[i] != NULL) {
//			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
//					"firstColumnsRowPointer[i] = "
//							+ firstColumnsRowPointer[i]->toString());
//		} else {
//			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
//					"firstColumnsRowPointer[i]=NULL");
//		}
//
//		if (lastColumnsRowPointer[i] != NULL) {
//			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
//					"lastColumnsRowPointer[i] = "
//							+ lastColumnsRowPointer[i]->toString());
//		} else {
//			utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
//					"lastColumnsRowPointer[i]=NULL");
//		}
//
//	}
//
//}

void Qlp::getSortedVariableConstraints(unsigned int index, std::vector<data::Constraint*>& constraints, std::vector<data::Constraint*>& cLessThanOrEqual, std::vector<data::Constraint*>& cGreaterThanOrEqual, std::vector<data::Constraint*>& cEquals, std::vector<data::Constraint*>& cRest) {

	cLessThanOrEqual.clear();
	cGreaterThanOrEqual.clear();
	cEquals.clear();
	cRest.clear();

	bool exists = false;
	if (this->getVar(index).getQuantifier() == data::QpVar::exists)
		exists = true;

	Constraint* constraint;
	QpRhs::RatioSign ratioSign;

	constraint = this->firstRow;
	while (constraint) {
		constraint->normalizeConstraint();
		if (constraint->pLastRowPoint->columnIndex != index) {
			cRest.push_back(constraint);
		} else {

			constraints.push_back(constraint);

			ratioSign = constraint->getRhsRatioSign();
			if (ratioSign == data::QpRhs::smallerThanOrEqual) {
				cLessThanOrEqual.push_back(constraint);
			} else if (ratioSign == data::QpRhs::greaterThanOrEqual) {
				cGreaterThanOrEqual.push_back(constraint);
			} else {
				cEquals.push_back(constraint);
			}
		}
		constraint = constraint->nextConstraint;
	}
}

std::string Qlp::getVariableStringByIndex(unsigned int index) const {
	std::string begin("");
	const data::QpVar* v;
	v = &this->getVarConst(index);
	begin += "\t ";
	begin += this->getVariableNameByIndex(index);
	begin += " [";
	if (!v->isBoundedBelow())
		begin += "-inf";
	else
		begin += v->getLowerBound().toString();
	begin += ",";
	if (!v->isBoundedAbove())
		begin += "+inf";
	else
		begin += v->getUpperBound().toString();
	begin += "] ";
	if (v->getQuantifier() == data::QpVar::exists) {
		begin += "E";
	} else if (v->getQuantifier() == data::QpVar::all) {
		begin += "A";
	} else {
		begin += "R";
	}
	return begin;
}

//void Qlp::printVariables() const {
//	for (unsigned int i = 0, size = this->getVariableCount(); i < size; i++)
//		utils::Logger::globalLog(utils::LOG_INFO, LOG_TAG,
//				this->getVariableStringByIndex(i));
//
//}

}
