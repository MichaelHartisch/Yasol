/*
*
* Solver: Constraint.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Datastructures/qlp/Constraint.hpp"
#include "Datastructures/qlp/Qlp.hpp"
namespace data {

std::string Constraint::LOG_TAG = "Constraint";

Constraint::Constraint(Qlp& p) :
		qlp(p), rhs(0, QpRhs::smallerThanOrEqual), elements(0), nextConstraint(NULL), previousConstraint(NULL), pFirstRowPoint(NULL), pLastRowPoint(NULL) {
}

Constraint::Constraint(Qlp& p, const data::QpRhs& r, Constraint* pCons) :
		qlp(p), rhs(r), elements(0), nextConstraint(NULL), previousConstraint(pCons), pFirstRowPoint(NULL), pLastRowPoint(NULL) {
	if (pCons)
		pCons->nextConstraint = this;
}

Constraint::~Constraint() {
	if (previousConstraint && nextConstraint) {
		nextConstraint->previousConstraint = this->previousConstraint;
		previousConstraint->nextConstraint = this->nextConstraint;
	} else if (previousConstraint && !nextConstraint) {
		previousConstraint->nextConstraint = NULL;
	} else if (!previousConstraint && nextConstraint) {
		nextConstraint->previousConstraint = NULL;
	}
	this->nextConstraint = this->previousConstraint = NULL;
}

bool Constraint::isEmpty() const {
	return !(pLastRowPoint && pFirstRowPoint);
}

unsigned int Constraint::getElementCount() const {
	return elements;
}

void Constraint::setRhsValue(const QpNum& v) {
	rhs.setValue(v);
}

const QpNum& Constraint::getRhsValue() const {
	return rhs.getValue();
}

void Constraint::setResponsibility(int r) {
    rhs.setResponsibility(r);
}
    
const int Constraint::getResponsibility() const {
    return rhs.getResponsibility();
}
    
void Constraint::setRhsRatioSign(QpRhs::RatioSign r) {
	rhs.setRatioSign(r);
}

QpRhs::RatioSign Constraint::getRhsRatioSign() const {
	return rhs.getRatioSign();
}

void Constraint::reverseRhsRatioSign() {
	rhs.reverseRatioSign();
}

const data::QpRhs& Constraint::getRhs() const {
	return rhs;
}

bool Constraint::isFirstConstraint() const {
	return ((previousConstraint == NULL) ? true : false);
}

bool Constraint::isLastConstraint() const {
	return ((nextConstraint == NULL) ? true : false);
}

bool Constraint::operator==(const data::Constraint& c) const {
	if ((this->elements != c.elements) || (this->rhs != c.rhs))
		return false;
	QlpMatCoeff* thisTmpPoint = pFirstRowPoint;
	QlpMatCoeff* rTmpPoint = c.pFirstRowPoint;
	while (thisTmpPoint && rTmpPoint) {
		if (!thisTmpPoint->equals(*rTmpPoint))
			return false;
		thisTmpPoint = thisTmpPoint->pNextRowPoint;
		rTmpPoint = rTmpPoint->pNextRowPoint;
	}
	return true;
}

//-------------------------------- Operations on Constraints ----------------->
void Constraint::removeConstraint() {
	if (!isEmpty()) {
		if (pFirstRowPoint == pLastRowPoint) {
			pFirstRowPoint->deleteFromCol();
		} else {
			do {
				pFirstRowPoint->deleteFromCol();
				pFirstRowPoint = pFirstRowPoint->pNextRowPoint;
				delete pFirstRowPoint->pPreviousRowPoint;
			} while (!pFirstRowPoint->isLastRowElement());
		} // TODO this could be done much better
		pLastRowPoint->deleteFromCol();
		delete pLastRowPoint;
	}

	if (nextConstraint) {
		nextConstraint->previousConstraint = this->previousConstraint;
	} else {
		qlp.setLastRowPtr(previousConstraint);
	}
	if (previousConstraint) {
		previousConstraint->nextConstraint = this->nextConstraint;
	} else {
		qlp.setFirstRowPtr(nextConstraint);
	}
	delete this;
}

void Constraint::multiplyConstraint(const data::QpNum& multiplier) {
	rhs *= multiplier;
	QlpMatCoeff* tmp = pFirstRowPoint;
	while (tmp) {
		tmp->value *= multiplier;
		tmp = tmp->pNextRowPoint;
	}
}

void Constraint::devideConstraint(const data::QpNum& divisor) {
	rhs /= divisor;
	QlpMatCoeff* tmp = pFirstRowPoint;
	while (tmp) {
		tmp->value /= divisor;
		tmp = tmp->pNextRowPoint;
	}
}

Constraint& Constraint::operator +(const data::Constraint& c) const {
	Constraint& r = (Constraint&) qlp.createRhsConstraint(rhs.getRatioSign(), rhs.getValue() - c.getRhsValue(), rhs.getResponsibility());
	QlpMatCoeff* tmp = pFirstRowPoint;
	while (tmp) {
		r.createConstraintElement(tmp->columnIndex, tmp->value);
		tmp = tmp->pNextRowPoint;
	}
	tmp = c.pFirstRowPoint;
	while (tmp) {
		r.createConstraintElement(tmp->columnIndex, tmp->value);
		tmp = tmp->pNextRowPoint;
	}
	return r;
}

Constraint& Constraint::operator -(const data::Constraint& c) const {
	Constraint& r = (Constraint&) qlp.createRhsConstraint(rhs.getRatioSign(), rhs.getValue() - c.getRhsValue(), rhs.getResponsibility());
	QlpMatCoeff* tmp = pFirstRowPoint;
	while (tmp) {
		r.createConstraintElement(tmp->columnIndex, tmp->value);
		tmp = tmp->pNextRowPoint;
	}
	tmp = c.pFirstRowPoint;
	while (tmp) {
		r.createConstraintElement(tmp->columnIndex, (-1) * tmp->value.asDouble());
		tmp = tmp->pNextRowPoint;
	}
	return r;
}

Constraint& Constraint::mergeConstraints(const data::Constraint& c) const {

	if (rhs.getRatioSign() == data::QpRhs::equal || c.getRhsRatioSign() == data::QpRhs::equal)
		throw utils::QlpSolverException("Can't merge constraints with equality sign in this method");

	data::QpNum f(rhs.getValue());
	f -= c.getRhsValue();
	//utils::Logger::globalLog(utils::LOG_INFO,LOG_TAG,"New Rhs Input: "+f.toString());
	Constraint& r = qlp.createRhsConstraint(rhs.getRatioSign(), f, rhs.getResponsibility());

	//utils::Logger::globalLog(utils::LOG_INFO,LOG_TAG,"New Constraint before adding: "+r.toString());

	QlpMatCoeff* currTmp = pFirstRowPoint;
	QlpMatCoeff* mergeTmp = c.pFirstRowPoint;

	if (currTmp && mergeTmp)
		while (1) {
			if (currTmp->columnIndex == mergeTmp->columnIndex) {
				//the points have the same index, thus we merge their coefficient
				//and create a new point if the new coeff is bigger than 0
				const data::QpNum& v1 = currTmp->value;
				const data::QpNum& v2 = mergeTmp->value;
				if ((v1 - v2) != 0.0)
					r.createConstraintElement(currTmp->columnIndex, v1 - v2);
				if (((currTmp = currTmp->pNextRowPoint) == NULL) | ((mergeTmp = mergeTmp->pNextRowPoint) == NULL))
					break;
			} else if (currTmp->columnIndex < mergeTmp->columnIndex) {
				//this points has a lower index than the current lowest index of the corresponding merge row
				r.createConstraintElement(currTmp->columnIndex, currTmp->value);
				if ((currTmp = currTmp->pNextRowPoint) == NULL)
					break;
			} else if (currTmp->columnIndex > mergeTmp->columnIndex) {
				//this point has a higher index than the lowest index of the corresponding merge row
				r.createConstraintElement(mergeTmp->columnIndex, mergeTmp->value * data::QpNum(-1.0));
				if ((mergeTmp = mergeTmp->pNextRowPoint) == NULL)
					break;
			}
		}

	if (currTmp) {
		while (currTmp) {
			r.createConstraintElement(currTmp->columnIndex, currTmp->value);
			currTmp = currTmp->pNextRowPoint;
		}
	} else if (mergeTmp) {
		while (mergeTmp) {
			r.createConstraintElement(mergeTmp->columnIndex, (mergeTmp->value) * data::QpNum(-1.0));
			mergeTmp = mergeTmp->pNextRowPoint;
		}
	}
	//utils::Logger::globalLog(utils::LOG_INFO,LOG_TAG,"New Constraint after adding: "+r.toString());
	return r;

}

void Constraint::normalizeConstraint() {
	if (this->elements){
		if(this->pLastRowPoint->value!=1.0){
		this->devideConstraint(this->pLastRowPoint->value);
		}
	}
}


void Constraint::createConstraintElements(const std::vector<data::IndexedElement>& c){

	if(!c.size())return;

	bool ordered = true;
	for(unsigned int i = 1; i < c.size(); i++){
		if(c[i].index<c[i-1].index){
			ordered=false;
			break;
		}
	}

	if(!ordered){
	for(unsigned int i = 0; i < c.size(); i++){
		this->createConstraintElement(c[i].index,c[i].value);
	}
	return;
	}

	QlpMatCoeff *tmp,*nP = new QlpMatCoeff(c[0].value,*this, c[0].index);
	this->pFirstRowPoint=this->pLastRowPoint=nP;
	++elements;
	for(unsigned int i = 1; i < c.size();i++){
		nP = new QlpMatCoeff(c[i].value,*this, c[i].index);
		// the new point becomes the last row point
		pLastRowPoint->pNextRowPoint = nP;
		// that former last row point is the previous point of the new one
		nP->pPreviousRowPoint = pLastRowPoint;
		// point p is the last new QlpMatCoeff int this row
		nP->pNextRowPoint = NULL;
		// this is the last point that was added
		pLastRowPoint = nP;
		// establish connections to the other Points in this col
		if (qlp.getLastColumnsRowPtr(c[i].index) != NULL) {

			if ((nP->pPreviousColPoint = qlp.getLastColumnsPointPtr(c[i].index)) == NULL) {
							throw utils::QlpSolverException("Calling Method -> Create Point");
			}

//			if ((nP->pPreviousColPoint = qlp.getLastColumnsRowPtr(c[i].index)->getPointPtrByIndex(c[i].index)) == NULL) {
//				throw utils::QlpSolverException("Calling Method -> Create Point");
//			}
			nP->pPreviousColPoint->pNextColPoint = nP;
		} else {
			// add this point's column to be the first new row for a variable
			qlp.setFirstColumnsRowPtr(*nP);
			//qlp.setFirstColumnsRowPtr(nP->columnIndex, nP->pRow);
		}
		// add this point's column to be the last new column for a variable
		qlp.setLastColumnsRowPtr(*nP);
		//qlp.setLastColumnsRowPtr(nP->columnIndex, nP->pRow);
		++elements;
	}
}
//
void Constraint::createConstraintElements(const Constraint& c){

	if(c.isEmpty())return;

	QlpMatCoeff* tmp;
	tmp=c.pFirstRowPoint;
	QlpMatCoeff* nP = new QlpMatCoeff(tmp->value,*this, tmp->columnIndex);
	this->pFirstRowPoint=this->pLastRowPoint=nP;
	++elements;
	tmp=tmp->pNextRowPoint;
	while(tmp!=NULL){
		unsigned int index = tmp->columnIndex;
		nP = new QlpMatCoeff(tmp->value,*this, index);
		// the new point becomes the last row point
		pLastRowPoint->pNextRowPoint = nP;
		// that former last row point is the previous point of the new one
		nP->pPreviousRowPoint = pLastRowPoint;
		// point p is the last new QlpMatCoeff int this row
		nP->pNextRowPoint = NULL;
		// this is the last point that was added
		pLastRowPoint = nP;
		tmp=tmp->pNextRowPoint;
		// establish connections to the other Points in this col
		if (qlp.getLastColumnsRowPtr(index) != NULL) {
			if ((nP->pPreviousColPoint = qlp.getLastColumnsPointPtr(index)) == NULL) {
				throw utils::QlpSolverException("Calling Method -> Create Point");
			}
			nP->pPreviousColPoint->pNextColPoint = nP;
		} else {
			// add this point's column to be the first new row for a variable
			qlp.setFirstColumnsRowPtr(*nP);
					//qlp.setFirstColumnsRowPtr(nP->columnIndex, nP->pRow);
				}
		// add this point's column to be the last new column for a variable
		qlp.setLastColumnsRowPtr(*nP);
				//qlp.setLastColumnsRowPtr(nP->columnIndex, nP->pRow);

		++elements;
	}
}


//--------------- Operations on single elements of the left-hand side of a constraint ------------------->
void Constraint::createConstraintElement(unsigned int index, const data::QpNum& value) {

	// check if a row point for this index does already exist
	QlpMatCoeff& p = *new QlpMatCoeff(value, *this, index);
	QlpMatCoeff* tmpCurr;
	QlpMatCoeff* tmpNext;

	// create connections to other points in this row
	if (pFirstRowPoint == NULL && pLastRowPoint == NULL) {
		// this is the first, and thus, also the last point of this row
		pFirstRowPoint = pLastRowPoint = &p;
		// this row point neither has a next point
		p.pNextRowPoint = NULL;
		// nor has it a previous point
		p.pPreviousRowPoint = NULL;
	} else if (pFirstRowPoint->columnIndex > index) {
		// the first row point has a larger index than the new point, thus, the
		// new point becomes the first row element
		pFirstRowPoint->pPreviousRowPoint = &p;
		// that former first row point is the next point of the new one
		p.pNextRowPoint = pFirstRowPoint;
		// point p is the first new QlpMatCoeff in this row
		p.pPreviousRowPoint = NULL;
		// this is the first row point
		pFirstRowPoint = &p;
	} else if (pLastRowPoint->columnIndex < index) {
		// the new point has a higher index than the last row point, thus,
		// the new point becomes the last row point
		pLastRowPoint->pNextRowPoint = &p;
		// that former last row point is the previous point of the new one
		p.pPreviousRowPoint = pLastRowPoint;
		// point p is the last new QlpMatCoeff int this row
		p.pNextRowPoint = NULL;
		// this is the last point that was added
		pLastRowPoint = &p;
	} else {
		if (this->getPointPtrByIndex(index)) {
				addConstraintElement(index, value);
				return;
			}
		// Now we are in the messy situation that we have to look where this point
		// must be sorted in
		tmpCurr = pFirstRowPoint;
		tmpNext = tmpCurr->pNextRowPoint;
		do {
			if (tmpCurr->columnIndex < index && tmpNext->columnIndex > index) {
				tmpCurr->pNextRowPoint = &p;
				p.pPreviousRowPoint = tmpCurr;
				p.pNextRowPoint = tmpNext;
				tmpNext->pPreviousRowPoint = &p;
				break;
			}
			tmpCurr = tmpNext;
			tmpNext = tmpCurr->pNextRowPoint;
		} while (!tmpCurr->isLastRowElement());
	}

	// establish connections to the other Points in this col
	if (qlp.getLastColumnsRowPtr(index) != NULL) {
		if ((p.pPreviousColPoint = qlp.getLastColumnsPointPtr(index)) == NULL) {
			throw utils::QlpSolverException("Calling Method -> Create Point");

		}
		p.pPreviousColPoint->pNextColPoint = &p;
	} else {
		// add this point's column to be the first new row for a variable
		qlp.setFirstColumnsRowPtr(p);
				//qlp.setFirstColumnsRowPtr(nP->columnIndex, nP->pRow);

	}
	// add this point's column to be the last new column for a variable
	qlp.setLastColumnsRowPtr(p);
			//qlp.setLastColumnsRowPtr(nP->columnIndex, nP->pRow);

	++elements;
}

data::QpNum Constraint::getConstraintElementCoeff(unsigned int index) const {
	if (getPointPtrByIndex(index))
		return getPointPtrByIndex(index)->value;
	return data::QpNum(true);
}

void Constraint::deleteConstraintElement(unsigned int index) {
	QlpMatCoeff* tmp = getPointPtrByIndex(index);
	if (!tmp)
		throw utils::QlpSolverException(utils::ToolBox::convertToString(index));
	if (pFirstRowPoint->columnIndex == index)
		pFirstRowPoint = tmp->pNextRowPoint;
	if (pLastRowPoint->columnIndex == index)
		pLastRowPoint = tmp->pPreviousRowPoint;
	tmp->deleteFromRowAndCol();
	--elements;
}

void Constraint::addConstraintElement(unsigned int index, const data::QpNum& val) {
	if (getPointPtrByIndex(index))
		getPointPtrByIndex(index)->value += val;
	else
		throw utils::DataStructureException("Row: Row Element does not exist");
}

void Constraint::eliminateConstraintElement(unsigned int index, const data::QpNum& val) {
	QlpMatCoeff* p = getPointPtrByIndex(index);
	if (!p)
		throw utils::DataStructureException("Row: Row Element does not exist");
	// insert the value and normalize the inequality
	this->rhs -= (p->value * val);
	// delete the point form the entire qlp
	p->deleteFromRowAndCol();
	--elements;
}

void Constraint::eliminateConstraintElement(unsigned int index, const data::Constraint& c) {
	// get the point to remove
	QlpMatCoeff* p = getPointPtrByIndex(index);
	if (!p)
		throw utils::DataStructureException("Row: Row Element does not exist");

	// get the coeff since we have to multiply each constraint element
	// of the inserted constraint with this value
	const data::QpNum coeff = p->value;
	// delete the point from the qlp
	p->deleteFromRowAndCol();
	--elements;
	this->rhs -= (coeff * c.getRhsValue());
	const Constraint* r = &c;

	QlpMatCoeff* tmpNewPoint = r->pFirstRowPoint;
	QlpMatCoeff* tmpOldPoint;
	for (unsigned int i = 0; i < c.getElementCount() - 1; ++i) {

		// id the old point does not exist create a new one, otherwise modify the existing
		if (!(tmpOldPoint = getPointPtrByIndex(tmpNewPoint->columnIndex))) {
			createConstraintElement(tmpNewPoint->columnIndex, tmpNewPoint->value * coeff);
		} else {
			tmpOldPoint->value = tmpOldPoint->value - (tmpNewPoint->value * coeff);
			if (tmpOldPoint->value.isZero())
				tmpOldPoint->deleteFromRowAndCol();
		}
		tmpNewPoint = tmpNewPoint->pNextRowPoint;
	}
}

int Constraint::getFirstConstraintElementIndex() const {
	if (pFirstRowPoint)
		return pFirstRowPoint->columnIndex;
	return -1;
}

int Constraint::getLastConstraintElementIndex() const {
	if (pLastRowPoint)
		return pLastRowPoint->columnIndex;
	return -1;
}

int Constraint::getNextConstraintElementIndex(unsigned int index) const {
	QlpMatCoeff* p = this->getPointPtrByIndex(index);
	if (p && !p->isLastRowElement())
		return p->pNextRowPoint->columnIndex;
	return -1;
}

int Constraint::getPreviousConstraintElementIndex(unsigned int index) const {
	QlpMatCoeff* p = this->getPointPtrByIndex(index);
	if (p && !p->isFirstRowElement())
		return p->pPreviousRowPoint->columnIndex;
	return -1;
}

std::vector<data::IndexedElement> Constraint::getElements() const {
	std::vector<data::IndexedElement> elemVec;
	QlpMatCoeff* tmp = this->pFirstRowPoint;
	while (tmp) {
		elemVec.push_back(data::IndexedElement(tmp->columnIndex, tmp->value));
		tmp = tmp->pNextRowPoint;
	}
	return elemVec;
}

void Constraint::getElements(std::vector<data::QpNum>& v) const {
	QlpMatCoeff* tmp = this->pFirstRowPoint;
	while (tmp) {
		v[tmp->columnIndex] = tmp->value;
		tmp = tmp->pNextRowPoint;
	}
}

void Constraint::getElementsSparse(std::vector<data::IndexedElement>& v) const {
	v.clear();
	QlpMatCoeff* tmp = this->pFirstRowPoint;
	while (tmp) {
		v.push_back(data::IndexedElement(tmp->columnIndex,tmp->value));
		tmp = tmp->pNextRowPoint;
	}
}

void Constraint::setElements(const std::vector<data::IndexedElement>& vec) {
	for (unsigned int i = 0; i < vec.size(); i++)
		this->createConstraintElement(vec[i].index, vec[i].value);

}

std::string Constraint::toString() const {
	std::string begin(" ");
	if (!isEmpty()) {
		QlpMatCoeff* tmp = pFirstRowPoint;
		while (!tmp->isLastRowElement()) {
			if (tmp != pFirstRowPoint && tmp->value > 0.0)
				begin += " + ";
			begin += this->Point2String(tmp);
			tmp = tmp->pNextRowPoint;
		}
		if (tmp != pFirstRowPoint && tmp->value > 0.0)
			begin += " + ";
		begin += this->Point2String(pLastRowPoint);
		begin += " ";
	}
	begin += rhs.toString();
	data::QpRhs::ratioSignString(rhs.getRatioSign());
	return begin;
}

const std::string Constraint::Point2String(const QlpMatCoeff*p) const {
	std::string str1(" ");
	str1 += p->value.toString();
	str1 += this->qlp.getVariableNameByIndex(p->columnIndex);	//" * x";
	return str1;
}

std::string Constraint::getDoubleString() const {
	std::string begin(" ");
	if (!isEmpty()) {
		QlpMatCoeff* tmp = pFirstRowPoint;
		while (!tmp->isLastRowElement()) {
			if (tmp != pFirstRowPoint && tmp->value > 0.0)
				begin += " + ";
			begin += tmp->value.toString();
			begin += qlp.getVariableNameByIndex(tmp->columnIndex);
			tmp = tmp->pNextRowPoint;
		}
		if (tmp != pFirstRowPoint && tmp->value > 0.0)
			begin += " + ";
		begin += pLastRowPoint->value.toString();
		begin += qlp.getVariableNameByIndex(pLastRowPoint->columnIndex);
		begin += " ";
	}
	begin += data::QpRhs::ratioSignString(rhs.getRatioSign());
	begin += " ";
	begin += rhs.getValue().toString();
	return begin;
}

}
