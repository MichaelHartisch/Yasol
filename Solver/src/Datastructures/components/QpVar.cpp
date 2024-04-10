/*
*
* Solver: QpVar.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Datastructures/components/QpVar.hpp"
namespace data {

QpVar::QpVar():name("no_name"), index(0), vRange(), vDistr(), q(exists), n(real){
	vRange.push_back(data::QpNum(true));
	vRange.push_back(data::QpNum(false));
}

QpVar::QpVar(const std::string& n, unsigned int i, const QpNum& l,
		const QpNum& u, NumberSystem ns, Quantifier q) :
		name(n), index(i), vRange(), vDistr(), q(q), n(ns) {
	vRange.push_back(l);
	vRange.push_back(u);
}

QpVar::QpVar(const std::string& n, unsigned int i,
		const std::vector<data::QpNum>& vr,
		const std::vector<data::QpRational>& vd, NumberSystem ns) :
		name(n), index(i), vRange(), vDistr(), q(random), n(ns) {
	this->setVariableRange(vr, vd);
}

QpVar::~QpVar() {
}

QpVar::QpVar(const QpVar& qvar) :
		name(qvar.name), index(qvar.index), vRange(qvar.vRange), vDistr(
				qvar.vDistr), q(qvar.q), n(qvar.n) {
}

QpVar& QpVar::operator=(const QpVar& qvar) {
	if (this == &qvar)
		return *this;
	this->name = qvar.name;
	this->index = qvar.index;
	this->vRange = qvar.vRange;
	this->vDistr = qvar.vDistr;
	this->q = qvar.q;
	return *this;
}

void QpVar::setName(std::string& n) {
	this->name = n;
}

const std::string& QpVar::getName() const {
	return this->name;
}

void QpVar::setIndex(unsigned int i) {
	this->index = i;
}

unsigned int QpVar::getIndex() const {
	return this->index;
}

void QpVar::setLowerBound(const QpNum& lb) {
	if (this->n == data::QpVar::binaries && (lb < 0 || lb > 1))
		throw utils::DataStructureException(
				"void QpVar::setLowerBound(const QpNum& lb) --> only [0,1] bounds allowed for binary variables.");
	vRange[0] = lb;
}

const QpNum& QpVar::getLowerBound() const {
	return vRange[0];
}

void QpVar::setUpperBound(const QpNum& ub) {
	if (this->n == data::QpVar::binaries && (ub < 0 || ub > 1))
		throw utils::DataStructureException(
				"void QpVar::setLowerBound(const QpNum& ub) --> only [0,1] bounds allowed for binary variables.");
	vRange[vRange.size() - 1] = ub;
}

const QpNum& QpVar::getUpperBound() const {
	return vRange[vRange.size() - 1];
}

void QpVar::setBounds(const QpNum& lb, const QpNum& ub) {
	if (vRange.size() > 2)
		vRange.resize(2, 0);
	this->setLowerBound(lb);
	this->setUpperBound(ub);
}

std::pair<QpNum,QpNum> QpVar::getBounds() const {
	return std::make_pair(vRange[0], vRange[vRange.size() - 1]);
}

void QpVar::setQuantifier(Quantifier q) {
	this->q = q;
}

QpVar::Quantifier QpVar::getQuantifier() const {
	return this->q;
}

void QpVar::setNumberType(NumberSystem n) {
	this->n = n;
	if(n==data::QpVar::binaries){
		if(this->vRange[0]<0)this->vRange[0].setZero();
		if(this->vRange[vRange.size()-1]<0)this->vRange[vRange.size()-1]=1;
	}
}
QpVar::NumberSystem QpVar::getNumberSystem() const {
	return this->n;
}

bool QpVar::operator ==(const QpVar& qvar) const {
	return ((this->index == qvar.index) && (this->vRange == qvar.vRange)
			&& (this->vDistr == qvar.vDistr) && (this->q == qvar.q)
			&& (this->n == qvar.n));
}

bool QpVar::operator !=(const QpVar& qvar) const {
	return !this->operator ==(qvar);
}

bool QpVar::isFree() const {
	return (!this->isBoundedBelow() && !this->isBoundedBelow());
}

bool QpVar::isBoundedAbove() const {
	return !vRange[vRange.size() - 1].isMaxInf();
}

bool QpVar::isBoundedBelow() const {
	return !vRange[0].isMinInf();
}

void QpVar::setFree() {
	if(this->n==data::QpVar::binaries)
		throw utils::DataStructureException("QpVar::setUnboundedBelow() --> not allowed for binary variables");
	this->setUnboundedBelow();
	this->setUnboundedAbove();
}

void QpVar::setUnboundedAbove() {
	if(this->n==data::QpVar::binaries)
		throw utils::DataStructureException("QpVar::setUnboundedAbove() --> not allowed for binary variables");
	vRange[vRange.size() - 1].setMaxInf();
}

void QpVar::setUnboundedBelow() {
	if(this->n==data::QpVar::binaries)
		throw utils::DataStructureException("QpVar::setUnboundedBelow() --> not allowed for binary variables");
	vRange[0].setMinInf();
}

void QpVar::setVariableRange(const std::vector<data::QpNum>& range,
		const std::vector<data::QpRational>& distr) {
	this->q = QpVar::random;
	this->vRange = range;
	this->vDistr = distr;
}

const std::vector<data::QpNum>& QpVar::getVariableRange() const {
	return this->vRange;
}

const std::vector<data::QpRational>& QpVar::getVariableDistribution() const {
	return this->vDistr;
}

std::string QpVar::toString() const {
	std::string s(name);
	s += " (i = ";
	s += utils::ToolBox::convertToString(this->index);
	s += " )[";
	s += data::QpNum::vecToString(vRange);
	s += "," + data::QpNum::vecToString(vDistr);
	s += "] ";
	if (this->q == exists) {
		s += "[exists,";
	} else if (this->q == all) {
		s += "[all,";
	} else {
		s += "[random,";
	}

	if (this->n == data::QpVar::real) {
		s += "real";
	} else if (this->n == data::QpVar::generals) {
		s += "general";
	} else {
		s += "binary";
	}
	if(this->vRange[0]==this->vRange[this->vRange.size()-1])
		s+=",fixed";
	s+="]";

	return s;
}

std::string QpVar::vecToString(const std::vector<data::QpVar>& vars) {
	std::string s("Variables [");
	s += utils::ToolBox::convertToString(vars.size());
	s += "]:\n";
	for (unsigned int i = 0; i < vars.size(); i++) {
		s += "\t\t" + vars[i].toString();
		s += "\n";
	}
	return s;
}

std::string QpVar::numberString(const QpVar::NumberSystem& number) {
	switch (number) {
	case QpVar::real:
		return "REALS   ";
	case QpVar::generals:
		return "GENERALS";
	case QpVar::binaries:
		return "BINARIES";
	default:
		throw "Unsupported Number Format";
	}
}

QpVar::NumberSystem QpVar::stringNumber(const std::string& number) {
	if (number == "REALS")
		return QpVar::real;
	else if (number == "GENERALS")
		return QpVar::generals;
	else if (number == "BINARIES")
		return QpVar::binaries;
	else
		throw "Unsupported String";
}

std::string QpVar::quantifierString(const QpVar::Quantifier& quant) {
	switch (quant) {
	case QpVar::exists:
		return "EXISTS";
	case QpVar::all:
		return "ALL   ";
	case QpVar::random:
		return "RANDOM";
	default:
		throw "Unsupported Quantifier";
	}
}

QpVar::Quantifier QpVar::stringQuantifier(const std::string& quant) {
	if (quant == "EXISTS")
		return QpVar::exists;
	else if (quant == "ALL   ")
		return QpVar::all;
	else if (quant == "RANDOM")
		return QpVar::random;
	else
		throw "Unsupported String";
}

}
