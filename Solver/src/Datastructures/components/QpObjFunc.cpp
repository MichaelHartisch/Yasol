/*
*
* Solver: QpObjFunc.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Datastructures/components/QpObjFunc.hpp"
namespace data {

QpObjFunc::QpObjFunc(Objective o, data::QpNum f) :
	obj(o), offset(f), coeffs() {
}

QpObjFunc::QpObjFunc(QpObjFunc::Objective obj,
		const std::vector<QpNum>& elemVec, data::QpNum f) :
	obj(obj), offset(f), coeffs(elemVec) {
}

QpObjFunc::QpObjFunc(QpObjFunc::Objective obj, const std::vector<
		data::IndexedElement>&elemVec, unsigned int size, data::QpNum f) :
	obj(obj), offset(f), coeffs(size, 0) {
	for (unsigned int i = 0; i < elemVec.size(); i++) {
		if (!elemVec[i].value.isZero()){
			this->coeffs[elemVec[i].index] = elemVec[i].value;
		}
	}
}

QpObjFunc::~QpObjFunc() {
}

QpObjFunc::QpObjFunc(const QpObjFunc& o) :
	obj(o.obj), offset(o.offset), coeffs(o.coeffs) {
}

QpObjFunc& QpObjFunc::operator=(const QpObjFunc& o) {
	if (this == &o)
		return *this;
	this->obj = o.obj;
	this->offset = o.offset;
	this->coeffs = o.coeffs;
	return *this;
}

const QpNum& QpObjFunc::operator [](unsigned int i) const {
	if (i >= this->coeffs.size()){
			throw utils::DataStructureException(
					"QpObjFunc::operator [](unsigned int i) --> (i >= this->coeffs.size())");
	}
	return this->coeffs[i];
}

QpNum QpObjFunc::getObjFuncValue(const std::vector<data::QpNum>& vec) const {
	QpNum num;
	if (vec.size() != this->coeffs.size()){
		throw utils::DataStructureException(
				"QpObjFunc () --> vec.size()!=coeffs.size()");
	}
	for (unsigned int i = 0, size = vec.size(); i < size; i++) {
		num += vec[i] * coeffs[i];
	}
	return num;
}

QpNum QpObjFunc::getObjFuncValue(const std::vector<data::IndexedElement>& ieVec) const {
	QpNum num;
	unsigned int index = 0;
	for (unsigned int i = 0, size = ieVec.size(); i < size; i++) {
		if((index = ieVec[i].index)>=coeffs.size()){
			throw utils::DataStructureException(
							"QpObjFunc () --> ieVec[i].index>=coeffs.size()");
		}
		num += ieVec[i].value * coeffs[index];
	}
	return num;
}

void QpObjFunc::setObjective(Objective o) {
	this->obj = o;
}
QpObjFunc::Objective QpObjFunc::getObjective() const {
	return this->obj;
}

void QpObjFunc::setOffset(const data::QpNum& f) {
	this->offset = f;
}

const data::QpNum& QpObjFunc::getOffset() const {
	return this->offset;
}

bool QpObjFunc::operator ==(const QpObjFunc& o) const {
	if (this->obj != o.obj)
		return false;
	if (this->offset != o.offset)
		return false;
	if (this->coeffs != o.coeffs)
		return false;
	return true;
}

bool QpObjFunc::operator !=(const QpObjFunc& o) const {
	return !this->operator ==(o);
}

void QpObjFunc::reverseObjectiveFuntion() {
	obj == min ? obj = max : obj = min;
	this->offset *= -1.0;
	for (unsigned int i = 0; i < this->coeffs.size(); i++) {
		if (!this->coeffs[i].isZero()){
			this->coeffs[i]*=-1.0;
		}
	}
}

void QpObjFunc::setObjElement(unsigned int index, const QpNum& value) {
	if (index >= this->coeffs.size()){
		throw utils::DataStructureException("QpObjFunc::setObjElement(unsigned int index, const QpNum& value) --> (index >= this->coeffs.size())");
	}
	this->coeffs[index] = value;
}

void QpObjFunc::setObjElements(const std::vector<QpNum>& elems){
	this->coeffs = elems;
}

void QpObjFunc::setObjElements(const std::vector<data::IndexedElement>& elems, unsigned int size){
	this->coeffs.clear();
	this->coeffs.resize(size);
	for(unsigned int i = 0; i < elems.size(); i++){
		this->setObjElement(elems[i].index,elems[i].value);
	}
}


std::vector<data::IndexedElement> QpObjFunc::getObjectiveElementsSparse() const {
	std::vector<data::IndexedElement> vec;
	for (unsigned int i = 0, size = this->coeffs.size(); i < size; i++) {
		if (!coeffs[i].isZero()){
			vec.push_back(data::IndexedElement(i, coeffs[i]));
		}
	}
	return vec;
}

const std::vector<QpNum>& QpObjFunc::getObjectiveElementsDense() const {
	return this->coeffs;
}

void QpObjFunc::setSize(unsigned int size) {
	this->coeffs.resize(size, 0);
}

unsigned int QpObjFunc::getSize() const {
	return this->coeffs.size();
}

void QpObjFunc::clearElements() {
	for (unsigned int i = 0; i < coeffs.size(); i++) {
		coeffs[i].setZero();
	}
}

void QpObjFunc::clear() {
	this->obj = min;
	this->offset.setZero();
	this->clearElements();
}

std::string QpObjFunc::toString() const {
	std::string s("Obj[");
	s += utils::ToolBox::convertToString(this->coeffs.size());
	s += "] --> ";
	this->obj == max ? s += "max " : s += "min ";
	if (!offset.isZero())
		s += offset.toString();
	s += data::QpNum::vecToString(this->coeffs);
	return s;
}

std::string QpObjFunc::objectiveString(QpObjFunc::Objective obj) {
	switch (obj) {
	case QpObjFunc::min:
		return "MINIMIZE";
	case QpObjFunc::max:
		return "MAXIMIZE";
	default:
		throw "Unsupported Objective";
	}
}

}
