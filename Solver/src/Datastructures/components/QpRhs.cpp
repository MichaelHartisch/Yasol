/*
*
* Solver: QpRhs.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Datastructures/components/QpRhs.hpp"
namespace data {

QpRhs::QpRhs(const QpNum& val, RatioSign ratio, data::QpRhs::Responsibility resp) :
    r(ratio), f(val), p(resp) {
}

QpRhs::QpRhs(double val, RatioSign ratio, int resp) :
    r(ratio), f(val), p((data::QpRhs::Responsibility)resp) {
}
    
QpRhs::QpRhs(double val, RatioSign ratio) :
    r(ratio), f(val) {
        p = data::QpRhs::UNDEF;
}
    
QpRhs::QpRhs(const QpNum& f, RatioSign ratio) :
	r(ratio), f(f) {
        p = data::QpRhs::UNDEF;
}

QpRhs::QpRhs(RatioSign ratio,double val) :
	r(ratio), f(val) {
        p = data::QpRhs::UNDEF;
}

QpRhs::QpRhs(RatioSign ratio,const QpNum& f) :
	r(ratio), f(f) {
        p = data::QpRhs::UNDEF;
}

QpRhs::QpRhs(const QpRhs& rhs) :
	r(rhs.r), f(rhs.f), p(rhs.p) {
}

QpRhs::~QpRhs() {
}

QpRhs& QpRhs::operator =(const QpRhs& rhs) {
	if (this == &rhs)
		return *this;
	this->r = rhs.r;
    this->f = rhs.f;
    this->p = rhs.p;
	return *this;
}

//-----------------Compound Assignment Operators with one QpNum and a Double--------------
QpRhs& QpRhs::operator +=(const QpNum& val) {
	this->f += val;
	return *this;
}

QpRhs& QpRhs::operator -=(const QpNum& val) {
	this->f -= val;
	return *this;
}

QpRhs& QpRhs::operator /=(const QpNum& val) {
	if (val.operator <(0.0))
		this->reverseRatioSign();
	this->f /= val;
	return *this;
}

QpRhs& QpRhs::operator *=(const QpNum& val) {
	if (val.operator <(0.0))
		this->reverseRatioSign();
	this->f *= val;
	return *this;
}

//------------------Comparison Operators between two Fractions---------------------------
bool QpRhs::operator ==(const QpRhs& rhs) const {
	return ((this->r == rhs.r) && (this->f == rhs.f));
}

bool QpRhs::operator !=(const QpRhs& rhs) const {
	return ((this->r != rhs.r) || (this->f != rhs.f));
}

void QpRhs::set(RatioSign ratio,const data::QpNum& f){
	r=ratio;
	this->f=f;
}

void QpRhs::setRatioSign(RatioSign ratio) {
	r = ratio;
}

const QpRhs::RatioSign& QpRhs::getRatioSign() const {
	return r;
}

void QpRhs::setValue(const QpNum& f) {
	this->f = f;
}

const QpNum& QpRhs::getValue() const {
	return this->f;
}

void QpRhs::setResponsibility(int p) {
    this->p = (data::QpRhs::Responsibility)p;
}
const int QpRhs::getResponsibility() const {
    return (int)p;
}

//----------------------------Helper Functions--------------------------------------------
std::string QpRhs::toString() const {
	std::string str1("");
	if (r == QpRhs::equal) {
		str1 += "=";
	} else if (r == QpRhs::smallerThanOrEqual) {
		str1 += "<=";
	} else {
		str1 += ">=";
	}
	return (str1 += this->f.toString());
}

std::string QpRhs::vecToString(const std::vector<data::QpRhs>& rhs) {
	std::string s("Rhs[");
	s+=utils::ToolBox::convertToString(rhs.size());
	s+="]:\n";
	for (unsigned int i = 0; i < rhs.size(); i++) {
		s += rhs[i].toString();
		s += "\n";
	}
	return s;
}

void QpRhs::reverseRatioSign() {
	if (this->r != QpRhs::equal)
		r == QpRhs::smallerThanOrEqual ? r = QpRhs::greaterThanOrEqual : r
				= QpRhs::smallerThanOrEqual;
}

std::string QpRhs::ratioSignString(const QpRhs::RatioSign& sign) {
	switch (sign) {
	case QpRhs::smallerThanOrEqual:
		return "<=";
	case QpRhs::greaterThanOrEqual:
		return ">=";
	case QpRhs::equal:
		return "=";
	default:
		throw "Unsupported Ratio Sign";
	}
}

QpRhs::RatioSign QpRhs::stringRatioSign(const std::string& sign) {
	if (sign == "<=")
		return QpRhs::smallerThanOrEqual;
	else if (sign == ">=")
		return QpRhs::greaterThanOrEqual;
	else if (sign == "=")
		return QpRhs::equal;
	else
		throw "Unsupported String";
}
}
