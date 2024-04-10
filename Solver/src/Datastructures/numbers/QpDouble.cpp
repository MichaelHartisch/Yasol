/*
*
* Solver: QpDouble.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Datastructures/numbers/QpDouble.hpp"
#include "Datastructures/numbers/QpRational.hpp"
#include <math.h>
#include "Utilities/ToolBox.hpp"
//#include <boost/lexical_cast.hpp>
namespace data {
////--------------------------------Simple and Complex Constructors and Destrcro------------------
QpDouble::QpDouble(int n, int d) :
		value(((double) n) / d), status(Num) {
}

QpDouble::QpDouble(bool minInf) :
		value(0), status(minInf ? MinInf : MaxInf) {
}

QpDouble::QpDouble(const QpDouble& val) :
		value(val.value), status(val.status) {
}

QpDouble::QpDouble(const QpRational& val) :
		value(val.asDouble()), status((QpDouble::Status) val.getStatus()) {
}

QpDouble::QpDouble(const double& val) :
		value(val), status(Num) {
}

QpDouble::QpDouble(const std::string& s) :
                value(0), status(Num) {
	    mpq_class n = utils::ToolBox::stringToMpqClass(s);
	    this->value = (double)n.getNominator() / (double)n.getDenominator();
	    this->status = Num;
        /*RatInf<mpq_class> n = utils::ToolBox::stringToTORationalInf(s);
        if(n.isInf){
                n.value.getNominator() < 0 ? status = MinInf : status = MaxInf ;
        }else{
                this->value = (double)n.value.getNominator() / (double)n.value.getDenominator();
        }*/
}

//----------------Assignment Operators-----------------------------------------
QpDouble& QpDouble::operator =(const double& val) {
	this->value = val;
	this->status = Num;
	return *this;
}

QpDouble& QpDouble::operator =(const QpDouble& val) {
	this->value = val.value;
	this->status = val.status;
	return *this;
}

QpDouble& QpDouble::operator =(const QpRational& val) {
	this->value = val.asDouble();
	this->status = (QpDouble::Status) val.getStatus();
	return *this;
}

QpDouble& QpDouble::operator =(const std::string& s) {
        MyRatNum mrn(utils::ToolBox::stringToMpqClass(s));
        value = (double)mrn.getNominator() / (double)mrn.getDenominator();
        this->status = Num;
        return *this;
}

//----------------Binary Arithmetic Operators with two Fractions-------------------
QpDouble QpDouble::operator+(const QpDouble& frac) const {
	return QpDouble(*this) += frac;
}

QpDouble QpDouble::operator-(const QpDouble& frac) const {
	return QpDouble(*this) -= frac;
}

QpDouble QpDouble::operator*(const QpDouble& frac) const {
	return QpDouble(*this) *= frac;
}

QpDouble QpDouble::operator/(const QpDouble& frac) const {
	return QpDouble(*this) /= frac;
}

//-----------------Compound Assignment Operators with two Fractions-------------------------
QpDouble& QpDouble::operator +=(const data::QpDouble& val) {
	if (DEBUG_NUM && (this->status != Num || val.status != Num)){
		throw utils::DataStructureException("QpDouble::operator +=(const data::QpDouble& val) --> (this->status != Num || frac.status != Num)");
	}
	this->value += (val.value);
	return *this;
}

QpDouble& QpDouble::operator -=(const data::QpDouble& val) {
	if (DEBUG_NUM && (this->status != Num || val.status != Num))
		throw utils::DataStructureException("QpDouble::operator -=(const data::QpDouble& val) --> (this->status != Num || frac.status != Num)");
	this->value -= (val.value);
	return *this;
}

QpDouble& QpDouble::operator /=(const data::QpDouble& val) {
	if (DEBUG_NUM && (this->status != Num || val.status != Num))
		throw utils::DataStructureException("QpDouble::operator /=(const data::QpDouble& val) --> (this->status != Num || frac.status != Num)");
	this->value /= (val.value);
	return *this;
}

QpDouble& QpDouble::operator *=(const data::QpDouble& val) {
	if (DEBUG_NUM && ((this->status != Num) || (val.status != Num)))
		throw utils::DataStructureException("QpDouble::operator *=(const data::QpDouble& val) --> (this->status != Num || frac.status != Num)");
	this->value *= (val.value);
	return *this;
}

//------------------Comparison Operators between two Fractions---------------------------
bool QpDouble::operator ==(const QpDouble& f) const {
	if (this == &f)
		return true;
	if (this->status == Num && f.status == Num) {
		if (!NUM_CAUT)
			return this->value == f.value;
		return fabs(this->value - f.value) < DOUBLE_EPSILON;
	} else {
		return this->status == f.status;
	}
}

bool QpDouble::operator !=(const QpDouble& f) const {
	return !this->operator ==(f);
}

bool QpDouble::operator <(const QpDouble& f) const {
	if (this == &f)
		return false;
	if (this->status == Num && f.status == Num) {
		if (!NUM_CAUT)
			return this->value < f.value;
		return (f.value - this->value) > DOUBLE_EPSILON; //*((std::fabs(this->value) < std::fabs(f.value) ? std::fabs(f.value) : std::fabs(this->value)) * DOUBLE_EPSILON);
	}
	if (this->status == MaxInf || f.status == MinInf)
		return false;
	if (this->status == MinInf || f.status == MaxInf)
		return true;
	throw utils::DataStructureException("bool QpDouble::operator <(const QpDouble& f) const --> Unsupported Comparison");
}

bool QpDouble::operator <=(const QpDouble& f) const {
	return ((this->operator <(f.value)) || (this->operator ==(f.value)));
}

bool QpDouble::operator>(const QpDouble& f) const {
	if (this == &f)
		return false;
	if (this->status == Num && f.status == Num) {
		if (!NUM_CAUT)
			return this->value > f.value;
		return (this->value - f.value) > DOUBLE_EPSILON; //((std::fabs(this->value) < std::fabs(f.value) ? std::fabs(f.value) : std::fabs(this->value)) * DOUBLE_EPSILON);
	}
	if (this->status == MinInf || f.status == MaxInf)
		return false;
	if (this->status == MaxInf || f.status == MinInf)
		return true;
	throw utils::DataStructureException("bool QpDouble::operator <(const QpDouble& f) const --> Unsupported Comparison");
}

bool QpDouble::operator >=(const QpDouble& f) const {
	return ((this->operator >(f)) || (this->operator ==(f)));
}

bool QpDouble::isZero() const {
	if (status != Num)
		return false;
	if (!NUM_CAUT)
		return this->value == 0;
	return fabs(this->value) < DOUBLE_EPSILON;
}

void QpDouble::setZero() {
	this->status = Num;
	this->value = 0;
}

bool QpDouble::isMaxInf() const {
	return (this->status == MaxInf);
}

bool QpDouble::isMinInf() const {
	return (this->status == MinInf);
}

void QpDouble::setMaxInf() {
	this->status = MaxInf;
	this->value = 0;
}

void QpDouble::setMinInf() {
	this->status = MinInf;
	this->value = 0;
}

QpDouble::Status QpDouble::getStatus() const {
	return this->status;
}

double QpDouble::asDouble() const {
	if (this->status == Num)
		return this->value;
	if (this->status == MaxInf)
		return MAX_QPDOUBLE_INF;
	if (this->status == MinInf)
		return MIN_QPDOUBLE_INF;
	throw utils::DataStructureException("double QpDouble::asDouble() const --> Unsupported Status.");
}

void QpDouble::cutDecimals(unsigned int prec){
	 double v[] = { 1, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e12, 1e13, 1e14, 1e15 };  // mgl. verl�ngern
	 value = floor(value * v[prec] + 0.5) / v[prec];
}

//----------------------------Helper Functions--------------------------------------------
std::string QpDouble::toString(unsigned int acc) const {
	std::string str1("");
	if (this->status == MaxInf) {
		str1 += "+inf";
	} else if (this->status == MinInf) {
		str1 += "-inf";
	} else {
		str1 += utils::ToolBox::convertToString(this->value, acc);
	}
	return str1;
}
}
