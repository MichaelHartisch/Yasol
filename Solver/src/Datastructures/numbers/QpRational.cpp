/*
*
* Solver: QpRational.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Datastructures/numbers/QpRational.hpp"
#include "Datastructures/numbers/QpDouble.hpp"
#include <math.h>
#include "Utilities/ToolBox.hpp"
#include <sstream>
namespace data {
////--------------------------------Simple and Complex Constructors and Destrcro------------------
QpRational::QpRational(int n, int d) :
		value(n, d), status(Num) {
}

QpRational::QpRational(bool minInf) :
		value(0), status(minInf ? MinInf : MaxInf) {
}

QpRational::QpRational(const QpRational& f) :
		value(f.value), status(f.status) {
}

QpRational::QpRational(const double& d) :
		value(d), status(Num) {
}

QpRational::QpRational(const std::string& s) :
                value(0), status(Num) {
	    mpq_class n = utils::ToolBox::stringToMpqClass(s);
	    this->value = n;
	    this->status = Num;
}

//----------------Assignment Operators-----------------------------------------
QpRational& QpRational::operator =(const double& d) {
	this->value = d;
	this->status = Num;
	return *this;
}

QpRational& QpRational::operator =(const QpRational& f) {
	if (this == &f)
		return *this;
	this->value = f.value;
	this->status = f.status;
	return *this;
}

//----------------Binary Arithmetic Operators with two Fractions-------------------
QpRational QpRational::operator+(const QpRational& frac) const {
	return QpRational(*this) += frac;
}

QpRational QpRational::operator-(const QpRational& frac) const {
	return QpRational(*this) -= frac;
}

QpRational QpRational::operator*(const QpRational& frac) const {
	return QpRational(*this) *= frac;
}

QpRational QpRational::operator/(const QpRational& frac) const {
	return QpRational(*this) /= frac;
}

//-----------------Compound Assignment Operators with two Fractions-------------------------
QpRational& QpRational::operator +=(const data::QpRational& frac) {
	if (DEBUG_NUM && (this->status != Num || frac.status != Num))
		throw utils::DataStructureException("QpRational::operator += --> (this->status != Num || frac.status != Num)");
	this->value.operator +=(frac.value);
	return *this;
}

QpRational& QpRational::operator -=(const data::QpRational& frac) {
	if (DEBUG_NUM && (this->status != Num || frac.status != Num))
		throw utils::DataStructureException("QpRational::operator -=  --> (this->status != Num || frac.status != Num)");
	this->value.operator -=(frac.value);
	return *this;
}

QpRational& QpRational::operator /=(const data::QpRational& frac) {
	if (DEBUG_NUM && (this->status != Num || frac.status != Num))
		throw utils::DataStructureException("QpRational::operator /=  --> (this->status != Num || frac.status != Num)");
	this->value.operator /=(frac.value);
	return *this;
}

QpRational& QpRational::operator *=(const data::QpRational& frac) {
	if (DEBUG_NUM && (this->status != Num || frac.status != Num))
		throw utils::DataStructureException("QpRational::operator *=  --> (this->status != Num || frac.status != Num)");
	this->value.operator *=(frac.value);
	return *this;
}

//------------------Comparison Operators between two Fractions---------------------------
bool QpRational::operator ==(const QpRational& f) const {
	if (this == &f)
		return true;
	if (this->status == Num && f.status == Num) {
		return this->value == f.value;
	} else {
		return this->status == f.status;
	}
}

bool QpRational::operator !=(const QpRational& f) const {
	return !this->operator==(f);
}

bool QpRational::operator <(const QpRational& f) const {
	if (this == &f)
		return false;
	if (this->status == Num && f.status == Num)
		return this->value < f.value;
	if (this->status == MaxInf || f.status == MinInf)
		return false;
	if (this->status == MinInf || f.status == MaxInf)
		return true;
	throw utils::DataStructureException("bool QpRational::operator>(const QpRational& f) const");
}

bool QpRational::operator <=(const QpRational& f) const {
	return ((this->operator <(f)) || (this->operator ==(f)));
}

bool QpRational::operator>(const QpRational& f) const {
	if (this == &f)
		return false;
	if (this -> status == Num || f.status == Num)
		return this->value > f.value;
	if (this->status == MinInf || f.status == MaxInf)
		return false;
	if (this->status == MaxInf || f.status == MinInf)
		return true;
	throw utils::DataStructureException("bool QpRational::operator>(const QpRational& f) const: ");
}

bool QpRational::operator >=(const QpRational& f) const {
	return ((this->operator >(f)) || (this->operator ==(f)));
}

bool QpRational::isZero() const {
	if (status != Num)
		return false;
	return this->value == 0;
}

void QpRational::setZero() {
	this->value = 0;
	this->status = Num;
}

bool QpRational::isMaxInf() const {
	return (this->status == MaxInf);
}

bool QpRational::isMinInf() const {
	return (this->status == MinInf);
}

void QpRational::setMaxInf() {
	this->value = 0;
	this->status = MaxInf;
}

void QpRational::setMinInf() {
	this->value = 0;
	this->status = MinInf;
}

QpRational::Status QpRational::getStatus() const {
	return this->status;
}

double QpRational::asDouble() const {
	if (this->status == Num) {
		MyRatNum mrn(value);
		return ((double)mrn.getNominator() / (double)mrn.getDenominator());
	}
	if (this->status == MaxInf)
		return MAX_QPDOUBLE_INF;
	else if (this->status == MinInf)
		return MIN_QPDOUBLE_INF;
	throw utils::DataStructureException("double QpRational::asDouble() const --> unsupported status.");
}

//----------------------------Helper Functions--------------------------------------------
std::string QpRational::toString(unsigned int acc) const {
        std::string str1;
        if (this->status == MaxInf) {
                str1 += "+inf";
        } else if (this->status == MinInf) {
                str1 += "-inf";
        } else if (DISPLAY_RATIONAL) {
                std::stringstream temp1,temp2;
                MyRatNum mrn(value);
                temp1<<mrn.getNominator();
                str1 += temp1.str();
                str1 += "/";
                temp1<<mrn.getDenominator();
                str1 += temp1.str();
        } else {
                MyRatNum mrn(value);
                double d = (double)mrn.getNominator() / (double)mrn.getDenominator();
                str1 += utils::ToolBox::convertToString(d, TO_STRING_ACCURACY);
        }
        return str1;
}

}
