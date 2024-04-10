/*
*
* Solver: QpRational.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QPRATIONAL_HPP_
#define QPRATIONAL_HPP_
#include "Settings/Settings.hpp"
//#include <gmpxx.h>
#include "Datastructures/numbers/myRatNum.hpp"
namespace data {
class QpDouble;
class QpRational {
public:

	typedef enum {
		MinInf, Num, MaxInf
	} Status;

	explicit QpRational(int = 0, int = 1);
	explicit QpRational(bool minInf);

	QpRational(const double&);
	QpRational(const QpRational&);
	QpRational(const QpDouble&);
	QpRational(const mpq_class&);
	QpRational(const std::string&);

	QpRational& operator =(const double&);
	QpRational& operator =(const mpq_class&);
	QpRational& operator =(const QpRational&);
	QpRational& operator =(const QpDouble&);
	QpRational& operator =(const std::string&);

	QpRational operator +(const QpRational&) const;
	QpRational operator -(const QpRational&) const;
	QpRational operator *(const QpRational&) const;
	QpRational operator /(const QpRational&) const;

	QpRational& operator +=(const QpRational&);
	QpRational& operator -=(const QpRational&);
	QpRational& operator *=(const QpRational&);
	QpRational& operator /=(const QpRational&);

	bool operator ==(const QpRational&) const;
	bool operator !=(const QpRational&) const;

	bool operator <(const QpRational&) const;
	bool operator <=(const QpRational&) const;
	bool operator >(const QpRational&) const;
	bool operator >=(const QpRational&) const;

	bool isZero() const;
	void setZero();
	bool isMaxInf() const;
	void setMaxInf();
	bool isMinInf() const;
	void setMinInf();

	Status getStatus()const;
	double asDouble() const;
	void cutDecimals(unsigned int){}
	mpq_class asMpqClass() const;

	//operator double() { return this->asDouble(); }

	std::string toString(unsigned int acc=TO_STRING_ACCURACY) const;

private:
	mpq_class value;
	Status status;
};
}
#endif
/* QPRATIONAL_HPP_ */
