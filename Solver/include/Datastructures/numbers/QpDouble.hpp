/*
*
* Solver: QpDouble.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QPDOUBLE_HPP_
#define QPDOUBLE_HPP_
#include "Settings/Settings.hpp"
//#include <gmpxx.h>
#include "Datastructures/numbers/myRatNum.hpp"
namespace data {
class QpRational;
class QpDouble {
public:

	typedef enum {
		MinInf, Num, MaxInf
	} Status;

	//--------------------------------Simple and Complex Constructors and Destrcro------------------
	explicit QpDouble(int = 0, int = 1);
	explicit QpDouble(bool minInf);

	QpDouble(const double&);
	QpDouble(const QpDouble&);
	explicit QpDouble(const QpRational&);
	explicit QpDouble(const std::string&);

	QpDouble& operator =(const double&);
	QpDouble& operator =(const QpDouble&);
	QpDouble& operator =(const QpRational&);
	QpDouble& operator =(const std::string&);

	QpDouble operator +(const QpDouble&) const;
	QpDouble operator -(const QpDouble&) const;
	QpDouble operator *(const QpDouble&) const;
	QpDouble operator /(const QpDouble&) const;

	QpDouble& operator +=(const QpDouble&);
	QpDouble& operator -=(const QpDouble&);
	QpDouble& operator *=(const QpDouble&);
	QpDouble& operator /=(const QpDouble&);

	bool operator ==(const QpDouble&) const;
	bool operator !=(const QpDouble&) const;

	bool operator <(const QpDouble&) const;
	bool operator <=(const QpDouble&) const;
	bool operator >(const QpDouble&) const;
	bool operator >=(const QpDouble&) const;

	bool isZero() const;
	void setZero();
	bool isMaxInf() const;
	bool isMinInf() const;
	void setMaxInf();
	void setMinInf();

	Status getStatus()const;

	double asDouble() const;
	void cutDecimals(unsigned int prec);
	//operator double() { return this->asDouble(); }
	std::string toString(unsigned int acc=TO_STRING_ACCURACY) const;

protected:

private:
	double value;
	Status status;
};
}
#endif
/* FRACTION_HPP_ */
