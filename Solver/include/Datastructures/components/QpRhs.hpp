/*
*
* Solver: QpRhs.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef RHS_HPP_
#define RHS_HPP_
#include "Datastructures/numbers/Numbers.hpp"
namespace data {
class QpRhs {
public:

    typedef enum {
        smallerThanOrEqual, equal, greaterThanOrEqual
    } RatioSign;
    typedef enum {
        EXISTENTIAL, UNIVERSAL, UNDEF
    } Responsibility;

    QpRhs(const QpNum& val, RatioSign ratio, data::QpRhs::Responsibility resp);
    QpRhs(double val, RatioSign ratio, int resp);
	QpRhs(double, RatioSign);
	QpRhs(const QpNum&, RatioSign);
	QpRhs(RatioSign = smallerThanOrEqual, double=0);
	QpRhs(RatioSign,const QpNum&);
	~QpRhs();

	QpRhs(const QpRhs&);
	QpRhs& operator =(const QpRhs&);

	QpRhs& operator +=(const QpNum&);
	QpRhs& operator -=(const QpNum&);
	QpRhs& operator *=(const QpNum&);
	QpRhs& operator /=(const QpNum&);

	bool operator ==(const QpRhs&) const;
	bool operator !=(const QpRhs&) const;

	void set(RatioSign ratio,const data::QpNum& f);
	void setRatioSign(RatioSign ratio);
	const RatioSign& getRatioSign() const;
	void reverseRatioSign();

    void setValue(const QpNum& f);
    const QpNum& getValue() const;
    void setResponsibility(int p);
    const int getResponsibility() const;

	std::string toString() const;
	static std::string vecToString(const std::vector<data::QpRhs>&);

	static std::string ratioSignString(const QpRhs::RatioSign& sign);
	static RatioSign stringRatioSign(const std::string& sign);

protected:
private:
	RatioSign r;
	QpNum f;
    Responsibility p;
};
}
#endif /*RHS_HPP_*/
