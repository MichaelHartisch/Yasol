/*
*
* Solver: QpObjFunc.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef OBJECTIVEFUNCTION_HPP_
#define OBJECTIVEFUNCTION_HPP_
//#include "Datastructures/Datastructures.hpp"
#include "Datastructures/numbers/Numbers.hpp"

namespace data {
class QpObjFunc {
public:

	typedef enum {
		max, min
	} Objective;

	//constructors and destructor
	QpObjFunc(Objective=min,data::QpNum=0.0);
	QpObjFunc(Objective, const std::vector<QpNum>&,data::QpNum=0.0);
	QpObjFunc(Objective, const std::vector<data::IndexedElement>&, unsigned int,data::QpNum=0.0);
	~QpObjFunc();

	//copy constructor and assignment operator
	QpObjFunc(const QpObjFunc&);
	QpObjFunc& operator=(const QpObjFunc&);

	//get coefficient at position i
	const QpNum& operator [](unsigned int)const;

	//get objective function value for a given variable assignment
	QpNum getObjFuncValue(const std::vector<data::QpNum>&)const;
	QpNum getObjFuncValue(const std::vector<data::IndexedElement>&)const;

	//comparison operators
	bool operator ==(const QpObjFunc&) const;
	bool operator !=(const QpObjFunc&) const;

	//setter and getter for objective
	void setObjective(Objective);
	Objective getObjective() const;

	//setter and getter for offset
	void setOffset(const data::QpNum&);
	const data::QpNum& getOffset() const;

	//getter for coefficients
	std::vector<data::IndexedElement> getObjectiveElementsSparse() const;
	const std::vector<QpNum>& getObjectiveElementsDense() const;

	//setter for coefficients
	void setObjElement(unsigned int, const QpNum&);
	void setObjElements(const std::vector<data::IndexedElement>&, unsigned int);
	void setObjElements(const std::vector<QpNum>&);

	//setter and getter for size
	void setSize(unsigned int);
	unsigned int getSize() const;

	//reverse objective, offset and coefficients multiplied by -1.0
	void reverseObjectiveFuntion();
	//all objective function coefficients are set to zero
	void clearElements();
	//objective is set to min, offset to zero, all coefficients to zero
	void clear();

	std::string toString() const;
	static std::string objectiveString(QpObjFunc::Objective);

protected:
private:
	Objective obj;
	data::QpNum offset;
	std::vector<data::QpNum> coeffs;
};

}
#endif /*OBJECTIVEFUNCTION_HPP_*/
