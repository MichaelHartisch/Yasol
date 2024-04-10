/*
*
* Solver: QlpMatCoeff.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef POINT_HPP_
#define POINT_HPP_

//#include "Settings.hpp"
#include "Datastructures/numbers/QpNum.hpp"

namespace data {

class Constraint;

class QlpMatCoeff {

	friend class QlpConverter;
	friend class Constraint;
	friend class Qlp;

protected:

	QlpMatCoeff(const data::QpNum& f, Constraint& row, unsigned int col) :
		value(f), columnIndex(col), pRow(&row), pNextRowPoint(NULL), pPreviousRowPoint(
				NULL), pNextColPoint(NULL), pPreviousColPoint(NULL) {
	}
	~QlpMatCoeff(){}

	// Decouples this point from points  it is connected to in its row.
	// The binding to points from other rows are hold
	void deleteFromRow();
	// Decouples this point from points in different cols it is connected to.
	// The binding to points within the row are hold.
	void deleteFromCol();
	// Deletes the point from the matrix and adapts the pointers
	// of all points it is connected to, to ensure consistency
	void deleteFromRowAndCol();

	inline const bool isFirstRowElement() const{
		return ((pPreviousRowPoint == NULL) ? true : false);
	}
	inline const bool isLastRowElement() const {
		return ((pNextRowPoint == NULL) ? true : false);
	}
	inline const bool isFirstColElement() const {
		return ((pPreviousColPoint == NULL) ? true : false);
	}
	inline const bool isLastColElement() const {
		return ((pNextColPoint == NULL) ? true : false);
	}
	inline const bool equals(const QlpMatCoeff& p) const {
		return (value == p.value && columnIndex == p.columnIndex);
	}
	//Hidden and not implemented
	QlpMatCoeff(const QlpMatCoeff&);
	QlpMatCoeff& operator=(const QlpMatCoeff&);
	bool operator==(const QlpMatCoeff&) const;

	// Coefficient
	data::QpNum value;
	// Index of the variable this QlpMatCoeff object corresponds to
	unsigned int columnIndex;
	// Pointer to the Row object that holds information about the restriction
	Constraint *pRow;
	// Pointer to the next entry  a_ij+* in the corresponding row that has a coefficient not equal to zero
	QlpMatCoeff *pNextRowPoint;
	// Pointer to the previous entry a_ij-* in the corresponding row that has a coefficient not equal to zero
	QlpMatCoeff *pPreviousRowPoint;
	// Pointer to the point in a next row where the point has a coefficient not equal to zero (a_i+*j != 0)
	QlpMatCoeff *pNextColPoint;
	// Pointer to the point in a previous row where the point has an coefficient not equal to zero (a_i-*j != 0)
	QlpMatCoeff *pPreviousColPoint;
};
}
#endif /* POINT_HPP_ */
