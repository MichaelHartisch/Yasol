/*
*
* Solver: Constraint.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef CONSTRAINT_HPP_
#define CONSTRAINT_HPP_

//#include "Settings.hpp"
#include "Datastructures/components/QpRhs.hpp"
#include "Datastructures/qlp/QlpMatCoeff.hpp"

namespace data {

class Qlp;

class Constraint {

	friend class QlpConverter;
	friend class QlpMatCoeff;
	friend class Qlp;

public:
    typedef enum {
        EXISTENTIAL, UNIVERSAL, UNDEF
    } Responsibility;

    // Returns if this Constraints has any elements on the left side
	bool isEmpty() const;
	//void setElementCount(unsigned int);
	unsigned int getElementCount() const;
	// Returns id this is the first row of the qlp
	bool isFirstConstraint() const;
	// Returns iof this is the last row of the qlp
	bool isLastConstraint() const;
	//Compares two Constraints with respect to their values, not their references
	bool operator==(const Constraint&) const;
    //Returns the responsibility for the constraint
    const int getResponsibility() const;
    void setResponsibility(int p);

    // Set Rhs value
	void setRhsValue(const QpNum&);
	// Get Rhs Value
	const QpNum& getRhsValue() const;
	// Set Ratio Sign
	void setRhsRatioSign(QpRhs::RatioSign);
	// Get Ratio Sign
	QpRhs::RatioSign getRhsRatioSign() const;
	// Reverse the Ratio Sign (= stays unchanged)
	void reverseRhsRatioSign();
	// Get reference to the rhs object
	const data::QpRhs& getRhs() const;

    // Get pointer to previous constraint or NULL if this is the first one
	inline Constraint* getPreviousConstraintPtr() {
		return previousConstraint;
	}
	// Get pointer to constant previous constraint or NULL if this is the first one
	inline const Constraint* getPreviousConstraintPtr() const {
		return previousConstraint;
	}
	// Get pointer to previous constraint or NULL if this is the last one
	inline const Constraint* getNextConstraintPtr() const {
		return nextConstraint;
	}
	// Get pointer to next constraint or NULL if this is the last one
	inline Constraint* getNextConstraintPtr() {
		return nextConstraint;
	}

	// Removes this constraint from
	void removeConstraint();
	// Multiplies the entire constraint
	void multiplyConstraint(const data::QpNum&);
	// Divides the entire constraint
	void devideConstraint(const data::QpNum&);
	// Creates a new constraint by summing up two existing constraints
	Constraint& operator +(const Constraint&) const;
	// Creates a new Constraint by subtracting two existing constraints
	Constraint& operator -(const Constraint&) const;
	// Eliminated the variable with the given index from this
	Constraint& mergeConstraints(const Constraint&) const;
	// Divide all Elements of constraint by coeff of highest index variable on lhs
	void normalizeConstraint();

	//---------------------- Methods for lhs elements --------------------------------->
	// Returns the index of the first element of a constraint with a coefficient
	// not equal to zero
	int getFirstConstraintElementIndex() const;
	// Returns the index of the last element of a constraint with a coefficient
	// not equal to zero
	int getLastConstraintElementIndex() const;
	// Returns the index of the following constraint element of a given element or -1
	// if it does not exist
	int getNextConstraintElementIndex(unsigned int) const;
	// Returns the index of the previous constraint element of a given element or -1
	// if it does not exist
	int getPreviousConstraintElementIndex(unsigned int) const;
	// Creates a new Constraint element
	void createConstraintElement(unsigned int, const data::QpNum&);
	//
	void createConstraintElements(const Constraint& c);
	//!!!!Assumes the variables to be ordered
	void createConstraintElements(const std::vector<data::IndexedElement>& c);

	// Returns the coefficient of the constraint entry with the given index
	data::QpNum getConstraintElementCoeff(unsigned int) const;
	// Deletes a constraint element
	void deleteConstraintElement(unsigned int);
	// Adds the value to the constraint element with the given index
	void addConstraintElement(unsigned int, const data::QpNum&);
	// Eliminates the constraint element with the given index by setting its value
	// and normalizing the resulting constraint
	void eliminateConstraintElement(unsigned int, const data::QpNum&);
	// Eliminates the constraint with the given index by inserting an constraint
	// with an equality sign and the index variable on the left sight
	void eliminateConstraintElement(unsigned int, const Constraint&);

	// Get all element values as vector of indexed elements <index,value>
	std::vector<data::IndexedElement> getElements() const;
	void getElements(std::vector<data::QpNum>&) const;
	void getElementsSparse(std::vector<data::IndexedElement>&) const;
	void setElements(const std::vector<data::IndexedElement>&);

	std::string toString() const;
	std::string getDoubleString() const;

	const std::string Point2String(const QlpMatCoeff*p) const;

protected:

	Constraint(Qlp&);
	Constraint(Qlp&, const data::QpRhs&, Constraint*);
	~Constraint();

	inline QlpMatCoeff* getPointPtrByIndex(unsigned int index) const {
		if (!(pFirstRowPoint || pLastRowPoint))
			return NULL;
		if(pLastRowPoint->columnIndex < index || pFirstRowPoint->columnIndex>index)
			return NULL;
		QlpMatCoeff* tmp = pFirstRowPoint;
		while (tmp) {
			if (tmp->columnIndex == index)
				return tmp;
			tmp = tmp->pNextRowPoint;
		}
		return NULL;
	}

	//Hidden
	Constraint(const Constraint&);
	Constraint& operator=(const Constraint&);
	//bool operator==(const Constraint&) const;

private:
    
	// Log string for debug output
	static std::string LOG_TAG;
	// Reference to qlp this constraint belongs to
	Qlp& qlp;
	// The Objective function
	data::QpRhs rhs;
	// The number of elements on the left hand side of this inequality
	unsigned int elements;
	// Pointer to the next Constraint
	Constraint* nextConstraint;
	// Pointer to the previous Constraint
	Constraint* previousConstraint;
	// Pointer to first point that corresponds to this restriction
	QlpMatCoeff* pFirstRowPoint;
	// Pointer to first point that corresponds to this restriction
	QlpMatCoeff* pLastRowPoint;
};
}

#endif /* CONSTRAINT_HPP_ */
