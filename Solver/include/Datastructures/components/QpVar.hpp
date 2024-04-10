/*
*
* Solver: QpVar.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QUANTIFIEDVARIABLE_HPP_
#define QUANTIFIEDVARIABLE_HPP_
#include "Datastructures/numbers/Numbers.hpp"
#include "cstring"

namespace data {
class QpVar {
public:

	//variable quantifier
	typedef enum {
		exists, all, random
	} Quantifier;

	//variable number system
	typedef enum {
		real, generals, binaries
	} NumberSystem;

	QpVar();

	//constructor for existential or universal variable
	QpVar(const std::string& name, unsigned int index, const QpNum& lb,
			const QpNum& ub, NumberSystem numberSystem = real,
			Quantifier quantifier = exists);

	//constructor for existential or universal variables
	QpVar(const std::string& name, unsigned int index, const std::vector<data::QpNum>& vRange,
				const std::vector<data::QpRational>& vDistr, NumberSystem ns);
	//destructor
	~QpVar();

	//copy oonstructor
	QpVar(const QpVar&);
	//assignment operator
	QpVar& operator=(const QpVar&);

	//set and get variable name
	void setName(std::string& n);
	const std::string& getName() const;

	//set and get the index of this variable
	void setIndex(unsigned int i);
	unsigned int getIndex() const;

	//set and get the lower/upper bound of this variable
	void setLowerBound(const QpNum& lb);
	const QpNum& getLowerBound() const;
	void setUpperBound(const QpNum& ub);
	const QpNum& getUpperBound() const;

	//set and get upper and lower bounds of this variable
	void setBounds(const QpNum&, const QpNum&);
	std::pair<QpNum,QpNum> getBounds() const;

	// set and get the quantifier of this variable
	void setQuantifier(Quantifier q);
	Quantifier getQuantifier() const;

	//set and get the number system type of this variable
	void setNumberType(NumberSystem n);
	NumberSystem getNumberSystem() const;

	//comparison operators
	bool operator ==(const QpVar&) const;
	bool operator !=(const QpVar&) const;

	// return if this is a free variable -inf <= x <= +inf
	bool isFree() const;
	// return if variable os bounded from above x <= ub
	bool isBoundedAbove() const;
	// return if variable is bounded from below lb <= x
	bool isBoundedBelow() const;
	// sets this variable to be free
	void setFree();
	// sets this variable to be unbounded from above
	void setUnboundedAbove();
	// sets this variable to be unbounded from below
	void setUnboundedBelow();

	//set the variable range of a random variable
	void setVariableRange(const std::vector<data::QpNum>& range,
			const std::vector<data::QpRational>& distr);
	// get the variable range of the variable, except for random variable the size is two
	const std::vector<data::QpNum>& getVariableRange() const;
	// get the variable range distribution of the variable, except for random variable the size is zero
	const std::vector<data::QpRational>& getVariableDistribution() const;

	// Convert a quantified variable object to a readable string
	std::string toString() const;

	// Convert a vector of quantified variable objects to a readable string
	static std::string vecToString(const std::vector<data::QpVar>&);
	// Converts the given quantifier to a string.
	static std::string quantifierString(const QpVar::Quantifier&);
	// Converts the given string to a quantifier.
	static QpVar::Quantifier stringQuantifier(const std::string&);
	// Converts the given quantifier to a string.
	static std::string numberString(const QpVar::NumberSystem&);
	// Converts the given string to a quantifier.
	static QpVar::NumberSystem stringNumber(const std::string&);

	// Get QpVar by name from a vector of QpVer objects, the first occurence will be passed back
	static QpVar& getQpVarByName(const std::string& name, std::vector<
			data::QpVar>& vec) {
		for (unsigned int i = 0, size = vec.size(); i < size; i++) {
			if (!std::strcmp(vec[i].getName().c_str(), name.c_str()))
				return vec[i];
		}
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable", "Name: "
				+ name);
		utils::Logger::globalLog(utils::LOG_INFO, "writeTexTable",
				"Variables: " + data::QpVar::vecToString(vec));
		throw utils::DataStructureException(
				"getQpVarByName --> variable does not exist: "+name);
	}
protected:
private:
	std::string name;
	unsigned int index;
	std::vector<data::QpNum> vRange;
	std::vector<data::QpRational> vDistr;
	Quantifier q;
	NumberSystem n;
};
}

#endif /*QUANTIFIEDVARIABLE_HPP_*/
