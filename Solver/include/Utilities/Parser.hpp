/*
*
* Solver: Parser.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef PARSER_HPP_
#define PARSER_HPP_
#include "Datastructures/Datastructures.hpp"
//#include <boost/tr1/unordered_map>
#include <unordered_map>
//using namespace std::tr1;
//using namespace std;

namespace utils {
class Parser {
public:
	static void createQlp(const std::string&, data::Qlp&);
	static void createQlp(const std::string&, data::QpObjFunc&,std::vector<data::QpVar>&,data::QpSparseMatrix&,std::vector<data::QpRhs>&);
protected:
	// Constructor, takes path to qlp
	Parser(const std::string&);
	// prints all lines that were read from file
	void printLines(std::list<std::string>&);
	// sorts all line to specific string arrays for seperate parts of the Qlp like bounds, quantfiers and so on
	void sortLines();
	// creates all columns for the Qlp
	void createQlpColumns();
	// creates all columns for the Qlp
	void createQlpColumns(std::vector<data::QpVar>&,data::QpObjFunc&);
	// creates all rows of the Qlp
	void createRows();
	// creates all rows of the Qlp
	void createRows(std::vector<data::QpVar>&,data::QpSparseMatrix&,std::vector<data::QpRhs>&);
	// creates all rows of the Qlp
	void createRows(std::vector<data::QpVar>& vars, data::QpSparseMatrix& lhsIE, std::vector<data::QpRhs>& rhs,std::string Q);
	// add the quantifiers
	void addQuantifiers(std::vector<data::QpVar>&);
	// restrict variables to be integral (generals or binaries)
	void addIntegrals(std::vector<data::QpVar>&);
	// add the quantifiers
	void addBounds(std::vector<data::QpVar>&);
	// create the order of the variables
	void createVariableOrder(std::vector<data::QpVar>&);
	// add the objective function to the qlp
	void addObjectiveFunction(std::vector<data::QpVar>&,data::QpObjFunc&);
	// fill Qlp objects
	void createQlp(data::Qlp&);
	// Get the index of a QpVar by its name
	int getVariableIndexByName(const std::string& varName) const;
	// Get QpVar by name from a vector of QpVer objects
	data::QpVar& getVarByName(const std::string& name, std::vector<
			data::QpVar>& vec)const;
	// Quick First check of used scientific Notation
	bool CheckScientificNotation(data::QpNum& coefficient, data::QpNum& exponent, const int Mode);
	bool BoundIsInf(std::string& bound);
        void CheckDigits(std::string& Coefficient,bool IsScientific=false);
	void clear(){
		objectiveFunction.clear();
		constraints.clear();
		u_constraints.clear();
		bounds.clear();
		generals.clear();
		binaries.clear();
		all.clear();
		exists.clear();
		random.clear();
		order.clear();
	}
private:
	static std::string LOG_TAG;
	static const int INIT = 1;
	static const int OBJECTIVE_FUNCTION = 2;
	static const int CONSTRAINTS = 3;
	static const int BOUNDS = 4;
	static const int GENERALS = 5;
	static const int BINARIES = 6;
	static const int QUANTORS = 7;
	static const int ALL = 8;
	static const int EXISTS = 9;
	static const int RANDOM = 10;
	static const int ORDER = 11;
 	static const int UNCERTAINTY = 12;
	static const int ST_RHS = 13;
	static const int ST_COEF = 14;
	static const int ST_BOUND = 15;
	data::QpObjFunc::Objective objective;
	std::list<std::string> qlp;
	std::string objectiveFunction;
	std::list<std::string> constraints;
	std::list<std::string> u_constraints;
	std::list<std::string> bounds;
	std::string generals;
	std::string binaries;
	std::string all;
	std::string exists;
	std::string random;
	std::string order;

	// Single Qlp parts
	data::QpObjFunc _obj;
	std::vector<data::QpVar> _vars;
	std::vector<data::QpRhs> _rhs;
	std::vector<std::vector<data::IndexedElement> > _lhsIE;

	std::unordered_map<std::string, unsigned int> nameToVariableIndex;

};
}

#endif /* PARSER_HPP_ */
