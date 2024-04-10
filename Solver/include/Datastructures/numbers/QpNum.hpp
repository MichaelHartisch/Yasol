/*
*
* Solver: QpNum.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QPNUM_HPP_
#define QPNUM_HPP_
#include "Utilities/ToolBox.hpp"
#include "Datastructures/numbers/QpDouble.hpp"
#include "Datastructures/numbers/QpRational.hpp"
namespace data {
class QpNum {
public:

	explicit QpNum(int n = 0, int d = 1) :
		value(n, d) {
	}

	explicit QpNum(bool minInf) :
		value(minInf) {
	}

	QpNum(const QpNum& frac) :
		value(frac.value) {
	}

	QpNum(const double& val) :
		value(val) {
	}

	/*QpNum(const mpq_class& d) :
		value(d) {
	}*/

	QpNum(const std::string& s) :
		value(s) {
	}

	QpNum& operator =(const double& val) {
		this->value = val;
		return *this;
	}

	/*QpNum& operator =(const mpq_class& d) {
		this->value = d;
		return *this;
	}*/

	QpNum& operator =(const QpNum& frac) {
		this->value = frac.value;
		return *this;
	}

	QpNum& operator =(const std::string& s) {
		this->value = s;
		return *this;
	}

	QpNum operator+(const QpNum& frac) const {
		return QpNum(*this) += frac;
	}

	QpNum operator-(const QpNum& frac) const {
		return QpNum(*this) -= frac;
	}

	QpNum operator*(const QpNum& frac) const {
		return QpNum(*this) *= frac;
	}

	QpNum operator/(const QpNum& frac) const {
		return QpNum(*this) /= frac;
	}

	QpNum& operator +=(const data::QpNum& frac) {
		this->value += frac.value;
		return *this;
	}

	QpNum& operator -=(const data::QpNum& frac) {
		this->value -= frac.value;
		return *this;
	}

	QpNum& operator /=(const data::QpNum& frac) {
		this->value /= frac.value;
		return *this;
	}

	QpNum& operator *=(const data::QpNum& frac) {
		this->value *= frac.value;
		return *this;
	}

	bool operator ==(const QpNum& f) const {
		return this->value.operator ==(f.value);
	}

	bool operator !=(const QpNum& f) const {
		return this->value.operator !=(f.value);
	}

	bool operator <(const QpNum& f) const {
		return this->value.operator <(f.value);
	}

	bool operator <=(const QpNum& f) const {
		return this->value.operator <=(f.value);
	}

	bool operator>(const QpNum& f) const {
		return this->value.operator >(f.value);
	}

	bool operator >=(const QpNum& f) const {
		return this->value.operator >=(f.value);
	}

	bool isZero() const {
		return this->value.isZero();
	}

	void setZero() {
		this->value.setZero();
	}

	bool isMaxInf() const {
		return this->value.isMaxInf();
	}

	bool isMinInf() const {
		return this->value.isMinInf();
	}

	void setMaxInf() {
		this->value.setMaxInf();
	}

	void setMinInf() {
		this->value.setMinInf();
	}

	double asDouble() const {
		return this->value.asDouble();
	}

	void cutDecimals(unsigned int prec){
		this->value.cutDecimals(prec);
	}

	/*mpq_class asMpqClass() const {
		return this->value.asMpqClass();
	}*/

	//operator double() { return this->value.operator double(); }

	std::string toString(unsigned int acc = TO_STRING_ACCURACY) const {
		return this->value.toString(acc);
	}

	//------------------------------- CLASS METHODS (some vector operations)-------------------------------->
	template<class T>
	static void stringVecToQpNumVec(std::string line, std::vector<T>& vec) {
		vec.clear();
		size_t pos = line.find('[');
		line = line.substr(pos + 1);
		while ((pos = line.find(',')) != std::string::npos) {
			vec.push_back(line.substr(0, pos));
			line = line.substr(pos + 1);
		}
		pos = line.find(']');
		vec.push_back(line.substr(0, pos));
	}

	template<class T>
	static bool equalQpNumVectors(const std::vector<T>& vec1,
			const std::vector<T>& vec2) {
		if (vec1.size() != vec2.size())
			return false;
		for (unsigned int i = 0, size = vec1.size(); i < size; i++) {
			if (vec1[i] != vec2[i]) {
				return false;
			}
		}
		return true;
	}

	template<class T>
	static bool equalQpNumVectorNonZeros(const std::vector<T>& vec1,
			const std::vector<T>& vec2) {
		if (vec1.size() != vec2.size())
			return false;
		unsigned int vec1NZ = 0, vec2NZ = 0;
		for (unsigned int i = 0, size = vec1.size(); i < size; i++) {
			if (vec1[i].isZero())vec1NZ++;
			if (vec2[i].isZero())vec2NZ++;
		}
		return (vec1NZ==vec2NZ);
	}

	template<class T>
	static std::string vecToString(const std::vector<T>& vec) {
		std::string s("[" + utils::ToolBox::convertToString(vec.size()));
		s += "]<";
		for (unsigned int i = 0; i < vec.size(); i++) {
			s += vec[i].toString();
			if (i + 1 < vec.size())
				s += ", ";
		}
		s += ">";
		return s;
	}

	template<class T>
	static std::string vecToStringSparse(const std::vector<T>& vec) {
		std::string s("[" + utils::ToolBox::convertToString(vec.size()));
		s += "]<";
		bool b = false;
		for (unsigned int i = 0; i < vec.size(); i++) {
			if (!vec[i].isZero()) {
				if (b)
					s += ", ";
				s += "(";
				s += utils::ToolBox::convertToString(i);
				s += ",";
				s += vec[i].toString();
				s += ")";
				b = true;
			}
		}
		s += ">";
		return s;
	}

	template<class T>
	static std::string matrixToString(const std::vector<std::vector<T> >& vec) {
		std::string s("\n-->\n");
		for (unsigned int i = 0; i < vec.size(); i++) {
			s += vecToString(vec[i]);
			s += "\n";
		}
		s += "<--";
		return s;
	}


	template<class T>
	static void  removeSmallValuesFromVec(std::vector<T>& vec, double abs = DOUBLE_EPSILON) {

#ifdef EXACT_ARITHMETIC
		return;
#endif

		for(unsigned int i = 0;  i < vec.size(); i++){
			if(!vec[i].isZero()&&std::fabs(vec[i].asDouble())<abs){
						vec[i].setZero();
			}
		}
	}



	template<class T>
	static std::vector<double> vecToDoubleVec(const std::vector<T>& vec) {
		std::vector<double> tmpVec;
		for (unsigned int i = 0; i < vec.size(); i++)
			tmpVec.push_back(vec[i].asDouble());
		return tmpVec;
	}

	template<class T>
	static T sumQpRationalVec(const std::vector<T>& vec) {
		T sum(0);
		for (unsigned int i = 0; i < vec.size(); i++) {
			sum += vec[i];
		}
		return sum;
	}

	template<class T>
	static T prodQpRationalVec(const std::vector<T>& vec) {
		T sum(1);
		for (unsigned int i = 0; i < vec.size(); i++) {
			sum *= vec[i];
		}
		return sum;
	}

	friend std::ostream& operator<<(std::ostream &cout, const QpNum& obj) {
		std::cout << obj.toString() << "\n";
		return cout;
	}

private:
#ifndef EXACT_ARITHMETIC
	data::QpDouble value;
#endif
#ifdef EXACT_ARITHMETIC
	data::QpRational value;
#endif
};

}
#endif
