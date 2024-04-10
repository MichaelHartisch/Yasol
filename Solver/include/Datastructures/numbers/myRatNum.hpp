/*
*
* Solver: myRatNum.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef MYRATNUM_HPP_
#define MYRATNUM_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <list>
#include <vector>
#include <sstream>

#ifdef WINDOWS
//#include <stdint.h>
#include <unistd.h>
#endif
class MyRatNum {
public:

	MyRatNum(int n = 0, int d = 1) :
		value(n, d) {
	}

	/*MyRatNum(bool minInf) :
		value(minInf) {
	}*/

	MyRatNum(const MyRatNum& frac) :
		value(frac.value) {
	}

	/*MyRatNum(const double& val) :
		value(val) {
	}

	MyRatNum(const std::string& s) :
		value(s) {
	}*/

	/*MyRatNum& operator =(const double& val) {
		this->value = val;
		return *this;
	}*/

	MyRatNum& operator =(const MyRatNum& frac) {
		this->value = frac.value;
		return *this;
	}

	MyRatNum& operator =(const int x) {
		this->value.first = x;
		this->value.second = 1;
		return *this;
	}

	/*MyRatNum& operator =(const std::string& s) {
		this->value = s;
		return *this;
	}*/

	MyRatNum operator+(const MyRatNum& frac) const {
		return MyRatNum(*this) += frac;
	}

	MyRatNum operator-(const MyRatNum& frac) const {
		return MyRatNum(*this) -= frac;
	}

	MyRatNum operator*(const MyRatNum& frac) const {
		return MyRatNum(*this) *= frac;
	}

	MyRatNum operator/(const MyRatNum& frac) const {
		return MyRatNum(*this) /= frac;
	}

	MyRatNum& operator +=(const MyRatNum& frac) {
		this->value.first = this->value.first * frac.value.second + this->value.second * frac.value.first;
		this->value.second = this->value.second * frac.value.second;
		canonicalize();
		return *this;
	}

	MyRatNum& operator -=(const MyRatNum& frac) {
		this->value.first = this->value.first * frac.value.second - this->value.second * frac.value.first;
		this->value.second = this->value.second * frac.value.second;
		canonicalize();
		return *this;
	}

	MyRatNum& operator /=(const MyRatNum& frac) {
		//std::cerr << "Div: " << value.first << " / " << value.second << " : " << frac.value.first << " / " << frac.value.second << std::endl;
		this->value.first *= frac.value.second;
		this->value.second *= frac.value.first;
		//std::cerr << "DivRes: " << value.first << " / " << value.second << std::endl;
		canonicalize();
		//std::cerr << "DivRes-C: " << value.first << " / " << value.second << std::endl;
		return *this;
	}

	MyRatNum& operator *=(const MyRatNum& frac) {
		this->value.first *= frac.value.first;
		this->value.second *= frac.value.second;
		canonicalize();
		return *this;
	}

	bool operator ==(const MyRatNum& f) const {
		return this->value.first ==(f.value.first) && this->value.second ==(f.value.second);
	}

	bool operator !=(const MyRatNum& f) const {
		return this->value.first !=(f.value.first) || this->value.second !=(f.value.second);
	}

	bool operator <(const MyRatNum& f) const {
		MyRatNum temp;
		temp = *this - f;
		if (temp.getNominator() < 0) return true;
		else return false;
		//return this->value.operator <(f.value);
	}

	bool operator <=(const MyRatNum& f) const {
		MyRatNum temp;
		temp = *this - f;
		if (temp.getNominator() <= 0) return true;
		else return false;
		//return this->value.operator <=(f.value);
	}

	bool operator>(const MyRatNum& f) const {
		MyRatNum temp;
		temp = *this - f;
		if (temp.getNominator() > 0) return true;
		else return false;
		//return this->value.operator >(f.value);
	}

	bool operator >=(const MyRatNum& f) const {
		MyRatNum temp;
		temp = *this - f;
		if (temp.getNominator() >= 0) return true;
		else return false;
		//return this->value.operator >=(f.value);
	}

	/*bool isZero() const {
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

	//operator double() { return this->value.operator double(); }

	std::string toString(unsigned int acc = TO_STRING_ACCURACY) const {
		return this->value.toString(acc);
	}*/

	//------------------------------- CLASS METHODS (some vector operations)-------------------------------->
	template<class T>
	static void stringVecToMyRatNumVec(std::string line, std::vector<T>& vec) {
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
	static bool equalMyRatNumVectors(const std::vector<T>& vec1,
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
	static bool equalMyRatNumVectorNonZeros(const std::vector<T>& vec1,
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
		std::string s("[" + std::string(vec.size()));
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
		std::string s("[" + std::string(vec.size()));
		s += "]<";
		bool b = false;
		for (unsigned int i = 0; i < vec.size(); i++) {
			if (!vec[i].isZero()) {
				if (b)
					s += ", ";
				s += "(";
			    std::stringstream temp;
			    temp<<i;
				s += temp.str();
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


	/*template<class T>
	static void  removeSmallValuesFromVec(std::vector<T>& vec, double abs = DOUBLE_EPSILON) {

#ifdef EXACT_ARITHMETIC
		return;
#endif

		for(unsigned int i = 0;  i < vec.size(); i++){
			if(!vec[i].isZero()&&std::fabs(vec[i].asDouble())<abs){
						vec[i].setZero();
			}
		}
	}*/



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

	friend std::ostream& operator<<(std::ostream &cout, const MyRatNum& obj) {
		std::cout << obj.value.first << " / " << obj.value.second << "\n";
		return cout;
	}

	int64_t ggt(int64_t a, int64_t b) {
		if (a == 0) return b;
		else if (b==0) return a;
		if (a < 0) a = -a;
		if (b < 0) b = -b;
		if (a < b) {
			int64_t t = a; a = b; b = t;
		}
		return ggt(b,(a % b));
	}

	void canonicalize() {
		int64_t g=ggt(value.first, value.second);
		//std::cerr << "value vor ggt: " << value.first << " / " << value.second << std::endl;
		//std::cerr << "ggt=" << g << std::endl;
		value.first /= g;
		value.second /= g;
		//std::cerr << "value nach ggt: " << value.first << " / " << value.second << std::endl;
		//std::cerr << "ggt=" << g << std::endl;
		if (value.second < 0) {
			value.second *= (-1);
			value.first *= (-1);
		}
		//std::cerr << "value am Ende: " << value.first << " / " << value.second << std::endl;

	}

	int64_t getNominator() { return value.first; }
	int64_t getDenominator() { return value.second; }

private:
	std::pair<int64_t, int64_t> value;
};

typedef MyRatNum mpq_class;
typedef int mpz_class;



#endif /* MYRATNUM_HPP_ */
