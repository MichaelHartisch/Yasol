/*
*
* Solver: Numbers.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef NUMBERS_HPP_
#define NUMBERS_HPP_
#include "Datastructures/numbers/QpNum.hpp"
#include "Datastructures/numbers/QpDouble.hpp"
#include "Datastructures/numbers/QpRational.hpp"
namespace data {

template<class T>
struct QpMatrix {
	typedef std::vector<std::vector<T> > Type;
};
struct IndexedElement;
typedef QpMatrix<data::IndexedElement>::Type QpSparseMatrix;
typedef QpMatrix<data::QpRational>::Type QpRationalMatrix;
typedef QpMatrix<data::QpNum>::Type QpNumMatrix;
typedef QpMatrix<data::QpDouble>::Type QpDoubleMatrix;
typedef QpMatrix<data::QpRational>::Type QpRationalMatrix;

struct QpCSMatrix {
	std::vector<QpNum> elements;
	std::vector<int> indices;
	std::vector<int> starts;
	std::vector<int> lengths;
};

struct IndexedElement {
	unsigned int index;
	data::QpNum value;
	IndexedElement(unsigned int i = 0, const data::QpNum& f = 0) :
		index(i), value(f) {
	}
	IndexedElement(const IndexedElement& ie) :
		index(ie.index), value(ie.value) {
	}

	void set(const IndexedElement& ie){
			this->index = ie.index;
			this->value = ie.value;
	}

	void set(unsigned int i, const data::QpNum& v){
		this->index = i;
		this->value = v;
	}

	IndexedElement& operator =(const IndexedElement& ie) {
		this->index = ie.index;
		this->value = ie.value;
		return *this;
	}
	bool operator ==(const IndexedElement& i) const {
		return (this->index == i.index && this->value == i.value);;
	}
	bool operator !=(const IndexedElement& i) const {
		return (this->index != i.index || this->value != i.value);
	}
	bool operator <(const IndexedElement& i) const {
		return (this->index < i.index);
	}

	std::string toString() const{
		return ((value.toString()) += "*x") += utils::ToolBox::convertToString(
					index);
	}
};

static void getSubMatrix(const QpSparseMatrix& source, QpSparseMatrix& target, unsigned int fRow, unsigned int tRow, unsigned int fCol, unsigned int tCol){
	if(fRow>tRow||fCol>tCol)
		throw utils::DataStructureException("getSubMatrix(const QpSparseMatrix& source ... ) --> fRow>tRow||fCol>tCol");
	if(source.size()<=tRow)
		throw utils::DataStructureException("getSubMatrix(const QpSparseMatrix& source ... ) --> source.size()<=tRow");
	target.clear();
	unsigned int newRows = 1+ tRow - fRow;
	for(unsigned int i = 0; i < newRows;i++){
		target.push_back(std::vector<data::IndexedElement>());
	}
	data::IndexedElement ie;
	for(unsigned int i = fRow, k = 0; i <= tRow;i++,k++){
		for(unsigned int j = 0; j < source[i].size();j++){
			if((ie=source[i][j]).index>=fCol && ie.index<=tCol){
				target[k].push_back(data::IndexedElement(ie.index-fCol,ie.value));
			}
		}
	}
	//TODO push back as many rows as diif 1+tRow-fRow
}

static bool operator ==(const std::vector<IndexedElement>& v1,const std::vector<IndexedElement>& v2){
	if(v1.size()!=v2.size())return false;
	for(unsigned int i = 0; i < v1.size();i++){
		if(v1[i]!=v2[i])return false;
	}
	return true;
}

static bool operator !=(const std::vector<IndexedElement>& v1,const std::vector<IndexedElement>& v2){
	return !operator ==(v1,v2);
}

static std::string indexedElementToString(const IndexedElement& i) {
	return ((i.value.toString()) += "*x") += utils::ToolBox::convertToString(
			i.index);
}

static std::string indexedElementVecToString(
		const std::vector<IndexedElement>& vec) {
	std::string s("[");
	for (unsigned int i = 0; i < vec.size(); i++) {
		s += indexedElementToString(vec[i]);
		if (i != vec.size() - 1)
			s += "+";
	}
	s += "]";
	return s;
}

static std::string indexedElementMatrixToString(
		const std::vector<std::vector<IndexedElement> >& m) {
	std::string s("[");
	for (unsigned int i = 0; i < m.size(); i++) {
		s += "\n\t" + indexedElementVecToString(m[i]);
	}
	s += "]";
	return s;
}

static void quickSortIeVec(std::vector<data::IndexedElement>& ieVec, int left, int right) {
      int i = left, j = right;
      data::IndexedElement tmp;
      int pivot = ieVec[(left + right) / 2].index;

      /* partition */
      while (i <= j) {

    	  while(ieVec[i].index < pivot)
    		  i++;
    	  while(ieVec[i].index > pivot)
    		  j--;

    	  if (i <= j) {
                  tmp = ieVec[i];
                  ieVec[i] = ieVec[j];
                  ieVec[j] = tmp;
                  i++;
                  j--;
            }
      };
      /* recursion */
      if (left < j)
            data::quickSortIeVec(ieVec, left, j);
      if (i < right)
            data::quickSortIeVec(ieVec, i, right);
}

static int partition(std::vector<data::IndexedElement> & vec, int left, int right, int pivotIndex)
{
    int pivot = vec[pivotIndex].index;
    do
    {
        while (vec[left].index < pivot) left++;
        while (vec[right].index > pivot) right--;
        if (left < right && vec[left].index != vec[right].index)
        {
        	data::IndexedElement tmp = vec[left];
        	vec[left]=vec[right];
            vec[right]=tmp;
        }
        else
        {
            return right;
        }
    }
    while (left < right);
    return right;
}

static void quicksort(std::vector<data::IndexedElement> & vec, int left, int right)
{
    if (left < right)
    {
        int pivot = (left + right) / 2; // middle
        int pivotNew = partition(vec, left, right, pivot);
        quicksort(vec, left, pivotNew - 1);
        quicksort(vec, pivotNew + 1, right);
    }
}
}

#endif /* NUMBERS_HPP_ */
