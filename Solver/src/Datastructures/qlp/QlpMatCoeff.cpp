/*
*
* Solver: QlpMatCoeff.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Datastructures/qlp/QlpMatCoeff.hpp"
#include "Datastructures/qlp/Qlp.hpp"

namespace data {

void QlpMatCoeff::deleteFromRow() {
	if (pNextRowPoint != NULL) {
		pNextRowPoint->pPreviousRowPoint = this->pPreviousRowPoint;
	} else {
		pRow->pLastRowPoint=this->pPreviousRowPoint;
	}
	if (pPreviousRowPoint != NULL) {
		pPreviousRowPoint->pNextRowPoint = this->pNextRowPoint;
	} else {
		pRow->pFirstRowPoint=this->pNextRowPoint;
	}
}

void QlpMatCoeff::deleteFromCol() {
	// check if it is an inner QlpMatCoeff or the last Col QlpMatCoeff
	if (pNextColPoint != NULL) {
		pNextColPoint->pPreviousColPoint = this->pPreviousColPoint;
	} else {
		if(pPreviousColPoint!=NULL){
			pRow->qlp.lastColumnsRowPointer[columnIndex] = pPreviousColPoint;
		}else{
			pRow->qlp.lastColumnsRowPointer[columnIndex] = NULL;
		}
	}

	// check if it is an inner QlpMatCoeff or the first Col QlpMatCoeff
	if (pPreviousColPoint != NULL) {
		pPreviousColPoint->pNextColPoint = this->pNextColPoint;
	} else {
		if(pNextColPoint!=NULL){
			pRow->qlp.firstColumnsRowPointer[columnIndex] = pNextColPoint;
		}else{
			pRow->qlp.firstColumnsRowPointer[columnIndex] = NULL;
		}
	}
}

// deltes this point from the hole Qlp
void QlpMatCoeff::deleteFromRowAndCol() {
	deleteFromRow();
	deleteFromCol();
	delete this;
}
}
