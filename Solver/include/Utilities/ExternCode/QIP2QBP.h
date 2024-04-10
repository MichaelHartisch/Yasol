/*
*
* Solver: QIP2QBP.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef QIP2QBP_HPP_
#define QIP2QBP_HPP_

#include "QLPReader.hpp"

#include "Utilities/Parser.hpp"
#include "Settings/Settings.hpp"

#include "Datastructures/Datastructures.hpp"
#include "Algorithms/Algorithms.hpp"



using namespace std;

class QIP2QBP
{
public:
	QIP2QBP(void);
	~QIP2QBP(void);
	data::Qlp qlp;
	data::Qlp qbp;
	std::vector<const data::QpVar *> existsVariablesQlp;
	std::vector<int> ExistBlockLength;
	std::vector<int> ExistBlockStart;
	std::vector<int> existIndexByIndex;
	std::vector<int> firstIndexBByIndex;
	std::vector<int> indexDByIndex;
	int NumberExistBlock;
	vector<int> numsBinaries;
	bool onlySum; // =1, wenn alle q_t im Modell nur als Bin�rdarstellung vorkommen sollen.
	int totalNumBinariesQbp;
	int totalNumDsQlp;

	void readQIP(string path);
	void prepareBs();
	int countBs(int);
	void makeQBP();
	void calcExistBlocks();
	void initExistIndexByIndex();
	void initFirstIndexBByIndex();
	void initIndexDByIndex();


	// gibt den Index der i-ten Variable des QLPs in numBinaries an (d.h. gibt aus, um die wievielte Existvariable es sich handelt)
	int getIndexExist(int i);
	// gibt den Index der ersten Binärvariable (b_1i) zur i-ten Existvariable im QBP an (d.h. nicht i-te Variable im QLP, sondern i-te Existvariable)
	int getFirstIndexB(int i);
	int getFirstIndexB_init(int i);
	// gibt den Index der i-ten Variable des QLPs im QBP an (nur für Forall-Variablen)
	int getIndexD(int i);
	// gibt bei onlySum=false für Existvariablen den Index der i-ten Variable des QLPs im QBP aus (Existvariablen aus QLP werden in den letzten Existblock des QBPs geschoben)
	int getIndexP(int i);
	void print();
};

#endif /* QIP2QBP_HPP_ */
