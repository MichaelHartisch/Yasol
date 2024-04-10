/*
*
* Solver: QlpStageRelaxer.cpp -- Copyright (c) 2010-2017 Jan Wolf
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


#include "Utilities/QlpStageRelaxer.hpp"


namespace utils {
std::string QlpStageRelaxer::LOG_TAG = "QlpStageSolver";
bool QlpStageRelaxer::OLD = true;

/* H�lt QlpSplitter
 * extrahier EXIT und UNIV Matrix
 * L�d sparse LP in QpExtSolCplex --> Univ Part schon nach hinten verschoben und eingesetzt (orig merken)
 * Bekommt Tiefe und pointer auf exitVec und univVec,
 * schaut dann welcher Teil ausgewertet werden muss und �ndert recht Seite
 */

QlpStageRelaxer::QlpStageRelaxer(const data::Qlp& q) : qlp(q),origVarLBs(),origVarUBs(),relaxers(),currVarLBs(),currVarUBs(),externSolvers(){
	this->initialize();
}

QlpStageRelaxer::~QlpStageRelaxer(){
	for(unsigned int i = 0; i < relaxers.size();i++){
		if(relaxers[i])	delete relaxers[i];
	}
	for(unsigned int i = 0; i < externSolvers.size();i++){
			if(externSolvers[i])	delete externSolvers[i];
	}
}

data::QpNum QlpStageRelaxer::solveStage(unsigned int stage){
	data::QpNum val;
	//TODO
	return val;
}

bool QlpStageRelaxer::solveStage(unsigned int stage, data::QpNum& retVal){
	if(OLD){
		retVal = this->relaxers[stage]->getUpperBound();
	}else{
		retVal = this->solveStage(stage);
	}
	return !retVal.isMaxInf();
}

void QlpStageRelaxer::setVariableFixation(unsigned int stage, unsigned int i, const data::QpNum& v){
	if(OLD){
	this->relaxers[stage]->extSol->setVarLB(i,v);
	this->relaxers[stage]->extSol->setVarUB(i,v);
	this->relaxers[stage]->ubComp=false;
	}
	this->currVarLBs[i]=v;
	this->currVarUBs[i]=v;

}

void QlpStageRelaxer::detachFixation(unsigned int stage, unsigned int i){
	if(OLD){
	this->relaxers[stage]->extSol->setVarLB(i,origVarLBs[i]);
	this->relaxers[stage]->extSol->setVarUB(i,origVarUBs[i]);
	this->relaxers[stage]->ubComp=false;
	}
	this->currVarLBs[i]=origVarLBs[i];
	this->currVarUBs[i]=origVarUBs[i];
}

void QlpStageRelaxer::initialize(){
	utils::QlpSplitter qlpSplitter(qlp);
	qlpSplitter.initSplitter(algorithm::Algorithm::WORST_CASE);
	this->origVarLBs = std::vector<data::QpNum>(this->qlp.getVariableCount(),0.0);
	this->origVarUBs = std::vector<data::QpNum>(this->qlp.getVariableCount(),0.0);
	this->currVarLBs = std::vector<data::QpNum>(this->qlp.getVariableCount(),0.0);
	this->currVarUBs = std::vector<data::QpNum>(this->qlp.getVariableCount(),0.0);
	std::vector<data::QpVar*> vVec = qlp.getVariableVector();

	for(unsigned int i = 0; i < origVarLBs.size();i++){
			this->origVarLBs[i]=vVec[i]->getLowerBound();
			this->origVarUBs[i]=vVec[i]->getUpperBound();
	}

	if(OLD){
	for(unsigned int i = 0, size = qlp.getStageCount(); i < size;i++){
		if(i>0){
			for(unsigned int j = 0; j < qlpSplitter.existVarIndAtDepth[i].first;j++){
				qlp.setObjectiveFunctionElement(j,0);
				vVec[j]->setQuantifier(data::QpVar::exists);
			}
		}
		this->relaxers.push_back(new utils::QlpRelaxer(qlp,true));
	}}else{

	}
}

}



