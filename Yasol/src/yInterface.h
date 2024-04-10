/*
*
* Yasol: yInterface.h -- Copyright (c) 2012-2017 Ulf Lorenz
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

#ifndef YINTERFACE_H_
#define YINTERFACE_H_

#include <queue>
#include <vector>
#include <set>

#include "QBPSolver.h"
#include "MipSolvers.h"
#define EXIST 0
#define UNIV  1

#define NORMAL_MODE      0
#define RESTRICTION_MODE 1
#define RELAXATION_MODE  2

const int maxUnivVars = 0;

struct node {
        int index;
        int i;
        bool onStack;
        int lowlink;
};

class IntInfo {
public:
   int index;
   int pt2leader;
   data::QpNum org_lb;
   data::QpNum org_ub;
   int org_ind;
   int bitcnt;
   std::string name;
   double tmp_x;
   int number;
};

class Component {
	std::vector< std::vector<int> > &varsOfComponents;
	int n;
	int m;
	std::vector< std::vector<int> > &conVec;
	std::vector<const data::Constraint *> &rows;
	int8_t *assigns;
	std::vector<bool> rowInd;
	std::vector<int> col;
	std::vector<int> d;
	std::vector<int> l;
	std::vector<int> f;
	std::vector<int> Pi;
	std::vector<int> Q;
	std::vector<int> stack;
	std::vector<data::IndexedElement> lhs;
	int TIME;
	const int WHITE = 0;
	const int GREY  = 1;
	const int BLACK = 2;
public:
	Component(int nn, int rs,std::vector< std::vector<int> > &vOfCs, std::vector< std::vector<int> > &cV, std::vector<const data::Constraint *> &rws, int8_t *as) :
	          n(nn), m(rs), varsOfComponents(vOfCs), conVec(cV), rows(rws), assigns(as) {
		rowInd.resize(nn);
		for (int i = 0; i < nn;i++) rowInd[i] = false;
		col.resize(nn);
		for (int i = 0; i < nn;i++) col[i] = WHITE;
		d.resize(nn);
		for (int i = 0; i < nn;i++) d[i] = -1;
		l.resize(nn);
		for (int i = 0; i < nn;i++) l[i] = -1;
		f.resize(nn);
		for (int i = 0; i < nn;i++) f[i] = -1;
		Pi.resize(nn);
		for (int i = 0; i < nn;i++) Pi[i] = -1;
		TIME = 0;
	}
	bool isBiConnected(int c) {  // under construction
		int v = varsOfComponents[c][0];
		for (int i = 0; i < n;i++) col[i] = WHITE;
		for (int i = 0; i < n;i++) Pi[i] = -1;
		TIME = 0;
		cerr << "Component " << c << " deals with the following Variables: ";
		for (int i = 0; i < varsOfComponents[c].size();i++) {
			cerr << varsOfComponents[c][i] << " ";
		}
		cerr << endl;
		MDFS_visit_iterativ(v);
		for (int i = 0; i < varsOfComponents[c].size();i++) {
			if (i == 0) {
				int rootCnt=0;
				for (int j = 0; j < varsOfComponents[c].size();j++)
					if (Pi[varsOfComponents[c][i]] == v) rootCnt++;
				if (rootCnt > 1) cerr << "y" << varsOfComponents[c][i] << " ist AP" << endl;
			} else {
				if (l[varsOfComponents[c][i]] < d[varsOfComponents[c][i]]) ;
				else cerr << "x" << varsOfComponents[c][i] << " ist AP" << endl;
			}
		}
		return true;
	}
	void MDFS_visit_iterativ(int u) {
		stack.clear();
		stack.push_back(u);
		TIME++;
	    d[u] = TIME;
	    l[u] = TIME;
		TIME++;
		cerr << "start with x" << u << " at TIME=" << TIME << endl;
		while (stack.size()>0) {
			int var = stack.back();
			//stack.pop_back();
			if (col[var] == WHITE) {
				col[var] = GREY;
				for (int j = 0; j < conVec[var].size();j++) {
					lhs.clear();
					if (rowInd[conVec[var][j]] == true) continue;
					lhs = rows[conVec[var][j]]->getElements();//conVec[j]->getElements();
					rowInd[conVec[var][j]] = true;
					for (int jj = 0; jj < lhs.size();jj++) {
						if (assigns[lhs[jj].index] != 2) continue;
						if (col[lhs[jj].index] == WHITE) {
							if (jj == 0) Pi[lhs[jj].index] = var;
							else Pi[lhs[jj].index] = lhs[jj-1].index;
							d[lhs[jj].index] = TIME;
							stack.push_back(lhs[jj].index);
							TIME++;
							if (lhs.size()>2) l[lhs[jj].index] = l[var];
							else l[lhs[jj].index] = d[lhs[jj].index];
						} else if (col[lhs[jj].index] == GREY && Pi[var] != lhs[jj].index) {
							l[var] = (l[var]<d[lhs[jj].index] ? l[var] : d[lhs[jj].index]);
						}
					}
				}
			} else if (col[var] == GREY) {
				col[var] = BLACK;
				TIME++;
				f[var] = TIME;
				if (l[Pi[var]] > l[var]) l[Pi[var]] = l[var];
				stack.pop_back();
			} else if (col[var] == BLACK) {
				if (l[Pi[var]] > l[var]) l[Pi[var]] = l[var];
				stack.pop_back();
			}
		}
	}

	void MDFS_visit2(int u) {
		col[u] = GREY;
		d[u] = TIME;
		l[u] = d[u];
		TIME = TIME + 1;
		for (int j = 0; j < conVec[u].size();j++) {
			lhs.clear();
			lhs = rows[conVec[u][j]]->getElements();
			for (int jj = 0; jj < lhs.size();jj++) {
				if (assigns[lhs[jj].index] != 2) continue;
				if (col[lhs[jj].index] == WHITE) {
					Pi[lhs[jj].index] = u;
					MDFS_visit2(lhs[jj].index);
					l[u] = (l[u] < l[lhs[jj].index] ? l[u] : l[lhs[jj].index]);
				}
				if (col[lhs[jj].index] == GREY && Pi[u] != lhs[jj].index) {
					l[u] = (l[u]<d[lhs[jj].index] ? l[u] : d[lhs[jj].index]);
				}
			}
		}
		col[u] = BLACK;
		f[u] = TIME;
		TIME = TIME + 1;
	}
	void NDFS_visit(int u, int c) {
		col[u] = GREY;
		for (int j = 0; j < conVec[u].size();j++) {
			lhs.clear();
			lhs = rows[conVec[u][j]]->getElements();
			for (int jj = 0; jj < lhs.size();jj++) {
				if (assigns[lhs[jj].index] != 2) continue;
				if (col[lhs[jj].index] == WHITE) {
					Q.push_back(lhs[jj].index);
					NDFS_visit(lhs[jj].index,c);
					if (l[lhs[jj].index] >= d[u]) { // u ist Artikulationspunkt
						cerr << "x" << u << " ist ein AP der Komponente " << c << endl;
						//solange Knoten von Q entfernen und ausgeben, bis v ausgegeben.
						//Knoten u auch mit ausgeben.
					}
				}
			}
		}
		col[u] = BLACK;
	}
};

class yInterface {
public:
	struct set_compare {
	    bool operator() (const std::pair<int,int>& l, const std::pair<int,int>& r) const{
	        return l.first < r.first;
	    }
	};
private:
	int yParamMaxHashSize;

	int processNo;
	int info_level;
  bool has_basis;
   	//data::Qlp qlp,dep,qlpRelax,qlptmp;
	utils::QlpStageSolver *QlpStSolve;
    coef_t result;
    time_t timeout;
    std::vector<data::QpNum> reducedCost;
    extSol::QpExternSolver::QpExtSolBase basis;
    const coef_t n_infinity = -1;
    const coef_t p_infinity = 1;
    std::vector<bool> alive;
    std::vector<int> sortcols;
    std::vector< std::pair< std::vector< std::pair<int,double> >,int > > SOSconstraints;
    std::set<std::pair<int,int>, set_compare > SOSvars;
    std::set<int> SOStabus;
    double LPoffset;
    std::vector<node> stack_S;
    std::vector< std::pair< std::vector<data::IndexedElement>,double > > RelaxationBuffer;
    bool gzipExists;
    std::string qdimacs2QipPath;

    void setQdimacs2QIPPath(std::string &str) {
       qdimacs2QipPath = str;
    }

    void yReadIniFile();

public:
	void yReadInput(int) ;
   	data::Qlp qlp,dep,qlpRelax, qlptmp;
	QBPSolver *qbp;
#ifndef NO_CGL
        CBCSolver *cbc;
#endif
        bool objInverted;
	const int GMI = 1;
	const int Cover = 2;
	const int LaP = 4;
	const int UserCut = 8;
	const int CGL_LIB = 16 | 32 | 64 | 128 | 256 | 512 | 1024 | 2048 | 4096;
    void ySetProcessNo(int pno) { processNo = pno; qbp->ySetProcessNo(pno); }
    int64_t getNumberOfDecisions() { return qbp->getNumberOfDecisions(); }
    int64_t getNumberOfPropagationSteps() { return qbp->getNumberOfPropagationSteps(); }
    int64_t getNumberOfQlpCalls() { return qbp->getNumberOfQlpCalls(); }
    int64_t getNumberOfLearntConstraints() { return qbp->getNumberOfLearntConstraints(); }
    int getmRows() { return qbp->mRows(); }
    int getnVars() { return qbp->nVars(); }
    time_t getTimeout() { return timeout; }
    coef_t getNinf() { return qbp->getNegativeInfinity(); }

	int confl_var;
	CRef confl, confl_partner;

	std::vector< IntInfo > integers;
	std::vector< int > integers_lim;
	double offset_shift=0.0;

	void writeRealName(int i) {
	  int leader = integers[i].pt2leader;
	  int leader_index = integers[leader].index;
	  std::string &name = integers[leader].name;
	  std::cout << "variable x" << i << ": name is " << name << std::endl; 
	}

	void updateNode(int nodeID) {
	  qbp->MCTS.updateFatherScore(nodeID, true);
	  std::cerr << "new father value = " << qbp->MCTS.nodes[ qbp->MCTS.nodes[nodeID].fatherID ].minmax_bnd << std::endl;
	}

	double ipow2(int x) {
		assert(x >= 0);
		double r=1.0;
		if (x <= 0) return 1.0;
		for (int i = 0; i < x;i++) r = r*2.0;
		return r;
	}

	void Reduce(std::string inputfile){

		//Inhalt von Eingabedatei in Datenstrukturen übernehmen
		std::vector<const data::QpRhs *> rhsVec=qlp.getRhsVecConst();
		std::vector<const data::Constraint *> conVec=qlp.getConstraintVecConst();
		data::QpObjFunc zielfkt=qlp.getObjectiveFunction();
		std::vector<data::IndexedElement> zielwerte=zielfkt.getObjectiveElementsSparse();
		int m=conVec.size(), n=getnVars();

		//Existenz- und All-Constraints sowie rhsVecs seperat abspeichern
		std::vector<const data::Constraint *> conVecE;
		std::vector<const data::Constraint *> conVecU;
		std::vector<const data::QpRhs *> rhsVecE;
		std::vector<const data::QpRhs *> rhsVecU;

		std::vector<bool> IsInE(n,false);
		std::vector<bool> IsInA(n,false);

		for(int i=0;i<conVec.size();i++){
			if(rhsVec[i]->getResponsibility()==0){ //Dann Existenzspieler-Constraint
				conVecE.push_back(conVec[i]);
				rhsVecE.push_back(rhsVec[i]);
			}
			if(rhsVec[i]->getResponsibility()==1){ //Dann Allspieler-Constraint
				conVecU.push_back(conVec[i]);
				rhsVecU.push_back(rhsVec[i]);
			}
		}

		//Anzahl Existenzconstraints
		int anzahlECon=0;
		for(int e=0;e<m;e++){
			if(rhsVec[e]->getResponsibility()==0){
				anzahlECon++;
			}
		}
		//Anzahl Allconstraints
		int anzahlUCon=m-anzahlECon;
		assert(anzahlECon==conVecE.size());
		assert(anzahlUCon==conVecU.size());

		//alle Variablen in einem Vektor in der richtigen (angegebenen) Reihenfolge!
		std::vector<data::QpVar *> vecVariablen=qlp.getVariableVector();

		//Existenxvariablen in seperaten Vektor speichern
		std::vector<const data::QpVar *> vecE=qlp.getVariableVectorByQuantifierConst(data::QpVar::exists);

		//Allvariablen in seperaten Vektor speichern
		std::vector<const data::QpVar *> vecU=qlp.getVariableVectorByQuantifierConst(data::QpVar::all);

		cerr << "Start building the reduced QIP." << endl;
		cerr << "INFO:" << endl;
		cerr << "==================================================================" << endl;
		cerr << "Variables with names 'p', 't_$INDEX' and 'CheckAllC$INDEX_InS' are added. Please ensure that your QIPID does not contain variables with such names!" <<endl;
		cerr << "Continuous Variables are ONLY allowed in a FINAL existential stage!" << endl;
		cerr << "Continuous variables are NOT allowed to be present within universal constraints!" <<endl;
		cerr << "The orientation of the objective in the reduced QIP will always be 'MINIMIZE'" << endl;
		cerr << "The first variable block should be an existential variable block." << endl;
		cerr << "The final variable block should be a existential /*universal*/ variable block." << endl;
                cerr << "==================================================================" << endl;

		//Existenz- oder All-Zugehörigkeit des ersten Blocks abfragen, getQuantifier gibt 0=E, 1=A
		int ersterBlock;
		if(vecVariablen[0]->getQuantifier()==1){//1=All; 0=Existenz
			cerr << "WARNING: First variable block is a universal variable block." << endl;
			cerr << "Resulting QIP may contain errors! \n" << endl; 
			ersterBlock=1;
		}
		else{
			ersterBlock=0;
		}

		bool FinalBlockIsUniversal=true;
		//Existenz- oder All-Zugehörigkeit des letzten Blocks abfragen, getQuantifier gibt 0=E, 1=A
		if(vecVariablen[vecVariablen.size()-1]->getQuantifier()==0){//1=All; 0=Existenz
			FinalBlockIsUniversal=false;
			cerr << "WARNING: Final variable block is an existential variable block." << endl;
			cerr << "An artificial final universal block is added.\n" << endl;
		}

		//Anzahl Blöcke insgesamt, ermittelt über Block der letzten Variable in vecVariablen
		int numBlocks=qbp->getBlock(vecVariablen[vecVariablen.size()-1]->getIndex());

		//Anzahl Existenzblöcke
		int numEBlocks=qlp.getStageCount(); //gibt Anzahl E-Blöcke, daraus U-Blöcke ableiten

		//Anzahl Allblöcke
		int numABlocks=numBlocks-numEBlocks;

		//mTilde, nur einmal (kein Vektor)
		double mTilde=0;
		int Switch=1;
		//cerr << zielfkt.getObjective() << endl;
		assert(zielfkt.getObjective()==data::QpObjFunc::Objective::min);
		if (zielfkt.getObjective()==data::QpObjFunc::Objective::max)
			Switch=-1;
		for(int a=0;a<zielwerte.size();a++){
			if(zielwerte[a].value.asDouble()<0)
				mTilde=mTilde+Switch*zielwerte[a].value.asDouble()*(qbp->getLowerBound(zielwerte[a].index)-qbp->getUpperBound(zielwerte[a].index));
			else 	
				mTilde=mTilde+Switch*zielwerte[a].value.asDouble()*(qbp->getUpperBound(zielwerte[a].index)-qbp->getLowerBound(zielwerte[a].index));
			}
		mTilde=mTilde+1;
		mTilde=Switch=mTilde;
		//0<R_LCD<=1 für Integer
		/*
		Weg über euklidischen Algorithmus führt leider bei ungleich vielen Nachkommastellen
		der Koeffizienten zu Ungenauigkeit
		Alternativer Ansatz -> Frage: wie klein kann es werden? 0.6 , 0.876 , => 0.001 kleinstmöglich
		 */
		int rhsR, anzahlNachkommastellen, globalAnzahlNachkommastellen;
		double nenner, ergebnisKehrwert;
		std::vector<double> rLCDVec;

		//L berechnen, bezieht sich auf All-System
		int lK;
		std::vector<int> lKVec;
		for(int c=0;c<anzahlUCon;c++){
			//____________________________________________________________________________
			//Calculate R^LCD
			//____________________________________________________________________________

			globalAnzahlNachkommastellen=0;
			std::vector<data::IndexedElement> lhs=conVecU[c]->getElements();
			std::vector<double> keineInt;
			for(int l=0;l<lhs.size();l++){
				if(lhs[l].value.asDouble()-floor(lhs[l].value.asDouble())!=0){
					keineInt.push_back(abs(lhs[l].value.asDouble()));

				}
			}
			if(rhsVecU[c]->getValue().asDouble()-floor(rhsVecU[c]->getValue().asDouble())!=0){
				keineInt.push_back(abs(rhsVecU[c]->getValue().asDouble()));
			}
			if(keineInt.size()==0){
				rLCDVec.push_back(1);
			}
			else{
				nenner=1.0;
				for(int o=0;o<keineInt.size();o++){
					anzahlNachkommastellen=0;
					while(keineInt[o]-floor(keineInt[o]) !=0){
						keineInt[o]=10*keineInt[o];
						anzahlNachkommastellen ++;
					}
					if(anzahlNachkommastellen>globalAnzahlNachkommastellen){
						globalAnzahlNachkommastellen=anzahlNachkommastellen;
					}
				}
				for(int a=0;a<globalAnzahlNachkommastellen;a++){
					nenner=nenner*10;
				}
				ergebnisKehrwert=1.0/(double)nenner;
				rLCDVec.push_back(ergebnisKehrwert);
			}

			//____________________________________________________________________________
			//Calculate L
			//____________________________________________________________________________

			if (conVecU[c]->getRhs().getRatioSign() == data::QpRhs::smallerThanOrEqual || conVecU[c]->getRhs().getRatioSign() == data::QpRhs::equal )
			{
				lK=0;
				for(int d=0;d<lhs.size();d++){
					IsInA[lhs[d].index]=true;
					if(lhs[d].value.asDouble()<0){
						lK=lK+(lhs[d].value.asDouble()*qbp->getUpperBound(lhs[d].index));
					}
					else{
						lK=lK+(lhs[d].value.asDouble()*qbp->getLowerBound(lhs[d].index));
					}
				}
				//cerr << c << "   " << lK <<endl;
				lKVec.push_back(lK);

			}

			if (conVecU[c]->getRhs().getRatioSign() == data::QpRhs::greaterThanOrEqual || conVecU[c]->getRhs().getRatioSign() == data::QpRhs::equal )
			{
				lK=0;
				for(int d=0;d<lhs.size();d++){
					IsInA[lhs[d].index]=true;
					if(lhs[d].value.asDouble()>0){
						lK=lK-(lhs[d].value.asDouble()*qbp->getUpperBound(lhs[d].index));
					}
					else{
						lK=lK-(lhs[d].value.asDouble()*qbp->getLowerBound(lhs[d].index));
					}
				}
				//cerr << c << "   " << lK <<endl;
				lKVec.push_back(lK);
			}
		}

		//M berechnen, bezieht sich auf Existenz-System
		int mK;
		std::vector<int> mKVec;
		for(int e=0;e<anzahlECon;e++){
			std::vector<data::IndexedElement> lhs=conVecE[e]->getElements();
			if (conVecE[e]->getRhs().getRatioSign() == data::QpRhs::smallerThanOrEqual || conVecE[e]->getRhs().getRatioSign() == data::QpRhs::equal ){
				mK=0;
				for(int f=0;f<lhs.size();f++){
					IsInE[lhs[f].index]=true;
					if(lhs[f].value.asDouble()<0){
						mK=mK+(lhs[f].value.asDouble()*qbp->getLowerBound(lhs[f].index));
					}
					else{
						mK=mK+(lhs[f].value.asDouble()*qbp->getUpperBound(lhs[f].index));
					}
				}
				mK=mK-rhsVecE[e]->getValue().asDouble();
				//cerr << e << " "<<mK << endl;

				if(mK>=0){
					mKVec.push_back(mK);
				}
				else{
					mK=0;
					mKVec.push_back(mK);
				}
			}

			// NO(!) ElseIf
			if (conVecE[e]->getRhs().getRatioSign() == data::QpRhs::greaterThanOrEqual || conVecE[e]->getRhs().getRatioSign() == data::QpRhs::equal){
				mK=0;
				for(int f=0;f<lhs.size();f++){
					IsInE[lhs[f].index]=true;
					if(lhs[f].value.asDouble()>0){
						mK=mK-(lhs[f].value.asDouble()*qbp->getLowerBound(lhs[f].index));
					}
					else{
						mK=mK-(lhs[f].value.asDouble()*qbp->getUpperBound(lhs[f].index));
					}
				}
				mK=mK+rhsVecE[e]->getValue().asDouble();
				//cerr << e << " "<<mK << endl;
				if(mK>0){
					mKVec.push_back(mK);
				}
				else{
					mK=0;
					mKVec.push_back(mK);
				}
			}


			//	cerr << "Error: Constraint not of type Ax<=b: " << endl;
		//		assert(0);
		}


		/*mit euklidischem Algorithmus
int rhsR, anzahlNachkommastellen, globalAnzahlNachkommastellen, grGemTeiler;
double nenner, ergebnisKehrwert;
std::vector<double> rLCDVec;
for(int g=0;g<conVecU.size();g++){
	globalAnzahlNachkommastellen=0;
	grGemTeiler=0;
	std::vector<double> keineInt;
	std::vector<data::IndexedElement> lhs=conVecU[g]->getElements();
	for(int l=0;l<lhs.size();l++){
		if(lhs[l].value.asDouble()-floor(lhs[l].value.asDouble())!=0){
			keineInt.push_back(abs(lhs[l].value.asDouble()));
		}
	}
	if(rhsVecU[g]->getValue().asDouble()-floor(rhsVec[g]->getValue().asDouble())!=0){
		keineInt.push_back(abs(rhsVec[g]->getValue().asDouble()));
	}
	if(keineInt.size()==0){
		rLCDVec.push_back(1);

	}
	else{
		nenner=1.0;
		for(int o=0;o<keineInt.size();o++){
			anzahlNachkommastellen=0;
			while(keineInt[o]-floor(keineInt[o]) !=0){
				keineInt[o]=10*keineInt[o];
				anzahlNachkommastellen ++;
			}
			if(anzahlNachkommastellen>globalAnzahlNachkommastellen){
				globalAnzahlNachkommastellen=anzahlNachkommastellen;
			}

		}
		int r=1, y, z, ggt1;
        y=keineInt[0]; //Bestimmung des ggT's von den ersten beiden Werten
        z=keineInt[1];
        if(z>y){ //Weil beim Euklid gelten muss: m>n
        	z=y;
        	y=keineInt[1];
        }
        while(r!=0) { //Euklid
        	r=y%z;
        	y=z;
        	if(r==0){
            	ggt1=z;
        	}
            z=r;
        }
        for(int q=2;q<=keineInt.size();q++){ //Fortlaufend wird jeweils ggT(ggT(a,b),c) bestimmt
            if(ggt1<keineInt[q]){ //Es muss gelten m>n
                z=ggt1;
                y=keineInt[q];
            }
            else{
                y=ggt1;
                z=keineInt[q];
            }
            while(r!=0){ //Euklid
                r=y%z;
                y=z;
                if(r==0){
                	ggt1=z;
                	r=1;
                }
                z=r;
            }
        }
		for(int a=0;a<globalAnzahlNachkommastellen;a++){
			nenner=nenner*10;
		}
		ergebnisKehrwert=abs(ggt1)/nenner;
		rLCDVec.push_back(ergebnisKehrwert);
	}
}
		 */
		//Ausgabe QIP

		//Schreibwerkzeug, Datei erstellen und öffnen
		fstream f;
		string Outputfile=inputfile+".reduced";
		f.open(Outputfile, ios::out);
		//Zielfunktion
		if (zielfkt.getObjective()==data::QpObjFunc::Objective::min)
			f<<"MINIMIZE"<<endl;
		else f<<"MAXIMIZE"<<endl;
		for(int h=0;h<zielwerte.size();h++){
		    if(zielwerte[h].value.asDouble()>=0)
			f<<"+";
		    else f<<"-";
		    f<<fabs(zielwerte[h].value.asDouble())<<vecVariablen[zielwerte[h].index]->getName()<<" ";
		}
		if(mTilde>=0)
			f<<"-"<<fabs(mTilde)<<"p"<<endl;
		else 
			f<<"+"<<fabs(mTilde)<<"p"<<endl;
		//Restriktionen
		f<<"SUBJECT TO"<<endl;

		//Vektor mit v anstatt x erstellen
		string varBez;
		std::vector<string> vVec;
		for(int i=0;i<vecVariablen.size();i++){
			varBez="v_"+vecVariablen[i]->getName();
			vVec.push_back(varBez);
		}

		//Teil1
		int indE;
		for (int b=0;b<numEBlocks;b++){
			if(ersterBlock==0){
				indE=2*b+1;
			}
			else{
				indE=2*(b+1);
			}
			std::vector<int>::iterator it = mKVec.begin();
			for(int a=0;a<conVecE.size();a++){
				std::vector<data::IndexedElement> lhs=conVecE[a]->getElements();

				if (conVecE[a]->getRhs().getRatioSign() == data::QpRhs::smallerThanOrEqual || conVecE[a]->getRhs().getRatioSign() == data::QpRhs::equal){
					for(int c=0;c<lhs.size();c++){
						if(lhs[c].value.asDouble()>=0){
							f<<"+";
						}
						if(qbp->getBlock(lhs[c].index)<=indE) f<<lhs[c].value.asDouble()<<vecVariablen[lhs[c].index]->getName()<<" ";
						else f<<lhs[c].value.asDouble()<<vVec[lhs[c].index]<<"_"<<indE<<" ";
					}
					if(indE-1>0 && *it!=0){
						f<<-*it<<"t"<<indE-1;
					}
					it++;
					f<<" <= "<<rhsVecE[a]->getValue().asDouble()<<endl;
				}
				if (conVecE[a]->getRhs().getRatioSign() == data::QpRhs::greaterThanOrEqual || conVecE[a]->getRhs().getRatioSign() == data::QpRhs::equal){
					for(int c=0;c<lhs.size();c++){
						if(lhs[c].value.asDouble()<=0){
							f<<"+";
						}
						else  f<<"-";
						if(qbp->getBlock(lhs[c].index)<=indE) f<<fabs(lhs[c].value.asDouble())<<vecVariablen[lhs[c].index]->getName()<<" ";
						else f<<fabs(lhs[c].value.asDouble())<<vVec[lhs[c].index]<<"_"<<indE<<" ";
					}
					if(indE-1>0 && -*it!=0){
						f<<-*it<<"t"<<indE-1;
					}
					it++;
					f<<" <= "<<-rhsVecE[a]->getValue().asDouble()<<endl;
				}
			}
			assert(it==mKVec.end());
		}

		//Teil2
		std::vector<int>::iterator it = mKVec.begin();

		for(int d=0;d<conVecE.size();d++){

			std::vector<data::IndexedElement> lhs=conVecE[d]->getElements();

			if (conVecE[d]->getRhs().getRatioSign() == data::QpRhs::smallerThanOrEqual || conVecE[d]->getRhs().getRatioSign() == data::QpRhs::equal){
				for (int e=0;e<lhs.size();e++){
					if(lhs[e].value.asDouble()>=0){
						f<<"+";
					}
					f<<lhs[e].value.asDouble()<<vecVariablen[lhs[e].index]->getName()<<" ";
				}
				if(-*it!=0){
					f<<-*it<<"p ";
				}
				it++;
				f<<"<= "<< (rhsVecE[d]->getValue().asDouble())<<endl;
			}
			if (conVecE[d]->getRhs().getRatioSign() == data::QpRhs::greaterThanOrEqual || conVecE[d]->getRhs().getRatioSign() == data::QpRhs::equal){
				for(int c=0;c<lhs.size();c++){
					if(lhs[c].value.asDouble()<0){
						f<<"+";
					}
					else  f<<"-";
					f<<abs(lhs[c].value.asDouble())<<vecVariablen[lhs[c].index]->getName()<<" ";
				}
				if(-*it!=0){
					f<<-*it<<"p ";
				}
				it++;
				f<<"<= " << (-rhsVecE[d]->getValue().asDouble())<<endl;
			}
		}
		assert(it==mKVec.end());


		//Teil3
		int indU;
		for(int b=0;b<numABlocks+(!FinalBlockIsUniversal);b++){
			if(ersterBlock==1){
				indU=2*b+1;
			}
			else{
				indU=2*(b+1);
			}
			std::vector<int>::iterator it = lKVec.begin();
			for(int a=0;a<conVecU.size();a++){
				std::vector<data::IndexedElement> lhs=conVecU[a]->getElements();
				if (conVecU[a]->getRhs().getRatioSign() == data::QpRhs::smallerThanOrEqual || conVecU[a]->getRhs().getRatioSign() == data::QpRhs::equal){
					for(int c=0;c<lhs.size();c++){
						if(-lhs[c].value.asDouble()>0) f<<"+";
						if(qbp->getBlock(lhs[c].index)<=indU) f<<-lhs[c].value.asDouble()<<vecVariablen[lhs[c].index]->getName()<<" ";
						else f<<-lhs[c].value.asDouble()<<vVec[lhs[c].index]<<"_"<<indU<<" ";

					}
					if(-(*it-rhsVecU[a]->getValue().asDouble()-rLCDVec[a])!=0){
						if(-(*it-rhsVecU[a]->getValue().asDouble()-rLCDVec[a])>0){
							f<<"+";
						}
						f<<-(*it-rhsVecU[a]->getValue().asDouble()-rLCDVec[a])<<"CheckAllC"<<a+1<<"_InS"<<indU<<" ";
					}
					f<<"<= "<<-*it<<endl;
					it++;
				}

				if (conVecU[a]->getRhs().getRatioSign() == data::QpRhs::greaterThanOrEqual || conVecU[a]->getRhs().getRatioSign() == data::QpRhs::equal){
					for(int c=0;c<lhs.size();c++){
						if(lhs[c].value.asDouble()>0) f<<"+";
						if(qbp->getBlock(lhs[c].index)<=indU) f<<lhs[c].value.asDouble()<<vecVariablen[lhs[c].index]->getName()<<" ";
						else f<<lhs[c].value.asDouble()<<vVec[lhs[c].index]<<"_"<<indU<<" ";

					}
					if(-(*it+rhsVecU[a]->getValue().asDouble()-rLCDVec[a])!=0){
						if(-(*it+rhsVecU[a]->getValue().asDouble()-rLCDVec[a])>0){
							f<<"+";
						}
						f<<-(*it+rhsVecU[a]->getValue().asDouble()-rLCDVec[a])<<"CheckAllC"<<a+1<<"_InS"<<indU<<" ";
					}
					f<<"<= "<<-*it<<endl;
					it++;
				}
			}
			assert(it==lKVec.end());
		}

		//Teil4
		f<<"+p";
		for(int a=0;a<numABlocks;a++){
			if(ersterBlock==1){
				indU=2*a+1;
			}
			else{
				indU=2*(a+1);
			}
			f<<" -t"<<indU;
		}
		if(!FinalBlockIsUniversal)
			f << " -t_final";
		f<<" <= "<<0<<endl;

		//Teil5
		for(int a=0;a<numABlocks+!FinalBlockIsUniversal;a++){
			if(ersterBlock==1){
				indU=2*a+1;
			}
			else{
				indU=2*(a+1);
			}
			if(a==numABlocks)
				f << " +t_final";
			else
				f<<"+t"<<indU;
			if(indU-2>0){
				f<<" -t"<<indU-2;
			}
			for(int b=1;b<=anzahlUCon;b++){
				f<<" -CheckAllC"<<b<<"_InS"<<indU;
			}
			f<<" <= "<<0<<endl;
		}

		//orderVec sortieren in finalOrderVec, x-Variablen hinzufügen, Bounds in Datei schreiben und Behilfs-Vektoren erstellen
		f<<"BOUNDS"<<endl;
		std::vector<string> finalOrderVec, finalExVec, finalAllVec, finalContVec;

		bool UniversalBlock=(bool)ersterBlock;
		bool SetAllXes =true;
		assert(UniversalBlock==vecVariablen[0]->getQuantifier());
		int BlockStart=0;
		int Next=0;
		for(int i=0;i<numBlocks;i++){
			if(vecVariablen[BlockStart]->getNumberSystem() ==data::QpVar::NumberSystem::real){
				cerr << "Tried to reduce QIP with universal system having non-integer variable " << vecVariablen[BlockStart]->getName() << endl;
				assert(i==numBlocks-1);
			}
			assert(BlockStart==0 || vecVariablen[BlockStart]->getQuantifier()!=UniversalBlock);
			if(BlockStart!=0){ //1=All; 0=Existenz
				UniversalBlock=!UniversalBlock;
			}
			Next=0;

			//The original variables
			while (BlockStart+Next<n &&UniversalBlock==vecVariablen[BlockStart+Next]->getQuantifier()){
				if((qbp->getLowerBound(BlockStart+Next)==0 && qbp->getUpperBound(BlockStart+Next)==1))finalOrderVec.push_back(vecVariablen[BlockStart+Next]->getName());
				else finalContVec.push_back(vecVariablen[BlockStart+Next]->getName());
				f<<qbp->getLowerBound(BlockStart+Next)<<" <= "<<vecVariablen[BlockStart+Next]->getName()<<" <= "<<qbp->getUpperBound(BlockStart+Next)<<endl;
				assert((qbp->getLowerBound(BlockStart+Next)==0 && qbp->getUpperBound(BlockStart+Next)==1) || (qbp->getQuantifier(BlockStart+Next)==EXIST && vecVariablen[BlockStart+Next]->getNumberSystem() ==data::QpVar::NumberSystem::real));
				if (vecVariablen[BlockStart+Next]->getQuantifier()==UNIV)
					finalAllVec.push_back(vecVariablen[BlockStart+Next]->getName());
				else finalExVec.push_back(vecVariablen[BlockStart+Next]->getName());
				Next++;
			}

			BlockStart=BlockStart+Next;
			Next=0;
			while (BlockStart+Next<n){
				int ind= BlockStart+Next;
				if((IsInA[ind] && UniversalBlock)||(IsInE[ind] & !UniversalBlock)){
					if((qbp->getLowerBound(BlockStart+Next)==0 && qbp->getUpperBound(BlockStart+Next)==1)) finalOrderVec.push_back(vVec[ind]+"_"+to_string(i+1));
					else finalContVec.push_back(vVec[ind]+"_"+to_string(i+1));
					f<<qbp->getLowerBound(ind)<<" <= "<<vVec[ind]+"_"+to_string(i+1)<<" <= "<<qbp->getUpperBound(ind)<<endl;
					if (UniversalBlock) finalAllVec.push_back(vVec[ind]+"_"+to_string(i+1));
					else finalExVec.push_back(vVec[ind]+"_"+to_string(i+1));
				}
				Next++;
			}
			if(UniversalBlock){
				//y-Variablen
				for(int con=0;con<rhsVecU.size();con++){
					string yVar="CheckAllC"+to_string(con+1)+"_InS"+to_string(i+1);
					finalOrderVec.push_back(yVar);
					f<<"0"<<" <= "<<yVar<<" <= "<<"1"<<endl;
					finalExVec.push_back(yVar);
				}
				//t-Variablen
				string tVar ="t"+to_string(i+1);
				finalOrderVec.push_back(tVar);
				f<<"0"<<" <= "<<tVar<<" <= "<<"1"<<endl;
				finalExVec.push_back(tVar);
			}
			if(i==numBlocks-1 && !FinalBlockIsUniversal){
				for(int con=0;con<rhsVecU.size();con++){
					string yVar="CheckAllC"+to_string(con+1)+"_InS"+to_string(i+2);
					finalOrderVec.push_back(yVar);
					f<<"0"<<" <= "<<yVar<<" <= "<<"1"<<endl;
					finalExVec.push_back(yVar);
				}

				assert(!UniversalBlock);
				string tVar ="t_final";
				finalOrderVec.push_back(tVar);
				f<<"0"<<" <= "<<tVar<<" <= "<<"1"<<endl;
				finalExVec.push_back(tVar);
			}
		}

		//p hinzufügen
		string pVar ="p";
		finalOrderVec.push_back(pVar);
		f<<"0"<<" <= "<<pVar<<" <= "<<"1"<<endl;
		finalExVec.push_back(pVar);

		/*//GENERAL
		f<<"GENERAL"<<endl;
		for(int i=0;i<finalGenVVec.size();i++){
			f<<finalGenVVec[i]<<" ";
		}
		f<<endl;*/

		//BINARY
		f<<"BINARY"<<endl;
		for(int j=0;j<finalOrderVec.size();j++){
			f<<finalOrderVec[j]<<" ";
		}
		f<<endl;

		//EXISTS
		f<<"EXISTS"<<endl;
		//for(int i=0;i<vecE.size();i++){
		//	f<<vecE[i]->getName()<<" ";
		//}
		for (int b=0;b<finalExVec.size();b++){
			f<<finalExVec[b]<<" ";
		}
		for (int b=0;b<finalContVec.size();b++){
                        f<<finalContVec[b]<<" ";
                }
		f<<endl;

		//ALL
		f<<"ALL"<<endl;
		//for(int j=0;j<vecU.size();j++){
		//	f<<vecU[j]->getName()<<" ";
		//}
		for (int b=0;b<finalAllVec.size();b++){
			f<<finalAllVec[b]<<" ";
		}
		f<<endl;

		//ORDER
		f<<"ORDER"<<endl;
		for(int i=0;i<finalOrderVec.size();i++){
			f<<finalOrderVec[i]<<" ";
		}
		for (int b=0;b<finalContVec.size();b++){
                        f<<finalContVec[b]<<" ";
                }
		f<<endl;

		//END
		f<<"END"<<endl;


		f.close();
		 cerr << "Reduced QIP written to " << Outputfile << endl;
		//cin.get();
	}
    int UcompV(std::vector<data::IndexedElement > &U, std::vector<data::IndexedElement > &V, data::IndexedElement rhs_u, data::IndexedElement rhs_v) // prepared to be "<="
    {
        if (U.size() == 0 && V.size() == 0) {
            assert(rhs_u.value.asDouble() >= 0.0);
            assert(rhs_v.value.asDouble() >= 0.0);
            return 0;
        }
        if (U.size() == 0) {
            assert(rhs_u.value.asDouble() >= 0.0);
            if (V[0].value.asDouble() < 0.0) return 1;
            else if (V[0].value.asDouble() > 0.0) return -1;
        }
        if (V.size() == 0) {
            assert(rhs_v.value.asDouble() >= 0.0);
            if (U[0].value.asDouble() < 0.0) return -1;
            else if (U[0].value.asDouble() > 0.0) return 1;
        }
        int i = 0; int j = 0;
        while (i < U.size() && j < V.size()) {
            double nextUcoef = 0.0;
            double nextVcoef = 0.0;
            if (U[i].index < V[j].index) {
                nextUcoef = U[i].value.asDouble();
            } else if (U[i].index > V[j].index) {
                nextVcoef = V[j].value.asDouble();
            } else { // equality
                nextUcoef = U[i].value.asDouble();
                nextVcoef = V[j].value.asDouble();
            }

            if (nextUcoef < nextVcoef -  1e-12) return -1;
            else if (nextUcoef > nextVcoef + 1e-12) return 1;
            else {
                i++; j++;
            }
        }
        if (U.size() < V.size()) {
            int iv = U.size();
            while (iv < V.size() && fabs( V[iv].value.asDouble() ) < LP_EPS ) iv++;
            if (iv >= V.size()) return 0;
            if (V[iv].value.asDouble() > 0.0) return -1;
            else if (V[iv].value.asDouble() < 0.0) return 1;
        }
        if (U.size() > V.size()) {
            int iu = V.size();
            while (iu < U.size() && fabs( U[iu].value.asDouble() ) <LP_EPS ) iu++;
            if (iu >= U.size()) return 0;
            if (U[iu].value.asDouble() > 0.0) return 1;
            else if (U[iu].value.asDouble() < 0.0) return -1;
        }
        if (rhs_u.value.asDouble() > rhs_v.value.asDouble()) return -1;
        if (rhs_u.value.asDouble() < rhs_v.value.asDouble()) return 1;
        return 0;//x.value.asDouble() > y.value.asDouble();//false;//x.index < y.index;
    }

    bool UltV(std::vector<data::IndexedElement > &U, std::vector<data::IndexedElement > &V, data::IndexedElement rhs_u, data::IndexedElement rhs_v)
    {
        if (XcompY(U,V,rhs_u,rhs_v) == -1) return true;
        return false;
    }
    
    bool UgtV(std::vector<data::IndexedElement > &U, std::vector<data::IndexedElement > &V, data::IndexedElement rhs_u, data::IndexedElement rhs_v)
    {
        if (XcompY(U,V,rhs_u,rhs_v) == 1) return true;
        return false;
    }
    
#ifdef NEXT_DOM
    int UcompV(std::vector<data::IndexedElement > &U, std::vector<data::IndexedElement > &V, data::IndexedElement rhs_u, data::IndexedElement rhs_v) // prepared to be "<="
    {
        if (U.size() == 0 && V.size() == 0) {
            //per definition no subsumption and no domination
            return 0;
        }
        // now assumed: (<U> >0 or <V> >0 ) and if one of them == 0 both rhs >= 0
        if (U.size() == 0) {
            assert(rhs_u >= 0.0);
            // then U redundant anyway
            return 1;
        }
        if (V.size() == 0) {
            assert(rhs_v >= 0.0);
            // then V redundant anyway
            return -1;
        }
        int i = 0; int j = 0;
        bool usharper = false;
        bool vsahrper = false;
        while (i < U.size() || j < V.size()) {
            double nextUcoef = 0.0;
            double nextVcoef = 0.0;
            if (U[i].index < V[j].index) {
                nextUcoef = U[i].value.asDouble();
            } else if (X[i].index < Y[j].index) {
                nextVcoef = V[j].value.asDouble();
            } else { // equality
                nextUcoef = U[i].value.asDouble();
                nextVcoef = V[j].value.asDouble();
            }
            
            if (rhs_u < rhs_v) {  // u sharper
                if (nextUcoef < nextVcoef) return 0;
            } else if (rhs_u > rhs_v) { // v sharper
                if (nextUcoef > nextVcoef) return 0;
            } else { //equality
                if (nextUcoef > nextVcoef) usharper = true;
                if (nextUcoef < nextVcoef) vsharper = true;
                if (usharper && vsharper) return 0;
            }
            
            if (X[i].index < Y[j].index) { if (X[i].value.asDouble() < 0.0) return 1; else return -1; }
            if (X[i].index > Y[j].index) { if (Y[i].value.asDouble() < 0.0) return -1; else return 1; }
            if (X[i].index == Y[j].index) {
                if (X[i].value.asDouble() < Y[j].value.asDouble() - 1e-12) {
                    return 1;
                }
                else if (X[i].value.asDouble() > Y[j].value.asDouble() + 1e-12) {
                    return -1;
                }
                i++; j++;
            }
        }
        if (Y.size() < X.size()) return -1;
        if (Y.size() > X.size()) return 1;
        if (x.value.asDouble() > y.value.asDouble()) return -1;
        if (x.value.asDouble() < y.value.asDouble()) return 1;
        return 0;//x.value.asDouble() > y.value.asDouble();//false;//x.index < y.index;
    }
#endif
    int XcompY(std::vector<data::IndexedElement > &X, std::vector<data::IndexedElement > &Y, data::IndexedElement x, data::IndexedElement y)
	{
		if (X.size() == 0 && Y.size() == 0) {
			//cerr << "X." << X.size() << "<=" << "Y." << Y.size() << endl;
			return 0;
		}
		if (X.size() == 0) {
			//cerr << "X." << X.size() << "<=" << "Y" << Y.size() << endl;
            if (Y[0].value.asDouble() > 0.0) return 1;
			else if (Y[0].value.asDouble() < 0.0) return -1;
		}
		if (Y.size() == 0) {
			//cerr << "X" << X.size() << "<=" << "Y." << Y.size() << endl;
            if (X[0].value.asDouble() > 0.0) return -1;
            else if (X[0].value.asDouble() < 0.0) return 1;
		}
		int i = 0; int j = 0;
		while (i < X.size() && j < Y.size()) {
			if (X[i].index < Y[j].index) { if (X[i].value.asDouble() < 0.0) return 1; else return -1; }
			if (X[i].index > Y[j].index) { if (Y[i].value.asDouble() < 0.0) return -1; else return 1; }
			if (X[i].index == Y[j].index) {
				if (X[i].value.asDouble() < Y[j].value.asDouble() - 1e-12) {
					return 1;
				}
				else if (X[i].value.asDouble() > Y[j].value.asDouble() + 1e-12) {
					return -1;
				}
				i++; j++;
			}
		}
        if (Y.size() < X.size()) {
            int ix = Y.size();
            while (ix < X.size() && fabs( X[ix].value.asDouble() ) < LP_EPS ) ix++;
            if (ix >= X.size()) return 0;
            if (X[ix].value.asDouble() > 0.0) return -1;
            else if (X[ix].value.asDouble() < 0.0) return 1;
        }
        if (Y.size() > X.size()) {
            int iy = X.size();
            while (iy < Y.size() && fabs( Y[iy].value.asDouble() ) < LP_EPS ) iy++;
            if (iy >= Y.size()) return 0;
            if (Y[iy].value.asDouble() > 0.0) return 1;
            else if (Y[iy].value.asDouble() < 0.0) return -1;
        }
		if (x.value.asDouble() > y.value.asDouble()) return -1;
		if (x.value.asDouble() < y.value.asDouble()) return 1;
		return 0;//x.value.asDouble() > y.value.asDouble();//false;//x.index < y.index;
	}

	bool XltY(std::vector<data::IndexedElement > &X, std::vector<data::IndexedElement > &Y, data::IndexedElement x, data::IndexedElement y)
	{
		if (XcompY(X,Y,x,y) == -1) return true;
		return false;
	}

	bool XgtY(std::vector<data::IndexedElement > &X, std::vector<data::IndexedElement > &Y, data::IndexedElement x, data::IndexedElement y)
	{
		if (XcompY(X,Y,x,y) == 1) return true;
		return false;
	}

	void yIselectionSort_Cols(int *colsort, int size, int maxLPstage)
	{
	    int     i, j, best_i;
	    int     tmp;
		extSol::QpExternSolver& eS = QlpStSolve->getExternSolver(maxLPstage);

	    for (i = 0; i < size-1; i++){
	        best_i = i;
	        for (j = i+1; j < size; j++){
	        	if ( XltY(*(eS.getCol_snapshot(colsort[j])), *(eS.getCol_snapshot(colsort[best_i])),
	        			  eS.getObj_snapshot(colsort[j]), eS.getObj_snapshot(colsort[best_i])) )
	                best_i = j;
	        }
	        tmp = colsort[i]; colsort[i] = colsort[best_i]; colsort[best_i] = tmp;
	    }
	}

	void yIsort(int* colsort, int size, int maxLPstage)
	{
		extSol::QpExternSolver& eS = QlpStSolve->getExternSolver(maxLPstage);
	    if (size <= 15)
	        yIselectionSort_Cols(colsort, size, maxLPstage);

	    else{
	        int         pivot_ix = size / 2;
	        int         pivot_col = colsort[pivot_ix];
	        int         tmp;
	        int         i = -1;
	        int         j = size;

	        for(;;){
	            do i++; while( XltY(*(eS.getCol_snapshot(colsort[i])), *(eS.getCol_snapshot(pivot_col)),
	            			        eS.getObj_snapshot(colsort[i]), eS.getObj_snapshot(pivot_col)) );
	            //do j--; while( XgtY(*(eS.getCol_snapshot(colsort[j/*pivot_ix*/])), *(eS.getCol_snapshot(pivot_col)),
    			//                    eS.getObj_snapshot(colsort[/*pivot_ix*/j]), eS.getObj_snapshot(pivot_col)) );
	            do j--; while( XltY(*(eS.getCol_snapshot(pivot_col)), *(eS.getCol_snapshot(colsort[j/*pivot_ix*/])),
    			                    eS.getObj_snapshot(pivot_col), eS.getObj_snapshot(colsort[j/*pivot_ix*/])) );

	            if (i >= j) break;

	            tmp = colsort[i]; colsort[i] = colsort[j]; colsort[j] = tmp;
	        }

	        yIsort(colsort    , i     , maxLPstage);
	        yIsort(&colsort[i], size-i, maxLPstage);
	    }
	}

	void sortCols(int maxLPstage) {
		int n = getnVars();
		sortcols.resize(n);
		for (int i = 0; i < n; i++) sortcols[i] = i;
		yIsort( sortcols.data(), sortcols.size(), maxLPstage);
        std::vector<int> tmp(sortcols.size());
        if(0)for (int i = 0; i < sortcols.size();i++) {
        	int x = sortcols[i];
        	tmp[i] = x;//.push_back(x);
        }
        if(0)for (int i = 0; i < tmp.size();i++)
        	sortcols[tmp[i]] = i;
	}

	/*
	int XdomY(extSol::QpExternSolver& externSolver, std::vector<data::IndexedElement > *X, std::vector<data::IndexedElement > *Y, data::IndexedElement x, data::IndexedElement y, bool useObj, std::vector<data::QpRhs> *RHSs)
	{
		// X and Y must be sorted with ascending index
		// was ist mit <=,>=,== ? ---> sense must be <= or == and obj minimize
		bool XdomY=true;
		bool YdomX=true;
		int n = getnVars();

		if (useObj) {
			if (x.index < n && y.index < n) {
				if (x.value.asDouble() < y.value.asDouble() - 1e-12) XdomY = false;
				else if (x.value.asDouble() > y.value.asDouble() + 1e-12) YdomX = false;
			} else if (x.index < n) {
				YdomX = false;
			} else if (y.index < n) {
				XdomY = false;
			} else useObj = false;
		}

		int i = 0; int j = 0;
		while (i < (*X).size() && j < (*Y).size()) {
			while (i < (*X).size() && j < (*Y).size() && (*X)[i].index < (*Y)[j].index ) { i++; YdomX = false; }
			if (i >= (*X).size() && j < (*Y).size()) XdomY = false;
			while (i < (*X).size() && j < (*Y).size() && (*X)[i].index > (*Y)[j].index ) { j++; XdomY = false; }
			if (i < (*X).size() && j >= (*Y).size()) YdomX = false;
			if (i >= (*X).size() || j >= (*Y).size()) {
				;
			} else if ((*X)[i].index == (*Y)[j].index) {
				assert((*RHSs)[(*X)[i].index].getRatioSign() != data::QpRhs::RatioSign::greaterThanOrEqual);
				if ((*X)[i].value.asDouble() < (*Y)[j].value.asDouble() - 1e-12) {
					XdomY = false;
					if ((*RHSs)[(*X)[i].index].getRatioSign() == data::QpRhs::RatioSign::equal) return 3;
				}
				else if ((*X)[i].value.asDouble() > (*Y)[j].value.asDouble() + 1e-12) {
					YdomX = false;
					if ((*RHSs)[(*X)[i].index].getRatioSign() == data::QpRhs::RatioSign::equal) return 3;
				}
				i++; j++;
				if (i >= (*X).size() && j < (*Y).size()) XdomY = false;
				if (i < (*X).size() && j >= (*Y).size()) YdomX = false;
			}
			if (!XdomY && !YdomX) break;
		}
		if (XdomY && YdomX) return 0;
		if (XdomY) return 1;
		if (YdomX) return -1;
		return 2;
	}
	*/

	int XdomY(extSol::QpExternSolver& externSolver, std::vector<data::IndexedElement > *X, std::vector<data::IndexedElement > *Y, data::IndexedElement x,
			data::IndexedElement y, bool useObj, std::vector<data::QpRhs> *RHSs, int colx, int coly)
	{
		// X and Y must be sorted with ascending index
		// was ist mit <=,>=,== ? ---> sense must be <= or == and obj minimize
		bool XdomY=true;
		bool YdomX=true;
		int n = getnVars();
		double xVal, yVal;
		int xInd, yInd;

		if (0&&colx == 90 && coly == 91) {
			cerr << "Rows with x91:";
			for (int i = 0; i < (*X).size();i++)
				cerr << (*X)[i].value.asDouble() << "u" << (*X)[i].index << " ";
			cerr << endl;
			cerr << "Rows with x92:";
			for (int i = 0; i < (*Y).size();i++)
				cerr << (*Y)[i].value.asDouble() << "u" << (*Y)[i].index << " ";
			cerr << endl;
		}
		if (0&&colx == 91 && coly == 90) {
			cerr << "Rows with x92:";
			for (int i = 0; i < (*X).size();i++)
				cerr << (*X)[i].value.asDouble() << "u" << (*X)[i].index << " ";
			cerr << endl;
			cerr << "Rows with x91:";
			for (int i = 0; i < (*Y).size();i++)
				cerr << (*Y)[i].value.asDouble() << "u" << (*Y)[i].index << " ";
			cerr << endl;
		}

		if (useObj) {
			xInd = x.index; yInd = y.index;
			if (x.index < n && y.index < n) {
				xVal = x.value.asDouble(); yVal = y.value.asDouble();
				if (xVal < yVal) XdomY = false;
				else if (xVal > yVal) YdomX = false;
			} else if (x.index < n) {
				xVal = x.value.asDouble(); yVal = 0.0;
				if (xVal < yVal) XdomY = false;
				else if (xVal > yVal) YdomX = false;
			} else if (y.index < n) {
				xVal = 0.0; yVal = y.value.asDouble();
				if (xVal < yVal) XdomY = false;
				else if (xVal > yVal) YdomX = false;
			} else useObj = false;
		}

		int i = 0; int j = 0;
		bool eqOc = false;
		while (i < (*X).size() || j < (*Y).size()) {
			if (i < (*X).size() && (*RHSs)[(*X)[i].index].getRatioSign() == data::QpRhs::RatioSign::equal) eqOc = true;
			if (j < (*Y).size() && (*RHSs)[(*Y)[j].index].getRatioSign() == data::QpRhs::RatioSign::equal) eqOc = true;
			if (i < (*X).size() && j < (*Y).size()) {
				if ((*X)[i].index < (*Y)[j].index) {
					xVal = (*X)[i].value.asDouble();
					yVal = 0.0;
					if (xVal < yVal) XdomY = false;
					else if (xVal > yVal) YdomX = false;
					i++;
				} else if ((*X)[i].index > (*Y)[j].index) {
					xVal = 0.0;
					yVal = (*Y)[j].value.asDouble();
					if (xVal < yVal) XdomY = false;
					else if (xVal > yVal) YdomX = false;
					j++;
				} else if ((*X)[i].index == (*Y)[j].index) {
					if (0&&colx == 90 && coly == 91) {
						cerr << (*X)[i].value.asDouble() << "x90 und " << (*Y)[j].value.asDouble() << "x91" << endl;
					}
					if (0&&colx == 91 && coly == 90) {
						cerr << (*X)[i].value.asDouble() << "x91 und " << (*Y)[j].value.asDouble() << "x90" << endl;
					}
					assert((*RHSs)[(*X)[i].index].getRatioSign() != data::QpRhs::RatioSign::greaterThanOrEqual);
					if ((*X)[i].value.asDouble() < (*Y)[j].value.asDouble() ) {
						XdomY = false;
					}
					else if ((*X)[i].value.asDouble() > (*Y)[j].value.asDouble() ) {
						YdomX = false;
					}
					i++; j++;
				}
			} else if (i < (*X).size()) {
				xVal = (*X)[i].value.asDouble();
				yVal = 0.0;
				if (xVal < yVal) XdomY = false;
				else if (xVal > yVal) YdomX = false;
				i++;
			} else if (j < (*Y).size()) {
				xVal = 0.0;
				yVal = (*Y)[j].value.asDouble();
				if (xVal < yVal) XdomY = false;
				else if (xVal > yVal) YdomX = false;
				j++;
			}
			if (!XdomY && !YdomX) break;
		}

		if (XdomY && YdomX) return 0;
		if (eqOc) return 3;
		if (XdomY) return 1;
		if (YdomX) return -1;
		return 2;
	}

	void findSymmetries(utils::QlpStageSolver& QlpStSolve, int maxLPstage, bool useObj, std::vector< std::pair<int,int> > &clist, int *types, int *block, int8_t * assigns, int *eas) {
#ifndef FIND_BUG
	  clist.clear();
	  //return;
#endif //#ifdef FIND_BUG
	  sortCols(maxLPstage);
	  clist.clear();
	  if(qbp->getUniversalConstraintsExist()) return;
	  //for (int i = 0; i < sortcols.size();i++) cerr << " " << sortcols[i];
	  //cerr << endl;
	  int n = getnVars();
	  int minBlock = n+10;
	  for (int i = 0; i < n;i++)
	    if (block[i] < minBlock) minBlock = block[i];
	  time_t start = time(NULL);
	  time_t ini_time = qbp->getInitilizationTime();
	  CliqueManager *pCM = qbp->getClipueManager();
	  int outputs=0;
	  int outputsSYMM=0;

	  //cerr << "ENTER FIND SYMM (time(NULL) - ini_time):" << (time(NULL) - ini_time) << endl;

	  int latestHitI=-1;
	  for (int i = 0/*1*/; i < n/*10*/ || (i < n /*&& time(NULL) - start < (time(NULL) - ini_time) * 0.2 + 5*/);i++) {
	    //cerr << "LOOP " << k << " ENTER FIND SYMM (time(NULL) - ini_time):" << (time(NULL) - ini_time) << " and time(NULL) - start:" << time(NULL) - start << endl;
	    if (latestHitI < i - 1000 - (int)sqrt((double)n)) break;
	    int latestHitIplusK=i+1-1;
	    for (int k = 1; k < n-i/*-k*/;k++) {
	      //int j = i+k;
	      int iplusk = i+k;
	      if(i==iplusk) continue;
	      if (latestHitIplusK < iplusk - 100 - (int)sqrt((double)n)) break;
	      if (QlpStSolve.getExternSolver(maxLPstage).getCol_snapshot(sortcols[i])->size() == 0) continue;
	      if (QlpStSolve.getExternSolver(maxLPstage).getCol_snapshot(sortcols[iplusk])->size() == 0) continue;
	      //if (block[sortcols[i]] != minBlock || block[sortcols[iplusk]] != minBlock) continue;
	      if (block[sortcols[i]] != block[sortcols[iplusk]]) continue;
	      if (types[sortcols[i]] != BINARY || types[sortcols[iplusk]] != BINARY) continue;
	      if (eas[sortcols[i]] == UNIV || eas[sortcols[iplusk]] == UNIV) continue;
	      if (assigns[sortcols[i]] != 2 || assigns[sortcols[iplusk]] != 2) continue;
	      //if (QlpStSolve.getExternSolver(maxLPstage).getCol_snapshot(sortcols[i])->size() > n-2) continue;
	      int xdy = XdomY( QlpStSolve.getExternSolver(maxLPstage),
			       QlpStSolve.getExternSolver(maxLPstage).getCol_snapshot(sortcols[i]),QlpStSolve.getExternSolver(maxLPstage).getCol_snapshot(sortcols[iplusk]),
			       QlpStSolve.getExternSolver(maxLPstage).getObj_snapshot(sortcols[i]),QlpStSolve.getExternSolver(maxLPstage).getObj_snapshot(sortcols[iplusk]),
			       useObj, QlpStSolve.getExternSolver(maxLPstage).getRowRhs_snapshot(), sortcols[i], sortcols[iplusk]);
	      if (/*k == 1 &&*/ xdy==0) {
		if (sortcols[i] >= sortcols[iplusk]) continue;
		latestHitI = i;
		latestHitIplusK = iplusk;
		//if (outputsSYMM < 5) cerr << "info SYMMETRIE FOUND!! " << i << "," << iplusk << " | " << sortcols[i] << ", " << sortcols[iplusk] << endl;
		//else if (outputsSYMM==5) cerr << " ... " << endl;
		outputsSYMM++;
		int j;
		bool foundFixing = false;
		j = pCM->FirstAdjacentInConflictGraph(2*sortcols[i]); // Fall x_i <= x_j
		if (j >= 0) {
		  //cerr << "j=" << j << endl;
		  int kk = pCM->NextAdjacentInConflictGraph(j);
		  while (kk >= 0) {
		    //cerr << "j=" << j << " kk=" << kk << " a(j)=" << CM.getAdjacent(j) << " a(kk)=" << CM.getAdjacent(kk) << endl;
		    int v1 = sortcols[i];
		    int v2 = pCM->getAdjacent(kk);
		    v2 /= 2;
		    assert(types[v2] != /*CONTINUOUS*/5000);
		    if (v2 == sortcols[iplusk] /*&& block[sortcols[i]] == minBlock && block[sortcols[iplusk]] == minBlock*/) {
		      //cerr << "found " << v1 << " <-> " << v2 << " -- " << solution[v1].asDouble() << "," << solution[v2].asDouble() << endl;
		      //cerr << "Bilde: (1-x" << v1 << ") + ";
		      if (pCM->getAdjacent(kk) & 1) {
			//cerr << "x" << v2 << " >=1 " << endl;
			//cerr << "Clause schon vorhanden 1" << endl;
		      } else {
			//cerr << "(1-x" << v2 << ") >= 1" << endl;
			if (0&&info_level >= -5) cerr << "Fall 1b set x" << sortcols[i] << " auf 0" << endl;
			clist.push_back(make_pair(-sortcols[i]-1,0));
			foundFixing = true;
		      }
		    }
		    kk = pCM->NextAdjacentInConflictGraph(kk);
		  }
		}
		j = pCM->FirstAdjacentInConflictGraph(2*sortcols[i]+1);  // Fall x_i <= x_j
		if (j >= 0) {
		  int kk = pCM->NextAdjacentInConflictGraph(j);
		  while (kk >= 0) {
		    int v1 = sortcols[i];
		    int v2 = pCM->getAdjacent(kk);
		    v2 /= 2;
		    assert(types[v2] != 5000/*CONTINUOUS*/);
		    if (v2 == sortcols[iplusk]  /*&& block[sortcols[i]] == minBlock && block[sortcols[iplusk]] == minBlock*/) {
		      //cerr << "found " << v1 << " <-> " << v2 << " -- " << solution[v1].asDouble() << "," << solution[v2].asDouble() << endl;
		      //cerr << "Bilde: x" << v1 << " + ";
		      if (pCM->getAdjacent(kk) & 1) {
			//cerr << "x" << v2 << " >=1 " << endl;
			//cerr << "set x" << sortcols[iplusk] << " auf 1" << endl;
			clist.push_back(make_pair(-sortcols[iplusk]-1,1));
			foundFixing = true;
		      } else {
			//cerr << "(1-x" << v2 << ") >= 1" << endl;
			//cerr << "Fall2b set x" << sortcols[i] << " = x" << sortcols[iplusk] << endl;
			clist.push_back(make_pair(sortcols[i],-sortcols[iplusk]-1));
			foundFixing = true;
		      }
		    }
		    kk = pCM->NextAdjacentInConflictGraph(kk);
		  }
		}
		if (foundFixing==true) {
		  //clist.push_back(make_pair(sortcols[i],sortcols[iplusk]));
		} else {
		  foundFixing = false;
		  j = pCM->FirstAdjacentInConflictGraph(2*sortcols[iplusk]); // Fall x_i >= x_j
		  if (j >= 0) {
		    //cerr << "j=" << j << endl;
		    int kk = pCM->NextAdjacentInConflictGraph(j);
		    while (kk >= 0) {
		      //cerr << "j=" << j << " kk=" << kk << " a(j)=" << CM.getAdjacent(j) << " a(kk)=" << CM.getAdjacent(kk) << endl;
		      int v1 = sortcols[iplusk];
		      int v2 = pCM->getAdjacent(kk);
		      v2 /= 2;
		      assert(types[v2] != /*CONTINUOUS*/5000);
		      if (v2 == sortcols[i]  /*&& block[sortcols[iplusk]] == minBlock && block[sortcols[i]] == minBlock*/) {
			//cerr << "found " << v1 << " <-> " << v2 << " -- " << solution[v1].asDouble() << "," << solution[v2].asDouble() << endl;
			//cerr << "Bilde: (1-x" << v1 << ") + ";
			if (pCM->getAdjacent(kk) & 1) {
			  //cerr << "x" << v2 << " >=1 " << endl;
			  //cerr << "Fall3a set x" << sortcols[i] << " = x" << sortcols[iplusk] << endl;
			  clist.push_back(make_pair(sortcols[i],-sortcols[iplusk]-1));
			  foundFixing = true;
			} else { //f�hrt zu keiner neuen Erkenntnis
			  //cerr << "(1-x" << v2 << ") >= 1" << endl;
			  //cerr << "Fall3b set x" << sortcols[iplusk] << " auf 0" << endl;
			  //cerr << "weil: SYMMETRIE FOUND!! i=" << i << ",k=" << k << " | x" << sortcols[i] << " <= x" << sortcols[iplusk] << endl;
			  foundFixing = true;
			  clist.push_back(make_pair(-sortcols[iplusk]-1,0));
			}
		      }
		      kk = pCM->NextAdjacentInConflictGraph(kk);
		    }
		  }
		  j = pCM->FirstAdjacentInConflictGraph(2*sortcols[iplusk]+1);// Fall x_i >= x_j
		  if (j >= 0) {
		    int kk = pCM->NextAdjacentInConflictGraph(j);
		    while (kk >= 0) {
		      int v1 = sortcols[iplusk];
		      int v2 = pCM->getAdjacent(kk);
		      v2 /= 2;
		      assert(types[v2] != 5000/*CONTINUOUS*/);
		      if (v2 == sortcols[i]  /*&& block[sortcols[iplusk]] == minBlock && block[sortcols[i]] == minBlock*/) {
			//cerr << "found " << v1 << " <-> " << v2 << " -- " << solution[v1].asDouble() << "," << solution[v2].asDouble() << endl;
			//cerr << "Bilde: x" << v1 << " + ";
			if (pCM->getAdjacent(kk) & 1) {
			  //cerr << "x" << v2 << " >=1 " << endl;
			  //cerr << "set x" << sortcols[i] << " auf 1" << endl;
			  clist.push_back(make_pair(-sortcols[i]-1,1));
			  foundFixing = true;
			} else {
			  //cerr << "(1-x" << v2 << ") >= 1" << endl;
			  //cerr << "Clause schon vorhanden 4" << endl;
			}
		      }
		      kk = pCM->NextAdjacentInConflictGraph(kk);
		    }
		  }
		  if (!foundFixing)
		    clist.push_back(make_pair(sortcols[i],sortcols[iplusk]));
		}
	      }
	      int sokoI=-10;
	      int sokoIK=-10;
	      if (xdy==-1) {
		//cerr << "SYMMERR" << endl;
		//if (info_level > 5) cerr << "DOMINANCE TYPE 1!! " << i << "," << iplusk << " | " << sortcols[i] << ", " << sortcols[iplusk] << endl;
		//clist.push_back(make_pair(sortcols[iplusk],sortcols[i]));
		sokoI = sortcols[iplusk];
		sokoIK = sortcols[i];
	      }
	      if (xdy==1) {
		sokoI = sortcols[i];
		sokoIK = sortcols[iplusk];
	      }
	      if (xdy==1 || xdy==-1) {
		latestHitI = i;
		latestHitIplusK = iplusk;
		//if (outputs < 5) cerr << "info DOMINANCE TYPE " << k << "!! " << i << "," << iplusk << " | " << sokoI << ", " << sokoIK << endl;
		//else if (outputs==5) cerr << " ... " << endl;
		outputs++;
		//clist.push_back(make_pair(sokoI,sokoIK));
		int j;
		bool foundFixing=false;
		j = pCM->FirstAdjacentInConflictGraph(2*sokoI); // Fall x_i <= x_j
		if (j >= 0) {
		  //cerr << "j=" << j << endl;
		  int kk = pCM->NextAdjacentInConflictGraph(j);
		  while (kk >= 0) {
		    //cerr << "j=" << j << " kk=" << kk << " a(j)=" << CM.getAdjacent(j) << " a(kk)=" << CM.getAdjacent(kk) << endl;
		    int v1 = sokoI;
		    int v2 = pCM->getAdjacent(kk);
		    v2 /= 2;
		    assert(types[v2] != /*CONTINUOUS*/5000);
		    if (v2 == sokoIK  /*&& block[sokoI] == minBlock && block[sokoIK] == minBlock*/) {
		      //cerr << "found " << v1 << " <-> " << v2 << " -- " << solution[v1].asDouble() << "," << solution[v2].asDouble() << endl;
		      //cerr << "Bilde: (1-x" << v1 << ") + ";
		      if (pCM->getAdjacent(kk) & 1) {
			//cerr << "x" << v2 << " >=1 " << endl;
			//cerr << "Clause schon vorhanden 1" << endl;
		      } else {
			//cerr << "(1-x" << v2 << ") >= 1" << endl;
			//cerr << "Fall 5b set x" << sokoI << " auf 0" << endl;
			//cerr << "weil: DOMINANCE FOUND!! i=" << i << ",k=" << k << " | x" << sokoI << " <= x" << sokoIK << endl;
			clist.push_back(make_pair(-sokoI-1,0));
			foundFixing = true;
		      }
		    }
		    kk = pCM->NextAdjacentInConflictGraph(kk);
		  }
		}
		j = pCM->FirstAdjacentInConflictGraph(2*sokoI+1);  // Fall x_i <= x_j
		if (j >= 0) {
		  int kk = pCM->NextAdjacentInConflictGraph(j);
		  while (kk >= 0) {
		    int v1 = sokoI;
		    int v2 = pCM->getAdjacent(kk);
		    v2 /= 2;
		    assert(types[v2] != 5000/*CONTINUOUS*/);
		    if (v2 == sokoIK /*&& block[sokoI] == minBlock && block[sokoIK] == minBlock*/) {
		      //cerr << "found " << v1 << " <-> " << v2 << " -- " << solution[v1].asDouble() << "," << solution[v2].asDouble() << endl;
		      //cerr << "Bilde: x" << v1 << " + ";
		      if (pCM->getAdjacent(kk) & 1) {
			//cerr << "x" << v2 << " >=1 " << endl;
			//cerr << "set x" << sokoIK << " auf 1" << endl;
			clist.push_back(make_pair(-sokoIK-1,1));
			foundFixing = true;
		      } else {
			//cerr << "(1-x" << v2 << ") >= 1" << endl;
			//cerr << "set x" << sokoI << " = x" << sokoIK << endl;
			clist.push_back(make_pair(sokoI,-sokoIK-1));
			foundFixing = true;
		      }
		    }
		    kk = pCM->NextAdjacentInConflictGraph(kk);
		  }
		}
		static int cntcl=0;
		if (foundFixing == false /*&& block[sokoI] == minBlock && block[sokoIK] == minBlock*/) {
		  clist.push_back(make_pair(sokoI,sokoIK));
		  /*
		    std::vector<data::IndexedElement > *X = QlpStSolve.getExternSolver(maxLPstage).getCol_snapshot(sortcols[i]);
		    std::vector<data::IndexedElement > *Y = QlpStSolve.getExternSolver(maxLPstage).getCol_snapshot(sortcols[iplusk]);
		    cout << "?x" << sortcols[i]+1 << " <= x" << sortcols[iplusk]+1 << endl;
		    cerr << "x" << sortcols[i] << " <= x" << sortcols[iplusk] << "; " << i << "," << iplusk << endl;
		    cerr << "Rows with xi:";
		    for (int i = 0; i < (*X).size();i++)
		    cerr << (*X)[i].value.asDouble() << "u" << (*X)[i].index << " ";
		    cerr << endl;
		    cerr << "Rows with xiplusk:";
		    for (int i = 0; i < (*Y).size();i++)
		    cerr << (*Y)[i].value.asDouble() << "u" << (*Y)[i].index << " ";
		    cerr << endl;
		  */
		  cntcl++;
		}
	      }
	    }
	  }
	  if (info_level >= 2) cerr << "Minblock=" << minBlock << endl;
	}

	double getLPoffset() {
		return LPoffset;
	}
	bool getIsInSOSvars(int x_i) {
		if (SOSvars.count(std::pair<int,int>(x_i,0)) > 0) return true;
		return false;
	}
	bool adaptSolution(std::vector<data::QpNum>& solution, int *types, extbool* assigns) {
		std::vector<int> SOSvarsSparse;
		for (int z = 0; z < solution.size() && z < getnVars();z++) {
			if (types[z]==0 /*BINARY*/) {
				if (SOSvars.count(std::pair<int,int>(z,0)) > 0) {
					if (assigns[z]!=2 && fabs((double)assigns[z]-solution[z].asDouble()) > 0.5) cerr << "SOSvar: " << (int)assigns[z] << "," << solution[z].asDouble() << endl;
					//std::cerr << "x" << x_i << " already in SOSvars." << std::endl;
					solution[z] = -10.0;
					SOSvarsSparse.push_back(z);
				} else if (assigns[z]!=2) {
					if (fabs((double)assigns[z]-solution[z].asDouble()) > 0.1) {
						cerr << "Error: assigned but fractional: " << (int)assigns[z] << "," << solution[z].asDouble() << endl;
						return false;
					}
					solution[z] = (double)assigns[z];
				}
			}
		}
		return true;

		bool uneval_exists=false;
		int round=0;
		int success = 0;
		do {
			round++;
			uneval_exists=false;
			for (int i = 0; i < SOSvarsSparse.size();i++) {
				assert(types[SOSvarsSparse[i]]==0 /*BINARY*/);
				if (solution[SOSvarsSparse[i]] < -0.5) {
					//suche zeile L in SOSconstraints

				    // binary search
					int v = SOSvarsSparse[i];
					if (round > solution.size()) cerr << "einzusetzende Variable ist x" << v << endl;
					int jj=SOSconstraints.size()-1,ii=0;
					int match=-1;
					while (jj >= ii) {
						int m = (ii + jj)>>1;
						//if (v==24) cerr << ii <<"," << jj <<  "," << m << "," << SOSconstraints[m].second<< endl;
						int w = fabs(SOSconstraints[m].second);
						if (w==v && SOSconstraints[m].second >= 0) {
							match = m;
							break;
						} else {
							if (ii==jj) break;
							if (v>w) {
								ii=m+1;
							} else if (v<w) {
								jj=m-1;
							}
						}
					}
					assert(match != -1);
					if (match!=-1) {
						//versuch x_SOSvarsSparse[i] aus L herauszurechnen.
						std::vector< std::pair<int,double> > &sosLHS = SOSconstraints[match].first;
						double result = 0.0;
						if (round > solution.size()) cerr << "Term und Werte: ";
						for (int ii = 0; ii < sosLHS.size();ii++) {
							if (sosLHS[ii].first > solution.size()) {
								if (round > solution.size()) cerr << sosLHS[ii].second << "x" << sosLHS[ii].first << "(.)" << " ";
								result = result + sosLHS[ii].second;
							} else if (solution[sosLHS[ii].first].asDouble() >= -1 ) {
								if (round > solution.size()) cerr << sosLHS[ii].second << "x" << sosLHS[ii].first << "(" << solution[sosLHS[ii].first].asDouble() << ")" << " ";
								result = result + solution[sosLHS[ii].first].asDouble()*sosLHS[ii].second;
							} else {
								uneval_exists = true;
								break;
							}
						}
						//wenn das geht: in solution eintragen
						//sonst: uneval_exists=true und next variable
						if (uneval_exists == false) {
							success++;
							if (-result < -1e-9 || -result > 1+1e-9) cerr << "RANGE: " << result << endl;
							assert(-result >= -1e-9 && -result <= 1+1e-9);
							solution[v] = -result;
							if (round > solution.size()) cerr << "MATCH!! solution[]=" << solution[match].asDouble() << endl;
						}
					}
				}
			}
			if (uneval_exists && round > solution.size()) {
				cerr << "Round " << round << " und success=" << success << endl;
				cerr << "#sos vars: " << SOSvarsSparse.size() << endl;


				char a;
				cin >> a;
			}
		} while (uneval_exists);
		//for (int ii=0;ii < solution.size();ii++) {
		//	cerr << "x" << ii << "=" << solution[ii].asDouble() << ", ";
		//}
		//cerr << endl;
		return true;
	}

	int64_t ggt(int64_t a, int64_t b);
	bool exactAvail(std::vector<data::IndexedElement> &table_lhs, std::vector<data::IndexedElement> &lhs, data::QpRhs table_rhs, data::QpRhs rhs);
	void simplifyCoefficients_01(std::vector<data::IndexedElement> &lhs, std::vector<data::IndexedElement> &ilhs, data::QpRhs &rhs, data::QpRhs &org_rhs, int8_t *ba);
	void analyzeAndChangeIntRow(std::vector<data::IndexedElement> &lhs, std::vector<data::IndexedElement> &ilhs, data::QpRhs &rhs, data::QpRhs &org_rhs, int8_t *ba);
	bool extractSOS(std::vector<data::IndexedElement> &org_lhs, data::QpRhs &org_rhs, int c) {return false;}
	void addNegativeOfSosC2ToSosC1(std::vector< std::pair<int,double> > &c1, std::vector< std::pair<int,double> > &c2, int x_i);
	bool checkStabilityOfMult(std::vector< std::pair<int,double> > &c1, std::vector< std::pair<int,double> > &c2, int x_i) ;
	void replaceSOSvars( int );
	void replaceSOSvars( std::vector< std::vector<data::IndexedElement> > &LHSs, std::vector<data::QpRhs> &RHSs, int c );
	void replaceSOSvarsInObj( std::vector<data::IndexedElement> &OBJ, data::QpRhs &OBJBND , int cntVars);
	bool preprocessConstraint(std::vector<int> &v_ids, std::vector<data::IndexedElement> &LHS_chg, data::QpRhs &RHS_chg,
			data::Qlp &qmip, utils::QlpStageSolver& QlpStSolve, HCTable * hct, std::pair<coef_t,uint64_t> &hp,
			int8_t *binary_assignments_, coef_t rhsbnd, int *types,
			int maxLPstage, std::vector<std::pair<int,double> > &cpropQ, coef_t *lowerBounds, coef_t *upperBounds,
			bool feasPhase, std::vector< std::pair<int,int> > &clist, int *block, int *eas);
#ifndef FIND_BUG
	void updateConstraints(data::Qlp &qmip, utils::QlpStageSolver& QlpStSolve, int8_t *ba, coef_t rhsbnd, int *, int,
			std::vector<std::pair<int,double> > &cpropQ, coef_t *lowerBounds, coef_t *upperBounds,
			bool feasPhase, std::vector< std::pair<int,int> > &clist, int *block, int *eas);
#else
	void updateConstraints(data::Qlp &qmip, utils::QlpStageSolver& QlpStSolve, int8_t *ba, coef_t rhsbnd, int *, int);
#endif
	void updateContBounds(data::Qlp &qmip, utils::QlpStageSolver& QlpStSolve, int8_t *binary_assignments_,
			int *types, int maxLPstage, std::vector< std::vector<data::IndexedElement> > &LHSs, std::vector<data::QpRhs> &RHSs , int cntVars,
			coef_t *lowerBounds, coef_t *upperBounds, std::vector<std::pair<int,double> > &cpropQ);
	void findComponents(data::Qlp &qmip, int8_t *assigns, int *components, std::vector< std::vector<int> > &varsOfComponents);
	void strongConnect(CliqueManager &CM, std::vector<node> &V, int8_t *assigns, int *components, int *index, node v, int *c_identifier, int *types);
	void findStrongComponents(CliqueManager &CM, int8_t *assigns, int *components, int *c_identifier, int *types);
	int makeBinary(data::Qlp &qlp, data::Qlp &qbp);

	algorithm::Algorithm::QlpSolution solveDEP(bool uBin = false) {
		data::Qlp qlp2;
                int ifl = getInfoLevel();
		std::cerr << "solveDEP InfoLevel is = " << getInfoLevel() << std::endl;
                if (ifl < 0) ifl = -ifl;
		if (ifl % 2 == 0 || uBin) {
		  std::cerr << "USE BINARIZATION." << endl;
		  integers.clear();
		  integers_lim.clear();
		  offset_shift = 0.0;
		  int SimplificationFailed=makeBinary(/*orgQlp, orgQlpBinarized*/qlp,qlp2);
		  algorithm::Qlp2Lp DEP(qlp2);
		  return DEP.solveQlp(algorithm::Algorithm::WORST_CASE);
                } else {
		  std::cerr << "NO BINAIZATION." << endl;
		  algorithm::Qlp2Lp DEP(qlp);
		  return DEP.solveQlp(algorithm::Algorithm::WORST_CASE);
		}

    }
    algorithm::Algorithm::QlpSolution resolveDEP(int8_t *assigns) {
		data::Qlp qlp2;
		int cnt=0;
                int ifl = getInfoLevel();
		std::cerr << "resolve DEP InfoLevel is = " << getInfoLevel() << std::endl;

                if (ifl < 0) ifl = -ifl;
		if (ifl % 2 == 0) {
		  std::cerr << "USE BINARIZATION." << endl;
		  integers.clear();
		  integers_lim.clear();
		  offset_shift = 0.0;
		  int SimplificationFailed=makeBinary(/*orgQlp, orgQlpBinarized*/qlp,qlp2);
                } else {
		  std::cerr << "NO BINARIZATION." << endl;
		  qlp = qlp2;
		}

		for (int z=0; z < qbp->nVars();z++)
				if (assigns[z] != 2) {
						   cnt++;
						   qlp.getVariableVector().at(z)->setLowerBound((double)assigns[z]);
						   qlp.getVariableVector().at(z)->setUpperBound((double)assigns[z]);
				} else if (qbp->isFixed(z)) {
						   qlp.getVariableVector().at(z)->setLowerBound((double)qbp->getFixed(z));
						   qlp.getVariableVector().at(z)->setUpperBound((double)qbp->getFixed(z));
				} else {
						   qlp.getVariableVector().at(z)->setLowerBound(0.0);
						   qlp.getVariableVector().at(z)->setUpperBound(1.0);
				}
		cerr << "es wurden " << cnt << " Varablen gestzt" << endl;
		algorithm::Qlp2Lp DEP(qlp2);
		return DEP.solveQlp(algorithm::Algorithm::WORST_CASE);
    }
    void dummyDEP() {
		algorithm::Qlp2Lp DEP(qlp);
        
        data::QpObjFunc obj;
        std::vector<data::QpVar> vars;
        data::QpSparseMatrix matrix;
        std::vector<data::QpRhs> rhs;
        
        utils::QlpConverter::convertToLP(qlp, obj, vars, matrix, rhs, utils::QlpConverter::COMPACT_VIEW, utils::QlpConverter::WORST);
    }
	void saveDepSolution(std::vector<data::QpNum>  &solu) {
		qbp->saveDepSolution(solu);
	}

	void setInfoLevel(int l) { info_level = l; }
	int getInfoLevel() { return info_level; }

	void setTimeout(time_t t) {
		timeout = t;
		qbp->setTimeout(t);
	}

    bool isOnTrack(ca_vec<int> &optSol, std::vector<data::IndexedElement>  &LHS, data::QpRhs &RHS) {
    	//if (USE_TRACKON == 0) return false;
    	return false;
    	cerr << ".";
    	double sum=0.0;
    	for (int i = 0;i < LHS.size();i++) {
    		int v = LHS[i].index;
    		double d = LHS[i].value.asDouble();
    		//if (d < -1e10 || d > 1e10 || optSol[v] > 2 || optSol[v] < 0) return true;
    		sum = sum + d * optSol[v];
    		if (d!=d) {
    			cerr << "d is NAN! Index=" << v << endl;
    			return true;
    		}
    		if (sum!=sum) {
    			cerr << "sum is NAN! Index=" << v << endl;
    			return true;
    		}
    	}
    	if (RHS.getRatioSign() == data::QpRhs::RatioSign::equal && (sum > RHS.getValue().asDouble()+LP_EPS || sum < RHS.getValue().asDouble()-LP_EPS)) {
    		cerr << "is " << sum << " == " << RHS.getValue().asDouble() << " ?" << endl;
    		return true;
    	} else if (RHS.getRatioSign() == data::QpRhs::RatioSign::smallerThanOrEqual && sum > RHS.getValue().asDouble()) {
    		cerr << "is " << sum << " <= " << RHS.getValue().asDouble() << " ?" << endl;
    		return true;
    	} else if (RHS.getRatioSign() == data::QpRhs::RatioSign::greaterThanOrEqual && sum < RHS.getValue().asDouble()){
    		cerr << "is " << sum << " >= " << RHS.getValue().asDouble() << " ?" << endl;
    		return true;
    	}
    	return false;
    }

    void yInit(data::Qlp& orgQlp,std::string inputfile="") {

      qlp = orgQlp;
      qlpRelax = orgQlp;

      utils::QlpConverter::relaxQlpNumberSystem(qlpRelax);

      std::vector<data::QpVar*> vars = qlpRelax.getVariableVectorByQuantifier(data::QpVar::all);
      std::vector<int> l_cU;
      if (vars.size() > maxUnivVars) {
	for(unsigned int i = 0; i < vars.size() - maxUnivVars;i++){
	  vars[i]->setQuantifier(data::QpVar::exists);
	  l_cU.push_back(vars[i]->getIndex());
	}
      }

      //QlpStSolve = new utils::QlpStageSolver(qlpRelax,true,false);
      qbp = new QBPSolver(qlpRelax);
#ifndef NO_CGL
      cbc = new CBCSolver();
#endif
      qbp->setyIF(this);
      qbp->determineFixedUniversalVars(l_cU); // not used, gives the possiblity to restrict the number of universal variables and leads to a relaxation
      yReadInput(yParamMaxHashSize);
      //cerr << "Start Ini AllSolver" << endl;
      qbp->AllSolver->init(orgQlp, data::QpRhs::Responsibility::UNIVERSAL);
      qbp->ExistSolver->init(orgQlp, data::QpRhs::Responsibility::EXISTENTIAL);
       if(qbp->getShowInfo()) cerr << "info: initialized AllSolver" << endl;
      if (qbp->nVars() > 40000 || qbp->getHasObjective() == false) {
	 if(qbp->getShowInfo()) std::cerr << "info: auto-deactivation of Cgl-Options" << std::endl;
	qbp->setUseCglRootCuts(false);
      }
#ifndef NO_CGL
      cbc->CutGenSolver->init(orgQlp, data::QpRhs::Responsibility::EXISTENTIAL);
      if (qbp->getInfoLevel()>-8 && qbp->getUseCglRootCuts()) {
	if(qbp->getShowInfo())cerr << "info: Initialized CBCSolver" << endl;
      }
      if (qbp->getInfoLevel()>-8) cerr << "Row Count is " << cbc->CutGenSolver->getRowCount() <<endl;
  //    cbc->CutGenSolver->solve();
//      cerr << "Solution at root: " << cbc->CutGenSolver->getObjValue() << endl;
//      cbc->CutGenSolver->CreateCuts();
#endif
      qbp->DepManagerInitGraph();

      if(qbp->getShowInfo()) cerr << "info: begin scan dependencies" << endl;
      qbp->DepManScanConstraints();
      if(qbp->getShowInfo()) cerr << "info: end scan dependencies" << endl;
      for (int z=0; z < qbp->nVars();z++) {
	qbp->insertVar2PriorityQueue(z);
      }
      qbp->DepManInitFillRate();
      qbp->setInputFileName(inputfile);
    }

	int yInit(std::string inputfile){
		data::Qlp orgQlp;
		data::Qlp orgQlpBinarized;
		utils::Parser::createQlp(inputfile, orgQlp/*orgQlpBinarized*/);
		if( orgQlp.getObjective() == data::QpObjFunc::max ){
		  orgQlp.reverseObjFunc();
		  objInverted = true;
		} else objInverted = false;
		int SimplificationFailed=makeBinary(orgQlp, orgQlpBinarized);
		if(SimplificationFailed) return SimplificationFailed;
		yInit( orgQlpBinarized,inputfile);
		this->qbp->setFinalOffset( orgQlp.getObjFuncOffset().asDouble());
		return 0;
	}

	bool nodeValueKnown(coef_t &r) {
		return false;
	}
	void preRec() {
		alive.push_back(true);
	}
	bool levelAlive() {
		bool x = alive[alive.size()-1];
		return x;
	}
	void postRec() {
		alive.pop_back();
	}
	bool deepImplicationAvailable() {
		return qbp->deepImplicationAvailable();
	}
	bool selectVariable(int &v) {
		for (int i = 0; i < getnVars();i++)
			if (qbp->getCurrentVariableAssignment(i) == extbool_Undef) {
				v = i;
				return true;
			}
		return false;
	}
	int choosePolarity() {
		return 0;
	}
	bool assignVariable(int pol, int var) {
		int oob = qbp->assign(var, pol, qbp->getTrailSize(),CRef_Undef, false);
		if (oob == ASSIGN_OK) return true;
		return false;
	}
	void killLevelsUntil(int &level, ValueConstraintPair &out_vcp) {
		QBPSolver::ValueConstraintPair VCP;
		//cerr << "go back to level " << level << " and imply " << out_vcp.v/2 << " to " << 1-(out_vcp.v&1) << endl;
		VCP.v = out_vcp.v;
		VCP.cr = out_vcp.cr;
		VCP.pos = out_vcp.pos;
		for (int i = /*qbp->decisionLevel();*/alive.size()-1; i > level-2;i--) {
			alive[i] = false;
			//cerr << "set dead" << i << endl;
		}
		qbp->addImplication(VCP);
	}
	bool propagate() {
		//if (qbp->propQ.size()==1) cerr << "eas=" << qbp->eas[qbp->propQ[0].v/2] << endl;
		bool r = qbp->propagate(confl, confl_var, confl_partner, false, false, 1000000);
		return r;
	}
	void write_nodeinfo(int nodeID){
	  qbp->write_nodeinfo(nodeID);
	}
	void write_successors(int nodeID) {
	  qbp->write_successors(nodeID);
	}

	bool analyzeImplicationGraph(ValueConstraintPair& out_info, int &target_level) {
		QBPSolver::ValueConstraintPair out_vcp;
		int out_target_dec_level;
		ca_vec<CoeVar>   out_learnt;

		if (qbp->getQuantifier(confl_var) == EXIST) {
			//out_vcp.v = -1;
			out_vcp.pos = -1;
			//out_vcp.cr = -1;
			//cerr << "TL=" << out_target_dec_level;
			if (qbp->analyze(confl, confl_var, confl_partner, out_learnt, out_target_dec_level, out_vcp) && out_vcp.pos != -1) {
				target_level = out_target_dec_level;
				out_info.cr = out_vcp.cr;
				out_info.pos = out_vcp.pos;
				out_info.v = out_vcp.v;
				//if (qbp->propQlimiter[out_vcp.v] <= 0) {
				//   qbp->propQ.push(out_vcp);
				//   qbp->propQlimiter[out_vcp.v] = qbp->propQ.size();
				//} else qbp->propQ[qbp->propQlimiter[out_vcp.v]-1] = out_vcp;
				return true;
			} else return false;
		}  else {
			//out_vcp.v = -1;
			out_vcp.pos = -1;
			//out_vcp.cr = -1;
			if (qbp->analyze4All(confl, confl_var, out_learnt, out_target_dec_level, out_vcp) && out_vcp.pos != -1) {
				target_level = out_target_dec_level;
				out_info.cr = out_vcp.cr;
				out_info.pos = out_vcp.pos;
				out_info.v = out_vcp.v;
				//if (qbp->propQlimiter[out_vcp.v] <= 0) {
				//   qbp->propQ.push(out_vcp);
				//   qbp->propQlimiter[out_vcp.v] = qbp->propQ.size();
				//} else qbp->propQ[qbp->propQlimiter[out_vcp.v]-1] = out_vcp;
				return true;
			} else return false;
		}
	}

	void unassignVariablesUntil(int v) {
		int wv;
		if (qbp->getTrailSize()>0) wv = qbp->getTrailsLast();
		else wv = -1;
		//cerr << "Trail (" << v << "," << qbp->decisionLevel()-1 << "):";
		for (int z = 0; z < qbp->getTrailSize();z++) {
			//cerr << qbp->trail[z] << (qbp->vardata[qbp->trail[z]].reason!=CRef_Undef ?"i-" :"-") << qbp->vardata[qbp->trail[z]].level << " ";
		}
		//cerr << endl;
		while (wv > -1 && wv != v) {
			qbp->unassign(wv, false);
			if (qbp->getTrailSize()>0) wv = qbp->getTrailsLast();
			else wv = -1;
		}
		assert(wv == v);
		qbp->unassign(wv, false);
	}
	coef_t yMySolve() {
		coef_t score;
		coef_t value;
		coef_t r;
		int v;
		ValueConstraintPair VCP;
		int target_level;

		if (nodeValueKnown(r)) return r;
		if (!selectVariable(v)) return p_infinity;
		if (qbp->getQuantifier(v)==EXIST) score = n_infinity;
		else                              score = p_infinity;
//cerr << "select " << v << " in " << qbp->decisionLevel() << endl;
		int p = choosePolarity();
		for (int i = 0; i <= 1;i++) {
			if (!levelAlive()) {
				//cerr << "leave dead level " << qbp->decisionLevel() << endl;
				//for (int o=0;o<alive.size();o++) cerr << "o=" << o << "A=" << alive[o];
				//cerr << endl;
				return score;
			}
			qbp->increaseDecisionLevel();
			bool success = assignVariable(!p ? i : 1-i ,v);
			//cerr << "assign " << v << "=" << (!p ? i : 1-i) << endl;
 			if (!success) continue;
			do {
				//if (qbp->propQ.size()>0) cerr << "prop1" << endl;
				success = propagate();
				//if (qbp->propQ.size()>0) cerr << "prop2" << endl;
				//else cerr << "S=" << success << endl;
				//cerr << "Trail after P(" << v << "," << qbp->decisionLevel()-1 << "):";
				for (int z = 0; z < qbp->getTrailSize();z++) {
					//cerr << qbp->trail[z] << (qbp->vardata[qbp->trail[z]].reason!=CRef_Undef ?"i-" :"-") << qbp->vardata[qbp->trail[z]].level << " ";
				}
				//cerr << endl;

				if (!success) {
					if (analyzeImplicationGraph(VCP, target_level)) {
						killLevelsUntil(target_level, VCP);
						//(v);
						//qbp->decreaseDecisionLevel();
						//return score;
						break;
					}
				} else {
					preRec();
				    value = yMySolve();
				    postRec();
					if (qbp->getQuantifier(v) == EXIST) score = max(score, value);
					if (qbp->getQuantifier(v) == UNIV ) score = min(score, value);
				}
			} while (levelAlive() && deepImplicationAvailable());
			//cerr << "unassign " << v << endl;
			unassignVariablesUntil(v);
			qbp->decreaseDecisionLevel();
	  	}
		//cerr << "leave level " << qbp->decisionLevel() << " with var " << v << endl;
		//cerr << "Trail3 (" << v << "," << qbp->decisionLevel()-1 << "):";
		for (int z = 0; z < qbp->getTrailSize();z++) {
			//cerr << qbp->trail[z] << (qbp->vardata[qbp->trail[z]].reason!=CRef_Undef ?"i-" :"-") << qbp->vardata[qbp->trail[z]].level << " ";
		}
		//cerr << "---------->" << score << endl;
		if (!levelAlive()) {
			//cerr << "leave deadDD level " << qbp->decisionLevel() << endl;
			return score;
		}
		return score;
	}

	inline bool isZero(double x, double epsZero = 1e-20) {
		return (x >= -epsZero && fabs(x) <= epsZero);
	}
	inline bool isOne(double x, double epsZero = 1e-20) {
		return (x >= 1.0-epsZero && fabs(x) <= 1.0+epsZero);
	}

	coef_t ySolve(time_t ini_time, int z) {
#ifdef NEWSOLVE
		preRec();
	    qbp->increaseDecisionLevel();
	    if (yMySolve() > 0) cerr << "has solution" << endl;
	    else cerr << "has no solution" << endl;
	    postRec();
	    qbp->decreaseDecisionLevel();
	    return 0.0;
#endif
		coef_t r;
		HTable *HT = qbp->getHTable();
		std::vector<int> pt2leaders;
		std::vector<int> bitcnts;
		qbp->setInitilizationTime(ini_time);
	    qbp->increaseDecisionLevel();
	    r = qbp->search(0, (void*)this);
		for (int i = 0; i < qbp->nVars(); i++) {
			bitcnts.push_back(integers[i].bitcnt);
			pt2leaders.push_back(integers[i].pt2leader);
		}
	    if ((processNo & 1) == 0) {
	    	if (fabs(qbp->getGlobalScore()) < -qbp->getNegativeInfinity()) {
				double min_slack = qbp->finalSolutionCheck(bitcnts,pt2leaders,objInverted);
				for (int i = 0; i < qbp->nVars(); i++) {
					//cerr << "x" << i << "(" << qbp->getType(i) << "," << qbp->getBlock(i) << "," << integers[i].bitcnt << ") ";
					if (qbp->getType(i) == 5000 /*CONTINUOUS*/) {
						if (qbp->getLowerBound(i) == qbp->getUpperBound(i)) {
							cerr << qbp->getLowerBound(i) << " ";
						} else if (qbp->getBlock(i) == 1) cerr << qbp->getFirstStageSolutionValue(i) << " ";
							   else cerr << "c" <<  qbp->getBlock(i) << " ";
					} else if (qbp->getBlock(i) == 1){
						if (integers[i].bitcnt > 1) {
							//cerr << "Var " << i<< " hat " << integers[i].bitcnt << "bits. Leader:" << integers[i].pt2leader;
							if (integers[i].pt2leader == i) {
								int res = 0;
								for (int z = 0;z < integers[i].bitcnt;z++) {
									res = 2*res + qbp->getFirstStageSolutionValue(i+z);
									//cerr << " | " << qbp->getFirstStageSolutionValue(i+z) << " | ";
								}
								cerr << res << " ";
								i+=integers[i].bitcnt;
								i--;
							} else cerr << "Error: missing bit" << endl;
						} else
							cerr << qbp->getFirstStageSolutionValue(i) << " ";
					} else {
						cerr << "s" << qbp->getBlock(i) << " ";
					}
				}
				cerr << endl;
	    	}
	    }
	    while (qbp->getVarPriorityQueueSize() > 0) {
	    	qbp->extractVarPriorityQueueMinimum();
	    }
		HT->~HTable();
	    if (info_level >= 0) return r;
	    else return (r > (coef_t)0 ? (coef_t)1 : (coef_t)0);
	}


	coef_t ySolve(time_t ini_time) {
    int pass = 1;
    int search_mode = /*NORMAL_MODE;*/RELAXATION_MODE;//RELAXATION_MODE;//RESTRICTION_MODE;//NORMAL_MODE;
    if(qbp->getUniversalConstraintsExist()){
	search_mode = NORMAL_MODE;
    }
    std::vector<data::IndexedElement> restrictlhs; 
    double restrictrhs;
    int widthFactor = 1;
#ifdef NEWSOLVE
    preRec();
    qbp->increaseDecisionLevel();
    if (yMySolve() > 0) cerr << "has solution" << endl;
    else cerr << "has no solution" << endl;
    postRec();
    qbp->decreaseDecisionLevel();
    return 0.0;
#endif
    coef_t r;
    std::vector<double> fstStSol;
    double global_score=qbp->getNegativeInfinity();
    double global_dual_bound=-qbp->getDontKnowValue();
    double alpha = global_score;// - fabs(global_score)*0.0001 - 0.00001;
    double beta = qbp->getDontKnowValue() / 2.0;//-qbp->getDontKnowValue();
    fstStSol.clear();
    qbp->setPhase(true);
    restrictlhs.clear();
    restrictrhs = 0.0;
    if(qbp->getShowInfo()){
    cerr << "learnDualCuts=" << qbp->getLearnDualCuts() << endl;
    cerr << "useGMI=" << qbp->getUseGMI() << endl;
    cerr << "showInfo=" << qbp->getShowInfo() << endl;
    cerr << "showWarning=" << qbp->getShowWarning() << endl;
    cerr << "showError=" << qbp->getShowError() << endl;
    cerr << "useCover=" << qbp->getUseCover() << endl;
    cerr << "usePump=" << qbp->getUsePump() << endl;
    cerr << "useMiniSeach=" << qbp->getUseMiniSearch() << endl;
    cerr << "useShadow=" << qbp->getUseShadow() << endl;
    cerr << "useLazyLP=" << qbp->getUseLazyLP() << endl;
    cerr << "useCglRootCuts=" << qbp->getUseCglRootCuts() << endl;
    cerr << "useCglRootPreprocess=" << qbp->getUseCglRootPreprocess() << endl;
    cerr << "reduceStrongBranching=" << qbp->getReduceStrongBranching() << endl;
    cerr << "useMcts=" << qbp->getUseMcts() << endl;
    cerr << "useFarJumping=" << qbp->getUseBendersBackJump() << endl;
    cerr << "useFastFix=" << qbp->getUseFastFix() << endl;
    cerr << "useImplications=" << qbp->getUseImplications() << endl;
    cerr << "useLimitedLP=" << qbp->getUseLimitedLP() << endl;
    cerr << "useStrongBranching=" << qbp->getUseStrongBranching() << endl;
    cerr << "useEarlyBackjump=" << qbp->getUseEarlyBackjump() << endl;
    cerr << "useBestLevelExtraction=" << qbp->getUseBestLevelExtraction() << endl;
    cerr << "useUniversalBackjump=" << qbp->getUseUniversalBackjump() << endl;
    cerr << "maxBaCLevel=" << qbp->getMaxBaCLevel() << endl;
    cerr << "useAlphabeta=" << qbp->getUseAlphabeta() << endl;
    cerr << "useScout=" << qbp->getUseScout() << endl;
    cerr << "useMonotones=" << qbp->getUseMonotones() << endl;
    cerr << "isSimplyRestricted=" << qbp->getIsSimplyRestricted() << endl;
    cerr << "writeOutputFile=" << qbp->getWriteOutputFile() << endl;
    cerr << "maintainPv=" << qbp->getMaintainPv() << endl;
    cerr << "gzipExists=" << gzipExists << endl;
    cerr << "qdimacs2qipConverter: " << qdimacs2QipPath << endl;
    }
    qbp->setInitilizationTime(ini_time);
   
    { 
      int maxLPStage = qbp->getMaxLPStage();
      utils::QlpStageSolver *QlpStSolve = qbp->getStageSolverPt();
      int numConstraints = qbp->getNumberOfConstraints();//(*QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot()).size();
      if (numConstraints > 5) {
	search_mode = NORMAL_MODE;
        if(qbp->getShowInfo()){
	  cerr << "info switch to normal mode" << endl;
        }
        cerr << "Yasol ready"<<endl;
      }
    }
    double res = qbp->searchInitialization(0, (void*)this);
    if (res < 0) return getNinf();
    do {
      HTable *HT = qbp->getHTable();
      qbp->setInitilizationTime(ini_time);
      qbp->increaseDecisionLevel();
      assert(qbp->decisionLevel() == 1);

      if (search_mode == NORMAL_MODE) {
	beta = -qbp->getDontKnowValue();
	alpha = qbp->getNegativeInfinity();
	if (qbp->getInfoLevel() > -8) std::cerr << "in NORMAL_MODE alpha=" << alpha << " beta=" << beta << std::endl;
	qbp->setPhase(true);
	qbp->setUseFstSTSolFirst(false);
	restrictlhs.clear();
	restrictrhs=0.0;
        //qbp->searchInitialization(0, (void*)this);
	r = qbp->search(0, (void*)this, search_mode, restrictlhs, restrictrhs, fstStSol, global_score, global_dual_bound, alpha, beta);
	break;
	//if (qbp->getInfoLevel() >= 0) return r;
	//else return (r > (coef_t)0 ? (coef_t)1 : (coef_t)0);
      } else if (search_mode == RESTRICTION_MODE) {
	double dual_bound;
	if (global_score > qbp->getDontKnowValue()) {
	  beta = -qbp->getDontKnowValue();
	  alpha = global_score - fabs(global_score)*0.001 - 0.001;
	  dual_bound = -qbp->getDontKnowValue();
	  qbp->setPhase(false);
	  qbp->setUseFstSTSolFirst(true);
	  //qbp->setFirstStageSolution(fstStSol);

	  restrictlhs.clear();
	  restrictrhs=0.0;
	  for (int g=0; g < fstStSol.size();g++) {
            if (0&&qbp->getType(g) != 0 /*i.e. not BINARY*/) continue;
	    data::IndexedElement e;
	    if (isZero(fstStSol[g])) {
	      e.index = g;
	      e.value = -1.0;
	      restrictlhs.push_back(e);
	    } else if (isOne(fstStSol[g])) {
	      e.index = g;
	      e.value = 1.0;
	      restrictrhs = restrictrhs + 1.0;
	      restrictlhs.push_back(e);
	    }
	  }
	  restrictrhs = restrictrhs - widthFactor * 5; // * (int)log2((double)qbp->nVars());//30.0;
	  cerr << "try LOCAL SEARCH with width " << widthFactor * 5 /* * (int)log2((double)qbp->nVars())*/ << endl;
	  if (widthFactor * 5 /* * (int)log2((double)qbp->nVars())*/ >= qbp->nVars()) pass = -1;
	} else {
	  beta = qbp->getDontKnowValue() / 2.0;//-qbp->getDontKnowValue();
	  alpha = qbp->getNegativeInfinity();
	  dual_bound=-qbp->getDontKnowValue();
	  qbp->setPhase(true);
	  qbp->setUseFstSTSolFirst(false);
	  if (restrictlhs.size() > 0) {
	    restrictrhs = restrictrhs - 20;
	    cerr << "try EXTENDED ROUNDING" << endl;
	  } else cerr << "try ROUNDING" << endl;
	}
	r = qbp->search(0, (void*)this, search_mode, restrictlhs, restrictrhs, fstStSol, global_score, dual_bound, alpha, beta);
	global_dual_bound = qbp->getGlobalDualBound();
	qbp->getFirstStageSolution(fstStSol);
	if (global_score < qbp->getGlobalScore()) {
	  global_score = qbp->getGlobalScore();
	  widthFactor = 1;
	} else {
	  widthFactor = widthFactor * 2;
	}
      } else if (RELAXATION_MODE) {
	qlptmp = qlp;
	const std::vector<data::QpNum> finalObjVec = qlptmp.getObjectiveFunctionValues();
	std::vector<const data::QpRhs *> rhsVec = qlptmp.getRhsVecConst();
	std::vector<const data::Constraint *> conVec = qlptmp.getConstraintVecConst();

	int numConstraints = rhsVec.size();
	if (numConstraints <= 1 && !(numConstraints == 1 && rhsVec[0]->getRatioSign() == data::QpRhs::equal)) {
	  if(qbp->getShowInfo()) std::cerr << "info: relaxation search not possible. Flip to standard mode." << endl;
	  //qbp->searchInitialization(0, (void*)this);
	  return qbp->search(0, (void*)this, search_mode, restrictlhs, restrictrhs, fstStSol, global_score, global_dual_bound, qbp->getNegativeInfinity(), -qbp->getDontKnowValue());
	}
	int i = 0;
	//std::cerr << "SAVE :" << std::endl;
	for (int row = 0; row < numConstraints;row++) {
	  RelaxationBuffer.resize(i+1);
	  RelaxationBuffer[i].first  = conVec[row]->getElements();
	  RelaxationBuffer[i].second = rhsVec[row]->getValue().asDouble();
	  data::QpRhs::RatioSign rs = rhsVec[row]->getRatioSign();
	  //std::cerr << row << "," << i << " : ";
	  //for (int ii = 0; ii < RelaxationBuffer[i].first.size();ii++) {
	  // std::cerr << RelaxationBuffer[i].first[ii].value.asDouble() << "y" << RelaxationBuffer[i].first[ii].index << " + ";
	  //}
	  //std::cerr << " 0 >= " << RelaxationBuffer[i].second << std::endl;

	  if (rs == data::QpRhs::equal || rs == data::QpRhs::smallerThanOrEqual) {
	    for (int ii = 0; ii < RelaxationBuffer[i].first.size();ii++) {
	      RelaxationBuffer[i].first[ii].value = -RelaxationBuffer[i].first[ii].value.asDouble();
	    }
	    RelaxationBuffer[i].second = -RelaxationBuffer[i].second;
	  }
	  i++;
	  if (rs == data::QpRhs::equal) {
	    RelaxationBuffer.resize(i+1);
	    RelaxationBuffer[i].first  = conVec[row]->getElements();
	    RelaxationBuffer[i].second = rhsVec[row]->getValue().asDouble();
	    i++;
	  }
	}

	HT->~HTable();
	if(qbp->getShowInfo()) std::cerr << "Info: ready for relaxation search " << endl;
        std::string StoreInputFilename=qbp->getInputFileName();
	delete qbp;
	qlp.deleteAllRows();
	qlpRelax.deleteAllRows();
	yInit( qlp, StoreInputFilename); //added inputfile	
        for (int i = 0; i < qbp->nVars();i++) {
            data::QpNum zero = 0.0;
	    qlp.setObjectiveFunctionElement(i, zero);
        }
	ySetProcessNo(processNo);
	utils::QlpStageSolver *QlpStSolve = qbp->getStageSolverPt();
        res = qbp->searchInitialization(0, (void*)this);
	if (res < 0) return getNinf();
	qbp->increaseDecisionLevel();
	qbp->setInitilizationTime(ini_time);
	qbp->setPhase(true);

	std::vector<data::IndexedElement> remLhs;
	double remRhs;

	// last row to constraints, second last row to objective
	int lastRow = RelaxationBuffer.size()-1;
	// 
	// last row is put into consraints in following while loop
	//
	// now the objective
	lastRow = RelaxationBuffer.size()-2;

	remLhs.clear();
	remRhs = 0.0;
	for (int i =0;i < RelaxationBuffer[lastRow].first.size();i++) {
	  //std::cerr << RelaxationBuffer[lastRow].first[i].value.asDouble() << "x" << RelaxationBuffer[lastRow].first[i].index << " + ";
	  remLhs.push_back(RelaxationBuffer[lastRow].first[i]);
	}
        remRhs = RelaxationBuffer[lastRow].second;
	//std::cerr << " 0 >= " << remRhs << std::endl;
        for (int i = 0; i < qbp->nVars();i++) {
            data::QpNum zero = 0.0;
	    qlp.setObjectiveFunctionElement(i, zero);
	    qlpRelax.setObjectiveFunctionElement(i, zero);
	    QlpStSolve->changeObjFuncCoeff(qbp->getMaxLPStage(), i, zero);
        }
	//std::cerr << "Minimze ";
        for (int i = 0; i < remLhs.size();i++) {
	  int index = remLhs[i].index;
	  double value = remLhs[i].value.asDouble();
            data::QpNum e = -value;
	    qlp.setObjectiveFunctionElement(index, e);
	    qlpRelax.setObjectiveFunctionElement(index, e);
	    QlpStSolve->changeObjFuncCoeff(qbp->getMaxLPStage(), index, e);
	    //std::cerr << e.asDouble() << "x" << index << " + ";
        }
	//std::cerr << " and beta=" <<  -remRhs << std::endl; 
	beta = -remRhs;
	alpha = qbp->getNegativeInfinity();
	
	while (RelaxationBuffer.size() > 0) {
	  // now put last row to qlp. => qlp has correct information, also in case of re-inits. ASSUMPTION: QlpStageSolver works with copy
	  int lastRow = RelaxationBuffer.size()-1;
	  remLhs.clear();
	  remRhs = 0.0;
	  for (int i =0;i < RelaxationBuffer[lastRow].first.size();i++) {
	    remLhs.push_back(RelaxationBuffer[lastRow].first[i]);
	  }
	  remRhs = RelaxationBuffer[lastRow].second;
	  //std::cerr << "add constraint at row " << lastRow << ": ";
	  //for (int i =0;i < remLhs.size();i++) {
	  //  std::cerr << remLhs[i].value.asDouble() << "x" << remLhs[i].index << " + ";
	  //}
	  //std::cerr << " 0 >= " << remRhs << std::endl;

	  data::QpRhs rhs;
	  rhs.setValue(remRhs);
	  rhs.setRatioSign(data::QpRhs::greaterThanOrEqual);
	  data::Constraint& c = qlp.createRhsConstraint(rhs);
	  c.setElements(remLhs);
	  data::Constraint& cc = qlpRelax.createRhsConstraint(rhs);
	  cc.setElements(remLhs);

	  qlptmp = qlp;

	  //std::cerr << "READY for relaxation search again" << endl;
	  std::string StoreInputFilename=qbp->getInputFileName();
	  delete qbp;
	  qlp = qlptmp;
	  //qlpRelax = qlp;
	  yInit( qlp, StoreInputFilename);//added inputfile	
	  ySetProcessNo(processNo);
	  utils::QlpStageSolver *QlpStSolve = qbp->getStageSolverPt();
	  res = qbp->searchInitialization(0, (void*)this);
	  if (res < 0) return getNinf();
	  qbp->increaseDecisionLevel();
	  qbp->setInitilizationTime(ini_time);
	  if (/*global_score < qbp->getDontKnowValue()*/pass==1) {
	    qbp->setPhase(true);
	    qbp->setUseFstSTSolFirst(false);
	  } else {
	    qbp->setPhase(false);
	    qbp->setUseFstSTSolFirst(true);
	  }

	  //std::cerr << "pre IN QLP" << std::endl;
	  int maxLPStage = qbp->getMaxLPStage();
	  vector<data::QpRhs> rhsVec;
	  rhsVec = qlp.getRhsVec( );
	  std::vector<data::Constraint*> lhss = qlp.getConstraintVec(); 
	  for (int i = 0; i < rhsVec.size();i++) {
	    double rhstmp = rhsVec.at( i ).getValue().asDouble();
	    std::vector<data::IndexedElement> rowtmp = lhss[i]->getElements();
	    
	    QlpStSolve->getExternSolver( maxLPStage ).addLProw_snapshot(rowtmp, rhsVec[i]);
	    if (0) {
	      double lhs = 0.0, rhs = rhstmp;
	      cerr << "real C-"<< i << ": ";
	      for (int h = 0; h < rowtmp.size();h++) {
		if (fstStSol.size() > rowtmp[h].index) lhs = lhs + rowtmp[h].value.asDouble() * fstStSol[rowtmp[h].index];
		cerr << rowtmp[h].value.asDouble() << "x" << rowtmp[h].index << " + ";
	      }
	      if (rhsVec.at( i ).getRatioSign() == data::QpRhs::greaterThanOrEqual)
		cerr << " 0 >= " << rhstmp << endl; 
	      else if (rhsVec.at( i ).getRatioSign() == data::QpRhs::smallerThanOrEqual)
		cerr << " 0 <= " << rhstmp << endl; 
	      else if (rhsVec.at( i ).getRatioSign() == data::QpRhs::equal)
		cerr << " 0 == " << rhstmp << endl; 
	      if (lhs < rhs) {
		std::cerr << "ERROR: lhs=" << lhs << " < " << rhs << "=rhs" << std::endl;
	      }
	    }
	  }




	  remLhs.clear();
	  remRhs = 0.0;
	  //std::cerr << "THIS SHOULD WORK:";
	  //for (int i = 0; i < fstStSol.size(); i++) {
	  //  std::cerr << " " << fstStSol[i];
	  //}
	  std::cerr << std::endl;
	  if (pass < 0)
	    qbp->setOutputSupport(false);
	  else
	    qbp->setOutputSupport(true);
	  r = qbp->search(0, (void*)this, search_mode, remLhs,  remRhs, fstStSol, global_score, global_dual_bound, alpha, beta);
	  if (qbp->getShowInfo())std::cerr << "info: return " << r << " gs=" << global_score << " beta=" << beta << std::endl;
	  if (pass < 0) break;
	      if (0) {
		std::cerr << "post IN QLP" << std::endl;
		int maxLPStage = qbp->getMaxLPStage();
		vector<data::QpRhs> rhsVec;
		rhsVec = qlp.getRhsVec( );
		std::vector<data::Constraint*> lhss = qlp.getConstraintVec(); 
		for (int i = 0; i < rhsVec.size();i++) {
		  double rhstmp = rhsVec.at( i ).getValue().asDouble();
		  //finde raus welche zeilen
		  std::vector<data::IndexedElement> rowtmp = lhss[i]->getElements();
		  //qlp.getRowLhs( i, rowtmp );
		  cerr << "real C-"<< i << ": ";
		  double lhs = 0.0, rhs = rhstmp;
		  cerr << "real C-"<< i << ": ";
		  for (int h = 0; h < rowtmp.size();h++) {
		    if (fstStSol.size() > rowtmp[h].index) lhs = lhs + rowtmp[h].value.asDouble() * fstStSol[rowtmp[h].index];
		    cerr << rowtmp[h].value.asDouble() << "x" << rowtmp[h].index << " + ";
		  }

		  if (rhsVec.at( i ).getRatioSign() == data::QpRhs::greaterThanOrEqual)
		    cerr << " 0 >= " << rhstmp << endl; 
		  else if (rhsVec.at( i ).getRatioSign() == data::QpRhs::smallerThanOrEqual)
		    cerr << " 0 <= " << rhstmp << endl; 
		  else if (rhsVec.at( i ).getRatioSign() == data::QpRhs::equal)
		    cerr << " 0 == " << rhstmp << endl; 
		  if (lhs < rhs) {
		    std::cerr << "ERROR: lhs=" << lhs << " < " << rhs << "=rhs" << std::endl;
		  }
		}
		std::cerr << "MINMIZE ";
		for (int ii = 0; ii < qbp->nVars();ii++) {
		  data::QpNum zero = qlp.getObjectiveFunctionElement(ii);
		  cerr << zero.asDouble() << "x" << ii << " + ";
		}
		std::cerr << std::endl;

	      }
	      if (0) {
		std::cerr << "Ext Sol cont" << std::endl;
		int maxLPStage = qbp->getMaxLPStage();
		vector<data::QpRhs> rhsVec;
		QlpStSolve->getExternSolver( maxLPStage ).getRhs( rhsVec );
		for (int i = 0; i < QlpStSolve->getExternSolver( maxLPStage ).getRowCount();i++) {
		  double rhstmp = rhsVec.at( i ).getValue().asDouble();
		  //finde raus welche zeilen
		  std::vector<data::IndexedElement> rowtmp;
		  QlpStSolve->getExternSolver( maxLPStage ).getRowLhs( i, rowtmp );
		  cerr << "real C-"<< i << ": ";
		  for (int h = 0; h < rowtmp.size();h++) {
		    cerr << rowtmp[h].value.asDouble() << "x" << rowtmp[h].index << " + ";
		  }
		  if (rhsVec.at( i ).getRatioSign() == data::QpRhs::greaterThanOrEqual)
		    cerr << " 0 >= " << rhstmp << endl; 
		  else if (rhsVec.at( i ).getRatioSign() == data::QpRhs::smallerThanOrEqual)
		    cerr << " 0 <= " << rhstmp << endl; 
		  else if (rhsVec.at( i ).getRatioSign() == data::QpRhs::equal)
		    cerr << " 0 == " << rhstmp << endl; 
		}
	      }
	      {
		//std::cerr << "im SNAPSHOT" << std::endl;
		int maxLPStage = qbp->getMaxLPStage();
		//vector<data::QpRhs> rhsVec;
		//QlpStSolve->getExternSolver( maxLPStage ).getRhs( rhsVec );
		for (int i = 0; i < (*QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot()).size();i++) {
		  double rhstmp = (*QlpStSolve->getExternSolver( maxLPStage ).getRowRhs_snapshot())[i].getValue().asDouble();//rhsVec.at( i ).getValue().asDouble();
		   //double rhstmp = rhsVec.at( i ).getValue().asDouble();
		  //finde raus welche zeilen
		  //std::vector<data::IndexedElement> rowtmp;
		  std::vector<data::IndexedElement> &rowtmp = *QlpStSolve->getExternSolver( maxLPStage ).getRowLhs_snapshot(i);
		   //QlpStSolve->getExternSolver( maxLPStage ).getRowLhs( i, rowtmp );
		  //cerr << "real C-"<< i << ": ";
		  //for (int h = 0; h < rowtmp.size();h++) {
		  //  cerr << rowtmp[h].value.asDouble() << "x" << rowtmp[h].index << " + ";
		  //}
		  //cerr << " <==> " << rhstmp << endl;
		}
	      }

	   

	  if (r > beta) {
	    //std::cerr << "RETURN as r" << r << ">" <<  beta << "=beta" << std::endl;
	    return qbp->getNegativeInfinity();
	  }
	  qbp->getFirstStageSolution(fstStSol);
	  if (0&&global_dual_bound > qbp->getGlobalDualBound()) {
	    global_dual_bound = qbp->getGlobalDualBound();
	  } 
	  //std::cerr << "RELAXATION run finished." << std::endl;
	  // Ok. Last row no more needed.
	  RelaxationBuffer[RelaxationBuffer.size()-1].first.clear();
	  RelaxationBuffer.pop_back();
	  
	  // add next constraint being the objective
	  for (int i = 0; i < qbp->nVars();i++) {
            data::QpNum zero = 0.0;
	    qlp.setObjectiveFunctionElement(i, zero);
	    qlpRelax.setObjectiveFunctionElement(i, zero);
	    QlpStSolve->changeObjFuncCoeff(qbp->getMaxLPStage(), i, zero);
	  }
	  lastRow = RelaxationBuffer.size()-2;
	  if (qbp->getShowInfo()) std::cerr << " Remaining iterations " << lastRow << std::endl;
	  //std::cerr << "LAST ROW " <<  lastRow << std::endl;
	  if (lastRow < 0) {
	    //std::cerr << "FINALE!!!!!" << std::endl;
	    // add original objective
	    //std::cerr << "Minimze III";
	    for (int i = 0; i < finalObjVec.size();i++) {
	      int index = i;
	      double value = finalObjVec[i].asDouble();
	      data::QpNum e = value;
	      qlp.setObjectiveFunctionElement(index, e);
	      qlpRelax.setObjectiveFunctionElement(index, e);
	      QlpStSolve->changeObjFuncCoeff(qbp->getMaxLPStage(), index, e);
	      //std::cerr << e.asDouble() << "x" << index << " + ";
	    }
	    //std::cerr << " and beta=" <<  -qbp->getDontKnowValue() << std::endl; 
	    
	    {
	      int maxLPStage = qbp->getMaxLPStage();
	      vector<data::QpRhs> rhsVec;
	      rhsVec = qlp.getRhsVec( );
	      std::vector<data::Constraint*> lhss = qlp.getConstraintVec(); 
	      for (int i = 0; i < rhsVec.size();i++) {
		double rhstmp = rhsVec.at( i ).getValue().asDouble();
		std::vector<data::IndexedElement> rowtmp = lhss[i]->getElements();
		
		QlpStSolve->getExternSolver( maxLPStage ).addLProw_snapshot(rowtmp, rhsVec[i]);
	      }
	    }
	    
	    
	    alpha = qbp->getDontKnowValue();
	    beta = -qbp->getDontKnowValue();
	    //remLhs.clear();
	    //remRhs = 0.0;
	    //r = qbp->search(0, (void*)this, search_mode, remLhs, remRhs, fstStSol, global_score, global_dual_bound, alpha, -qbp->getDontKnowValue());
	    pass = -1;
	    //break;
	  } else {
	    remLhs.clear();
	    remRhs=0.0;
	    for (int i =0;i < RelaxationBuffer[lastRow].first.size();i++) {
	      remLhs.push_back(RelaxationBuffer[lastRow].first[i]);
	    }
	    remRhs = RelaxationBuffer[lastRow].second;
	    //std::cerr << "Minimze II";
	    for (int i = 0; i < remLhs.size();i++) {
	      int index = remLhs[i].index;
	      double value = remLhs[i].value.asDouble();
	      data::QpNum e = -value;
	      qlp.setObjectiveFunctionElement(index, e);
	      qlpRelax.setObjectiveFunctionElement(index, e);
	      QlpStSolve->changeObjFuncCoeff(qbp->getMaxLPStage(), index, e);
	      //std::cerr << e.asDouble() << "x" << index << " + ";
	    }
	    //std::cerr << " and beta=" <<  -remRhs << std::endl; 
	    alpha = qbp->getDontKnowValue();
	    beta = -remRhs;
	    pass++;
	  }
	}
      }

      if (/*pass == 10 ||*/ pass < 0) break;

      while (qbp->getVarPriorityQueueSize() > 0) {
	qbp->extractVarPriorityQueueMinimum();
      }
      HT->~HTable();
      pass++;
      std::cerr << "RESTART " << pass-1 << endl;
      std::string StoreInputFilename=qbp->getInputFileName();
      delete qbp;
      yInit( qlp, StoreInputFilename); // add inputfile	
      ySetProcessNo(processNo);

    } while (1);
    std::vector<int> pt2leaders;
    std::vector<int> bitcnts;
    for (int i = 0; i < qbp->nVars(); i++) {
      bitcnts.push_back(integers[i].bitcnt);
      pt2leaders.push_back(integers[i].pt2leader);
    }
    if ((processNo & 1) == 0) {
      if (fabs(qbp->getGlobalScore()) < -qbp->getNegativeInfinity()) {
	double min_slack = qbp->finalSolutionCheck(bitcnts,pt2leaders,objInverted);
	for (int i = 0; i < qbp->nVars(); i++) {
	  //cerr << "x" << i << "(" << qbp->getType(i) << "," << qbp->getBlock(i) << "," << integers[i].bitcnt << ") ";
	  if (qbp->getType(i) == 5000 /*CONTINUOUS*/) {
	    if (qbp->getLowerBound(i) == qbp->getUpperBound(i)) {
	      cerr << qbp->getLowerBound(i) << " ";
	    } else if (qbp->getBlock(i) == 1) cerr << qbp->getFirstStageSolutionValue(i) << " ";
	    else cerr << "c" <<  qbp->getBlock(i) << " ";
	  } else if (qbp->getBlock(i) == 1){
	    if (integers[i].bitcnt > 1) {
	      //cerr << "Var " << i<< " hat " << integers[i].bitcnt << "bits. Leader:" << integers[i].pt2leader;
	      if (integers[i].pt2leader == i) {
		int res = 0;
		for (int z = 0;z < integers[i].bitcnt;z++) {
		  res = 2*res + qbp->getFirstStageSolutionValue(i+z);
		  //cerr << " | " << qbp->getFirstStageSolutionValue(i+z) << " | ";
		}
		cerr << res << " ";
		i+=integers[i].bitcnt;
		i--;
	      } else cerr << "Error: missing bit" << endl;
	    } else
	      cerr << qbp->getFirstStageSolutionValue(i) << " ";
	  } else {
	    cerr << "s" << qbp->getBlock(i) << " ";
	  }
	}
	cerr << endl;
      }
    }

    if (qbp->getInfoLevel() >= -OUT_OFFSET) return r;
    else return (r > (coef_t)0 ? (coef_t)1 : (coef_t)0);
  }

	yInterface() {
		yParamMaxHashSize = 10000000;
		timeout = time(NULL) + 1000000000;
		info_level = 1;
		qbp = 0;
	    result = (coef_t)0;
	    qbp = 0;
	    QlpStSolve = 0;
	    LPoffset = 0.0;
	    confl = CRef_Undef;
	    confl_partner = CRef_Undef;
	    offset_shift = 0.0;
	    confl_var = 0;
	}

	double computeCutRatio(vector< pair<unsigned int, double> >& cut);
	double finalSolutionCheck();

#ifndef FIND_BUG
	virtual int GenerateCutAndBranchCuts( extSol::QpExternSolver& externSolver, vector< vector< data::IndexedElement > > &listOfCutsLhs,
			       vector< data::QpNum > &listOfCutsRhs, vector<int> &listOfCutsVars,
					      int treedepth, int currentBlock, bool &global_valid, std::vector<unsigned int> &candidates, int cuttype, int* types, int8_t *assigns, unsigned int initime, int* solu, int *fixs, int *blcks, int orgN, double intLB);
#else
	virtual int GenerateCutAndBranchCuts( extSol::QpExternSolver& externSolver, vector< vector< data::IndexedElement > > &listOfCutsLhs,
			       vector< data::QpNum > &listOfCutsRhs,
					      int treedepth, int currentBlock, bool &global_valid, std::vector<unsigned int> &candidates, int cuttype, int* types, int8_t *assigns, unsigned int initime, int orgN, double intLB);
#endif
	virtual double getReducedCostDj(int ix, bool &l, bool &u, bool &loUp);
	virtual void getRCandB(extSol::QpExternSolver& externSolver);
	virtual int isReducedCostFixed( double z, double lpval, int ix, double solx, double &rc);
	virtual void setBranchingDecision(int &pick, int &left, int &right);
	virtual void moveUp(coef_t& value, coef_t uBound, int status);
	virtual void moveDown(int toLevel, int pick, int val, int val_ix);

	std::vector<double> getFirstStageSolution();
	~yInterface() {
	  //std::cerr << "delete CBC" << std::endl;
#ifndef NO_CGL
	  delete cbc;
#endif
	  //std::cerr << "delete QBP" << std::endl;
	  delete qbp;
	  //delete QlpStSolve;
	}
};


#endif /* YINTERFACE_H_ */
