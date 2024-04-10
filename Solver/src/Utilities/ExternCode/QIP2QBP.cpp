/*
*
* Solver: QIP2QBP.cpp -- Copyright (c) 2010-2017 Jan Wolf
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

#include "Utilities/ExternCode/QIP2QBP.h"
//#include <boost/algorithm/string.hpp>


QIP2QBP::QIP2QBP(void) : qlp(), qbp(), existsVariablesQlp(), ExistBlockLength(),
ExistBlockStart(),existIndexByIndex(), firstIndexBByIndex(), indexDByIndex(), NumberExistBlock(),numsBinaries(), onlySum(true),
totalNumBinariesQbp(),totalNumDsQlp(){
}

QIP2QBP::~QIP2QBP(void) {
}

void QIP2QBP::readQIP(string path) {
	std::string inputFile(path);
	utils::Parser::createQlp(inputFile, qlp);

}
void QIP2QBP::makeQBP() {

	prepareBs();
	calcExistBlocks();
	initExistIndexByIndex();
	initIndexDByIndex();
	initFirstIndexBByIndex();

	int pos = 0;
	int index = 0;

	cout << "Variablen anlegen / kopieren" << endl;

	std::vector<data::QpVar*> vVec = qlp.getVariableVector();

	//Variablen anlegen / kopieren
	for (unsigned int i = 0, size = qlp.getVariableCount(); i < size; i++) {
		if (vVec[i]->getQuantifier() == data::QpVar::exists) {
			for (int j = 0; j < numsBinaries.at(pos); j++) {
				stringstream s;
				s << "b_" << index;
				data::QpVar var = data::QpVar(s.str(), index, 0, 1,
						data::QpVar::binaries, data::QpVar::exists);
				qbp.createVariable(var);
				index++;
			}
			pos++;
		} else {
			stringstream s;
			s << "D_" << index;
			if (vVec[i]->getQuantifier()
					== data::QpVar::random) {
				cout << "Random-Variablen im Problem" << endl;
			}
			//cout<<qlp.getVariableByIndex(i).getNumberSystem()<<endl;
			data::QpVar var = data::QpVar(s.str(), index,
					vVec[i]->getLowerBound(),
					vVec[i]->getUpperBound(),
					vVec[i]->getNumberSystem(),
					vVec[i]->getQuantifier());
			qbp.createVariable(var);
			index++;
		}

	}

	if (!onlySum) {
		for (unsigned int i = 0; i < qlp.getVariableVector().size(); i++) {
			if (vVec[i]->getQuantifier()
					== data::QpVar::exists) {
				stringstream s;
				s << "p_" << index;
				data::QpVar var = data::QpVar(s.str(), index,
						vVec[i]->getLowerBound(),
						vVec[i]->getUpperBound(),
						data::QpVar::real, data::QpVar::exists);
				qbp.createVariable(var);
				index++;
			}
		}
	}

	cout << "Zielfunktion anlegen / kopieren" << endl;
	//Zielfunktion anlegen / kopieren
	std::vector<data::IndexedElement> objCoefQlp =
			qlp.getObjectiveFunction().getObjectiveElementsSparse();

	data::QpNum value;

	qbp.setObjective(qlp.getObjective());

	for (unsigned int i = 0; i < objCoefQlp.size(); i++) {
		data::QpVar* var = &qlp.getVariableByIndex(objCoefQlp.at(i).index);

		if (var->getQuantifier() == data::QpVar::exists && onlySum) {

			// Zielfunktionskoeffizienten für b_tj anlegen
			for (int j = 0; j < numsBinaries.at(getIndexExist(var->getIndex())); j++) {
				index = getFirstIndexB(getIndexExist(var->getIndex())) + j;
				(value = objCoefQlp.at(i).value.asDouble()) *= pow(2.0, j);
				qbp.setObjectiveFunctionElement(index, value);

			}
		} else {
			//Zielfunktionskoeffizienten kopieren
			if (var->getQuantifier() == data::QpVar::exists) {
				index = getIndexP(var->getIndex());
			} else {
				index = getIndexD(var->getIndex());
			}
			value = objCoefQlp.at(i).value.asDouble();
			qbp.setObjectiveFunctionElement(index, value);
		}
	}

	//Methode zum Debuggen, gibt alle Variablen mit index und Name aus
	//print();

	std::vector<data::Constraint const *> cVec;
	std::vector<data::IndexedElement> ieVec;

	cout << "Constraints anlegen / kopieren" << endl;
	//Constraints anlegen / kopieren
	if (onlySum) {

		cVec = qlp.getConstraintVecConst();

		// Constraints kopieren mit Ersetzen der Existvariablen durch Binärdarstellung
		for (unsigned int i = 0; i < cVec.size(); i++) {

			data::Constraint& c = qbp.createRhsConstraint(cVec[i]->getRhs());

			//c.setRhsRatioSign(cVec[i]->getRhsRatioSign());

			unsigned int index;
			data::QpNum value;

			cVec[i]->getElementsSparse(ieVec);


			for (unsigned int k = 0; k
					< ieVec.size(); k++) {
				data::QpVar* var =
						&qlp.getVariableByIndex(ieVec[k].index);
				if (var->getQuantifier() == data::QpVar::exists) {
					//p_t bzw. ersetzt durch b_it
					for (int j = 0; j < numsBinaries.at(getIndexExist(
							var->getIndex())); j++) {
						index = getFirstIndexB(getIndexExist(var->getIndex()))
								+ j;
						value = pow(2.0, j) * ieVec[k].value.asDouble();
						c.createConstraintElement(index, value);
					}
				} else {
					//D_t
					index = getIndexD(var->getIndex());
					value = ieVec[k].value;
					c.createConstraintElement(index, value);
				}
			}
		}

		//Bounds der Existvariablen auf Binärdarstellung übertragen:
		for (unsigned int i = 0; i < qlp.getVariableCount(); i++) {
			data::QpVar* var = &qlp.getVariableByIndex(i);
			if(var->getQuantifier() == data::QpVar::exists){

				// Obere Schranke
				data::Constraint& upper =
					qbp.createRhsConstraint(data::QpRhs::smallerThanOrEqual, var->getUpperBound()); //(data::QpRhs::equal,rhs-QPnum)
				//Untere Schranke
				data::Constraint& lower =
								qbp.createRhsConstraint(data::QpRhs::greaterThanOrEqual, var->getLowerBound()); //(data::QpRhs::equal,rhs-QPnum)
				int index;
				data::QpNum value;
				// b_it
				for (int j = 0; j < numsBinaries.at(getIndexExist(i)); j++) {
					index = getFirstIndexB(getIndexExist(i)) + j;
					value = pow(2.0, j);
					upper.createConstraintElement(index, value);
					lower.createConstraintElement(index, value);
				}
			}
		}

	} else {

		throw utils::QlpSolverException("WRONG TURN");

		// Alle Constraints aus qlp in qbp kopieren
		for (unsigned int i = 0; i < qlp.getConstraintVec().size(); i++) {
			data::Constraint& c = qbp.createRhsConstraint(
					qlp.getConstraintVec().at(i)->getRhs());
			c.setRhsRatioSign(qlp.getConstraintVec().at(i)->getRhsRatioSign());
			unsigned int index;
			data::QpNum value;
			for (unsigned int j = 0; j
					< qlp.getConstraintVec().at(i)->getElementCount(); j++) {
				data::QpVar* var =
						&qlp.getVariableByIndex(
								qlp.getConstraintVec().at(i)->getElements().at(
										j).index);
				if (var->getQuantifier() == data::QpVar::exists) {
					index = getIndexP(var->getIndex());
				} else {
					index = getIndexD(var->getIndex());
				}
				value = qlp.getConstraintVec().at(i)->getElements().at(
						j).value;
				c.createConstraintElement(index, value);
			}

		}

		//  Neue anlegen q_t = sum_i b_it * 2^i
		for (unsigned int i = 0; i < existsVariablesQlp.size(); i++) {
			data::Constraint& c =
					qbp.createRhsConstraint(data::QpRhs::equal, 0); //(data::QpRhs::equal,rhs-QPnum)

			int index;
			data::QpNum value;
			// b_it
			for (int j = 0; j < numsBinaries.at(i); j++) {
				index = getFirstIndexB(i) + j;
				value = pow(2.0, j);
				c.createConstraintElement(index, value);
			}
			// ursprüngliche Existvariable q_t
			index = totalNumBinariesQbp + totalNumDsQlp + i;
			value = -1;
			c.createConstraintElement(index, value);
		}

	}

}

void QIP2QBP::prepareBs() {
	totalNumBinariesQbp = 0;
	existsVariablesQlp = qlp.getVariableVectorByQuantifierConst(
			data::QpVar::exists);

	// Strukturierung der Bin�rvariablen:
	for (unsigned int i = 0; i < existsVariablesQlp.size(); i++) {
		//cout << "ub: " << existsVariablesQlp.at(i)->getUpperBound().asDouble() <<endl;

		if(existsVariablesQlp.at(i)->getUpperBound().asDouble() >=1000000000000 || existsVariablesQlp.at(i)->getUpperBound().asDouble() < 0){
			cout<<"unbeschränkte Variable: " << existsVariablesQlp.at(i)->getName() <<endl;
			exit(0);
		}

		numsBinaries.push_back(countBs(
				existsVariablesQlp.at(i)->getUpperBound().asDouble()));

		//cout << "i: " << i << "---> " << countBs(existsVariablesQlp.at(i)->getUpperBound().asDouble())<<endl;
	}

	// Anzahl der ben�tigten Bin�rvariablen bestimmen:
	for (unsigned int i = 0; i < numsBinaries.size(); i++) {
		totalNumBinariesQbp = totalNumBinariesQbp + numsBinaries.at(i);
	}
}

// Z�hlt, wie viele Bin�rvariablen zur Darstellung einer Existenzvariable mit der �bergegebenen oberen Schranke ben�tigt werden. 
int QIP2QBP::countBs(int boundExVar) {
	if (boundExVar == 0) {
		return 1;
	}
	int i = 0;
	int j = 0;
	while (j < boundExVar) {
		j = j + pow(2.0, i);
		i++;
	}
	return i;
}

void QIP2QBP::initExistIndexByIndex(){

	int i = 0;
		for (int j = 0; j < qlp.getVariableCount(); j++) {
			if (qlp.getVariableByIndex(j).getQuantifier() == data::QpVar::exists) {
				this->existIndexByIndex.push_back(i);
				i++;
			}else{
				this->existIndexByIndex.push_back(-1);
			}
		}
}

void QIP2QBP::initFirstIndexBByIndex(){
	for (unsigned int j = 0; j < qlp.getVariableCount(); j++) {
		this->firstIndexBByIndex.push_back(this->getFirstIndexB_init(j));
	}
}

void QIP2QBP::initIndexDByIndex(){
	int indexReturn = 0;
	int i = 0;
	for (int j = 0; j < qlp.getVariableCount(); j++) {
		if (qlp.getVariableByIndex(j).getQuantifier() == data::QpVar::exists) {
			indexReturn = indexReturn + numsBinaries.at(i);
			i++;
		} else {
			indexReturn++;
		}
		this->indexDByIndex.push_back(indexReturn);
	}
}


void QIP2QBP::calcExistBlocks() {
	int elementNumber = 0;
	NumberExistBlock = 1;
	ExistBlockStart.push_back(0);
	for (unsigned int i = 0; i < existsVariablesQlp.size() - 1; i++) {
		if (existsVariablesQlp.at(i)->getIndex()
				!= existsVariablesQlp.at(i + 1)->getIndex() - 1) {
			ExistBlockLength.push_back(elementNumber);
			ExistBlockStart.push_back(existsVariablesQlp.at(i + 1)->getIndex());
			elementNumber = 0;
			NumberExistBlock++;
		}
		elementNumber++;
	}
	ExistBlockLength.push_back(elementNumber);

	//Berechnung totalNumDs
	totalNumDsQlp = qlp.getVariableCount() - existsVariablesQlp.size();
}

int QIP2QBP::getIndexExist(int index) {
	return this->existIndexByIndex[index]; //precomputed
//	int i = 0;
//	for (int j = 0; j < index; j++) {
//		if (qlp.getVariableByIndex(j).getQuantifier() == data::QpVar::exists) {
//			i++;
//		}
//	}
//	return i;
}

int QIP2QBP::getFirstIndexB(int index) {
	return this->firstIndexBByIndex[index]; //precomputed
//	int i = 0;
//	int indexReturn = 0;
//	for (unsigned int j = 0; j < qlp.getVariableCount(); j++) {
//		if (index > i && qlp.getVariableByIndex(j).getQuantifier()
//				== data::QpVar::exists) {
//			indexReturn = indexReturn + numsBinaries.at(i);
//			i++;
//		} else if (index == i && qlp.getVariableByIndex(j).getQuantifier()
//				== data::QpVar::exists) {
//			i++;
//		} else if (index >= i) {
//			indexReturn++;
//		}
//	}
//	return indexReturn;
}

int QIP2QBP::getFirstIndexB_init(int index) {

	int i = 0;
	int indexReturn = 0;

	for (unsigned int j = 0; j < qlp.getVariableCount(); j++) {
		if (index > i && qlp.getVariableByIndex(j).getQuantifier()
				== data::QpVar::exists) {
			indexReturn = indexReturn + numsBinaries.at(i);
			i++;
		} else if (index == i && qlp.getVariableByIndex(j).getQuantifier()
				== data::QpVar::exists) {
			i++;
		} else if (index >= i) {
			indexReturn++;
		}
	}
	return indexReturn;
}


int QIP2QBP::getIndexD(int index) {
	return this->indexDByIndex[index];//precomputed
//	int indexReturn = 0;
//	int i = 0;
//	for (int j = 0; j < index; j++) {
//		if (qlp.getVariableByIndex(j).getQuantifier() == data::QpVar::exists) {
//			indexReturn = indexReturn + numsBinaries.at(i);
//			i++;
//		} else {
//			indexReturn++;
//		}
//	}
//	return indexReturn;
}

int QIP2QBP::getIndexP(int index) {
	int indexReturn = 0;
	indexReturn = totalNumBinariesQbp + totalNumDsQlp + getIndexExist(index);
	return indexReturn;
}

void QIP2QBP::print() {
	for (unsigned int i = 0; i < qbp.getVariableVector().size(); i++) {
		cout << "index " << qbp.getVariableVector().at(i)->getIndex() << endl;
		cout << "name " << qbp.getVariableVector().at(i)->getName() << endl;
	}
}

/*//Main-Methode für Consolen-Aufruf
int main( int argc, char *argv[] ) {

	if( argc < 3 ){
		cout << "Mindestens 2 Argumente benötigt. \n 1. Argument: Input-Dateipfad \n 2. Argument: 1, wenn alle Existvariablen q_t aus dem Eingabe-LP im Modell nur als Binärdarstellung vorkommen sollen." << endl;
		exit( EXIT_FAILURE );
	}

	string fileIn = argv[1];
	string umw = argv[2];

	//vector <string> fields;
    //boost::split( fields, fileIn, "." );
	//string dateiName = fields.at(0);


	stringstream st(fileIn);
	string dateiName;
	getline(st, dateiName, '.');


	string fileOutQBP = dateiName + "_QBP_" + umw +".qlp";
	string fileOutDEPBin = dateiName + "_BinDEP_" + umw +".lp";
	string fileOutDEP = dateiName + "_DEP" +".lp";
	string losung = dateiName +"_" + umw + "Loesung.txt";


	std::cout << "foobar" << std::endl;
	QIP2QBP q = QIP2QBP();
	q.onlySum = ( umw == "1" );

	q.qlp = QLPReader::read( fileIn.c_str(), "", "_" );
	//q.readQIP(fileIn);

	 //ofstream myfile1;
	 //myfile1.open( "/home/utz/Documents/Johannes/TEST.qlp" );   // Schreibt original QLP raus
	 //myfile1 << q.qlp.toQlpFileString(false);
	 //myfile1.close();


	//q.makeQBP();

	//qbp in Datei rausschreiben;
//	 ofstream myfile;
//	 myfile.open( fileOutQBP.c_str() );
//	 myfile << q.qbp.toQlpFileString(false);
//	 myfile.close();

	// cout << "QBP-Datei geschrieben" << endl;


	 // DEP schreiben
	 data::Qlp dep;
	 utils::QlpConverter::convertToLP(q.qlp,dep,utils::QlpConverter::COMPACT_VIEW,utils::QlpConverter::WORST);
	 utils::ToolBox::writeToFile(fileOutDEP,dep.toQlpFileString(true));


	 // Binäres DEP schreiben
//	 data::Qlp depBin;
//	 utils::QlpConverter::convertToLP(q.qbp,depBin,utils::QlpConverter::COMPACT_VIEW,utils::QlpConverter::WORST);
//	 utils::ToolBox::writeToFile(fileOutDEPBin,depBin.toQlpFileString(true));

	 cout << "DEP-Dateien geschrieben" << endl;




//	 //Lösung in Datei rausschreiben;
//	 ofstream myfile2;
//	 myfile2.open( losung.c_str() );
//	 myfile2 << "Lösung von " << fileIn.c_str()<<endl;
//	 cout << "Lösung von " << fileIn.c_str()<<endl;
//
//	 {
//		 algorithm::Algorithm::QlpSolution s = algorithm::Qlp2Lp( q.qbp ).solveQlp( algorithm::Algorithm::WORST_CASE );
//		 myfile2 << "Solution status: " << s.getSolutionStatusString()<< endl;
//		 cout << "Solution status: " << s.getSolutionStatusString()<<endl;
//		 if( s.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL || s.getSolutionStatus() == extSol::QpExternSolver::NUM_BEST ){
//			 myfile2 << "Objective: " << s.getObjFunctionValue().asDouble() << endl;
//			 myfile2 << "Solution-Size: " << s.getSolutionVector().size()<< endl;
//			 myfile2 << "Solution: "<< endl;
//			 for(unsigned int i = 0; i < s.getSolutionVector().size();i++){
//				 myfile2 << i << " : " << s.getSolutionVector().at(i).asDouble()<< endl;
//			 }
//		 }
//	 }
//	 {
//		 algorithm::Algorithm::QlpSolution s = algorithm::Qlp2Lp( q.qlp ).solveQlp( algorithm::Algorithm::WORST_CASE );
//		 myfile2 << "----------------------------------------------------------------------------------------------"<< endl;
//		 myfile2 << "Org. Solution status: " << s.getSolutionStatusString()<< endl;
//
//		 if( s.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL || s.getSolutionStatus() == extSol::QpExternSolver::NUM_BEST ){
//			 myfile2 << "Org. Objective: " << s.getObjFunctionValue().asDouble()<< endl;
//			 myfile2 << "Org. Solution-Size: " << s.getSolutionVector().size() << endl;
//			 myfile2 << "Org. Solution: "<< endl;
//			 for(int i = 0; i < s.getSolutionVector().size();i++){
//				 myfile2 << i << " : " << s.getSolutionVector().at(i).asDouble()<< endl;
//			}
//		 }
//	 }
//	 myfile2.close();


	 cout << "Fertig" << endl;


	return EXIT_SUCCESS;
}


//Main-Methode für Debug
//int main() {
//
//	string fileIn = "/home/utz/Documents/Johannes/Original/modelQuantified2.lp";
//	string umw = "1";
//	string fileOut = fileIn + "QBP" + umw +".qlp";
//
//
//	std::cout << "foobar" << std::endl;
//	QIP2QBP q = QIP2QBP();
//	q.onlySum = ( umw == "1" );
//
//	q.qlp = QLPReader::read( fileIn.c_str(), "", "_" );
//
//	//q.readQIP(fileIn);
//	 ofstream myfile1;
//	 myfile1.open( "/home/utz/Documents/Johannes/TEST.qlp" );
//	 myfile1 << q.qlp.toQlpFileString(false);
//	 myfile1.close();
//
//	q.makeQBP();
//
//	//qbp in Datei rausschreiben;
//	 ofstream myfile;
//	 myfile.open( fileOut.c_str() );
//	 myfile << q.qbp.toQlpFileString(false);
//	 myfile.close();
//
//	 cout << "QBP-Datei geschrieben" << endl;
//
//	 {
//		 algorithm::Algorithm::QlpSolution s = algorithm::Qlp2Lp( q.qbp ).solveQlp( algorithm::Algorithm::WORST_CASE );
//		 cout << "Solution status: " << s.getSolutionStatusString() << endl;
//
//		 if( s.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL || s.getSolutionStatus() == extSol::QpExternSolver::NUM_BEST ){
//			 cout << "Objective: " << s.getObjFunctionValue().asDouble() << endl;
//			 cout << "Solution: " << endl;
//			 cout << "Solution-Size: " << s.getSolutionVector().size() << endl;
//			 for(int i = 0; i < s.getSolutionVector().size();i++){
//				 cout << i << " : " << s.getSolutionVector().at(i).asDouble()<< endl;
//			 }
//		 }
//	 }
//	 {
//		 algorithm::Algorithm::QlpSolution s = algorithm::Qlp2Lp( q.qlp ).solveQlp( algorithm::Algorithm::WORST_CASE );
//		 cout << "Org. Solution status: " << s.getSolutionStatusString() << endl;
//
//		 if( s.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL || s.getSolutionStatus() == extSol::QpExternSolver::NUM_BEST ){
//			 cout << "Org. Objective: " << s.getObjFunctionValue().asDouble() << endl;
//			 cout << "Org. Solution: " << endl;
//			 cout << "Org. Solution-Size: " << s.getSolutionVector().size() << endl;
//			 for(int i = 0; i < s.getSolutionVector().size();i++){
//					cout << i << " : " << s.getSolutionVector().at(i).asDouble()<< endl;
//			}
//		 }
//	 }
//	 cout << "Fertig" << endl;
//
//
//	return EXIT_SUCCESS;
//}*/
