/*
*
* Yasol: CutAdder.cpp -- Copyright (c) 2012-2014 Ulf Lorenz, Thomas Opfer
* Yasol: CutAdder.cpp -- Copyright (c) 2014-2017 Michael Hartisch, Ulf Lorenz, Thomas Opfer
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

#include "CutAdder.h"
//#define USE_NBD_CLP
//#undef USE_NBD_CPLEX_C

//#undef USE_NBD_CLP
//#define USE_NBD_CPLEX_C
#ifdef USE_NBD_HIGHS
#include "ExternSolvers/QpExtSolHighs.hpp"
#include "io/FilereaderMps.h"
#include "cmath"
#include <assert.h>
#endif
#ifdef USE_NBD_CLP
#include "ExternSolvers/QpExtSolCLP.hpp"
#include "ClpSimplex.hpp"
#include <CoinWarmStartBasis.hpp>
#include "CoinFactorization.hpp"
#include <ClpPrimalColumnSteepest.hpp>
#include <ClpDualRowSteepest.hpp>
#define LONG_MAX 0x7fffffffffffffffL
#endif
#ifdef USE_NBD_CPLEX_C
#include "ExternSolvers/QpExtSolCplexC.hpp"
#include "cmath"
#endif

using namespace std;

//#define TRACK_SOL
#define CCT_EPS 1e-10

const int BUG_FINDER=0;
const bool SHOW_COVER = (BUG_FINDER>=1?true:false);

//#define LIFTING

#ifdef LIFTING
#else
vector< pair< vector< pair<unsigned int, double> >, double > > CutAdder::getCoverCuts( extSol::QpExternSolver &extSolver, int *types, int8_t *assigns, int decLev, unsigned int initime, int * solu, int* fixs, int *blcks, int orgN ){
  //vector< pair< vector< pair<unsigned int, double> >, double > > CutAdder::getCoverCuts( CPXENVptr env, CPXCLPptr lp ){
		const double eps = 1e-9;
		time_t T0=time(NULL);
//cerr << "enter ";
		vector< pair< vector< pair<unsigned int, double> >, double > > cuts;

		if( extSolver.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL ||
				extSolver.getSolutionStatus() == extSol::QpExternSolver::NUM_BEST ||
				extSolver.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL_INFEAS){
	//		cout << endl << endl << "Beginn Tornister-Deckel-Reduzierungen" << endl << endl;

			static int fstM=-1;
			static int m_r = -1;
			if (fstM==-1 || decLev <= 1 || extSolver.getRowCount() < m_r) fstM = extSolver.getRowCount();
			unsigned int m = fstM;//extSolver.getRowCount();
			static bool f_t = false;
			if (decLev <= 1 || extSolver.getRowCount() < m_r) f_t=false;
			if (f_t == false) {
				f_t = true;
				m_r = m;
			}
			/*static*/ unsigned int n = extSolver.getVariableCount();
			//m = extSolver.getRowCount();
			if (decLev <= 1) {
				//n = extSolver.getVariableCount();
				m_r = extSolver.getRowCount();
			}
			if (m_r > extSolver.getRowCount()) m_r = extSolver.getRowCount();
			m = m_r;
            m = m_r = extSolver.getRowCount();
            
			cuts.reserve( 2*m+2 );

			//cerr << "enter:" << extSolver.getRowCount() << "," << m_r << "," << m << endl;

			{

				/*static*/ vector<double> xlpopt( n );
				/*static*/ vector<data::QpNum> colLower(n);
	            /*static*/ vector<data::QpNum> colUpper(n);
				/*static*/ vector<data::QpNum> objVals( n );

				extSolver.getValues( objVals );
				extSolver.getLB(colLower);
				extSolver.getUB(colUpper);
				for( unsigned int i = 0; i < n; ++i ){
					xlpopt[i] = objVals[i].asDouble();
				}

				/*static*/ vector<data::QpRhs> rhsVec( m_r );
				if (rhsVec.size() < m_r) rhsVec.resize(m_r);
#ifdef DEPRICATED_SLACKS
				///*static*/ vector<double> slacks(m_r);
				//if (slacks.size() < m_r) slacks.resize(m_r);
#ifdef USE_NBD_CPLEX_C
				int status = CPXXgetslack(
						*(CPXENVptr*)extSolver.getSolverEnv(),
						*(CPXLPptr*)extSolver.getSolverModel(),
						                slacks.data(), 0, m_r-1);
#elif USE_NBD_CLP
				//slacks.resize(this->getRowCount());
				const double *ra = ((ClpSimplex*)(extSolver.getSolverModel()))->getRowActivity();
				const double *lb = ((ClpSimplex*)(extSolver.getSolverModel()))->getRowLower();
				const double *ub = ((ClpSimplex*)(extSolver.getSolverModel()))->getRowUpper();
				for (int i = 0; i < m_r; ++i) {
				    slacks[i] = (lb[i] != -COIN_DBL_MAX ? lb[i] : ub[i])  - ra[i];
				}
#elif USE_NBD_IGHS
#endif
#endif
				extSolver.getRhs( rhsVec );
				extSolver.prepareMatrixRowForm();
				int cnt_succ=0;
				for (int ii=0;ii<1;ii++) {
				  int m = extSolver.getRowRhs_snapshot()->size();

					for( unsigned int i = 0; i < m; ++i ){
					  //if (i > m_r) break;
					  //if ((time(NULL)-initime)/10 < time(NULL) - T0) break;
						//if(cnt_succ > 1000+10*log2(n+m)) break;
						std::vector<data::IndexedElement> &rowtmp = *extSolver.getRowLhs_snapshot(i);
						double rhstmp = (*extSolver.getRowRhs_snapshot())[i].getValue().asDouble();//rhsVec.at( i ).getValue().asDouble();
						if (extSolver.getLazyStatus(i) == true || rowtmp.size() == 0) continue;
						if (extSolver.getStatus(i)) {
						  //cerr << "not C-"<< i << endl;
						  //for (int h = 0; h < rowtmp.size();h++) {
						  //  cerr << rowtmp[h].value.asDouble() << "x" << rowtmp[h].index << " + ";
						  //}
						  //cerr << " 0 <= " << rhstmp << endl; 
						  continue;
						}

						//extSolver.getRowLhs( i, rowtmp );
						if (rowtmp.size() == 0) continue;

						//if (rowtmp.size() > /*(double)n /*/ sqrt((double)n)) continue;

						double testRow = 0;

						for( unsigned int k = 0; k < rowtmp.size(); ++k ){
							testRow += xlpopt.at( rowtmp.at( k ).index ) * rowtmp.at( k ).value.asDouble();
						}

						if (fabs(testRow - rhstmp) > fabs(rhstmp)*(1e-9) + 1e-9) continue;

						//if( slacks[i] > abs(rhstmp)*(1e-9) + 1e-9) continue;
						//if (0&&abs(testRow - rhstmp) > rhstmp*(1e-1) + 1e-9 ) {
							//cerr << "+";
							//continue;
						//} //else cerr << "-";
						//if( rhsVec.at( i ).getRatioSign() == data::QpRhs::smallerThanOrEqual || rhsVec.at( i ).getRatioSign() == data::QpRhs::equal ){
						if( (*extSolver.getRowRhs_snapshot())[i].getRatioSign() == data::QpRhs::smallerThanOrEqual || (*extSolver.getRowRhs_snapshot())[i].getRatioSign() == data::QpRhs::equal ){
							vector<pair<double, unsigned int> > rowsparse;
							rowsparse.reserve( rowtmp.size() );

							for( unsigned int j = 0; j < rowtmp.size(); ++j ){
								rowsparse.push_back( make_pair( rowtmp.at( j ).value.asDouble(), rowtmp.at( j ).index ) );
							}

							vector<pair<unsigned int, double> > reslhs;
							double resrhs;
							if( createCover( xlpopt, rowsparse, rhstmp, CCT_EPS, reslhs, resrhs, types, assigns, decLev, colLower.data(), colUpper.data(), solu, fixs, blcks, orgN ) ){
								vector< pair<unsigned int, double> > cut;
								cut.reserve( n );

								// Ausgabetest ANFANG
								//cout << "Cover found." << endl;
								//cout << "ORG:";
								//for (int z=0; z < rowtmp.size();z++) {
								//	cout << rowtmp[z].value.asDouble() << "x" << rowtmp[z].index << "=" << xlpopt[rowtmp[z].index] << " + ";
								//}
								//cout << " <= " << rhstmp << endl;

								double testLP = 0;

								for( unsigned int k = 0; k < reslhs.size(); ++k ){
									testLP += xlpopt.at( reslhs.at( k ).first ) * reslhs.at( k ).second;
									cut.push_back( make_pair( reslhs.at( k ).first, -reslhs.at( k ).second ) );
									//cout << "x"<<reslhs.at( k ).first << "=" <<reslhs.at( k ).second << " + ";
								}
								//cout << " <= " << resrhs << endl;
								if( testLP > resrhs /*+ eps*/ && cut.size() > 0 ){
									//cout << "LP ok: " << testLP << " > " << resrhs << endl;
									cuts.push_back( make_pair( cut, -resrhs ) );
									cnt_succ++;
								} else {
									//cout << "!!!!! LP wird NICHT abgeschnitten: " << testLP << " <= " << resrhs << endl;
								}
								// Ausgabetest ENDE
							}
						}

						//if( rhsVec.at( i ).getRatioSign() == data::QpRhs::greaterThanOrEqual || rhsVec.at( i ).getRatioSign() == data::QpRhs::equal ){
						if( (*extSolver.getRowRhs_snapshot())[i].getRatioSign() == data::QpRhs::greaterThanOrEqual || (*extSolver.getRowRhs_snapshot())[i].getRatioSign() == data::QpRhs::equal ){
							// Koeffizienten/RHS negieren und selber Code wie oben.

							//double rhstmp = -rhsVec.at( i ).getValue().asDouble();
						        double rhstmp = -(*extSolver.getRowRhs_snapshot())[i].getValue().asDouble();
							vector<pair<double, unsigned int> > rowsparse;
							rowsparse.reserve( rowtmp.size() );

							for( unsigned int j = 0; j < rowtmp.size(); ++j ){
								rowsparse.push_back( make_pair( -rowtmp.at( j ).value.asDouble(), rowtmp.at( j ).index ) );
							}


							vector<pair<unsigned int, double> > reslhs;
							double resrhs;
							if( createCover( xlpopt, rowsparse, rhstmp, CCT_EPS, reslhs, resrhs, types, assigns, decLev, colLower.data(), colUpper.data(), solu, fixs, blcks, orgN ) ){
								vector< pair<unsigned int, double> > cut;
								cut.reserve( n );

								// Ausgabetest ANFANG
		//						cout << "Cover found." << endl;

								double testLP = 0;

								for( unsigned int k = 0; k < reslhs.size(); ++k ){
									testLP += xlpopt.at( reslhs.at( k ).first ) * reslhs.at( k ).second;
									cut.push_back( make_pair( reslhs.at( k ).first, -reslhs.at( k ).second ) );
								}

		//						cout << " <= " << coverRhs << endl;
								double tmp = -resrhs;
								if( testLP > resrhs /*+ eps*/ && cut.size() > 0 /*&&
										cleanCut( cut, cut.size(), tmp, objVals,
												colLower, colUpper, n, types) */) {
		//							cout << "LP ok: " << testLP << " > " << coverRhs << endl;
									cuts.push_back( make_pair( cut, tmp ) );
									cnt_succ++;
								} else {
		//							cout << "!!!!! LP wird NICHT abgeschnitten: " << testLP << " <= " << coverRhs << endl;
								}
								// Ausgabetest ENDE

							}
						}

					}
				}

			}


	//		cout << endl << endl << "Ende Tornister-Deckel-Reduzierungen" << endl << endl;


		} else {
			if (info_level >= 2) cerr << "LP status inappropriate for Cover cut generation." << extSolver.getSolutionStatus() << endl;
		}


		return cuts;

	}


#define OneEM12 1e-12
#define OneEM9 1e-9
#define OneEM7 1e-7
#define OneEM6 1e-6
#define OneEM5 1e-5
#define COVER_RHS_EPS (1e-9 * org_rowsparse.size() + 1e-7)

// expects sense = "<="
bool CutAdder::createCover( vector<double>& xlpopt, const vector<pair<double, unsigned int> >& org_rowsparse, const double org_rhs, const double eps, vector<pair<unsigned int, double> >& reslhs, double& resrhs, int *types, int8_t *assigns, int decLev, data::QpNum* colLower, data::QpNum* colUpper, int * solu, int *fixs, int *blcks , int orgN){
  //bool CutAdder::createCover( vector<double>& xlpopt, const vector<pair<double, unsigned int> >& rowsparse, double rhs, const double eps, vector<pair<unsigned int, double> >& reslhs, double& resrhs, int *types, int8_t *assigns, int decLev, data::QpNum* colLower, data::QpNum* colUpper ){
#ifndef FIND_BUG
  const unsigned int nbig = xlpopt.size();
  //bool CutAdder::createCover( const vector<double>& xlpopt, const vector<pair<double, unsigned int> >& rowsparse, double rhs, vector<unsigned int>& resPos, vector<unsigned int>& resNeg, const double eps ){

  const unsigned int n = xlpopt.size();

  double rhs = org_rhs + COVER_RHS_EPS;

  // Komplementieren

  static vector<double> xlpotmp(n);// = xlpopt;
  static vector<double> xlpotmp2(n);// = xlpopt;
  if (n > xlpotmp.size()) xlpotmp.resize(n);
  if (n > xlpotmp2.size()) xlpotmp2.resize(n);
  //for (int z=0;z<xlpopt.size();z++) xlpotmp.push_back(xlpopt[z]);
  //vector<bool> complemented( n, false );
  double rhstmp = rhs;
  bool containsReals=false;
  bool containsInts=false;
  bool printout=true;
  //if (org_rowsparse.size() > 10)
  if (BUG_FINDER == 0) 
     printout = false;
  bool containsNoFrac = true;
  vector<pair<double, unsigned int> > realPart;
  vector<pair<double, unsigned int> > rowsparse;
  int minblock=n+2;
  int maxblock=0;
  double orgLhs = 0.0;
		
  for( unsigned int j = 0; j < org_rowsparse.size(); ++j ){
    if (org_rowsparse[j].second >= orgN) return false;
    orgLhs = orgLhs + org_rowsparse[j].first * xlpopt[ org_rowsparse[j].second ];
    if (assigns[org_rowsparse.at( j ).second] == 2 /*extboolUndef*/) {
      if (blcks[org_rowsparse.at( j ).second] > maxblock) maxblock = blcks[org_rowsparse.at( j ).second];
      if (blcks[org_rowsparse.at( j ).second] < minblock) minblock = blcks[org_rowsparse.at( j ).second];
    }
    if (types[org_rowsparse.at( j ).second] != 0 /*not BINARY*/) {
      realPart.push_back(org_rowsparse[j]);
      containsReals=true;
      rowsparse.push_back(org_rowsparse[j]);
      /* ganz frueher: rechte Seite konservativ anpassen:
	 if (rowsparse.at( j ).first >= 0) {
	 rhstmp = rhstmp - colLower[rowsparse.at( j ).second].asDouble()*rowsparse.at( j ).first;
	 } else {
	 rhstmp = rhstmp - colUpper[rowsparse.at( j ).second].asDouble()*rowsparse.at( j ).first;
	 }
	 xlpotmp[rowsparse.at( j ).second] = xlpopt[rowsparse.at( j ).second];
	 xlpopt[rowsparse.at( j ).second] = 0.0;
      */
    } else {
      containsInts=true;
      rowsparse.push_back(org_rowsparse[j]);
    }
    if (types[org_rowsparse.at( j ).second]==0 && xlpopt.at( org_rowsparse.at( j ).second ) > OneEM12 && xlpopt.at( org_rowsparse.at( j ).second ) < 1.0-OneEM12) {
      containsNoFrac = false;
    }
  }
        
  //if (orgLhs < rhs - 0.0001 - rhs * 0.01) return false;
  if (rowsparse.size() <= ((containsReals && containsInts) ? 2 : 1) || (containsReals && !containsInts) || (containsNoFrac)) {
    //cerr << "ups, no more a constraint in covercut generation" << endl;
    return false;
  }

  if (printout /*&& containsReals && containsInts*/) {
    double solusum_double=0.0;
    cerr << "try to get cut from row ";
    for( unsigned int j = 0; j < rowsparse.size(); ++j ){
      solusum_double = solusum_double + solu[rowsparse[j].second]*rowsparse[j].first;
      cerr << rowsparse.at( j ).first << (types[rowsparse.at( j ).second] == 0?"x":"y") << rowsparse.at( j ).second << "(" << xlpopt[rowsparse.at( j ).second] << ")";
      if (types[rowsparse.at( j ).second] != 0) {
	cerr << "[" << colLower[rowsparse.at( j ).second].asDouble() << "," << colUpper[rowsparse.at( j ).second].asDouble() << "]";
      } else {
	cerr << "[" << solu[rowsparse.at( j ).second] << "," << xlpopt[rowsparse.at( j ).second] << "]";
      }

      cerr << " + ";
    }
        
    cerr << " 0 <= " << rhs << " >= solusum=" << solusum_double  << " len=" << org_rowsparse.size() << endl;
  }
        
  for( unsigned int j = 0; j < rowsparse.size(); ++j ){
    if (rowsparse[j].first > -OneEM9/*12*/ && rowsparse[j].first < OneEM9/*12*/) {
      cerr << "Warning: Cover coefficient < OneEM12." << endl;
      //return false;
    }
    assert(rowsparse[j].first <= -OneEM12 || rowsparse[j].first >= OneEM12);
    //assert(types[rowsparse[j].second] == 0);
  }

  //complement variables
  double shift=0.0;
  int countOrgInts=0;
  for( unsigned int j = 0; j < rowsparse.size(); ++j ){
    xlpotmp2[rowsparse.at( j ).second] = xlpopt[rowsparse.at( j ).second];
    if (types[rowsparse.at( j ).second] == 0 /*BINARY*/)
      countOrgInts++;
    if( rowsparse.at( j ).first < 0) {
      if (types[rowsparse.at( j ).second] == 0 /*BINARY*/){
	unsigned int ind = rowsparse.at( j ).second;
	xlpopt.at( ind ) = 1 - xlpopt.at( ind );
	//complemented.at( ind ) = true;
	rhstmp -= rowsparse.at( j ).first;
      } else {
	unsigned int ind = rowsparse.at( j ).second;
	shift = 	-colLower[rowsparse.at( j ).second].asDouble();;
	xlpopt.at( ind ) = xlpopt.at( ind ) + shift;
	rhstmp += shift*rowsparse.at( j ).first;
      }
    } else if (types[rowsparse.at( j ).second] != 0 /*not BINARY*/) {
      unsigned int ind = rowsparse.at( j ).second;
      shift = 	colUpper[rowsparse.at( j ).second].asDouble();
      rowsparse.at( j ).first = -rowsparse.at( j ).first;
      xlpopt.at( ind ) = -xlpopt.at( ind ) + shift;
      rhstmp += shift*(-rowsparse.at( j ).first);
    }
    xlpotmp[rowsparse.at( j ).second] = xlpopt[rowsparse.at( j ).second];
  }

  if (printout /*&& containsReals && containsInts*/) {
    cerr << "try to get cut from modified row ";
    for( unsigned int j = 0; j < rowsparse.size(); ++j ){
      cerr << rowsparse.at( j ).first << (types[rowsparse.at( j ).second] == 0?"x":"y") << rowsparse.at( j ).second << "(" << xlpopt[rowsparse.at( j ).second] << ")";
      if (types[rowsparse.at( j ).second] != 0) {
	cerr << "[" << colLower[rowsparse.at( j ).second].asDouble()+shift << "," << colUpper[rowsparse.at( j ).second].asDouble()+shift << "]";
      }
      cerr << " + ";
    }
            
    cerr << " 0 <= " << rhstmp << " len=" << org_rowsparse.size() << endl;
  }

  // Cover berechnen

  /*static*/ vector<pair<double, unsigned int> > result;
  /*static*/ vector<pair<double, unsigned int> > S;
  /*static*/ vector<pair<double, unsigned int> > remS;
  /*static*/ vector<pair<double, unsigned int> > K_minus_S;

  //if (result.capacity() < rowsparse.size() ) result.reserve( rowsparse.size() );
  result.clear();
  S.clear();
  remS.clear();
  K_minus_S.clear();
  static vector<double> rowdense( n );
  if (n > rowdense.size()) rowdense.resize(n);

  for( unsigned int i = 0; i < rowsparse.size(); ++i ){

    const unsigned int ind = rowsparse.at( i ).second;
    if (types[ind] != 0 /*not binary*/) continue;

    rowdense.at( ind ) = rowsparse.at( i ).first;
    if( xlpopt.at( ind ) > 0.0 ){
      result.push_back( make_pair( xlpopt[ind] /*rowsparse[i].first*/, ind ) );
    } else {
      K_minus_S.push_back( make_pair( fabs(rowdense[ind]), ind) ); // oder abs?
    }
  }

  sort( result.begin(), result.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return p1.first > p2.first; } );

  /*static*/ vector<unsigned int> res;
  res.clear();
  //for (int u = 0; u < realPart.size();u++) {
  //    assert(xlpopt[realPart[u].second] >= 0);
  //    rhstmp = rhstmp + fabs(realPart[u].first) * xlpopt[realPart[u].second];
  //}
        
  double sum = 0;
  double sumSolu=0.0;
  bool contEx=false;
  for( unsigned int j = 0; j < result.size(); ++j ){
    assert(types[result[j].second] == 0);
    if (types[result[j].second] != 0 /*BINARY */) {
      contEx=true;
      continue;
    }
    sum += fabs(rowdense.at( result[j].second )); // oder abs?
#ifdef TRACK_SOL
    if (solu[result[j].second] == 1 && rowdense.at( result[j].second ) > 0)
      sumSolu += fabs(rowdense.at( result[j].second ));
    else if (solu[result[j].second] == 0 && rowdense.at( result[j].second ) < 0)
      sumSolu += fabs(rowdense.at( result[j].second ));
#endif
    if( sum > rhstmp +OneEM12 + fabs(rhstmp)*(OneEM12) ){ // oder abs?
      //cerr << abs(rhstmp)*(1.0/(OneEM12)) << endl;
      //if (res.capacity() < j) res.reserve( j );
      for( unsigned int k = 0; k <= j; ++k ){
	if (types[result[k].second] != 0 /*BINARY */) continue;
	res.push_back( result.at( k ).second );
      }
      for( unsigned int k = j+1; k < result.size(); k++ ){
	if (types[result[j].second] != 0 /*BINARY */) continue;
	assert(types[result[j].second] == 0);
	if (assigns[result[j].second] == 2 /*&& xlpopt[result[j].second] > 0.00001*/) {
	  K_minus_S.push_back( make_pair( fabs(rowdense[result.at( k ).second]), result.at( k ).second )); // oder abs?
	}
      }
      res.clear();
      break;
    }
  }
  if (printout) cerr << "CONT EXISTS:" << contEx << endl;
  if (printout /*&& containsReals && containsInts*/) {
    cerr << "sum=" << sum << " und rhs=" << rhstmp << " und eps=" << eps << endl;
    cerr << "sumSolu=" << sumSolu << " und rhs=" << rhstmp << " und eps=" << eps << endl;
    for (int j = 0; j < res.size(); j++) {
      cerr << " x" << res[j] << "(" << rowdense[res[j]] << "," << solu[res[j]]  << ")";
    }
    cerr << endl;
  }
        
  // Schritt 3 im Paper
  if(res.size() > 0){
    bool removed = false;
    vector<bool> remove( res.size(), false );
    for( int j = res.size() - 1; j >= 0; --j ){

      // TODO macht das wirklich Sinn??
      if( 0&&xlpopt.at( res.at( j ) ) >= 1 - 1e-20 ){
	break;
      }

      if( sum - fabs(rowdense.at( res.at( j ) ) ) > rhstmp +OneEM7/*12*/ + abs(rhstmp)*(OneEM9/*12*/)){ // oder abs?
	if (printout) cerr << "rausnehmen von x" << res[j] << " mit Gewicht " << rowdense[res[j]] << " Rest:" << (sum - fabs(rowdense.at( res.at( j ) ) ))-rhstmp << endl;
	removed = true;
	remove.at( j ) = true;
	sum -= fabs(rowdense.at( res.at( j ) ) ); // oder abs?
      }
    }

    if( removed ){
      // Etwas overkill, dafür verstehbar
      vector<unsigned int> oldRes = res;
      //cerr << "A:o.size()=" << oldRes.size() << " res.size=" << res.size() << endl;
      res.clear();
      //cerr << "B:o.size()=" << oldRes.size() << " res.size=" << res.size() << endl;
      res.reserve( oldRes.size() );
      //cerr << "C:o.size()=" << oldRes.size() << " res.size=" << res.size() << endl;
      for( unsigned int i = 0; i < oldRes.size(); ++i ){
	if( !remove.at( i ) ){
	  res.push_back( oldRes.at( i ) );
	} else {
	  K_minus_S.push_back( make_pair( fabs(rowdense[oldRes.at( i )]), oldRes.at( i ) )); // oder abs?
	}
      }
      /*
	cerr << "after DEL sum=" << sum << " und rhs=" << rhstmp << " und eps=" << eps << endl;
	for (int j = 0; j < res.size(); j++) {
	cerr << " x" << res[j];
	}
	cerr << endl;
      */
    }
  }
  //cerr << "rt=" << rhstmp;

  if( !res.size() || res.size() == countOrgInts){
    for( unsigned int j = 0; j < rowsparse.size(); ++j ){
      if( org_rowsparse.at( j ).first < 0) {
	if (types[rowsparse.at( j ).second] == 0 /*BINARY*/){
	  unsigned int ind = rowsparse.at( j ).second;
	  xlpopt.at( ind ) = 1 - xlpopt.at( ind );
	  //complemented.at( ind ) = true;
	} else {
	  unsigned int ind = rowsparse.at( j ).second;
	  shift = 	-colLower[rowsparse.at( j ).second].asDouble();;
	  xlpopt.at( ind ) = xlpopt.at( ind ) - shift;
	}
      } else if (types[rowsparse.at( j ).second] != 0 /*not BINARY*/) {
	unsigned int ind = rowsparse.at( j ).second;
	shift = 	colUpper[rowsparse.at( j ).second].asDouble();
	xlpopt.at( ind ) = xlpopt.at( ind ) - shift;
	xlpopt.at( ind ) = -xlpopt.at( ind );
	// is same as: xlpopt.at( ind ) = -xlpopt.at( ind ) + shift;
	rowsparse.at( j ).first = -rowsparse.at( j ).first;
      }
    }
    return false;
  }

  if (printout && containsReals && containsInts) {
    cerr << "Minimum Cuts:" << endl;
    for (int j = 0; j < res.size(); j++) {
      cerr << " x" << res[j] << "(" << rowdense[res[j]] << ")";
    }
    cerr << endl;
  }
#ifdef TRACK_SOL
  cerr << "C1:" << endl;
  int solusum=0.0;
  for (unsigned int i = 0; i < res.size(); ++i ){
    cerr << " x" << res[i] << "(" << (solu[res[i]]) << ")";
    solusum++;
  }
  cerr << " + 0 <= " << res.size() -1 << " cf." << solusum << endl;
#endif
  reslhs.clear();
  resrhs = (double)res.size()-1.0;
  for (int i = 0; i < res.size();i++) {
    S.push_back( make_pair( 1.0, res.at( i ) ));
    remS.push_back( make_pair( 1.0, res.at( i ) ));
  }
  {
    double testLP = 0.0;
    double sumSolu=0.0;
#ifdef TRACK_SOL
    for( unsigned int k = 0; k < S.size(); ++k ){
      testLP += xlpopt.at( S.at( k ).second ) * S.at( k ).first;
      if (solu[S.at( k ).second] == 1 && rowdense.at( S.at( k ).second ) > 0)
	sumSolu += S[k].first;//fabs(rowdense.at( res[j] ));
      else if (solu[S.at( k ).second] == 0 && rowdense.at( S.at( k ).second ) < 0)
	sumSolu += S[k].first;//fabs(rowdense.at( S.at( k ).second ));
    }

    if (1||sumSolu > resrhs) {
      cerr << "?";
      cerr << "?" << sumSolu << "<=" << resrhs <<";" << endl;
    }
#endif
  }

  if (SHOW_COVER) {
		  
    cerr << "rowsparse und rhs=" << resrhs << " und eps=" << eps << endl;
    for (int j = 0; j < org_rowsparse.size(); j++) {
      cerr << org_rowsparse[j].first << "x" << org_rowsparse[j].second << " + ";
    }
    cerr << " <= "<< rhs << endl;
    cerr << "COVER:";
    for (int j = 0; j < S.size(); j++) {
      cerr << S[j].first << "x" << S[j].second << " + ";
    }
    cerr << " <= "<< resrhs << endl;
    cerr << " END "<< endl;
  }

  // lifting
  // Cover ist in S, Koeffizientenmenge K ist im Original rowsparse; K-S ist in K_minus_S
  // sortiere ggfs. K-S, kleine Koeffizienten zuerst. Code funktioniert nur, wenn so sortiert wird. Warum?

  double ttestLP = 0.0;
  for( unsigned int k = 0; k < S.size(); ++k ){
    ttestLP += xlpopt.at( S.at( k ).second ) * S.at( k ).first;
  }
  //if (testLP > resrhs) cerr << "+";
  //if (containsReals) cerr << "+++";
  bool printcut = true;//false;
  if (!containsReals /*&& ttestLP <= resrhs*/) {
#define USE_LIFTING_I
    if (decLev < (int)sqrt((double)n)) {

      sort( K_minus_S.begin(), K_minus_S.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return p1.first*xlpotmp[p1.second] > p2.first*xlpotmp[p2.second]; } );
      /*
	for (int j = 0; j < K_minus_S.size(); j++) {
	cerr << K_minus_S[j].first << "x" << K_minus_S[j].second << " + ";
	}
	cerr << " ascend order check "<< endl;
      */

      // for (int k = |K-S-1|; k >= 0; k--)

      sort( S.begin(), S.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return p1.first / fabs(rowdense[p1.second]) > p2.first / fabs(rowdense[p2.second]); } ); // oder abs?
      int remSz = S.size();
#ifdef USE_LIFTING_I
      for (int k = 0; ttestLP <= resrhs && k < K_minus_S.size();k++) {
	if (assigns[K_minus_S[k].second] != 2) continue;
	if (xlpopt[K_minus_S[k].second] > 1.0 - OneEM9) continue;
	if (xlpopt[K_minus_S[k].second] <= OneEM9) continue;
	//for (int o=0;o<S.size();o++)
	//assert(K_minus_S[k].second != S[o].second);
	//   sortiere S gem�� Quotient fk/|aik|, gro�e zuerst, fk ist der Koeffizient in S, aik der im Original
	//for (int o=0;o<S.size()-1;o++) {
	//cerr << S[o].first / abs(rowdense[S[o].second]) << "|" << S[o+1].first / abs(rowdense[S[o+1].second]) << " ";
	//assert(S[o].first / abs(rowdense[S[o].second]) >= S[o+1].first / abs(rowdense[S[o+1].second]) );
	//}
	//cerr << endl;
	//sort( S.begin(), S.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return p1.first / abs(rowdense[p1.second]) > p2.first / abs(rowdense[p2.second]); } );
	/*
	  for (int j = 0; j < S.size(); j++) {
	  cerr << S[j].first << "x" << S[j].second << ":" << S[j].first / abs(rowdense[S[j].second])<< ":" <<  S[j].first << " + ";
	  }
	  cerr << " ascend order check "<< endl;
	*/

	//   a0 = cut.rhs

	double a0 = rhstmp;

	//   ak = Koeffizient von koeff[K-S [k]]

	double ak = fabs(rowdense[K_minus_S[k].second]); // oder abs?

	//   rhs = a0-ak

	double lrhs = a0-ak;

	//   volle_anteile = 0

	double complete = 0.0;

	//   for (int r = 0; r < S;r++)
	//      if (koeff[S[r]] <= rhs)
	//          rhs -= koeff[S[r]]
	//          volle_anteile ++
	//      else
	//          anteil = rhs / koeff[S[r]]
	//          fk = |S|-1-volle_anteile - floor(anteil)
	//          if (fk>0) fuege fk*xk zu S hinzu und nimm xk aus K-S raus
	//          break

	int rs=S.size();
	for (int r = 0; r < rs/*S.size()*/;r++) {
	  assert(S[r].first < -OneEM12 || S[r].first > OneEM12);
	  if (fabs(rowdense[S[r].second]) <= lrhs) { // oder abs?
	    lrhs = lrhs - fabs(rowdense[S[r].second]); // oder abs?
	    complete = complete + S[r].first;//1.0;
	  } else {
	    double part = lrhs;// / fabs(rowdense[S[r].second]); // oder abs?
	    int fk = resrhs - complete;// - floor(part);//- floor(part*S[r].first);
	    if (fk > 0+eps) {
	      std::pair<double, unsigned int> P( fk, K_minus_S[k].second );
	      if (0&&P.first / abs(rowdense[P.second]) <=
		  S[S.size()-1].first / abs(rowdense[S[S.size()-1].second] - 0*OneEM12 - 0*(OneEM12)*abs(S[S.size()-1].first / abs(rowdense[S[S.size()-1].second])))) {
		S.push_back(make_pair( (double)fk, K_minus_S[k].second ));
	      } else {
		assert(S.size() > 0);
		int ix = S.size()-1;
		S.push_back(make_pair(-1.0,0));
		if (P.first / fabs(rowdense[P.second]) <= S[ix].first / fabs(rowdense[S[ix].second])) {
		  S[ix+1].first = P.first;
		  S[ix+1].second = P.second;
#ifdef TRACK_SOL
		  cerr << "add " << fk << "x" << P.second << endl;
#endif
		} else {
		  while (ix > 0 && P.first / fabs(rowdense[P.second]) > S[ix].first / fabs(rowdense[S[ix].second])) { // oder abs?
		    //cerr << P.first / abs(rowdense[P.second]) << "*" << S[ix].first / abs(rowdense[S[ix].second]) << "*" << ix << " ";
		    S[ix+1].first = S[ix].first;
		    S[ix+1].second = S[ix].second;
		    ix--;
		    if (P.first / fabs(rowdense[P.second]) <= S[ix].first / fabs(rowdense[S[ix].second])) {  // oder abs?
		      ix++;
		      break;
		    }
		  }
		  if (ix == 0) {
		    if (S[ix].first < 0.0 || P.first / fabs(rowdense[P.second]) <= S[ix].first / fabs(rowdense[S[ix].second])) {
		      S[ix+1].first = P.first;
		      S[ix+1].second = P.second;
		    } else {
		      S[ix].first = P.first;
		      S[ix].second = P.second;
		    }
		  } else {
		    S[ix].first = P.first;
		    S[ix].second = P.second;
		  }
		  //cerr << "B:ix=" << ix << " P=" << P.first << "x"<< P.second << " resrhs="<< resrhs<< endl;
		  //cerr << "S[ix].first=" << S[ix].first << " P.first / fabs(rowdense[P.second]) <= S[ix].first / fabs(rowdense[S[ix].second])=" << (int)(P.first / fabs(rowdense[P.second]) <= S[ix].first / fabs(rowdense[S[ix].second])) << " S[ix+1].first=" << S[ix+1].first << endl;
		}
		assert(S[S.size()-1].first < -OneEM12 || S[S.size()-1].first > OneEM12);
#ifdef TRACK_SOL
		printcut = true;
#endif
	      }
	      /*
		for (int j = 0; j < S.size(); j++) {
		cerr << S[j].first << "x" << S[j].second << ":" << S[j].first / abs(rowdense[S[j].second])<< ":" <<  S[j].first << " + ";
		}
		cerr << " added last " << part << "," << complete << "," << rhstmp<< endl;
	      */
	      //resrhs = S.size()-1;
	    }
	    break;
	  }
	}
      }
#endif
      double test_binwert_min=0.0;
      for( unsigned int kk = 0; kk < resrhs+1; ++kk ){
	test_binwert_min = test_binwert_min + S[kk].first*fabs(rowdense[S[kk].second]);
      }
      if (test_binwert_min <= rhstmp) do {
	  sort( S.begin(), S.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return fabs(rowdense[p1.second])*p1.first < p2.first*fabs(rowdense[p2.second]); } );
	  test_binwert_min=0.0;
	  for( unsigned int kk = 0; kk < resrhs+1; ++kk ){
	    test_binwert_min = test_binwert_min + S[kk].first*fabs(rowdense[S[kk].second]);
	  }
	  if (test_binwert_min <= rhstmp) {
	    S[0] = S[S.size()-1];
	    S.pop_back();
	    if (info_level >= 2) cerr << " sBy1;";
	  }
	} while (test_binwert_min <= rhstmp  && S.size() > resrhs+1);
      bool isDoubleContained=false;
      for( unsigned int k = 0; k < S.size()-1; ++k ){
	if (S.at( k ).second == S.at( k+1 ).second) {
	  if (info_level >= 2) cerr << "Warning: cover cut lifting" << endl;
	  //isDoubleContained = true;
	  //break;
	  if (k+1 == S.size()-1) S.pop_back();
	  else {
	    for (int kk = k+1; kk < S.size()-1;kk++) {
	      S[kk] = S[kk+1];
	    }
	    S.pop_back();
	  }
	}
      }
      for( unsigned int k = 0; k < S.size()-1; ++k ){
	if (S.at( k ).second == S.at( k+1 ).second) {
	  if (info_level >= 2) cerr << "Error: cover cut lifting" << endl;
	  isDoubleContained = true;
	  break;
	}
      }
      if (isDoubleContained || test_binwert_min <= rhstmp || S.size() < 2) {
	S.clear();
	res.clear();
	for (int o = 0; o < remS.size();o++) {
	  S.push_back(remS[o]);
	  res.push_back(remS[o].second);
	}
	resrhs = (double)res.size()-1.0;
	if (info_level >= 2) cerr << "lost lifting." << endl;
#ifdef TRACK_SOL
	cerr << "B binwert =" << test_binwert_min << " >? rhstmp =" << rhstmp << endl;
#endif
      }
      double testLP = 0.0;
      double sumSolu=0.0;
#ifdef TRACK_SOL
      for( unsigned int k = 0; k < S.size(); ++k ){
	testLP += xlpopt.at( S.at( k ).second ) * S.at( k ).first;
	if (solu[S.at( k ).second] == 1 && rowdense.at( S.at( k ).second ) > 0)
	  sumSolu += S[k].first;//fabs(rowdense.at( res[j] ));
	else if (solu[S.at( k ).second] == 0 && rowdense.at( S.at( k ).second ) < 0)
	  sumSolu += S[k].first;//fabs(rowdense.at( S.at( k ).second ));
      }

      if (testLP > resrhs) {
	cerr << "*";
	printcut = true;
	cerr << "*" << sumSolu << "<=" << resrhs <<";";
      }
#endif
      testLP = 0.0;
      for( unsigned int k = 0; k < S.size(); ++k ){
	testLP += xlpopt.at( S.at( k ).second ) * S.at( k ).first;
      }
      /*else cerr << "-";*/
      static int thecount=0;
      if (SHOW_COVER) {
	cerr << "SECOND rowsparse und rhs=" << resrhs << " und eps=" << eps << endl;
	for (int j = 0; j < org_rowsparse.size(); j++) {
	  cerr << org_rowsparse[j].first << "x" << org_rowsparse[j].second << " + ";
	}
	cerr << " <= "<< rhs << endl;
	cerr << "END-COVER:";
	for (int j = 0; j < S.size(); j++) {
	  cerr << S[j].first << "x" << S[j].second << " + ";
	}
	cerr << " <= "<< resrhs << endl;
	cerr << " END "<< endl;
      }


      if (minblock == maxblock && /*thecount < 5000 &&*/ testLP <= resrhs - OneEM12 - fabs(resrhs)*OneEM12 && decLev <= (int)sqrt((double)n)) {
	// look for a 1,k-cut 1K cut 1k cut 1-k cut
	// old minimal cover still in res
	// find max |aij| of S
	S.clear();
	res.clear();
	for (int o = 0; o < remS.size();o++) {
	  S.push_back(remS[o]);
	  res.push_back(remS[o].second);
	}
	resrhs = (double)res.size()-1.0;

	double max_a = -1.0;
	int max_a_index = -1;
	for (int o = 0; o < res.size();o++)
	  if (fabs(rowdense[res[o]]) > max_a) {
	    max_a = fabs(rowdense[res[o]]);
	    max_a_index = o;
	  }
	double max2_a = -1.0;
	int max2_a_index = -1;
	for (int o = 0; o < res.size();o++) {
	  if (fabs(rowdense[res[o]]) > max2_a && o != max_a_index) {
	    max2_a = fabs(rowdense[res[o]]);
	    max2_a_index = o;
	  }
	}
	if (max2_a != max_a) {
	  int t = res[max_a_index];
	  // K_minus_S m�sste noch ok sein
	  /*static*/ vector<pair<double, unsigned int> > Sstar;
	  Sstar.clear();
	  for (int o = 0; o < res.size();o++) {
	    if (fabs(rowdense[res[o]]) < max_a) {
	      Sstar.push_back(make_pair(1.0, res[o]));
	    }
	  }
	  double c1=0.0;
	  bool tin=false;
	  for (int i = 0; i < res.size(); i++) {
	    c1 = c1 + fabs(rowdense[res[i]]); // oder abs?
	    if (res[i] == t) tin = true;
	  }
	  //rhstmp = rhstmp + OneEM12 + fabs(rhstmp)*(OneEM12); // oder abs?
	  assert(tin==true);
	  if (info_level >= 2) if (c1 <= rhstmp) cerr << "Error: c1=" << c1 << " <= rhstmp=" << rhstmp << endl;
	  if (c1 <= rhstmp || c1-fabs(rowdense[t]) > rhstmp) {
	    cerr << "Warning in cover cut generation. c1:" << c1 << " rhstmp=" << rhstmp << " c1-fabs(row...):" << c1-fabs(rowdense[t]) << endl;
	    rhstmp = rhstmp + OneEM12 + fabs(rhstmp)*(OneEM12);
	    if (c1 <= rhstmp || c1-fabs(rowdense[t]) > rhstmp) {
	      cerr << "Error in cover cut generation. c1:" << c1 << " rhstmp=" << rhstmp << " c1-fabs(row...):" << c1-fabs(rowdense[t]) << endl;
	      return false;
	    }
	    //return false;
	  }
	  assert(c1 > rhstmp);
	  assert(c1-fabs(rowdense[t]) <= rhstmp); // oder abs?
	  for (int i = 0; i < res.size(); i++) {
	    if (c1-fabs(rowdense[res[i]]) > rhstmp && info_level >= 2) cerr << "c1=" << c1-fabs(rowdense[res[i]]) << " rhstmp=" << rhstmp << endl; // oder abs?
	    rhstmp = rhstmp + OneEM12 + fabs(rhstmp)*(OneEM12); // oder abs?
	    if (c1-fabs(rowdense[res[i]]) > rhstmp) {
	      cerr << "Error in cover cut generation. c1:" << c1 << " rhstmp=" << rhstmp << " c1-fabs(row...):" << c1-fabs(rowdense[res[i]]) << endl;
	      return false;
	    } 
	    assert(c1-fabs(rowdense[res[i]]) <= rhstmp); // oder abs?
	  }

	  bool found1k=false;
	  if (Sstar.size() == res.size()-1) {
	    //cerr << "O";
	    //else cerr << "(" << max_a << "," << Sstar.size() << "," << res.size() << ")";
	    double k = res.size()-1;
	    double check_aij=0.0;
	    double check_ws=0.0;
	    double a0 = rhstmp;
	    //cerr << "a0=" << a0;
	    for (int o = 0; o < Sstar.size();o++) {
	      check_aij = check_aij + fabs(rowdense[Sstar[o].second]); // oder abs?
	      check_ws = check_ws + xlpopt.at( Sstar.at( o ).second );
	    }
	    //cerr << "::check_aij=" << check_aij << endl;
	    //cerr << "::max2_a=" << max2_a << endl;
	    assert(check_aij <= a0);
	    //cerr << "caij=" << check_aij << " a0=" << a0 << " t=" << check_aij + abs(rowdense[t]) << endl;
	    //cerr << "cws=" << check_ws << " k=" << k << " check_es+t=" << check_ws+xlpopt[t]<< endl;
	    for (int o = 0; o < K_minus_S.size();o++) {
	      if (xlpopt[K_minus_S.at( o ).second] < OneEM9) continue;
	      if (fabs(rowdense[K_minus_S[o].second]) > max2_a) continue; // oder abs?
	      if (check_aij + fabs(rowdense[K_minus_S[o].second]) <= a0 - OneEM9 - (OneEM9)*fabs(a0)) { // oder abs?
		//check dass K_minus_S[o].second nicht in Sstar
		for (int l=0; l < Sstar.size();l++) {
		  assert(Sstar[l].second != K_minus_S[o].second);
		}
		// Sstar = Sstar union {K_minus_S[o].second]}
		Sstar.push_back(make_pair(1.0, K_minus_S[o].second));
		double r = Sstar.size();
		assert(check_aij + fabs(rowdense[K_minus_S[o].second]) <= rhstmp); // oder abs?
		// (r-k+1)*x_t + Sum_(j in T(r)) x_j <= r
		// T(r) : Summe aller Teilmengen von S* der Gr��e r
		if ((r-k+1)*xlpopt[t] + (check_ws+xlpopt.at( K_minus_S.at( o ).second )) > r-OneEM9) {
#ifdef TRACK_SOL
		  cerr << "F:" << (Sstar.size()-k+1)*xlpopt[t] + (check_ws+xlpopt.at( K_minus_S.at( o ).second )) << " > " << Sstar.size() << endl;
#endif
		  //cerr << "G:" <<
		  /*
		    cerr << "show OLD cut: ";
		    for (int jj=0; jj < S.size();jj++) {
		    cerr << S[jj].first << "x" << S[jj].second << " + ";
		    sum = sum + abs(rowdense[S[jj].second]);
		    }
		    cerr << " 0 <= " << resrhs << endl;
		    cerr << "show OLD cover: ";
		    double rr=0.0;
		    for (int jj=0; jj < S.size();jj++) {
		    rr = rr + abs(rowdense[S[jj].second]);
		    cerr << abs(rowdense[S[jj].second]) <<  " + ";
		    }
		    cerr << " 0 = " << rr << " <=? " << rhstmp << endl;
		  */
		  S.clear();
		  for (int l=0;l < Sstar.size();l++) {
		    S.push_back(make_pair(Sstar[l].first,Sstar[l].second));
		  }
		  sort( S.begin(), S.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return fabs(rowdense[p1.second])*p1.first < p2.first*fabs(rowdense[p2.second]); } );
		  double test_binwert_min=0.0;
		  double test_binwert_max=0.0;
		  for( unsigned int kk = 0; kk < k/*r*//*Sstar.size()*/; ++kk ){
		    test_binwert_min = test_binwert_min + fabs(rowdense[S[kk].second]);
		  }
		  for( unsigned int kk = 0; kk < k/*r*//*Sstar.size()*/; ++kk ){
		    int kkk = S.size()-1-kk;
		    test_binwert_max = test_binwert_max + fabs(rowdense[S[kkk].second]);
		  }
		  if (0&&test_binwert_max > a0) {
#ifdef TRACK_SOL
		    cerr << "binwert=" << test_binwert_max << " > a0 =" << a0 << endl;
#endif
		    break;
		  }
		  if (test_binwert_min + fabs(rowdense[t]) <= a0) {
#ifdef TRACK_SOL
		    cerr << "binwert+{t}=" << test_binwert_min+fabs(rowdense[t]) << " <= a0 =" << a0 << endl;
#endif
		    break;
		  }
		  found1k = true;
		  S.push_back(make_pair((double)(r-k+1.0),t));
		  assert(r-k+1.0 > OneEM12 || r-k+1.0 < -OneEM12);
		  resrhs = r;//Sstar.size();
#ifdef TRACK_SOL
		  double sum=0.0;
		  cerr << "show cover: ";
		  for (int jj=0; jj < remS.size();jj++) {
		    if (rowdense[remS[jj].second] >= 0)
		      cerr << remS[jj].first << "x" << remS[jj].second <<"(" << solu[remS[jj].second] << ") + ";
		    else
		      cerr << remS[jj].first << "x" << remS[jj].second <<"(" << 1.0-solu[remS[jj].second] << ") + ";
		    sum = sum + fabs(rowdense[remS[jj].second]); // oder abs?
		  }
		  cerr << " 0 <= " << remS.size()-1.0 << endl;
		  cerr << "show 1k sum: ";
		  double rr=0.0;
		  double xisum = 0.0;
		  for (int jj=0; jj < S.size();jj++) {
		    rr = rr + abs(rowdense[S[jj].second]);
		    xisum = xisum + xlpopt.at( S.at( jj ).second ) * S.at( jj ).first;
		    cerr << abs(rowdense[S[jj].second]) <<  " + ";
		  }
		  cerr << " 0 = " << rr << " <=? " << rhstmp << endl;
		  double binwert=0.0;
		  double flowwert = 0.0;
		  for( unsigned int k = 0; k < S.size(); ++k ){
		    flowwert = flowwert + xlpopt[ S.at( k ).second ] * S.at( k ).first;
		    if (rowdense[S.at( k ).second] >= 0) {
		      binwert = binwert + solu[ S.at( k ).second ] * S.at( k ).first;
		    } else {
		      binwert = binwert + (1.0 - solu[ S.at( k ).second ]) * S.at( k ).first;
		    }
		  }
		  cerr << "show cut: ";
		  for (int jj=0; jj < S.size();jj++) {
		    if (rowdense[S[jj].second] >= 0)
		      cerr << S[jj].first << "x" << S[jj].second << "(" << solu[S[jj].second]<< "," << xlpopt[S[jj].second] << ")" << " + ";
		    else
		      cerr << S[jj].first << "x" << S[jj].second << "(" << 1.0-solu[S[jj].second]<< "," << xlpopt[S[jj].second] << ")" << " + ";
		  }
		  cerr << " 0 <=? " << r << endl;
		  cerr << "binwert=" << binwert << " und flowwert=" << flowwert << endl;

		  if (rr > rhstmp) {
		    cerr << "aijsumr=" << rr << " > " << "rhs=" << rhstmp << endl;
		    cerr << "xijsum=" << xisum << " > " << "rhs=" << r << endl;
		    cerr << "check+aj+t=" << check_aij + abs(rowdense[K_minus_S[o].second])+abs(rowdense[t]) << " >? " << a0 << endl;
		  }
#endif
		  //thecount++;
		  //assert(sum <= rhstmp);
		  double testLP = 0.0;
		  for( unsigned int k = 0; k < S.size(); ++k ){
		    testLP += xlpopt.at( S.at( k ).second ) * S.at( k ).first;
		  }
#define USE_LIFTING_II
#ifdef USE_LIFTING_II
		  ///// Lifting II
		  if (decLev <= (int)/*log2*/sqrt((double)n) /*testLP <= resrhs*/) {
		    double remRhs = resrhs;
		    remS.clear();
		    for (int z = 0; z < S.size(); z++) {
		      remS.push_back(S[z]);
		    }
		    sort( S.begin(), S.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return p1.first / fabs(rowdense[p1.second]) > p2.first / fabs(rowdense[p2.second]); } ); // oder abs?
		    /*static*/ std::vector<bool> isInS(n+1);
		    for (int kk = 0; kk < K_minus_S.size();kk++) {
		      isInS[K_minus_S[kk].second] = false;
		    }
		    for (int kk = 0; kk < S.size();kk++) {
		      isInS[S[kk].second] = true;
		    }
		    for (int kk = 0; kk < K_minus_S.size();kk++) {
		      if(isInS[K_minus_S[kk].second]) {
			//K_minus_S[kk] = K_minus_S[K_minus_S.size()-1];
			//K_minus_S.pop_back();
			K_minus_S[kk].first = 0.0;
		      }
		    }
		    for (int kk = 0; kk < K_minus_S.size();kk++) {
		      if (assigns[K_minus_S[kk].second] != 2) continue;
		      if (xlpopt[K_minus_S[kk].second] > 1.0 - OneEM9) continue;
		      if (K_minus_S[kk].first == 0.0) continue;
		      //for (int o=0;o<S.size();o++)
		      //assert(K_minus_S[k].second != S[o].second);
		      //   sortiere S gem�� Quotient fk/|aik|, gro�e zuerst, fk ist der Koeffizient in S, aik der im Original
		      //for (int o=0;o<S.size()-1;o++) {
		      //cerr << S[o].first / abs(rowdense[S[o].second]) << "|" << S[o+1].first / abs(rowdense[S[o+1].second]) << " ";
		      //assert(S[o].first / abs(rowdense[S[o].second]) >= S[o+1].first / abs(rowdense[S[o+1].second]) );
		      //}
		      //cerr << endl;
		      //sort( S.begin(), S.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return p1.first / abs(rowdense[p1.second]) > p2.first / abs(rowdense[p2.second]); } );
		      /*
			for (int j = 0; j < S.size(); j++) {
			cerr << S[j].first << "x" << S[j].second << ":" << S[j].first / abs(rowdense[S[j].second])<< ":" <<  S[j].first << " + ";
			}
			cerr << " ascend order check "<< endl;
		      */

		      //   a0 = cut.rhs

		      double a0 = rhstmp;

		      //   ak = Koeffizient von koeff[K-S [k]]

		      double ak = fabs(rowdense[K_minus_S[kk].second]); // oder abs?

		      //   rhs = a0-ak

		      double lrhs = a0-ak;

		      //   volle_anteile = 0

		      double complete = 0.0;

		      //   for (int r = 0; r < S;r++)
		      //      if (koeff[S[r]] <= rhs)
		      //          rhs -= koeff[S[r]]
		      //          volle_anteile ++
		      //      else
		      //          anteil = rhs / koeff[S[r]]
		      //          fk = |S|-1-volle_anteile - floor(anteil)
		      //          if (fk>0) fuege fk*xk zu S hinzu und nimm xk aus K-S raus
		      //          break

		      int remSsize = S.size();
		      for (int r = 0; r < remSsize;r++) {

			assert(S[r].first < -OneEM12 || S[r].first > OneEM12);
			if (fabs(rowdense[S[r].second]) <= lrhs) { // oder abs?
			  lrhs = lrhs - fabs(rowdense[S[r].second]); // oder abs?
			  complete = complete + S[r].first;//1.0;
			} else {
			  double part = lrhs;// / fabs(rowdense[S[r].second]); // oder abs?
			  int fk = resrhs - complete;// - floor(part);//- floor(part*S[r].first);
			  if (fk > 0+eps) {
			    std::pair<double, unsigned int> P( fk, K_minus_S[kk].second );
			    if (0&&P.first / abs(rowdense[P.second]) <=
				S[S.size()-1].first / abs(rowdense[S[S.size()-1].second] - 0*OneEM12 - 0*(OneEM12)*abs(S[S.size()-1].first / abs(rowdense[S[S.size()-1].second])))) {
			      S.push_back(make_pair( (double)fk, K_minus_S[kk].second ));
			    } else {
			      assert(S.size() > 0);
			      int ix = S.size()-1;
			      S.push_back(make_pair(-1.0,0));
			      if (P.first / fabs(rowdense[P.second]) <= S[ix].first / fabs(rowdense[S[ix].second])) {
				S[ix+1].first = P.first;
				S[ix+1].second = P.second;
			      } else {
				while (ix > 0 && P.first / fabs(rowdense[P.second]) > S[ix].first / fabs(rowdense[S[ix].second])) { // oder abs?
				  //cerr << P.first / abs(rowdense[P.second]) << "*" << S[ix].first / abs(rowdense[S[ix].second]) << "*" << ix << " ";
				  S[ix+1].first = S[ix].first;
				  S[ix+1].second = S[ix].second;
				  ix--;
				  if (P.first / fabs(rowdense[P.second]) <= S[ix].first / fabs(rowdense[S[ix].second])) {  // oder abs?
				    ix++;
				    break;
				  }
				}
				if (ix == 0) {
				  if (S[ix].first < 0.0 || P.first / fabs(rowdense[P.second]) <= S[ix].first / fabs(rowdense[S[ix].second])) {
				    S[ix+1].first = P.first;
				    S[ix+1].second = P.second;
				  } else {
				    S[ix].first = P.first;
				    S[ix].second = P.second;
				  }
				} else {
				  S[ix].first = P.first;
				  S[ix].second = P.second;
				}
				//cerr << "A:ix=" << ix << " P=" << P.first << "x"<< P.second << " resrhs="<< resrhs<< endl;
				//cerr << "S[ix].first=" << S[ix].first << " P.first / fabs(rowdense[P.second]) <= S[ix].first / fabs(rowdense[S[ix].second])=" << (int)(P.first / fabs(rowdense[P.second]) <= S[ix].first / fabs(rowdense[S[ix].second])) << " S[ix+1].first=" << S[ix+1].first << endl;
			      }
			      assert(S[S.size()-1].first < -OneEM12 || S[S.size()-1].first > OneEM12);
			      //printcut = true;
			    }
			    /*
			      for (int j = 0; j < S.size(); j++) {
			      cerr << S[j].first << "x" << S[j].second << ":" << S[j].first / abs(rowdense[S[j].second])<< ":" <<  S[j].first << " + ";
			      }
			      cerr << " added last " << part << "," << complete << "," << rhstmp<< endl;
			    */
			    //resrhs = S.size()-1;
			  }
			  break;
			}
		      }
		    }
		    //sort( S.begin(), S.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return fabs(rowdense[p1.second])*p1.first < p2.first*fabs(rowdense[p2.second]); } );
		    double test_binwert_min=0.0;
		    for( unsigned int kk = 0; kk < resrhs+1; ++kk ){
		      test_binwert_min = test_binwert_min + S[kk].first*fabs(rowdense[S[kk].second]);
		    }
		    if (test_binwert_min <= rhstmp || S.size() < 2) {
		      S.clear();
		      res.clear();
		      found1k = false;
		      for (int o = 0; o < remS.size();o++) {
			S.push_back(remS[o]);
			res.push_back(remS[o].second);
		      }
		      resrhs = remRhs;
		      if (info_level >= 2) cerr << "lost lifting II." << endl;
#ifdef TRACK_SOL
		      cerr << "C binwert =" << test_binwert_min << " >? rhstmp =" << rhstmp << endl;
#endif
		    }
		  }
		  ///// end of lifting
		  /*static*/ std::vector<bool> isInS(n+1);

		  //cerr << "ENTER 1k LIFT II" << endl;
		  for (int kk=0;kk<S.size();kk++) {
		    isInS[S[kk].second] = false; 
		  }
		  for (int kk=0;kk<S.size();kk++) {
		    if (isInS[S[kk].second]) {
		      cerr<< "Error: 1-k cut with double element. repairing ..." << endl;
		      S[kk].first = 0.0;
		      S.clear();
		      res.clear();
		      found1k = false;
		      break;
		    } else isInS[S[kk].second] = true; 
		  }
		  sort( S.begin(), S.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return p1.first  > p2.first; } );
		  if (S.size() > 1)
		    assert(S[0] >= S[S.size()-1]);
		  while(S.size() > 0 && S[S.size()-1].first == 0.0) 
		    S.pop_back();
		  if (S.size() == 0) {
		    cerr << "Warning: 1-k cuts empty." << endl;
		    //assert(0);
		  }
#endif
		  break;
		} else {
		  check_aij = check_aij + fabs(rowdense[K_minus_S[o].second]); // oder abs?
		  //cerr << "check_aij=" << check_aij << endl;
		  check_ws = check_ws + xlpopt.at( K_minus_S.at( o ).second );
		  //if (o == K_minus_S.size()-1) cerr << "N";
		}
	      }
	    }
	  }
	  double ttestLP = 0.0;
	  for( unsigned int k = 0; k < S.size(); ++k ){
	    ttestLP += xlpopt.at( S.at( k ).second ) * S.at( k ).first;
	  }

	  if (ttestLP > resrhs) {
	    if (info_level >= 2) cerr << "/";
	    char a;
	    //cin >> a;
#ifdef TRACK_SOL
	    printcut = true;
#endif
	  }
	  if (!found1k) {
	    //cerr << "NO 1K!" << endl;
	    resrhs = 0.0;
	    S.clear();
	  } else {
	    resrhs = resrhs + OneEM6 + fabs(resrhs)*OneEM7;
	  }
	} else if (0){
	  resrhs = 1.0;
	  S.clear();
	}
      } else {
	printcut = true;//false;
	//resrhs = 1.0;
	//S.clear();
	//else cerr << "G";
      }
    } else if (0) {
      resrhs = 1.0;
      S.clear();
    }
    double testLP = 0.0;
    sumSolu=0.0;
#ifdef TRACK_SOL
    for( unsigned int k = 0; k < S.size(); ++k ){
      testLP += xlpopt.at( S.at( k ).second ) * S.at( k ).first;
      if (solu[S.at( k ).second] == 1 && rowdense.at( S.at( k ).second ) > 0)
	sumSolu += S[k].first;//fabs(rowdense.at( res[j] ));
      else if (solu[S.at( k ).second] == 0 && rowdense.at( S.at( k ).second ) < 0)
	sumSolu += S[k].first;//fabs(rowdense.at( S.at( k ).second ));
    }
    cerr << "#" << sumSolu << "<=" << resrhs <<";";
#endif
    int cntnegs=0;
    if (SHOW_COVER) {
      cerr << "THIRD rowsparse und rhs=" << resrhs << " und eps=" << eps << endl;
      for (int j = 0; j < org_rowsparse.size(); j++) {
	cerr << org_rowsparse[j].first << "x" << org_rowsparse[j].second << " + ";
      }
      cerr << " <= "<< rhs << endl;
      cerr << "FINAL-COVER:";
      for (int j = 0; j < S.size(); j++) {
	cerr << S[j].first << "x" << S[j].second << " + ";
      }
      cerr << " <= "<< resrhs << endl;
      cerr << " END "<< endl;
    }

    for( unsigned int i = 0; i < S.size(); ++i ){
      //if (S[i].first > -OneEM12 && S[i].first < OneEM12) continue;
      assert(types[S[i].second] == 0);
      if( rowdense[S[i].second] < 0 && types[S[i].second] == 0 /*BINARY*/ ){
	resrhs -= S[i].first;
	cntnegs++;
	reslhs.push_back( make_pair( S[i].second, -S[i].first));
      } else {
	reslhs.push_back( make_pair( S[i].second, S[i].first ));
      }
    }
#ifdef TRACK_SOL
    cerr << "# TO BE TURNED:" << cntnegs << endl;
#endif
  } else
    if (containsReals) {  //cerr << "resrhs=" << resrhs;
      // Cover ist in S, Koeffizientenmenge K ist im Original rowsparse; K-S ist in K_minus_S
      // sortiere ggfs. K-S, kleine Koeffizienten zuerst. Code funktioniert nur, wenn so sortiert wird. Warum?
      // nur S und K werden im Folgenden benoetigt
            
      // Step 1: extend cover
      double contOverhead = 0.0;
      for( unsigned int j = 0; j < rowsparse.size(); ++j ){
	if( types[rowsparse.at( j ).second] != 0 /*not BINARY*/){
	  unsigned int ind = rowsparse.at( j ).second;
	  double contOverhead = contOverhead + xlpopt.at( ind ) * rowsparse.at( j ).first;
	}
      }
            
      sort( K_minus_S.begin(), K_minus_S.end(), []( pair<double, unsigned int> p1, pair<double, unsigned int> p2 ){ return p1.first*xlpotmp2[p1.second] > p2.first*xlpotmp2[p2.second]; } );
                
      double curhs=rhstmp;
      if (0&&contOverhead < 0) {
	curhs -= contOverhead;
	rhstmp -= contOverhead;
	for (int i = 0; i < S.size();i++) {
	  curhs = curhs - fabs(rowdense[S[i].second]);
	  //cerr << " a_[S_" << i << "]=" << "a_[" << S[i].second <<  "]=" << rowdense[S[i].second];
	}
	if (curhs >= 0)
	  for (int u = 0;u < K_minus_S.size();u++) {
	    if (xlpotmp[K_minus_S[u].second] <= OneEM9) continue;
	    assert(xlpotmp[K_minus_S[u].second] >= 0);
	    if (fabs(rowdense[K_minus_S[u].second]) <= curhs) {
	      S.push_back(K_minus_S[u]);
	      curhs -= fabs(rowdense[K_minus_S[u].second]);
	    } else {
	      S.push_back(K_minus_S[u]);
	      curhs -= fabs(rowdense[K_minus_S[u].second]);
	      break;
	    }
	  }
      }
            
      // Step 2: build lambda, resrhs and reslhs
            
      double lambda=0.0;
      for (int i = 0; i < S.size();i++) {
	lambda = lambda + fabs(rowdense[S[i].second]);
	//cerr << " a_[S_" << i << "]=" << "a_[" << S[i].second <<  "]=" << rowdense[S[i].second];
      }
      if (info_level > 1) if (fabs(lambda - sum) > OneEM6)
			    cerr << " lambda=" << lambda << " sum=" << sum << endl;
      //assert(fabs(lambda - sum) <= 1e-6 );
      lambda = lambda - rhstmp;
            
      resrhs=0.0;
      for (int i = 0; i < S.size();i++) {
	S[i].first = (fabs(rowdense[S[i].second]) < lambda ? fabs(rowdense[S[i].second]) : lambda);
	resrhs = resrhs + S[i].first;
	if (rowdense.at( S[i].second ) < 0) {
	  resrhs -= S[i].first;
	  S[i].first *= -1.0;
	}
      }
      resrhs = resrhs - lambda + OneEM9;
            
      reslhs.clear();
      for (int i = 0; i < S.size();i++) {
	reslhs.push_back( make_pair( S[i].second, S[i].first ));
      }
            
      if (printout) cerr << " lambda=" << lambda << " neuesResrhs=" << resrhs << " und rhs =" << rhs << endl;

    }

  //if (testLP > resrhs) cerr << "G";
  //else cerr << "S";


  // Rueckkomplementieren


#ifdef OLD_REVERT
  for( unsigned int j = 0; j < rowsparse.size(); ++j ){
    if( rowsparse.at( j ).first < 0 && types[rowsparse.at( j ).second] == 0 /*BINARY*/){
      unsigned int ind = rowsparse.at( j ).second;
      xlpopt.at( ind ) = 1 - xlpopt.at( ind );
      //complemented.at( ind ) = false;
    }
  }
#endif
  for( unsigned int j = 0; j < rowsparse.size(); ++j ){
    if( org_rowsparse.at( j ).first < 0) {
      if (types[rowsparse.at( j ).second] == 0 /*BINARY*/){
	unsigned int ind = rowsparse.at( j ).second;
	xlpopt.at( ind ) = 1 - xlpopt.at( ind );
	//complemented.at( ind ) = true;
	if (rowsparse.at( j ).second > xlpopt.size() || rowsparse.at( j ).second > xlpotmp2.size()) {
	  cerr << j << " " << rowsparse.at( j ).second << " " << xlpopt.size() << " " << xlpotmp2.size() << endl;
	}
	if (fabs(xlpopt.at( rowsparse.at( j ).second ) - xlpotmp2.at( rowsparse.at( j ).second )) > OneEM12) {
	  cerr << "Fall 1: xlpopt=" << xlpopt.at( rowsparse.at( j ).second ) << " xlpotmp=" << xlpotmp2.at( rowsparse.at( j ).second ) << endl;
	}
      } else {
	unsigned int ind = rowsparse.at( j ).second;
	shift = 	-colLower[rowsparse.at( j ).second].asDouble();;
	xlpopt.at( ind ) = xlpopt.at( ind ) - shift;
	if (fabs(xlpopt.at( rowsparse.at( j ).second ) - xlpotmp2.at( rowsparse.at( j ).second )) > OneEM12) {
	  cerr << "Fall 2:xlpopt=" << xlpopt.at( rowsparse.at( j ).second ) << " xlpotmp=" << xlpotmp2.at( rowsparse.at( j ).second ) << endl;
	}
      }
    } else if (types[rowsparse.at( j ).second] != 0 /*not BINARY*/) {
      unsigned int ind = rowsparse.at( j ).second;
      shift = 	colUpper[rowsparse.at( j ).second].asDouble();
      xlpopt.at( ind ) = xlpopt.at( ind ) - shift;
      xlpopt.at( ind ) = -xlpopt.at( ind );
      // is same as: xlpopt.at( ind ) = -xlpopt.at( ind ) + shift;
      rowsparse.at( j ).first = -rowsparse.at( j ).first;
      if (fabs(xlpopt.at( rowsparse.at( j ).second ) - xlpotmp2.at( rowsparse.at( j ).second )) > OneEM5) {
	cerr << "Fall 3: xlpopt=" << xlpopt.at( rowsparse.at( j ).second ) << " xlpotmp=" << xlpotmp2.at( rowsparse.at( j ).second ) << endl;
      }
    }
    assert( fabs(xlpopt.at( rowsparse.at( j ).second ) - xlpotmp2.at( rowsparse.at( j ).second )) <= OneEM5);
  }
        
  sumSolu=0.0;
#ifdef TRACK_SOL
  for( unsigned int k = 0; k < reslhs.size(); ++k ){
    if (solu[reslhs.at( k ).first] == 1)
      sumSolu += reslhs[k].second;//fabs(rowdense.at( res[j] ));
  }
  cerr << "%" << sumSolu << "<=" << resrhs <<";";
  if (sumSolu > resrhs) printcut = true;
#endif
  if (realPart.size()>0) {
    if (printout) cerr << "jup, a mixed cut:" << realPart.size() << "," << reslhs.size() << endl;
    for( unsigned int j = 0; j < realPart.size(); j++ ){
      reslhs.push_back(std::make_pair(realPart[j].second,realPart[j].first));
      if (realPart.at( j ).first < 0) {
	double shift = -colLower[realPart.at( j ).second].asDouble();
	resrhs = resrhs - shift * realPart.at( j ).first;
	//unsigned int ind = rowsparse.at( j ).second;
	//shift = 	-colLower[rowsparse.at( j ).second].asDouble();;
	//xlpopt.at( ind ) = xlpopt.at( ind ) + shift;
	//rhstmp += shift*rowsparse.at( j ).first;
      } else {
	shift = 	colUpper[realPart.at( j ).second].asDouble();
	resrhs += shift*(realPart.at( j ).first);
	//unsigned int ind = rowsparse.at( j ).second;
	//shift = 	colUpper[rowsparse.at( j ).second].asDouble();
	//rowsparse.at( j ).first = -rowsparse.at( j ).first;
	//xlpopt.at( ind ) = -xlpopt.at( ind ) + shift;
	//rhstmp += shift*(-rowsparse.at( j ).first);
      }
    }
    if (printout) for (int j = 0; j < reslhs.size(); j++) {
	cerr << reslhs[j].second << " x" << reslhs[j].first << "(" << xlpopt[reslhs[j].first] << ") + ";
      }
    if (printout) cerr << " + 0 <= " << resrhs << endl;
    double test_binwert_min=0.0;
    for( unsigned int kk = 0; kk < reslhs.size(); ++kk ){
      test_binwert_min = test_binwert_min + reslhs[kk].second*xlpopt[reslhs[kk].first];
    }
    if (test_binwert_min <= resrhs) {
      reslhs.clear();
      resrhs = 0;
      return false;
    }
    //if (test_binwert_min > resrhs) {
    //    cerr << "WOW";
    //}
  }
  if (printcut) {
    bool isOnTrack=true;
    if(0)for (int zz = 0; zz < n; zz++) {
	if (assigns[zz] != 2/*extbool_Undef*/) {
	  if (assigns[zz] != solu[zz]) isOnTrack = false;
	} else if (fixs[zz] != 2) {
	  if (fixs[zz] != solu[zz]) isOnTrack = false;
	}
      }
    double binwert=0.0;
    double flowwert = 0.0;
    for( unsigned int k = 0; k < reslhs.size(); ++k ){
      binwert = binwert + solu[ reslhs.at( k ).first ] * reslhs.at( k ).second;
      flowwert = flowwert + xlpopt[ reslhs.at( k ).first ] * reslhs.at( k ).second;
    }
    if (isOnTrack && binwert > resrhs) {
      cerr << "show cut: ";
      for (int jj=0; jj < reslhs.size();jj++) {
	cerr << reslhs[jj].second << "x" << reslhs[jj].first << "(" << xlpopt[reslhs[jj].first] << ")" << " + ";
      }
      cerr << " 0 >= " << resrhs << endl;
      cerr << "Flow: " << flowwert << "  ? <= " << resrhs << endl;
      cerr << "ERROR: " << binwert << " not <= " << resrhs << endl;
      char a;
      cin >> a;
    }
  }
  /*
    for( unsigned int j = 0; j < rowsparse.size(); ++j ){
    if (types[rowsparse.at( j ).second] != 0) {  0=BINARY
    xlpopt[rowsparse.at( j ).second] = xlpotmp[rowsparse.at( j ).second];
    }
    }
  */
  //for (int z=0;z<xlpopt.size();z++) xlpopt[z] = xlpotmp[z];
  return true;
}

#else //FIND_BUG
#endif //FIND_BUG

#endif

#define GMI4
#ifdef GMI4
	vector< pair< vector< pair<unsigned int, double> >, double > > CutAdder::getGMICuts( extSol::QpExternSolver &extSolver, vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev, unsigned int initime, std::vector<int> &listOfCutsVars ){
		//Accuracy
		const double eps = 1e-9;
		//Starting time of GMI generation
		unsigned int gmistim = time(NULL);

		// Vector for the final cuts
		/*static*/ vector< pair< vector< pair<unsigned int, double> >, double > > cuts(candidateVariables.size(), make_pair(vector< pair<unsigned int, double> >(), 0));
		cuts.clear();
		listOfCutsVars.clear();

		// We may need to call ClpModel::status()

		// Unfortunately the status of Clp is hard coded ...
		// -1 - did not run
		//  0 - optimal
		//  1 - primal infeasible
		//  2 - dual infeasible
		//  3 - stopped on iterations or time
		//  4 - stopped due to errors
		//  5 - stopped by event handler

		if (extSolver.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL){
			const unsigned int m = extSolver.getRowCount();
			const unsigned int n = extSolver.getVariableCount();

			//Get solution values of the primal variables
			/*static*/ vector<data::QpNum> objVals(n);
			extSolver.getValues(objVals);

			//Vector for lower bounds; Default is 0, for binary Vars
			std::vector<data::QpNum> lbVec(n);
			extSolver.getLB(lbVec);

			/*static*/ vector<double> l(m + n, 0);
			if (l.capacity() < m + n) { l.reserve(m + n); l.resize(m + n); }
			for (int z = 0; z<n; z++) l[z] = lbVec[z].asDouble();

			//Vector for upper bounds; Default is 1, for binary Vars
			std::vector<data::QpNum> ubVec(n);
			extSolver.getUB(ubVec);
			/*static*/ vector<double> u(m + n, 1);
			if (u.capacity() < m + n) { u.reserve(m + n); u.resize(m + n); }
			for (int z = 0; z< n; z++) u[z] = ubVec[z].asDouble();

			//for(int i=0;i<n;i++)
				//cerr << i << " " << (int)assigns[i] << " " << lbVec[i].asDouble() << " " << ubVec[i].asDouble()<<endl;
			//cin.get();
			//Vector for original RHS of the constraint system= Right hand side vector b; Used for Substitution Slacks and Calculation of A_B^{-1}b
			/*static*/ vector<double> RHS(m, 0);
			if (RHS.capacity() < m) { RHS.reserve(m ); RHS.resize(m); }

			// Vector (Size n) for all non-basic-variables: N[i]=the i-th non-basic-variable
			/*static*/ vector<int> N(n);
			if (N.capacity() < n) { N.reserve( n); N.resize( n); }

			// NInv[i]=number of non-basic-variable; and -1 if i is basic
			/*static*/ std::vector<int> Ninv(m + n, -1);
			if (Ninv.capacity() < m + n) { Ninv.reserve(m + n); Ninv.resize(m + n); }

			// Vector (SIze m) for all Basic-Variables: B[i]=the i-th basic-variable
			/*static*/ std::vector<int> B(m);
			if (B.capacity() < m) { B.reserve(m); B.resize(m); }

			// BInv[i]=number of basic-variable; and -1 if i is non-basic
			/*static*/ std::vector<int> Binv(m + n, -1);
			if (Binv.capacity() < m + n) { Binv.reserve(m + n); Binv.resize(m + n); }
			for (int z = 0; z< m + n; z++) Binv[z] = -1;

			// Information for each variable (normal and slack)
			/*static*/ vector<int> basestat(n + m);
			if (basestat.capacity() < m + n) { basestat.reserve(m + n); basestat.resize(m + n); }

#ifdef USE_NBD_CLP
			ClpSimplex *M;
			M = (ClpSimplex*)extSolver.getSolverModel();
			M->dual(0,7);
#endif
#ifdef USE_NBD_HIGHS_X
			Highs *M;
			M = (Highs*)extSolver.getSolverEnv();
			M->setOptionValue("solver", kSimplexString);
			M->setOptionValue("simplex_strategy", kSimplexStrategyDual);
			M->run();
#endif


			/*static*/ extSol::QpExternSolver::QpExtSolBase base;
			extSolver.getBase(base);
			if (base.variables.size() < n) {
			  cerr << "Error: got no base in GMI" << endl;
			  return cuts;
			}

			/*static*/ std::vector<data::QpRhs> rhsVec;
			extSolver.getRhs(rhsVec);



			for (unsigned int i = 0; i < m; ++i){
				// For a slack variable either the upper, or lower bound is zero. Which one is used depends on ratiosign of corresponding row
				l.at(n + i) = 0;
				u.at(n + i) = 0;
				RHS.at(i) = rhsVec.at(i).getValue().asDouble();	//It is simply the b_i entry of row i
				if (rhsVec.at(i).getRatioSign() == data::QpRhs::greaterThanOrEqual){

					// Switch Lower/upper bound, because cplex only returns "AtLower" of Slacks, even if they are <=0
					// For CLP this should not do any harm, since CLP set it right in the first place
					if (base.constraints.at(i) == extSol::QpExternSolver::AtLower){
						base.constraints.at(i) = extSol::QpExternSolver::AtUpper;
					}
				}
			}

			// Make sure that the "AtLower", and "AtUpper" is set correctly for binaries
			for (unsigned int i = 0; i < n; ++i){
				if (types[i]==0 && abs(objVals[i].asDouble() - 1) < 0.5 && base.variables.at(i) == extSol::QpExternSolver::AtLower){
					base.variables.at(i) = extSol::QpExternSolver::AtUpper;
				}
				if (types[i]==0 &&abs(objVals[i].asDouble() - 1) > 0.5 && base.variables.at(i) == extSol::QpExternSolver::AtUpper){
					base.variables.at(i) = extSol::QpExternSolver::AtLower;

				}
			}

			// Save the (changed) "AtLower", "AtUpper" and "Basis" information in bastestat for both primal and slack variables
			for (unsigned int i = 0; i < n; ++i)
				basestat.at(i) = base.variables.at(i);

			for (unsigned int i = 0; i < m; ++i)
				basestat.at(n + i) = base.constraints.at(i);
			//Find out which variables are non-basic
			unsigned int ntmp = 0;
			for (unsigned int i = 0; i < n + m; ++i){
				if (basestat.at(i) != extSol::QpExternSolver::Basic){
					N.at(ntmp) = i;
					Ninv.at(i) = ntmp;
					++ntmp;
				}
				else
					Ninv.at(i) = -1;
			}
			//std::cerr << "regular n=" << n << " m=" << m << " nbasics=" << ntmp << std::endl;
			if (ntmp != n)
				throw runtime_error("Invalid basis size.");

#ifdef USE_NBD_CLP
			// GET THE ORDER OF THE BASIC VARIABLES

			//static std::vector<int> basis(m);
			//if (basis.capacity()<m) { basis.reserve(m); basis.resize(m); }
			//basis.clear();

            if (M->problemStatus() || M->algorithm()>=0) {
	      //std::cerr << "Error: GMI in trouble: Stat:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
	      //return cuts;
	      std::cerr << "Warning: GMI in trouble: Stat:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
                M->setAlgorithm(-1);
                M->dual(0,7);
                
                if (M->problemStatus() || M->algorithm()>=0) {
                    std::cerr << "Warning A: GMI in trouble: Stat:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
                    M->allSlackBasis(true);   // reset basis
                    M->dual(0,7);
                }
                if (M->problemStatus() || M->algorithm()>=0) {
                    std::cerr << "Warning B: GMI in trouble: Stat:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
                    M->allSlackBasis(true);   // reset basis
                    M->primal(0,7);
                }
                if (M->problemStatus() || M->algorithm()>=0) {
                    std::cerr << "Warning C: GMI in trouble: Stat:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
                    M->allSlackBasis(true);   // reset basis
                    M->primal(0,3);
                } else std::cerr << "info: GMI ok" << std::endl;
            


                for (int z = 0; z< m + n; z++) Binv[z] = -1;
                // Information for each variable (normal and slack)
                ///*static*/ vector<int> basestat(n + m);
                if (basestat.capacity() < m + n) { basestat.reserve(m + n); basestat.resize(m + n); }
                ///*static*/ extSol::QpExternSolver::QpExtSolBase base;
		{
		  //ClpSimplex *M;
		  //M = (ClpSimplex*)extSolver.getSolverModel();
		  //M->setAlgorithm(-1);
		  //M->dual(0,7);
		}
                extSolver.getBase(base);
		if (base.variables.size() < n) {
		  cerr << "Error: got no base in gmi ii" << endl;
		  return cuts;
		}

                ///*static*/ std::vector<data::QpRhs> rhsVec;
                extSolver.getRhs(rhsVec);
                for (unsigned int i = 0; i < m; ++i){
                    // For a slack variable either the upper, or lower bound is zero. Which one is used depends on ratiosign of corresponding row
                    l.at(n + i) = 0;
                    u.at(n + i) = 0;
                    RHS.at(i) = rhsVec.at(i).getValue().asDouble();	//It is simply the b_i entry of row i
                    if (rhsVec.at(i).getRatioSign() == data::QpRhs::greaterThanOrEqual){
                        
                        // Switch Lower/upper bound, because cplex only returns "AtLower" of Slacks, even if they are <=0
                        // For CLP this should not do any harm, since CLP set it right in the first place
                        if (base.constraints.at(i) == extSol::QpExternSolver::AtLower){
                            base.constraints.at(i) = extSol::QpExternSolver::AtUpper;
                        }
                    }
                }
                // Make sure that the "AtLower", and "AtUpper" is set correctly for binaries
                for (unsigned int i = 0; i < n; ++i){
                    if (types[i]==0 && abs(objVals[i].asDouble() - 1) < 0.5 && base.variables.at(i) == extSol::QpExternSolver::AtLower){
                        base.variables.at(i) = extSol::QpExternSolver::AtUpper;
                    }
                    if (types[i]==0 &&abs(objVals[i].asDouble() - 1) > 0.5 && base.variables.at(i) == extSol::QpExternSolver::AtUpper){
                        base.variables.at(i) = extSol::QpExternSolver::AtLower;
                        
                    }
                }
                // Save the (changed) "AtLower", "AtUpper" and "Basis" information in bastestat for both primal and slack variables
                for (unsigned int i = 0; i < n; ++i)
                    basestat.at(i) = base.variables.at(i);
                for (unsigned int i = 0; i < m; ++i)
                    basestat.at(n + i) = base.constraints.at(i);
                //Find out which variables are non-basic
                ntmp = 0;
                for (unsigned int i = 0; i < n + m; ++i){
                    if (basestat.at(i) != extSol::QpExternSolver::Basic){
                        N.at(ntmp) = i;
                        Ninv.at(i) = ntmp;
                        ++ntmp;
                    }
                    else
                        Ninv.at(i) = -1;
                }

		//std::cerr << "n=" << n << " m=" << m << " nbasics=" << ntmp << std::endl;
                
                if (M->problemStatus() || M->algorithm()>=0) {
                    std::cerr << "Warning II: GMI in trouble: Stat:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
		    return cuts;
                    M->allSlackBasis(true);   // reset basis
                    //extSolver.prepareMatrixRowForm();
                    M->primal(0,7);
                    //M->initialSolve();
                    //ClpDualRowSteepest dualSteepest;
                    //M->setDualRowPivotAlgorithm(dualSteepest);
                    //M->allSlackBasis(true);   // reset basis
                    //M->dual(0,1);
                   if (M->problemStatus() || M->algorithm()>=0) {
                      std::cerr << "Error: GMI not allowed." << std::endl;
                      return cuts;
                   }
                }
            }

            //CoinWarmStartBasis *wsb = M->getBasis();
			std::vector<int> head(m);
			if (head.capacity()<m) { head.reserve(m); head.resize(m); }
			//std::cerr << "G";

			M->getBasics(&(head.data()[0]));
			//std::cerr << "g";
			/*for (int i = 0; i < n; i++) {
				if (wsb->getStructStatus(i) == CoinWarmStartBasis::basic) {
					basis.push_back(i);
				}
			}
			for (int i = 0; i < m; i++) {
				if (wsb->getArtifStatus(i) == CoinWarmStartBasis::basic) {
					basis.push_back(n + i);
				}
			}

			//wsb->print();
			delete wsb;
			*/
			for (unsigned int i = 0; i < m; ++i){
							B.at(i) = head.at(i);
							if (Ninv[B.at(i)]!=-1) {
							  std::cerr << "Error in get GMI cuts. Ninv=" << Ninv[B.at(i)] << ", B[i]=" << B[i] << std::endl;
								return cuts;
							}
							assert(Ninv[B.at(i)]==-1);
							Binv.at(B.at(i)) = i;
			}
#endif
#ifdef USE_NBD_HIGHS
	Highs *highsPt = (Highs*)extSolver.getSolverEnv();
	Highs &highs = *highsPt;
	const HighsInfo& info = highs.getInfo();
	//const unsigned int m = getRowCount();
	//const unsigned int n = getVariableCount();

	//std::vector<int> B(m);
	//std::vector<int> Binv(m + n, -1);
	//std::vector<int> N(n);

	//std::vector<int> Ninv(m + n, -1);
	std::vector<int> basis(m);
	basis.clear();
	//lhs.clear();
	const HighsBasis& highs_basis = highs.getBasis();
	const bool has_basis = info.basis_validity;

	if (!has_basis) {
	  basis.clear();
	  std::cerr << "Exception caught getting base. No basis in GMI. " << std::endl;
	  return cuts;
	}

	std::vector<HighsInt> head(m);
	if (head.capacity()<m) { head.reserve(m); head.resize(m); }
	/* getBasicVariablesArray:
	 * Entries are indices of columns if in [0,num_col), and entries in [num_col,
	 * num_col+num_row) are (num_col+row_index).
	 */
	/* getBasicVariables:
	 * Non-negative entries are indices of columns,
	 * and negative entries are -(row_index+1).
	 */
	//const HighsInt* head2 = highs.getBasicVariablesArray();
	highs.getBasicVariables(&(head.data()[0]));

	/*for (unsigned int i = 0; i < m; ++i){
	  B.at(i) = head[i];
	  if (Ninv[B.at(i)]!=-1) {
	    std::cerr << "Error in get GMI cuts. Ninv=" << Ninv[B.at(i)] << ", B[i]=" << B[i] << std::endl;
	    return cuts;
	  }
	  assert(Ninv[B.at(i)]==-1);
	  Binv.at(B.at(i)) = i;
	}*/
	for (unsigned int i = 0; i < m; ++i){
	  if (head[i] >= 0){
	    //assert(head[i] == head2[i]);
	    B.at(i) = head[i];
	  }
	  else {
	    //assert(n-(head[i]+1) == head2[i]);
	    B.at(i) = n - (head[i] + 1);    
	  }
	  Binv.at(B.at(i)) = i;
	}
	
#endif
#ifdef USE_NBD_CPLEX_C
			// GET THE ORDER OF THE BASIC VARIABLES
			std::vector<CPXDIM> head(m);
			if (head.capacity()<m) { head.reserve(m); head.resize(m); }
			/*
			* An array. The array contains the indices of the variables in the resident basis,
			* where basic slacks are specified by the negative of the corresponding row index minus 1 (one);
			* that is, -rowindex - 1.
			* The array must be of length at least equal to the number of rows in the LP problem object.
			*/
			if (CPXXgetbhead(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), head.data(), NULL)){
				throw utils::ExternSolverException("Unable to obtain basis.");
			}

			for (unsigned int i = 0; i < m; ++i){
				if (head.at(i) >= 0){
					B.at(i) = head.at(i);
				}
				else {
				        B.at(i) = n - (head.at(i) + 1);    
				}
				Binv.at(B.at(i)) = i;
			}
#endif
			std::vector< std::vector< data::IndexedElement >  > allrows2;
#ifdef USE_NBD_CLP
					//Needed later...
			extSolver.prepareMatrixRowForm();
			for (int i = 0; i < m;i++) {
			  std::vector< data::IndexedElement > T;
			  allrows2.push_back(T);
			  extSolver.getRowLhs/*_prep*//*_rows*/(i, allrows2[i]);
			}
#endif
#ifdef USE_NBD_HIGHS
			//Needed later...
			extSolver.prepareMatrixRowForm();
			for (int i = 0; i < m;i++) {
			  std::vector< data::IndexedElement > T;
			  allrows2.push_back(T);
			  extSolver.getRowLhs(i, allrows2[i]);
			  /*for (int j = 0 ; j < allrows2[i].size();j++)
			    cerr << allrows2[i][j].value.asDouble() << "y" << allrows2[i][j].index << " ";
			  cerr << endl;
			  */
			}
#endif
#ifdef USE_NBD_CPLEX_C
	      {
		  //void QpExtSolCplexC::getRowLhs(unsigned int ri, std::vector<data::IndexedElement>& lhs) {
		  int vars = n;

		  static std::vector<double> rowtmp(vars);
		  static std::vector<CPXDIM> rowindtmp(vars);
		  /*static*/ std::vector<CPXNNZ> rowstarts(m);
		  CPXNNZ size = 0, foo = 0, surplus = 0;
		  int missing = CPXXgetrows(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), &size, rowstarts.data(), rowindtmp.data(),rowtmp.data(), 0, &surplus, 0, m-1);
		  assert(surplus < 0);
		  missing = -surplus;
		  rowtmp.resize(missing);
		  rowindtmp.resize(missing);

		  if ((missing=CPXXgetrows(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), &size, rowstarts.data(), rowindtmp.data(), rowtmp.data(), rowtmp.size(), &surplus, 0, m-1))) {
		  std::cerr << "Error: Extern Solver denies row."<< "all " << surplus << " " << missing << std::endl;
		    return cuts;
		  }
		  int crow=0;
		  int i =0;
		  while (i < size) {
		      std::vector< data::IndexedElement > T;
		      allrows2.push_back(T);
		      do {
		         data::IndexedElement e;
                         e.index = rowindtmp[i];
                         e.value = rowtmp[i];
		         allrows2[crow].push_back(e);
			 i++;
		      } while ((crow < rowstarts.size()-1) ? i < rowstarts[crow+1] : i < size);
		      crow++;
		  }

	   }
#endif


			// Try creating GMI-cut for each candidate
			for (unsigned int iind = 0; iind < candidateVariables.size(); ++iind){
			  //if ((double)(time(NULL) - gmistim) > 0.02*(double)(time(NULL) - initime)) break;
			  //if(cuts.size()>sqrt(n)) break;
				unsigned int ind = candidateVariables[iind];

				//If variable is too close to being integer: continue
				if (Binv[ind] == -1 || abs(objVals.at(ind).asDouble() - floor(objVals.at(ind).asDouble() + 0.5)) <= eps){
                    if (Binv[ind] == -1 && abs(objVals.at(ind).asDouble() - floor(objVals.at(ind).asDouble() + 0.5)) > eps)
                        cerr << "Warning: GMI-cut failed." << endl;
				}
				else {
					assert(Binv[ind] != -1);
					//Vector for the corresponding row of the (solved) simplex tableau. Equal to A_B^{-1}A
#ifdef USE_NBD_HIGHS
					std::vector<data::QpNum> SimplexRow(n + m);
#endif
#ifdef USE_NBD_CPLEX_C
					std::vector<data::QpNum> SimplexRow(n + m);
#endif
#ifdef USE_NBD_CLP
					std::vector<double> SimplexRow(n + m, 0);
#endif
					for (int i=0; i < n+m;i++)
					  SimplexRow[i] = 0.0;

					//Vector of the corresponding inverse-basic-part-matrix: A_B^{-1}
					std::vector<double> BInvRow(m);

					// Not sure why two different vectors are needed...
					/*static*/ vector<double> gmicoeff(n);
					/*static*/ vector<data::QpNum> gmicoeff2;
					gmicoeff2.resize(n);
					//extSolver.getRowLhsOfTableauByColumn( ind, gmicoeff2 );
#ifdef USE_NBD_HIGHS
					//std::vector<double> BInvRow(m);
					HighsInt row_num_nz=0;
					/*HVector row_ep;
					row_ep.setup(m);

					highs.getBasisInverseRowSparse(Binv[ind], row_ep);
					for (int i=0;i<m;i++) {
					  BInvRow[i] = 0.0;
					}
					for (int i=0; i < row_ep.count;i++) {
					  BInvRow[row_ep.index[i]] = row_ep.array[row_ep.index[i]];
					}*/
					highs.getBasisInverseRow(Binv[ind], BInvRow.data(), &row_num_nz,NULL);
					for (int i=0;i<m;i++)
                                          BInvRow[i] = -BInvRow[i];
					extSolver.getBinvArow(Binv[ind], SimplexRow);
					
					for (unsigned int i = 0; i < m; ++i){
						if (Ninv[n + i] != -1)
							SimplexRow[n + i] = BInvRow[i];
					}

					// Get the coefficients of the non-basic-variables of the row of interest
					for (unsigned int i = 0; i < n; ++i){
					  gmicoeff2[i] = SimplexRow[N[i]].asDouble();
					}
#endif
#ifdef USE_NBD_CPLEX_C

					if (CPXXbinvrow(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), Binv[ind], BInvRow.data())){
						throw utils::ExternSolverException("CPXXbinvrow error");
					}


					extSolver.getBinvArow(Binv[ind], SimplexRow);
					//This only returns A_B^{-1}A for primal variables; Slack part is still missing
					//if (CPXXbinvarow(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), Binv[ind], SimplexRow.data())){
					//	throw utils::ExternSolverException("CPXXbinvarow error");
					//}


					//complete simplex row to (A_B^{-1}(A   I))= (A_B^{-1}A    A_B^{-1}), for slack columns
					for (unsigned int i = 0; i < m; ++i){
						if (Ninv[n + i] != -1)
							SimplexRow[n + i] = BInvRow[i];
					}

					//ONCE USED FOR CHECKING CPLEX_SLACK PROPERTIES
					/*static vector<double> slacks(m);
					if (slacks.size() < m) slacks.resize(m);
					int status = CPXXgetslack(
						*(CPXENVptr*)extSolver.getSolverEnv(),
						*(CPXLPptr*)extSolver.getSolverModel(),
						slacks.data(), 0, m - 1);
					for (int i = 0; i < m; i++)if (BInvRow[i] == 0){
						cout << "SLACK-VALS " << BInvRow[i] << " " << rhsVec.at(i).getRatioSign() << " " << slacks.at(i) << endl;
					}
					cin.get();
					*/


					// Get the coefficients of the non-basic-variables of the row of interest
					for (unsigned int i = 0; i < n; ++i){
					  gmicoeff2[i] = SimplexRow[N[i]].asDouble();
					}

					//Each entry of basics should be zero, unless it is the candidate variable Just for checking
					/*for (unsigned int i = 0; i < m; ++i){
						if (B[i] != ind)
							assert(SimplexRow[B[i]] == 0);

						if (B[i] ==ind)
							assert(SimplexRow[B[i]]==1);
					}*/

#endif
#ifdef USE_NBD_CLP
					//get complete simplex row (A_B^{-1}(A   I))= (A_B^{-1}A    A_B^{-1}) including slack columns

					M->getBInvARow(Binv[ind], &(SimplexRow.data()[0]), &(SimplexRow.data()[n]));
					M->getBInvRow(Binv[ind], &(BInvRow.data()[0]));

					//Needed later...
					//extSolver.prepareMatrixRowForm();

					gmicoeff2.resize(n);

					for (unsigned int i = 0; i < m; ++i){
						if (Ninv[n + i] != -1)
							SimplexRow[n + i] = BInvRow[i];
					}
					/*for (unsigned int i = 0; i < m; ++i){
						assert(SimplexRow[n+i] == BInvRow[i]);
						if (B[i] != ind){
							assert(SimplexRow[B[i]] == 0);
						}
					}*/
					for (unsigned int i = 0; i < n; ++i)
						gmicoeff2[i] = SimplexRow[N[i]];
#endif

					if (gmicoeff2.size() != n){
						continue;
						//throw runtime_error( "unexpected" );
					}

					//gmicoeff stores the simplex row entries of the non-basic variables
					for (unsigned int i = 0; i < n; ++i)
						gmicoeff[i] = gmicoeff2[i].asDouble();

					// Generate GMI Cut

					double ai0 = 0;
                    double maxRhs=0.0;
					//ai0=A_B^-1b
                    for (int i = 0; i < m; i++) {
                        if (fabs(RHS.at(i))>maxRhs) maxRhs = fabs(RHS.at(i));
#ifdef USE_NBD_HIGHS
			ai0 += SimplexRow.at(n+i).asDouble()*RHS[i];
#endif
#ifdef USE_NBD_CLP
			ai0 += SimplexRow.at(n+i)*RHS[i];
#endif
#ifdef USE_NBD_CPLEX
			ai0 += SimplexRow.at(n+i).asDouble()*RHS[i];
#endif
                    }

					//cerr<< "Bef "<<ai0<< " "<<objVals.at(ind).asDouble()<<endl;

					//THE FOLLOWING STEPS ARE NEEDED, SINCE GMI CUTS ARE CREATED FOR Non-BASICS IN R+ AT LOWER BOUND 0
					// Transform each non-basic-variable into a non-negative variable at it's lower bound with value 0
					for (unsigned int k = 0; k < n; ++k){
						unsigned int nind = N[k];
						if (basestat.at(nind) == extSol::QpExternSolver::AtLower){
							// Lower Bound
							// If X is at lower bound: replace X=LB+Xnew
							ai0 -= gmicoeff[k] * l[nind];
						}
						else if (basestat.at(nind) == extSol::QpExternSolver::AtUpper) {
							// Upper Bound
							// If X is at upper bound: replace X=UB-Xnew
							ai0 -= gmicoeff[k] * u[nind];
							gmicoeff[k] *= -1;

						}

					}
					//cerr<< "A "<<ai0<< " "<<objVals.at(ind).asDouble()<<endl;
					// BETTER EPSILON???
					//ai0 = objVals.at(ind).asDouble();
					if (0&&abs(ai0- objVals.at(ind).asDouble())>/*1e-12 * maxRhs +*/ eps) {
						cerr << "Error: in GMI: " << abs(ai0- objVals.at(ind).asDouble()) << endl;
                        continue;
						break;
					}
					//assert(abs(ai0- objVals.at(ind).asDouble())<=eps);

					// Maybe objVals is more accurate!? Calculating ai0 in the lines before is just used for checking; Could be omitted to safe time
                    if (ai0<1e-9 || objVals.at(ind).asDouble() < 1e-9) continue;
		    //ai0 = objVals.at(ind).asDouble();
                    if (ai0 < 0) continue;
					assert(ai0 >= 0);

					// Cut from row
					vector<double> cutcoeff(n);
					double fi0 = ai0 - floor(ai0);
					if (abs(fi0) <= eps){
						// k-cut separation impossible
						continue;
					}

					// Here the well-known GMI formula is finally used
					for (unsigned int k = 0; k < n; ++k){
						unsigned int nind = N[k];
						double & aij = gmicoeff[k];
						if (nind < n &&  types[nind]==0){
							// Integer nonbasic Variables
							double fij = aij - floor(aij);
							if (fij <= fi0){
								cutcoeff[k] = fij;
							}
							else {
								cutcoeff[k] = fi0*(1 - fij) / (double)(1 - fi0);
							}
						}
						else {
							// Real nonbasic variables
							if (aij >= 0){
								cutcoeff[k] = aij;
							}
							else {
								cutcoeff[k] = fi0*(-aij) / (double)(1 - fi0);
							}
						}
					}

					// undo variable substititions
					double beta = fi0;// -OneEM12 - abs(fi0)*OneEM12;
					for (unsigned int k = 0; k < n; ++k){
						//get the k-th non-basic variable, which is then nind
						unsigned int nind = N[k];
						if (basestat.at(nind) == extSol::QpExternSolver::AtLower){
							// Lower Bound
							// Re-Substitute Xnew=X-LB;
							beta += cutcoeff[k] * l[nind];
						}
						else if (basestat.at(nind) == extSol::QpExternSolver::AtUpper) {
							// Upper Bound
							// Re-Substitute Xnew=UB-X;
							beta -= cutcoeff[k] * u[nind];
							cutcoeff[k] *= -1;
						}
					}

					vector<double> cut(n, 0);
					for (int i = 0; i < n; i++) cut[i] = 0; //Just to be sure
					{
					  std::vector<data::IndexedElement> *rowtmp = NULL;
					  for (unsigned int k = 0; k < n; ++k){
					    unsigned int nind = N[k];
					    if (nind < n){
					      //If variable is primal
					      cut[nind] += cutcoeff[k];
					    }
					    else {
					      //If slack: substitute slack in terms of the corresponding row ax+s=b
					      if (cutcoeff[k] != 0){ //Or too close to zero....
						const unsigned int row = nind - n;
						rowtmp = &allrows2[row];
						//extSolver.getRowLhs/*_rows*/(row, rowtmp);
						
						for (unsigned int l = 0; l < rowtmp->size(); ++l){
						  cut[rowtmp->at(l).index] -= cutcoeff[k] * rowtmp->at(l).value.asDouble();
						}
						beta -= cutcoeff[k] * RHS.at(row);
					      }
					    }
					  }
					  if (rowtmp ==NULL || rowtmp->size() == 0){
					    break;
					  }
					}
					
					vector<pair<unsigned int, double> > cutvec;
					double max_c = 0.0;
					for( unsigned int i = 0; i < n; ++i ){
					  if( abs( cut.at( i ) ) > max_c ){
					    max_c = abs( cut.at( i ));
					  }
					}
					for( unsigned int i = 0; i < n; ++i ){
					  if (abs( cut.at( i )) > 1e-20) {
					    if( abs( cut.at( i ) ) > max_c * 1e-15 || types[i] != 0 ){
					      cutvec.push_back( make_pair( i, cut.at( i ) ) );
					    } else if (cut[i]<0) beta = beta - cut.at( i );
					  }
					}

					/*
					//if (lhs.size() > 0 && computeCutRatio(lhs) > 10000000) {
					std::sort(cutvec.begin(),cutvec.end(),[](std::pair<int, double> p1, std::pair<int,double> p2){ return fabs(p1.second) > fabs(p2.second); });
					double minDelta = fabs(beta)*0.001 + 0.01;
					//cerr << "1: first elem:" << lhs[0].value.asDouble() << " last:" << lhs[lhs.size()-1].value.asDouble() << endl;
					//if(0)cerr << "-d-";lhs.clear(); rhs = 0.0; }
					
					double llhs=0.0;
					double lrhs=beta;
					std::vector<data::QpNum> ubs;
					double mini=1e300;
					double maxi=-1e300;
					for (int h=0;h<cutvec.size();h++) {
					  if (fabs(cutvec[h].second) < mini) mini=fabs(cutvec[h].second);
					  if (fabs(cutvec[h].second) > maxi) maxi=fabs(cutvec[h].second);
					  llhs = llhs + cutvec[h].second*(objVals[cutvec[h].first].asDouble());
					}
					while (cutvec.size() > 0 && fabs(cutvec[0].second)/fabs(cutvec[cutvec.size()-1].second) > 100000) {
					  if (llhs+minDelta+fabs(cutvec[cutvec.size()-1].second) < beta) {
					    beta = beta - fabs(cutvec[cutvec.size()-1].second);
					    cutvec.pop_back();
					  } else break;
					}
					//cerr << "2: first elem:" << lhs[0].value.asDouble() << " last:" << lhs[lhs.size()-1].value.asDouble() << endl;
					//cerr << "stage " << stage << " lhs=" << llhs << " < " << lrhs << "=rhs --- min=" << mini << " max=" << maxi << endl; 

					*/
					
					if (cutvec.size() > 0 /*&&
										  cleanCut(
										  cutvec, cutvec.size(), beta, objVals,
										  colLower, colUpper, n, types)*/) {
						//cuts.at( iind ) = make_pair( cutvec, beta );
						cuts.push_back(make_pair(cutvec, beta));
						listOfCutsVars.push_back(candidateVariables[iind]);
					}
					/*if (cutvec.size() > 0 && cleanCut(cutvec, cutvec.size(), beta, objValues, n, rhsVec) ) {
					OsiRowCut rc;
					rc.setRow(cutNz, cutIndex, cutElem);
					rc.setLb(-param.getINFINIT());
					rc.setUb(cutRhs);
					cs.insert(rc);
					}*/
				}
			}
		}
		else cerr << "LP status inappropriate for GMI cut generation. " << extSolver.getSolutionStatus() << endl;

		return cuts;
	}

#endif

#define LPCuts
#ifdef LPCuts


vector< pair< vector< pair<unsigned int, double> >, double > > CutAdder::getLPCuts(extSol::QpExternSolver &extSolver, vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev, unsigned int initime, std::vector<int> &listOfCutsVars){

  //Output vector "cuts"
  /*static*/ vector< pair< vector< pair<unsigned int, double> >, double > > cuts(candidateVariables.size(), make_pair(vector< pair<unsigned int, double> >(), 0));
  cuts.clear();
  listOfCutsVars.clear();
  if (info_level >= 2) 
    cerr <<"Get L&P Cuts"<<endl;
  double epsNull = 1e-7;//1e-9;

  if (extSolver.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL){

    unsigned int m = extSolver.getRowCount();
    const unsigned int n = extSolver.getVariableCount();

    vector<data::QpNum> x_star(n);
    extSolver.getValues(x_star);
    vector< pair<int, double> > sortedCand;
    data::QpNum dummyZero = 0.0;
    data::QpNum dummyOne = 1.0;

    //Dummy vector with UNsigned entries (1,...,m)
    vector<unsigned int> uRowCount(m);
    //Dummy vector with entries (1,...,m)
    vector<int> RowCount(m);


    // Sort candidate variables ascending regarding their solution value
    for (int i = 0; i < candidateVariables.size(); i++){
      if (1 - x_star[candidateVariables[i]].asDouble()>epsNull && x_star[candidateVariables[i]].asDouble() > epsNull)
	sortedCand.push_back(make_pair(candidateVariables[i], x_star[candidateVariables[i]].asDouble()));
    }
    sort(sortedCand.begin(), sortedCand.end(), [](pair<int, double> p1, pair<int, double> p2){ return p1.second < p2.second; });

    // Store LP entries: rhs, obj, bounds, etc.

    //RightHandSide vector, with value and sense
    std::vector<data::QpRhs> SaveRhsVec(m);

    //RightHandSide vector, only value ATTENTION: Initialized LATE
    std::vector<data::QpNum> SaveRhsVals(m);
    extSolver.getRhs(SaveRhsVec);

    std::vector<data::QpNum> SaveLbVec(n);
    extSolver.getLB(SaveLbVec);

    std::vector<data::QpNum> SaveUbVec(n);
    extSolver.getUB(SaveUbVec);

    int SaveObjSense = 0;

    //righthand side minus row activity (for not(!) ranged rows)
    vector<double> SaveSlacks(m);
    std::vector<double> SaveObjCoeffs(n);

    //Vector for changing sense of row to "Ranged"
    vector<char> SenseRow(m);
#ifdef USE_NBD_CLP
    // Get the Clp-Model
    ClpSimplex *M;
    M = (ClpSimplex*)extSolver.getSolverModel();
    
    // Get the LHS-activity of all rows
    double * ActivityCLP=NULL;
    ActivityCLP=(double*)M->getRowActivity();
    assert(ActivityCLP!=NULL);

    // Save the objective coefficients
    double * SaveObjCoeffsCLP =NULL;
    SaveObjCoeffsCLP=(double *) M->getObjCoefficients();
    assert(SaveObjCoeffsCLP!=NULL);
    for (int i=0;i<n;i++)
      SaveObjCoeffs[i]=SaveObjCoeffsCLP[i];

    //Save sense of optimization
    SaveObjSense=M->getObjSense();

    //Get the Slack Values of each row by calculating RHS-LHS

    for (int row=0;row<m;row++){
      SaveSlacks[row]=SaveRhsVec[row].getValue().asDouble()- ActivityCLP[row];
      //cerr << SaveSlacks[row] << endl;
      /*if(SaveRhsVec.at(row).getRatioSign()== data::QpRhs::smallerThanOrEqual&&SaveSlacks[row]<=-epsNull)
	cerr << "Slack negative "<< SaveSlacks[row]<<endl;
	//assert(SaveSlacks[row]>=-epsNull);
	else if (SaveRhsVec.at(row).getRatioSign()== data::QpRhs::greaterThanOrEqual&& SaveSlacks[row]>=epsNull)
	cerr << "Slack positive "<< SaveSlacks[row]<<endl;
	///assert(SaveSlacks[row]<=epsNull);
	else if (SaveRhsVec.at(row).getRatioSign()== data::QpRhs::equal && abs(SaveSlacks[row])>=epsNull)
	cerr << "Slack not 0 "<< SaveSlacks[row]<<endl;

	//assert(abs(SaveSlacks[row])<=epsNull);
	*/
    }
#endif
#ifdef USE_NBD_CPLEX_C
    int statusSlack = CPXXgetslack(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), SaveSlacks.data(), 0, m - 1);
    int statusObj = CPXXgetobj(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), SaveObjCoeffs.data(), 0, n - 1);
    SaveObjSense = CPXXgetobjsen(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel());
    if (statusSlack+statusObj!=0)
      throw utils::ExternSolverException("CPXXgetslack or CPXXgetobj error");

#endif

    // Create MLP instance for variable k:
    //max y_k(-x*_k)
    //s.t.: Ay = b*x_star[k] (f�r bin�re Programme)
    //0<=y<=x*

    // Clear objective function, set LB to Zero, set UB to x_star
    for (int i = 0; i < n; i++){
      extSolver.changeObjFuncCoeff(i, dummyZero);
      extSolver.setVarLB(i, dummyZero);
      data::QpNum dummy = x_star[i];
      extSolver.setVarUB(i, dummy);
    }

    for (int j = 0; j < m; j++){
      SenseRow[j] = 'R';
      RowCount[j] = j;
      uRowCount[j] = j;
    }


#ifdef USE_NBD_CPLEX_C
    CPXXchgobjsen(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), CPX_MAX);
    CPXXchgsense(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), m, RowCount.data(), SenseRow.data());
#endif
#ifdef USE_NBD_CLP
    //Set sense to maximization
    M->setOptimizationDirection(-1);
    // Setting the sense of the rows to ranged in upcomming USE_NBD_CLP block
#endif
    //---------------------------------------------------------------------------------------------------
    // Start generating cuts for each candidate variable

    for (int k_ind = 0; k_ind < sortedCand.size() /*&& k_ind < 10*/; k_ind++){
      const int k = sortedCand[k_ind].first;
      assert(k < n);
      assert(sortedCand[k_ind].second == x_star[k].asDouble());

      //Set new objective function
      extSolver.changeObjFuncCoeff(k, dummyOne);


      //Set Lower and Upper Bound of constraints
      //Ranged constraints have the form:-   RHS+Range <= aTx <= RHS, if Range is negative
      //-RHS <= aTx <= RHS+Range, if Range is positive
      bool DangerousSlack = false;
      vector<double> NewRange(m);
      std::vector<data::QpNum> NewRhs(m);
      if (k_ind != 0 && sortedCand[k_ind - 1].second == sortedCand[k_ind].second) goto SameRHS;

      for (int j = 0; j < m; j++){

	NewRhs[j] = x_star[k] * SaveRhsVec[j].getValue();
	//double RHS_VAL = x_star[k].asDouble() * SaveRhsVec[j].getValue().asDouble();
	//CPXXchgcoef(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), j, -1, RHS_VAL);
	switch (SaveRhsVec.at(j).getRatioSign()){
	case data::QpRhs::smallerThanOrEqual:
	  //cerr << "S";
	  // Constraint j will have the form: RHS_VAL-SaveSlack[j] <= A_j*x <= RHS_VAL
	  if (SaveSlacks[j] >= 0){
	    NewRange[j] = -SaveSlacks[j];
	    //CPXXchgcoef(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), j, -2, -SaveSlacks[j]);
	    break;
	  }
	  else if (SaveSlacks[j] >= -epsNull){
	    NewRange[j] = 0;
	    //CPXXchgcoef(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), j, -2, 0);
	    break;
	  }
	  else{
	    DangerousSlack = true;
	    if (info_level >= 2) cerr << "Warning: Negative Slack: " << SaveSlacks[j] << endl;
	    break;
	  }

	case data::QpRhs::equal:
	  //cerr << "EQ";
	  // Constraint j will have the form: RHS_VAL <= A_j*x <= RHS_VAL
	  if (abs(SaveSlacks[j]) <= epsNull){
	    NewRange[j] =0;

	    //CPXXchgcoef(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), j, -2, 0);
	    break;
	  }
	  else{
	    DangerousSlack = true;
	    if (info_level >= 2) cerr << "Warning: Non-Zero Slack: " << SaveSlacks[j] << endl;
	    break;
	  }

	case data::QpRhs::greaterThanOrEqual:
	  //cerr << "G";
	  // Constraint j will have the form: RHS_VAL <= A_j*x <= RHS_VAL-SaveSlack[j]
	  if (SaveSlacks[j] <= 0){
	    NewRange[j] = -SaveSlacks[j];
	    //CPXXchgcoef(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), j, -2, -SaveSlacks[j]);
	    break;
	  }
	  else if (SaveSlacks[j] <= epsNull){
	    NewRange[j] = 0;
	    //CPXXchgcoef(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), j, -2, 0);
	    break;
	  }
	  else{
	    DangerousSlack = true;
	    if (info_level >= 2) cerr << "Warning: Positive Slack: " << SaveSlacks[j] << endl;
	    break;
	  }
	default:
	  throw runtime_error("Unexpected sense.");
	  break;
	}
	if (DangerousSlack) goto NoCuts;

	//Possible: extSolver.changeRhsElements might be faster!?
      }
#ifdef USE_NBD_CPLEX_C
      if(CPXXchgrngval(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), m, RowCount.data(), NewRange.data()))
	throw utils::ExternSolverException("CPXXchgrngval error");
      extSolver.changeRhsElements(uRowCount, NewRhs);

#endif
#ifdef USE_NBD_CLP
      // for each row: set new bounds
      for (int j = 0; j < m; j++){
	//cerr << NewRhs[j].asDouble()<< " " <<NewRange[j]<<endl;
	switch (SaveRhsVec.at(j).getRatioSign()){
	case data::QpRhs::smallerThanOrEqual:

	  assert(NewRhs[j].asDouble()+NewRange[j]<=NewRhs[j].asDouble());
	  M->setRowBounds(j,NewRhs[j].asDouble()+NewRange[j],NewRhs[j].asDouble());
	  break;

	case data::QpRhs::equal:
	  M->setRowBounds(j,NewRhs[j].asDouble(),NewRhs[j].asDouble());
	  break;

	case data::QpRhs::greaterThanOrEqual:
	  assert(NewRhs[j].asDouble()+NewRange[j]>=NewRhs[j].asDouble());
	  M->setRowBounds(j,NewRhs[j].asDouble(),NewRhs[j].asDouble()+NewRange[j]);
	  break;
	default:
	  throw runtime_error("Unexpected sense.");
	  break;
	}
      }
#endif

      //RUN

      // Solve the MLP
    SameRHS:
      extSolver.solve(100000, 100000);
      if (info_level >= 2) cerr << "$";// << (extSolver.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL) << " " << extSolver.getObjValue() - x_star[k] << endl;


      // If optimal solution found and objval (minus constant) is negative
      if (extSolver.getSolutionStatus() == extSol::QpExternSolver::OPTIMAL&& extSolver.getObjValue() - x_star[k] < 0){ //or <-epsNull
	vector<double> Slacks(m);
#ifdef USE_NBD_CPLEX_C
	int status = CPXXgetslack(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), Slacks.data(), 0, m - 1);
	assert(status == 0);
#endif
#ifdef USE_NBD_CLP
	//Get the new Slacks
	double * NewActivityCLP=NULL;
	NewActivityCLP=(double*)M->getRowActivity();
	assert(NewActivityCLP!=NULL);

	for (int row=0;row<m;row++){
	  Slacks[row]=M->getRowUpper()[row]- NewActivityCLP[row];
        if (Slacks[row]<-1e-3/*-epsNull*/) {
            cerr << "Error: Numerical instability at cut generation: " << Slacks[row] << endl;
            goto NoCuts;
        }
	  //assert(Slacks[row]>=-epsNull);
	}
#endif
	// basic/lower-bound/upper-bound information for variables and constraints/slacks
	/*static*/ vector<int> basestat(n + m);
	if (basestat.capacity() < n + m) { basestat.reserve(n + m); basestat.resize(n + m); }
	// Basis und Schranken der Schlupfvariablen auslesen
	/*static*/ extSol::QpExternSolver::QpExtSolBase base;
	extSolver.getBase(base);
	if (base.variables.size() < n) {
	  cerr << "Error: got no base" << endl;
	  goto NoCuts;
	}
	for (unsigned int i = 0; i < n; ++i)
	  basestat.at(i) = base.variables.at(i);

	for (unsigned int j = 0; j < m; j++)
	  basestat.at(n + j) = base.constraints.at(j);


	// Basic variables
	/*static*/ std::vector<int> B(m);
	if (B.capacity() < m) { B.reserve(m); B.resize(m); }

	// Binv(i) -> i is  Binv(i)-th basic variable; -1 if non-basic
	/*static*/ std::vector<int> Binv(n + m, -1);
	if (Binv.capacity() < n + m) { Binv.reserve(n + m); Binv.resize(n + m); }

	//Non-Basic variables
	/*static*/ std::vector<int> NN(n);
	if (NN.capacity() <n) { NN.reserve(n); NN.resize(n); }

	// Ninv(i) -> i is  Ninv(i)-th non-basic variable; -1 if basic
	/*static*/ std::vector<int> Ninv(n + m, -1);
	if (Ninv.capacity() < n + m) { Ninv.reserve(n + m); Ninv.resize(n + m); }
	for (int z = 0; z < n + m; z++){
	  Ninv[z] = -1;
	  Binv[z] = -1;
	}

	/*
	 * An array. The array contains the indices of the variables in the resident basis,
	 * where basic slacks are specified by -rowindex - 1.
	 */
#ifdef USE_NBD_CPLEX_C
	std::vector<CPXDIM> head(m);
	if (head.capacity()<m) { head.reserve(m); head.resize(m); }

	if (CPXXgetbhead(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), head.data(), NULL)){
	  throw utils::ExternSolverException("Unable to obtain basis.");
	}

	//detect and store basic variables
	for (unsigned int i = 0; i < m; ++i){
	  if (head.at(i) >= 0)
	    B[i] = head.at(i);
	  else
	    B[i] = n - (head.at(i) + 1);

	  Binv[B[i]] = i;
	}
#endif
#ifdef USE_NBD_CLP
          ClpSimplex *M;
          M = (ClpSimplex*)extSolver.getSolverModel();
          
          if (M->problemStatus() || M->algorithm()>=0) {
	    //std::cerr << "Warning: GMI in trouble LP: Stat:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
	    //M->setAlgorithm(-1);
	    // M->dual(0,7);

	      if (M->problemStatus() || M->algorithm()>=0) {
		//std::cerr << "Warning: GMI in trouble: LP:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
              M->setAlgorithm(-1);
              M->dual(0,7);
  
                if (M->problemStatus() || M->algorithm()>=0) {
		  std::cerr << "Warning A: GMI in trouble: LP:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
		  M->allSlackBasis(true);   // reset basis
		  M->dual(0,7);
                }
                if (M->problemStatus() || M->algorithm()>=0) {
		  std::cerr << "Warning B: GMI in trouble: LP:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
		  M->allSlackBasis(true);   // reset basis
		  M->primal(0,7);
                }
                if (M->problemStatus() || M->algorithm()>=0) {
		  std::cerr << "Warning C: GMI in trouble: LP:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
		  M->allSlackBasis(true);   // reset basis
		  M->primal(0,3);
                }
	      }	
              


	      if (M->problemStatus() || M->algorithm()>=0) 
                std::cerr << "Error: GMI in trouble LP: Stat:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
              
              ///*static*/ extSol::QpExternSolver::QpExtSolBase base;
              extSolver.getBase(base);
	      if (base.variables.size() < n) {
		cerr << "Error: got no base" << endl;
		goto NoCuts;
	      }

              for (unsigned int i = 0; i < n; ++i)
                  basestat.at(i) = base.variables.at(i);
              
              for (unsigned int j = 0; j < m; j++)
                  basestat.at(n + j) = base.constraints.at(j);
              
              
              // Basic variables
              ///*static*/ std::vector<int> B(m);
              if (B.capacity() < m) { B.reserve(m); B.resize(m); }
              
              // Binv(i) -> i is  Binv(i)-th basic variable; -1 if non-basic
              ///*static*/ std::vector<int> Binv(n + m, -1);
              if (Binv.capacity() < n + m) { Binv.reserve(n + m); Binv.resize(n + m); }
              
              //Non-Basic variables
              ///*static*/ std::vector<int> NN(n);
              if (NN.capacity() <n) { NN.reserve(n); NN.resize(n); }
              
              // Ninv(i) -> i is  Ninv(i)-th non-basic variable; -1 if basic
              ///*static*/ std::vector<int> Ninv(n + m, -1);
              if (Ninv.capacity() < n + m) { Ninv.reserve(n + m); Ninv.resize(n + m); }
              for (int z = 0; z < n + m; z++){
                  Ninv[z] = -1;
                  Binv[z] = -1;
              }

              
              if (M->problemStatus() || M->algorithm()>=0) {
                  std::cerr << "Warning LP II: GMI in trouble: Stat:" << M->problemStatus() << " Alg:" << M->algorithm() << std::endl;
                  M->allSlackBasis(true);   // reset basis
                  //extSolver.prepareMatrixRowForm();
                  
                  M->initialSolve();
                  //ClpDualRowSteepest dualSteepest;
                  //M->setDualRowPivotAlgorithm(dualSteepest);
                  //M->allSlackBasis(true);   // reset basis
                  //M->dual(0,1);
                  if (M->problemStatus() || M->algorithm()>=0) {
                      std::cerr << "Error: LP not allowed." << std::endl;
		      goto NoCuts;//return cuts;
                  }
              }
          }
          
          //CoinWarmStartBasis *wsb = M->getBasis();
          std::vector<int> head(m);
          if (head.capacity()<m) { head.reserve(m); head.resize(m); }
          //std::cerr << "H";
          M->dual(0,7);
	M->getBasics(&(head.data()[0]));
	for (unsigned int i = 0; i < m; ++i){
	  B.at(i) = head.at(i);
	  Binv.at(B.at(i)) = i;
	}
#endif

	int ntmp = 0;
	for (unsigned int i = 0; i <n + m; i++){
	  if (Binv[i] == -1){
	    NN[ntmp] = i;
	    Ninv[i] = ntmp;
	    ntmp++;
	  }
	}
	assert(ntmp == n);
	if (Binv[k] != -1){//Only if the candidate variable is basic (should always be true)
	  // Vector equal to the row of the basis inverse regarding variable k
	  std::vector<double> Lambda(m);
	  std::vector<double> s_vec(n);
	  std::vector<double> t_vec(n);
	  std::vector<double> Tableau_row(n + m, 0);

#ifdef USE_NBD_CPLEX_C
	  if (CPXXbinvrow(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), Binv[k], Lambda.data()))
	    throw utils::ExternSolverException("CPXXbinvrow error");

	  if (CPXXbinvarow(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), Binv[k], Tableau_row.data()))
	    throw utils::ExternSolverException("CPXXbinvarow error");

	  for (unsigned int i = 0; i < m; ++i){
	    assert(Tableau_row[n + i] == 0);
	    //// logische Variablen
	    if (Ninv[n + i] != -1){
	      Tableau_row[n + i] = -Lambda[i];

	      /*if (SaveRhsVec.at(i).getRatioSign() == data::QpRhs::greaterThanOrEqual)
		Tableau_row[n + i] = -Lambda[i];
		else if (SaveRhsVec.at(i).getRatioSign() == data::QpRhs::smallerThanOrEqual)
		Tableau_row[n + i] = -Lambda[i];
		else {
		if (SaveRhsVec.at(i).getRatioSign() != data::QpRhs::equal) goto NextCandidate;
		Tableau_row[n + i] = -Lambda[i];
		}*/
	    }
	  }
#endif
#ifdef USE_NBD_CLP
	  //AUFPASSEN!!!!
	  //Get the whole row of the tableau
	  M->getBInvARow(Binv[k], &(Tableau_row.data()[0]), &(Tableau_row.data()[n]));
	  M->getBInvRow(Binv[k], &(Lambda.data()[0]));
	  for (unsigned int i = 0; i < m; ++i){

	    //assert(Tableau_row[n + i] == 0);
	    //// logische Variablen
	    //cerr << i << " " << Tableau_row[n+i] << endl;
	    if (Ninv[n + i] != -1){
	      assert(Tableau_row[n + i] == Lambda[i]);
	      Tableau_row[n + i] = -Lambda[i];
	    }

	    // HERE WE HAVE TO CHEAT A BIT, SINCE CLP's getBase DOES SOMETHING UNEXPECTED...
	    if (abs(Tableau_row[n+i])<epsNull) Tableau_row[n+i]=0;
	    else if (Tableau_row[n+i]>0)basestat.at(n+i)=extSol::QpExternSolver::AtLower;
	    else basestat.at(n+i)=extSol::QpExternSolver::AtUpper;


	  }
	  for (int i=0;i<n;i++){
	    if (abs(Tableau_row[i])<epsNull) Tableau_row[i]=0;
	    else if (Tableau_row[i]>0)basestat.at(i)=extSol::QpExternSolver::AtLower;
	    else basestat.at(i)=extSol::QpExternSolver::AtUpper;

	    if(M->getColumnStatus(i)==5){
	      if (Tableau_row[i]>=0)basestat.at(i)=0;
	      else basestat.at(i)=2;
	    }

	    //if(M->getColumnStatus(i)==2){
	    //Tableau_row[i] *=-1;
	    //if(M->getColumnStatus(i)==3) basestat.at(i)=2;
	    //cerr << i << " STATSCOL " <<basestat.at(i) << " "<< (int)M->getColumnStatus(i)<<endl;
	    //}
	  }
	  //cin.get();


	  extSolver.prepareMatrixRowForm();
#endif
	  std::vector<double> u_vec(m);
	  std::vector<double> v_vec(m);

	  for (int j = 0; j < m; j++){
	    if (Ninv[n + j] == -1){//i ist Basis-Variable
	      assert(Binv[n + j] != -1);
	      u_vec[j] = 0;
	      v_vec[j] = 0;
	    }
	    else{//j is non basic
	      if (basestat.at(n + j) == extSol::QpExternSolver::AtLower){
		if (Tableau_row[n + j] > 0){// This should be the case
		  v_vec[j] = Tableau_row[n + j];
		  u_vec[j] = 0;
		}
		else if (abs(Tableau_row[n + j]) < epsNull){// Counts as zero
		  u_vec[j] = 0;
		  v_vec[j] = 0;
		}
		else{   //Too dangerous
		  if (info_level >= 2)  cerr << Tableau_row[n + j] << "Warning: Tableau Entry negative!!" << endl;// 1e-10
		  //cin.get();
		  goto NextCandidate;
		}

	      }
	      else if (basestat.at(n + j) == extSol::QpExternSolver::AtUpper){
		if (Tableau_row[n + j] < 0){// This should be the case!
		  v_vec[j] = 0;
		  u_vec[j] = -Tableau_row[n + j];
		}
		else if (abs(Tableau_row[n + j]) < epsNull){// Counts as Zero
		  u_vec[j] = 0;
		  v_vec[j] = 0;
		}
		else{
		  if (info_level >= 2) cerr << Tableau_row[n + j] << "Warning: Tableau Entry positive!!" << endl;// Too dangerous
		  //cin.get();
		  goto NextCandidate;
		}
	      }
	      else{
		u_vec[j] = 0;
		v_vec[j] = 0;
	      }
	    }
	  }

	  /*for (int t = 0; t < m; t++){
	    if ((Lambda[t] != 0 || (u_vec[t] - v_vec[t]) != 0) && abs(Lambda[t] -(u_vec[t] - v_vec[t]))>1e-9){
	    cerr << "Lambda " << Lambda[t] << " " << u_vec[t] - v_vec[t] << endl;
	    goto NextCandidate;
	    }
	    }*/


        for (int j = 0; j < n; j++){
            if (Ninv[j] == -1){//i ist Basis-Variable
                assert(Binv[j] != -1);
                s_vec[j] = 0;
                t_vec[j] = 0;
            }
            else{//j ist Nicht-Basis-Variable
                if (basestat.at(j) == extSol::QpExternSolver::AtUpper){
                    if (Tableau_row[j] <= 0){
                        s_vec[j] = -Tableau_row[j];
                        t_vec[j] = 0;
                    }
                    else if (abs(Tableau_row[j]) < epsNull){
                        s_vec[j] = 0;
                        t_vec[j] = 0;
                    }
                    else{
                        if (info_level >= 2) cerr << Tableau_row[j] << "Warning: Tableau entry positive" << endl;
                        goto NextCandidate;
                    }
                }
                else{
                    if (basestat.at(j) != extSol::QpExternSolver::AtLower) {
                        cerr << "Error: LP cut, undefined status." << endl;
                        goto NextCandidate;
                    }
                    assert(basestat.at(j) == extSol::QpExternSolver::AtLower);
                    if (Tableau_row[j] >= 0){
                        s_vec[j] = 0;
                        t_vec[j] = Tableau_row[j];
                    }
                    else if (abs(Tableau_row[j]) < epsNull){
                        s_vec[j] = 0;
                        t_vec[j] = 0;
                    }
                    else
                    {
                        if (info_level >= 2)  cerr << Tableau_row[j] << "Warning: Tableau Entry negative" << endl;
                        goto NextCandidate;
                    }
                    
                }
                
            }
        }

	  //Siehe Beweis Lemma 3 Bonami: Erforderlich? Nicht-Basis-Variablen, die nicht nicht 0 sind fehlen vielleicht?
	  double b_Lambda = 0;
	  for (int s = 0; s < m; s++){
	    b_Lambda += Lambda[s] * SaveRhsVec[s].getValue().asDouble();

	    //if (SaveRhsVec.at(s).getRatioSign() == data::QpRhs::greaterThanOrEqual)
	    //b_Lambda += Lambda[s] * SaveRhsVec[s].getValue().asDouble();
	    //else if (SaveRhsVec.at(s).getRatioSign() == data::QpRhs::smallerThanOrEqual)
	    //b_Lambda += Lambda[s] * SaveRhsVec[s].getValue().asDouble();//+, since Lambda AND the RHS-value are turned around
	    //else {
	    //if (SaveRhsVec.at(s).getRatioSign() != data::QpRhs::equal) goto NextCandidate;
	    //b_Lambda += Lambda[s] * SaveRhsVec[s].getValue().asDouble();//SAVE
	    //}
	  }

	  // Generate Cut s'x +(b'Lambda-1)x_k>=0
	  vector<double> cut(n, 0);
	  double cut_rhs = 0;

	  double u_0 = -b_Lambda + 1;
	  double v_0 = 1 - u_0;
	  std::vector<double> uTA(n);
	  std::vector<double> vTA(n);
	  std::vector<double> em(n);


	  if (abs(u_0 - (1 - x_star[k].asDouble())) > epsNull){
	    //Important??????
	    //goto NextCandidate;
	  }
	  if (u_0 >= 1 || u_0 <= 0){
	    goto NextCandidate;
	  }
	  assert(s_vec[k] == 0);

	  for (int i = 0; i < n; i++){
	    cut[i] = s_vec[i];
	    uTA[i] = 0;
	    vTA[i] = 0;
	    em[i] = 0;
	  }
	  cut[k] = b_Lambda - 1;

	  for (int j = 0; j < m; j++){
	    if (u_vec[j] > 0 || v_vec[j] > 0){
	      assert(Ninv[n + j] != -1);
	      std::vector<data::IndexedElement> A_row;
	      extSolver.getRowLhs(j, A_row);
	      if (A_row.size() <= 0) goto NextCandidate;
	      for (int l = 0; l < A_row.size(); ++l){ //-1, da A=(A,I)
		em[A_row[l].index] += Lambda[j] * A_row[l].value.asDouble();
		if (SaveRhsVec.at(j).getRatioSign() == data::QpRhs::greaterThanOrEqual){
		  uTA[A_row[l].index] += u_vec[j] * A_row[l].value.asDouble();
		  vTA[A_row[l].index] += v_vec[j] * A_row[l].value.asDouble();
		}
		else if (SaveRhsVec.at(j).getRatioSign() == data::QpRhs::smallerThanOrEqual){
		  uTA[A_row[l].index] -= v_vec[j] * A_row[l].value.asDouble();
		  vTA[A_row[l].index] -= u_vec[j] * A_row[l].value.asDouble();
		}
		else {
		  if (SaveRhsVec.at(j).getRatioSign() != data::QpRhs::equal) goto NextCandidate;
		  uTA[A_row[l].index] += u_vec[j] * A_row[l].value.asDouble();
		  vTA[A_row[l].index] += v_vec[j] * A_row[l].value.asDouble();
		}

		//if (A_row.at(l).index == k || types[j] != 0){
		if (SaveRhsVec.at(j).getRatioSign() == data::QpRhs::greaterThanOrEqual)
		  cut[A_row.at(l).index] += u_vec[j] * A_row.at(l).value.asDouble();
		else if (SaveRhsVec.at(j).getRatioSign() == data::QpRhs::smallerThanOrEqual)
		  cut[A_row.at(l).index] -= v_vec[j] * A_row.at(l).value.asDouble();
		else {
		  if (SaveRhsVec.at(j).getRatioSign() != data::QpRhs::equal) goto NextCandidate;
		  cut[A_row.at(l).index] += u_vec[j] * A_row.at(l).value.asDouble();

		}
		//}
	      }
	      if (SaveRhsVec.at(j).getRatioSign() == data::QpRhs::greaterThanOrEqual)
		cut_rhs += u_vec[j] * (SaveRhsVec[j].getValue().asDouble());
	      else if (SaveRhsVec.at(j).getRatioSign() == data::QpRhs::smallerThanOrEqual)
		cut_rhs -= v_vec[j] * (SaveRhsVec[j].getValue().asDouble());
	      else {
		if (SaveRhsVec.at(j).getRatioSign() != data::QpRhs::equal) goto NextCandidate;
		cut_rhs += u_vec[j] * (SaveRhsVec[j].getValue().asDouble());
	      }
	    }
	  }

	  for (int j = 0; j<n; j++){
	    if (j != k && types[j] == 0){
	      //if (cut[j] - max(uTA[j], vTA[j]) >= epsNull)
	      //cerr << cut[j] << " " << uTA[j] << " " << vTA[j] << endl;
	      //assert(cut[j] - max(uTA[j], vTA[j])<epsNull);
              if (abs(em[j] - (uTA[j] - vTA[j])) >= epsNull) {
		  std::cerr << "Warning: abs(em[j] - (uTA[j] - vTA[j])) < epsNull. Under control" << std::endl;
                  goto NoCuts;
	      }
	      assert(abs(em[j] - (uTA[j] - vTA[j])) < epsNull);
	      double Left = uTA[j] - u_0*floor(em[j]);
	      double Right = vTA[j] + v_0*ceil(em[j]);
	      if (Left < Right) cut[j] = Left;
	      else cut[j] = Right;
	    }
	  }

	  // Old (spearated) Strengthening; Now already happens before

	  //double u_0 = -b_Lambda + 1;
	  //double v_0 = 1 - u_0;
	  //if (abs(u_0 - (1 - x_star[k].asDouble()))>epsNull){
	  ////cerr << "U0 Error " << u_0 << " " << 1 - x_star[k].asDouble() << " " << x_star[k].asDouble() << endl;
	  //goto NextCandidate;
	  //}
	  ////WICHTIG!??!?!? u_0==1-x_star[k] nur manchmal...
	  //if (u_0 >= 1 || u_0 <= 0){
	  ////vector<data::QpNum> y_star(n);
	  ////extSolver.getValues(y_star);
	  ////cerr << "u0 " << u_0 << " " << (1 - x_star[k].asDouble()) << " " << x_star[k] <<  endl;
	  //goto NextCandidate;
	  //}
	  ////assert(u_0<1 && u_0>0);
	  //std::vector<double> uTA(n);
	  //std::vector<double> vTA(n);
	  //std::vector<double> em(n);
	  //for (int j = 0; j < n; j++){
	  //uTA[j] = 0;
	  //vTA[j] = 0;
	  //em[j] = 0;
	  //}
	  //for (int row = 0; row < m; row++){
	  //if (Lambda[row] != 0 || u_vec[row] != 0 || v_vec[row] != 0){//Lambda m�sste reichen!?
	  //std::vector<data::IndexedElement> Arow;
	  //extSolver.getRowLhs(row, Arow);
	  //for (int j = 0; j < Arow.size(); j++){
	  //em[Arow[j].index] += Lambda[row] * Arow[j].value.asDouble();
	  //if (SaveRhsVec.at(row).getRatioSign() == data::QpRhs::greaterThanOrEqual){
	  //uTA[Arow[j].index] += u_vec[row] * Arow[j].value.asDouble();
	  //vTA[Arow[j].index] += v_vec[row] * Arow[j].value.asDouble();
	  //}
	  //else if(SaveRhsVec.at(row).getRatioSign() == data::QpRhs::smallerThanOrEqual){
	  //uTA[Arow[j].index] -= v_vec[row] * Arow[j].value.asDouble();
	  //vTA[Arow[j].index] -= u_vec[row] * Arow[j].value.asDouble();
	  //}
	  //else {
	  //if (SaveRhsVec.at(row).getRatioSign() != data::QpRhs::equal) goto NextCandidate;
	  //uTA[Arow[j].index] -= v_vec[row] * Arow[j].value.asDouble();
	  //vTA[Arow[j].index] -= u_vec[row] * Arow[j].value.asDouble();

	  //}
	  //}
	  //}
	  //}

	  //for (int j = 0; j<n; j++){
	  //if (j != k && types[j] == 0){         sNull)
	  //cerr << cut[j] << " " << uTA[j] << " " << vTA[j] << endl;
	  //assert(cut[j] - max(uTA[j], vTA[j])<epsNull);
	  //assert(abs(em[j] - (uTA[j] - vTA[j])) < epsNull);
	  //double Left = uTA[j] - u_0*floor(em[j]);
	  //double Right = vTA[j] + v_0*ceil(em[j]);
	  //if (Left < Right) cut[j] = Left;
	  //else cut[j] = Right;
	  //}
	  //}

	  vector<pair<unsigned int, double> > cutvec;
	  double max_c = 0;
	  for (unsigned int i = 0; i < n; ++i){
	    if (abs(cut.at(i)) > max_c){
	      max_c = abs(cut.at(i));
	    }
	  }
	  if (max_c == 0){
	    if (info_level >= 2) cerr << max_c << " max_c" << endl;
	    goto NextCandidate;
	  }

	  for (unsigned int i = 0; i < n; ++i){
	    if (abs(cut[i]) > 1e-20) {
	      if (abs(cut.at(i)) > max_c * 1e-9 || types[i] != 0)
		cutvec.push_back(make_pair(i, cut[i]));
	      else{
		//Not so sure if this is ALWAYS correct... especially for non-binaries
		cut_rhs += fmax(0, -cut.at(i));// -= cut.at(i);//-=abs(cut[i]); //
	      }
	    }
	  }

	  //double LHS_Val = 0;
	  //for (int v = 0; v < cutvec.size(); v++)
	  //LHS_Val += cutvec[v].second *x_star[cutvec[v].first].asDouble();

	  //cut_rhs += -OneEM12 - abs(cut_rhs)*OneEM12;
	  if (cutvec.size() > 0) { // && LHS_Val < cut_rhs)
	    cuts.push_back(make_pair(cutvec, cut_rhs));
	    listOfCutsVars.push_back(candidateVariables[k_ind]);
	  }
	  //else cerr << "No Cut " << LHS_Val << " " << cut_rhs << endl;
	}
	else{
	  // Theoretically not possible: Bonami, Lemma 3
	  if (info_level >= 2) cerr << "Candidate variable not basic in MLP" << endl;
	}

      }
    NextCandidate:
      extSolver.changeObjFuncCoeff(k, dummyZero);
      //if(cuts.size()>sqrt(n)) break;
    }

  NoCuts:
    // LP zur�ck bauen
#ifdef USE_NBD_CPLEX_C
    for (int j = 0; j < m; j++){
      SaveRhsVals[j] = SaveRhsVec[j].getValue();
      switch (SaveRhsVec.at(j).getRatioSign()){
      case data::QpRhs::smallerThanOrEqual:
	SenseRow[j] = 'L';
	//CPXXchgcoef(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), j, -1, SaveRhsVals[j].asDouble());

	break;
      case data::QpRhs::equal:
	SenseRow[j] = 'E';
	break;
      case data::QpRhs::greaterThanOrEqual:
	SenseRow[j] = 'G';
	break;
      default:
	throw runtime_error("Unexpected sense.");
	break;
      }
    }
    //Sense und RHS
    CPXXchgsense(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), m, RowCount.data(), SenseRow.data());

    extSolver.changeRhsElements(uRowCount, SaveRhsVals);
    CPXXchgobjsen(*(CPXENVptr*)extSolver.getSolverEnv(), *(CPXLPptr*)extSolver.getSolverModel(), SaveObjSense);

#endif
#ifdef USE_NBD_CLP
    for (int j = 0; j < m; j++){
      switch (SaveRhsVec.at(j).getRatioSign()){
      case data::QpRhs::smallerThanOrEqual:
	M->setRowBounds(j,- COIN_DBL_MAX,SaveRhsVec[j].getValue().asDouble());
	break;

      case data::QpRhs::equal:
	M->setRowBounds(j,SaveRhsVec[j].getValue().asDouble(),SaveRhsVec[j].getValue().asDouble());
	break;

      case data::QpRhs::greaterThanOrEqual:
	M->setRowBounds(j,SaveRhsVec[j].getValue().asDouble(),COIN_DBL_MAX);
	break;

      default:
	throw runtime_error("Unexpected sense.");
	break;
      }
    }
    M->setOptimizationDirection(SaveObjSense); 
#endif
    for (int i = 0; i < n; i++){
      extSolver.setVarLB(i, SaveLbVec[i]);
      extSolver.setVarUB(i, SaveUbVec[i]);
      extSolver.changeObjFuncCoeff(i, SaveObjCoeffs[i]);
    }
    extSolver.solve(100000, 100000);

  }
  //cuts.clear();
  if (info_level >= 2) cerr << cuts.size() << " LP Cuts created" << endl;
  return cuts;
}

#endif


