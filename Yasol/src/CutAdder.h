/*
*
* Yasol: CutAdder.h -- Copyright (c) 2012-2014 Ulf Lorenz, Thomas Opfer
* Yasol: CutAdder.h -- Copyright (c) 2014-2017 Michael Hartisch, Ulf Lorenz, Thomas Opfer
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

#ifndef CUTADDER_H_
#define CUTADDER_H_
#ifdef WINDOWS
#include <time.h>
#endif
#include "ExternSolvers/QpExternSolver.hpp"


#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <algorithm>
#include <functional>
#include <unordered_map>

#include <assert.h>

#define DEBUG_OUT_EQ 0

class Param_c {
public:

};

#define fabs(x) ((x) >= 0 ? (x) : -(x))
class CutAdder {
public:
  static std::vector<std::pair<int, double>> globalConditionalLowerBounds;
  static std::vector<std::pair<int, double>> globalConditionalUpperBounds;

  static void printRowEQ(const std::vector<std::pair<unsigned int, double>>& lhs, double rhs, int *solu, int n, const double * xlp, const double* xlpExtr, int *types, bool invertRS=false) {
	  double lhs_lp=0.0;
	  double lhs_opt=0.0;
	  //if (!DEBUG_OUT_EQ) return;
	  //if (lhs.size() > 20) return;
	  for( unsigned int j = 0; j < lhs.size(); ++j ){
	    if (lhs[j].first < n) {
	      if (types[lhs[j].first] == 0)
		std::cerr << lhs[j].second << "x" << lhs[j].first << "(" << solu[lhs.at( j ).first]<< ") + ";
	      else
		std::cerr << lhs[j].second << "y" << lhs[j].first << "(" << solu[lhs.at( j ).first]<< ") + ";
	      lhs_lp = lhs_lp + lhs.at( j ).second * xlp[lhs.at( j ).first];
	      lhs_opt = lhs_opt + lhs.at( j ).second * solu[lhs.at( j ).first];
	    } else {
	      std::cerr << lhs[j].second << "Y" << lhs[j].first << "(" << xlpExtr[lhs.at( j ).first-n]<< ") + ";
	      lhs_lp = lhs_lp + lhs.at( j ).second * xlpExtr[lhs.at( j ).first-n];
	      lhs_opt = lhs_opt + lhs.at( j ).second * solu[lhs.at( j ).first];
	    }
	  }
	  std::cerr << "==" << lhs_opt << (invertRS ? " >= " : " <= ") << rhs << std::endl;
	  double lhs_o=0.0;
	  return;
	  std::vector<double> slacks;
	  if (solu != 0) {
	    for( unsigned int j = 0; j < lhs.size(); ++j ){
	      if (lhs[j].first < n)
		lhs_o = lhs_o + lhs.at( j ).second * solu[lhs.at( j ).first];
	      else {
		slacks.push_back(lhs.at( j ).second);
	      }
	    }
	    assert(slacks.size()<=1);
	    double mult=1.0;
	    if (slacks.size()==1 && slacks[0] < 0.0)
	      mult=-1.0;
	    if (mult*lhs_o > mult*rhs) {
	      for( unsigned int j = 0; j < lhs.size(); ++j ){
		if (lhs[j].first < n) {
		  std::cerr << "x" << lhs.at( j ).first << " solu:" << solu[lhs.at( j ).first] << " coef=" << lhs.at( j ).second << std::endl;
		  //lhs_o = lhs_o + lhs.at( j ).second * solu[lhs.at( j ).first];
		}
	      }
	      std::cerr << "lhs_o=" << lhs_o << " <!= " << rhs << std::endl;
	      assert(0);
	    }
	  }
	  return;
        }
  static double debugCheckSol(const std::vector<std::pair<unsigned int, double>>& lhs, double rhs, int *solu, int n, bool invert=false) {
          if (solu==0) return -1e100;
	  double lhs_o=0.0;
	  for( unsigned int j = 0; j < lhs.size(); ++j ){
	    if (lhs.at( j ).first<n) 
	      lhs_o = lhs_o + lhs.at( j ).second * solu[lhs.at( j ).first];
	    else
	      lhs_o = lhs_o + lhs.at( j ).second * 0;
	  }
	  if (invert) {
	    lhs_o = -lhs_o;
	    rhs = -rhs;
	  }
	  if (lhs_o > rhs + 1e-7) {
	    for( unsigned int j = 0; j < lhs.size(); ++j ){
	      if (lhs.at( j ).first<n) 	      
		std::cerr << "x" << lhs.at( j ).first << " solu:" << solu[lhs.at( j ).first] << " coef=" << lhs.at( j ).second << std::endl;
	      else
		;		
	    }
	    std::cerr << "lhs_o=" << lhs_o << " <!= " << rhs << std::endl;
	    assert(0);
	  }
	  return lhs_o;
        }
        static bool isCutViolated(const std::vector<std::pair<unsigned int, double>>& cut, double rhs, const std::vector<double>& xlpopt) {
	  double act = 0.0;
	  double const tol = 1e-6;
	  for (int i=0; i < cut.size();i++) {
	    unsigned int j = cut[i].first;
	    double coef = cut[i].second;
	    act += coef * xlpopt[j];
	  }
	  return act > rhs + 10 * tol + tol*fabs(fabs(rhs)+fabs(act));
	}
        static uint64_t approximateDenominator(double frac, double tol, int maxIter = 1000) {
	  double x = frac;
	  uint64_t h1 = 1, h2 = 0;
	  uint64_t k1 = 0, k2 = 1;

	  for (int i = 0; i < maxIter; ++i) {
	    int a = (int)std::floor(x);
	    uint64_t h = a * h1 + h2;
	    uint64_t k = a * k1 + k2;
	    if (std::abs(h / (double)k - frac) < tol)
	      return k;
	    h2 = h1; h1 = h;
	    k2 = k1; k1 = k;
	    x = 1.0 / (x - a);
	    if (k > 1e12) break; // Notfalls abbrechen
	  }
	  return 1;
	}
        // Hauptfunktion: finde gemeinsame Skalierung
        static double integralScale(const std::vector<double>& vals,
				    double deltadown = 1e-6,
				    double deltaup = 1e-13) {
	  if (vals.empty()) return 0.0;

	  // Finde kleinsten und größten Absolutwert
	  auto minmax = std::minmax_element(vals.begin(), vals.end(),
					    [](double a, double b) { return std::abs(a) < std::abs(b); });
	  double minval = std::abs(*minmax.first);
	  double maxval = std::abs(*minmax.second);

	  int expshift = 0;
	  if (minval > deltaup) std::frexp(minval, &expshift);
	  expshift = std::max(-expshift, 0) + 3;

	  int expMaxVal = 0;
	  std::frexp(maxval, &expMaxVal);
	  expMaxVal = std::min(expMaxVal, 32);
	  if (expMaxVal + expshift > 32) expshift = 32 - expMaxVal;

	  uint64_t denom = uint64_t{75} << expshift;
	  uint64_t startDenom = denom;

	  // Skaliere erstes Element und prüfe Rundung
	  double val = vals[0] * denom;
	  double downval = std::floor(val + deltaup);
	  double frac = val - downval;

	  if (frac > deltadown) {
	    denom *= approximateDenominator(frac, deltaup);
	    val = vals[0] * denom;
	    downval = std::floor(val + deltaup);
	    frac = val - downval;
	    if (frac > deltadown) return 0.0;
	  }

	  int64_t currgcd = std::llabs((int64_t)downval);

	  for (size_t i = 1; i < vals.size(); ++i) {
	    val = vals[i] * denom;
	    downval = std::floor(val + deltaup);
	    frac = val - downval;

	    if (frac > deltadown) {
	      denom = startDenom;
	      val = vals[i] * denom;
	      frac = val - std::floor(val);
	      denom *= approximateDenominator(frac, deltaup);
	      val = vals[i] * denom;
	      downval = std::floor(val + deltaup);
	      frac = val - downval;
	      if (frac > deltadown) return 0.0;
	    }

	    if (currgcd != 1)
	      currgcd = Gcd(currgcd, (int64_t)downval);
	  }

	  return denom / (double)currgcd;
	}
	static int64_t Gcd(int64_t a, int64_t b) {
		  int64_t c;
		  if ( a < 0 ) a = - a;
		  if ( b < 0 ) b = - b;
		  if ( a < b ) { c = a; a = b; b = c; }
		  while ( b != 0 ) {
		    c = a % b; a = b; b = c;
		  }
		  return(a);
	}
    static double getMirCoefficient(double d, double f, double eps) {
      double delta = d - floor(d) - f;
      if (delta > eps)
          return floor(d) + delta / (1 - f);
      else
          return floor(d);
    }
    static bool gotoLowerSubst(const double aj, const double xlp,  const double lb,  const double ub)
    {
        return xlp - lb < ub - xlp;
    }
    
    static double coeffNorm2(const std::vector<std::pair<unsigned int, double> >& lhs) {
        long double s = 0.0L;
        for (const auto& kv : lhs) {
            const long double a = kv.second;
            s += a * a;
        }
        const long double n = std::sqrt(s);
        return static_cast<double>(n);
    }
    
    // Aktivität a^T x
    static inline double activity(const std::vector<std::pair<unsigned int, double> >& lhs, const double* x, int n) {
      double s = 0.0;
      for (const auto& kv : lhs) {
	const unsigned j = kv.first;
	const long double a = kv.second;
	const long double xj = (x && j < (unsigned)n) ? x[j] : 0.0L;
	s += a * xj;
      }
      return s;
    }
    static void computeCmirViolationAndEfficacy(const std::vector<std::pair<unsigned int,double>>& lhs_cMIR, double rhs_cMIR, double sCoef, double sStar, const double* xlp, int n, double& viol, double& effi, double eps) {
      const double act = activity(lhs_cMIR, xlp, n);
      double v = (act - sCoef * sStar) - rhs_cMIR;   // ≤-Ungleichung
      if (v < eps) v = 0.0;

      const double lhsNorm = coeffNorm2(lhs_cMIR);
      const double n2 = std::sqrt( lhsNorm*lhsNorm + sCoef*sCoef );
      effi = (n2 > eps) ? (v / n2) : 0.0;
      viol = v / n2;
    }

  static void buildcMirInequality(double delta, double rhs,const std::vector<std::pair<unsigned int,double>>& rowIntpart, const data::QpNum* colUpperBound, const std::vector<bool>& isComplemented, std::vector<std::pair<unsigned int,double>>& lhs_cMIR,double& rhs_cMIR,double& sCoef,int n) {
    const double eps = 1e-9;
    lhs_cMIR.clear();
    
    const double scalrhs = rhs / delta;
    double f = scalrhs - std::floor(scalrhs);
    //if (f < eps) f = eps;
    //if (f > 1.0 - eps) f = 1.0 - eps;
    rhs_cMIR = std::floor(scalrhs);

    assert(isComplemented.size() == rowIntpart.size());
    
    for (int i = 0; i < rowIntpart.size(); ++i) {
      const unsigned int j = rowIntpart[i].first;
      const double a = rowIntpart[i].second;
      const bool jIsCompl = isComplemented[i];
        
      if (!jIsCompl) {
	const double mc = getMirCoefficient(a / delta, f, 1e-6);
	if (std::abs(mc) > 1e-12) lhs_cMIR.emplace_back(j, mc);
      } else {
	const double mc = getMirCoefficient(-a / delta, f, 1e-6);
	if (std::abs(mc) > 1e-12) lhs_cMIR.emplace_back(j, -mc);
	const double Uj = (j < (unsigned)n && colUpperBound) ? colUpperBound[j].asDouble() : 0.0;
	rhs_cMIR -= mc * Uj;
      }
    }
    
    // s-Koeffizient
    sCoef = 1.0 / (delta * (1.0 - f));
  }



  static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > > *getGMICuts( extSol::QpExternSolver &extSolver, std::vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev, unsigned int initime, std::vector<int> &listOfCutsVars , int*, int*, int*, int*, int);

  // MIRsmartEQ
  static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >* getMIRsmartEQ( extSol::QpExternSolver &extSolver, int *types, int8_t *assigns, int decLev, unsigned int initime, int * solu, int* fixs, int *blcks, int *eass, int orgN ) {
    unsigned int m = extSolver.getRowCount();
    unsigned int n = extSolver.getVariableCount();
    
    extSolver.prepareMatrixRowForm();
    std::vector<data::QpRhs> rhsVec( m );
    std::vector< std::vector< data::IndexedElement >  > allrows;
    bool success = getAllRows(extSolver,allrows);
    extSolver.getRhs( rhsVec );
    for (int i=0; i < allrows.size();i++) {
        if (rhsVec[i].getRatioSign() == data::QpRhs::greaterThanOrEqual ) {
	    for (int ii=0;ii<allrows[i].size();ii++) {
	      int idx=allrows[i][ii].index;
	      double coeff=allrows[i][ii].value.asDouble();
	      allrows[i][ii].value = -coeff;
	    }
	    rhsVec[i].setRatioSign(data::QpRhs::smallerThanOrEqual);
	    rhsVec[i].setValue(-rhsVec[i].getValue().asDouble());
	}
    }
    return getMIRsmartCore( extSolver, allrows, rhsVec, types, assigns, decLev, initime, solu, fixs, blcks, eass, orgN, 2, false);    
  }
  
  static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >* getMIRsmartCore( extSol::QpExternSolver &extSolver, std::vector<std::vector<data::IndexedElement>> &allRows, std::vector<data::QpRhs> &rhsVec, int *types, int8_t *assigns, int decLev, unsigned int initime, int * solu, int* fixs, int *blcks, int *eass, int orgN, int maxAggr, bool fromCmirize );
  static bool selectRowToAggregateEQ( const std::vector<std::vector<data::IndexedElement>> &allRows, const std::vector< std::pair<unsigned int, double> >& rowAggregated, const data::QpNum* colUpperBound, const data::QpNum* colLowerBound, const std::vector< std::pair<unsigned int, double> >& setRowsAggregated, const double* xlp, int& rowSelected, int& colSelected, std::vector<std::pair<int, double>>& tightenedLoBnds,std::vector<std::pair<int, double>> &tightenedUpBnds, int *types, int numCols , int orgN); 

  static void refreshVariableBndConstraintsEQ(std::vector<std::pair<int, double>>& tightenedLoBnds,std::vector<std::pair<int, double>> &tightenedUpBnds, const std::vector< std::vector< data::IndexedElement >  >& allRows, const std::vector<data::QpRhs>& rhsVec, int * types, int numCols, int orgN );

  static inline void boundSubstituteIntegerEQ(int j, double a_j, std::vector<std::pair<unsigned int,double>>& intAcc, std::vector<int>& intAccDense, const int* types){
    if (intAccDense[j] == -1) {
      intAccDense[j] = intAcc.size();
      intAcc.emplace_back(j, a_j);
    } else {
      intAcc[intAccDense[j]].second += a_j;
    }
  }
  static bool boundSubstituteContinuousEQ(int indCol, double coefCol, const double* xlp, const double* xlpSlacks, const data::QpNum* colUpperBound, const data::QpNum* colLowerBound, const std::vector<std::pair<int,double>>& tightenedLoBnds, const std::vector<std::pair<int,double>>& tightenedUpBnds, std::vector<std::pair<unsigned int,double>>& intAcc, std::vector<int>& intAccIndicesDense, std::vector<std::pair<unsigned int,double>> &contVariablesInS, std::vector<int>& contVariablesInSDense, double& rhsAfterSubst, double& sStar, int numCols, int *types, double infinity=1e15, double eps=1e-6){
    //assert(indCol < orgN || indCol >= numCols);
    if (indCol >= numCols) { //slack variable
      if (coefCol < -eps) {
	contVariablesInS.emplace_back(indCol, coefCol);
	// sStar += (-a_j) * (x - LB); hier x = xlpSlacks[indCol-numCols], LB=0
	sStar -= coefCol * xlpSlacks[indCol - numCols];
      }
      return true;
    }
    const double LoB = (tightenedLoBnds.size()>indCol&&tightenedLoBnds[indCol].first >= 0) ? tightenedLoBnds[indCol].second * xlp[tightenedLoBnds[indCol].first]
                                                 : colLowerBound[indCol].asDouble();
    const double UpB = (tightenedUpBnds.size()>indCol&&tightenedUpBnds[indCol].first >= 0) ? tightenedUpBnds[indCol].second * xlp[tightenedUpBnds[indCol].first]
                                                 : colUpperBound[indCol].asDouble();
    // Wenn beide Schranken "unendlich": kein cmir Element formbar 
    if (LoB <= -infinity && UpB >= infinity) return false;  

    if (gotoLowerSubst(coefCol, xlp[indCol], LoB, UpB)) {
        // tightenedLoBnds[.] -> Binder in Integer-Teil, sonst RHS-Shift
      if (tightenedLoBnds.size()>indCol&&tightenedLoBnds[indCol].first >= 0) {
	  if (fabs(coefCol * tightenedLoBnds[indCol].second) > 1e-12) {
	    if (intAccIndicesDense[tightenedLoBnds[indCol].first] == -1) {
	      intAccIndicesDense[tightenedLoBnds[indCol].first] = intAcc.size();
	      intAcc.emplace_back(tightenedLoBnds[indCol].first, coefCol * tightenedLoBnds[indCol].second);
	    } else {
	      intAcc[intAccIndicesDense[tightenedLoBnds[indCol].first]].second += coefCol * tightenedLoBnds[indCol].second;
	    }
	  }
        } else {
	  rhsAfterSubst -= coefCol * LoB;
        }
        // S-Term bei a_j < 0: +a_j contVariablesInS
	// Update sStar
        if (coefCol < -eps) {
	  if (fabs(coefCol) > 1e-12) {
	    if (contVariablesInSDense[indCol] == -1) {
	      contVariablesInSDense[indCol] = contVariablesInS.size();
	      contVariablesInS.emplace_back(indCol, coefCol);
	    } else {
	      contVariablesInS[contVariablesInSDense[indCol]].second = coefCol;
	    }
	  }
	  sStar -= coefCol * (xlp[indCol] - LoB);
        }
    } else {

      if (tightenedUpBnds.size()>indCol&&tightenedUpBnds[indCol].first >= 0 ) {
	  if (fabs(coefCol * tightenedUpBnds[indCol].second) > 1e-12) {
	    if (intAccIndicesDense[tightenedUpBnds[indCol].first] == -1) {
	      intAccIndicesDense[tightenedUpBnds[indCol].first] = intAcc.size();
	      intAcc.emplace_back(tightenedUpBnds[indCol].first, coefCol * tightenedUpBnds[indCol].second);
	    } else {
	      intAcc[intAccIndicesDense[tightenedUpBnds[indCol].first]].second += coefCol * tightenedUpBnds[indCol].second;
	    }
	  }
	}
	else {
	  rhsAfterSubst -= coefCol * UpB;
	}
	// Update sStar
	if (coefCol > eps) {
	  if (fabs(coefCol) > 1e-12) {
	    if (contVariablesInSDense[indCol] == -1) {
	      contVariablesInSDense[indCol] = contVariablesInS.size();
	      contVariablesInS.emplace_back(indCol, -coefCol);
	    } else {
	      contVariablesInS[contVariablesInSDense[indCol]].second = -coefCol;
	    }
	  }
	  sStar += coefCol * (UpB - xlp[indCol]);
	}
    }

    return true;
}
  static bool substituteBoundsAndSplitIntS_EQ( const std::vector<std::pair<unsigned int, double> >& rowIn, const double* xlp, const double* xlpExtra, const data::QpNum* colUpperBound, const data::QpNum* colLowerBound, std::vector<std::pair<unsigned int, double> >& rowIntPart, double& rhsAfterSubst, double& sStar, std::vector<std::pair<unsigned int, double> >& contVariablesInS, int numCols, int *types, const std::vector<std::pair<int, double>>& tightenedLoBnds, const std::vector<std::pair<int, double>>& tightenedUpBnds, int numSlacks ) {
  double infinity = 1e15;
  double eps=1e-6;
  int j;
  assert(rowIntPart.size()==0);
  assert(contVariablesInS.size()==0);

  std::vector<int> rowIntPartIndicesDense(numCols+numSlacks,-1);
  for (int ii=0;ii<rowIntPart.size();ii++)
    rowIntPartIndicesDense[rowIntPart[ii].first] = ii;
  std::vector<int> contVariablesInSDense(numCols+numSlacks,-1);
  for (int ii=0;ii<contVariablesInS.size();ii++)
    contVariablesInSDense[contVariablesInS[ii].first] = ii;

  for ( j = 0; j < rowIn.size(); ++j) {
    // get index and coefficient of column j in the given row
    const int indCol = rowIn[j].first;
    const double coefCol = rowIn[j].second;

    // if variable is mini-DEP variable, but not model variable in original QMIP
    //if (indCol >= orgN && indCol < numCols) return false;

    // if the lower bound is equal to the upper bound, remove variable
    if ( (indCol < numCols) &&
	 (fabs(colLowerBound[indCol].asDouble() - colUpperBound[indCol].asDouble()) < 1e-10) ) {
      rhsAfterSubst -= coefCol * colLowerBound[indCol].asDouble();
      continue;
    }

    // |a|≈0 -> relax on cost of rhs
    if (fabs(coefCol) < eps) {
      if (coefCol<0.0)
	rhsAfterSubst -= coefCol*colUpperBound[indCol].asDouble();
      else
	rhsAfterSubst -= coefCol*colLowerBound[indCol].asDouble();
      continue;
    }
    // set the coefficients of the integer variables
    if ( (indCol < numCols)  && (types[indCol]==0 /* BINARY */) ) {
      boundSubstituteIntegerEQ(indCol, coefCol, rowIntPart, rowIntPartIndicesDense, types);
      continue;
    }

    if (!boundSubstituteContinuousEQ(
            indCol, coefCol, xlp, xlpExtra, colUpperBound, colLowerBound, tightenedLoBnds, tightenedUpBnds,
            rowIntPart, rowIntPartIndicesDense,
            contVariablesInS, contVariablesInSDense,
            rhsAfterSubst,
            sStar,
            numCols, types)) {
      sStar = 0.0;
      rhsAfterSubst=0.0;
      rowIntPart.clear();
      contVariablesInS.clear();
      return false; 
    }
  }
    
  if (contVariablesInS.size() == 0) return false; 
  if (rowIntPart.size() == 0) return false;       
  for ( j = 0; j < rowIntPart.size(); ++j) {
    int indCol = rowIntPart[j].first;
    // if the coefficient is zero, disregard
    if (fabs(rowIntPart[j].second) < eps) continue;
    // if the lower bound is not zero, then return false
    if (fabs(colLowerBound[indCol].asDouble()) > eps) {
      return false;
    }
  }
  return true;
}
static double improvedUbOf(
    int j,
    const data::QpNum* globalUB,
    const std::vector<std::pair<int,double>>& tightenedUB,
    const double* xlp
){
    if (j >= 0 && j < (int)tightenedUB.size()) {
        const int idx = tightenedUB[j].first;
        const double c = tightenedUB[j].second;
        if (idx >= 0) return c * xlp[idx];
    }
    return globalUB[j].asDouble();
}

static double improvedLbOf(
    int j,
    const data::QpNum* globalLB,
    const std::vector<std::pair<int,double>>& tightenedLB,
    const double* xlp
){
    if (j >= 0 && j < (int)tightenedLB.size()) {
        const int idx = tightenedLB[j].first;
        const double c = tightenedLB[j].second;
        if (idx >= 0) return c * xlp[idx];
    }
    return globalLB[j].asDouble();
}
// re-formulate the best cut with the model variables  , i.e. without slacks
// expectation: bestCutIndicesDense has size (n + numSlacks) and is initialisied with -1.
static inline bool expandSTermThroughSlacks(
    const std::vector<std::pair<unsigned,double>>& sStarTerms, // (varIdx, coefInS)
    double sCoefBestCut,
    int n, int numSlacks,
    const std::vector<std::vector<data::IndexedElement>>& allRows,
    const std::vector<data::QpRhs>& rhsVec,
    const data::QpNum*                                    colUpperBound,
    const data::QpNum*                                    colLowerBound,
    const std::vector<std::pair<int,double>>&             tightenedLoBnds,
    const std::vector<std::pair<int,double>>&             tightenedUpBnds,
    const int* listRowsAggregated,
    const std::vector<std::pair<unsigned,double>>& inRowLhs, 
    const std::vector<double>& inRowDense, 
    std::vector<std::pair<unsigned,double>>& bestCutLhs, 
    std::vector<int>& bestCutIndicesDense,
    double& bestCutRhs,
    const double *xlp,
    std::pair<std::vector<std::pair<unsigned int,double>>,double>& cMirCut,    
    double tol = 1e-12
){
    for (int j = 0; j < sStarTerms.size(); ++j) {
        int varix = sStarTerms[j].first;
        double varco = sStarTerms[j].second;
        
        if (varix < n) {  // variable is model variable
            double LoB = improvedLbOf(varix, colLowerBound, tightenedLoBnds, xlp);
            double UpB = improvedUbOf(varix, colUpperBound, tightenedUpBnds,xlp);
	    std::pair<int, double> VLB = std::make_pair(-1,0.0);
            std::pair<int, double> VUB = std::make_pair(-1,0.0);
	    if (tightenedLoBnds.size() > varix) VLB = tightenedLoBnds[varix];
	    if (tightenedUpBnds.size() > varix) VLB = tightenedUpBnds[varix];
            // Select the bound substitution
            if (gotoLowerSubst(inRowDense[varix], xlp[varix], LoB, UpB)) {
	        if (fabs(sCoefBestCut * varco) > 1e-12) {		      
		  if (bestCutIndicesDense[varix] == -1) {
		    bestCutIndicesDense[varix] = bestCutLhs.size();
		    bestCutLhs.emplace_back(varix, sCoefBestCut * varco);
		  } else {
		    bestCutLhs[bestCutIndicesDense[varix]].second = sCoefBestCut * varco;
		  }		      
		}
                if (VLB.first >= 0) {
                    int indVLB = tightenedLoBnds[varix].first; 
		    if (fabs(- sCoefBestCut * varco * tightenedLoBnds[varix].second) > 1e-12) {
		      if (bestCutIndicesDense[indVLB] == -1) {
			bestCutIndicesDense[indVLB] = bestCutLhs.size();
			bestCutLhs.emplace_back(indVLB, - sCoefBestCut * varco * VLB.second);
		      } else {
			bestCutLhs[bestCutIndicesDense[indVLB]].second += (- sCoefBestCut * varco * VLB.second);
		      }		      
		    }
                } else {
		    bestCutRhs += sCoefBestCut * varco * colLowerBound[varix].asDouble();
                }
            } else {
	        if (fabs(- sCoefBestCut * varco) > 1e-12) {		      
		  if (bestCutIndicesDense[varix] == -1) {
		    bestCutIndicesDense[varix] = bestCutLhs.size();
		    bestCutLhs.emplace_back(varix, - sCoefBestCut * varco);
		  } else {
		    bestCutLhs[bestCutIndicesDense[varix]].second = - sCoefBestCut * varco;
		  }		      
		}
                if (VUB.first >= 0) {
                    int indVUB = VUB.first;
		    if (fabs(sCoefBestCut * varco * VUB.second) > 1e-12) {
		      if (bestCutIndicesDense[indVUB] == -1) {
			bestCutIndicesDense[indVUB] = bestCutLhs.size();
			bestCutLhs.emplace_back(indVUB, sCoefBestCut * varco * VUB.second);
		      } else {
			bestCutLhs[bestCutIndicesDense[indVUB]].second += sCoefBestCut * varco * VUB.second;
		      }		      
		    }
                } else {
		  bestCutRhs -= sCoefBestCut * varco * colUpperBound[varix].asDouble();
                }
            }
        }
        else {  // variable is slack variable
            // in this case the LB = 0 and the UB = infinity
            const int iRow = listRowsAggregated[varix - n];
            const double mult = (rhsVec.at(iRow).getRatioSign()==data::QpRhs::smallerThanOrEqual)
                                ? (- sCoefBestCut * varco)
                                : (  sCoefBestCut * varco);
            bestCutRhs += rhsVec.at(iRow).getValue().asDouble() * mult;
	    const std::vector<data::IndexedElement>& row = allRows[iRow];
            int nElements = row.size();
            for (int i=0;i<row.size();i++) {
	      if (fabs(row[i].value.asDouble()*mult) > 1e-12) {
		if (bestCutIndicesDense[row[i].index] == -1) {
		  bestCutIndicesDense[row[i].index] = bestCutLhs.size();
		  bestCutLhs.emplace_back(row[i].index, row[i].value.asDouble()*mult);
		} else {
		  bestCutLhs[bestCutIndicesDense[row[i].index]].second += row[i].value.asDouble()*mult;
		}		      
	      }
	    }
        }
    }
    // Check the violation of the cut with the original variables.
    double cutRHS = bestCutRhs;
    double violation = 0.0;
    double normCut = 0.0;
    double largest=0.0;

    for (int j = 0; j < bestCutLhs.size(); ++j) {
      double value = bestCutLhs[j].second;
	if (largest < fabs(value))
	  largest = fabs(value);
    }
    double testValue=fmax(1.0e-6*largest,1.0e-12);
    for (int j = 0; j < bestCutLhs.size(); ++j) {
        int column = bestCutLhs[j].first;
        double value = bestCutLhs[j].second;
        if (fabs(value)>testValue) {
            violation += value * xlp[column];
            normCut += value * value;
        } else if (value) {
	    bestCutLhs[j].second = 0.0;
            if (value>0.0) {
	      double modification = value*colLowerBound[column].asDouble();
	      if (colLowerBound[column].asDouble()>0.0) {
                    modification=0.0;
                }
                cutRHS -= modification;
            } else {
	      double modification = value*colUpperBound[column].asDouble();
	      if (colUpperBound[column].asDouble()<0.0) {
                    modification=0.0;
                }
                cutRHS -= modification;
            }
        }
    }
    for (int j=0; j < bestCutLhs.size();j++)
      if (bestCutLhs[j].second == 0.0) {
	bestCutLhs[j] = bestCutLhs[bestCutLhs.size()-1];
	bestCutLhs.pop_back();
	j--;
      }
    violation -= cutRHS;
    violation /= sqrt(normCut);
    
    if ( violation > tol ) {
        cMirCut.first.clear();
	for (int j=0; j < bestCutLhs.size();j++)
	  cMirCut.first.emplace_back(bestCutLhs[j]);
	cMirCut.second = cutRHS;
        bestCutLhs.clear();
        return true;
    } else {
        bestCutLhs.clear();
    }    
    return false;
}
static bool cMirSeparationEQ_core(
    const double*                                         xlp,
    double                                                sStar,
    const data::QpNum*                                    colUpperBound,
    const data::QpNum*                                    colLowerBound,
    const std::vector<std::pair<unsigned,double>>&        intPart,           
    const double&                                         rhsAfterSubst,   
    const std::vector<std::pair<unsigned,double>>&        sTerms,
    double&                                               sCoefBestCut,
    std::vector<std::pair<unsigned,double>>&              cMirCutLhs,
    double&                                               cMirCutRhs,
    const std::vector<std::pair<int,double>>&             tightenedLoBnds,
    const std::vector<std::pair<int,double>>&             tightenedUpBnds,
    int                                                   n,
    int*                                                  types,
    double                                                eps = 1e-9,
    double                                                infinity = 1e15
){
    // Vorbedingungen
    if (sTerms.empty() || intPart.empty()) return false;
    for (const auto& t : intPart) {
        const int var = (int)t.first;
        if (var < n && fabs(improvedLbOf(var, colLowerBound, tightenedLoBnds, xlp)) > eps)
            return false; // Int-LB muss 0 sein
    }

    // ---------- Set isComplemented und complement of that set, notComplemented aufbauen ----------
    std::vector<bool> isComplemented(intPart.size(), false);
    std::vector<std::pair<unsigned int,double>> notComplemented;  // (posInIntPart, |x - U/2|)
    notComplemented.reserve(intPart.size());
    double rhsUse = rhsAfterSubst;

    // --- isComplemented und notComplemented aufbauen + Numerator anpassen ---
    for (int ii = 0; ii < (int)intPart.size(); ++ii) {
      const unsigned var = intPart[ii].first;
      const double   aj  = intPart[ii].second;
      if (var >= (unsigned)n) continue;
      
      const double Uj = improvedUbOf((int)var, colUpperBound, tightenedUpBnds, xlp);
      if (!(Uj < infinity)) continue;           // kein UB → bleibt in intPart
      
      const double xj = xlp[var];
      if (xj >= 0.5 * Uj) {
        isComplemented[ii] = true;            // var ∈ isComplemented
        rhsUse -= aj * Uj;                    // <<--- WICHTIG: Numerator anpassen
      } else if (xj > eps && xj < Uj - eps) {
        notComplemented.emplace_back(ii, fabs(xj - 0.5 * Uj));
      }
    }

    if (notComplemented.size() > 0) {
      std::sort(notComplemented.begin(), notComplemented.end(), [](std::pair<int, double> e1, std::pair<int, double> e2) {
	  return e1.second < e2.second;});
    }

    // ---------- Beste Basis-Ungleichung über alle Pivots ----------
    auto try_build_and_score = [&](double delta,
                               double rhsUse,               
                               const std::vector<bool>& complMask,
                               std::vector<std::pair<unsigned,double>>& outLHS,
                               double& outRHS, double& outS, double& outViol, double& outEffi)
			       {
				 std::vector<std::pair<unsigned,double>> lhs_cMIR;
				 double rhs_cMIR = 0.0, sCoef = 0.0;
				 buildcMirInequality(delta, rhsUse, intPart, colUpperBound, complMask,
						     lhs_cMIR, rhs_cMIR, sCoef, n);
				 double viol = 0.0, effi = 0.0;
				 computeCmirViolationAndEfficacy(lhs_cMIR, rhs_cMIR, sCoef, sStar, xlp, n, viol, effi, eps);
				 outLHS.swap(lhs_cMIR);
				 outRHS = rhs_cMIR; outS = sCoef; outViol = viol; outEffi = effi;
				 return (viol > eps);
			       };
    // --- Pivot-Suche (Basis) mit rhsUse ---
    bool found = false;
    double bestDelta = 0.0, bestEffi = 0.0, bestViol = 0.0, bestS = 0.0, bestRHS = 0.0;
    std::vector<std::pair<unsigned,double>> bestLHS;
    std::vector<bool> bestCompl = isComplemented;

    for (int ii=0; ii< intPart.size();ii++) {
      std::pair<int,double> pivot = intPart[ii];
      const unsigned pivVar = pivot.first;
      const double   delta  = pivot.second;
      if (pivVar >= (unsigned)n) continue;
      if (types && types[pivVar] != 0) continue;
      if (delta <= eps) continue;
      const double U_piv = improvedUbOf((int)pivVar, colUpperBound, tightenedUpBnds, xlp);
      if (!(U_piv < infinity)) continue;

      std::vector<std::pair<unsigned,double>> candLHS;
      double candRHS=0.0, candS=0.0, viol=0.0, effi=0.0;
      if (!try_build_and_score(delta, rhsUse, isComplemented, candLHS, candRHS, candS, viol, effi)) continue;

      if (!found || effi > bestEffi) {
        found = true; bestEffi = effi; bestViol = viol; bestDelta = delta;
        bestLHS = candLHS; bestRHS = candRHS; bestS = candS; bestCompl = isComplemented;
      }
      if (bestViol > eps && ii > 1000) break;
    }
    if (!found) return false;

    // --- Delta-Skalierungen 1/8... 1 (weiter mit rhsUse) ---
    for (double m=0.125; m < 1.01; m = 2*m) {
      const double deltaScaled = bestDelta * m; 
      if (deltaScaled <= eps) break;

      std::vector<std::pair<unsigned,double>> candLHS;
      double candRHS=0.0, candS=0.0, viol=0.0, effi=0.0;
      if (!try_build_and_score(deltaScaled, rhsUse, bestCompl, candLHS, candRHS, candS, viol, effi)) continue;

      if (effi > bestEffi) {
	bestEffi = effi; bestViol = viol; bestDelta = deltaScaled;
	bestLHS = candLHS; bestRHS = candRHS; bestS = candS;
      }
    }

    // --- notComplemented „move-to-isComplemented“ mit localRhsUse ---
    std::vector<bool> currCompl = bestCompl;
    double localRhsUse = rhsUse;
    for (int ii=0; ii < notComplemented.size();ii++) {
      std::pair<unsigned int, double> q = notComplemented[ii];
      const int pos = (int)q.first;
      if (pos < 0 || pos >= (int)intPart.size() || currCompl[pos]) continue;

      const unsigned var = intPart[pos].first;
      const double   aj  = intPart[pos].second;
      const double   Uj  = improvedUbOf((int)var, colUpperBound, tightenedUpBnds, xlp);
      if (!(Uj < infinity)) continue;

      const double rhsCand = localRhsUse - aj * Uj;     // lokale Anpassung
      currCompl[pos] = true;

      std::vector<std::pair<unsigned,double>> candLHS;
      double candRHS=0.0, candS=0.0, viol=0.0, effi=0.0;
      if (try_build_and_score(bestDelta, rhsCand, currCompl, candLHS, candRHS, candS, viol, effi) && effi > bestEffi) {
	bestEffi = effi; bestViol = viol; bestLHS = candLHS; bestRHS = candRHS; bestS = candS;
	localRhsUse = rhsCand;                         // behalten
      } else {
	currCompl[pos] = false;                       // revert
      }
      if (bestViol > eps && ii > 1000) break;
    }

    // --- Outputs setzen ---
    cMirCutLhs   = std::move(bestLHS);
    cMirCutRhs   = bestRHS;
    sCoefBestCut = bestS;
    return true;
}

  static bool cMirSeparationEQ( const std::vector<std::vector<data::IndexedElement>>& allRows, const std::vector<std::pair<unsigned int, double>>& rowAggregated, const int* listRowsAggregated, const std::vector<data::QpRhs> & rhsVec, const double* xlp, const data::QpNum* colUpperBound, const data::QpNum* colLowerBound, double& rowAggregatedRhs, std::pair< std::vector<std::pair<unsigned int, double>>,double > & cMirCut, const std::vector<std::pair<int, double>>& tightenedLoBnds, const std::vector<std::pair<int, double>>& tightenedUpBnds , int n, int *solu, double *xlpxtra, int *types, int *eass, bool isSave, int maxAggr)
{
  double sStar=0.0;
  int numSlacks=maxAggr;
  std::vector<std::pair<unsigned int, double>> intPartOfInputRow;
  double rhsAfterSubst = rowAggregatedRhs;
  std::vector<std::pair<unsigned int, double>> sTerms;
  std::vector<std::pair<unsigned int, double>> bestCut;
  double rhsBestCut = 0.0;
    // call bound substitution heuristic  
    if (!substituteBoundsAndSplitIntS_EQ( rowAggregated, xlp, xlpxtra, colUpperBound, colLowerBound, intPartOfInputRow, rhsAfterSubst, sStar, sTerms, n, types, tightenedLoBnds, tightenedUpBnds, numSlacks ))
      return false;

    std::vector<double> rowAggregatedDense(n+numSlacks,0.0);
    for (int ii=0;ii<rowAggregated.size();ii++) {
      rowAggregatedDense[rowAggregated[ii].first] = rowAggregated[ii].second;
    }
    std::vector<int> bestCutIndicesDense(n+numSlacks,-1);
    double sCoefBestCut = 0.0;
    
    if (!cMirSeparationEQ_core(xlp, sStar, colUpperBound, colLowerBound, intPartOfInputRow, rhsAfterSubst, sTerms, sCoefBestCut, /*lhsout, rhsout*/bestCut,rhsBestCut, tightenedLoBnds, tightenedUpBnds, n, types, 1e-9))
      return false;

    // write the best cut found with the model variables
    for (int ii=0;ii<n;ii++)
      bestCutIndicesDense[ii]=-1;
    for (int ii = 0; ii < bestCut.size();ii++) {
      bestCutIndicesDense[bestCut[ii].first] = ii;
    }
    return expandSTermThroughSlacks(
         sTerms,
         sCoefBestCut,
         n, numSlacks,
         allRows, rhsVec,
	 colUpperBound,
	 colLowerBound,
	 tightenedLoBnds,
	 tightenedUpBnds, listRowsAggregated,
         /*in-row*/  rowAggregated, rowAggregatedDense,
         /*LHS/RHS*/ bestCut, bestCutIndicesDense, rhsBestCut, xlp, cMirCut,
	 1e-4
    );
  }


static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >* CMIRizeEQ( extSol::QpExternSolver &extSolver, const std::vector<std::pair<unsigned,double>>& invec, double inrhs, int *types, int8_t *assigns, int decLev, unsigned int initime, int * solu, int* fixs, int *blcks, int *eass, int orgN, char senseChar,
    const std::vector<std::pair<int,double>>*          tightenedLoBndsOpt = nullptr,
    const std::vector<std::pair<int,double>>*          tightenedUpBndsOpt = nullptr
		      ){
    unsigned int m = extSolver.getRowCount();
    unsigned int n = extSolver.getVariableCount();
    
    std::vector<data::QpRhs> rhsVec( 1 );
    std::vector< std::vector< data::IndexedElement >  > allrows(1);
    allrows.resize(1);
    allrows[0].resize(invec.size());
    rhsVec.resize(1);
    rhsVec[0].setValue(senseChar=='>'?-inrhs:inrhs);
    rhsVec[0].setRatioSign(data::QpRhs::smallerThanOrEqual);
    for (int i=0; i < invec.size();i++) {
      data::IndexedElement e;
      e.index = invec[i].first;
      if (senseChar=='>') e.value = -invec[i].second;
      else                e.value = invec[i].second;
      allrows[0][i] = e;
    }

    return getMIRsmartCore( extSolver, allrows, rhsVec, types, assigns, decLev, initime, solu, fixs, blcks, eass, orgN, 1, true);    
}
// end of MIRsmart
  
static void aggregate2RowsInFirst(std::vector<std::pair<unsigned int, double>> &aggRow /*must be sorted by ascending index*/, double &rhs_aggRow, const std::vector<std::pair<unsigned int, double>> &b, double rhs_b, double w, double tol, int pivot);  
  static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >* getCoverCuts( extSol::QpExternSolver &extSolver, int *types, int8_t *assigns, int decLev, unsigned int, int *, int *, int*, int*, int );
	static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >* getLPCuts(extSol::QpExternSolver &extSolver, std::vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev, unsigned int initime, std::vector<int> &listOfCutsVars);

	static void setOrgM(int m) { orgM = m; }
	static int getOrgM() { return orgM; }
private:
	static int orgM;
	static Param_c param;
	static const int info_level = 1;
	CutAdder();
	static inline bool isZero(double x, double epsZero = 1e-20) {
		return (fabs(x) <= epsZero);
	}
        static bool getAllRows(extSol::QpExternSolver &extSolver, std::vector< std::vector< data::IndexedElement >  > &allrows2);

	static bool createCover( std::vector<double>& xlpopt, const std::vector<std::pair<double, unsigned int> >& rowsparse, double rhs, const double eps, std::vector<std::pair<unsigned int, double> >& reslhs, double& resrhs, int *types,  int8_t *assigns, int decLev, data::QpNum*, data::QpNum*, int *solu, int*fixs, int *blcks, int orgN );
};

#endif /* CUTADDER_H_ */
