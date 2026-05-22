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
#include <unordered_set>

#include <assert.h>

#include "GFkSolve.h"

#define DEBUG_OUT_EQ 0

class Param_c {
public:

};

#define fabs(x) ((x) >= 0 ? (x) : -(x))
class CutAdder {
public:
  static std::vector<std::pair<int, double>> globalConditionalLowerBounds;
  static std::vector<std::pair<int, double>> globalConditionalUpperBounds;
  static std::vector<data::QpNum> rootLBs;
  static std::vector<data::QpNum> rootUBs;  
  static std::vector<uint8_t> boundTypes;
  enum BoundType : uint8_t {
			    yNone       = 0,
			    ySimpleLb   = 1,  // z = x - lb
			    ySimpleUb   = 2,  // z = ub - x
			    yVariableLb = 3,  // z = x - (coef*y)
			    yVariableUb = 4   // z = (coef*y) - x
  };
  enum class AggregationMode {
			      RowWiseCMIR,
			      HeuristicAggregation,
			      ModK_GFk,
			      FromCMIRize
  };

  struct CoverGenState {
    std::vector<int> cover;
    double coverweight = 0.0;
    double lambda = 0.0;

    std::vector<double> upper;
    std::vector<double> solval;
    std::vector<uint8_t> complementation;
    std::vector<uint8_t> isintegral;

    std::vector<double> vals;   // coeffs of current row (mutable)
    std::vector<int>    inds;   // indices of current row
    double rhs = 0.0;
    bool integralSupport = false;
    bool integralCoefficients = false;
    int rowlen = 0;
    double initialScale = 1.0;

    std::vector<int> integerinds;
    std::vector<double> deltas;

    double feastol = 1e-9;
    double epsilon = 1e-12;
  };

  // Returns a random float 0 <= x < 1. Seed must never be 0.
  static inline double drand(double& seed) {
    seed *= 1389796;
    int q = (int)(seed / 2147483647);
    seed -= (double)q * 2147483647;
    return seed / 2147483647; }
  
  // Returns a random integer 0 <= x < size. Seed must never be 0.
  static inline uint64_t irand(double& seed, uint64_t size) {
    return (uint64_t)(drand(seed) * size); }

  static bool getAllRows_snapshot(extSol::QpExternSolver &extSolver, std::vector< std::vector< data::IndexedElement >  > &allrows2, std::vector<data::QpRhs>&rhsVec) {
    unsigned int n = extSolver.getVariableCount();
    int m = extSolver.getRowRhs_snapshot()->size();
    allrows2.clear();
    rhsVec.clear();

    for( unsigned int i = 0; i < m; ++i ){
      std::vector<data::IndexedElement> &rowtmp = *extSolver.getRowLhs_snapshot(i);
      data::QpRhs rhs = (*extSolver.getRowRhs_snapshot())[i];
      if (extSolver.getLazyStatus(i) == true || rowtmp.size() == 0) continue;
      if (0&&extSolver.getStatus(i)) {
	continue;
      }
      if (rowtmp.size() == 0) continue;
      allrows2.push_back(rowtmp);
      rhsVec.push_back(rhs);
    }
    return true;
  }
  struct CoverPrep {
    std::vector<int> inds;                 // indices (0..n-1)
    std::vector<double> vals;              // coeffs >= 0
    std::vector<double> upper;             // U_i = ub-lb (finite)
    std::vector<double> solval;            // shifted/complemented LP value in [0,U]
    std::vector<uint8_t> isintegral;       // 1 if integer (types[i]==0)
    std::vector<uint8_t> complementation;  // 1 if complemented
    double rhs = 0.0;                      // transformed rhs
  };

  // returns false if row not suitable for cover-separation
  static bool preprocessForCoverFromRow(
					const std::vector<std::pair<unsigned int,double>>& row, // sparse
					double rhs,
					const data::QpRhs::RatioSign sense,
					const data::QpNum* lb,
					const data::QpNum* ub,
					const double* xlp,
					int n,                 // number of original vars
					int orgN,              // your "bigX" threshold (>=orgN skip)
					const int* types,      // 0==integer
					const int8_t* assigns, // 2==free
					const int* fixs,       // 2==free
					CoverPrep& out,
					bool &conIn,
					double inf = 1e15,
					double eps = 1e-12)
  {
    out.inds.clear(); out.vals.clear(); out.upper.clear();
    out.solval.clear(); out.isintegral.clear(); out.complementation.clear();
    out.rhs = rhs;
    conIn=false;

    // 1) bring to <= by possibly multiplying with -1
    double mult = 1.0;
    if (sense == data::QpRhs::greaterThanOrEqual) mult = -1.0;
    // equality: leave as-is (cover on == is tricky; for now we reject)
    if (sense == data::QpRhs::equal) return false;

    // 2) build list of candidate items (only original vars < n)
    for (const auto& p : row) {
      const int idx = (int)p.first;
      //if (eass[idx]==1) continue;
      if (idx >= n) continue;            // ignore slacks n+a
      if (idx < n && idx >= orgN) return false;     // bigX present -> skip (keep it strict)

      double a = mult * p.second;
      //if (fabs(a) <= eps) continue;

      // only integer vars for first version
      if (types[idx] != 0) {
	// allow continuous only if fixed
	const bool fixedByAssign = (assigns[idx] != 2);
	const bool fixedByFixs   = (fixs[idx] == 0 || fixs[idx] == 1);
	const bool fixedByBounds = (fabs(ub[idx].asDouble() - lb[idx].asDouble()) < 1e-7);
	if (!(fixedByAssign || fixedByFixs || fixedByBounds)) conIn=true;//return false;
	// if it is fixed, it can be substituted away by shifting rhs:
	// treat it like constant -> incorporate directly and skip as item
	/*
	const double x = xlp[idx];
	out.rhs -= a * x;
	continue;
	*/
      }

      
      const double lbi = lb[idx].asDouble();
      const double ubi = ub[idx].asDouble();
      /*
      double lbi=0.0;
      double ubi=0.0;
      if (types[idx] != 0) {
	lbi = lb[idx].asDouble();
	ubi = ub[idx].asDouble(); 
      } else {
	lbi = 0.0; ubi = 1.0;
	if (CutAdder::rootLBs.size() == n) {
	  lbi = CutAdder::rootLBs[idx].asDouble();
	  ubi = CutAdder::rootUBs[idx].asDouble();
	}
      }
      */
      if (ubi >= inf/10) return false; // need finite U for cover
      const double Ui = ubi - lbi;
      if (0&&Ui <= eps) {
	// fixed integer -> substitute
	const double x = xlp[idx];
	out.rhs -= a * x;
	continue;
      }

      // shift x = lb + x', so rhs -= a*lb
      out.rhs -= a * lbi;
      double xsh = xlp[idx] - lbi;

      // make coefficient nonnegative by complementing if needed
      uint8_t comp = 0;
      if (a < 0.0) {
	comp = 1;
	// a*(lb + x') with a<0 -> after complement:
	// a*x' = a*(Ui - x'') = a*Ui + (-a)*x''
	out.rhs -= a * Ui;
	a = -a;
	xsh = Ui - xsh;
      }

      // clamp xsh (numerical)
      if (xsh < 0.0) xsh = 0.0;
      if (xsh > Ui)  xsh = Ui;

      out.inds.push_back(idx);
      out.vals.push_back(a);
      out.upper.push_back(Ui);
      out.solval.push_back(xsh);
      out.isintegral.push_back(types[idx] == 0 ? 1 : 0);
      //out.isintegral.push_back(1);
      out.complementation.push_back(comp);
    }

    // rhs must be meaningful
    if (out.inds.empty()) return false;

    // also ensure rhs is not negative huge etc.
    if (!std::isfinite(out.rhs)) return false;

    return true;
  }
  // Returns true if a cover was found.
  // Outputs:
  //  - cover: positions (indices into prep arrays, not variable indices!)
  //  - coverweight: sum a_i * upper_i over cover
  //  - lambda: coverweight - rhs
  static bool determineCoverHiGHSlike(
				      const CoverPrep& prep,
				      std::vector<int>& cover,
				      double& coverweight,
				      double& lambda,
				      const data::QpNum* colLowerBound,
				      const data::QpNum* colUpperBound,
				      bool lpSol = true,
				      double feastol = 1e-9,
				      double eps = 1e-12)
  {
    cover.clear();
    coverweight = 0.0;
    lambda = 0.0;
    const double cover_eps=0.0;//1e-9;

    const int len = (int)prep.inds.size();
    if (len == 0) return false;
    int numInt = 0;
    for (int i=0;i<len;++i) if (prep.isintegral[i] && prep.upper[i] > eps) ++numInt;
    if (numInt < 2) return false;

    // We need something to cover: rhs must be < sum a_i * upper_i to have any chance.
    double maxWeight = 0.0;
    for (int i = 0; i < len; ++i) {
      if (!prep.isintegral[i]) continue;
      if (prep.upper[i] <= eps) continue;
      maxWeight += prep.vals[i] * prep.upper[i];
    }
    if (maxWeight <= prep.rhs + feastol + fabs(fabs(prep.rhs)+fabs(maxWeight))*feastol) return false;

    // Candidate score: HiGHS uses LP solution information; we mimic a robust variant.
    // Prefer variables that are "active" in the LP solution and have large a_i*upper_i.
    std::vector<int> order;
    order.reserve(len);
    for (int i = 0; i < len; ++i) {
      if (!prep.isintegral[i]) continue;
      if (prep.upper[i] <= eps) continue;
      if (colUpperBound[prep.inds[i]].asDouble()-colLowerBound[prep.inds[i]].asDouble() < 1e-8 && prep.upper[i] >= 1.0-1e-8)
	continue;
      order.push_back(i);
    }

    auto score = [&](int i) -> double {
		   const double ai = prep.vals[i];
		   const double Ui = prep.upper[i];
		   assert(fabs(1.0-Ui) < 0.01);
		   const double xi = prep.solval[i];
		   // "lpSol" tries to pick vars with larger xi; but also weight with ai*Ui.
		   // This mirrors the HiGHS idea: likely cover items should be strong and present in LP.
		   return (lpSol ? (xi /*+ 1e-6*/) : 1.0) * (ai * Ui);
		 };

    std::sort(order.begin(), order.end(),
	      [&](int a, int b) { return score(a) > score(b); });

    // Greedy cover: add items until coverweight > rhs
    long double coverweight2 = 0.0;
    for (int i = 0; i < len; ++i) {
      if (!prep.isintegral[i]) continue;
      if (prep.upper[i] <= eps) continue;
      if (colUpperBound[prep.inds[i]].asDouble()-colLowerBound[prep.inds[i]].asDouble() >= 1e-8 ||
	  prep.upper[i] < 1.0-1e-8)
	continue;
      const double ai = prep.vals[i];
      const double Ui = prep.upper[i];
      const double w  = ai * Ui;
      if (w <= eps) continue;
      cover.push_back(i);
      coverweight2 += w;
      if (coverweight2 > prep.rhs /*+ feastol*/+ fabs(fabs(coverweight2)+fabs(prep.rhs))*cover_eps + cover_eps) break;
    }
    
    for (int pos : order) {
      const double ai = prep.vals[pos];
      const double Ui = prep.upper[pos];
      const double w  = ai * Ui;
      if (w <= eps) continue;
      cover.push_back(pos);
      coverweight2 += w;
      if (coverweight2 > prep.rhs /*+ feastol*/+ fabs(fabs(coverweight2)+fabs(prep.rhs))*cover_eps + cover_eps) break;
    }

    if (coverweight2 <= prep.rhs /*+ feastol*/+ fabs(fabs(coverweight2)+fabs(prep.rhs))*cover_eps + cover_eps) {
      cover.clear();
      coverweight = 0.0;
      return false;
    }

    // Minimalize the cover: try removing items while keeping it a cover.
    // Classic minimal cover: remove items with smallest weight first.
    std::sort(cover.begin(), cover.end(), [&](int a, int b) {
					    return prep.vals[a] * prep.upper[a] < prep.vals[b] * prep.upper[b];
					  });

    for (int it = 0; it < (int)cover.size(); /*increment inside*/) {
      const int pos = cover[it];
      const double w = prep.vals[pos] * prep.upper[pos];
      if (coverweight2 - w > prep.rhs + cover_eps/*feastol*/ + fabs(fabs(coverweight2)+fabs(prep.rhs))*cover_eps+cover_eps) {
	coverweight2 -= w;
	cover.erase(cover.begin() + it);
	continue;
      }
      ++it;
    }

    coverweight = (double)coverweight2;
    lambda = coverweight - prep.rhs;
    if (lambda <= feastol) return false;

    return !cover.empty();
  }
  static bool separateLiftedKnapsackCoverHiGHSlike(
						   const CoverPrep& prep,
						   std::vector<int> cover,               // positions in prep arrays (copy, will be sorted)
						   const double lambda,
						   std::vector<double>& alpha,           // output: coeffs on prep positions
						   double& beta,                         // output: rhs on prep vars
						   const data::QpNum* lb,
						   const data::QpNum* ub,
						   int *eass,
						   int *blcks,
						   int orgN,
						   const double feastol = 1e-9,
						   const double eps = 1e-12)
  {
    const int rowlen = (int)prep.inds.size();
    const int coversize = (int)cover.size();
    if (coversize <= 1) return false;
    if (rowlen == 0) return false;
    if (!(lambda > feastol)) return false;
    //return false;
    // local copy of coefficients (HiGHS modifies vals in place)
    std::vector<double> vals = prep.vals;

    // sort cover by descending vals
    std::sort(cover.begin(), cover.end(),
	      [&](int a, int b) { return vals[a] > vals[b]; });

    // compute abar (HiGHS: balancing using lambda)
    long double abartmp = vals[cover[0]];
    double sigma = lambda;
    for (int i = 1; i < coversize; ++i) {
      const long double delta = abartmp - vals[cover[i]];
      const long double kdelta = (long double)i * delta;
      if (kdelta < sigma) {
	abartmp = vals[cover[i]];
	sigma -= kdelta;
      } else {
	abartmp -= sigma * (1.0 / (double)i);
	sigma = 0.0;
	break;
      }
    }
    if (sigma > 0.0) abartmp = prep.rhs / (long double)coversize;
    const long double abar = abartmp;

    // build S[] and coverflag[] (HiGHS)
    std::vector<double> S(coversize, 0.0);
    std::vector<int8_t> coverflag(rowlen, 0);  // default 0 for non-cover
    long double sum = 0.0;
    int cplussize = 0;

    /*
    std::cerr << "ORIGINAL COVER: ";
    for (int i = 0; i < coversize; ++i) {
      std::cerr << " " << prep.inds[cover[i]]; 
    }
    std::cerr << std::endl;
    */
    for (int i = 0; i < coversize; ++i) {
      sum += /*std::*/fmin(abar, vals[cover[i]]);
      S[i] = (double)sum;

      if (vals[cover[i]] > abar + feastol) {
	++cplussize;
	coverflag[cover[i]] = 1;   // C+
      } else {
	coverflag[cover[i]] = -1;  // C0
      }
    }

    bool halfintegral = false;

    // lifting function g(z) (HiGHS)
    auto g = [&](double z) -> double {
	       const long double hfrac = z / abar;
	       double coef = 0.0;
	       const double cover_lift_eps = 0.0;

	       int h = (int)std::floor(hfrac + 0.5);
	       if (abar>0.1 && h != 0 &&
		   std::abs(hfrac - (double)h) * /*std::*/fmax(1.0, abar) <= eps &&
		   h <= cplussize - 1) {
		 halfintegral = true;
		 coef = 0.5;
	       }

	       h = std::max(h - 1, 0);
	       for (; h < coversize; ++h) {
		 if (z <= S[h] + feastol) break;
	       }

	       return (coef + (double)h)-cover_lift_eps*(double)coversize*fabs(coef + (double)h);
	     };

    // HiGHS: rhs = coversize - 1, vals transformed to lifted coefficients
    beta = (double)coversize - 1.0;
    /*
    std::cerr << "coversize="<<coversize<<" cplussize="<<cplussize
	      <<" abar="<<abar<<" lambda="<<lambda<<" beta="<<beta<<"\n";
    for (int h=0; h<std::min(coversize, 10); ++h)
      std::cerr << "S["<<h<<"]="<<S[h]<<"\n";
    */
    
    int maxBlock=-1;
    int minBlock=blcks[prep.inds[0]];
    for (int i = 0; i < rowlen; ++i) {
      assert(ub[prep.inds[i]].asDouble()>=lb[prep.inds[i]].asDouble());
      //if (ub[prep.inds[i]].asDouble()-lb[prep.inds[i]].asDouble() < 1e-8)
      //continue;
      if (blcks[prep.inds[i]] > maxBlock) maxBlock = blcks[prep.inds[i]];
      else if (blcks[prep.inds[i]] < minBlock) minBlock = blcks[prep.inds[i]];
    }
    for (int i = 0; i < rowlen; ++i) {
      if (vals[i] == 0.0) continue;
      //assert(ub[prep.inds[i]].asDouble()-lb[prep.inds[i]].asDouble() > 1e-8);
      if (coverflag[i] == 0 && (/*blcks[prep.inds[i]] <*/minBlock != maxBlock || prep.inds[i] >= orgN || eass[prep.inds[i]] == 1 || /*xlp[prep.inds[i]]<=1e-4 || xlp[prep.inds[i]]>=1.0-1e-4 ||*/ ub[prep.inds[i]].asDouble()-lb[prep.inds[i]].asDouble() < 1e-8 )) {
	vals[i] = 0.0;
      } else if (coverflag[i] == -1) {
	vals[i] = 1.0;
      } else {
	int z=vals[i];
	//std::cerr << "z="<<z<<" g(z)="<<g(z)<<" (S1="<<S[1]<<")\n";
	vals[i] = fmax(g(vals[i]),0.0);
      }
    }

    if (halfintegral) {
      beta *= 2.0;
      for (int i = 0; i < rowlen; ++i) vals[i] *= 2.0;
    }
    beta = beta;

    alpha.swap(vals);
    return true;
  }
  // Port von HiGHS: separateLiftedMixedBinaryCover()
  // Erwartet: prep ist already knapsack-normalisiert (vals>=0, rhs entsprechend),
  // cover sind Positionen in prep-Arrays, lambda = coverweight - rhs (>0)
  static bool separateLiftedMixedBinaryCoverHiGHSlike(
						      const CoverPrep& prep,
						      std::vector<int> cover,               // copy; will be sorted
						      const double lambda,
						      std::vector<double>& alpha,           // output coeffs on prep positions
						      double& beta,                         // output rhs on prep vars
						      int *eass,
						      int *blcks,
						      int orgN,
						      const double epsilon = 1e-12)
  {
    const int rowlen = (int)prep.inds.size();
    const int coversize = (int)cover.size();
    if (coversize == 0) return false;
    if (!(lambda > 0.0)) return false;
    //return false;
    
    // local copy of coefficients (HiGHS modifies vals in place)
    std::vector<double> vals = prep.vals;

    // coverflag
    std::vector<uint8_t> coverflag(rowlen, 0);
    for (int i = 0; i < coversize; ++i) coverflag[cover[i]] = 1;

    // sort cover by descending vals
    std::sort(cover.begin(), cover.end(),
	      [&](int a, int b) { return vals[a] > vals[b]; });

    // Build S over the subset of cover items with vals > lambda (within epsilon)
    std::vector<double> S(coversize, 0.0);
    int p = coversize;
    double sum = 0.0;

    for (int i = 0; i < coversize; ++i) {
      if (vals[cover[i]] - lambda <= epsilon) {
	p = i;
	break;
      }
      sum += vals[cover[i]];
      S[i] = sum;
    }
    if (p == 0) return false;

    // Lifting function phi(a)
    auto phi = [&](double a) -> double {
		 for (int i = 0; i < p; ++i) {
		   if (a <= S[i] - lambda) return (double)((long double)i * lambda);
		   if (a <= S[i])         return (double)((long double)(i + 1) * lambda + (a - (long double)S[i]));
		 }
		 return (double)((long double)p * lambda + ((long double)a - S[p - 1]));
	       };

    // In HiGHS: rhs = -lambda; then add min(vals_i, lambda) for cover members
    double rhs = -lambda;
    int maxBlock = -1;
    int minBlock = blcks[prep.inds[0]];
    for (int i = 0; i < rowlen; ++i) 
      if (blcks[prep.inds[i]] > maxBlock) maxBlock = blcks[prep.inds[i]];
      else if (blcks[prep.inds[i]] < minBlock) minBlock = blcks[prep.inds[i]];
    // Only keep integral vars (binary in your world). Nonintegral with vals>=0 are dropped.
    for (int i = 0; i < rowlen; ++i) {
      if (!prep.isintegral[i]) {
	// In your prep: vals[i] >= 0 typically -> can set to 0 safely
	if (vals[i] < 0) return false; // would break validity assumptions
	vals[i] = 0.0;
	continue;
      }

      if (coverflag[i]) {
	vals[i] = std::min(vals[i], lambda);
	rhs += vals[i];
      } else if (maxBlock != minBlock) {
	vals[i] = 0.0;
      } else {
	vals[i] = phi(vals[i]);
      }
    }

    alpha.swap(vals);
    beta = rhs;
    //std::cerr << "MIXEDBIN!" << std::endl;
    return true;
  }
  static void combineDuplicates(
				std::vector<std::pair<unsigned,double>>& row,
				const double dropTol = 1e-12)
  {
    std::unordered_map<unsigned, double> acc;
    acc.reserve(row.size() * 2);

    for (const auto& [idx, coef] : row) {
      if (fabs(coef) <= dropTol) continue;
      acc[idx] += coef;
    }

    row.clear();
    row.reserve(acc.size());

    for (const auto& kv : acc) {
      const unsigned idx = kv.first;
      const double coef  = kv.second;
      if (fabs(coef) > dropTol) row.emplace_back(idx, coef);
    }

    // optional: sort by index for determinism (nice for debugging & hashing)
    std::sort(row.begin(), row.end(),
	      [](std::pair<unsigned,double>& a, std::pair<unsigned,double>& b){ return a.first < b.first; });
  }
  static bool untransformAndGenerateCutEQ(
					// --- Input: already aggregated row and rhs (your snippet uses these) ---
					const std::vector<std::pair<unsigned, double>>& rowAggregated,
					double rhsAggregatedValue,

					// --- Input: original row + rhs used for preprocess (you already pass these) ---
					const std::vector<std::pair<unsigned, double>>& rowToUse,
					const data::QpNum& rhsOfRowToUse,

					// --- LP solution (main + extra for "idx>=n" terms in the ax-check) ---
					const std::vector<double>& xlpopt,
					const std::vector<double>& xlpExtra,

					// --- bounds (need .asDouble() inside) ---
					const std::vector<data::QpNum>& colLowerBound,
					const std::vector<data::QpNum>& colUpperBound,

					// --- dimensions and var meta ---
					unsigned n,
					unsigned orgN,
					int *types,
				        int8_t *assigns,
					int *fixs,
					int *eass,
					int *blcks,

					// --- output: your cut container (same structure you showed) ---
					std::vector<std::pair<std::vector<std::pair<unsigned, double>>, double>>& cuts,

					// --- knobs ---
					bool debug_output = false,
					double feastol = 1e-9,
					double eps = 1e-12,
					double violTol = 1e-7)
  {
    bool comparisonOk=false;
    bool prepOk=false;
    bool coverOk=false;

    CoverPrep prep;

    // 1) consistency check s = rhs - a*x >= 0 (on aggregated row)
    double ax = 0.0;
    for (auto [idx,coef] : rowAggregated) {
      if ((int)idx < (int)n) ax += coef * xlpopt[idx];
      else ax += coef * xlpExtra[(int)idx - (int)n];
    }
    const double s = rhsAggregatedValue - ax;
    if (s < -violTol) {
      comparisonOk=false;
    } else {
      comparisonOk=true;
    }

    if (!comparisonOk) return false;

    // 2) preprocess for cover
    double minSlack =  1e100;
    double maxSlack = -1e100;
    int cntSlack = 0;
    for (auto [idx,coef] : rowToUse) {
      if ((int)idx >= (int)n) { // slack/extra
	cntSlack++;
	minSlack = std::min(minSlack, coef);
	maxSlack = std::max(maxSlack, coef);
      }
    }
    //if (cntSlack) std::cerr << "[COVER] slacks="<<cntSlack<<" minCoef="<<minSlack<<" maxCoef="<<maxSlack<<"\n";
    //assert(minSlack >= 0);
    bool conIn=false;
    if (minSlack >= 0.0) 
      prepOk = preprocessForCoverFromRow(
				       rowToUse, rhsOfRowToUse.asDouble(), data::QpRhs::smallerThanOrEqual,
				       colLowerBound.data(), colUpperBound.data(),
				       xlpopt.data(), n, orgN, types, assigns, fixs, prep, conIn);

    if (!prepOk) return false;

    // 3) determine cover
    std::vector<int> cover;
    double coverweight = 0.0;
    double lambda = 0.0;

    const bool hasCover = determineCoverHiGHSlike(
						  prep, cover, coverweight, lambda, colLowerBound.data(), colUpperBound.data(),
						  /*lpSol=*/true, /*feastol=*/feastol, /*eps=*/eps);

    coverOk = hasCover;
    if (!coverOk) return false;

    // 4) decide which lifting: knapsack vs mixed-binary (HiGHS-like)
    bool hasContinuous = conIn;
    for (int pos = 0; pos < (int)prep.inds.size(); ++pos) {
      if (!prep.isintegral[pos] && fabs(prep.vals[pos]) > eps) {
	hasContinuous = true;
	break;
      }
    }

    std::vector<double> alpha;
    double beta = 0.0;
    bool liftOk = false;

    if (!hasContinuous) {
      liftOk = separateLiftedKnapsackCoverHiGHSlike(prep, cover, lambda, alpha, beta, colLowerBound.data(), colUpperBound.data(), eass, blcks, orgN, feastol, eps);
    } else {
      liftOk = separateLiftedMixedBinaryCoverHiGHSlike(prep, cover, lambda, alpha, beta, eass, blcks, orgN, eps);
    }

    if (!liftOk) return false;

    // 5) backtransform prep cut -> original cutSparse: sum c*x <= rhsCut
    std::vector<std::pair<unsigned, double>> cutSparse;
    cutSparse.reserve(prep.inds.size());
    double rhsCut = beta;

    for (int pos = 0; pos < (int)prep.inds.size(); ++pos) {
      const double a = alpha[pos];
      if (a == 0.0) continue;

      const unsigned col = (unsigned)prep.inds[pos];
      const double lb = colLowerBound[col].asDouble();
      const double ub = colUpperBound[col].asDouble();

      if (prep.complementation[pos] == 0) {
	cutSparse.push_back({col, +a});
	rhsCut += a * lb;
      } else {
	cutSparse.push_back({col, -a});
	rhsCut += a * ub;
      }
    }

    // After you built cutSparse (in the "mod-k z-space") and rhsCut:
    //
    // interpret each term as coefficient * z_j and rewrite z_j in terms of x/y.

    /*
    std::vector<std::pair<unsigned,double>> cutFinal;
    cutFinal.reserve(cutSparse.size() + 8); // maybe add some y-terms

    double rhsFinal = rhsCut;

    // Helper: add (idx,coef) accumulating duplicates (optional; simplest: push and combine later)
    auto pushTerm = [&](unsigned idx, double coef) {
		      if (fabs(coef) <= 1e-18) return;
		      cutFinal.emplace_back(idx, coef);
		    };

    for (auto [col, coefZ] : cutSparse) {
      // col is an x-index (0..n-1) but might represent z depending on boundTypes
      assert(col < n);
      const uint8_t bt = CutAdder::boundTypes[col];

      if (bt == CutAdder::yNone) {
	// z == x (no substitution)
	pushTerm(col, coefZ);
	continue;
      }

      const double lb = colLowerBound[col].asDouble();
      const double ub = colUpperBound[col].asDouble();

      if (bt == CutAdder::ySimpleLb) {
	// z = x - lb  =>  coefZ*z = coefZ*x - coefZ*lb
	pushTerm(col, coefZ);
	rhsFinal += coefZ * lb;  // move (-coefZ*lb) to rhs of <= :  coefZ*x <= rhs + coefZ*lb
	continue;
      }

      if (bt == CutAdder::ySimpleUb) {
	// z = ub - x  =>  coefZ*z = coefZ*ub - coefZ*x
	pushTerm(col, -coefZ);
	rhsFinal += coefZ * ub;
	continue;
      }

      if (bt == CutAdder::yVariableLb) {
	// z = x - (coefY * y)   with (yIndex, coefY) stored in globalConditionalLowerBounds[col]
	const int yIndex = CutAdder::globalConditionalLowerBounds[col].first;
	const double coefY = CutAdder::globalConditionalLowerBounds[col].second;
	// coefZ*z = coefZ*x - coefZ*coefY*y
	pushTerm(col, coefZ);
	if (yIndex >= 0) pushTerm((unsigned)yIndex, -coefZ * coefY);
	else {
	  // if no yIndex, treat as simple lb? (depends on your conventions)
	  // safest: reject
	  return false;
	}
	continue;
      }

      if (bt == CutAdder::yVariableUb) {
	// z = (coefY * y) - x   with (yIndex, coefY) stored in globalConditionalUpperBounds[col]
	const int yIndex = CutAdder::globalConditionalUpperBounds[col].first;
	const double coefY = CutAdder::globalConditionalUpperBounds[col].second;
	// coefZ*z = coefZ*coefY*y - coefZ*x
	pushTerm(col, -coefZ);
	if (yIndex >= 0) pushTerm((unsigned)yIndex, +coefZ * coefY);
	else return false;
	continue;
      }

      // unknown type
      return false;
    }
    combineDuplicates(cutFinal, 1e-12);

    // Optionally: combine duplicates in cutFinal here (recommended)
    cutSparse.swap(cutFinal);
    rhsCut = rhsFinal;
    */
    // 6) violation check (<= cut)
    double act = 0.0;
    double largest = 0.0;
    double smallest = 1e100;
    for (auto [j,c] : cutSparse) {
      largest=fmax(largest,fabs(c));
      smallest=fmin(smallest,fabs(c));
      act += c * xlpopt[j];
    }
    const double viol = act - rhsCut;

    if (viol <= violTol) return false;
    if (largest>1.0e8*smallest||largest>1.0e7||smallest<1.0e-5) return false;
			

    if (debug_output) {
      std::cerr << "transformAndGenerateCutEQ: added cover cut, viol=" << viol
		<< " hasContinuous=" << hasContinuous
		<< " |cut|=" << cutSparse.size()
		<< " lambda=" << lambda << "\n";
    }

    // 7) push into your cuts container (you store >= form)
    cuts.emplace_back();
    for (int ii = 0; ii < (int)cutSparse.size(); ++ii) {
      cuts.back().first.emplace_back(cutSparse[ii].first, -cutSparse[ii].second);
    }
    cuts.back().second = -rhsCut;

    return true;
  }

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
    
    static double coeffNorm2Approx(const std::vector<std::pair<unsigned int, double> >& lhsO) {
        double s = 0.0;
	int samplesize=lhsO.size();
	std::vector<std::pair<unsigned int, double> > lhs;
	if (samplesize > 100) samplesize = 100;
	for (int i=0;i < samplesize/*loglen1*/;i++) {
	  int pos;
	  if (samplesize == lhsO.size()) pos = i; 
	  else pos = irand(CutAdder::random_seed,lhsO.size());
	  lhs.push_back(lhsO[pos]);
	}
        for (const auto& kv : lhs) {
            const long double a = kv.second;
            s += a * a;
        }
	if (samplesize<lhsO.size())
	  s = s *(double)lhsO.size() / (double)samplesize;
        return std::sqrt(s);
    }
    
    // Aktivität a^T x
    static inline double activityApprox(const std::vector<std::pair<unsigned int, double> >& lhsO, const double* x, int n) {
      double s = 0.0;
      int samplesize=lhsO.size();
      std::vector<std::pair<unsigned int, double> > lhs;
      if (samplesize > 100) samplesize = 100;
      for (int i=0;i < samplesize/*loglen1*/;i++) {
	int pos;
	if (samplesize == lhsO.size()) pos = i; 
	else pos = irand(CutAdder::random_seed,lhsO.size());
	lhs.push_back(lhsO[pos]);
      }
      for (const auto& kv : lhs) {
	const unsigned j = kv.first;
	const long double a = kv.second;
	const long double xj = (x && j < (unsigned)n) ? x[j] : 0.0L;
	s += a * xj;
      }
      if (samplesize<lhsO.size())
	s = s * (double)lhsO.size() / (double)samplesize;
      return s;
    }
    // Aktivität a^T x
  static inline double activityApprox(const std::vector<data::IndexedElement >& lhsO, const double* x, int n) {
      double s = 0.0;
      int samplesize=lhsO.size();
      std::vector<data::IndexedElement > lhs;
      if (samplesize > 100) samplesize = 100;
      for (int i=0;i < samplesize/*loglen1*/;i++) {
	int pos;
	if (samplesize == lhsO.size()) pos = i; 
	else pos = irand(CutAdder::random_seed,lhsO.size());
	lhs.push_back(lhsO[pos]);
      }
      for (const auto& kv : lhs) {
	const unsigned j = kv.index;
	const long double a = kv.value.asDouble();
	const long double xj = (x && j < (unsigned)n) ? x[j] : 0.0L;
	s += a * xj;
      }
      if (samplesize<lhsO.size())
	s = s * (double)lhsO.size() / (double)samplesize;
      return s;
    }

    // Für <=-Constraints: slack = rhs(row) - a^T x
    // Rückgabe robust (clamp) in [0, rhs], falls rhs endlich ist.
    static inline double slackValueAtXlp(
				       int row,
				       const std::vector<std::vector<data::IndexedElement>>& allrows,
				       const std::vector<double>& xlp,
				       const std::vector<data::QpRhs>& rhsVec)
    {
      // Slackbounds für <=: [0, rhs]
      const double ub = rhsVec[row].getValue().asDouble();   // slackUpper
      const double lb = 0.0;                                // slackLower
      
      // Wenn ub "inf" sein kann, hier ggf. anders behandeln.
      // In deinem Setup sagst du: klar 0..INF, aber rhs ist i.d.R. endlich.
      // Wir clampen nur, wenn ub sinnvoll ist.
      const double act = activityApprox(allrows[row], xlp.data(), xlp.size());
      double s = ub - act;
      
      // Robustheit: clamp
      if (s < lb) s = lb;
      if (s > ub) s = ub;
      
      return s;
    }
  
    static void computeCmirViolationAndEfficacy(const std::vector<std::pair<unsigned int,double>>& lhs_cMIR, double rhs_cMIR, double sCoef, double sStar, const double* xlp, int n, double& viol, double& effi, double eps) {
      const double act = activityApprox(lhs_cMIR, xlp, n);
      double v = (act - sCoef * sStar) - rhs_cMIR;   // ≤-Ungleichung
      if (v < eps) v = 0.0;

      const double lhsNorm = coeffNorm2Approx(lhs_cMIR);
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

    struct VecHash {
    size_t operator()(std::vector<GFkSolveClass::SolutionEntry> const& v) const {
      size_t h = v.size();
      for (auto const& e : v) {
	h ^= std::hash<int>()(e.index + 1315423911 * e.weight);
      }
      return h;
    }
  };	  
  struct VecEq {
    bool operator()(const std::vector<GFkSolveClass::SolutionEntry>& a,
		    const std::vector<GFkSolveClass::SolutionEntry>& b) const {
      if (a.size() != b.size()) return false;
      for (size_t i = 0; i < a.size(); ++i) {
	if (a[i].index  != b[i].index)  return false;
	if (a[i].weight != b[i].weight) return false;
      }
      return true;
    }
  };

  static inline int64_t nearestInteger(double x) {
    return static_cast<int64_t>(std::llround(x));
  }
  struct EqTransformForModk {
    // Integer part (für GF(k)-System)
    std::vector<int>    intInds;   // nur binär/int/später ggf. slack-int
    std::vector<double> intVals;   // Koeffizienten (noch unskaliert)
    double rhsAfterSubst;

    // S-Term (für späteren Cut-Bau / expandSTermThroughSlacks)
    std::vector<std::pair<unsigned int,double>> sTerms;  // kann auch Slack-Indizes enthalten
    double sStar;                                // Wert in LP-Lösung (optional)

    // Für “solval”-ähnlichen Filter im Row-Loop
    std::vector<double> intDistance; // distance-to-chosen-bound pro intInd (optional)
    bool ok;
  };
  static inline EqTransformForModk transformEqForModkRow(
							    int row,
							    const std::vector<std::vector<data::IndexedElement>>& allrows,
							    const std::vector<data::QpRhs>& rhsVec,
							    const std::vector<double>& xlp,
							    const std::vector<double>& xlpxtra,
							    const std::vector<data::QpNum>& colLowerBound,
							    const std::vector<data::QpNum>& colUpperBound,
							    const std::vector<std::pair<int, double>>& gclbs,
							    const std::vector<std::pair<int, double>>& gcubs,
							    int *types,
							    const int *eass,
							    double feastol,
							    int n,
							    bool integralPositive)  // analog zu deinem Flag
  {
    EqTransformForModk out;

    // 1) Ausgangszeile (<=) holen
    //    Annahme (wie in deinem Code): rhsVec[row] ist rowUpper
    double rhs = rhsVec[row].getValue().asDouble();

    // 2) Als EQ interpretieren: a x + s = rhs
    //    Slack-Index ist n + row. (Das passt zu deinen Slack-Hilfsfunktionen)
    //    Wir fügen den Slack hier NUR in sTerms ein, nicht in intPart.
    //    (Für GF(k) willst du üblicherweise nur "echte" Integer-Variablen.)
    std::vector<std::pair<unsigned int,double>> eqRow;
    for (int ii=0; ii < allrows[row].size();ii++)
      eqRow.emplace_back(allrows[row][ii].index,allrows[row][ii].value.asDouble());
    // //eqRow.emplace_back(n + row, 1.0);

    // 3) Jetzt Boundsubstitution + Split in Integer-Part und S-Term
    //    -> Wir nutzen deine EQ-Routine.
    //    ACHTUNG: Ich kann die exakte Signatur in deiner Datei hier nicht
    //    100% wortgleich garantieren, aber inhaltlich ist es genau dieser Call.
    std::vector<std::pair<unsigned int,double>> rowIntPart;
    std::vector<std::pair<unsigned int,double>> contInS;
    double rhsAfterSubst = rhs;
    double sStar = 0.0;

    // Pseudocode: ersetze durch deine echte Signatur / Parameterreihenfolge
    bool ok = substituteBoundsAndSplitIntS_EQ(
					      /*row=*/eqRow,
					      /*xlp=*/xlp.data(), xlpxtra.data(),
					      /*colUpperBound etc...*/colUpperBound.data(),colLowerBound.data(),
					      /*outputs*/ rowIntPart, rhsAfterSubst, sStar, contInS,
					      n, types, gclbs, gcubs, xlpxtra.size()
					      );

    if (!ok) {
      out.ok = false;
      return out;
    }

    out.rhsAfterSubst = rhsAfterSubst;
    out.sStar = sStar;
    out.sTerms = std::move(contInS);  // inkl. Slack-Anteil, falls er in S bleibt

    // 4) GF(k)-Part: nur Binär/Integer (bei dir: types[col]==0 ist binär)
    //    Dazu bauen wir intInds/intVals und ersetzen "solval" durch Distanz.
    for (const auto& [col, a] : rowIntPart) {
      if (col < 0) continue;
      if (col >= n) {
	// Slack im intPart? im Regelfall ignorieren für GF(k)
	// (kann man später anders machen, aber erstmal robust)
	continue;
      }
      if (types[col] != 0) continue; // nur binär in deinem Setting

      // Distance-to-chosen-bound als solval-Ersatz:
      // Wir nehmen tightened bounds (gclbs/gcubs), weil "VLB/VUB ist automatisch drin".
      //double lb = gclbs[col];
      //double ub = gcubs[col];
      auto VLB = gclbs[col];
      double lb = (VLB.first >= 0) ? VLB.second * xlp[VLB.first]
	: colLowerBound[col].asDouble();
      
      auto VUB = gcubs[col];
      double ub = (VUB.first >= 0) ? VUB.second * xlp[VUB.first]
	: colUpperBound[col].asDouble();

      // Für Binär sollte das [0,1] sein; bei Fixierungen ggf. lb==ub.
      // Wir wählen die Bound-Seite:
      // - Wenn integralPositive: lenke negative Koeffizienten über "ub" etc.
      //   (du kannst hier exakt deine bisherige boundType-Logik nachbauen)
      bool chooseLb = true;
      if (integralPositive) {
	// grobe Default-Regel: positive Koeff -> wähle lb, negative -> wähle ub
	chooseLb = (a >= 0.0);
      } else {
	// "nearest-bound" (stabiler Default)
	double dlb = fabs(xlp[col] - lb);
	double dub = fabs(ub - xlp[col]);
	chooseLb = (dlb <= dub);
      }

      // Distanz wie in selectRowToAggregateEQ:
      double dist = fmin(xlp[col] - lb, ub - xlp[col]);
      
      // numerische Bereinigung
      //if (dist < 0.0 && dist > -10*feastol) dist = 0.0;
      if (dist < 0.0) dist = 0.0;

      out.intInds.push_back(col);
      out.intVals.push_back(a);
      out.intDistance.push_back(dist);
    }

    out.ok = true;
    return out;
  }

// helper: current bounds exactly like selectRowToAggregateEQ
static inline void getCurrentBoundsEQStyle(
					   int col,
					   const std::vector<std::pair<int,double>>& gclbs,
					   const std::vector<std::pair<int,double>>& gcubs,
					   const std::vector<data::QpNum>& colLowerBound,
					   const std::vector<data::QpNum>& colUpperBound,
					   const std::vector<double>& xlp,
					   double& LoB,
					   double& UpB)
  {
    const auto VLB = gclbs[col];
    LoB = (VLB.first >= 0) ? (VLB.second * xlp[VLB.first]) : colLowerBound[col].asDouble();

    const auto VUB = gcubs[col];
    UpB = (VUB.first >= 0) ? (VUB.second * xlp[VUB.first]) : colUpperBound[col].asDouble();
  }


  static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >*getMODKcuts( extSol::QpExternSolver &extSolver, std::vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev, unsigned int initime, std::vector<int> &listOfCutsVars, int*solu, int *fixs, int* blcks, int* eass, int orgN );
  static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >*getMODKcutsProgressive( extSol::QpExternSolver &extSolver, std::vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev, unsigned int initime, std::vector<int> &listOfCutsVars, int*solu, int *fixs, int* blcks, int* eass, int orgN )
 {
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
    return getMIRsmartFlexiCore( extSolver, allrows, rhsVec, types, assigns, decLev, initime, solu, fixs, blcks, eass, orgN, 2, CutAdder::AggregationMode::ModK_GFk);    
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
    //return getMIRsmartFlexiCore( extSolver, allrows, rhsVec, types, assigns, decLev, initime, solu, fixs, blcks, eass, orgN, 2, CutAdder::AggregationMode::HeuristicAggregation);
    return getMIRsmartCore( extSolver, allrows, rhsVec, types, assigns, decLev, initime, solu, fixs, blcks, eass, orgN, 2, false);    
  }

  static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >* getMIRsmartCore( extSol::QpExternSolver &extSolver, std::vector<std::vector<data::IndexedElement>> &allRows, std::vector<data::QpRhs> &rhsVec, int *types, int8_t *assigns, int decLev, unsigned int initime, int * solu, int* fixs, int *blcks, int *eass, int orgN, int maxAggr, bool fromCmirize );
  static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >* getMIRsmartFlexiCore( extSol::QpExternSolver &extSolver, std::vector<std::vector<data::IndexedElement>> &allRows, std::vector<data::QpRhs> &rhsVec, int *types, int8_t *assigns, int decLev, unsigned int initime, int * solu, int* fixs, int *blcks, int *eass, int orgN, int maxAggr, CutAdder::AggregationMode aggMode );
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
    if (fabs(coefCol) < eps && indCol < numCols) {
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
	    if (tightenedUpBnds.size() > varix) VUB = tightenedUpBnds[varix];
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
            const int idx = listRowsAggregated[varix - n];
            const double mult = (rhsVec.at(idx).getRatioSign()==data::QpRhs::smallerThanOrEqual)
                                ? (- sCoefBestCut * varco)
                                : (  sCoefBestCut * varco);
            bestCutRhs += rhsVec.at(idx).getValue().asDouble() * mult;
	    const std::vector<data::IndexedElement>& row = allRows[idx];
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

  /*
  for (int ii = 0; ii < (int)rowAggregated.size(); ++ii) {
    const unsigned idx = rowAggregated[ii].first;
    if (idx >= (unsigned)(n + numSlacks)) {
      std::cerr << "CMIR dense overflow: idx=" << idx
		<< " n=" << n << " numSlacks=" << numSlacks
		<< " (n+numSlacks=" << (n+numSlacks) << ")\n";
      abort();
    }
  }
  */
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

    return getMIRsmartFlexiCore( extSolver, allrows, rhsVec, types, assigns, decLev, initime, solu, fixs, blcks, eass, orgN, 1, /*true*/CutAdder::AggregationMode::FromCMIRize);    
}
// end of MIRsmart
  
  static void aggregate2RowsInFirst(std::vector<std::pair<unsigned int, double>> &aggRow /*must be sorted by ascending index*/, double &rhs_aggRow, const std::vector<std::pair<unsigned int, double>> &b, double rhs_b, double w, double tol, int pivot);

  static bool determineCover(CoverGenState& st, bool lpSol=true);
  static void separateLiftedKnapsackCover(CoverGenState& st, /*out cut*/...);
  static bool separateLiftedMixedBinaryCover(CoverGenState& st, /*out cut*/...);
  static bool separateLiftedMixedIntegerCover(CoverGenState& st, /*out cut*/...);
  
  static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >* getCoverCuts( extSol::QpExternSolver &extSolver, int *types, int8_t *assigns, int decLev, unsigned int, int *, int *, int*, int*, int );
	static std::vector< std::pair< std::vector< std::pair<unsigned int, double> >, double > >* getLPCuts(extSol::QpExternSolver &extSolver, std::vector<unsigned int> candidateVariables, int *types, int8_t *assigns, int decLev, unsigned int initime, std::vector<int> &listOfCutsVars);

	static void setOrgM(int m) { orgM = m; }
	static int getOrgM() { return orgM; }
private:
	static int orgM;
	static Param_c param;
        static double random_seed;
	static const int info_level = 1;
        CutAdder();
	static inline bool isZero(double x, double epsZero = 1e-20) {
		return (fabs(x) <= epsZero);
	}
        static bool getAllRows(extSol::QpExternSolver &extSolver, std::vector< std::vector< data::IndexedElement >  > &allrows2);

	static bool createCover( std::vector<double>& xlpopt, const std::vector<std::pair<double, unsigned int> >& rowsparse, double rhs, const double eps, std::vector<std::pair<unsigned int, double> >& reslhs, double& resrhs, int *types,  int8_t *assigns, int decLev, data::QpNum*, data::QpNum*, int *solu, int*fixs, int *blcks, int orgN );
};

#endif /* CUTADDER_H_ */
