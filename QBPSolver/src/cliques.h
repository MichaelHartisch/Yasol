/*
*
* Yasol: cliques.h -- Copyright (c) 2012-2017 Ulf Lorenz
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

#include "multiset_plus.h"
#include <string>
#include <algorithm>

#include "Settings/Settings.hpp"
#include "Algorithms/Algorithms.hpp"
#include "Datastructures/Datastructures.hpp"
#include "Utilities/Parser.hpp"
#include "DataStructures.h"
#include "HashTable.h"

#define USE_CLICK_EXTR_OUT 0


class ValueConstraintPair {
public:
	 CRef cr;
	 Var  v;
	 int  pos;
	 ValueConstraintPair(CRef lcr, Var lv, int lpos) : cr(lcr), v(lv), pos(lpos) {}
	 ValueConstraintPair() {}
};

class CIpair {
public:
	int32_t constraint;
	int32_t index;
	CIpair(int32_t c, int32_t i) {
		constraint=c;index=i;
	}
	CIpair() { constraint = index = -1; }
};

class GUBelement {
public:
	Var x;
	CIpair N,S;
	int32_t W,E;
	GUBelement(int v) {
		x = v;
		N.constraint = S.constraint = -1;
		W = E = -1;
	}
};

class _GUBGraph {
	ca_vec<CIpair> head_list;
	//ca_vec< ca_vec<GUBelement> > graph;
	ca_vec<bool> sharp;
public:
	void InitGUBGraph(int n) {
		head_list.clear();
		for (int i=0; i < n;i++)
			head_list.push(CIpair(-1,-1));
	}
	/*void AddGUBconstraint(ca_vec<CoeVar> &constraint, coef_t rhs) {
		graph.push();
		sharp.push(false);
		int ci = graph.size()-1;
		for (int i = 0; i < constraint.size();i++) {
			graph[ci].push(GUBelement(constraint[i].x));
			if (i >= 1) {
				graph[ci][i-1].E = i;
				graph[ci][i].W = i-1;
			}
			graph[ci][i].S  = head_list[graph[ci][i].x];
			head_list[graph[ci][i].x] = CIpair(ci,i);
			if (graph[ci][i].S.constraint > -1) graph[ graph[ci][i].S.constraint ][ graph[ci][i].S.index ].N = CIpair(ci,i);
		}
		return;
		cerr << "*";
	}
	bool assignAndScanGUBgraph(int var, int val, ca_vec<ValueConstraintPair> &propQ, ca_vec<int> &propQlimiter,
			                    ca_vec<int> &litList, ca_vec<extbool> &assigns) {
        int lit;
		int olit = var+var;
		if (val == 0) olit |= 1;
		int row,row2,col;
		int ix1,ix2,conflict=-1;
		litList.clear();
		litList.push(olit);
		while (litList.size() > 0) {
			lit = litList[litList.size()-1];
			litList.pop();
			for (row = head_list[lit].constraint, col = head_list[lit].index;
					row != -1;
					row2 = row, row = graph[row][col].S.constraint, col = graph[row2][col].S.index) {
				if (sharp[row]) {
					//conflict. Learn!
				}
				sharp[row] = true;
				for (int left = col - 1; left != -1; left = graph[row][left].W) {
					if (graph[row][left].x != olit && assigns[graph[row][left].x>>1] == extbool_Undef && propQlimiter[graph[row][left].x ^ 1] <= 0) {
						propQ.push( ValueConstraintPair(CRef_Undef, graph[row][left].x ^ 1, left));
						propQlimiter[graph[row][left].x ^ 1] = propQ.size();
						if (propQlimiter[graph[row][left].x] > 0) {
							conflict = row;
							ix1 = propQlimiter[graph[row][left].x ^ 1] - 1;
							ix2 = propQlimiter[graph[row][left].x] - 1;
							break;
						}
						litList.push(graph[row][left].x ^ 1);
					}
				}
				if (conflict < 0) {
					for (int right = col + 1; right != -1; right = graph[row][right].E) {
						if (graph[row][right].x != olit && assigns[graph[row][right].x>>1] == extbool_Undef && propQlimiter[graph[row][right].x ^ 1] <= 0) {
							propQ.push(ValueConstraintPair(CRef_Undef, graph[row][right].x ^ 1, right));
							propQlimiter[graph[row][right].x ^ 1] =
									propQ.size();
							if (propQlimiter[graph[row][right].x] > 0) {
								conflict = row;
								ix1 = propQlimiter[graph[row][right].x ^ 1] - 1;
								ix2 = propQlimiter[graph[row][right].x] - 1;
								break;
							}
							litList.push(graph[row][right].x ^ 1);
						}
					}
				}
				if (conflict >= 0) {
					// conflict. Learn!
				}
			}
		}
	}
    bool isGUBConstraint(ca_vec<CoeVar> &constraint, coef_t rhs) {
    	float cnt_neg=0.0;
    	for (int i=0; i < constraint.size();i++) {
    		if (constraint[i].coef >= 1.0+COEF_EPS || constraint[i].coef <= 1.0-COEF_EPS) return false;
    		if (sign(constraint[i])) cnt_neg = cnt_neg + 1.0;
    	}
    	if (rhs >= 1.0+cnt_neg+COEF_EPS || rhs <= 1.0+cnt_neg-COEF_EPS) return false;
    	return true;
    }*/
};

class ConflictGraphEdge {
public:
	union { struct { uint32_t v1,v2; } i32; uint64_t v12; } e;
	uint32_t getv1() { return e.i32.v1; }
	uint32_t getv2() { return e.i32.v2; }
	uint64_t gete()  { return e.v12; }
	ConflictGraphEdge(int32_t v1, uint32_t v2) { e.i32.v1 = v1; e.i32.v2 = v2; }
	ConflictGraphEdge(uint64_t edge) { e.v12 = edge; }
};

class _ConflictGraph {
public:
	ca_vec<int32_t> last;
	ca_vec<int32_t> adj;
	ca_vec<int32_t> next;
	int FirstAdjacent(int i) {
	  if (i > -1)
		return last[i];
	  else return -1;
	}
	int NextAdjacent(int j) {
		if (j > -1) return next[j];
		else return -1;
	}
	void AddEdge(int i, int j) {
		adj.push(j);
		next.push(last[i]);
		last[i] = adj.size()-1;
		adj.push(i);
		next.push(last[j]);
		last[j] = adj.size()-1;
	}
	void DelLastEdge(int &i, int &j) {
	        if (adj.size() < 2) {
		  assert (adj.size() == 0);
		  i = j = -1;
		  return ;
	        }     
	        // deleted is edge (i,j)
	        j = adj[adj.size()-2];
	        i = adj[adj.size()-1];
	        last[ j  ] = next[ last[j]  ];
		last[ i  ] = next[ last[i]  ];
		next.pop();
		next.pop();
		adj.pop();
		adj.pop();
	}
	void InitConflictGraph(int n) {
		adj.clear();
		next.clear();
		last.clear();
		for (int i=0; i < n; i++)
			last.push(-1);
	}
	int getSize() {
		return last.size();
	}
};

/*class _EdgeContainer : public AVLCntContainer<uint64_t>
{
public:
	NodeType<uint64_t>* iter;
	bool isEmpty() {
		if (this->size == 0) return true;
		return false;
	}
	uint64_t getFstRetData() {
		iter = this->begin();
		return iter->data;
	}
	bool isFinished() {
		if (iter == this->end()) return true;
		return false;
	}
	uint64_t getNxtRetData() {
		iter = this->next;
		return iter->data;
	}
}*/

class CliqueManager {
	int nVars;
    _GUBGraph GUBGraph;
	_ConflictGraph ConflictGraph;
	AVLCntContainer<uint64_t> EdgeContainer;
    // int Size() { return size; }
    // void Insert(T v, int q) { _Insert(first, v, q); }
    // void Remove(T v, int q) { _Remove(first,v,q); }
    // void Insert(T v) { _Insert(first,v,1); }
    // void Remove(T v) { _Remove(first,v,1); }
    // int IsContained(T v) { return _IsContained(first, v); }
    // void Clear() { _Clear(first); }
    // int IsEmpty() { if (first) return 0; else return 1; }
public:
	struct VarData { CRef reason; int level; };
	ca_vec<Cli_Table> cliques;
	CliqueManager() {

	}
    void CliqueManagerInitGraph(int n) {
       nVars=n;
       ConflictGraph.InitConflictGraph(n+n);
       GUBGraph.InitGUBGraph(n+n);
    }

    int getConflictGraphSize() {
    	return ConflictGraph.getSize();
    }

    void AddEdge(uint32_t l1, uint32_t l2) {
    	if (!EdgeIsInContainer(l1, l2)) {
    		//std::cerr << "E2C(" << l1 << "," << l2 << ")";
    		AddEdge2Container(l1,l2);
    		//assert(EdgeIsInContainer(l1, l2));
    		ConflictGraph.AddEdge(l1, l2);
    		/*bool isIn = false;
    		int next = FirstAdjacentInConflictGraph(l1);
    		while (next >= 0) {
    			int node = getAdjacent(next);
    			if (node == l2) isIn = true;
    			next = NextAdjacentInConflictGraph(next);
    		}
    		assert(isIn);
    		isIn = false;
    		next = FirstAdjacentInConflictGraph(l2);
    		while (next >= 0) {
    			int node = getAdjacent(next);
    			if (node == l1) isIn = true;
    			next = NextAdjacentInConflictGraph(next);
    		}
    		assert(isIn);
    		*/
    	}
    }
    void AddEdge2Container(uint32_t l1, uint32_t l2) {
    	ConflictGraphEdge ce(l1,l2);
    	EdgeContainer.Insert(ce.gete());
    }
    void DelEdgeFromContainer(uint32_t &l1, uint32_t &l2) {
        int i,j;
	ConflictGraph.DelLastEdge(i,j);
	l1 = (uint32_t)i;
	l2 = (uint32_t)j;
    	ConflictGraphEdge ce(l1,l2);
    	if (EdgeContainer.IsContained(ce.gete())) {
	  EdgeContainer.Remove(ce.gete());
	  assert(!EdgeContainer.IsContained(ce.gete()));
	}
    }

    bool EdgeIsInContainer(uint32_t l1, uint32_t l2) {
    	ConflictGraphEdge ce(l1,l2);
    	return EdgeContainer.IsContained(ce.gete());
    }
    int getAdjacent(int idx) {
    	return ConflictGraph.adj[idx];
    }
	int FirstAdjacentInConflictGraph(int i) {
		return ConflictGraph.FirstAdjacent(i);
	}
	int NextAdjacentInConflictGraph(int j) {
		return ConflictGraph.NextAdjacent(j);
	}

    /*void addGUBConstraint(ca_vec<CoeVar> &constraint, coef_t rhs) {
        if (GUBGraph.isGUBConstraint(constraint, rhs)) {
        	GUBGraph.AddGUBconstraint(constraint, rhs);
        }
    }
    bool isGUBConstraint(ca_vec<CoeVar> &constraint, coef_t rhs) {
    	float cnt_neg=0.0;
    	static int cnt=0;
    	for (int i=0; i < constraint.size();i++) {
    		if (constraint[i].coef >= 1.0+COEF_EPS || constraint[i].coef <= 1.0-COEF_EPS) return false;
    		if (sign(constraint[i])) cnt_neg = cnt_neg + 1.0;
    	}
    	if (rhs >= 1.0+cnt_neg+COEF_EPS || rhs <= 1.0+cnt_neg-COEF_EPS) return false;
    	cnt++;
    	//
    	//cerr << cnt << ": ";
    	//for (int i=0; i < constraint.size();i++) {
    	//	cerr << (sign(constraint[i]) ? "-" : "") << constraint[i].coef << "x" << i << " ";
    	//}
    	//cerr << rhs << endl;
    	return true;
    }
    */
    void extractCliqueFromLPconstraint(ca_vec<CoeVar> &constraint_lhs, coef_t rhs, ca_vec<CoeVar> &clique, ca_vec<CoeVar> &zero_forced) {
    	// convert from >= to <=
    	//SearchOrderLt SOL;
    	int i;
    	if (USE_CLICK_EXTR_OUT) cerr << "enter extraction" << endl;
    	//static ca_vec<CoeVar> constraint_lhs;
    	//constraint_lhs.clear();
    	//for (int j = 0; j < c_lhs.size();j++) constraint_lhs.push(c_lhs[j]);
    	clique.clear();
    	zero_forced.clear();
    	rhs = -rhs;
    	for (i = 0; i < constraint_lhs.size();i++) {
    		coef_t k = constraint_lhs[i].coef;
    		assert(k >= 0.0);
    		if (!sign(constraint_lhs[i])) rhs = rhs + k;
    	}
    	// convert a copy of constraint to knapsack
    	// constraint_lhs itself contains the knapsack constraint. sign(x_i) == 0 means x_i is replaced by (1-x'_i)

    	// sort it non-increasing --> alle Koeffizenten sind positiv, wegen Konvertierung
    	// MUST BE SORTED BEVORE!!! sort(constraint_lhs,SOL);
    	// solange Koeffizent > rechte Seite --> zugeh�rige Variable = 0
    	for (i = 0; i < constraint_lhs.size();i++) {
    		if (constraint_lhs[i].coef > rhs) zero_forced.push(constraint_lhs[i]);
    		else break;
    	}
    	// solange x_i == x_i+1 == 1 => gr��er rechte Seite: geh�ren beide zu Clique.
    	bool lastHasVal=false;
    	CoeVar last;
    	for ( ; i < constraint_lhs.size()-1;i++) {
    		CoeVar Z;
    		if (constraint_lhs[i].coef + constraint_lhs[i+1].coef > rhs) {
    			Z.x = constraint_lhs[i].x^1;
    			Z.coef = 1.0;
    			clique.push(Z);
    			lastHasVal = true;
    			last.x = constraint_lhs[i+1].x^1;
    			last.coef = 1.0;
    		}
    		else break;
    	}
    	if (lastHasVal) {
    		clique.push(last);
    	}
    	if (USE_CLICK_EXTR_OUT) cerr << "leave extraction" << endl;
    }

    void addNativeClique(ca_vec<CoeVar> &C, int numVars) {
    	Cli_Table *CT = new Cli_Table(numVars, (C.size() * 3) / 2);
    	cliques.push(*CT);
    	for (int i = 0; i < C.size();i++) {
    		uint64_t H = CT->getHashConstant(var(C[i])+var(C[i]));
    		CT->setEntry(H, C[i]);
    	}
    }

    void extractImplis(ca_vec<CoeVar> &propQ, int v, int p, int n,ca_vec<extbool> &assigns) {
    	propQ.clear();
    	std::vector<bool> used(n+n,false);
    	int next = FirstAdjacentInConflictGraph(v+v+(1-p));
    	used[v+v+(1-p)] = true;
    	while (next >=0) {
    		int node = getAdjacent(next);
    		if (used[node] == false && assigns[node>>1] == extbool_Undef) {
    		    used[node] = true;
    		    propQ.push(mkCoeVar(node>>1,1.0,(node&1) ? true : false));
    		    //dfsExtractImplis(propQ, node, used,assigns);
    		    dfsExtractImplis(propQ, node^1, used,assigns);
    		}
    		next = NextAdjacentInConflictGraph(next);
    	}
    }
    void dfsExtractImplis(ca_vec<CoeVar> &propQ, int node, std::vector<bool> &used,ca_vec<extbool> &assigns) {
    	if (used[node] == true) return;
    	int next = FirstAdjacentInConflictGraph(node);
    	used[node] = true;
    	while (next >=0) {
    		int node = getAdjacent(next);
    		if (used[node] == false && assigns[node>>1] == extbool_Undef) {
    		    used[node] = true;
    		    propQ.push(mkCoeVar(node>>1,1.0,(node&1) ? true : false));
    		    //dfsExtractImplis(propQ, node, used,assigns);
    		    dfsExtractImplis(propQ, node^1, used,assigns);
    		}
    		next = NextAdjacentInConflictGraph(next);
    	}
    }
    int forcedByConflictTable(int v, int n,ca_vec<extbool> &assigns, int &level, VarData *vd, int &sigvar) {
    	int p=0;
    	int mustnot0=false;
    	int mustnot1=false;
    	level = n+5;
    	int next = FirstAdjacentInConflictGraph(v+v+(1-p));
    	while (next >=0) {
    		int node = getAdjacent(next);
    		if (assigns[node>>1] != extbool_Undef) {
    			if (node&1) { // Verbot hat Vorzeichen gesetzt => ist not-ed
    				if (assigns[node>>1] == 0) {
    					mustnot0 = true;
    	    			if (vd[node>>1].level < level) {
    	    				level = vd[node>>1].level;
    	    				sigvar = node;
    	    			}
    				}
    			} else {
    				if (assigns[node>>1] == 1) {
    					mustnot0 = true;
    	    			if (vd[node>>1].level < level) {
    	    				level = vd[node>>1].level;
    	    				sigvar = node;
    	    			}
    				}
    			}
    		}
    		next = NextAdjacentInConflictGraph(next);
    	}
    	p=1;
    	next = FirstAdjacentInConflictGraph(v+v+(1-p));
    	while (next >=0) {
    		int node = getAdjacent(next);
    		if (assigns[node>>1] != extbool_Undef) {
    			if (node&1) { // Verbot hat Vorzeichen gesetzt => ist not-ed
    				if (assigns[node>>1] == 0) {
    					mustnot1 = true;
    	    			if (vd[node>>1].level < level) {
    	    				level = vd[node>>1].level;
    	    				sigvar = node;
    	    			}
    				}
    			} else {
    				if (assigns[node>>1] == 1) {
    					mustnot1 = true;
    	    			if (vd[node>>1].level < level) {
    	    				level = vd[node>>1].level;
    	    				sigvar = node;
    	    			}
    				}
    			}
    		}
    		next = NextAdjacentInConflictGraph(next);
    	}
    	if (mustnot1 && mustnot0) return 4;
    	if (mustnot1 == true) return 0;
    	if (mustnot0 == true) return 1;
    	return extbool_Undef;
    }

    bool checkTheGraph(ca_vec<int> &optSol) {
    	return true;
    	for (int z=0; z < optSol.size();z++) {
    		assert(optSol[z] == 0 || optSol[z] == 1);
    		if (optSol[z] == 1) {
            	int next = FirstAdjacentInConflictGraph(z+z);
            	while (next >=0) {
            		int node = getAdjacent(next);
            		int realVal = optSol[node>>1];
            		int notAllowed = ((node&1) ? 0 : 1);
            		if (realVal == notAllowed) return false;
            		next = NextAdjacentInConflictGraph(next);
            	}
    		} else {
            	int next = FirstAdjacentInConflictGraph(z+z+1);
            	while (next >=0) {
            		int node = getAdjacent(next);
            		int realVal = optSol[node>>1];
            		int notAllowed = ((node&1) ? 0 : 1);
              		if (realVal == notAllowed) return false;
            		next = NextAdjacentInConflictGraph(next);
            	}
    		}
    	}
    	return true;
    }
};
