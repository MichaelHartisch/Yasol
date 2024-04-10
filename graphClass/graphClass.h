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

#ifndef GRAPHCLASS_H
#define GRAPHCLASS_H

#include "multiset_plus.h"
#include <string>
#include <math.h>
#include <algorithm>

class GeneralGraphEdge {
 public:
  union { struct { uint32_t v1,v2; } i32; uint64_t v12; } e;
  uint32_t getv1() { return e.i32.v1; }
  uint32_t getv2() { return e.i32.v2; }
  uint64_t gete()  { return e.v12; }
  GeneralGraphEdge(int32_t v1, uint32_t v2) { e.i32.v1 = v1; e.i32.v2 = v2; }
  GeneralGraphEdge(uint64_t edge) { e.v12 = edge; }
};

class _GeneralGraph {
 public:
  ca_vec<int> listStart;
  ca_vec<int> adj;
  ca_vec<int> next;
  ca_vec<int> prev;
  int nextFree;
  int edgeCnt;
  bool out;
  _GeneralGraph() {
    InitGraph(10,false);
  }
  int FirstAdjacent(int i) {
    if (i > -1) {
      if (i >= getNodeCnt()) {
	for (int j = getNodeCnt(); j <= i;j++)
	  listStart.push(-1);
      } 
      return listStart[i];
    } else return -1;
  }
  int NextAdjacent(int j) {
    if (j > -1) return next[j];
    else return -1;
  }
  int PrevAdjacent(int j) {
    if (j != listStart[j]) return prev[j];
    else return -1;
  }
  int getNewIndex4Adj() {
    if (nextFree == -1) {
      adj.push(-1);
      next.push(-1);
      prev.push(-1);
      return prev.size()-1;
    } else {
      int ix = nextFree;
      nextFree = prev[nextFree];
      return ix;
    }
  }
  void AddEdge(int i, int j) {
    int ix_i = getNewIndex4Adj();
    adj[ix_i] = i;
    next[ix_i] = listStart[j];
    //the following if (listStart[j] >= 0) is added because runtimt error. Why was this "forgotten"? -> unclear
    if (listStart[j] >= 0) prev[listStart[j]] = ix_i;
    listStart[j] = ix_i;
    prev[listStart[j]] = -1;
    edgeCnt++;
  }

  int findAddress(int i, int j) {
    int ix = listStart[i];
    while (adj[ix] != j) {
      ix = next[ix];
      if (ix < 0) return -1; //is not contained
    }
    return ix;
  }
  void DeleteEdge(int i, int j) {
    if (out) cerr << "try to delete directed edge (" << i << "," << j << ")" << endl; 
    int i_ix = findAddress(i,j);
    if (i_ix < 0) {
      if (out) cerr << "Warning: Edge " << " does not exist." << endl;
      return;
    }
    int j_ix = i_ix + 1;
    for (int h = 1; h < 2;h++) {
      int ix;
      int v;
      edgeCnt--;
      if (h == 0) { ix = j_ix; v = j; }
      else        { ix = i_ix; v = i; }
      if (next[ix] < 0 && prev[ix] < 0) { // weder prev noch next
	if (nextFree < 0) nextFree = ix;
	else {
	  prev[ix] = nextFree;
	  nextFree = ix;
	}
	listStart[v] = -1;
      } else if (prev[ix] < 0) { //kein prev, aber next
	prev[next[ix]] = -1;
	listStart[v] = next[ix];
      } else if (next[ix] < 0) { //kein next, aber prev
	next[prev[ix]] = -1;
      } else {
	next[prev[ix]] = next[ix];
	prev[next[ix]] = prev[ix];
      }
    }
  }

  void InitGraph(int n, bool lout) {
    nextFree = -1;
    edgeCnt = 0;
    out = lout;
    adj.clear();
    next.clear();
    listStart.clear();
    for (int i=0; i < n; i++)
      listStart.push(-1);
  }
  int getNodeCnt() {
    return listStart.size();
  }
  int getEdgeCnt() {
    return edgeCnt;
  }
  int getEdgeSpace() {
    return adj.size();
  }

  bool getEdge(int e, int &i, int &j) {
    if (e >= 0 && e < getEdgeSpace()) {
      bool found = false;

      for (int jjj=0;!found && jjj<getNodeCnt();jjj++) {
	int jj;
	jj = FirstAdjacent(jjj);
	if (jj >= 0) {
	  int w = adj[jj];
	  int ex = -1;
	  //cerr << "1) found (" << jjj << "," << w << ") " << endl;                  
	  ex = findAddress(jjj,w);
	  if (ex == e) {
	    if (out) cerr << "FOUND EDGE " << e << endl;
	    found = true;
	    i = jjj;
	    j = w;
	    break;
	  }   
	  int kk = NextAdjacent(jj);
	  while (kk >= 0) {
	    int w = adj[kk];
	    //cerr << "2) found (" << jjj << "," << w << ") " << endl;
	    ex = findAddress(jjj,w);
	    if (ex == e) {
	      if (out) cerr << "FOUND EDGE " << e << endl;
	      found = true;
	      i = jjj;
	      j = w;
	      break;
	    }   
	    kk = NextAdjacent(kk);
	  }
	}
      }
      if (found) return true;
      else return false;
    }
    return false;
  }
};

class GraphManager {
  int numNodes;
  bool out;
  _GeneralGraph GeneralGraph;
  AVLCntContainer<uint64_t> EdgeContainer;
  // int Size() { return size; }
  // void Insert(T v, int q) { _Insert(first, v, q); }
  // void Remove(T v, int q) { _Remove(first,v,q); }
  // void Insert(T v) { _Insert(first,v,1); }
  // void Remove(T v) { _Remove(first,v,1); }
  // int IsContained(T v) { return _IsContained(first, v); }
  // void Clear() { _Clear(first); }
  // int IsEmpty() { if (first) return 0; else return 1; }

  double random_seed;
  // Returns a random double 0 <= x < 1. Seed must never be 0.
  static inline double drand(double& seed) {
    seed *= 1389796;
    int q = (int)(seed / 2147483647);
    seed -= (double)q * 2147483647;
    return seed / 2147483647; }

  // Returns a random integer 0 <= x < size. Seed must never be 0.
  static inline int irand(double& seed, int size) {
    return (int)(drand(seed) * size); }

 public:
  GraphManager() {
    random_seed = 2375.0; //(double)time(NULL);
  }
  void GraphManagerInitGraph(int n, bool lout) {
    numNodes=n;
    out = lout;
    GeneralGraph.InitGraph(n,lout);
  }

  int getNodeCnt() {
    return GeneralGraph.getNodeCnt();
  }

  bool AddEdge(int l1, int l2) {
    GeneralGraphEdge ce(l1,l2);
    EdgeContainer.IsContained(ce.gete());
    if (!EdgeContainer.IsContained(ce.gete())) {
      //std::cerr << "E2C(" << l1 << "," << l2 << ")";
      EdgeContainer.Insert(ce.gete());
      GeneralGraph.AddEdge(l2, l1);
      return true;
    } else return false;
  }

  bool DeleteEdge(int l1, int l2) {
    GeneralGraphEdge ce(l1,l2);
    EdgeContainer.IsContained(ce.gete());
    if (EdgeContainer.IsContained(ce.gete())) {
      //std::cerr << "E2C(" << l1 << "," << l2 << ")";
      EdgeContainer.Remove(ce.gete());
      GeneralGraph.DeleteEdge(l1, l2);
      return true;
    } else return false;
  }

  bool getEdgeByIndex(int e, int &i, int &j) {
     return GeneralGraph.getEdge(e,i,j);
  }

  int getAdjacent(int idx) {
    return GeneralGraph.adj[idx];
  }
  int FirstAdjacentInGraph(int i) {
    return GeneralGraph.FirstAdjacent(i);
  }
  int NextAdjacentInGraph(int j) {
    return GeneralGraph.NextAdjacent(j);
  }

  void printGraph() {
    if (!out) return;
    cerr << endl;
    for (int jjj=0;jjj<getNodeCnt();jjj++) {
      int j;
      j = FirstAdjacentInGraph(jjj);
      if (j >= 0) {
	cerr << "Node: " << jjj << " Neighbours: ";
	int w = getAdjacent(j);
	cerr << w << " ";                   
	int kk = NextAdjacentInGraph(j);
	while (kk >= 0) {
	  int w = getAdjacent(kk);
	  cerr << w << " ";
	  kk = NextAdjacentInGraph(kk);
	}
	cerr << endl;
      }
    }
  }

  bool checkConsistency() {
    if (out) cerr << "There are " << GeneralGraph.getNodeCnt() << " nodes." << endl;
    if (out) cerr << "There are " << GeneralGraph.getEdgeCnt() << " edges." << endl;
    if (out) cerr << "There currently is space for  " << GeneralGraph.getEdgeSpace() << " edges." << endl;
    return true;
  }


  void initTestGraph(int n, bool lout) {
    int numEdges=0;
    GraphManager *GM = new GraphManager();
    GM->GraphManagerInitGraph(n, lout);
    int last_i=-1;
    int last_j=-1;
    for (int k = 0; k < (int)((double)n*sqrt(n));k++) {
      int i = irand(random_seed,n);
      int j = irand(random_seed,n);
      if (i != j) {
	if (GM->AddEdge(i,j)) {
	  numEdges++;
	  if (lout) cerr << "Add directed edge (" << i << "," << j << ")" << endl;
	  last_i = i;
	  last_j = j;
	}
	if (GM->AddEdge(i,j)) {
	  numEdges++;
	  cerr << "Error: 2nd attempt should fail: Add directed edge (" << i << "," << j << ")" << endl;
	}
      }
    }    
    GM->printGraph();
    if (last_i > -1) {
      if (lout) cerr << "Delete last Edge (" << last_i << "," << last_j << ")";
      GM->DeleteEdge(last_i, last_j);
    } else if (out) cerr << "Warning: no last edge to be deleted" << endl;
    for (int k=0;k<(int)sqrt(n);k++) {
      int i,j;
      if (lout) cerr << "There are " << GM->GeneralGraph.getEdgeCnt() << " directed edges in the graph. Reserved are " << GM->GeneralGraph.getEdgeSpace() << endl;
      int e = irand(random_seed,GM->GeneralGraph.getEdgeSpace());
      bool success = GM->getEdgeByIndex(e,i,j);
      if (success) {
	if (lout) cerr << "Delete Edge (" << i << "," << j << ")";
	GM->DeleteEdge(i, j);
      } else {
	if (lout) cerr << "Edge " << e <<  " does not exist." << endl; 
      }
    }
    GM->checkConsistency();
    if (lout) GM->printGraph();
  }

};
#endif
