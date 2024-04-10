/*
*
* Yasol: dependencies.h -- Copyright (c) 2012-2017 Ulf Lorenz
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

#ifndef DEPENDENCIES_H_
#define DEPENDENCIES_H_

#include "multiset_plus.h"
#include "UnionFind.h"

using namespace std;
class DepManGraphEdge {
public:
	static const int deleted=1;
	union { struct { unsigned int v1 : 31,v2 : 31; unsigned int flags : 2; } node; uint64_t edge; } e;
	uint32_t getv1() { return e.node.v1; }
	uint32_t getv2() { return e.node.v2; }
	uint32_t getflags() { return e.node.flags; }
	uint64_t gete()  { return e.edge; }
	DepManGraphEdge(int32_t v1, uint32_t v2) { e.node.v1 = v1; e.node.v2 = v2; e.node.flags = 0; }
	DepManGraphEdge(uint64_t edge) { e.edge = edge; }
};

class _DependencyGraph {
public:
	ca_vec<int32_t> last;
	ca_vec<int32_t> adj;
	ca_vec<int32_t> next;
	int FirstAdjacent(int i) {
		return last[i];
	}
	int NextAdjacent(int j) {
		if (j > -1) return next[j];
		else return -1;
	}
	bool EdgeExists(int n1, int n2) {
	    for (int j=FirstAdjacent(n1);j>=0;j=NextAdjacent(j)) {
	    	int32_t content = adj[j];
	    	if (content == n2) return true;
	    }
	    return false;
	}
	void AddEdge(int i, int j) {
		if (EdgeExists(i,j)) return;
		adj.push(j);
		next.push(last[i]);
		last[i] = adj.size()-1;
		/*adj.push(i);
		next.push(last[j]);
		last[j] = adj.size()-1;*/
	}
	void DeleteEdge(int i, int j) {
		int searched_jj = -1;
		int last_jj = -1;
		int last_but_one_jj = -1;
    	for (int jj=FirstAdjacent(i);jj>=0;jj=NextAdjacent(jj)) {
    		last_but_one_jj = last_jj;
    		last_jj = jj;
    		if (adj[jj] == j) searched_jj = jj;
    	}
    	assert(searched_jj >=0 && last_jj >= 0);
    	adj[searched_jj] = adj[last_jj];
		if (last_but_one_jj >= 0) next[last_but_one_jj] = -1;
	}
	void ChangeEdge(int i, int j, int d) {
		// very special procedure; (i,d) is replaced by (i,j)
    	for (int jj=FirstAdjacent(i);jj>=0;jj=NextAdjacent(jj)) {
    		int32_t content = adj[jj];
    		if (content == d) {
    			if (!EdgeExists(i,j)) adj[jj] = j;
    			else DeleteEdge(i,d);
    			return;
    		}
    	}
    	assert(0);
	}
	void InitDependencyGraph(int n) {
		for (int i=0; i < n; i++)
			last.push(-1);
	}
};

class DependManager {
	static const int dmN  = 0;
	static const int dmUB = 1;
	static const int dmLB = 2;
	const ca_vec<int> &block;
	const ca_vec<int> &eas;
	ca_vec<int> depths;
	ca_vec<int32_t> comp;
	ca_vec<AVLCntContainer<uint64_t> > H_i;
	ca_vec<AVLCntContainer<uint64_t> > H_ip1;
	AVLCntContainer<uint64_t> mixedEdges;
	int mode;
	int nVars;
	_DependencyGraph DependencyGraph;
	_DependencyGraph AntiGraph;
	//AVLCntContainer<uint64_t> EdgeContainer;
	UnionFind Partition;
    // int Size() { return size; }
    // void Insert(T v, int q) { _Insert(first, v, q); }
    // void Remove(T v, int q) { _Remove(first,v,q); }
    // void Insert(T v) { _Insert(first,v,1); }
    // void Remove(T v) { _Remove(first,v,1); }
    // int IsContained(T v) { return _IsContained(first, v); }
    // void Clear() { _Clear(first); }
    // int IsEmpty() { if (first) return 0; else return 1; }
public:
	ca_vec<int32_t> fillrate;
	static int compareBlock(const void *a, const void *b) {
		std::pair<int, int>  *pa = (std::pair<int, int> *)a;
		std::pair<int, int>  *pb = (std::pair<int, int> *)b;
		if (pa->second < pb->second)
			return -1;
		if (pa->second == pb->second) {
			if (pa->first < pb->first)
				return -1;
			if (pa->first == pb->first)
				return 0;
			if (pa->first > pb->first)
				return 1;
		}
		if (pa->second > pb->second)
			return 1;
		return 0;
	}
	static int compareInt(const void *a, const void *b) {
		int  *pa = (int *)a;
		int  *pb = (int *)b;
		if (*pa < *pb)
			return -1;
		if (*pa == *pb)
		    return 0;
		if (*pa > *pb)
			return 1;
		return 0;
	}
	bool dependency(Var x, Var y) {
		if (block[y] <= block[x]) return false;
		return true;
		return YDependsOnX(x,y);
		//return EdgeIsInContainer(x,y);
		// original: return (block[x] < block[y]);
	}
	bool equiv(Var x, Var y) {
		if (block[x] == block[y]) return true;
		return false;
		if (block[x] < block[y] && YDependsOnX(x,y)) return false;
		if (block[x] > block[y] && YDependsOnX(y,x)) return false;
		//if (EdgeIsInContainer(x,y) || EdgeIsInContainer(y,x)) return false;
		return true;
		// original: return (block[x] == block[y]);
	}

	DependManager(ca_vec<int> &blo, ca_vec<int> &ea) : block(blo), eas(ea) {
		mode = dmN;
		nVars = 0;
	}
	void setMode(int m) { mode = m; }
	void setNVars(int n) { nVars = n; }
	void PartitionInit(int n) {
		Partition.UnionFindInit(n);
		for (int t=0;t<n;t++) {
			comp.push(-1);
			H_i.push();
			H_ip1.push();
		}
	}
	void rekSetComponents(int c, int v) {
		comp[v] = c;
    	for (int j=DependencyGraph.FirstAdjacent(v);j>=0;j=DependencyGraph.NextAdjacent(j)) {
    		int32_t content = DependencyGraph.adj[j];
     		if (comp[content] == -1) {
    			if (Partition.FindRepresenter(content) != Partition.FindRepresenter(v)) rekSetComponents(c, content);
    			else comp[content] = c;
    		}
    	}
    	for (int j=AntiGraph.FirstAdjacent(v);j>=0;j=AntiGraph.NextAdjacent(j)) {
    		int32_t content = AntiGraph.adj[j];
     		if (comp[content] == -1) {
    			if (Partition.FindRepresenter(content) != Partition.FindRepresenter(v)) rekSetComponents(c, content);
    			else comp[content] = c;
    		}
    	}
	}
	void setComponents() {
		for (int t = 0; t < nVars; t++) {
			if (comp[t] == -1) {
				cerr << "comp set " << t << ":" <<  Partition.FindRepresenter(t) << endl;
				if (Partition.FindRepresenter(t) != -1 && Partition.FindRepresenter(t) != t) {
					//assert(Partition.FindRepresenter(t) < t);
					comp[t] = comp[Partition.FindRepresenter(t)];
	    		} else rekSetComponents(t, t);
			}
		}
	}

	int findRoot(int noderep) {
		int j;
		while (AntiGraph.FirstAdjacent(noderep) >= 0) {
		   j = AntiGraph.FirstAdjacent(noderep);
		   assert(AntiGraph.NextAdjacent(j) < 0);
		   noderep = Partition.FindRepresenter(AntiGraph.adj[j]);
		}
		return noderep;
	}

	void mergePredTrees(int tree_node) {
    	int j = AntiGraph.FirstAdjacent(tree_node);
		cerr << "merge A " << j << " " << AntiGraph.adj[j] << endl;
    	if (j < 0) return;
		cerr << "merge B" << endl;
    	int pred1 = AntiGraph.adj[j];
    	j = AntiGraph.NextAdjacent(j);
        if (j < 0) return;
		cerr << "merge C" << endl;

    	int pred2 = AntiGraph.adj[j];
    	assert(pred1!=tree_node);
    	assert(pred2!=tree_node);
    	j = AntiGraph.NextAdjacent(j);
    	assert(j<0);
    	cerr << "in mergePredTrees mit " << pred1 << "," << pred2 << endl;

    	if (block[pred1] == block[pred2]) {
    		int node, node_to_be_deleted;
    		Partition.Union(pred1,pred2);
    		if (Partition.FindRepresenter(pred1) == pred1 && Partition.FindRepresenter(pred2) != pred2) {
    			cerr << " unify " << pred1 << " und " << pred2 << " zu " << pred1 << endl;
    			node = pred1;
    			node_to_be_deleted = pred2;
    		} else if (Partition.FindRepresenter(pred2) == pred2 && Partition.FindRepresenter(pred1) != pred1) {
    			cerr << " unify " << pred1 << " und " << pred2 << " zu " << pred2 << endl;
    			node = pred2;
    			node_to_be_deleted = pred1;
    		} else {
    			cerr << pred1 << " " << pred2 << " " << Partition.FindRepresenter(pred2)  << " " << Partition.FindRepresenter(pred2) << endl;
    			assert(0);
    		}
        	for (int j=DependencyGraph.FirstAdjacent(node_to_be_deleted);j>=0;j=DependencyGraph.NextAdjacent(j)) {
        		int32_t content = DependencyGraph.adj[j];
        		if (Partition.FindRepresenter(node)!=Partition.FindRepresenter(content))
        		   DependencyGraph.AddEdge(Partition.FindRepresenter(node),Partition.FindRepresenter(content));
    			cerr << Partition.FindRepresenter(node) << " ---> " << Partition.FindRepresenter(content) << endl;
        		if (Partition.FindRepresenter(node)!=Partition.FindRepresenter(content))
    			    AntiGraph.ChangeEdge(Partition.FindRepresenter(content),Partition.FindRepresenter(node), node_to_be_deleted);
    			cerr << Partition.FindRepresenter(node) << " <--- " << Partition.FindRepresenter(content) << endl;
        	}
        	for (int j=AntiGraph.FirstAdjacent(node_to_be_deleted);j>=0;j=AntiGraph.NextAdjacent(j)) {
        		int32_t content = AntiGraph.adj[j];
        		if (Partition.FindRepresenter(node)!=Partition.FindRepresenter(content))
        		    DependencyGraph.ChangeEdge(Partition.FindRepresenter(content),Partition.FindRepresenter(node), node_to_be_deleted);
        		cerr << Partition.FindRepresenter(content) << " ---> " << Partition.FindRepresenter(node) << endl;
        		if (Partition.FindRepresenter(node)!=Partition.FindRepresenter(content))
    			    AntiGraph.AddEdge(Partition.FindRepresenter(node),Partition.FindRepresenter(content));
    			cerr << Partition.FindRepresenter(content) << " <--- " << Partition.FindRepresenter(node) << endl;
        	}
        	mergePredTrees(node);
    	} else if (block[pred1] < block[pred2]) {
    		DependencyGraph.ChangeEdge(pred1,pred2,tree_node);
			cerr << Partition.FindRepresenter(pred1) << " ---> " << Partition.FindRepresenter(pred2) << endl;
			AntiGraph.DeleteEdge(tree_node,pred1);
    		AntiGraph.AddEdge(pred2, pred1);
			cerr << Partition.FindRepresenter(pred1) << " <--- " << Partition.FindRepresenter(pred2) << endl;
    		mergePredTrees(pred2);
    	} else {
    		DependencyGraph.ChangeEdge(pred2,pred1,tree_node);
			cerr << Partition.FindRepresenter(pred2) << " ---> " << Partition.FindRepresenter(pred1) << endl;
			AntiGraph.DeleteEdge(tree_node,pred2);
    		AntiGraph.AddEdge(pred1, pred2);
			cerr << Partition.FindRepresenter(pred2) << " <--- " << Partition.FindRepresenter(pred1) << endl;
    		mergePredTrees(pred1);
    	}
    }

	int FindSmallestAnchestor(int startnode, int b) {
		int j;
		int leader=startnode, second=-1;
		if (AntiGraph.FirstAdjacent(startnode) >= 0)
		    cerr << "|Anchestor von " << startnode << ": " << Partition.FindRepresenter(AntiGraph.adj[AntiGraph.FirstAdjacent(startnode)]);
		else
		    cerr << "|Anchestor von " << startnode << ": -1";
        if (AntiGraph.FirstAdjacent(startnode)<0) return startnode;
		while (AntiGraph.FirstAdjacent(startnode) >= 0) {
		   j = AntiGraph.FirstAdjacent(startnode);
		   assert(AntiGraph.NextAdjacent(j) < 0);
		   startnode = Partition.FindRepresenter(AntiGraph.adj[j]);
		   second = leader;
		   leader = startnode;
		   if (block[leader] < b) {
			   assert(second>-1);
			   cerr << "smallest Anch=" << second << endl;
			   return second;
		   }
		}
		cerr << "smallest Anch=" << leader << endl;
		return leader;
	}

	void addSubtree(int node, std::vector<int> &unionSet) {
    	bool answer=true;
        static ca_vec<int>      ana_stack;
        static ca_vec<int>      ana_seen_stack;
        static ca_vec<int8_t>   seen;
        static bool first = true;
        if (first) for(int z=0;z<nVars;z++) {
        	seen.push(0);
        	first = false;
        }
	    ana_stack.clear();
	    ana_seen_stack.clear();
	    ana_stack.push( node );
	    seen[node] = 1;
	    ana_seen_stack.push(node);

        while (ana_stack.size() > 0) {
        	int node = ana_stack.last();
	    	ana_stack.pop();
	    	for (int j=DependencyGraph.FirstAdjacent(node);j>=0;j=DependencyGraph.NextAdjacent(j)) {
	    		int32_t content = DependencyGraph.adj[j];
	    		if (seen[content]) continue;
	    		unionSet.push_back(content);
	    	    ana_stack.push( content );
	    	    seen[content] = 1;
	    	    ana_seen_stack.push(content);
	    	}
	    }
        while(ana_seen_stack.size() > 0) {
        	seen[ana_seen_stack.last()] = 0;
        	ana_seen_stack.pop();
        }
	}


	void scanConstraints( ca_vec<CRef>  &constraints, ConstraintAllocator  &constraintallocator) {
		std::vector<std::pair<int, int> > VarBlock;
		PartitionInit(nVars);
		return;
		cerr << "start scanning deps" << endl;
		for (int zz=0; zz < constraints.size(); zz++) {
			VarBlock.clear();
			Constraint &c = constraintallocator[constraints[zz]];
			cerr << "C"<<zz<<",("<<c.size()<<"):" << endl;
			for (int ii = 0; ii < c.size(); ii++) {
				std::pair<int,int> p(var(c[ii]), block[var(c[ii])]);
				VarBlock.push_back(p);
				cerr << p.first << "(" << p.second << ") | ";
			}
			cerr << endl;
			cerr << "start sorting" << endl;
			std::qsort((void*) VarBlock.data(), VarBlock.size(), sizeof(std::pair<int,int>), compareBlock);
			cerr << "stop sorting, start simplifying" << endl;
			int t = 0;
			for (t=0;t<VarBlock.size();t++)
				if (eas[VarBlock[t].first]==EXIST) {
					break;
				}
			for (;t<VarBlock.size();t++) {
				if (t+1<VarBlock.size() && eas[VarBlock[t+1].first]==EXIST && block[VarBlock[t+1].first] == block[VarBlock[t].first]) {
					Partition.Union(VarBlock[t].first,VarBlock[t+1].first);
					VarBlock[t+1].first = VarBlock[t].first;
				} else if (t==VarBlock.size()-1) {
					// fertig
				} else if (eas[VarBlock[t+1].first]==UNIV) {
					cerr << VarBlock[t].first << "E(" << VarBlock[t].second << ") | ";
					t = t+1;
					while (t < VarBlock.size() && eas[VarBlock[t].first]==UNIV) {
						cerr << VarBlock[t].first << "A(" << VarBlock[t].second << ") | ";
						t++;
					}
				}
				cerr << VarBlock[t].first << "E(" << VarBlock[t].second << ") | ";
			}
			cerr << "stop simplifying, start compression" << endl;
		    if (VarBlock.size()>1) { // compress double occuring Vars
			    int target=0, kompress_start=0, kompress_end=0;
		    	for (kompress_end=0 ;kompress_end < VarBlock.size();) {
		    		kompress_start=kompress_end;
			    	while (kompress_end < VarBlock.size() &&
			    		   VarBlock[kompress_start].first==VarBlock[kompress_end].first) {
			    			   kompress_end++;
			    	}
		    		VarBlock[target].second = VarBlock[kompress_start].second;
		    		VarBlock[target].first = Partition.FindRepresenter(VarBlock[kompress_start].first);
		    		target++;
		    	}
		    	/*int end = VarBlock.size()-target;
		    	for (int i=0; i < end;i++)
		    		VarBlock.pop_back();*/
		    	while (VarBlock.size() > target) VarBlock.pop_back();
		    }
			cerr << "finish compression" << endl;
			//creiere tupel var/block, sortiere, laufe von vorn nach hinten: wenn v < h add dep
			for (int jj=0; jj < VarBlock.size(); jj++) {
				if (eas[VarBlock[jj].first]==UNIV) {
					continue;
				}
				for (int ii = jj+1; ii < VarBlock.size(); ii++) {
					if (eas[VarBlock[ii].first]==UNIV) continue;
					while (ii < VarBlock.size() && Partition.FindRepresenter(VarBlock[ii].first) == Partition.FindRepresenter(VarBlock[jj].first)) {
                        ii++;
					}
					if (ii >= VarBlock.size()) break;
					if (Partition.FindRepresenter(VarBlock[ii].first) < Partition.FindRepresenter(VarBlock[jj].first)) {
						std::pair<int, int> tmp;
						tmp = VarBlock[ii];
						VarBlock[ii] = VarBlock[jj];
						VarBlock[jj] = tmp;
					}
					if (Partition.FindRepresenter(VarBlock[ii].first) <= Partition.FindRepresenter(VarBlock[jj].first)) {
						cerr << VarBlock[ii].first  << " " << VarBlock[jj].first  << " " << Partition.FindRepresenter(VarBlock[ii].first) << " " << Partition.FindRepresenter(VarBlock[jj].first) << endl;
						cerr << "D"<<zz<<",("<<c.size()<<"):" << endl;
						for (int ii = 0; ii < c.size(); ii++) {
							std::pair<int,int> p(var(c[ii]), block[var(c[ii])]);
							VarBlock.push_back(p);
							cerr << Partition.FindRepresenter(p.first) << "(" << p.second << ") | ";
						}
						cerr << endl;
						assert(0);//ii++;//break;
					}
					else if (Partition.FindRepresenter(VarBlock[ii].first) > Partition.FindRepresenter(VarBlock[jj].first)) {
						   if (DependencyGraph.EdgeExists(Partition.FindRepresenter(VarBlock[jj].first),Partition.FindRepresenter(VarBlock[ii].first))) break;
						   DependencyGraph.AddEdge(Partition.FindRepresenter(VarBlock[jj].first),Partition.FindRepresenter(VarBlock[ii].first));
						   cerr << VarBlock[jj].first << " ---> " << VarBlock[ii].first << endl;
						   AntiGraph.AddEdge(Partition.FindRepresenter(VarBlock[ii].first),Partition.FindRepresenter(VarBlock[jj].first));
						   cerr << VarBlock[jj].first << " <--- " << VarBlock[ii].first << endl;
						   cerr << "main loop merge predtrees von " << Partition.FindRepresenter(VarBlock[ii].first) << "," << VarBlock[ii].first << endl;
                           mergePredTrees(Partition.FindRepresenter(VarBlock[ii].first));
                   		for (int i = 0; i < nVars; i++) {
                   			if (Partition.FindRepresenter(i) != i) continue;
                   	    	for (int j=DependencyGraph.FirstAdjacent(i);j>=0;j=DependencyGraph.NextAdjacent(j)) {
                   	    		int32_t content = DependencyGraph.adj[j];
                   	     	    cerr << "Edge: " << i << " ----> " << content << endl;
                   	    	}
                   		}
                           break;
					}
				}
			}
		}
		setComponents();
		cerr << endl << "components: ";
		for (int t=0; t < nVars;t++)
			cerr << " " << comp[t];
		cerr << endl;
		for (int zz=0; zz < constraints.size(); zz++) {
			VarBlock.clear();
			Constraint &c = constraintallocator[constraints[zz]];
			cerr << "C"<<zz<<",("<<c.size()<<"):" << endl;
			for (int ii = 0; ii < c.size(); ii++) {
				std::pair<int,int> p(var(c[ii]), block[var(c[ii])]);
				VarBlock.push_back(p);
			}
			std::qsort((void*) VarBlock.data(), VarBlock.size(), sizeof(std::pair<int,int>), compareBlock);
			for (int t=0;t<VarBlock.size();t++) {
				VarBlock[t].first = Partition.FindRepresenter(VarBlock[t].first);
				VarBlock[t].second = block[Partition.FindRepresenter(VarBlock[t].first)];
			}
			cerr << "stop simplifying, start compression" << endl;
		    if (VarBlock.size()>1) { // compress double occuring Vars
			    int target=0, kompress_start=0, kompress_end=0;
		    	for (kompress_end=0 ;kompress_end < VarBlock.size();) {
		    		kompress_start=kompress_end;
			    	while (kompress_end < VarBlock.size() &&
			    		   VarBlock[kompress_start].first==VarBlock[kompress_end].first) {
			    			   kompress_end++;
			    	}
		    		VarBlock[target].second = VarBlock[kompress_start].second;
		    		VarBlock[target].first = Partition.FindRepresenter(VarBlock[kompress_start].first);
		    		target++;
		    	}
		    	while (VarBlock.size() > target) VarBlock.pop_back();
		    }
			cerr << "finish compression" << endl;
			cerr << "follow Varblock. size=" << VarBlock.size() << endl;
			for (int jj=0; jj < VarBlock.size(); jj++) {
				cerr << VarBlock[jj].first << "E(" << VarBlock[jj].second << ") | ";
				int ii = jj;
				for (; ii < VarBlock.size() && eas[VarBlock[ii].first]==UNIV; ii++)
					;
				cerr << endl << "Find smallest Anchestor of " << VarBlock[ii].first << " ";
				int x = FindSmallestAnchestor(Partition.FindRepresenter(VarBlock[ii].first), block[VarBlock[jj].first]);
				cerr << endl << "Find smallest Anchestor(2) of " << VarBlock[ii].first << " ";
				int y = FindSmallestAnchestor(Partition.FindRepresenter(VarBlock[ii].first), block[VarBlock[jj].first]+1);
				cerr << endl;
				cerr << x << "?isContauined?" << H_i[Partition.FindRepresenter(VarBlock[jj].first)].IsContained(x) << endl;
				cerr << y << "?isContauined?" << H_ip1[Partition.FindRepresenter(VarBlock[jj].first)].IsContained(y) << endl;
				cerr << ii << " A " << jj << " " <<  VarBlock[jj].first << " " << VarBlock[ii].first << endl;
				if (!H_i[Partition.FindRepresenter(VarBlock[jj].first)].IsContained(x))
					H_i[Partition.FindRepresenter(VarBlock[jj].first)].Insert(x);
				cerr << ii << " B " << jj << " " <<  VarBlock[jj].first << " " << VarBlock[ii].first << endl;
				if (!H_ip1[Partition.FindRepresenter(VarBlock[jj].first)].IsContained(y))
					H_ip1[Partition.FindRepresenter(VarBlock[jj].first)].Insert(y);
				cerr << ii << " C " << jj << " " <<  VarBlock[jj].first << " " << VarBlock[ii].first << endl;

			}
			cerr << "finish constraint" << endl;
		}
		for (int i = 0, li = -1; i < nVars; i++) {
			while(li>=0 && i < nVars && Partition.FindRepresenter(i)==Partition.FindRepresenter(li) ) i++;
			li=i;
			if (i >= nVars) break;
			for (int j = 0, lj = -1; j < nVars; j++) {
				while(lj>=0 && j < nVars && Partition.FindRepresenter(j)==Partition.FindRepresenter(lj) ) j++;
				lj=j;
				if (j >= nVars) break;
				if (eas[Partition.FindRepresenter(j)] != eas[Partition.FindRepresenter(i)]) {
					std::vector<int> unionSet;
					for (NodeType<uint64_t>* iter = H_i[j].begin();iter != H_i[j].end(); iter = H_i[j].next()) {
						int node = iter->data;
						unionSet.push_back(node);
						addSubtree(node, unionSet);
					}
					for (NodeType<uint64_t>* iter = H_ip1[i].begin();iter != H_ip1[i].end(); iter = H_ip1[i].next()) {
						int node = iter->data;
						unionSet.push_back(node);
						addSubtree(node, unionSet);
					}
					std::qsort((void*) unionSet.data(), unionSet.size(), sizeof(int), compareInt);
					if (unionSet.size() > 0)
			            for (int t=0; t < unionSet.size()-1;t++) {
			        	    if (unionSet[t] == unionSet[t+1]) {
			        		    cerr << "DEPENDENCY!!!  " << i << " " << j << " " << t << " " << unionSet.size() << endl;
			        		    if (block[i] < block[j]) {
			        		    	cerr << j << " on " << i << endl;
			        		    	if (!mixedEdges.IsContained(DepManGraphEdge(i,j).gete() ))
			        		    		mixedEdges.Insert(DepManGraphEdge(i,j).gete());
			        		    } else if (block[i] > block[j]) {
			        		    	cerr << i << " on " << j << endl;
			        		    	if (!mixedEdges.IsContained(DepManGraphEdge(j,i).gete() ))
			        		    		mixedEdges.Insert(DepManGraphEdge(j,i).gete());
			        		    } else cerr << i << "<=>" << j << endl;
			        	    }
			            }
					cerr << "durch" << endl;
				}
				cerr << "durch2" << endl;
			}
			cerr << "durch3" << endl;
		}
		cerr << "durch4" << endl;
		for (NodeType<uint64_t>* iter = mixedEdges.begin();iter != mixedEdges.end(); iter = mixedEdges.next()) {
			   DepManGraphEdge e = DepManGraphEdge(iter->data);
			   DependencyGraph.AddEdge(Partition.FindRepresenter( e.getv1() ),Partition.FindRepresenter( e.getv2() ));
			   cerr << e.getv1() << " ---> " << e.getv2() << endl;
			   AntiGraph.AddEdge(Partition.FindRepresenter( e.getv2() ),Partition.FindRepresenter( e.getv1() ));
			   cerr << e.getv1() << " <--- " << e.getv2() << endl;
               mergePredTrees(Partition.FindRepresenter( e.getv2() ));
		}
		setComponents();
		cerr << endl;
		for (int t=0; t < nVars;t++)
			cerr << " " << comp[t];
		cerr << endl;
		for (int ii = 0; ii < nVars;ii++) {
			if (Partition.FindRepresenter(ii) == ii) {
				cerr << "H_I(" << ii << ") : " << H_i[Partition.FindRepresenter(ii)].Size() << endl;
				cerr << "H_Ip1(" << ii << ") : " << H_ip1[Partition.FindRepresenter(ii)].Size() << endl;
			}
		}
		for (int i = 0; i < nVars; i++) {
			if (Partition.FindRepresenter(i) != i) continue;
	    	for (int j=DependencyGraph.FirstAdjacent(i);j>=0;j=DependencyGraph.NextAdjacent(j)) {
	    		int32_t content = DependencyGraph.adj[j];
	     	    cerr << "Edge: " << i << " ----> " << content << endl;
	    	}
		}
		cerr << endl << "components: ";
		setComponents();
		for (int t=0; t < nVars;t++)
			cerr << " " << comp[t];
		cerr << endl;
    	fillrate.capacity(nVars+2);
    	depths.capacity(nVars+2);
    	for (int i = 0; i < nVars; i++) {
    		fillrate[i] = 0;
    		depths[i] = determineDepth(i);
    	}
		cerr << endl << "end scanning deps" << endl;
	}

    void DepManagerInitGraph(int n) {
       nVars=n;
       DependencyGraph.InitDependencyGraph(n);
       AntiGraph.InitDependencyGraph(n);
       //GUBGraph.InitGUBGraph(n+n);
    }

    void initFillRate(int n, ca_vec<extbool> &assigns) {
    	fillrate.capacity(n+2);
    	for (int i = 0; i < n; i++) {
    		fillrate[i] = 0;
    	}
		for (int i = 0; i < nVars; i++)
			if (assigns[i] == extbool_Undef)
				fillrate[Partition.FindRepresenter(i)]++;
		//cerr << "Fillrates: ";
		//for (int i = 0; i < nVars; i++)
		//	cerr << " " << fillrate[i];
    }
    void increaseFillrate(int x) {
    	fillrate[Partition.FindRepresenter(x)]++;
    }
    void decreaseFillrate(int x) {
    	fillrate[Partition.FindRepresenter(x)]--;
    	//assert(fillrate[Partition.FindRepresenter(x)] >= 0);
    }
    int sumAllRates() {
    	int z=0;
    	for (int i = 0; i < nVars;i++)
    		z+=fillrate[i];
    	return z;
    }
    bool isValid(int x) {
    	int j;
        if (AntiGraph.FirstAdjacent(Partition.FindRepresenter(x))<0) {
        	return true;
        } else
        	return RisValid(x);
    }
    bool RisValid(int x) {
    	int j;
        if (AntiGraph.FirstAdjacent(Partition.FindRepresenter(x))<0) {
 		   if (fillrate[Partition.FindRepresenter(x)] == 0) return true;
 		   else return false;
        }
		if (AntiGraph.FirstAdjacent(Partition.FindRepresenter(x)) >= 0) {
		   j = AntiGraph.FirstAdjacent(Partition.FindRepresenter(x));
		   assert(AntiGraph.NextAdjacent(j) < 0);
		   if (fillrate[Partition.FindRepresenter(x)] == 0) return RisValid(AntiGraph.adj[j]);
		   else return false;
		}
		return false;
    }
    bool YDependsOnX(int x, int y) {
    	if (block[x] >= block[y]) return false;
    	bool answer=false;
        static ca_vec<int>      ana_stack;
        static ca_vec<int>      ana_seen_stack;
        static ca_vec<int8_t>   seen;
        static bool first = true;
        if (first) for(int z=0;z<nVars;z++) {
        	seen.push(0);
        	first = false;
        }
        x = Partition.FindRepresenter(x);
        y = Partition.FindRepresenter(y);
	    /*ana_stack.clear();
	    ana_seen_stack.clear();
    	for (int j=DependencyGraph.FirstAdjacent(x);j>=0;j=DependencyGraph.NextAdjacent(j)) {
    		int32_t content = DependencyGraph.adj[j];
    		if (seen[content]) continue;
    	    ana_stack.push( content );
    	    seen[content] = 1;
    	    ana_seen_stack.push(content);
    	}

        while (ana_stack.size() > 0) {
        	int k = ana_stack.last();
        	//cerr << "*" << k;
	    	ana_stack.pop();
	    	if (k == y) {
	    		answer=true;
	    		//cerr << "D(" << x << "," << y << ") ";
	    		while (ana_stack.size()>0) ana_stack.pop();
	    		break;
	    	}
	    	for (int j=DependencyGraph.FirstAdjacent(k);j>=0;j=DependencyGraph.NextAdjacent(j)) {
	    		int32_t content = DependencyGraph.adj[j];
	    		if (seen[content]) continue;
	    	    ana_stack.push( content );
	    	    seen[content] = 1;
	    	    ana_seen_stack.push(content);
	    	}
	    }
        while(ana_seen_stack.size() > 0) {
        	seen[ana_seen_stack.last()] = 0;
        	ana_seen_stack.pop();
        }*/
		int j;
		int node = y;
        if (AntiGraph.FirstAdjacent(node)<0) return false;
		while (AntiGraph.FirstAdjacent(node) >= 0) {
		   j = AntiGraph.FirstAdjacent(node);
		   assert(AntiGraph.NextAdjacent(j) < 0);
		   node = Partition.FindRepresenter(AntiGraph.adj[j]);
		   if (node == x) return true;
		}
    	if (0&&x==12&&y==29) {
    		if (answer==false) cerr << "N(" << x << "," << y << ") ";
    		else cerr << "D(" << x << "," << y << ") ";
    		//return true;
    	}
    	return answer;
    }
    int determineDepth(int x) {
        x = Partition.FindRepresenter(x);
		int j;
		int depth=0;
        if (AntiGraph.FirstAdjacent(x)<0) return depth;
		while (AntiGraph.FirstAdjacent(x) >= 0) {
		   j = AntiGraph.FirstAdjacent(x);
		   assert(AntiGraph.NextAdjacent(j) < 0);
		   x = Partition.FindRepresenter(AntiGraph.adj[j]);
		   depth++;
		}
    	return depth;
    }
};

#endif /* DEPENDENCIES_H_ */
