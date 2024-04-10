/*
*
* Yasol: UnionFind.h -- Copyright (c) 2012-2017 Ulf Lorenz
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

#ifndef UNIONFIND_H_
#define UNIONFIND_H_

// textbook approach

class UnionFind {
	int N;
	std::vector<int> predecessor;
	std::vector<int> rank;
public:
	void UnionFindInit(int n) {
		N = n;
		for (int i=0; i < n;i++) {
			predecessor.push_back(-1);
			rank.push_back(1);
		}
	}
	int FindRepresenter(int x) {
		if (predecessor[x] == -1) return x;
		else {
			predecessor[x] = FindRepresenter(predecessor[x]);
			return predecessor[x];
		}
	}
	void UnionRepresenters(int x, int y) {
		if (x==y) return;
		if (rank[x] == rank[y]) rank[x]++;
		if (/*rank[x] < rank[y]*/x > y) {
			predecessor[x] = y;
		} else {
			predecessor[y] = x;
		}
	}
	void Union(int x, int y) {
		UnionRepresenters(FindRepresenter(x), FindRepresenter(y));
	}
};


#endif /* UNIONFIND_H_ */
