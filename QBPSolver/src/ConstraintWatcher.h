/*
*
* Yasol: ConstraintWatcher.h -- Copyright (c) 2012-2017 Ulf Lorenz
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


#ifndef CONSTRAINTWATCHER_H_
#define CONSTRAINTWATCHER_H_

#include "DataStructures.h"

#define TRACEVAR 0
#define TRACECON 392

static int cntt=0;

class ConstraintWatcher {
public:
	struct s_stack {
		int32_t lit;
		unsigned int occpt : 31;
		unsigned int active: 1;
	};
	int nVars;
	ca_vec<s_stack> stack;
    ca_vec<int32_t>  stack_lim;
    ca_vec<ca_vec<int32_t> > p_occlist;
    ca_vec<ca_vec<int32_t> > n_occlist;
    ca_vec<int32_t>  p_occpt;
    ca_vec<int32_t>  n_occpt;
    bool is_ready;

    bool isReady() { return is_ready; }

    void initConstraintWatcher(int n, ConstraintAllocator  &ca, ca_vec<CRef>  &constraints, int8_t *assigns, bool FeasOnly) {
    	nVars = n;
    	for (int i = 0; i < n; i++) {
    		p_occlist.push();
    		n_occlist.push();
    		p_occpt.push(-1);
    		n_occpt.push(-1);
    	}
    	stack_lim.push(0);
    	//printOccPts();
    	for (int i=constraints.size()-1; i >= 0;i--) {
    	//for (int i = 0; i<constraints.size(); i++) {
        	Constraint &c = ca[constraints[i]];
    		if (c.header.learnt != 0) {
    			cerr << "Error: Watched Constraint is learnt?" << endl;
    			continue;
    		}
        	c.header.watched = 1;
        	if (i>0 && c.saveFeas(assigns)) continue;
    		for (int j = 0; j < c.size();j++) {
    			//if (i == 828) cerr << "Var=" << var(c[j]) << endl;
    			if (sign(c[j])) {
    				if (n_occpt[var(c[j])] == -1) {
    		    		n_occpt[var(c[j])]=0;
    				}
    		        n_occlist[var(c[j])].push(i);
    		        //if (var(c[j]) == TRACEVAR) cerr << TRACEVAR << " ist neg in C" << i << endl;
    			} else {
    				if (p_occpt[var(c[j])] == -1) {
    		    		p_occpt[var(c[j])]=0;
    				}
    		        p_occlist[var(c[j])].push(i);
    		        //if (var(c[j]) == TRACEVAR) cerr << TRACEVAR << " ist pos in C" << i << endl;
    			}
    		}
    	}
    	is_ready = true;
    	//printOccPts();
    	//getchar();
    	//cerr << "getCwatcher of n" << TRACEVAR << ":" << getCWatcher(TRACEVAR*2+1) << endl;
    	//cerr << "getCwatcher of p" << TRACEVAR << ":" << getCWatcher(TRACEVAR*2) << endl;
    }
    int32_t getNextConstraint(CoeVar lit, ConstraintAllocator  &ca, ca_vec<CRef>  &constraints, int8_t *assigns) {
    	if (sign(lit)) {
    		do {
    			n_occpt[var(lit)]++;
    			if (n_occpt[var(lit)] >= n_occlist[var(lit)].size()) return -1;
    			int ci = n_occlist[var(lit)][n_occpt[var(lit)]];
    			Constraint &c=ca[constraints[ci]];
    			assert (!c.header.learnt);// break;
    			if (ci == 0 || !c.saveFeas(assigns)) return ci;
    		} while (1);
    	} else {
    		do {
    			p_occpt[var(lit)]++;
    			if (p_occpt[var(lit)] >= p_occlist[var(lit)].size()) return -1;
    			int ci = p_occlist[var(lit)][p_occpt[var(lit)]];
    			Constraint &c=ca[constraints[ci]];
    			assert (!c.header.learnt);// break;
    			if (ci == 0 || !c.saveFeas(assigns)) return ci;
    		} while (1);
    	}
    	return -1;
    }
    int getCWatcher(int lit) {
    	if (lit&1) {
    		if (n_occpt[lit >> 1] < 0) return -1;
			if (n_occpt[lit >> 1] >= n_occlist[lit>>1].size()) return -1;
    		return n_occlist[lit>>1][n_occpt[lit >> 1]];
    	} else {
    		if (p_occpt[lit >> 1] < 0) return -1;
			if (p_occpt[lit >> 1] >= p_occlist[lit>>1].size()) return -1;
    		return p_occlist[lit>>1][p_occpt[lit >> 1]];
    	}
    }

    void printOccPts() {
    	for (int z = 0; z < nVars;z++) cerr << p_occpt[z] << " ";
    	cerr << endl;
    	for (int z = 0; z < nVars;z++) cerr << n_occpt[z] << " ";
    	cerr << endl;
    }
    void assign(int variable, int signum, ConstraintAllocator  &ca, ca_vec<CRef>  &constraints, int8_t *assigns) {
    	//if (variable==TRACEVAR) cerr << "set " << TRACEVAR <<  " " << stack.size() << " " << stack_lim.size()<<endl;
    	cntt++;
    	int ci = getCWatcher(variable+variable+signum);
		if (ci==-1) {
			stack_lim.push(stack.size());
			return;
		}
		if (ci==0) {
			stack_lim.push(stack.size());
			return;
		}
    	while (ci != -1) {
			Constraint &c = ca[constraints[ci]];
			if (!c.saveFeas(assigns)) {
				stack_lim.push(stack.size());
				return;
			}

			//if (ci == TRACECON) cerr << TRACECON << " wird rausgenommen wegen var " << variable << "mit s=" << signum
			//		<<" und " << (int)assigns[variable] << endl;
			//if (cntt==454) cerr << "CONSTRAINT RAUS:" << ci << " weil x" << variable << "=" << 1-signum << endl;
			//if (cntt==454) c.print(c,assigns,false);
			for (int i = 0; i < c.size();i++) {
				int cii = getCWatcher(var(c[i])+var(c[i])+sign(c[i]));
				if (cii == ci) { // c ist im watcher von c[j] vorne
					//cerr << "ist vorn bei " << var(c[i])+var(c[i])+sign(c[i]) << " ->raus nit activity=" << c.header.watched << endl;
					//if (ci == TRACECON && variable == TRACEVAR) cerr << TRACECON << " wird wirklich rausgenommen" << endl;
					//assert (getCWatcher(var(c[i])+var(c[i])+sign(c[i])) != 872);
					s_stack E;
					E.lit = var(c[i])+var(c[i]) + sign(c[i]);
					E.active = (c.header.watched&1);
					if (0&&ci == TRACECON) {
						if (sign(c[i])) cerr << "n_occpt wird gesaved: v=" << var(c[i]) << endl;
						else cerr << "p_occpt wird gesaved: v="  << var(c[i]) << endl;
					}
					E.occpt = (sign(c[i]) ? n_occpt[var(c[i])] : p_occpt[var(c[i])]);
					getNextConstraint(c[i], ca, constraints, assigns);
					stack.push(E);
				}
			}
			//assert(c.header.watched==1);
			c.header.watched = 0;
			//break;
			ci = getCWatcher(variable+variable+signum);
    	}
		stack_lim.push(stack.size());
    }
    void unassign(int variable, ConstraintAllocator  &ca, ca_vec<CRef>  &constraints) {
    	stack_lim.pop();
    	if (0&&variable==TRACEVAR) {
    		cerr << "unset " << TRACEVAR << " " << stack.size() << " " << stack_lim.size() <<endl << "stacklim:";
    		for (int z=0; z < stack_lim.size();z++) cerr << stack_lim[z] << " " ;
    		cerr << endl;
    	}
    	int targetstp = stack_lim.last();
    	while (stack.size() > targetstp/*stack_lim[stack_lim.size()-1]*/) {
    		int a;
    		s_stack E;
    		E = stack.last();
    		stack.pop();
    		if (stack.size()==0) assert(targetstp==0);
    		if (E.lit & 1) { // hat vorzeichen
    			n_occpt[E.lit >> 1] = E.occpt;
    			a = E.active;
    			Constraint &c = ca[constraints[n_occlist[E.lit>>1][n_occpt[E.lit>>1]]]];
    			////assert(n_occlist[E.lit>>1][n_occpt[E.lit>>1]] != 0);
    			//if (cntt==454) cerr << "REIN:" << n_occlist[E.lit>>1][n_occpt[E.lit>>1]] << endl;
    			//if (cntt==454) cerr << " weil x" << (E.lit>>1) << "=" << 1-(E.lit&1) << " mit activity=" << a << endl;
    			c.header.watched = (E.active&1);
    			//if ((E.lit >> 1) == TRACEVAR) cerr << "back: v" << TRACEVAR << " mit c=" << n_occlist[E.lit >> 1][n_occpt[E.lit >> 1]] << endl;
    		} else {
    			p_occpt[E.lit >> 1] = E.occpt;
    			a = E.active;
    			Constraint &c = ca[constraints[p_occlist[E.lit>>1][p_occpt[E.lit>>1]]]];
    			////assert(p_occlist[E.lit>>1][p_occpt[E.lit>>1]] != 0);
    			//if (cntt==454) cerr << "REIN:" << p_occlist[E.lit>>1][p_occpt[E.lit>>1]] << endl;
    			//if (cntt==454) cerr << " weil x" << (E.lit>>1) << "=" << 1-(E.lit&1) << " mit activity=" << a << endl;
    			//cerr << "LIT" << E.lit << "geht rin mit activity=" << E.active << endl;
    			c.header.watched = (E.active&1);
    			//if ((E.lit >> 1) == TRACEVAR) cerr << "back: v" << TRACEVAR << " mit c=" << n_occlist[E.lit >> 1][n_occpt[E.lit >> 1]] << endl;
    		}
    	}
    	//if (variable==TRACEVAR) cerr << "unset stacksize=" << stack.size() << " mit " << stack_lim.size() << endl;
    	//if (cntt==454) cerr << "U ";
    	//if (cntt==454) printOccPts();
    }
    ConstraintWatcher() { is_ready = false; nVars = 0; }
};


#endif /* CONSTRAINTWATCHER_H_ */
