/*
*
* Yasol: HashTable.h -- Copyright (c) 2012-2017 Ulf Lorenz
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

#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include "DataStructures.h"
#define EXIST 0
#define UNIV  1
#define LB   1
#define UB   2
#define FIT  3
#define CONSTRAINT 4
#define ASSIGN_OK -1

using namespace std;
class HTentry {
public:
	uint64_t hash;
	coef_t value;
	int polvar;
	int32_t weight;
	int32_t lengthOfTrail;
	int32_t objective_iterations;
	int8_t ea;
	int8_t bound;

	int getVar() { return (polvar>>1); }
	void setPolVar(int v, int p) { polvar = v+v+p; }
	int getPol() { return (polvar & 1); }
	void setEntry(int64_t h,coef_t val,int pv,int32_t w, int8_t lea, int8_t b, int32_t l, int32_t oi) {
		hash = h;
		value = val;
		polvar = pv;
		weight = w;
		ea = lea;
		bound = b;
		lengthOfTrail = l;
		objective_iterations = oi;
	}
};

class HTable {
	double seed;
	HTentry* table;
	uint64_t *(hashconsts[2]);
	uint64_t hash;
	int size;
public:
	int getSize() { return size; }
	uint64_t getHashConstant(int lit) {
		return hashconsts[lit&1][lit>>1];
	}
	HTable( int num_vars, int s) {
        size = s;
		seed = 25.0;
		hash = 0;
		table = (HTentry*)realloc(NULL/*table*/,(size+4)*sizeof(HTentry));
		hashconsts[0] = (uint64_t*)realloc(NULL/*hashconsts[0]*/ ,(num_vars+2)*sizeof(uint64_t));
		hashconsts[1] = (uint64_t*)realloc(NULL/*hashconsts[1]*/ ,(num_vars+2)*sizeof(uint64_t));
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < num_vars; j++)
				hashconsts[i][j] = irand(seed,0xffffffffffffffff);

		for (int i = 0; i < size; i++) {
			table[i].hash = 0;
			table[i].bound = 0;
		}
	}

	~HTable() {
	        if (table) free(table);
		for (int i = 0; i < 2; i++)
		  if (hashconsts[i]) free(hashconsts[i]);
		hashconsts[0] = hashconsts[1] = 0;
		table = 0;
	}

	inline void assign(int var, int pol) {
		if (pol == extbool_Undef) {
		  cerr << "Passiert das?" << endl;
		}
		hash ^= hashconsts[pol][var];
	}

	inline void unassign(int var, int pol) {
		if (pol == extbool_Undef) {
		  cerr << "Warning: HT->unassign with problem." << endl;
		  return;
		}
		hash ^= hashconsts[pol][var];
	}

	bool getEntry(HTentry **hte, int32_t l) {
		int64_t index = (hash / 100000) % size;
		if (table[index].hash == 0) return false;
		if (table[index].hash == hash && table[index].lengthOfTrail == l) {
			*hte = &(table[index]);
			cerr << "H";
			return true;
		}
		return false;
	}

	bool setEntry(coef_t val,int pol, int var,int32_t w, int8_t lea, int8_t b, int32_t l, int32_t oi, coef_t dont_know, bool break_from_outside) {
        if (break_from_outside) return false;
		int64_t index = (hash / 100000) % size;
		if (!((b == FIT || b == LB) && val > dont_know))
			if (!((b == FIT || b == UB) && val < dont_know)) return false;
		if (table[index].bound == CONSTRAINT) return false;
		if (table[index].hash == hash && table[index].lengthOfTrail == l) {
		  if (/*table[index].weight < w ||*/ table[index].hash == 0) {
			  table[index].setEntry(hash, val,var+var+pol,w, lea, b, l, oi);
			  return true;
		  }
		} else if (table[index].lengthOfTrail >= l || table[index].hash == 0) {
			  table[index].setEntry(hash, val,var+var+pol,w, lea, b, l, oi);
			  return true;
		}
		//cerr << "M";
		return false;
	}

    bool SatConstraintExists(ca_vec<CoeVar> &ps) {
    	uint64_t hash=0;
    	for (int i=0; i < ps.size(); i++)
    		hash ^= hashconsts[sign(ps[i])][var(ps[i])];
		int64_t index = hash % size;
		if (table[index].bound == CONSTRAINT) {
			if (table[index].hash == hash && table[index].lengthOfTrail == ps.size()) {
	        	cerr << "SAT-Constraint exists!!!" << endl;
				return true;
			}
		}
    	return false;
    }

    bool SatConstraintAdd(ca_vec<CoeVar> &ps) {
    	uint64_t hash=0;
    	for (int i=0; i < ps.size(); i++)
    		hash ^= hashconsts[sign(ps[i])][var(ps[i])];
		int64_t index = hash % size;
		if (table[index].bound != CONSTRAINT) {
			table[index].hash = hash;
			table[index].lengthOfTrail = ps.size();
			return true;
		}
    	return false;
    }

    bool SatConstraintDelete(ca_vec<CoeVar> &ps) {
    	uint64_t hash=0;
    	for (int i=0; i < ps.size(); i++)
    		hash ^= hashconsts[sign(ps[i])][var(ps[i])];
		int64_t index = hash % size;
		if (table[index].bound == CONSTRAINT && table[index].hash == hash && table[index].lengthOfTrail == size) {
			table[index].hash = 0;
			table[index].lengthOfTrail = 0;
			return true;
		}
    	return false;
    }

    bool SatConstraintExists(CoeVar* data, int s) {
    	uint64_t hash=0;
    	for (int i=0; i < s; i++)
    		hash ^= hashconsts[sign(data[i])][var(data[i])];
		int64_t index = hash % size;
		if (table[index].bound == CONSTRAINT) {
			if (table[index].hash == hash && table[index].lengthOfTrail == s) return true;
		}
    	return false;
    }

    bool SatConstraintAdd(CoeVar* data, int s) {
    	uint64_t hash=0;
    	for (int i=0; i < s; i++)
    		hash ^= hashconsts[sign(data[i])][var(data[i])];
		int64_t index = hash % size;
		if (table[index].bound != CONSTRAINT) {
			table[index].hash = hash;
			table[index].lengthOfTrail = s;
			return true;
		}
    	return false;
    }

    bool SatConstraintDelete(CoeVar* data, int s) {
    	uint64_t hash=0;
    	for (int i=0; i < s; i++)
    		hash ^= hashconsts[sign(data[i])][var(data[i])];
		int64_t index = hash % size;
		if (table[index].bound == CONSTRAINT && table[index].hash == hash && table[index].lengthOfTrail == s) {
			table[index].hash = 0;
			table[index].lengthOfTrail = 0;
			return true;
		}
    	return false;
    }

    // Returns a random float 0 <= x < 1. Seed must never be 0.
    static inline double drand(double& seed) {
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;
        return seed / 2147483647; }

    // Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline uint64_t irand(double& seed, uint64_t size) {
        return (uint64_t)(drand(seed) * size); }

};

// ---- special hashtable for LP-cuts

class HTCutentry {
public:
	uint64_t hash;
	coef_t coefsum;
	int index;

	void setEntry(int64_t h,coef_t cs, int ind=0) {
		hash = h;
		coefsum = cs;
		index = ind;
	}
};

class HCTable {
	double seed;
	HTCutentry* table;
	uint64_t *(hashconsts[2]);
	uint64_t hash;
	int size;
        int num_vars;
public:
	uint64_t getHashConstant(int lit) {
		return hashconsts[lit&1][lit>>1];
	}
	int getSize() { return size; }
	HCTable( int nv, int s) {
                size = s;
                num_vars = nv;
		seed = (double)time(NULL);
		hash = 0;
		table = (HTCutentry*)realloc(NULL,(size+4)*sizeof(HTCutentry));
		hashconsts[0] = (uint64_t*)realloc(NULL,(num_vars+2)*sizeof(uint64_t));
		hashconsts[1] = (uint64_t*)realloc(NULL,(num_vars+2)*sizeof(uint64_t));
		for (int i = 0; i < 1; i++)
			for (int j = 0; j < num_vars; j++) {
				hashconsts[i][j] = irand(seed,0xffffffffffffffff);
				hashconsts[i+1][j] = hashconsts[i][j] ^ 0xffffffffffffffff;
				//cerr << "H" << i << j << "=" << ((void*)hashconsts[i][j]) << endl;
				//cerr << "H" << i+1 << j << "=" << ((void*)hashconsts[i+1][j]) << endl;
			}
		for (int i = 0; i < size; i++) {
			table[i].hash = 0;
		}
	}

	~HCTable() {
		free(table);
		for (int i = 0; i < 2; i++)
			free(hashconsts[i]);
	}

	void clear(int nv=0) {
	        if (nv > num_vars) {
		   num_vars = nv;
		   hashconsts[0] = (uint64_t*)realloc(hashconsts[0]/*NULL*/,(num_vars+2)*sizeof(uint64_t));
		   hashconsts[1] = (uint64_t*)realloc(hashconsts[1]/*NULL*/,(num_vars+2)*sizeof(uint64_t));
		   for (int i = 0; i < 1; i++)
			for (int j = 0; j < num_vars; j++) {
				hashconsts[i][j] = irand(seed,0xffffffffffffffff);
				hashconsts[i+1][j] = hashconsts[i][j] ^ 0xffffffffffffffff;
				//cerr << "H" << i << j << "=" << ((void*)hashconsts[i][j]) << endl;
				//cerr << "H" << i+1 << j << "=" << ((void*)hashconsts[i+1][j]) << endl;
			}
		}                
		for (int i = 0; i < size; i++) {
			table[i].hash = 0;
		}
	}

    std::pair<coef_t,uint64_t> computeHash(std::vector<data::IndexedElement> &lhs, double rhs) {
    	uint64_t h=0;
    	coef_t csum=0.0;
        static const uint64_t geqH = 0x99647;
   	    for (int i = 0; i < lhs.size();i++) {
    		if (lhs[i].value.asDouble() < 0.0) h ^= hashconsts[1][lhs[i].index];
    		else                               h ^= hashconsts[0][lhs[i].index];
    		csum += abs(lhs[i].value.asDouble());
    	}
    	csum += abs(rhs);
        h ^= geqH;
    	return std::pair<coef_t,uint64_t>(/*csum*/rhs,h);
    }

    std::pair<coef_t,uint64_t> computeHash(std::vector<data::IndexedElement> &lhs, double rhs, data::QpRhs::RatioSign rs) {
    	uint64_t h=0;
    	//coef_t csum=0.0;
    	static const uint64_t leqH = 0xfa99647;
    	static const uint64_t geqH = 0x99647;
    	static const uint64_t eqH  = 0xffffa99647;
    	for (int i = 0; i < lhs.size();i++) {
    		double coef = lhs[i].value.asDouble();
    		if (rs == 2) coef = -coef;
		//std::cerr << "S" << lhs[i].index << "," << num_vars << std::endl;
    		if (coef < 0.0) h ^= hashconsts[1][lhs[i].index];
    		else            h ^= hashconsts[0][lhs[i].index];
    		//csum += abs(lhs[i].value.asDouble());
    	}
    	if (rs == 2) {
    		rhs = -rhs;
    		rs = data::QpRhs::smallerThanOrEqual;
    	}
		if (rs == 1/*data::QpRhs::RatioSign::equal*/) {
			h ^= eqH;
		} else if (rs == 0/*data::QpRhs::RatioSign::smallerThanOrEqual*/) {
			h ^= leqH;
		} if (rs == 2) {
			h ^= geqH;
		}
    	return std::pair<coef_t,uint64_t>(/*csum*/rhs,h);
    }

    /*std::pair<coef_t,uint64_t> computeHash(std::vector<coef_t> &lhs_c, std::vector<unsigned int> &lhs_x, coef_t rhs) {
    	uint64_t h=0;
    	coef_t csum=0.0;
    	for (int i = 0; i < lhs_c.size();i++) {
    		if (lhs_c[i] < 0.0) h ^= hashconsts[1][lhs_x[i]];
    		else                h ^= hashconsts[0][lhs_x[i]];
    		csum += abs(lhs_c[i]);
    	}
    	csum += abs(rhs);
    	return std::pair<coef_t,uint64_t>(csum,h);
    }*/

	bool getEntry(HTCutentry **hte, uint64_t &hash, coef_t &csum) {
		//return false;
		int64_t index = hash % size;
		if (table[index].hash == 0) return false;
		if (table[index].hash == hash && abs(csum-table[index].coefsum) < 0.0001) {
			*hte = &(table[index]);
			//cerr << "H";// << ((void*)hash);
			return true;
		}
		return false;
	}

	bool getEntry(HTCutentry **hte, uint64_t &hash, coef_t &csum, data::QpRhs::RatioSign rs) {
		//return false;
		int64_t index = hash % size;
		if (table[index].hash == 0) return false;
		if (table[index].hash != hash) index = (index + 1) % size;
		if (table[index].hash == 0) return false;
		if (table[index].hash != hash) index = (index + 1) % size;
		if (table[index].hash == 0) return false;
		if (table[index].hash == hash && rs == 1/*data::QpRhs::RatioSign::equal*/ && abs(csum-table[index].coefsum) < 0.0001) {
			*hte = &(table[index]);
			//cerr << "H1";// << ((void*)hash);
			return true;
		} else if (table[index].hash == hash && rs == 0/*data::QpRhs::RatioSign::smallerThanOrEqual*/ && csum>=table[index].coefsum-0.0001) {
			*hte = &(table[index]);
			//cerr << "H2";// << ((void*)hash);
			return true;
		} if (table[index].hash == hash && rs == 2 && csum<=table[index].coefsum+0.0001) {
			*hte = &(table[index]);
			//cerr << "H3";// << ((void*)hash);
			return true;
		}
		return false;
	}

	bool setEntry(coef_t &csum, uint64_t &hash, int s=0) {
		int64_t index = hash % size;
		if (table[index].hash != 0) index = (index + 1) % size;
		if (table[index].hash != 0) index = (index + 1) % size;
		if (table[index].hash != 0) return false;
		table[index].setEntry(hash, csum, s);
		return false;
	}

	void delEntry(coef_t &csum, uint64_t &hash) {
		int64_t index = hash % size;
		if (table[index].hash == hash && abs(csum-table[index].coefsum) < 0.0001)
			table[index].setEntry(0, 0);
	}

	// Returns a random float 0 <= x < 1. Seed must never be 0.
    static inline double drand(double& seed) {
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;
        return seed / 2147483647; }

    // Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline uint64_t irand(double& seed, uint64_t size) {
        return (uint64_t)(drand(seed) * size); }
};

// ---- special hashtable for cliques

class Cli_entry {
public:
	uint64_t hash;
	CoeVar cova;

	void setEntry(int64_t h,CoeVar &cv) {
		hash = h;
		cova = cv;
	}
};

class Cli_Table {
	double seed;
	Cli_entry* table;
	uint64_t *(hashconsts[2]);
	bool hash_consts_exist = false;
	uint64_t hash;
	int size;
	CoeVar nonCoeVar;
public:
	uint64_t getHashConstant(int lit) {
		return hashconsts[lit&1][lit>>1];
	}
	Cli_Table( int num_vars, int s) {
        size = s;
		seed = (double)time(NULL);
		hash = 0;
		table = (Cli_entry*)realloc(NULL,(size+4)*sizeof(Cli_entry));
		nonCoeVar = mkCoeVar(-1,0.0,false);
		if (hash_consts_exist == false) {
			hash_consts_exist = true;
			hashconsts[0] = (uint64_t*)realloc(NULL,(num_vars+2)*sizeof(uint64_t));
			hashconsts[1] = (uint64_t*)realloc(NULL,(num_vars+2)*sizeof(uint64_t));
			for (int i = 0; i < 1; i++)
				for (int j = 0; j < num_vars; j++) {
					hashconsts[i][j] = irand(seed,0xffffffffffffffff);
					hashconsts[i+1][j] = hashconsts[i][j] ^ 0xffffffffffffffff;
					//cerr << "H" << i << j << "=" << ((void*)hashconsts[i][j]) << endl;
					//cerr << "H" << i+1 << j << "=" << ((void*)hashconsts[i+1][j]) << endl;
				}
		}
		for (int i = 0; i < size; i++) {
			table[i].hash = 0;
		}
	}

	~Cli_Table() {
		free(table);
		for (int i = 0; i < 2; i++) {
			if (hashconsts[i] != NULL) {
				free(hashconsts[i]);
				hashconsts[i] = NULL;
			}
		}
	}

    uint64_t computeHash(CoeVar &cv) {
    	uint64_t h=hashconsts[1][var(cv)];
    	return h;
    }

	bool isInClique(Cli_entry **hte, uint64_t &hash, int v) {
		//return false;
		int index = hash % size;
		if (table[index].hash == 0) return false;
		for (int i = 0; i < 5;i++)
			if (table[(index+i)%size].hash == hash) {
				index = (index+i)%size;
				break;
			}
		if (table[index].hash == hash && v == var(table[index].cova)) {
			*hte = &(table[index]);
			//cerr << "H";// << ((void*)hash);
			return true;
		}
		return false;
	}

	bool setEntry(uint64_t &hash, CoeVar &cv) {
		int index = hash % size;
		for (int i = 0; i < 5; i++) {
			if (table[(index+i)%size].hash == 0) {
				table[(index+i)%size].setEntry(hash, cv);
				return true;
			}
		}
		return false;
	}

	/*void delEntry(uint64_t &hash) {
		int64_t index = hash % size;
		if (table[index].hash == hash)
			table[index].setEntry(0, nonCoeVar);
		else if (table[index].hash != 0) {
			for (int i = 1; i < 5;i++)
				if (table[(index+i)%size].hash == 0 || table[(index+i)%size].hash == hash) {
					index = (index+i)%size;
					break;
				}
			if (table[index].hash == hash)
				table[index].setEntry(0, nonCoeVar);
		}
	}*/

	// Returns a random float 0 <= x < 1. Seed must never be 0.
    static inline double drand(double& seed) {
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;
        return seed / 2147483647; }

    // Returns a random integer 0 <= x < size. Seed must never be 0.
    static inline uint64_t irand(double& seed, uint64_t size) {
        return (uint64_t)(drand(seed) * size); }
};

#endif /* HASHTABLE_H_ */
