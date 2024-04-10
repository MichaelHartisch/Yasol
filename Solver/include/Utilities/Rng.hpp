/*
*
* Solver: Rng.hpp -- Copyright (c) 2010-2017 Jan Wolf
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

#ifndef RNG_HPP_
#define RNG_HPP_
/*
 Random number generator class
 =============================
 History:

 Created - Sarah "Voodoo Doll" White (2006/01/24)
 =============================
 Description:

 This class wraps the Mersenne Twister generator
 with a public interface that supports three common
 pseudorandom number requests:

 === Uniform deviate [0,1) ===
 Random rnd(seed);
 double r = rnd.uniform();

 === Uniform deviate [0,hi) ===
 Random rnd(seed);
 unsigned long r = rnd.uniform(hi);

 === Uniform deviate [lo,hi) ===
 Random rnd(seed);
 unsigned long r = rnd.uniform(lo, hi);

 seed, lo, and hi are user supplied values, with
 seed having a default setting of 1 for debugging
 and testing purposes.
 */
namespace utils {
class Rng {
	// Arbitrary constants that work well
	static const int N = 624;
	static const int M = 397;
	static const unsigned long MATRIX_A = 0x9908b0dfUL;
	static const unsigned long UPPER_MASK = 0x80000000UL;
	static const unsigned long LOWER_MASK = 0x7fffffffUL;
	static const unsigned long MAX = 0xffffffffUL;

public:
	Rng(unsigned long seed=1):next(0) {
		seedgen(seed);
	}

	// Return a uniform deviate in the range [0,1)
	double uniform() {
		return randgen() * (1.0 / (MAX + 1.0));
	}

	// Return a uniform deviate in the range [0,hi)
	int uniform(int hi){
		return (uniform() * hi);
	}

	// Return a uniform deviate in the range [lo,hi)
	int uniform(int lo, int hi) {
		int val;
		while ((val = lo + uniform(hi - lo + 1)) > hi || val < lo)
			;
		return val;
	}

private:

	void seedgen(unsigned long seed){
		x[0] = seed & MAX;
		for (int i = 1; i < N; i++) {
			x[i] = (1812433253UL * (x[i - 1] ^ (x[i - 1] >> 30)) + i);
			x[i] &= MAX;
		}
	}

	unsigned long randgen(){
		unsigned long rnd;
		// Refill the pool when exhausted
		if (next == N) {
			int a;
			for (int i = 0; i < N - 1; i++) {
				rnd = ((x[i] & UPPER_MASK) | (x[i + 1] & LOWER_MASK));
				a = (rnd & 0x1UL) ? MATRIX_A : 0x0UL;
				x[i] = x[(i + M) % N] ^ (rnd >> 1) ^ a;
			}

			rnd = ((x[N - 1] & UPPER_MASK) | (x[0] & LOWER_MASK));
			a = (rnd & 0x1UL) ? MATRIX_A : 0x0UL;
			x[N - 1] = x[M - 1] ^ (rnd >> 1) ^ a;
			next = 0; // Rewind index
		}

		rnd = x[next++]; // Grab the next number
		// Voodoo to improve distribution
		rnd ^= (rnd >> 11);
		rnd ^= (rnd << 7) & 0x9d2c5680UL;
		rnd ^= (rnd << 15) & 0xefc60000UL;
		rnd ^= (rnd >> 18);

		return rnd;
	}

	unsigned long x[N]; // Random number pool
	int next; // Current pool index
};
}
#endif /*RNG_HPP_*/
