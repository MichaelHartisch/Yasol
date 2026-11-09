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

#include <iostream>
using namespace std;

#define _AVL_MAIN
#include "cliques.h"

void CliqueManager::buildCliqueSepGraph(
    const std::vector<double>& xlp,
    const int* types,
    CliqueSepGraph& graph,
    double eps)
{
    graph.clear();

    if (types == nullptr)
        return;

    if ((int)xlp.size() < nVars)
        return;

    // ------------------------------------------------------------
    // Mapping from original CM literal numbers to local nodes.
    //
    // CM literal encoding:
    //
    //     2*i     : x_i = 1
    //     2*i + 1 : x_i = 0
    // ------------------------------------------------------------

    graph.localOfLiteral.assign(2 * nVars, -1);

    graph.literal.reserve(2 * nVars);
    graph.weight.reserve(2 * nVars);

    auto addLiteral = [&](int lit, double weight) {
        if (weight <= eps) {
            graph.zeroLiteral.push_back(lit);
            return;
        }
        const int local =
            (int)graph.literal.size();
        graph.localOfLiteral[lit] = local;
        graph.literal.push_back(lit);
        graph.weight.push_back(weight);
    };

    // ------------------------------------------------------------
    // First pass:
    // create local nodes for binary literals having positive
    // weight in the current LP solution.
    // ------------------------------------------------------------

    for (int i = 0; i < nVars; ++i) {

        // Yasol: types[i] == 0 means binary.
        if (types[i] != 0)
            continue;

        double x = xlp[i];

        // LP solutions may be microscopically outside [0,1].
        if (x < 0.0)
            x = 0.0;
        else if (x > 1.0)
            x = 1.0;

        addLiteral(2 * i,     x);
        addLiteral(2 * i + 1, 1.0 - x);
    }

    // ------------------------------------------------------------
    // Second pass:
    // copy the relevant part of the Yasol conflict graph.
    //
    // IMPORTANT:
    // Access ConflictGraph directly here. Do NOT use
    // EdgeIsInContainer(), FirstAdjacentInConflictGraph(), ...
    // inside the later clique search.
    // ------------------------------------------------------------

    graph.adj.resize(graph.literal.size());

    for (int u = 0;
         u < (int)graph.literal.size();
         ++u) {

        const int lit =
            graph.literal[u];

        for (int e = ConflictGraph.last[lit];
             e >= 0;
             e = ConflictGraph.next[e]) {

            const int otherLit =
                ConflictGraph.adj[e];

            if (otherLit < 0 ||
                otherLit >= (int)graph.localOfLiteral.size())
	      continue;

            const int v =
	      graph.localOfLiteral[otherLit];

            // The other endpoint was not selected for the local graph.
            if (v < 0)
	      continue;

            if (v == u)
	      continue;

            graph.adj[u].push_back(v);
        }
    }
    // ------------------------------------------------------------
    // Add the trivial conflict between the two literals of the
    // same binary variable:
    //
    //     x_i = 0   and   x_i = 1
    //
    // Yasol's ConflictGraph does not need to store these edges
    // explicitly, but for clique separation they belong to the
    // local conflict graph.
    // ------------------------------------------------------------

    int complementPairs   = 0;
    int complementAlready = 0;
    int complementAdded   = 0;

    for (int i = 0; i < nVars; ++i) {

      const int lit1 = 2 * i;
      const int lit0 = 2 * i + 1;

      const int u = graph.localOfLiteral[lit1];
      const int v = graph.localOfLiteral[lit0];

      if (u < 0 || v < 0)
        continue;

      ++complementPairs;

      bool alreadyThere = false;

      for (int w : graph.adj[u]) {
        if (w == v) {
	  alreadyThere = true;
	  break;
        }
      }

      if (alreadyThere) {

        ++complementAlready;

      } else {

        graph.adj[u].push_back(v);
        graph.adj[v].push_back(u);

        ++complementAdded;
      }
    }
    for (auto& a : graph.adj) {
      std::sort(a.begin(), a.end());
      a.erase(std::unique(a.begin(), a.end()), a.end());
    }
    // ------------------------------------------------------------
    // Build mapping for zero-weight literals:
    //
    // original CM literal -> index in zeroLiteral
    // ------------------------------------------------------------

    graph.zeroLocalOfLiteral.assign(2 * nVars, -1);

    for (int z = 0;
	 z < (int)graph.zeroLiteral.size();
	 ++z) {
      graph.zeroLocalOfLiteral[graph.zeroLiteral[z]] = z;
    }


    // ------------------------------------------------------------
    // Build neighbourhoods for zero-weight literals.
    //
    // zeroAdj[z]:
    //     neighbours among POSITIVE-weight local nodes.
    //
    // zeroZeroAdj[z]:
    //     neighbours among ZERO-weight nodes, using zero-local
    //     numbering.
    // ------------------------------------------------------------

    graph.zeroAdj.resize(graph.zeroLiteral.size());
    graph.zeroZeroAdj.resize(graph.zeroLiteral.size());

    for (int z = 0;
	 z < (int)graph.zeroLiteral.size();
	 ++z) {
      const int lit =
        graph.zeroLiteral[z];
      for (int e = ConflictGraph.last[lit];
	   e >= 0;
	   e = ConflictGraph.next[e]) {
        const int otherLit =
	  ConflictGraph.adj[e];
        if (otherLit < 0 ||
            otherLit >= 2 * nVars)
	  continue;

        // ----------------------------------------
        // zero -> positive-weight neighbour
        // ----------------------------------------

        const int v =
	  graph.localOfLiteral[otherLit];
        if (v >= 0)
	  graph.zeroAdj[z].push_back(v);

        // ----------------------------------------
        // zero -> zero-weight neighbour
        // ----------------------------------------

        const int zz =
	  graph.zeroLocalOfLiteral[otherLit];
        if (zz >= 0 && zz != z)
	  graph.zeroZeroAdj[z].push_back(zz);
      }
      // Clean positive neighbourhood.

      std::sort(
		graph.zeroAdj[z].begin(),
		graph.zeroAdj[z].end()
		);

      graph.zeroAdj[z].erase(
			     std::unique(
					 graph.zeroAdj[z].begin(),
					 graph.zeroAdj[z].end()
					 ),
			     graph.zeroAdj[z].end()
			     );


      // Clean zero neighbourhood.

      std::sort(
		graph.zeroZeroAdj[z].begin(),
		graph.zeroZeroAdj[z].end()
		);

      graph.zeroZeroAdj[z].erase(
				 std::unique(
					     graph.zeroZeroAdj[z].begin(),
					     graph.zeroZeroAdj[z].end()
					     ),
				 graph.zeroZeroAdj[z].end()
				 );
    }
    if (USE_CLIQUE_OUT) std::cerr
      << "CLIQUESEP COMPLEMENTS:"
      << " pairs=" << complementPairs
      << " already=" << complementAlready
      << " added=" << complementAdded
      << std::endl;
    long long zeroZeroDirectedEdges = 0;
    for (const auto& a : graph.zeroZeroAdj)
      zeroZeroDirectedEdges += (long long)a.size();

    if (USE_CLIQUE_OUT) std::cerr
      << "CLIQUESEP ZERO:"
      << " zeroLiterals=" << graph.zeroLiteral.size()
      << " zeroZeroEdges=" << zeroZeroDirectedEdges / 2
      << std::endl;
}

void CliqueManager::separateCliquesHighsLike(
					     const std::vector<double>& xlp,
					     const int* types,
					     std::vector<CliqueSepClique>& result,
					     double minWeight,
					     int maxCalls,
					     int maxCliques)
{
  result.clear();

  CliqueSepGraph graph;

  buildCliqueSepGraph(
		      xlp,
		      types,
		      graph
		      );

  std::vector<int> P;
  P.reserve(graph.literal.size());
  int maxDegree = 0;
  for (int u = 0; u < (int)graph.literal.size(); ++u) {
    const int degree = (int)graph.adj[u].size();
    if (degree == 0)
      continue;
    P.push_back(u);
    if (degree > maxDegree)
      maxDegree = degree;
  }
  std::sort(P.begin(), P.end(), [&](int a, int b) { return graph.weight[a] > graph.weight[b]; });
  if (USE_CLIQUE_OUT) std::cerr
    << "CLIQUESEP SEARCH:"
    << " candidates=" << P.size()
    << " maxDegree=" << maxDegree
    << std::endl;

  std::vector<int> R;
  std::vector<int> X;
  int calls = 0;
  bronKerboschHighsLike(
			graph,
			R,
			P,
			X,
			0.0,
			minWeight,
			maxCalls,
			maxCliques,
			calls,
			result
			);
  int totalAddedZero = 0;
  int maxAddedZero = 0;

  for (auto& clique : result) {

    // ------------------------------------------------------------
    // Convert the positive-weight clique literals back to
    // local positive graph nodes.
    // ------------------------------------------------------------

    std::vector<int> cliqueNodes;
    cliqueNodes.reserve(clique.literals.size());

    for (int lit : clique.literals) {

      if (lit < 0 ||
	  lit >= (int)graph.localOfLiteral.size())
	continue;

      const int v =
	graph.localOfLiteral[lit];

      if (v >= 0)
	cliqueNodes.push_back(v);
    }


    // ------------------------------------------------------------
    // Collect zero-weight literals that conflict with every
    // member of the existing clique.
    //
    // Note:
    // ConflictGraph does not explicitly contain the trivial
    // complement edge. Therefore literals belonging to the same
    // binary variable are handled explicitly.
    // ------------------------------------------------------------

    std::vector<int> candidates;

    for (int z = 0;
         z < (int)graph.zeroLiteral.size();
         ++z) {

      const int zeroLit =
	graph.zeroLiteral[z];

      bool fits = true;

      for (int cliqueLit : clique.literals) {

	const int v =
	  graph.localOfLiteral[cliqueLit];

	if (v < 0) {
	  fits = false;
	  break;
	}

	// Same binary variable:
	// zeroLit and cliqueLit are complementary literals.
	if ((zeroLit >> 1) == (cliqueLit >> 1))
	  continue;

	if (!std::binary_search(
				graph.zeroAdj[z].begin(),
				graph.zeroAdj[z].end(),
				v)) {

	  fits = false;
	  break;
	}
      }

      if (fits)
	candidates.push_back(z);
    }


    // ------------------------------------------------------------
    // Greedy maximal clique inside the zero-weight candidates.
    //
    // High zero-degree first tends to leave more possibilities
    // for subsequent additions.
    // ------------------------------------------------------------

    std::sort(
	      candidates.begin(),
	      candidates.end(),
	      [&](int a, int b) {
		return graph.zeroZeroAdj[a].size()
		  > graph.zeroZeroAdj[b].size();
	      }
	      );

    std::vector<int> selectedZero;

    for (int z : candidates) {

      bool fits = true;

      for (int zz : selectedZero) {

	if (!std::binary_search(
				graph.zeroZeroAdj[z].begin(),
				graph.zeroZeroAdj[z].end(),
				zz)) {

	  fits = false;
	  break;
	}
      }

      if (!fits)
	continue;

      selectedZero.push_back(z);
    }


    // ------------------------------------------------------------
    // Extend the result clique.
    // Weight stays unchanged because these literals have LP
    // weight zero.
    // ------------------------------------------------------------

    const int added =
      (int)selectedZero.size();

    for (int z : selectedZero)
      clique.literals.push_back(
				graph.zeroLiteral[z]
				);

    totalAddedZero += added;

    if (added > maxAddedZero)
      maxAddedZero = added;
  }

  if (USE_CLIQUE_OUT) std::cerr
    << "CLIQUESEP EXTEND:"
    << " totalAddedZero=" << totalAddedZero
    << " maxAddedZero=" << maxAddedZero
    << std::endl;
  if (USE_CLIQUE_OUT) std::cerr
    << "CLIQUESEP RESULT:"
    << " calls=" << calls
    << " cliques=" << result.size()
    << std::endl;

  for (int i = 0; i < (int)result.size() && i < 10; ++i) {

    if (USE_CLIQUE_OUT) std::cerr
      << "  clique " << i
      << " size=" << result[i].literals.size()
      << " weight=" << result[i].weight
      << std::endl;
  }
    
  long long numDirectedEdges = 0;
  for (const auto& a : graph.adj)
    numDirectedEdges += (long long)a.size();

  if (USE_CLIQUE_OUT) std::cerr
    << "CLIQUESEP:"
    << " nVars=" << nVars
    << " localNodes=" << graph.literal.size()
    << " localEdges=" << numDirectedEdges / 2
    << " minWeight=" << minWeight
    << " maxCalls=" << maxCalls
    << " maxCliques=" << maxCliques
    << std::endl;
}
void CliqueManager::bronKerboschHighsLike(
    const CliqueSepGraph& graph,
    std::vector<int>& R,
    std::vector<int> P,
    std::vector<int> X,
    double currentWeight,
    double minWeight,
    int maxCalls,
    int maxCliques,
    int& calls,
    std::vector<CliqueSepClique>& result)
{
    if (calls >= maxCalls)
        return;

    if ((int)result.size() >= maxCliques)
        return;

    ++calls;

    // ------------------------------------------------------------
    // Simple weight bound:
    //
    // Even if every node in P could be added, can we still
    // exceed minWeight?
    // ------------------------------------------------------------

    double possibleWeight = currentWeight;

    for (int v : P)
        possibleWeight += graph.weight[v];

    if (possibleWeight <= minWeight)
        return;

    // ------------------------------------------------------------
    // R is maximal iff neither P nor X contains an extension.
    // ------------------------------------------------------------

    if (P.empty() && X.empty()) {

        if (currentWeight > minWeight) {

            CliqueSepClique clique;
            clique.weight = currentWeight;
            clique.literals.reserve(R.size());

            for (int v : R)
                clique.literals.push_back(graph.literal[v]);

            result.push_back(std::move(clique));
        }

        return;
    }

    // ------------------------------------------------------------
    // Choose a Bron-Kerbosch pivot from P union X.
    //
    // We prefer the vertex having most neighbours inside P.
    // ------------------------------------------------------------

    int pivot = -1;
    int bestCount = -1;

    auto countNeighboursInP =
        [&](int u)
        {
            int cnt = 0;

            for (int v : P) {

                if (std::binary_search(
                        graph.adj[u].begin(),
                        graph.adj[u].end(),
                        v))
                    ++cnt;
            }

            return cnt;
        };

    for (int u : P) {

        const int cnt = countNeighboursInP(u);

        if (cnt > bestCount) {
            bestCount = cnt;
            pivot = u;
        }
    }

    for (int u : X) {

        const int cnt = countNeighboursInP(u);

        if (cnt > bestCount) {
            bestCount = cnt;
            pivot = u;
        }
    }

    // ------------------------------------------------------------
    // Branch only over P \ N(pivot).
    // ------------------------------------------------------------

    std::vector<int> branchVertices;
    branchVertices.reserve(P.size());

    for (int v : P) {

        const bool neighbourOfPivot =
            pivot >= 0 &&
            std::binary_search(
                graph.adj[pivot].begin(),
                graph.adj[pivot].end(),
                v);

        if (!neighbourOfPivot)
            branchVertices.push_back(v);
    }

    // Prefer high LP-weight literals first.
    std::sort(
        branchVertices.begin(),
        branchVertices.end(),
        [&](int a, int b) {
            return graph.weight[a] > graph.weight[b];
        }
    );

    for (int v : branchVertices) {

        if (calls >= maxCalls ||
            (int)result.size() >= maxCliques)
            return;

        // v may already have been moved P -> X by an earlier branch.
        auto itP =
            std::find(P.begin(), P.end(), v);

        if (itP == P.end())
            continue;

        std::vector<int> newP;
        std::vector<int> newX;

        newP.reserve(P.size());
        newX.reserve(X.size());

        for (int u : P) {
            if (std::binary_search(
                    graph.adj[v].begin(),
                    graph.adj[v].end(),
                    u))
                newP.push_back(u);
        }

        for (int u : X) {
            if (std::binary_search(
                    graph.adj[v].begin(),
                    graph.adj[v].end(),
                    u))
                newX.push_back(u);
        }

        R.push_back(v);

        bronKerboschHighsLike(
            graph,
            R,
            std::move(newP),
            std::move(newX),
            currentWeight + graph.weight[v],
            minWeight,
            maxCalls,
            maxCliques,
            calls,
            result
        );

        R.pop_back();

        // Standard Bron-Kerbosch:
        //
        //     P := P \ {v}
        //     X := X union {v}

        P.erase(itP);
        X.push_back(v);
    }
}
