#include "YasolSymmetry.h"
#include "QBPSolver.h"

void YasolSymmetry::clear() {
    graph.clear();
    varNode.clear();
    rowNode.clear();
    //permutations.clear();

    colorMap.clear();
    nextColor = 1;

    symPartition.cells.clear();
}

int YasolSymmetry::addNode(int color) {
    int id = (int) graph.size();
    graph.push_back(Node{color, {}});
    return id;
}

void YasolSymmetry::addEdge(int a, int b) {
    graph[a].adj.push_back(b);
    graph[b].adj.push_back(a);
}

long long YasolSymmetry::encodeDouble(double x) {
    const double scale = 1e9;
    return (long long) std::llround(x * scale);
}

int YasolSymmetry::getColorId(const std::vector<long long>& key) {
    auto it = colorMap.find(key);

    if (it != colorMap.end())
        return it->second;

    int id = nextColor++;
    colorMap[key] = id;
    return id;
}

int YasolSymmetry::makeVarColor(
    int type,
    int block,
    double lb,
    double ub,
    double obj
) {
    return getColorId({
        1,
        type,
        block,
        encodeDouble(lb),
        encodeDouble(ub),
        encodeDouble(obj)
    });
}

int YasolSymmetry::makeRowColor(
    double rhs,
    int nnz,
    double maxCoeff
) {
    return getColorId({
        2,
        encodeDouble(rhs),
        nnz,
        encodeDouble(maxCoeff)
    });
}

int YasolSymmetry::makeCoeffColor(double coeff) {
    return getColorId({
        3,
        encodeDouble(coeff)
    });
}

YasolSymmetryResult YasolSymmetry::run(QBPSolver& qbp) {
    int n = qbp.nVars();

    if (persistentNumVars != n) {
        persistentResult = YasolSymmetryResult();
        persistentNumVars = n;
    }

    YasolSymmetryResult& result = persistentResult;

    clear();

    if (n == 0 || qbp.getNumberOfConstraints() <= 1)
        return result;

    // Sicherheitsbudgets fuer grosse Instanzen.
    // Negative Werte wuerden das jeweilige Limit deaktivieren.
    groupSearch.maxSearchNodes = 100000;
    groupSearch.maxLeaves = 100000;
    groupSearch.maxDepth = 10000;
    groupSearch.maxGenerators = 3000;
    
    groupSearch.infoout = false;
    groupSearch.validateGenerators = true;

    groupSearch.loadFromQbp(qbp, 1e-9);
    groupSearch.initializeHashValues();

    if (!groupSearch.partitionRefinement())
        return result;

    groupSearch.removeFixPoints();
    groupSearch.initializeGroundSet();

    groupSearch.printOrbitCells();
    groupSearch.printSearchStats("after ground set");

    groupSearch.permutations.clear();
    groupSearch.storedAutomorphisms.clear();
    groupSearch.numStoredAutomorphisms = 0;
    groupSearch.searchAutomorphisms();
    //groupSearch.searchAutomorphismsMultiStart();
    groupSearch.printSearchStats("after searchAutomorphisms");
    
    /*
      Ergänzt sichere Generatoren aus verbliebenen Zellen.
      Durch addPermutationGenerator() sind sie validiert.
    */
    //groupSearch.extractTrivialGenerators();
    //groupSearch.printSearchStats("after extractTrivialGenerators");
    
    //groupSearch.extractOrbits(result);

    YasolSymmetryResult runResult;
    groupSearch.extractOrbits(runResult);
    mergeOrbitResults(result, runResult, qbp.nVars());

    result.orbitIdOfVar.assign(qbp.nVars(), -1);
    result.posInOrbitOfVar.assign(qbp.nVars(), -1);
 
    for (auto& orb : result.orbits) {
      std::sort(orb.begin(), orb.end());
    }
    for (int o = 0; o < (int)result.orbits.size(); ++o) {
      for (int p = 0; p < (int)result.orbits[o].size(); ++p) {
        int v = result.orbits[o][p];
	
        if (v >= 0 && v < qbp.nVars()) {
	  result.orbitIdOfVar[v] = o;
	  result.posInOrbitOfVar[v] = p;
        }
      }
    }
 
    result.success = !result.orbits.empty();
    result.partial = groupSearch.abortedByBudget;
 
    if (groupSearch.infoout) {
      std::cerr << "\n==== FINAL ORBITS ====\n" << "Number of candidate orbits = " << result.orbits.size() << "\n";
      for (int k = 0; k < (int)result.orbits.size(); ++k) {
	std::cerr << "Orbit " << k << ": ";
	for (int v : result.orbits[k])
	  std::cerr << v << " ";
	std::cerr << "\n";
      }
      std::cerr << "======================\n";
    }
    return result;
}

void YasolSymmetry::buildColoredGraph(QBPSolver& qbp) {

    int n = qbp.nVars();
    int m = qbp.getNumberOfConstraints();

    varNode.resize(n);
    rowNode.resize(m, -1);

    //-----------------------------------
    // Objective lesen: constraints[0]
    //-----------------------------------

    Constraint* obj = qbp.getConstraint(0);

    std::vector<double> objCoeff(n, 0.0);

    for (int k = 0; k < obj->size(); ++k) {
        int j = var((*obj)[k]);

        double a = (*obj)[k].coef;
        if (sign((*obj)[k]) == 1)
            a = -a;

        objCoeff[j] = a;
    }

    //-----------------------------------
    // Variablenknoten
    //-----------------------------------

    for (int j = 0; j < n; ++j) {
        int color = makeVarColor(
            qbp.getType(j),
            qbp.getBlock(j),
            qbp.getLowerBound(j),
            qbp.getUpperBound(j),
            objCoeff[j]
        );

        varNode[j] = addNode(color);
    }

    //-----------------------------------
    // Constraintknoten + Koeffizientenknoten
    //-----------------------------------

    for (int r = 1; r < m; ++r) {
        Constraint* c = qbp.getConstraint(r);

        //--------------------------------
        // max |a_ij| für Zeilenfarbe
        //--------------------------------

        double maxAbsCoeff = 0.0;

        for (int k = 0; k < c->size(); ++k) {
            double a = (*c)[k].coef;   // Betrag
            if (a > maxAbsCoeff)
                maxAbsCoeff = a;
        }

        rowNode[r] = addNode(
            makeRowColor(
                c->header.rhs,
                c->size(),
                maxAbsCoeff
            )
        );

        //--------------------------------
        // Koeffizientenknoten
        //--------------------------------

        for (int k = 0; k < c->size(); ++k) {
            int j = var((*c)[k]);

            double a = (*c)[k].coef;   // Betrag
            if (sign((*c)[k]) == 1)
                a = -a;

            int coeffNode = addNode(
                makeCoeffColor(a)
            );

            addEdge(varNode[j], coeffNode);
            addEdge(coeffNode, rowNode[r]);
        }
    }
}

void YasolSymmetry::initializePartition() {
    symPartition.cells.clear();

    std::map<int, std::vector<int>> byColor;

    for (int v = 0; v < (int)graph.size(); ++v) {
        byColor[graph[v].color].push_back(v);
    }

    for (auto& it : byColor) {
        symPartition.cells.push_back(std::move(it.second));
    }
}


void YasolSymmetry::detectAutomorphisms() {
  bool output = false;//true;
  initializePartition();
  if (output) {
    std::cerr << "partition cells = " << symPartition.cells.size() << "\n";
    for (int i=0; i < symPartition.cells.size(); i++) {
      if (symPartition.cells[i].size() > 1) {
        std::cerr << "cell " << i << ": ";
        for (int v : symPartition.cells[i]) {
	  std::cerr << v << " ";
        }
        std::cerr << "\n";
      }
    }
  }
  refinePartition();
  if (output) {
    std::cerr << "refined partition cells = " << symPartition.cells.size() << "\n";
    int nontrivVarCells = 0;
    std::vector<int> nodeToVar(graph.size(), -1);
    for (int j = 0; j < (int)varNode.size(); ++j)
      nodeToVar[varNode[j]] = j;

    for (int i = 0; i < (int)symPartition.cells.size(); ++i) {
      std::vector<int> vars;

      for (int node : symPartition.cells[i]) {
        if (nodeToVar[node] >= 0)
	  vars.push_back(nodeToVar[node]);
      }

      if (vars.size() >= 2) {
        nontrivVarCells++;

        std::cerr << "refined var cell " << i << ": ";
        for (int v : vars)
	  std::cerr << v << " ";
        std::cerr << "\n";
      }
    }
    std::cerr << "nontrivial refined variable cells = " << nontrivVarCells << "\n";
  }
  extractGenerators();
}

uint64_t YasolSymmetry::mix64(uint64_t x) const {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    x = x ^ (x >> 31);
    return x;
}

YasolSymmetry::Signature YasolSymmetry::computeSignature(int node) {
    uint64_t h1 = 0;
    uint64_t h2 = 0;

    for (int nb : graph[node].adj) {
        uint64_t z = mix64((uint64_t)graph[nb].color);
        h1 += z;
        h2 += z * z;
    }

    Signature s;
    s.ownColor = graph[node].color;
    s.degree = (int)graph[node].adj.size();
    s.h1 = h1;
    s.h2 = h2;

    return s;
}

void YasolSymmetry::refinePartition() {

    bool changed = true;

    while (changed) {

        changed = false;

        std::vector<std::vector<int>> newCells;
        newCells.reserve(symPartition.cells.size());

        for (auto& cell : symPartition.cells) {

            if (cell.size() <= 1) {
                newCells.push_back(std::move(cell));
                continue;
            }

            std::unordered_map<Signature,std::vector<int>,SignatureHash> split;

            split.reserve(cell.size() * 2);

            for (int node : cell) {
                split[computeSignature(node)].push_back(node);
            }

            if (split.size() > 1)
                changed = true;

            for (auto& p : split) {
                newCells.push_back(std::move(p.second));
            }
        }

        symPartition.cells = std::move(newCells);

        for (int c = 0; c < (int)symPartition.cells.size(); ++c) {
            for (int node : symPartition.cells[c]) {
                graph[node].color = c + 1;
            }
        }
    }
}

void YasolSymmetry::extractGenerators() {
    // erstmal leer
    //permutations.clear();
}

void YasolSymmetry::extractOrbits(YasolSymmetryResult& result) {
  //result.orbits.clear();

    std::vector<int> nodeToVar(graph.size(), -1);

    for (int j = 0; j < (int)varNode.size(); ++j) {
        nodeToVar[varNode[j]] = j;
    }

    for (const auto& cell : symPartition.cells) {
        std::vector<int> orbit;

        for (int node : cell) {
            int v = nodeToVar[node];

            if (v >= 0)
                orbit.push_back(v);
        }

        if (orbit.size() >= 2)
            result.orbits.push_back(std::move(orbit));
    }
}


void YasolGroupSearch::loadFromQbp(QBPSolver& qbp, double eps) {
    coeffColorMap.clear();
    nextCoeffColor = 1;
    numCol = qbp.nVars();

    // constraints[0] ist Objective, echte Constraints starten bei 1
    int qbpNumCons = qbp.getNumberOfConstraints();
    numRow = qbpNumCons - 1;

    numVertices = numCol + numRow;
    numActiveCols = numCol;

    fixedVertex.assign(numCol + numRow, 0);

    for (int j = 0; j < numCol; ++j) {
      if (qbp.isFixed(j) || qbp.getAssignment(j)!=2) {
        fixedVertex[j] = 1;
      }
    }

    Gstart.assign(numVertices + 1, 0);
    Gend.clear();
    Gedge.clear();

    vertexToCell.assign(numVertices, 0);
    currentPartition.clear();
    currentPartitionLinks.clear();

    //--------------------------------------
    // Objective-Koeffizienten sammeln
    //--------------------------------------

    std::vector<double> objCoeff(numCol, 0.0);

    Constraint* obj = qbp.getConstraint(0);

    for (int k = 0; k < obj->size(); ++k) {
        int j = var((*obj)[k]);

        double a = (*obj)[k].coef;
        if (sign((*obj)[k]) == 1)
            a = -a;

        objCoeff[j] = a;
    }

    //--------------------------------------
    // Zeilen- und Spaltenzählung
    //--------------------------------------

    std::vector<int> colSizes(numCol, 0);
    std::vector<int> rowSizes(numRow, 0);

    int numNz = 0;

    for (int r = 1; r < qbpNumCons; ++r) {
        Constraint* c = qbp.getConstraint(r);
        int gr = r - 1;

        rowSizes[gr] = c->size();
        numNz += c->size();

        for (int k = 0; k < c->size(); ++k) {
            int j = var((*c)[k]);
            colSizes[j]++;
        }
    }

    //--------------------------------------
    // Gstart für Spalten
    //--------------------------------------

    int offset = 0;

    for (int j = 0; j < numCol; ++j) {
        Gstart[j] = offset;
        offset += colSizes[j];
    }

    //--------------------------------------
    // Gstart für Zeilen
    //--------------------------------------

    for (int r = 0; r < numRow; ++r) {
        Gstart[numCol + r] = offset;
        offset += rowSizes[r];
    }

    Gstart[numVertices] = offset;

    // Jede Matrix-NZ erscheint zweimal: col->row und row->col
    Gedge.assign(2 * numNz, std::make_pair(-1, 0u));

    //--------------------------------------
    // temporäre Schreibpositionen
    //--------------------------------------

    std::vector<int> colWrite(numCol);
    std::vector<int> rowWrite(numRow);

    for (int j = 0; j < numCol; ++j)
        colWrite[j] = Gstart[j];

    for (int r = 0; r < numRow; ++r)
        rowWrite[r] = Gstart[numCol + r];

    //--------------------------------------
    // Kanten eintragen
    //--------------------------------------

    for (int r = 1; r < qbpNumCons; ++r) {
        Constraint* c = qbp.getConstraint(r);
        int gr = r - 1;
        int rowVertex = numCol + gr;

        for (int k = 0; k < c->size(); ++k) {
            int j = var((*c)[k]);

            double a = (*c)[k].coef;
            if (sign((*c)[k]) == 1)
                a = -a;

            uint32_t edgeColor =
                (uint32_t) makeCoeffColor(a);

            int pc = colWrite[j]++;
            Gedge[pc] = std::make_pair(rowVertex, edgeColor);

            int pr = rowWrite[gr]++;
            Gedge[pr] = std::make_pair(j, edgeColor);
        }
    }

    Gend.assign(Gstart.begin() + 1, Gstart.end());

    //--------------------------------------
    // Startzellen/Farben für Spalten
    //--------------------------------------

    std::map<std::vector<long long>, int> cellMap;
    int nextCell = 1;

    auto getCell = [&](const std::vector<long long>& key) -> int {
        auto it = cellMap.find(key);
        if (it != cellMap.end())
            return it->second;

        int id = nextCell++;
        cellMap[key] = id;
        return id;
    };

    for (int j = 0; j < numCol; ++j) {
        std::vector<long long> key;

        key.push_back(1); // column marker
        key.push_back(qbp.getType(j));
        key.push_back(qbp.getBlock(j));
        key.push_back(encodeDouble(qbp.getLowerBound(j)));
        key.push_back(encodeDouble(qbp.getUpperBound(j)));
        key.push_back(encodeDouble(objCoeff[j]));
        key.push_back(colSizes[j]);

        vertexToCell[j] = getCell(key);
    }

    //--------------------------------------
    // Startzellen/Farben für Zeilen
    //--------------------------------------

    for (int r = 1; r < qbpNumCons; ++r) {
        Constraint* c = qbp.getConstraint(r);
        int gr = r - 1;

        double maxAbsCoeff = 0.0;

        for (int k = 0; k < c->size(); ++k) {
            double a = (*c)[k].coef;
            if (a > maxAbsCoeff)
                maxAbsCoeff = a;
        }

        std::vector<long long> key;

        key.push_back(2); // row marker
        key.push_back(encodeDouble(c->header.rhs));
        key.push_back(c->size());
        key.push_back(encodeDouble(maxAbsCoeff));

        vertexToCell[numCol + gr] = getCell(key);
    }

    //--------------------------------------
    // currentPartition nach Zellen sortieren
    //--------------------------------------

    currentPartition.resize(numVertices);

    for (int v = 0; v < numVertices; ++v)
        currentPartition[v] = v;

    std::sort(
        currentPartition.begin(),
        currentPartition.end(),
        [&](int v1, int v2) {
            if (vertexToCell[v1] != vertexToCell[v2])
                return vertexToCell[v1] < vertexToCell[v2];
            return v1 < v2;
        }
    );

    //--------------------------------------
    // currentPartitionLinks + vertexToCell
    // auf Zellstart-Indizes umstellen
    //--------------------------------------

    currentPartitionLinks.resize(numVertices);

    int cellStart = 0;
    int cellNumber = vertexToCell[currentPartition[0]];

    for (int i = 0; i < numVertices; ++i) {
        int vertex = currentPartition[i];

        if (vertexToCell[vertex] != cellNumber) {
            currentPartitionLinks[cellStart] = i;
            cellStart = i;
            cellNumber = vertexToCell[vertex];
        }

        vertexToCell[vertex] = cellStart;
        currentPartitionLinks[i] = cellStart;
    }

    currentPartitionLinks[cellStart] = numVertices;

    bool output = false;//true;
    if (output) {
      std::cerr	<< "\n==== GROUP SEARCH ====\n";

      std::cerr << "numCol = " << numCol << "\n";

      std::cerr << "numRow = " << numRow << "\n";

      std::cerr << "numVertices = " << numVertices << "\n";

      std::cerr << "numEdges = " << Gedge.size() << "\n";

      std::cerr << "initial partition cells = ";

      int ncells = 0;

      for (int i = 0; i < numVertices;) {
	++ncells;
	i = currentPartitionLinks[i];
      }
      std::cerr << ncells << "\n";

      if (infoout) {
	for (int i = 0; i < std::min(numVertices,20);++i) {
	  std::cerr << "v " << i << " -> cell " << vertexToCell[i] << "\n";
	}
      }
    }    
}

long long YasolGroupSearch::encodeDouble(double x) const {
    const double scale = 1e9;
    return (long long) std::llround(x * scale);
}

int YasolGroupSearch::makeCoeffColor(double coeff) {
    std::vector<long long> key;
    key.push_back(encodeDouble(coeff));

    auto it = coeffColorMap.find(key);
    if (it != coeffColorMap.end())
        return it->second;

    int id = nextCoeffColor++;
    coeffColorMap[key] = id;
    return id;
}

int YasolGroupSearch::cellSize(int cell) const {
    return currentPartitionLinks[cell] - cell;
}

void YasolGroupSearch::markCellForRefinement(int cell) {
    if (cellSize(cell) <= 1)
        return;
    if ( cellInRefinementQueue[cell] )
        return;
    cellInRefinementQueue[cell] = true;
    refinementQueue.push_back(cell);
}

void YasolGroupSearch::initializeHashValues() {
    vertexHash.clear();
    cellInRefinementQueue.assign(numVertices,false);
    refinementQueue.clear();
    for (int v = 0; v < numVertices; v++) {
        uint32_t h = 0;
        int cell = vertexToCell[v];
        for (int e = Gstart[v]; e < Gend[v]; e++) {
            uint32_t c = (uint32_t) vertexToCell[ Gedge[e].first ];
            c ^= Gedge[e].second;
            h = h * 2654435761u + c;
        }
        vertexHash[v] = h;
        markCellForRefinement( cell );
    }
}

bool YasolGroupSearch::partitionRefinement() {
    while (!refinementQueue.empty()) {
        int cell = refinementQueue.back();
        refinementQueue.pop_back();

        if (cell < 0 || cell >= (int)cellInRefinementQueue.size())
            continue;

        cellInRefinementQueue[cell] = false;

        if (cell < 0 || cell >= numVertices)
            continue;

        if (cellSize(cell) <= 1)
            continue;

        int begin = cell;
        int end = currentPartitionLinks[cell];

        if (begin < 0 || end > numVertices || begin >= end)
            continue;

        std::sort(
            currentPartition.begin() + begin,
            currentPartition.begin() + end,
            [&](int a, int b) {
                uint32_t ha = 0;
                uint32_t hb = 0;

                auto ia = vertexHash.find(a);
                if (ia != vertexHash.end())
                    ha = ia->second;

                auto ib = vertexHash.find(b);
                if (ib != vertexHash.end())
                    hb = ib->second;

                if (ha != hb)
                    return ha < hb;

                return a < b;
            }
        );

        int cellStart = begin;

        uint32_t lastHash = 0;
        {
            int v0 = currentPartition[begin];
            auto it = vertexHash.find(v0);
            if (it != vertexHash.end())
                lastHash = it->second;
        }

        for (int i = begin + 1; i < end; ++i) {
            int v = currentPartition[i];

            uint32_t h = 0;
            auto it = vertexHash.find(v);
            if (it != vertexHash.end())
                h = it->second;

            if (h != lastHash) {
                if (!splitCell(cellStart, i)) {
                    // Pruning durch Prefix-Vergleich.
                    // Queue und Hashwerte sauber leeren.
                    for (int c : refinementQueue) {
                        if (c >= 0 &&
                            c < (int)cellInRefinementQueue.size())
                            cellInRefinementQueue[c] = false;
                    }

                    refinementQueue.clear();
                    vertexHash.clear();

                    return false;
                }

                cellStart = i;
                lastHash = h;
            }

            updateCellMembership(i, cellStart, true);
        }
    }

    vertexHash.clear();

    return true;
}

void YasolGroupSearch::removeFixPoints() {

    std::vector<int> oldPartition = std::move(currentPartition);
    std::vector<int> oldLinks = std::move(currentPartitionLinks);

    int oldNumVertices = numVertices;

    currentPartition.clear();

    //----------------------------------
    // 1. Nur relevante nichttriviale Restzellen kopieren.
    //
    // Fixierte Variablen bleiben im vollen Graphen erhalten,
    // werden aber aus der aktiven Suchpartition entfernt.
    //----------------------------------

    std::vector<int> newCellStarts;

    int skippedFixedVars = 0;
    int skippedSingletons = 0;

    for (int c = 0; c < oldNumVertices;) {
        int end = oldLinks[c];

        std::vector<int> kept;
        kept.reserve(end - c);

        for (int i = c; i < end; ++i) {
            int v = oldPartition[i];

            // Fixierte Variablen nicht mehr individualisieren.
            // Sie bleiben aber im vollen Graphen vorhanden und
            // werden in isPermutationAutomorphism() implizit festgehalten.
            if (v < numCol && isFixedVariableVertex(v)) {
                ++skippedFixedVars;
                continue;
            }

            kept.push_back(v);
        }

        if ((int)kept.size() > 1) {
            int newStart = (int)currentPartition.size();
            newCellStarts.push_back(newStart);

            for (int v : kept)
                currentPartition.push_back(v);
        } else {
            skippedSingletons += (int)kept.size();
        }

        c = end;
    }

    //----------------------------------
    // 2. Neue aktive Suchgröße setzen
    //----------------------------------

    numVertices = (int)currentPartition.size();

    //----------------------------------
    // 3. vertexToCell zurücksetzen
    //----------------------------------

    for (int i = 0; i < (int)vertexToCell.size(); ++i)
        vertexToCell[i] = -1;

    //----------------------------------
    // 4. Links neu aufbauen
    //----------------------------------

    currentPartitionLinks.assign(numVertices, -1);

    for (int idx = 0; idx < (int)newCellStarts.size(); ++idx) {
        int start = newCellStarts[idx];
        int end = (idx + 1 < (int)newCellStarts.size())
            ? newCellStarts[idx + 1]
            : numVertices;

        currentPartitionLinks[start] = end;

        for (int i = start; i < end; ++i) {
            int v = currentPartition[i];
            vertexToCell[v] = start;

            if (i != start)
                currentPartitionLinks[i] = start;
        }
    }

    //----------------------------------
    // 5. Aktive Spalten zählen
    //----------------------------------

    numActiveCols = 0;

    for (int v : currentPartition) {
        if (v < numCol)
            ++numActiveCols;
    }

    //----------------------------------
    // 6. Queue/Hashes leeren
    //----------------------------------

    vertexHash.clear();
    refinementQueue.clear();
    cellInRefinementQueue.assign(numVertices, false);

    if (infoout) {
        std::cerr << "removeFixPoints: skipped fixed vars = " << skippedFixedVars << ", skipped singleton leftovers = " << skippedSingletons << ", active vertices = " << numVertices << "\n";
    }
}

void YasolGroupSearch::initializeGroundSet() {
    vertexGroundSet.clear();
    vertexGroundSet.reserve(currentPartition.size());

    for (int v : currentPartition) {
        if (v < 0 || v >= numCol + numRow)
            continue;

        // Sollte nach removeFixPoints() eigentlich nicht mehr vorkommen,
        // aber als Sicherheitsgurt drinlassen.
        if (v < numCol && isFixedVariableVertex(v))
            continue;

        vertexGroundSet.push_back(v);
    }

    std::sort(vertexGroundSet.begin(), vertexGroundSet.end());

    vertexPosition.assign(numCol + numRow, -1);

    for (int i = 0; i < (int)vertexGroundSet.size(); ++i) {
        vertexPosition[vertexGroundSet[i]] = i;
    }

    // Online-Orbits fuer die Generatorsuche.
    // Diese DSU dient nur dazu, redundante Leaf-Permutationen nicht
    // auf das maxGenerators-Budget zaehlen zu lassen. extractOrbits()
    // rekonstruiert die finalen Orbits weiterhin aus den gespeicherten
    // Permutationen.
    orbitPartition.resize(numCol);
    orbitSize.assign(numCol, 1);
    for (int j = 0; j < numCol; ++j)
        orbitPartition[j] = j;
}

int YasolGroupSearch::getColumnOrbitRoot(int col) {
    if (col < 0 || col >= numCol)
        return -1;

    int root = orbitPartition[col];
    while (root != orbitPartition[root])
        root = orbitPartition[root];

    while (col != root) {
        int next = orbitPartition[col];
        orbitPartition[col] = root;
        col = next;
    }

    return root;
}

bool YasolGroupSearch::mergeColumnOrbits(int col1, int col2) {
    if (col1 < 0 || col1 >= numCol || col2 < 0 || col2 >= numCol)
        return false;

    int root1 = getColumnOrbitRoot(col1);
    int root2 = getColumnOrbitRoot(col2);

    if (root1 < 0 || root2 < 0 || root1 == root2)
        return false;

    if (orbitSize[root1] < orbitSize[root2])
        std::swap(root1, root2);

    orbitPartition[root2] = root1;
    orbitSize[root1] += orbitSize[root2];

    return true;
}

void YasolGroupSearch::extractTrivialGenerators() {
    for (int c = 0; c < numVertices;) {
        int end = currentPartitionLinks[c];
        int size = end - c;

        if (size >= 2) {
            int baseVertex = currentPartition[c];
            int p1 = vertexPosition[baseVertex];

            if (p1 >= 0) {
                for (int i = c + 1; i < end; ++i) {
                    int otherVertex = currentPartition[i];
                    int p2 = vertexPosition[otherVertex];

                    if (p2 < 0)
                        continue;

                    std::vector<int> perm = vertexGroundSet;
                    std::swap(perm[p1], perm[p2]);

                    addPermutationGenerator(std::move(perm), "trivial-cell");
                }
            }
        }
        c = end;
    }
}

void YasolGroupSearch::extractOrbits(YasolSymmetryResult& result) {
  //result.orbits.clear();

    std::vector<int> parent(numCol);
    std::vector<int> size(numCol, 1);

    for (int j = 0; j < numCol; ++j)
        parent[j] = j;

    auto find = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    auto unite = [&](int a, int b) {
        int ra = find(a);
        int rb = find(b);

        if (ra == rb)
            return;

        if (size[ra] < size[rb])
            std::swap(ra, rb);

        parent[rb] = ra;
        size[ra] += size[rb];
    };

    for (const auto& perm : permutations) {
        for (int p = 0; p < (int)perm.size(); ++p) {
            int v1 = vertexGroundSet[p];
            int v2 = perm[p];

            if (v1 < numCol && v2 < numCol && v1 != v2)
                unite(v1, v2);
        }
    }

    std::map<int, std::vector<int>> classes;

    for (int j = 0; j < numCol; ++j) {
        int r = find(j);
        classes[r].push_back(j);
    }

    for (auto& it : classes) {
        if (it.second.size() >= 2)
            result.orbits.push_back(std::move(it.second));
    }
}

void YasolGroupSearch::createNode() {
    SearchNode node;
    ++numSearchNodes;

    node.stackStart = (int)cellCreationStack.size();
    node.certificateEnd = (int)currNodeCertificate.size();
    node.targetCell = -1;
    node.lastDistinguished = -1;

    node.branchCandidates.clear();
    node.nextBranchCandidate = 0;

    nodeStack.push_back(std::move(node));
}

int YasolGroupSearch::selectTargetCell() {
    int i = 0;

    if (nodeStack.size() > 1)
        i = nodeStack[nodeStack.size() - 2].targetCell;

    while (i < numVertices) {
        if (cellSize(i) > 1)
            return i;
        ++i;
    }

    return -1;
}

bool YasolGroupSearch::splitCell(int cell, int splitPoint) {
    const int oldEnd = currentPartitionLinks[cell];

    if (splitPoint <= cell || splitPoint >= oldEnd)
        return false;

    const int leftSize = splitPoint - cell;
    const int rightSize = oldEnd - splitPoint;

    uint32_t hLeft = 0x12345678u;
    uint32_t hRight = 0x9abcdef0u;

    for (int i = cell; i < splitPoint; ++i) {
        int v = currentPartition[i];

        uint32_t hv = 0;
        auto it = vertexHash.find(v);
        if (it != vertexHash.end())
            hv = it->second;

        hLeft = mixHash32(
            hLeft,
            (uint32_t)v,
            hv
        );
    }

    for (int i = splitPoint; i < oldEnd; ++i) {
        int v = currentPartition[i];

        uint32_t hv = 0;
        auto it = vertexHash.find(v);
        if (it != vertexHash.end())
            hv = it->second;

        hRight = mixHash32(
            hRight,
            (uint32_t)v,
            hv
        );
    }

    uint32_t cert = 0xdeadbeefu;

    cert = mixHash32(
        cert,
        (uint32_t)cell,
        (uint32_t)splitPoint
    );

    cert = mixHash32(
        cert,
        (uint32_t)leftSize,
        hLeft
    );

    cert = mixHash32(
        cert,
        (uint32_t)rightSize,
        hRight
    );

    // Prefix-Pruning.
    // currNodeCertificate enthält cert hier noch NICHT.
    if (!firstLeafCertificate.empty()) {
        int pos = (int)currNodeCertificate.size();

        if (firstLeafPrefixLen == pos &&
            pos < (int)firstLeafCertificate.size() &&
            cert == firstLeafCertificate[pos]) {
            ++firstLeafPrefixLen;
        }

        if (bestLeafPrefixLen == pos &&
            pos < (int)bestLeafCertificate.size() &&
            cert == bestLeafCertificate[pos]) {
            ++bestLeafPrefixLen;
        }

        if (firstLeafPrefixLen <= pos &&
            bestLeafPrefixLen <= pos) {

            uint32_t diffVal;

            if (bestLeafPrefixLen == pos)
                diffVal = cert;
            else
                diffVal = currNodeCertificate[bestLeafPrefixLen];

            if (diffVal > bestLeafCertificate[bestLeafPrefixLen]) {
                ++numPruned;
                return false;
            }
        }
    }

    currentPartitionLinks[splitPoint] = oldEnd;
    currentPartitionLinks[cell] = splitPoint;

    cellCreationStack.push_back(splitPoint);
    currNodeCertificate.push_back(cert);

    return true;
}
bool YasolGroupSearch::updateCellMembership(
    int pos,
    int cell,
    bool markForRefinement
) {
    int vertex = currentPartition[pos];

    if (vertexToCell[vertex] == cell)
        return false;

    vertexToCell[vertex] = cell;

    if (pos != cell)
        currentPartitionLinks[pos] = cell;

    if (markForRefinement) {
        for (int e = Gstart[vertex]; e < Gend[vertex]; ++e) {
            int nb = Gedge[e].first;
            int nbCell = vertexToCell[nb];

            if (nbCell < 0)
                continue;

            if (cellSize(nbCell) <= 1)
                continue;

            // einfache Hash-Aktualisierung
            uint32_t h = vertexHash[nb];
            h = h * 2654435761u
                + (uint32_t)(cell ^ Gedge[e].second);
            vertexHash[nb] = h;

            markCellForRefinement(nbCell);
        }
    }

    return true;
}

void YasolGroupSearch::initializeBranchCandidates(SearchNode& node) {
    node.branchCandidates.clear();
    node.nextBranchCandidate = 0;

    const int cell = node.targetCell;

    if (cell < 0)
        return;

    const int end = currentPartitionLinks[cell];

    for (int i = cell; i < end; ++i)
        node.branchCandidates.push_back(currentPartition[i]);

    // Kandidaten mit hoeherem lokalen Struktur-Score zuerst.
    std::sort(node.branchCandidates.begin(), node.branchCandidates.end(),
              [&](int a, int b) {
                  const double sa = vertexBranchingScore(a);
                  const double sb = vertexBranchingScore(b);

                  if (sa != sb)
                      return sa > sb;

                  return a < b;
              });
}

int YasolGroupSearch::findVertexPositionInCell(int cell, int vertex) const {
    if (cell < 0)
        return -1;

    const int end = currentPartitionLinks[cell];

    for (int i = cell; i < end; ++i) {
        if (currentPartition[i] == vertex)
            return i;
    }

    return -1;
}

bool YasolGroupSearch::determineNextToDistinguish() {
    SearchNode& curr = nodeStack.back();

    distinguishCands.clear();

    if (curr.targetCell < 0)
        return false;

    if (curr.branchCandidates.empty())
        initializeBranchCandidates(curr);

    while (curr.nextBranchCandidate < (int)curr.branchCandidates.size()) {
        const int v = curr.branchCandidates[curr.nextBranchCandidate++];

        // Nach Backtracking kann die Zelle veraendert/restauriert sein.
        // Deshalb suchen wir die aktuelle Position des Kandidaten in der Target-Zelle.
        const int pos = findVertexPositionInCell(curr.targetCell, v);

        if (pos < 0)
            continue;

        // Beim ersten Kandidaten erlauben wir den Referenzast immer.
        // Danach nur noch Kandidaten, die nicht durch bereits gespeicherte
        // Automorphismen als redundant erkannt werden.
        if (curr.lastDistinguished == -1 || checkStoredAutomorphism(v)) {
            distinguishCands.push_back(&currentPartition[pos]);
            return true;
        }
    }

    return false;
}

bool YasolGroupSearch::distinguishVertex(int targetCell) {
    if (distinguishCands.size() != 1)
        return false;

    if (targetCell < 0)
        return false;

    int targetEnd = currentPartitionLinks[targetCell];

    if (targetEnd <= targetCell + 1)
        return false;

    int newCell = targetEnd - 1;

    std::swap(*distinguishCands[0], currentPartition[newCell]);

    nodeStack.back().lastDistinguished =
        currentPartition[newCell];

    if (!splitCell(targetCell, newCell))
        return false;

    updateCellMembership(newCell, newCell);

    return true;
}

std::vector<std::vector<int>> YasolGroupSearch::savePartition() {
    // Deshalb vorerst leerer Platzhalter.
    return std::vector<std::vector<int>>();
}

void YasolGroupSearch::restorePartition(const std::vector<std::vector<int>>& p) {
    (void)p;
    // Platzhalter. Echte Backtracking-Restauration bauen wir danach
    // über cellCreationStack/currNodeCertificate.
}

void YasolGroupSearch::searchAutomorphisms() {
    numSearchNodes = 0;
    numLeaves = 0;
    numPruned = 0;
    numGenerated = 0;
    numGraphMismatches = 0;
    numAutoSkipped = 0;
    abortedByBudget = false;
    numDepthPruned = 0;
    
    //permutations.clear();

    storedAutomorphisms.clear();
    numStoredAutomorphisms = 0;

    firstLeafCertificate.clear();
    bestLeafCertificate.clear();

    firstLeafOrder.clear();
    bestLeafOrder.clear();

    firstLeafGraph.clear();
    bestLeafGraph.clear();

    firstLeafPrefixLen = 0;
    bestLeafPrefixLen = 0;
    firstPathDepth = 0;
    bestPathDepth = 0;

    bestCertificate.clear();

    nodeStack.clear();
    currNodeCertificate.clear();
    cellCreationStack.clear();
    distinguishCands.clear();

    createNode();

    while (!nodeStack.empty()) {

      // Harte globale Budgets: Nodes, Leaves, Generatoren.
      if (searchBudgetExceeded()) {
        abortedByBudget = true;
        break;
      }
      
      // Tiefenlimit: nur aktuellen Ast abschneiden, nicht die ganze Suche.
      if (maxDepth >= 0 && (int)nodeStack.size() >= maxDepth) {
        ++numDepthPruned;
        ++numPruned;
	
        switchToNextNode((int)nodeStack.size() - 1);
        continue;
      }

      int target = selectTargetCell();

      if (target == -1) {
        processCurrentLeaf();
        switchToNextNode((int)nodeStack.size() - 1);
        continue;
      }

      SearchNode& curr = nodeStack.back();

      // Wenn fuer diesen Suchknoten noch keine Target-Zelle gesetzt ist,
      // setzen wir sie jetzt und initialisieren spaeter die sortierten Kandidaten.
      if (curr.targetCell == -1) {
        curr.targetCell = target;
      }

      if (!determineNextToDistinguish()) {
        switchToNextNode((int)nodeStack.size() - 1);
        continue;
      }

      if (!distinguishVertex(curr.targetCell)) {
        switchToNextNode((int)nodeStack.size() - 1);
        continue;
      }

      if (!partitionRefinement()) {
        switchToNextNode((int)nodeStack.size() - 1);
        continue;
      }
      createNode();
    }

    if (infoout) {
      std::cerr << "symmetry search stats:" << " nodes=" << numSearchNodes << " leaves=" << numLeaves << " pruned=" << numPruned << " graphMismatch=" << numGraphMismatches << " generators=" << numGenerated << " autoSkipped=" << numAutoSkipped << "\n";
    }
}

void YasolGroupSearch::descendWithDistinguishedVertex(int targetCell, int v) {
    (void)targetCell;
    (void)v;
    // Platzhalter.
}

std::vector<int> YasolGroupSearch::computeCertificate() {
    std::vector<int> cert;

    for (uint32_t x : currNodeCertificate)
        cert.push_back((int)x);

    return cert;
}

std::vector<int> YasolGroupSearch::extractLeafPermutation() {
    return vertexGroundSet;
}

std::vector<int> YasolGroupSearch::getCellVertices(int cell) const {
    std::vector<int> result;

    if (cell < 0 || cell >= numVertices)
        return result;

    int end = currentPartitionLinks[cell];

    for (int i = cell; i < end; ++i)
        result.push_back(currentPartition[i]);

    return result;
}

std::vector<int> YasolGroupSearch::extractLeafOrder() {
    std::vector<int> order;

    int c = 0;
    while (c < numVertices) {
        int e = currentPartitionLinks[c];

        for (int i = c; i < e; ++i)
            order.push_back(currentPartition[i]);

        c = e;
    }

    return order;
}

void YasolGroupSearch::switchToNextNode(int backtrackDepth) {
    int stackEnd = (int)cellCreationStack.size();

    if (backtrackDepth == 0) {
      if (infoout)
	std::cerr << "switchToNextNode depth zero cleanup\n";
      if (!cellCreationStack.empty())
	cleanupBacktrack(0);
      
      cellCreationStack.clear();
      currNodeCertificate.clear();
      nodeStack.clear();
      return;
    }

    nodeStack.resize(backtrackDepth);

    do {
        if (searchBudgetExceeded()) {
	  abortedByBudget = true;
	  nodeStack.clear();
	  break;
	}

        SearchNode& currNode = nodeStack.back();

        backtrack(currNode.stackStart, stackEnd);
        stackEnd = currNode.stackStart;

        if (firstPathDepth > (int)nodeStack.size())
            firstPathDepth = (int)nodeStack.size();

        if (bestPathDepth > (int)nodeStack.size())
            bestPathDepth = (int)nodeStack.size();

        if (firstLeafPrefixLen > currNode.certificateEnd)
            firstLeafPrefixLen = currNode.certificateEnd;

        if (bestLeafPrefixLen > currNode.certificateEnd)
            bestLeafPrefixLen = currNode.certificateEnd;

        currNodeCertificate.resize(currNode.certificateEnd);
	firstLeafPrefixLen =
	  std::min(firstLeafPrefixLen, currNode.certificateEnd);
	
	bestLeafPrefixLen =
	  std::min(bestLeafPrefixLen, currNode.certificateEnd);
	
        if (!determineNextToDistinguish()) {
            nodeStack.pop_back();
            continue;
        }

        cleanupBacktrack(stackEnd);

        int targetCell = currNode.targetCell;

        if (!distinguishVertex(targetCell)) {
            nodeStack.pop_back();
            continue;
        }

        if (!partitionRefinement()) {
            stackEnd = (int)cellCreationStack.size();
            continue;
        }

        createNode();
        break;

    } while (!nodeStack.empty());
}

void YasolGroupSearch::backtrack(
    int backtrackStackNewEnd,
    int backtrackStackEnd
) {
    for (int stackPos = backtrackStackEnd - 1;
         stackPos >= backtrackStackNewEnd;
         --stackPos) {

        int cell = cellCreationStack[stackPos];

        int newStart = getCellStart(cell - 1);
        int currEnd = currentPartitionLinks[cell];

        currentPartitionLinks[cell] = newStart;
        currentPartitionLinks[newStart] = currEnd;
    }
}

void YasolGroupSearch::cleanupBacktrack(int cellCreationStackPos) {
    for (int stackPos = (int)cellCreationStack.size() - 1;
         stackPos >= cellCreationStackPos;
         --stackPos) {

        int cell = cellCreationStack[stackPos];

        int cellStart = getCellStart(cell);
        int cellEnd = currentPartitionLinks[cellStart];

        for (int v = cell;
             v < cellEnd && vertexToCell[currentPartition[v]] == cell;
             ++v) {
            updateCellMembership(v, cellStart, false);
        }
    }

    cellCreationStack.resize(cellCreationStackPos);
}

int YasolGroupSearch::getCellStart(int pos) {
    int startPos = currentPartitionLinks[pos];

    if (startPos > pos)
        return pos;

    if (currentPartitionLinks[startPos] < startPos) {
        do {
            linkCompressionStack.push_back(pos);
            pos = startPos;
            startPos = currentPartitionLinks[startPos];
        } while (currentPartitionLinks[startPos] < startPos);

        do {
            currentPartitionLinks[linkCompressionStack.back()] = startPos;
            linkCompressionStack.pop_back();
        } while (!linkCompressionStack.empty());
    }

    return startPos;
}

void YasolGroupSearch::processCurrentLeaf() {
    std::vector<int> leafOrder = extractLeafOrder();
    ++numLeaves;

    if (firstLeafCertificate.empty()) {
        firstLeafCertificate = currNodeCertificate;
        bestLeafCertificate = currNodeCertificate;

        firstLeafOrder = leafOrder;
        bestLeafOrder = leafOrder;

        firstLeafGraph = dumpCurrentGraph();
        bestLeafGraph = firstLeafGraph;

        firstLeafPrefixLen = (int)currNodeCertificate.size();
        bestLeafPrefixLen = (int)currNodeCertificate.size();

        firstPathDepth = (int)nodeStack.size();
        bestPathDepth = (int)nodeStack.size();

        if (infoout) std::cerr
            << "first leaf cert size="
            << currNodeCertificate.size()
            << " graphSize="
            << firstLeafGraph.size()
            << "\n";

        return;
    }

// leaf-first wie bisher immer
addGeneratorFromLeafOrder(firstLeafOrder, leafOrder, "leaf-first");

// leaf-best nur, wenn es eine echte andere Referenz ist
// und der Zertifikatspfad plausibel ist
if (!bestLeafOrder.empty() &&
    bestLeafOrder.size() == leafOrder.size() &&
    bestLeafOrder != firstLeafOrder &&
    currNodeCertificate.size() >= bestLeafPrefixLen) {
    addGeneratorFromLeafOrder(bestLeafOrder, leafOrder, "leaf-best");
}

    // Best leaf aktualisieren.
    if (currNodeCertificate < bestLeafCertificate) {
        bestLeafCertificate = currNodeCertificate;
        bestLeafOrder = leafOrder;
        bestLeafPrefixLen = (int)currNodeCertificate.size();
        bestPathDepth = (int)nodeStack.size();
        bestLeafGraph = dumpCurrentGraph();

        if (infoout)
            std::cerr << "better leaf found\n";
    }
}

uint32_t YasolGroupSearch::getVertexHash(int vertex) const {
    auto it = vertexHash.find(vertex);

    if (it == vertexHash.end())
        return 0u;

    return it->second;
}

std::vector<uint64_t> YasolGroupSearch::dumpCurrentGraph() const {
    std::vector<uint64_t> graph;

    for (int v = 0; v < numCol; ++v) {
        int fromCell = vertexToCell[v];

        for (int e = Gstart[v]; e < Gend[v]; ++e) {
            int w = Gedge[e].first;
            uint32_t color = Gedge[e].second;

            int toCell = vertexToCell[w];

            graph.push_back(
                encodeGraphTriple(fromCell, toCell, color)
            );
        }

        for (int e = Gend[v]; e < Gstart[v + 1]; ++e) {
            int w = Gedge[e].first;
            uint32_t color = Gedge[e].second;

            int toCell = w; // bereits komprimierter Zellvertreter/Fixpunkt

            graph.push_back(
                encodeGraphTriple(fromCell, toCell, color)
            );
        }
    }

    std::sort(graph.begin(), graph.end());

    return graph;
}

bool
YasolGroupSearch::compareCurrentGraph(
    const std::vector<uint64_t>& other,
    int& wrongCell
) const {

    std::vector<uint64_t> current =
        dumpCurrentGraph();

    wrongCell = -1;

    static int printedMismatches = 0;
    const int maxPrintedMismatches = 20;

    if (current.size() != other.size()) {
        wrongCell = -2;

        if (printedMismatches < maxPrintedMismatches) {
            ++printedMismatches;

            if (infoout) std::cerr
                << "compareCurrentGraph size mismatch:"
                << " current.size=" << current.size()
                << " other.size=" << other.size()
                << "\n";
        }

        return false;
    }

    for (size_t i = 0;
         i < current.size();
         ++i) {

        if (current[i] != other[i]) {
            wrongCell = (int)i;

            if (printedMismatches < maxPrintedMismatches) {
                ++printedMismatches;

                if (infoout) std::cerr
                    << "compareCurrentGraph mismatch at dump index "
                    << i
                    << " current=" << current[i]
                    << " other=" << other[i]
                    << "\n";
            }

            return false;
        }
    }

    return true;
}

uint64_t YasolGroupSearch::encodeGraphTriple(
    int fromCell,
    int toCell,
    uint32_t edgeColor
) const {
    uint64_t a = (uint64_t)(uint32_t)fromCell;
    uint64_t b = (uint64_t)(uint32_t)toCell;
    uint64_t c = (uint64_t)edgeColor;

    uint64_t h = 1469598103934665603ULL;

    h ^= a;
    h *= 1099511628211ULL;

    h ^= b + 0x9e3779b97f4a7c15ULL;
    h *= 1099511628211ULL;

    h ^= c + 0xc2b2ae3d27d4eb4fULL;
    h *= 1099511628211ULL;

    return h;
}

bool YasolGroupSearch::addVerifiedLeafPermutation(
    std::vector<int>&& perm,
    const char* source
) {
    if ((int)perm.size() != (int)vertexGroundSet.size())
        return false;

    bool nontrivial = false;

    for (int i = 0; i < (int)perm.size(); ++i) {
        const int v = vertexGroundSet[i];

        if (perm[i] != v) {
            nontrivial = true;
            break;
        }
    }

    if (!nontrivial)
        return false;

    for (const auto& oldPerm : permutations) {
        if (oldPerm == perm)
            return false;
    }

    if (numStoredAutomorphisms < maxStoredAutomorphisms) {
        storedAutomorphisms.insert(
            storedAutomorphisms.end(),
            perm.begin(),
            perm.end()
        );

        ++numStoredAutomorphisms;
    }

    permutations.push_back(std::move(perm));
    ++numGenerated;

    return true;
}


bool YasolGroupSearch::addPermutationGenerator(
    std::vector<int>&& perm,
    const char* source
) {
    const int nFull = fullNumVertices();

    if ((int)perm.size() != (int)vertexGroundSet.size())
        return false;

    bool nontrivial = false;

    int movedTotal = 0;
    int movedCols = 0;
    int movedRows = 0;

    for (int i = 0; i < (int)perm.size(); ++i) {
      int v = vertexGroundSet[i];

      if (perm[i] == v)
        continue;

      nontrivial = true;
      ++movedTotal;

      if (v < numCol)
        ++movedCols;
      else
        ++movedRows;
    }

    if (!nontrivial)
      return false;
 
if (validateGenerators && !isPermutationAutomorphism(perm)) {
    ++numGraphMismatches;

    return false;
}
 if (0&&validateGenerators && !isPermutationAutomorphism(perm)) {
        ++numGraphMismatches;

        if (infoout) {
            std::cerr << "rejected non-automorphism generator from "
                      << source << ":";
            for (int p = 0; p < (int)perm.size(); ++p) {
                int v = vertexGroundSet[p];
                int w = perm[p];
                if (v != w)
                    std::cerr << " " << v << "->" << w;
            }
            std::cerr << "\n";
        }

        return false;
 }

    for (const auto& oldPerm : permutations) {
        if (oldPerm == perm)
            return false;
    }

    if (infoout) {
        std::cerr << "accepted automorphism generator from "
                  << source << ":";
        for (int p = 0; p < (int)perm.size(); ++p) {
            int v = vertexGroundSet[p];
            int w = perm[p];
            if (v != w)
	      std::cerr << " " << v << "->" << w;
        }
        std::cerr << "\n";
	if (infoout) {
	  int movedTotal = 0;
	  int movedCols = 0;
	  int movedRows = 0;

	  for (int i = 0; i < (int)perm.size(); ++i) {
	    if (perm[i] == vertexGroundSet[i])
	      continue;

	    ++movedTotal;

	    int v = vertexGroundSet[i];

	    if (v < numCol)
	      ++movedCols;
	    else
	      ++movedRows;
	  }

	  std::cerr << "accepted generator: movedTotal=" << movedTotal
		    << " movedCols=" << movedCols
		    << " movedRows=" << movedRows
		    << "\n";
	}
    }

    if (numStoredAutomorphisms < maxStoredAutomorphisms) {
      storedAutomorphisms.insert(
				 storedAutomorphisms.end(),
				 perm.begin(),
				 perm.end()
				 );
      ++numStoredAutomorphisms;
    }
    
    permutations.push_back(std::move(perm));
    ++numGenerated;

    return true;
}

bool YasolGroupSearch::isPermutationAutomorphism(
    const std::vector<int>& perm
) const {
    const int nFull = fullNumVertices();

    if ((int)perm.size() != (int)vertexGroundSet.size())
        return false;

    if ((int)Gstart.size() < nFull + 1)
        return false;

    std::vector<int> image(nFull);

    for (int v = 0; v < nFull; ++v)
        image[v] = v;

    std::vector<unsigned char> used(nFull, 0);

    for (int p = 0; p < (int)vertexGroundSet.size(); ++p) {
        int v = vertexGroundSet[p];
        int w = perm[p];

        if (v < 0 || v >= nFull)
            return false;

        if (w < 0 || w >= nFull)
            return false;

        int wp = vertexPosition[w];

        if (wp < 0 || wp >= (int)vertexGroundSet.size())
            return false;

        if (used[w])
            return false;

        used[w] = 1;
        image[v] = w;
    }

    auto edgeCode = [](int a, int b, uint32_t color) -> uint64_t {
        uint64_t x = 1469598103934665603ULL;

        x ^= (uint64_t)a + 0x9e3779b97f4a7c15ULL + (x << 6) + (x >> 2);
        x *= 1099511628211ULL;

        x ^= (uint64_t)b + 0x85ebca6b27d4eb2fULL + (x << 6) + (x >> 2);
        x *= 1099511628211ULL;

        x ^= (uint64_t)color + 0xc2b2ae3d27d4eb4fULL + (x << 6) + (x >> 2);
        x *= 1099511628211ULL;

        return x;
    };

    std::vector<uint64_t> original;
    std::vector<uint64_t> mapped;

    original.reserve(Gedge.size());
    mapped.reserve(Gedge.size());

    for (int v = 0; v < nFull; ++v) {
        for (int e = Gstart[v]; e < Gstart[v + 1]; ++e) {
            int w = Gedge[e].first;
            uint32_t color = Gedge[e].second;

            if (w < 0 || w >= nFull)
                return false;

            original.push_back(edgeCode(v, w, color));
            mapped.push_back(edgeCode(image[v], image[w], color));
        }
    }

    std::sort(original.begin(), original.end());
    std::sort(mapped.begin(), mapped.end());

    return original == mapped;
}

void YasolGroupSearch::printSearchStats(const char* where) const {
    const int fullNumVertices = numCol + numRow;

    if (infoout)  std::cerr
        << "\n==== GROUP SEARCH STATS: " << where << " ====\n"
        << "numCol             = " << numCol << "\n"
        << "numRow             = " << numRow << "\n"
        << "full graph vertices= " << fullNumVertices << "\n"
        << "active vertices    = " << numVertices << "\n"
        << "numEdges           = " << Gedge.size() << "\n"
        << "groundSet size     = " << vertexGroundSet.size() << "\n"
        << "nodes              = " << numSearchNodes << "\n"
        << "leaves             = " << numLeaves << "\n"
        << "pruned             = " << numPruned << "\n"
        << "graphMismatch      = " << numGraphMismatches << "\n"
        << "generators         = " << numGenerated << "\n"
        << "all permutations   = " << permutations.size() << "\n"
        << "stored automorphs  = " << numStoredAutomorphisms << "\n"
        << "maxStoredAutomorphs= " << maxStoredAutomorphisms << "\n"
        << "autoSkipped        = " << numAutoSkipped << "\n"
	<< "depthPruned        = " << numDepthPruned << "\n"
        << "abortedByBudget    = " << abortedByBudget << "\n"
        << "maxSearchNodes     = " << maxSearchNodes << "\n"
        << "maxLeaves          = " << maxLeaves << "\n"
        << "maxDepth           = " << maxDepth << "\n"
        << "maxGenerators      = " << maxGenerators << "\n"
	<< "===============================\n";
}

void YasolGroupSearch::printOrbitCells() const {
    if (!infoout) return; 
    std::cerr << "\n==== NONTRIVIAL CURRENT CELLS ====\n";

    int nontrivial = 0;

    for (int c = 0; c < numVertices;) {
        int end = currentPartitionLinks[c];
        int size = end - c;

        if (size > 1) {
            ++nontrivial;
	    if (infoout) {
	      std::cerr << "cell " << c << " : ";
	      
	      for (int i = c; i < end; ++i) {
                int v = currentPartition[i];
		
                if (v < numCol)
		  std::cerr << "x" << (v + 1) << " ";
                else
		  std::cerr << "r" << (v - numCol) << " ";
	      }
	      
	      std::cerr << "\n";
	    }
	}

        c = end;
    }

    std::cerr << "nontrivial cells = " << nontrivial << "\n";
    std::cerr << "===============================\n";
}

int YasolGroupSearch::fullNumVertices() const {
    return numCol + numRow;
}

bool YasolGroupSearch::addGeneratorFromLeafOrder(
    const std::vector<int>& referenceOrder,
    const std::vector<int>& leafOrder,
    const char* source
) {
    if (referenceOrder.size() != leafOrder.size()) {
        return false;
    }

    std::vector<int> perm = vertexGroundSet;
    bool nontrivial = false;
    int moved = 0;
    int mapped = 0;
    int badPosition = 0;

    for (int i = 0; i < (int)leafOrder.size(); ++i) {
        int a = referenceOrder[i];
        int b = leafOrder[i];

        if (a == b)
            continue;

        ++moved;

        int pa = -1;

        if (a >= 0 && a < (int)vertexPosition.size())
            pa = vertexPosition[a];

        if (pa >= 0 && pa < (int)perm.size()) {
            perm[pa] = b;
            nontrivial = true;
            ++mapped;
        } else {
            ++badPosition;
        }
    }

    if (!nontrivial) {
        return false;
    }

bool ok = addVerifiedLeafPermutation(std::move(perm), source);
return ok;
//bool ok = addPermutationGenerator(std::move(perm), source);
//return ok;
}

bool YasolGroupSearch::searchBudgetExceeded() const {
    if (maxSearchNodes >= 0 && numSearchNodes >= maxSearchNodes)
        return true;

    if (maxLeaves >= 0 && numLeaves >= maxLeaves)
        return true;

    if (maxGenerators >= 0 && numGenerated >= maxGenerators)
        return true;

    // Achtung:
    // maxDepth ist KEIN globales Suchbudget.
    // Eine erreichte Tiefengrenze soll nur den aktuellen Ast prunen,
    // nicht die komplette Suche abbrechen.
    return false;
}

double YasolGroupSearch::vertexBranchingScore(int vertex) const {
    double score = 0.0;

    if (vertex < 0 || vertex >= (int)Gstart.size())
        return score;

    const int deg = Gend[vertex] - Gstart[vertex];

    // 1. Normaler Grad.
    score += (double)deg;

    // 2. Verschiedene Kantenfarben.
    // Billig implementiert: Farben sammeln, sortieren, unique.
    std::vector<uint32_t> colors;
    colors.reserve(deg);

    int nonSingletonNeighborCells = 0;

    for (int e = Gstart[vertex]; e < Gend[vertex]; ++e) {
        const int nb = Gedge[e].first;
        const uint32_t col = Gedge[e].second;

        colors.push_back(col);

        if (nb >= 0 && nb < (int)vertexToCell.size()) {
            int nbCell = vertexToCell[nb];
            if (nbCell >= 0 && nbCell < numVertices && cellSize(nbCell) > 1)
                ++nonSingletonNeighborCells;
        }
    }

    std::sort(colors.begin(), colors.end());
    colors.erase(std::unique(colors.begin(), colors.end()), colors.end());

    score += 5.0 * (double)colors.size();

    // 3. Nachbarn in nichttrivialen Zellen:
    // solche Vertices koennen beim Individualisieren oft mehr Refinement ausloesen.
    score += 2.0 * (double)nonSingletonNeighborCells;

    // 4. Leichter Bonus fuer Variablen, die in der objective-Zeile haengen.
    // objective ist constraint[0], also Zeilenknoten numCol + 0.
    const int objectiveVertex = numCol;
    for (int e = Gstart[vertex]; e < Gend[vertex]; ++e) {
        if (Gedge[e].first == objectiveVertex) {
            score += 10.0;
            break;
        }
    }

    return score;
}

double YasolGroupSearch::cellBranchingScore(int cell) const {
    const int end = currentPartitionLinks[cell];
    const int size = end - cell;

    if (size <= 1)
        return -std::numeric_limits<double>::infinity();

    double avgVertexScore = 0.0;

    for (int i = cell; i < end; ++i)
        avgVertexScore += vertexBranchingScore(currentPartition[i]);

    avgVertexScore /= (double)size;

    // Wichtig: kleine Zellen bevorzugen.
    // Das reduziert die Branching-Breite.
    //
    // Gleichzeitig sollen strukturell wirksame Zellen bevorzugt werden.
    return 100000.0 / (double)size + avgVertexScore;
}

bool YasolGroupSearch::isFixedVariableVertex(int v) const {
    if (v < 0 || v >= numCol)
        return false;

    if (v >= (int)fixedVertex.size())
        return false;

    return fixedVertex[v] != 0;
}

bool YasolGroupSearch::checkStoredAutomorphism(int vertex) {
    if (numStoredAutomorphisms == 0)
        return true;

    if (vertex < 0 || vertex >= (int)vertexPosition.size())
        return true;

    const int pos = vertexPosition[vertex];

    if (pos < 0 || pos >= (int)vertexGroundSet.size())
        return true;

    const int groundSize = (int)vertexGroundSet.size();

    for (int a = 0; a < numStoredAutomorphisms; ++a) {
        const int* aut =
            storedAutomorphisms.data() + (size_t)a * groundSize;

        bool automorphismUseful = true;

        // Der gespeicherte Automorphismus darf nur benutzt werden,
        // wenn er die bisherigen Distinguished-Entscheidungen
        // des aktuellen Suchpfads fixiert.
        for (int j = (int)nodeStack.size() - 2; j >= firstPathDepth; --j) {
            if (j < 0 || j >= (int)nodeStack.size())
                continue;

            const int fixedVertex = nodeStack[j].lastDistinguished;

            if (fixedVertex < 0)
                continue;

            if (fixedVertex >= (int)vertexPosition.size())
                continue;

            const int fixedPos = vertexPosition[fixedVertex];

            if (fixedPos < 0 || fixedPos >= groundSize)
                continue;

            if (aut[fixedPos] != vertexGroundSet[fixedPos]) {
                automorphismUseful = false;
                break;
            }
        }

        if (!automorphismUseful)
            continue;

        const int image = aut[pos];

        // Orbit-Minimum-Idee: wenn ein im aktuellen Stabilizer zulaessiger
        // gespeicherter Automorphismus diesen Kandidaten auf einen kleineren
        // Kandidaten abbildet, brauchen wir diesen Branch nicht.
        if (image < vertex) {
            ++numAutoSkipped;
            return false;
        }
    }

    return true;
}

std::vector<int> YasolGroupSearch::collectRootTargetCells() const {
    std::vector<int> cells;

    int c = 0;

    while (c < numVertices) {
        int e = currentPartitionLinks[c];
        int size = e - c;

        if (size > 1) {
            bool hasCol = false;

            for (int i = c; i < e; ++i) {
                if (currentPartition[i] < numCol) {
                    hasCol = true;
                    break;
                }
            }

            if (hasCol)
                cells.push_back(c);
        }

        c = e;
    }

    std::sort(
        cells.begin(),
        cells.end(),
        [&](int a, int b) {
            int sa = currentPartitionLinks[a] - a;
            int sb = currentPartitionLinks[b] - b;

            if (sa != sb)
                return sa < sb;

            return a < b;
        }
    );

    return cells;
}

void YasolGroupSearch::searchAutomorphismsMultiStart() {
    std::vector<int> rootCells = collectRootTargetCells();

    if (rootCells.empty()) {
        searchAutomorphisms();
        return;
    }

    PartitionSnapshot rootState = saveState();

    int oldMaxSearchNodes = maxSearchNodes;
    int oldMaxLeaves = maxLeaves;

    const int perRootNodes = 2000;
    const int perRootLeaves = 2000;

    int rootsTried = 0;

    for (int rootCell : rootCells) {
        if (numGenerated >= maxGenerators)
            break;

        if (rootsTried >= 100)
            break;

        restoreState(rootState);

        forcedRootCell = rootCell;

        maxSearchNodes = perRootNodes;
        maxLeaves = perRootLeaves;

        numSearchNodes = 0;
        numLeaves = 0;
        numPruned = 0;
        abortedByBudget = false;

        searchAutomorphisms();

        ++rootsTried;

        if (infoout) std::cerr
            << "multistart rootCell=" << rootCell
            << " rootsTried=" << rootsTried
            << " generators=" << numGenerated
            << " storedPerms=" << permutations.size()
            << "\n";
    }

    forcedRootCell = -1;

    maxSearchNodes = oldMaxSearchNodes;
    maxLeaves = oldMaxLeaves;

    restoreState(rootState);
}

YasolGroupSearch::PartitionSnapshot
YasolGroupSearch::saveState() const {
    PartitionSnapshot s;

    s.numVertices = numVertices;
    s.numActiveCols = numActiveCols;

    s.currentPartition = currentPartition;
    s.currentPartitionLinks = currentPartitionLinks;
    s.vertexToCell = vertexToCell;
    s.vertexPosition = vertexPosition;
    s.vertexGroundSet = vertexGroundSet;

    s.cellCreationStack = cellCreationStack;
    s.cellInRefinementQueue = cellInRefinementQueue;
    s.refinementQueue = refinementQueue;
    s.currNodeCertificate = currNodeCertificate;

    s.firstLeafPrefixLen = firstLeafPrefixLen;
    s.bestLeafPrefixLen = bestLeafPrefixLen;
    s.firstPathDepth = firstPathDepth;
    s.bestPathDepth = bestPathDepth;

    return s;
}

void
YasolGroupSearch::restoreState(const PartitionSnapshot& s) {
    numVertices = s.numVertices;
    numActiveCols = s.numActiveCols;

    currentPartition = s.currentPartition;
    currentPartitionLinks = s.currentPartitionLinks;
    vertexToCell = s.vertexToCell;
    vertexPosition = s.vertexPosition;
    vertexGroundSet = s.vertexGroundSet;

    cellCreationStack = s.cellCreationStack;
    cellInRefinementQueue = s.cellInRefinementQueue;
    refinementQueue = s.refinementQueue;
    currNodeCertificate = s.currNodeCertificate;

    firstLeafPrefixLen = s.firstLeafPrefixLen;
    bestLeafPrefixLen = s.bestLeafPrefixLen;
    firstPathDepth = s.firstPathDepth;
    bestPathDepth = s.bestPathDepth;

    nodeStack.clear();
    distinguishCands.clear();
    linkCompressionStack.clear();
    vertexHash.clear();
}

void YasolSymmetry::mergeOrbitResults(
    YasolSymmetryResult& target,
    const YasolSymmetryResult& source,
    int nVars
) {
    std::vector<int> parent(nVars);
    std::vector<int> size(nVars, 1);

    for (int i = 0; i < nVars; ++i)
        parent[i] = i;

    auto find = [&](int x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    auto unite = [&](int a, int b) {
        int ra = find(a);
        int rb = find(b);

        if (ra == rb)
            return;

        if (size[ra] < size[rb])
            std::swap(ra, rb);

        parent[rb] = ra;
        size[ra] += size[rb];
    };

    auto addOrbit = [&](const std::vector<int>& orb) {
        if (orb.size() < 2)
            return;

        int first = -1;

        for (int v : orb) {
            if (v < 0 || v >= nVars)
                continue;

            if (first < 0)
                first = v;
            else
                unite(first, v);
        }
    };

    // alte Orbits behalten
    for (const auto& orb : target.orbits)
        addOrbit(orb);

    // neue Orbits dazumergen
    for (const auto& orb : source.orbits)
        addOrbit(orb);

    // target neu aufbauen
    std::map<int, std::vector<int>> classes;

    for (int v = 0; v < nVars; ++v) {
        int r = find(v);
        if (r != v || size[r] > 1) {
            classes[r].push_back(v);
        }
    }

    target.orbits.clear();

    for (auto& it : classes) {
        if (it.second.size() >= 2)
            target.orbits.push_back(std::move(it.second));
    }

    // orbitIdOfVar neu setzen
    target.orbitIdOfVar.assign(nVars, -1);
    target.posInOrbitOfVar.assign(nVars, -1);
 
    for (int o = 0; o < (int)target.orbits.size(); ++o) {
      for (int p = 0; p < (int)target.orbits[o].size(); ++p) {
        int v = target.orbits[o][p];

        if (v >= 0 && v < nVars) {
            target.orbitIdOfVar[v] = o;
            target.posInOrbitOfVar[v] = p;
        }
      }
    }

    target.success = target.success || source.success;
    target.partial = target.partial || source.partial;
}


#if 0
//Debugging
std::vector<YasolGroupSearch::DebugGraphEntry>
YasolGroupSearch::dumpCurrentGraphDebug() const {
    std::vector<DebugGraphEntry> out;

    for (int v = 0; v < fullNumVertices; ++v) {
        int fromCell = vertexToCell[v];

        for (int e = Gstart[v]; e < Gstart[v + 1]; ++e) {
            int w = Gedge[e].first;
            uint32_t color = Gedge[e].second;

            int toCell = vertexToCell[w];

            uint64_t packed = packGraphEntry(fromCell, toCell, color);

            out.push_back({fromCell, toCell, color, packed});
        }
    }

    std::sort(out.begin(), out.end(),
              [](const DebugGraphEntry& a,
                 const DebugGraphEntry& b) {
                  return a.packed < b.packed;
              });

    return out;
}
#endif
