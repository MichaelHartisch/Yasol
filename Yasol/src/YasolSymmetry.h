#ifndef YASOL_SYMMETRY
#define YASOL_SYMMETRY

#include <vector>
#include <map>
#include <cmath>
#include <unordered_map>
#include <cstdint>
#include <algorithm>
#include <limits>

class QBPSolver;

struct YasolSymmetryFixing {
    int var;
    int value; // 0 oder 1
};

struct YasolSymmetryResult {
    bool success = false;
    bool partial = false;

    // Nichttriviale Variablenorbits.
    // Die Einträge sind Original-Variablenindizes.
    std::vector<std::vector<int>> orbits;

    // orbitIdOfVar[v] == -1 bedeutet: Variable v liegt in keinem nichttrivialen Orbit.
    std::vector<int> orbitIdOfVar;
    std::vector<int> posInOrbitOfVar;
public:
    int getOrbitPos(int var) const {
      if (var < 0 || var >= (int)posInOrbitOfVar.size())
	return -1;      
      return posInOrbitOfVar[var];
    }
    int getOrbitID(int var) const {
        if (var < 0 || var >= (int)orbitIdOfVar.size())
            return -1;

        return orbitIdOfVar[var];
    }

    bool hasOrbit(int var) const {
        return getOrbitID(var) >= 0;
    }
};

class YasolGroupSearch {
public:
  struct PartitionSnapshot {
    int numVertices;
    int numActiveCols;

    std::vector<int> currentPartition;
    std::vector<int> currentPartitionLinks;
    std::vector<int> vertexToCell;
    std::vector<int> vertexPosition;
    std::vector<int> vertexGroundSet;

    std::vector<int> cellCreationStack;
    std::vector<unsigned char> cellInRefinementQueue;
    std::vector<int> refinementQueue;
    std::vector<uint32_t> currNodeCertificate;

    std::vector<int> nodeStackDummy; // nicht benutzen, nur falls du später erweitern willst

    int firstLeafPrefixLen;
    int bestLeafPrefixLen;
    int firstPathDepth;
    int bestPathDepth;
};

PartitionSnapshot saveState() const;
void restoreState(const PartitionSnapshot& s);
  
    static inline uint32_t mixHash32(uint32_t h, uint32_t a, uint32_t b) {
      uint32_t x = h;

      x ^= a + 0x9e3779b9u + (x << 6) + (x >> 2);
      x *= 0x85ebca6bu;
      
      x ^= b + 0xc2b2ae35u + (x << 6) + (x >> 2);
      x *= 0x27d4eb2du;
      
      x ^= x >> 15;
      x *= 0x846ca68bu;
      x ^= x >> 16;
      
      return x;
    }
  std::vector<uint32_t> firstLeafCertificate;
  std::vector<uint32_t> bestLeafCertificate;

  std::vector<int> firstLeafOrder;
  std::vector<int> bestLeafOrder;
  std::vector<int> linkCompressionStack;

  int firstLeafPrefixLen = 0;
  int bestLeafPrefixLen = 0;
  int firstPathDepth = 0;
  int bestPathDepth = 0;

  std::vector<int> storedAutomorphisms;
  int numStoredAutomorphisms = 0;
  int maxStoredAutomorphisms = 64;

  void loadFromQbp(QBPSolver& qbp, double eps);
  bool initializeDetection();
  void run(YasolSymmetryResult& result);
  int cellSize(int cell) const;
  void markCellForRefinement(int cell);
  void initializeHashValues();
  bool partitionRefinement();
  void removeFixPoints();

  void initializeGroundSet();
  void extractTrivialGenerators();
  void extractOrbits(YasolSymmetryResult& result);

  void createNode();
  int selectTargetCell();
  bool determineNextToDistinguish();
  bool splitCell(int cell, int splitPoint);
  bool updateCellMembership(int pos, int cell, bool markForRefinement = true);
  bool distinguishVertex(int targetCell);
  uint32_t getVertexHash(int vertex) const;
  
  //bool distinguishSpecificVertex(int targetCell, int vertex);
  std::vector<int> getCellVertices(int cell) const;
  std::vector<std::vector<int>> savePartition();
  void restorePartition(const std::vector<std::vector<int>>& p);
  void processCurrentLeaf();
  void searchAutomorphisms();
  void searchAutomorphismsMultiStart();
  // //void searchAutomorphismsRecursive(std::vector<int>& path,int depth);
  void descendWithDistinguishedVertex(int targetCell, int v);
  //void handleLeaf();
  std::vector<int> extractLeafOrder();
  std::vector<int> computeCertificate();
  std::vector<int> extractLeafPermutation();

  int getCellStart(int pos);

  void backtrack(int backtrackStackNewEnd, int backtrackStackEnd);
  void cleanupBacktrack(int cellCreationStackPos);
  void switchToNextNode(int backtrackDepth);
  
  //PartitionSnapshot saveState() const;
  //void restoreState(const PartitionSnapshot& s);
  
  void debugOneLevelSearch();
  //private:
  std::map<std::vector<long long>, int> coeffColorMap;
  int nextCoeffColor = 1;

  long long encodeDouble(double x) const;

  int makeCoeffColor(double coeff);

  int numCol = 0;
  int numRow = 0;
  int numVertices = 0;
  int numActiveCols = 0;
      
  std::vector<int> Gstart;
  std::vector<int> Gend;
  std::vector<std::pair<int, uint32_t>> Gedge;
      
  std::vector<int> currentPartition;
  std::vector<int> currentPartitionLinks;
  std::vector<int> vertexToCell;
  std::vector<int> vertexPosition;
  std::vector<int> vertexGroundSet;
  std::vector<std::vector<int>> permutations;
      
  std::vector<int> orbitPartition;
  std::vector<int> orbitSize;

  std::unordered_map<int,uint32_t> vertexHash;

  std::vector<uint8_t> cellInRefinementQueue;

  std::vector<int> refinementQueue;

  struct SearchNode {
    int stackStart;
    int certificateEnd;
    int targetCell;
    int lastDistinguished;

    // Geordnete Kandidatenliste fuer die aktuelle Target-Zelle.
    // Die Eintraege sind Original-Vertex-IDs aus currentPartition.
    std::vector<int> branchCandidates;
    int nextBranchCandidate = 0;
  };

  std::vector<SearchNode> nodeStack;
  std::vector<int> bestCertificate;
  std::vector<int> firstLeafPermutation;
  std::vector<std::vector<int>> generators;
  std::vector<int> cellCreationStack;
  std::vector<uint32_t> currNodeCertificate;
  std::vector<int*> distinguishCands;
  std::vector<char> fixedVertex;

  int numSearchNodes = 0;
  int numLeaves = 0;
  int numPruned = 0;
  int numGenerated = 0;
  int numGraphMismatches = 0;
  int numAutoSkipped = 0;

  // Harte Sicherheitsbudgets.
  // Negative Werte bedeuten: kein Limit.
  int maxSearchNodes = 100000; //20000
  int maxLeaves = 10000; // 5000
  int maxDepth = 100; //50
  int maxGenerators = 1000; //500
  int forcedRootCell = -1;
  int localNodeLimit = 2000;

  
  bool abortedByBudget = false;
  int numDepthPruned = 0;
  
  std::vector<uint64_t> firstLeafGraph;
  std::vector<uint64_t> bestLeafGraph;
  
  std::vector<uint64_t> dumpCurrentGraph() const;
  
  bool compareCurrentGraph( const std::vector<uint64_t>& other, int &wrongDumpIndex) const;
  uint64_t encodeGraphTriple(int fromCell, int toCell, uint32_t edgeColor) const;
  bool checkStoredAutomorphism(int vertex);
  bool infoout = false;
  bool validateGenerators = true;
  
  bool addPermutationGenerator( std::vector<int>&& perm, const char* source );
  bool addVerifiedLeafPermutation(std::vector<int>&& perm, const char* source);
  int getColumnOrbitRoot(int col);
  bool mergeColumnOrbits(int col1, int col2);
  bool addGeneratorFromLeafOrder( const std::vector<int>& referenceOrder, const std::vector<int>& leafOrder, const char* source
				  );
  bool isPermutationAutomorphism(const std::vector<int>& perm ) const;
  int fullNumVertices() const;
  
  void printSearchStats(const char* where) const;
  void printOrbitCells() const;

  bool searchBudgetExceeded() const;
  
  double cellBranchingScore(int cell) const;
  double vertexBranchingScore(int vertex) const;
  
  void initializeBranchCandidates(SearchNode& node);
  int findVertexPositionInCell(int cell, int vertex) const;
  
  bool isFixedVariableVertex(int v) const;

  //for Debugging
  struct DebugGraphEntry {
    int fromCell;
    int toCell;
    uint32_t color;
    uint64_t packed;
  };
  std::vector<int> collectRootTargetCells() const;
  std::vector<YasolGroupSearch::DebugGraphEntry> dumpCurrentGraphDebug() const;
};

class YasolSymmetry {
public:
    YasolSymmetryResult run(QBPSolver& qbp);

private:
    YasolGroupSearch groupSearch;

    struct Node {
        int color;
        std::vector<int> adj;
    };
    struct Partition {
        std::vector<std::vector<int>> cells;
    };

    std::vector<Node> graph;
    std::vector<int> varNode;
    std::vector<int> rowNode;

  //std::vector<std::vector<int>> permutations;

    Partition symPartition;
    YasolSymmetryResult persistentResult;
    int persistentNumVars = -1;
  

    std::map<std::vector<long long>, int> colorMap;
    int nextColor = 1;

    struct Signature {
      int ownColor;
      int degree;
      uint64_t h1;
      uint64_t h2;
      
      bool operator==(const Signature& o) const {
	return ownColor == o.ownColor &&
	  degree == o.degree &&
	  h1 == o.h1 &&
	  h2 == o.h2;
      }
    };

    struct SignatureHash {
      std::size_t operator()(const Signature& s) const {
	uint64_t h = 1469598103934665603ULL;
	h ^= (uint64_t)s.ownColor + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
	h ^= (uint64_t)s.degree   + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
	h ^= s.h1                 + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
	h ^= s.h2                 + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
	return (std::size_t)h;
      }
    };

    uint64_t mix64(uint64_t x) const;
    Signature computeSignature(int node);

    void clear();

    int addNode(int color);
    void addEdge(int a, int b);

    long long encodeDouble(double x);
    int getColorId(const std::vector<long long>& key);

    int makeVarColor(
        int type,
        int block,
        double lb,
        double ub,
        double obj
    );

    int makeRowColor(
        double rhs,
        int nnz,
	double maxCoeff
    );

    int makeCoeffColor(double coeff);

    void loadFromQbp(QBPSolver& qbp, double eps);
    void buildColoredGraph(QBPSolver& qbp);
    void detectAutomorphisms();

    void initializePartition();
    void refinePartition();
    void extractGenerators();

    void extractOrbits(YasolSymmetryResult& result);
    void mergeOrbitResults(YasolSymmetryResult& target, const YasolSymmetryResult& source, int nVars);
};

struct Partition {
    std::vector<int> color;
    std::vector<std::vector<int>> cells;
};

#endif
