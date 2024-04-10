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

#ifndef MCTSCLASS_H
#define MCTSCLASS_H

#define deb_o  0
#include "../graphClass/multiset_plus.h"
#include "../graphClass/graphClass.h"
#include <string>
#include<iomanip> 
//#include<conio.h> 
#include <math.h>
#include <algorithm>

//class GraphManager {
  //void GraphManagerInitGraph(int n, bool lout) ;
  //int getNodeCnt();
  //bool AddEdge(int l1, int l2);
  //bool DeleteEdge(int l1, int l2);
  //bool getEdgeByIndex(int e, int &i, int &j);
  //int getAdjacent(int idx);
  //int FirstAdjacentInGraph(int i);
  //int NextAdjacentInGraph(int j);
  //void printGraph();
  //bool checkConsistency();
  //void initTestGraph(int n, bool lout);
//};
#define EXIST 0
#define UNIV  1

class Node {
 public:
  const double EQ_EPS = 1e-10;
  int ID;
  int fatherID;
  double AVscore;
  double ExistScoreSum;
  double UnivScoreSum;
  int64_t ExistScoreCnt;
  int64_t UnivScoreCnt;
  bool innerNode;
  int who2move;
  double minmax_bnd;
  std::pair<int,int> potentialMoves;
  std::vector< std::pair< int,double > > DeltaSolution;
  double lowerBound;
  double upperBound;
  double pseudoScore;
  uint64_t pseudoScoreCnt;
  double Activity;
  bool isFullNode;
  int entryVar;
  int entryVal;
  int visits;
  bool isClosed;
  double LpVal;
  bool LPexists;
  bool isFeasible;

  bool gotEvalFromSucc;
  bool winnerIsExist;

  Node(int i, std::vector<double> &sol, std::vector<double> &sol_cmp, double lb, double ub, bool bndValid, double psScore, double Act, bool pseudoValid, bool activityValid, bool solutionValid, uint64_t & psCnt, int eVar, int eVal, int fID, int w2m, int i1, int i2, double mmbd) {
    if (!isFullNode) isFullNode = solutionValid;
    ID = i;
    fatherID = fID;
    entryVar = eVar;
    entryVal = eVal;
    ExistScoreSum = 0.0;
    UnivScoreSum = 0.0;
    ExistScoreCnt = 1;
    UnivScoreCnt = 1;
    winnerIsExist = true;
    innerNode = false;
    who2move = w2m;
    minmax_bnd = mmbd;
    assert(w2m==0);
    potentialMoves.first = i1;
    potentialMoves.second = i2;
    if (solutionValid) {
      DeltaSolution.clear();
      for (int z=0; z < sol.size();z++) {
	if (sol_cmp.size() != sol.size() || fabs(sol[z] - sol_cmp[z]) > EQ_EPS )
	  DeltaSolution.push_back(std::make_pair(z,sol[z]));
      }
    } else {
    }
    if (activityValid) {
      Activity = Act;
    }
    if (pseudoValid) {
      pseudoScore = psScore;
      pseudoScoreCnt = psCnt;
    }
    if (bndValid) {
      upperBound = ub;
      lowerBound = lb;
    }
    visits = 0;
    isClosed = false;
    LPexists = false;
    gotEvalFromSucc = false;
  }
  Node(int i, int fID, int w2m, int i1, int i2, double mmbd) {
    ID = i;
    fatherID = fID;
    ExistScoreSum = 0.0;
    UnivScoreSum = 0.0;
    ExistScoreCnt = 1;
    UnivScoreCnt = 1;
    winnerIsExist=true;
    DeltaSolution.clear();
    Activity = 0.0;
    pseudoScore = 0.0;
    pseudoScoreCnt = 1.0;
    upperBound = std::numeric_limits<double>::max();
    lowerBound = -std::numeric_limits<double>::max();
    entryVar = -1;
    entryVal = -1;
    visits = 0;
    isClosed = false;
    gotEvalFromSucc = false;
    innerNode=false;
    who2move = w2m;
    minmax_bnd = mmbd;
    potentialMoves.first = i1;
    potentialMoves.second = i2;
    LPexists = false;
  }
  void updateLabel(int var, int val) {
    entryVar = var;
    entryVal = val;
    //cerr << "Node " << ID << " has got y" << var << "=" << val << endl;
  } 
  void updateVisits(int x) {
    visits += x;
    gotEvalFromSucc = false;
    //cerr << "update node " << ID << " has " << visits << endl;
  } 
  void updateUpperBound(double newBnd) {
    upperBound = newBnd;
  }
  void updatePseudoCost(double *pCo, uint64_t *pCoCnt) {
    pseudoScore = *pCo;
    pseudoScoreCnt = *pCoCnt;
  }
  void updateActivity(double *ac) {
    Activity = *ac;
  }
};

class _MCTS {
 public:
  double dont_know;
  double n_infinity;
  double p_infinity;
  double currentSimVal;
  bool currentSimWinnerIsExist;
  bool currentSimValid;
  std::vector<int> blocks;
  std::vector<int> eas;
  std::vector<int> types;
  GraphManager GM;
  std::vector< Node > nodes;
  int maxBlock;
  time_t last_time_check;
  _MCTS() {
    currentSimVal = 0.0;
    currentSimValid=false;
    currentSimWinnerIsExist=true;
    maxBlock = 1;
    last_time_check = time(NULL);
  }
  bool rootExists() {
    if (nodes.size() > 0) return true;
    else return false;
  }
  void AddRoot() {
    int Vto=-1;
    for (int i=0;i < blocks.size();i++)
      if (blocks[i]==1) Vto=i;
      else break;
    assert(Vto>=0);
    Node nd(0,-1,eas[0],0,Vto,dont_know);
    assert(nodes.size()== 0);
    nodes.push_back(nd);
  }

  void updateBlockAndPotenitalMoves(int id, int exaV) {
    int pBlock = blocks[exaV];
    int Vto=-1;
    int Vfro=-1;
    if (nodes[id].who2move >= 0) return;
    for (int i=0;i < blocks.size();i++) {
      if (maxBlock < blocks[i])
	maxBlock = blocks[i];
    }
    for (int i=0;i < blocks.size();i++) {
      if (blocks[i] < pBlock) continue;
      if (Vfro == -1 && blocks[i] == pBlock) Vfro = Vto = i;
      else {
	if (blocks[i]==pBlock) Vto=i;
	else break;
      }
    }
    assert(Vto>=0);
    assert(Vfro>=0);
    assert(Vto>=Vfro);
    nodes[id].who2move = eas[exaV];
    //assert(nodes[id].who2move==0);
    nodes[id].potentialMoves.first = Vfro;
    nodes[id].potentialMoves.second = Vto;
  }

  void initBlocksAndEasAndValues(int *b, int *ea, int *t, int N, double dn, double ninf) {
    dont_know = dn;
    n_infinity = ninf;
    p_infinity = -ninf;
    blocks.clear();
    eas.clear();
    for (int i =0; i < N;i++)
      blocks.push_back(b[i]);
    for (int i =0; i < N;i++) {
      eas.push_back(ea[i]);
      //cerr << eas[i] << endl;
    }
    //cerr << endl;
    for (int i =0; i < N;i++)
      types.push_back(t[i]);
  }

  int findSucc(int id, int toVar, int toVal) {
    int j;
    j = GM.FirstAdjacentInGraph(id);
    if (j >= 0) {
      int w = GM.getAdjacent(j);
      //w first succ
      //cerr << "FOUND for ID=" << id << " x" << toVar << "=" << toVal << " :Succ=" << nodes[w].ID << " with y" << nodes[w].entryVar << "=" << nodes[w].entryVal << endl;
      if ( nodes[w].entryVar == toVar && nodes[w].entryVal == toVal) 
	return nodes[w].ID;
      int kk = GM.NextAdjacentInGraph(j);
      while (kk >= 0) {
	int w = GM.getAdjacent(kk);
	//cerr << "FOUND for ID=" << id << " x" << toVar << "=" << toVal << " :Succ=" << nodes[w].ID << " with y" << nodes[w].entryVar << "=" << nodes[w].entryVal << endl;
	// w next succ
	if ( nodes[w].entryVar == toVar && nodes[w].entryVal == toVal) 
	  return nodes[w].ID;
	kk = GM.NextAdjacentInGraph(kk);
      }
    }
    return -1;
  }
  void printNodeInfo(int nodeID, int pick, int pol,
		    double* pseudo0score,
		    double* pseudo1score,
		    uint64_t* pseudo0scoreCnt,
		    uint64_t* pseudo1scoreCnt,
		    double* pActivity,
		    double* nActivity,
		    int8_t *assigns,
		    int8_t *killers,
		    double upperBnd,
		    int n) {
    double score0 = (nodes[0].ExistScoreSum/nodes[0].ExistScoreCnt);

    cerr << "This is node " << nodeID << ", " << (nodes[nodeID].who2move==EXIST?"e":"u")<< " [" << nodes[nodeID].potentialMoves.first << ","<< nodes[nodeID].potentialMoves.second <<  "] and picked is x" << pick << "=" << pol << " mmE=" << nodes[nodeID].minmax_bnd<< endl;
    if (nodes[nodeID].visits > 2 && !nodes[nodeID].innerNode)
      assert(0); 
    //welche sons als Knoten vorhanden; node nr, eas, var, closed, minmaxest, stats;
    int j;
    j = GM.FirstAdjacentInGraph(nodeID);
    if (j >= 0) {
      int w = GM.getAdjacent(j);
      //w first succ


      int brother = w;
      double lossS = -((nodes[w].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[w].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
      if (nodes[w].entryVal == 0) brother++;
      else brother--;
      double loss2 = /*sqrt*/(fabs((pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])
				   * (pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar])));
      double loss3 = sqrt(fabs((pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])
				   * (pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar])));
      double loss = 1.0/((nodes[w].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[w].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
      if (nodes[w].UnivScoreCnt == 0 || nodes[w].UnivScoreSum /* + nodes[brother].UnivScoreSum*/ == 0.0) loss = 0.0;//1.0;
      if (pseudo1scoreCnt[nodes[w].entryVar] == 0 || pseudo0scoreCnt[nodes[w].entryVar] == 0) loss2 = 0;
      if (pseudo1scoreCnt[nodes[w].entryVar] == 0 || pseudo0scoreCnt[nodes[w].entryVar] == 0) loss3 = 0;
      if (nodes[w].UnivScoreCnt == 0 || nodes[w].UnivScoreSum == 0.0) lossS = 0;

      double bloss = 1.0/((nodes[brother].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[brother].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
      if (nodes[brother].UnivScoreCnt == 0 || nodes[brother].UnivScoreSum /* + nodes[brother].UnivScoreSum*/ == 0.0) bloss = 0.0;

      double score = (nodes[nodeID].ExistScoreSum/nodes[nodeID].ExistScoreCnt) * (1.0/nodes[w].ExistScoreCnt) 
	+ (nodes[w].ExistScoreSum/nodes[w].ExistScoreCnt)  * (1.0-1.0/nodes[w].ExistScoreCnt);    
      double score2 = -score + sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
      double score8 = score/nodes[0].visits+/*fabs(score)*loss*/lossS + loss2*(nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
      double score9 = score/nodes[0].visits+/*fabs(score)*loss*/loss + loss3*(nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
      double score10 = score/nodes[0].visits+/*fabs(score)*loss*/loss + (loss2/ (fmax(fabs(score),1.0) )) * (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
      double score11 = score/nodes[0].visits+/*fabs(score)*loss*/loss + (loss2/ (fmax(fabs(score),1.0) )) + (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
      double score12 = score/nodes[0].visits+/*fabs(score)*loss*/loss + 1e5 * (loss2/ (fmax(fabs(score),1.0) )) + (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
      double score13 = score/nodes[0].visits+/*fabs(score)*loss*/loss + 1e5 * (loss2/ (fmax(fabs(score),1.0) )) + (nodes[brother].isClosed?2.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
      double score19 = score/nodes[0].visits+(fabs(score0)+1000.0)*loss*bloss + loss + (nodes[brother].isClosed?200000.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*(0.001+loss));
      double score20 = score/nodes[0].visits+(fabs(score0)+1000.0)*loss*bloss + loss + (nodes[brother].isClosed?1.0:1.0)*(sqrt((double)nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*(0.001+loss));
      
      //double score2 = score*loss + (sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*/*loss*/(1.0/((nodes[w].UnivScoreSum+nodes[brother].UnivScoreSum) / (nodes[w].UnivScoreCnt+nodes[brother].UnivScoreCnt))));
      score10 = score20;

      cerr << "   FOUND x" << nodes[w].entryVar << "=" << nodes[w].entryVal << " :SuccID:" << nodes[w].ID << " closed:" << nodes[w].isClosed << nodes[brother].isClosed <<" mE=" << nodes[w].minmax_bnd << " LP:" << nodes[w].LpVal;
      cerr << setprecision(5) << " sc=" << score10 << " loss=" << loss << " b-loss=" << bloss << " stats:" << nodes[w].ExistScoreSum << " / " << nodes[w].UnivScoreSum << " // " << nodes[w].ExistScoreCnt << " / " << nodes[w].UnivScoreCnt << endl;
      //if ( nodes[w].entryVar == toVar && nodes[w].entryVal == toVal) 
      //return nodes[w].ID;
      int kk = GM.NextAdjacentInGraph(j);
      while (kk >= 0) {
	int w = GM.getAdjacent(kk);
	//cerr << "FOUND for ID=" << id << " x" << toVar << "=" << toVal << " :Succ=" << nodes[w].ID << " with y" << nodes[w].entryVar << "=" << nodes[w].entryVal << endl;
	// w next succ
	int brother = w;
	double lossS = -((nodes[w].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[w].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
	if (nodes[w].entryVal == 0) brother++;
	else brother--;
	double loss2 = /*sqrt*/(fabs((pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])
				     * (pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar])));
	double loss3 = sqrt(fabs((pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])
				     * (pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar])));
	double loss = 1.0/((nodes[w].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[w].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
	if (nodes[w].UnivScoreCnt == 0 || nodes[w].UnivScoreSum /* + nodes[brother].UnivScoreSum*/ == 0.0) loss = 0.0;//1.0;
	if (pseudo1scoreCnt[nodes[w].entryVar] == 0 || pseudo0scoreCnt[nodes[w].entryVar] == 0) loss2 = 0;
	if (pseudo1scoreCnt[nodes[w].entryVar] == 0 || pseudo0scoreCnt[nodes[w].entryVar] == 0) loss3 = 0;
	if (nodes[w].UnivScoreCnt == 0 || nodes[w].UnivScoreSum == 0.0) lossS = 0;

	double bloss = 1.0/((nodes[brother].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[brother].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
	if (nodes[brother].UnivScoreCnt == 0 || nodes[brother].UnivScoreSum /* + nodes[brother].UnivScoreSum*/ == 0.0) bloss = 0.0;


	double score = (nodes[nodeID].ExistScoreSum/nodes[nodeID].ExistScoreCnt) * (1.0/nodes[w].ExistScoreCnt) 
	  + (nodes[w].ExistScoreSum/nodes[w].ExistScoreCnt)  * (1.0-1.0/nodes[w].ExistScoreCnt);    
	double score2 = -score + sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	double score8 = score/nodes[0].visits+/*fabs(score)*loss*/lossS + loss2*(nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	double score9 = score/nodes[0].visits+/*fabs(score)*loss*/loss + loss3*(nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	double score10 = score/nodes[0].visits+/*fabs(score)*loss*/loss + (loss2/ (fmax(fabs(score),1.0) )) * (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	double score11 = score/nodes[0].visits+/*fabs(score)*loss*/loss + (loss2/ (fmax(fabs(score),1.0) )) + (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	double score12 = score/nodes[0].visits+/*fabs(score)*loss*/loss + 1e5 * (loss2/ (fmax(fabs(score),1.0) )) + (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	double score13 = score/nodes[0].visits+/*fabs(score)*loss*/loss + 1e5 * (loss2/ (fmax(fabs(score),1.0) )) + (nodes[brother].isClosed?2.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	double score19 = score/nodes[0].visits+(fabs(score0)+1000.0)*loss*bloss + loss + (nodes[brother].isClosed?200000.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*(0.001+loss));
	double score20 = score/nodes[0].visits+(fabs(score0)+1000.0)*loss*bloss + loss + (nodes[brother].isClosed?1.0:1.0)*(sqrt((double)nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*(0.001+loss));
	score10=score20;
	cerr << "   FOUND x" << nodes[w].entryVar << "=" << nodes[w].entryVal << " :SuccID:" << nodes[w].ID << " closed:" << nodes[w].isClosed << nodes[brother].isClosed << " mE=" << nodes[w].minmax_bnd << " LP:" << nodes[w].LpVal;
	cerr << setprecision(5) << " sc=" << score10 << " loss=" << loss << " b-loss=" << bloss << " stats:" << nodes[w].ExistScoreSum << " / " << nodes[w].UnivScoreSum << " // " << nodes[w].ExistScoreCnt << " / " << nodes[w].UnivScoreCnt << endl;
	//if ( nodes[w].entryVar == toVar && nodes[w].entryVal == toVal) 
	//return nodes[w].ID;
	kk = GM.NextAdjacentInGraph(kk);
      }
    }

  }

#define THESCORE score20 //score20  //ohne "nutz nur wenn auch bestS_pseudo definiert ist"
	  //score2: 46 - 695s // 138 - 8752s, 1F // fast0507 nein, mzzv42 nein
	  //score2: 46 - ???s // 137 - ???s, ? Absturz // fast0507 nein, mzzv42 nein ---- nur !feasPhase Zugsammeln
	  //score3: 
	  //score4: 
	  //score5: 
	  //score6: 
	  //score7: 46 - ???s // 138 - ????s, ?? // fast0507 nein, mzzv42 nein
	  //score7: 46 - 608s // 137 - 7571s, 1 Absturz // fast0507 nein, mzzv42 nein ---- nur !feasPhase Zugsammeln
	  //score8: 
	  //score9: 
	  //score10: 
	  //score11: 
	  //score12: 
	  //score13: 46 - 1438s // 137 - 6105s 1 Absturz // fast0507 nein, mzzv42 nein
	  //score20: 46 - 568 s // 137 - 10134s // fast0507 nein, mzzv42 ja

  void findBestSucc(int id, int &best_succ, int &best_dir, 
		    double* pseudo0score,
		    double* pseudo1score,
		    uint64_t* pseudo0scoreCnt,
		    uint64_t* pseudo1scoreCnt,
		    double* pActivity,
		    double* nActivity,
		    int8_t *assigns,
		    int8_t *killers,
		    double upperBnd,
		    int n) {
    int j;
    int best_w=-1;
    bool doO;
    double best_score = -std::numeric_limits<double>::max();
    best_succ = -1;
    
    if (id < 0 || !nodes[id].innerNode) return;
    //if (id != 0) return;
    if (nodes[id].isClosed) {
      best_succ = -2;
      return;
    }

    //if (nodes[0].visits < 50) return;

    best_dir = 0;
    if (last_time_check+30 < time(NULL)) {
      doO=true;
      last_time_check = time(NULL);
    } else 
      doO = false;

    //doO = true;

    if (id==0 && doO) {
      std::cerr << "MCTS node" << id << (nodes[id].innerNode?" is inner node." : "is a leaf.") << std::endl;
      std::cerr << "moves from" << nodes[id].potentialMoves.first << " to " << nodes[id].potentialMoves.second << " w2m:" << nodes[id].who2move << std::endl;
      std::cerr << "#possible successors:" << 2*(nodes[id].potentialMoves.second-nodes[id].potentialMoves.first) << std::endl;
      std::cerr << "LpVal exists:" << (nodes[id].LPexists ? " yes":" no") << std::endl;
    }

    int bestS_pseudo=-1;
    int bestS_inf=-1;
    int bestS_infII=-1;
    int bestS_ninf=-1;
    int bestS_act=-1;
    int bestS_pseudo_cnt=0;
    int bestS_inf_cnt=0;
    int bestS_infII_cnt=0;
    int bestS_ninf_cnt=0;
    int bestS_all_cnt=0;
    int bestS_latest=-1;
    int sumVisits = 0;
    double bestS_pseudo_score=n_infinity;
    double bestS_inf_score=n_infinity;
    double bestS_infII_score=n_infinity;
    double bestS_ninf_score=n_infinity;
    double bestS_act_score=n_infinity;

    double score0 = (nodes[0].ExistScoreSum/nodes[0].ExistScoreCnt);

    j = GM.FirstAdjacentInGraph(id);
    if (j >= 0) {
      int w = GM.getAdjacent(j);
      //w first succ

      if ( nodes[w].entryVar >= 0 && assigns[nodes[w].entryVar] == 2 && !nodes[w].isClosed) {
	sumVisits += nodes[w].visits;
	if (nodes[w].entryVal == 0) {
	  bestS_all_cnt++;
	  if (pseudo0scoreCnt[nodes[w].entryVar] > 0) bestS_pseudo_cnt++;
	} else if (nodes[w].entryVal == 1) {
	  bestS_all_cnt++;
	  if (pseudo1scoreCnt[nodes[w].entryVar] > 0) bestS_pseudo_cnt++;
	} else assert(0);
	bestS_ninf_cnt += nodes[w].ExistScoreCnt;
	bestS_inf_cnt  += nodes[w].UnivScoreCnt;
	bestS_latest = w;

	if (nodes[w].LPexists && pseudo0scoreCnt[nodes[w].entryVar] > 0 && pseudo1scoreCnt[nodes[w].entryVar] > 0) {
	  //std::cerr << nodes[w].LpVal << "," << pseudo0score[nodes[w].entryVar] << "," << pseudo1score[nodes[w].entryVar] << std::endl;
	  double loss = (/*nodes[w].LpVal -*/ (pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])) 
	    * (/*nodes[w].LpVal -*/ (pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar]));
	  if (loss > bestS_pseudo_score) {
	    bestS_pseudo_score = loss;
	    bestS_pseudo = w;
	  }
	}
	if (eas[nodes[w].entryVar] == EXIST) {
	  double score = (nodes[id].ExistScoreSum/nodes[id].ExistScoreCnt) * (1.0/nodes[w].ExistScoreCnt) 
	    + (nodes[w].ExistScoreSum/nodes[w].ExistScoreCnt)  * (1.0-1.0/nodes[w].ExistScoreCnt);    
	  double loss = ((pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])) 
	              * ((pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar]));
	  if (nodes[w].entryVal == 0) loss = pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar];
	  else loss = pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar];
	  //loss = 1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt);
	  //loss = 1.0;
	  //if (nodes[w].UnivScoreSum == 0) loss = 1.0;
	  double score2 = score + (sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)) * loss; 
	  if (score2 > bestS_ninf/*nodes[w].ExistScoreSum / nodes[w].ExistScoreCnt > bestS_ninf_score && nodes[w].ExistScoreCnt > 1*/) {
	    bestS_ninf_score = score2;//nodes[w].ExistScoreSum / nodes[w].ExistScoreCnt;
	    bestS_ninf = w;
	  }
	} else if (eas[nodes[w].entryVar] == UNIV) {
	  double score = (nodes[id].ExistScoreSum/nodes[id].ExistScoreCnt) * (1.0/nodes[w].ExistScoreCnt) 
	    + (nodes[w].ExistScoreSum/nodes[w].ExistScoreCnt)  * (1.0-1.0/nodes[w].ExistScoreCnt);    
	  double loss = 1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt);
	  if (nodes[w].UnivScoreSum == 0) loss = 1.0;
	  double score2 = -score + (sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*loss);
	  if (score2 > bestS_ninf_score) {
	    //if (bestS_ninf_score != n_infinity) cerr << "take " << w << " because " << score2 << " > " << bestS_ninf_score << endl;
	    bestS_ninf_score = score2;
	    //bestS_ninf_score = -nodes[w].ExistScoreSum / nodes[w].ExistScoreCnt + sqrt(nodes[id].visits) / (nodes[w].ExistScoreCnt+nodes[w].UnivScoreCnt-1);
	    //bestS_ninf_score = -nodes[w].ExistScoreSum / nodes[w].ExistScoreCnt;
	    bestS_ninf = w;
	  }
	}
	if (eas[nodes[w].entryVar] == UNIV) {
	  if (-nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt > bestS_inf_score && nodes[w].UnivScoreCnt > 1) {
	    if (killers[nodes[w].entryVar] != 0 && killers[nodes[w].entryVar] != 1) {
	      bestS_inf_score = 1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt) + 0.5*(sqrt(nodes[0].visits) / (nodes[w].visits+1));
	      bestS_inf = w;
	    } else if (killers[nodes[w].entryVar] == nodes[w].entryVal) {
	      bestS_inf_score = 100+1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt)*(sqrt(nodes[0].visits) / (nodes[w].visits+1));
	      bestS_inf = w;
	    } else {
	      bestS_inf_score = 100+1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt)*(sqrt(nodes[0].visits) / (nodes[w].visits+1));
	      bestS_inf = w;
	    }
	  }
	} else if (eas[nodes[w].entryVar] == EXIST) {
	  int brother = w;
	  double score = (nodes[id].ExistScoreSum/nodes[id].ExistScoreCnt) * (1.0/nodes[w].ExistScoreCnt) 
	    + (nodes[w].ExistScoreSum/nodes[w].ExistScoreCnt)  * (1.0-1.0/nodes[w].ExistScoreCnt);    
	  double loss = 1.0/((nodes[w].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[w].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
	  double lossS = -((nodes[w].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[w].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
	    if (nodes[w].entryVal == 0) brother++;
	    else brother--;
	    double bloss = 1.0/((nodes[brother].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[brother].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
	    if (nodes[brother].UnivScoreCnt == 0 || nodes[brother].UnivScoreSum /* + nodes[brother].UnivScoreSum*/ == 0.0) bloss = 0.0;//1.0;
	    //if (nodes[brother].isClosed) loss = loss * 2.0;
	    //loss = loss / (nodes[brother].UnivScoreSum / nodes[brother].UnivScoreCnt);
	    if (nodes[w].UnivScoreCnt == 0 || nodes[w].UnivScoreSum /* + nodes[brother].UnivScoreSum*/ == 0.0) loss = 0.0;//1.0;
	    if (nodes[w].UnivScoreCnt == 0 || nodes[w].UnivScoreSum == 0.0) lossS = 0;
	    double loss2 = /*sqrt*/(fabs((pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])
				   * (pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar])));
	    double loss3 = sqrt(fabs((pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])
				   * (pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar])));
	    if (pseudo1scoreCnt[nodes[w].entryVar] == 0 || pseudo0scoreCnt[nodes[w].entryVar] == 0) loss2 = 0;
	    if (pseudo1scoreCnt[nodes[w].entryVar] == 0 || pseudo0scoreCnt[nodes[w].entryVar] == 0) loss3 = 0;
	    double score2 = score+fabs(score)*loss + (nodes[brother].isClosed?200000.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*loss);
	    double score3 = /*score+fabs(score)*loss*/loss + loss2 + (nodes[brother].isClosed?2.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*loss);
	    double score4 = score+/*fabs(score)*loss*/loss2 + loss*(nodes[brother].isClosed?2.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score5 = score+/*fabs(score)*loss*/lossS + loss2*(nodes[brother].isClosed?2.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score6 = /*score+*//*fabs(score)*loss*/loss2 + loss*(nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score7 = score/nodes[0].visits+/*fabs(score)*loss*/loss2 + loss*(nodes[brother].isClosed?2.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score8 = score/nodes[0].visits+/*fabs(score)*loss*/lossS + loss2*(nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score9 = score/nodes[0].visits+/*fabs(score)*loss*/loss + loss3*(nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score10 = score/nodes[0].visits+/*fabs(score)*loss*/loss + (loss2/ (fmax(fabs(score),1.0) )) * (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score11 = score/nodes[0].visits+/*fabs(score)*loss*/loss + (loss2/ (fmax(fabs(score),1.0) )) + (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score12 = score/nodes[0].visits+/*fabs(score)*loss*/loss + 1e5 * (loss2/ (fmax(fabs(score),1.0) )) + (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score13 = score/nodes[0].visits+(fabs(score0)+10.0)*loss + 0.0 * (loss2/ (fmax(fabs(score),1.0) )) + (0.001+loss)*(nodes[brother].isClosed?2000000.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score19 = score/nodes[0].visits+(fabs(score0)+1000.0)*loss*bloss + loss + (nodes[brother].isClosed?200000.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*(0.001+loss));
	      double score20 = score/nodes[0].visits+(fabs(score0)+1000.0)*loss*bloss + loss + (nodes[brother].isClosed?1.0:1.0)*(sqrt((double)nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*(0.001+loss));
	    
	    //double score2 = score*loss + (sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*/*loss*/(1.0/((nodes[w].UnivScoreSum+nodes[brother].UnivScoreSum) / (nodes[w].UnivScoreCnt+nodes[brother].UnivScoreCnt))));
	  score2 = THESCORE;
	  if (score2 > bestS_inf_score) {
	    bestS_inf_score = score2;
	    bestS_inf = w;
	  }
	}
	if (pActivity[nodes[w].entryVar] + nActivity[nodes[w].entryVar] > bestS_act_score) {
	  bestS_act_score = pActivity[nodes[w].entryVar] + nActivity[nodes[w].entryVar];
	  bestS_act = w;
	}	
      }
      int kk = GM.NextAdjacentInGraph(j);
      while (kk >= 0) {
	int w = GM.getAdjacent(kk);
	// w next succ
	if ( nodes[w].entryVar >= 0 && assigns[nodes[w].entryVar] == 2  && !nodes[w].isClosed) {
	  sumVisits += nodes[w].visits;
	  if (nodes[w].entryVal == 0) {
	    bestS_all_cnt++;
	    if (pseudo0scoreCnt[nodes[w].entryVar] > 0) bestS_pseudo_cnt++;
	  } else if (nodes[w].entryVal == 1) {
	    bestS_all_cnt++;
	    if (pseudo1scoreCnt[nodes[w].entryVar] > 0) bestS_pseudo_cnt++;
	  } else assert(0);
	  bestS_ninf_cnt += nodes[w].ExistScoreCnt;
	  bestS_inf_cnt  += nodes[w].UnivScoreCnt;
	  bestS_latest = w;
	  if (nodes[w].LPexists && pseudo0scoreCnt[nodes[w].entryVar] > 2 && pseudo1scoreCnt[nodes[w].entryVar] > 2) {
	    //std::cerr << nodes[w].LpVal << "," << pseudo0score[nodes[w].entryVar] << "," << pseudo1score[nodes[w].entryVar] << std::endl;
	    double loss = (/*nodes[w].LpVal -*/ (pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])) 
	      * (/*nodes[w].LpVal -*/ (pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar]));
	    if (loss > bestS_pseudo_score) {
	      bestS_pseudo_score = loss;
	      bestS_pseudo = w;
	    }
	  }
	  if (eas[nodes[w].entryVar] == EXIST) {
	    double score = (nodes[id].ExistScoreSum/nodes[id].ExistScoreCnt) * (1.0/nodes[w].ExistScoreCnt) 
	      + (nodes[w].ExistScoreSum/nodes[w].ExistScoreCnt)  * (1.0-1.0/nodes[w].ExistScoreCnt);    
	    //double score2 = score + sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	    double loss = ((pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])) 
	                * ((pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar]));
	    if (nodes[w].entryVal == 0) loss = pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar];
	    else loss = pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar];
	    //loss = 1.0;
	    //loss = 1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt);
	    //if (nodes[w].UnivScoreSum == 0) loss = 1.0;
	    double score2 = score + (sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)) * loss; 

	    if (score2 > bestS_ninf/*nodes[w].ExistScoreSum / nodes[w].ExistScoreCnt > bestS_ninf_score && nodes[w].ExistScoreCnt > 1*/) {
	      bestS_ninf_score = score2;//nodes[w].ExistScoreSum / nodes[w].ExistScoreCnt;
	      bestS_ninf = w;
	    }
	  } else if (eas[nodes[w].entryVar] == UNIV) {
	    double score = (nodes[id].ExistScoreSum/nodes[id].ExistScoreCnt) * (1.0/nodes[w].ExistScoreCnt) 
	      + (nodes[w].ExistScoreSum/nodes[w].ExistScoreCnt)  * (1.0-1.0/nodes[w].ExistScoreCnt);    
	    double loss = 1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt);
	    if (nodes[w].UnivScoreSum == 0) loss = 1.0;
	    double score2 = -score + (sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*loss);
	    if (score2 > bestS_ninf_score) {
	      bestS_ninf_score = score2;
	      //bestS_ninf_score = -nodes[w].ExistScoreSum / nodes[w].ExistScoreCnt + sqrt(nodes[id].visits) / (nodes[w].ExistScoreCnt+nodes[w].UnivScoreCnt-1);
	      //if (bestS_ninf_score != n_infinity) cerr << "take " << w << " because " << score2 << " > " << bestS_ninf_score << endl;
	      bestS_ninf = w;
	    }
	  }
	  if (eas[nodes[w].entryVar] == UNIV) {
	    if (-nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt > bestS_inf_score && nodes[w].UnivScoreCnt > 1) {
	      if (killers[nodes[w].entryVar] != 0 && killers[nodes[w].entryVar] != 1) {
		bestS_inf_score = 1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt) + 0.5*(sqrt(nodes[0].visits) / (nodes[w].visits+1));
		bestS_inf = w;
	      } else if (killers[nodes[w].entryVar] == nodes[w].entryVal) {
		bestS_inf_score = 100+1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt)*(sqrt(nodes[0].visits) / (nodes[w].visits+1));
		bestS_inf = w;
	      } else {
		bestS_inf_score = 100+1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt)*(sqrt(nodes[0].visits) / (nodes[w].visits+1));
		bestS_inf = w;
	      }
	    }
	  } else if (eas[nodes[w].entryVar] == EXIST) {
	      int brother = w;
	      double score = (nodes[id].ExistScoreSum/nodes[id].ExistScoreCnt) * (1.0/nodes[w].ExistScoreCnt) 
		+ (nodes[w].ExistScoreSum/nodes[w].ExistScoreCnt)  * (1.0-1.0/nodes[w].ExistScoreCnt);    
	      double loss = 1.0/((nodes[w].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[w].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
	      double lossS = -((nodes[w].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[w].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
	      //double loss = 1.0/(nodes[w].UnivScoreSum / nodes[w].UnivScoreCnt);
	      if (nodes[w].entryVal == 0) brother++;
	      else brother--;
	      double bloss = 1.0/((nodes[brother].UnivScoreSum/*+nodes[brother].UnivScoreSum*/) / (nodes[brother].UnivScoreCnt/*+nodes[brother].UnivScoreCnt*/));
	      if (nodes[brother].UnivScoreCnt == 0 || nodes[brother].UnivScoreSum /* + nodes[brother].UnivScoreSum*/ == 0.0) bloss = 0.0;	      //if (nodes[brother].isClosed) loss = loss * 2.0;
	      //loss = loss / (nodes[brother].UnivScoreSum / nodes[brother].UnivScoreCnt);
	      if (nodes[w].UnivScoreCnt == 0 || nodes[w].UnivScoreSum /* + nodes[brother].UnivScoreSum*/ == 0.0) loss = 0.0;//1.0;

	      if (nodes[w].UnivScoreCnt == 0 || nodes[w].UnivScoreSum ==0.0) lossS = 0;
	      //double score2 = score*loss + (sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*/*loss*/(1.0/((nodes[w].UnivScoreSum+nodes[brother].UnivScoreSum) / (nodes[w].UnivScoreCnt+nodes[brother].UnivScoreCnt))));
	      double loss2 = /*sqrt*/(fabs((pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])
				   * (pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar])));
	      double loss3 = sqrt(fabs((pseudo0score[nodes[w].entryVar] / pseudo0scoreCnt[nodes[w].entryVar])
				   * (pseudo1score[nodes[w].entryVar] / pseudo1scoreCnt[nodes[w].entryVar])));
	      if (pseudo1scoreCnt[nodes[w].entryVar] == 0 || pseudo0scoreCnt[nodes[w].entryVar] == 0) loss2 = 0;
	      if (pseudo1scoreCnt[nodes[w].entryVar] == 0 || pseudo0scoreCnt[nodes[w].entryVar] == 0) loss3 = 0;
	      double score2 = score+fabs(score)*loss + (nodes[brother].isClosed?200000.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*loss);
	      double score3 = /*score+fabs(score)*loss*/loss + loss2 + (nodes[brother].isClosed?2.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*loss);
	      double score4 = /*score+fabs(score)*loss*/loss2 + loss*(nodes[brother].isClosed?2.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score5 = score+/*fabs(score)*loss*/lossS + loss2*(nodes[brother].isClosed?2.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score6 = /*score+*//*fabs(score)*loss*/loss2 + loss*(nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score7 = score/nodes[0].visits+/*fabs(score)*loss*/loss2 + loss*(nodes[brother].isClosed?2.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score8 = score/nodes[0].visits+/*fabs(score)*loss*/lossS + loss2*(nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score9 = score/nodes[0].visits+/*fabs(score)*loss*/loss + loss3*(nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score10 = score/nodes[0].visits+/*fabs(score)*loss*/loss + (loss2/ (fmax(fabs(score),1.0) )) * (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score11 = score/nodes[0].visits+/*fabs(score)*loss*/loss + (loss2/ (fmax(fabs(score),1.0) )) + (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score12 = score/nodes[0].visits+/*fabs(score)*loss*/loss + 1e5 * (loss2/ (fmax(fabs(score),1.0) )) + (nodes[brother].isClosed?2.0:1.0)*(log2(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score13 = score/nodes[0].visits+(fabs(score0)+10.0)*loss + 0.0 * (loss2/ (fmax(fabs(score),1.0) )) + (0.001+loss)*(nodes[brother].isClosed?2000000.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      double score19 = score/nodes[0].visits+(fabs(score0)+1000.0)*loss*bloss + loss + (nodes[brother].isClosed?200000.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*(0.001+loss));
	      double score20 = score/nodes[0].visits+(fabs(score0)+1000.0)*loss*bloss + loss + (nodes[brother].isClosed?1.0:1.0)*(sqrt((double)nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*(0.001+loss));
	      //double score13 = score/nodes[0].visits+fabs(score+10000.0)*loss + 0.0 * (loss2/ (fmax(fabs(score),1.0) )) + (0.001+loss)*(nodes[brother].isClosed?2000000.0:1.0)*(sqrt(nodes[0/*nodeID*/].visits)) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1);
	      //double score2 = score*loss + (sqrt(nodes[0/*nodeID*/].visits) / (nodes[w].visits/*ExistScoreCnt+nodes[w].UnivScoreCnt-1*/+1)*/*loss*/(1.0/((nodes[w].UnivScoreSum+nodes[brother].UnivScoreSum) / (nodes[w].UnivScoreCnt+nodes[brother].UnivScoreCnt))));
	      score2 = THESCORE;
	    if (score2 > bestS_inf_score) {
	      bestS_inf_score = score2;
	      bestS_inf = w;
	    }
	  }
	  
	  if (pActivity[nodes[w].entryVar] + nActivity[nodes[w].entryVar] > bestS_act_score) {
	    bestS_act_score = pActivity[nodes[w].entryVar] + nActivity[nodes[w].entryVar];
	    bestS_act = w;
	  }
	}
	kk = GM.NextAdjacentInGraph(kk);
      }
    }

    if (id==0 && doO) {
      std::cerr << "There are " << bestS_all_cnt << " successors." << std::endl;
      std::cerr << "#visits is " << nodes[id].visits << " sum over all visits = " << sumVisits << std::endl;
      if (bestS_pseudo>=0) std::cerr << "Pseusocosts point to x " << nodes[bestS_pseudo].entryVar << " with score:" << bestS_pseudo_score<< std::endl;
      if (bestS_ninf>=0) std::cerr << "Exist-stats point to x " << nodes[bestS_ninf].entryVar << " with score:" << bestS_ninf_score << " #" << nodes[bestS_ninf].ExistScoreCnt << std::endl;
      if (bestS_inf>=0) std::cerr << "Univ-stats  point to x " << nodes[bestS_inf].entryVar << " with score:" << -bestS_inf_score << " #" << nodes[bestS_inf].UnivScoreCnt << std::endl;
      if (bestS_act>=0) std::cerr << "Activites   point to x " << nodes[bestS_act].entryVar << " with score:" << bestS_act_score << std::endl;
    }

    if (0&&bestS_latest >= 0) {
      int brother = bestS_latest;
      if (nodes[bestS_latest].entryVal == 0) brother++;
      else brother--; 
      best_succ = nodes[bestS_latest].entryVar;
      if (nodes[bestS_latest].UnivScoreCnt < nodes[brother].UnivScoreCnt) 
	best_dir = nodes[brother].entryVal;
      else 
	best_dir = nodes[bestS_latest].entryVal;
    } else if (0&&bestS_pseudo >= 0 && nodes[id].visits < 10) {
	  best_succ = nodes[bestS_pseudo].entryVar;
	  best_dir = nodes[bestS_pseudo].entryVal;
    } else {
      if ((bestS_inf>=0&&bestS_ninf>=0 &&nodes[bestS_ninf].ExistScoreCnt < nodes[bestS_inf].UnivScoreCnt) ||
	  (bestS_inf>=0&&bestS_ninf<0)) {
	if (bestS_pseudo < 0) {
	  best_succ = -1;
	  return;
	}
	if (bestS_inf>=0) {
	  best_dir = nodes[bestS_inf].entryVal;
	  assert(!nodes[bestS_inf].isClosed);
	  best_succ = nodes[bestS_inf].entryVar;
	  int broth = bestS_inf;
	  double s0,s1;
	  if (nodes[bestS_inf].entryVal == 0) {
	    broth++;
	    s0 = nodes[bestS_inf].ExistScoreSum / nodes[bestS_inf].ExistScoreCnt;
	    s1 = nodes[broth].ExistScoreSum / nodes[broth].ExistScoreCnt;
	  } else {
	    broth--;
	    s1 = nodes[bestS_inf].ExistScoreSum / nodes[bestS_inf].ExistScoreCnt;
	    s0 = nodes[broth].ExistScoreSum / nodes[broth].ExistScoreCnt;
	  }
	  if (nodes[broth].ExistScoreCnt <= 1 || nodes[bestS_inf].ExistScoreCnt <= 1) s0 = s1 = 0.0;

	  /*
	  if (s0 + pseudo1score[nodes[bestS_inf].entryVar] / pseudo1scoreCnt[nodes[bestS_inf].entryVar] > s1 + pseudo0score[nodes[bestS_inf].entryVar] / pseudo0scoreCnt[nodes[bestS_inf].entryVar]) {
	    best_dir = 0;
	    if ((nodes[bestS_inf].entryVal == 0 && nodes[bestS_inf].isClosed) || (nodes[broth].entryVal == 0 && nodes[broth].isClosed)) 
	      best_dir = 1;
	  } else { 
	    best_dir = 1;
	    if ((nodes[bestS_inf].entryVal == 1 && nodes[bestS_inf].isClosed) || (nodes[broth].entryVal == 1 && nodes[broth].isClosed)) 
	      best_dir = 0;
	  }
	  */
	  if (0&&bestS_pseudo >= 0 && //(nodes[bestS_inf].visits < 10) && 
              (pseudo1score[nodes[bestS_pseudo].entryVar]+pseudo0score[nodes[bestS_pseudo].entryVar] > nodes[bestS_pseudo].visits * 1e-5 * fabs(nodes[/*bestS_pseudo*/id].ExistScoreSum / nodes[/*bestS_pseudo*/id].ExistScoreCnt))) {
	    best_succ = nodes[bestS_pseudo].entryVar;
	    int broth = bestS_pseudo;
	    if (nodes[bestS_pseudo].entryVal == 0) {
	      broth++;
	    } else {
	      broth--;
	    }
	    if (pseudo1score[nodes[bestS_pseudo].entryVar] / pseudo1scoreCnt[nodes[bestS_pseudo].entryVar] > pseudo0score[nodes[bestS_pseudo].entryVar] / pseudo0scoreCnt[nodes[bestS_pseudo].entryVar]) {
	      best_dir = 0;
	      if ((nodes[bestS_pseudo].entryVal == 0 && nodes[bestS_pseudo].isClosed) || (nodes[broth].entryVal == 0 && nodes[broth].isClosed)) 
		best_dir = 1;
	    } else { 
	      best_dir = 1;
	      if ((nodes[bestS_pseudo].entryVal == 1 && nodes[bestS_pseudo].isClosed) || (nodes[broth].entryVal == 1 && nodes[broth].isClosed)) 
		best_dir = 0;
	    }
	  } else if (0){
	    best_succ = nodes[bestS_act].entryVar;
	    if (pActivity[nodes[bestS_act].entryVar] > nActivity[nodes[bestS_act].entryVar]) 
	      best_dir = 1;
	    else 
	      best_dir = 0;
	  }
	  if (id==0 && doO) std::cerr << "chose univ side" << std::endl; 
	}
      } else {
	if (0&&bestS_pseudo < 0) {
	  best_succ = -1;
	  return;
	}
	if (bestS_ninf>=0) {
	  int brother = nodes[bestS_ninf].ID;
	  if (nodes[bestS_ninf].entryVal == 0) brother++;
	  else brother--;
	  assert(brother>=0);
	  assert(nodes[bestS_ninf].entryVar == nodes[brother].entryVar);
	  assert(!nodes[bestS_ninf].isClosed);
	  best_succ = nodes[bestS_ninf].entryVar;
	  best_dir = nodes[bestS_ninf].entryVal;
	  if (0&&eas[nodes[bestS_ninf].entryVar] == UNIV && nodes[bestS_ninf].visits > /*1 &&*/ nodes[brother].visits /*== 0*/ && !nodes[brother].isClosed) {
 
	    best_dir = nodes[brother].entryVal;
	  }
	  if (id==0 && doO) std::cerr << "chose exist side" << std::endl; 
	  //#define PSEUDO_DOM
#ifdef PSEUDO_DOM
	  if (bestS_pseudo >= 0 && //(nodes[bestS_inf].visits < 10) && 
              (pseudo1score[nodes[bestS_pseudo].entryVar]+pseudo0score[nodes[bestS_pseudo].entryVar] > nodes[bestS_pseudo].visits * 1e-5 * fabs(nodes[/*bestS_pseudo*/id].ExistScoreSum / nodes[/*bestS_pseudo*/id].ExistScoreCnt))) {
	    best_succ = nodes[bestS_pseudo].entryVar;
	    int broth = bestS_pseudo;
	    if (nodes[bestS_pseudo].entryVal == 0) {
	      broth++;
	    } else {
	      broth--;
	    }
	    if (pseudo1score[nodes[bestS_pseudo].entryVar] / pseudo1scoreCnt[nodes[bestS_pseudo].entryVar] > pseudo0score[nodes[bestS_pseudo].entryVar] / pseudo0scoreCnt[nodes[bestS_pseudo].entryVar]) {
	      best_dir = 0;
	      if ((nodes[bestS_pseudo].entryVal == 0 && nodes[bestS_pseudo].isClosed) || (nodes[broth].entryVal == 0 && nodes[broth].isClosed)) 
		best_dir = 1;
	    } else { 
	      best_dir = 1;
	      if ((nodes[bestS_pseudo].entryVal == 1 && nodes[bestS_pseudo].isClosed) || (nodes[broth].entryVal == 1 && nodes[broth].isClosed)) 
		best_dir = 0;
	    }
	  } else if (0){
	    best_succ = nodes[bestS_act].entryVar;
	    if (pActivity[nodes[bestS_act].entryVar] > nActivity[nodes[bestS_act].entryVar]) 
	      best_dir = 1;
	    else 
	      best_dir = 0;
	  }
#endif
	}
      }
    }
    return;
  }

  void updateVisits(int id, int x) {
    nodes[id].updateVisits(x);
    currentSimValid = false;
  } 

  bool gotEvalFromSucc(int id) {
    return nodes[id].gotEvalFromSucc;
  }

  void updateBounds(int nodeID, std::vector< std::pair< std::pair<double,double>, int > > &bndList, int n) {
    if (nodeID < 0) return;
    std::vector<int> bndListIdcs;
    for (int i=0; i < n;i++) {
      bndListIdcs.push_back(-1);
    }
    for (int i = 0; i < bndList.size();i++) {
      bndListIdcs[bndList[i].second] = i;
    }

    {
      int j;
      int id = nodeID;
      j = GM.FirstAdjacentInGraph(id);
      if (j >= 0) {
	int w = GM.getAdjacent(j);
	//w first succ
	int node = nodes[w].ID;
	int var = nodes[w].entryVar;
	int val = nodes[w].entryVal;
	if (bndListIdcs[var] >= 0) {
	  int i = bndListIdcs[var];
	  double nodeuBnd = (bndList[i].first.first > bndList[i].first.second ? bndList[i].first.first : bndList[i].first.second);
	  if (nodeuBnd < nodes[nodeID].upperBound) nodes[nodeID].upperBound = nodeuBnd;
	  if (nodeuBnd < nodes[node].upperBound) nodes[node].upperBound = nodeuBnd;
	}
	int kk = GM.NextAdjacentInGraph(j);
	while (kk >= 0) {
	  int w = GM.getAdjacent(kk);
	  // w next succ
	  int node = nodes[w].ID;
	  int var = nodes[w].entryVar;
	  int val = nodes[w].entryVal;
	  if (bndListIdcs[var] >= 0) {
	    int i = bndListIdcs[var];
	    double nodeuBnd = (bndList[i].first.first > bndList[i].first.second ? bndList[i].first.first : bndList[i].first.second);
	    if (nodeuBnd < nodes[nodeID].upperBound) nodes[nodeID].upperBound = nodeuBnd;
	    if (nodeuBnd < nodes[node].upperBound) nodes[node].upperBound = nodeuBnd;
	  }
	  kk = GM.NextAdjacentInGraph(kk);
	}
      }
    }
  }

  void addMIPevalExact(int nodeID, double LPval, bool isFeas) {
    addLPeval(nodeID, LPval, isFeas);
    if (!nodes[nodeID].isClosed) {
      if (isFeas)
	setClosed(nodeID, LPval, LPval);
      else 
	;//setClosed(nodeID, n_infinity, n_infinity);
    }
  }
  void addLPeval(int nodeID, double LPval, bool isFeas) {
    if (!nodes[nodeID].isClosed && (!nodes[nodeID].LPexists || !nodes[nodeID].isFeasible || nodes[nodeID].LpVal > LPval)) {
      nodes[nodeID].LpVal = LPval;
      nodes[nodeID].LPexists = true;
      nodes[nodeID].isFeasible = isFeas;
      if (!isFeas) {
	nodes[nodeID].isClosed = true;
	nodes[nodeID].LpVal = n_infinity;
	nodes[nodeID].minmax_bnd = n_infinity;
	nodes[nodeID].lowerBound = n_infinity;
	nodes[nodeID].upperBound = n_infinity;
      }
      assert(LPval >= n_infinity);
    }
  }
  bool isClosed(int nodeID, double &l, double &u) {
    if (nodes[nodeID].isClosed == true) {
      l = nodes[nodeID].lowerBound;
      u = nodes[nodeID].upperBound;
      return true;
    } else {
      l = nodes[nodeID].lowerBound;
      u = nodes[nodeID].upperBound;
      return false;
    }
    return false;
  }
  void setClosed(int nodeID, double l, double u) {
    //return;
    int brother = nodeID;
    if (nodes[nodeID].entryVal == 0) brother++;
    else brother--;
    if (nodes[nodeID].isClosed) return;
    //cerr << "close node " << nodeID << " WITH " << l << " father is " << nodes[nodeID].fatherID << " broth is " << brother;
    //if (brother >= 0) cerr << " val of broth " << nodes[brother].minmax_bnd << " isclosed:" << nodes[brother].isClosed << endl;
    //else cerr << endl;

    if (0&&nodeID == 0) {
      cerr << "set node " << nodeID << " closed. to "  << l << endl;
      assert(0); 
    }
    ////assert(nodeID != 96);
    nodes[nodeID].isClosed = true;
    nodes[nodeID].minmax_bnd = l;
    nodes[nodeID].lowerBound = l;
    nodes[nodeID].upperBound = u;
    if (!(l != dont_know && l >= n_infinity)) cerr << "set closed:" << nodeID << " l=" << l << " u=" << u << endl;
    assert(l != dont_know && l >= n_infinity);
  }
  void updateNodeScore(int nodeID) {
    assert(0);
    //cerr << "mache update bei node " << nodeID << endl;
    int f = nodes[nodeID].fatherID;
    //assert(nodes[nodeID].who2move == EXIST || nodes[nodeID].who2move == UNIV);
    assert(f < 0 || nodes[f].who2move == EXIST || nodes[f].who2move == UNIV);

    if (f >= 0) {
      if (nodes[f].minmax_bnd == n_infinity)
	setClosed(f, n_infinity, n_infinity);
      if (nodes[f].isClosed) return;

      nodes[f].AVscore = nodes[nodeID].AVscore;
      nodes[f].gotEvalFromSucc = true;
      nodes[f].winnerIsExist = nodes[nodeID].winnerIsExist;
      int father2nodeVar = nodes[nodeID].entryVar;
      int father2nodeVal = nodes[nodeID].entryVal;
      int brother = findSucc(f,father2nodeVar, 1-father2nodeVal);
      //assert(brother >= 0);
      //assert(nodes[brother].who2move >= 0);
      bool bothSuccClosed = false;
      if (nodes[f].who2move >= 0 && nodes[nodeID].isClosed && nodes[brother].isClosed) bothSuccClosed = true;
      if (bothSuccClosed) {
	assert(nodes[nodeID].minmax_bnd != dont_know && nodes[brother].minmax_bnd != dont_know);
	if (nodes[f].who2move == EXIST) {
	  nodes[f].minmax_bnd = fmax(nodes[nodeID].minmax_bnd, nodes[brother].minmax_bnd);
	  setClosed(f,nodes[f].minmax_bnd,nodes[f].minmax_bnd);
	} else {
	  if (nodes[f].minmax_bnd < fmin(nodes[nodeID].minmax_bnd, nodes[brother].minmax_bnd) || nodes[f].minmax_bnd == dont_know) 
	    nodes[f].minmax_bnd = fmin(nodes[nodeID].minmax_bnd, nodes[brother].minmax_bnd);
	  setClosed(f,nodes[f].minmax_bnd,nodes[f].minmax_bnd);
	}
	if(nodes[f].minmax_bnd == dont_know || nodes[f].lowerBound < n_infinity) {
	  cerr << "Error: closed mit Wert " << nodes[f].minmax_bnd << endl;
	  cerr << "Info: replcae value with -inf " << endl;
	  nodes[f].lowerBound = n_infinity;
	  assert(0);
	}
	assert(nodes[f].minmax_bnd != dont_know && nodes[f].lowerBound >= n_infinity);
      } else { 
	if (nodes[f].who2move == EXIST) {
	  if (nodes[nodeID].minmax_bnd == n_infinity || nodes[nodeID].minmax_bnd > dont_know) {
	    if (nodes[nodeID].minmax_bnd > nodes[f].minmax_bnd) {
	      nodes[f].minmax_bnd = nodes[nodeID].minmax_bnd;
	    }
	    if (nodes[nodeID].minmax_bnd == n_infinity) assert(nodes[nodeID].isClosed);
	  } 
	} else if (nodes[f].who2move == UNIV) {
	  if (nodes[nodeID].minmax_bnd == n_infinity || nodes[nodeID].minmax_bnd > dont_know) {
	    double minmax_bnd;
	    if (nodes[nodeID].minmax_bnd == n_infinity || nodes[brother].minmax_bnd == n_infinity) {
	      minmax_bnd = n_infinity;
	      setClosed(f,n_infinity,n_infinity);
	      if (nodes[brother].who2move < 0) {
		//cerr << "Warning: minmax yes, but fll node no" << endl;
	      }
	    } else if (nodes[nodeID].minmax_bnd == dont_know || nodes[brother].minmax_bnd == dont_know) {
	      if (nodes[nodeID].minmax_bnd == dont_know && nodes[brother].minmax_bnd == dont_know) {
		minmax_bnd = dont_know;
	      } else if (nodes[nodeID].minmax_bnd == dont_know) {
		minmax_bnd = dont_know;
	      } else if (nodes[brother].minmax_bnd == dont_know) {
		minmax_bnd = dont_know;
	      }
	    } else {
	      assert(nodes[nodeID].who2move >= 0 && nodes[brother].who2move >= 0);
	      minmax_bnd = fmin(nodes[nodeID].minmax_bnd, nodes[brother].minmax_bnd);
	    } 
	    if (nodes[nodeID].isClosed && nodes[nodeID].minmax_bnd < nodes[f].minmax_bnd) {
	      assert(nodes[nodeID].minmax_bnd == n_infinity || nodes[nodeID].minmax_bnd > dont_know);
	      nodes[f].minmax_bnd = nodes[nodeID].minmax_bnd;
	    } else {
	      if (nodes[f].minmax_bnd==dont_know || (minmax_bnd<nodes[f].minmax_bnd && minmax_bnd!=dont_know)) {
		nodes[f].minmax_bnd = minmax_bnd;
		if (nodes[f].minmax_bnd!=dont_know) nodes[f].lowerBound=nodes[f].minmax_bnd;
	      }
	      assert(!nodes[f].isClosed || nodes[f].minmax_bnd!=dont_know);
	      assert(nodes[f].minmax_bnd==dont_know || nodes[f].minmax_bnd >= n_infinity);
	      if (nodes[f].minmax_bnd > nodes[f].lowerBound) nodes[f].lowerBound=nodes[f].minmax_bnd;
	    }
	  } else {
	  }
	  if (nodes[nodeID].minmax_bnd>dont_know && nodes[brother].minmax_bnd>dont_know) {
	    //if (nodes[f].minmax_bnd == dont_know)
	    //  printNodeInfo(f, father2nodeVar, father2nodeVal); 
	    assert(nodes[f].lowerBound > dont_know); 
	  }
	} else assert(0);
      } 
    }
    if (!nodes[nodeID].winnerIsExist) {
      if (deb_o) std::cerr << "propagate univ value:" << nodes[nodeID].AVscore << std::endl;
      if(fabs(nodes[nodeID].AVscore) < fabs(dont_know)/2) {
	nodes[nodeID].UnivScoreSum = (2.0*nodes[nodeID].UnivScoreSum + 3.0*nodes[nodeID].AVscore)*0.2;
	nodes[nodeID].UnivScoreCnt = nodes[nodeID].UnivScoreCnt + 1;
      }
    } else {
      if (deb_o) std::cerr << "propagate exist value:" << nodes[nodeID].AVscore << std::endl;
      if (fabs(nodes[nodeID].AVscore) < fabs(dont_know)/2) {
	nodes[nodeID].ExistScoreSum = nodes[nodeID].ExistScoreSum + nodes[nodeID].AVscore;
	nodes[nodeID].ExistScoreCnt = nodes[nodeID].ExistScoreCnt + 1;
      } 
    }
  }

  void updateFatherScore(int nodeID, bool interactive = false) {
    int f = nodes[nodeID].fatherID;
    int b = nodeID; 
    if (nodes[nodeID].entryVal == 0) b++;
    else b--;
    if (f < 0) return;

    assert(nodes[f].who2move == EXIST || nodes[f].who2move == UNIV);
    if (nodes[f].isClosed) {
      std::cerr << "Warning: father is closed. f:" << nodes[f].minmax_bnd << " s0=" << nodes[nodeID].minmax_bnd << " s1=" << nodes[b].minmax_bnd << std::endl;
      return;
    }
    //assert(nodes[f].isClosed==false);
    bool bothSuccClosed = false;
    if (nodes[nodeID].minmax_bnd < n_infinity) nodes[nodeID].minmax_bnd = dont_know;
    if (nodes[b].minmax_bnd < n_infinity)      nodes[b].minmax_bnd = dont_know;

    if (nodes[nodeID].isClosed) {
      if (nodes[nodeID].lowerBound == dont_know || nodes[nodeID].upperBound == dont_know) {
	cerr << "Error: node " << nodeID << " is closed, but value unknown. who2move:" << nodes[nodeID].who2move << endl;
	nodes[nodeID].isClosed = false;
      }
    }
    if (nodes[b].isClosed) {
      if (nodes[b].lowerBound == dont_know || nodes[b].upperBound == dont_know) {
	cerr << "Error: node " << b << " is closed, but value unknown. who2move:" << nodes[b].who2move << endl;
	nodes[b].isClosed = false;
      }
    }
    if (nodes[f].who2move >= 0 && nodes[nodeID].isClosed && nodes[b].isClosed) bothSuccClosed = true;
    /*
      nodes[f].AVscore = nodes[nodeID].AVscore;
      nodes[f].gotEvalFromSucc = true;
      nodes[f].winnerIsExist = nodes[nodeID].winnerIsExist;
    */
    if (interactive) std::cerr << "bothSuccClosed=" << bothSuccClosed << std::endl;

    if (interactive) std::cerr << "detected son value:" << nodes[b].minmax_bnd << ", " << nodes[nodeID].minmax_bnd << std::endl;
    if (interactive) std::cerr << "who2move(father)=" << (nodes[f].who2move==EXIST?"Exist":"Univ") << std::endl;

    if (nodes[f].who2move == EXIST) {
      double mm = n_infinity;
      if (nodes[nodeID].minmax_bnd > mm) mm = nodes[nodeID].minmax_bnd;
      if (nodes[b].minmax_bnd > mm) mm = nodes[b].minmax_bnd;
      if (nodes[f].minmax_bnd == dont_know)  nodes[f].minmax_bnd = mm;
      else if (nodes[f].minmax_bnd < mm) nodes[f].minmax_bnd = mm;
      else mm = nodes[f].minmax_bnd;
      if (interactive) std::cerr << "mm=" << mm << "nodes[f].minmax_bnd=" << nodes[f].minmax_bnd << std::endl;
      if (bothSuccClosed) setClosed(f, mm, mm);
      else {
	//assert(nodes[f].minmax_bnd > n_infinity);
	if (nodes[f].minmax_bnd <= n_infinity && mm > n_infinity) std::cerr << "Warning: nodes[f].minmax_bnd <= n_infinity i.e. = " << nodes[f].minmax_bnd<< " current mm=" << mm << " son1=" << nodes[nodeID].minmax_bnd << "son2=" << nodes[b].minmax_bnd << std::endl;
      }
      if (nodes[f].minmax_bnd != dont_know) nodes[f].lowerBound = nodes[f].minmax_bnd;
      else nodes[f].lowerBound = n_infinity; 
    } else if (nodes[f].who2move == UNIV) {
      if (nodes[nodeID].minmax_bnd == n_infinity || nodes[nodeID].minmax_bnd > dont_know) {
	double minmax_bnd=dont_know;
	if (nodes[nodeID].minmax_bnd == n_infinity && nodes[nodeID].isClosed == false) {
	  nodes[nodeID].minmax_bnd = nodes[nodeID].lowerBound = dont_know;
	  nodes[nodeID].upperBound = -n_infinity;
	}
	if (nodes[nodeID].minmax_bnd == n_infinity || nodes[b].minmax_bnd == n_infinity) {
	  if (nodes[nodeID].minmax_bnd == n_infinity) assert(nodes[nodeID].isClosed);
	  if (nodes[  b   ].minmax_bnd == n_infinity) assert(nodes[b].isClosed);
	  minmax_bnd = n_infinity;
	  nodes[f].lowerBound = n_infinity; 
	  nodes[f].minmax_bnd = n_infinity; 
	  //setClosed(f,n_infinity,n_infinity);
	} else if (nodes[nodeID].minmax_bnd == dont_know || nodes[b].minmax_bnd == dont_know) {
	  if (nodes[nodeID].minmax_bnd == dont_know && nodes[b].minmax_bnd == dont_know) {
	    minmax_bnd = dont_know;
	  } else if (nodes[nodeID].minmax_bnd == dont_know) {
	    minmax_bnd = dont_know;
	    if (nodes[b].isClosed) {
	      nodes[f].upperBound = nodes[b].minmax_bnd; 
	    }
	  } else if (nodes[b].minmax_bnd == dont_know) {
	    minmax_bnd = dont_know;
	    if (nodes[nodeID].isClosed) {
	      nodes[f].upperBound = nodes[nodeID].minmax_bnd; 
	    }
	  }
	} else {
	  assert(nodes[nodeID].who2move >= 0 && nodes[b].who2move >= 0);
	  minmax_bnd = fmin(nodes[nodeID].minmax_bnd, nodes[b].minmax_bnd);
	}
	if (minmax_bnd != dont_know) {
	  nodes[f].minmax_bnd = minmax_bnd; 
	  //if (bothSuccClosed) setClosed(f, minmax_bnd, minmax_bnd);
	}

	if (nodes[nodeID].isClosed && nodes[nodeID].minmax_bnd < nodes[f].minmax_bnd) {
	  assert(0); //should be outdated
	  assert(nodes[nodeID].minmax_bnd == n_infinity || nodes[nodeID].minmax_bnd > dont_know);
	  nodes[f].minmax_bnd = nodes[nodeID].minmax_bnd;
	} else {
	  if (nodes[f].minmax_bnd==dont_know || (minmax_bnd<nodes[f].minmax_bnd && minmax_bnd!=dont_know)) {
	    assert(minmax_bnd==dont_know); //should be outdated   
	    nodes[f].minmax_bnd = minmax_bnd;
	    if (nodes[f].minmax_bnd!=dont_know) nodes[f].lowerBound=nodes[f].minmax_bnd;
	  }
	  assert(!nodes[f].isClosed || nodes[f].minmax_bnd!=dont_know);
	  assert(nodes[f].minmax_bnd==dont_know || nodes[f].minmax_bnd >= n_infinity);
	  if (nodes[f].minmax_bnd > nodes[f].lowerBound) {
	    //assert(0); //should be outdated   
	    nodes[f].lowerBound=nodes[f].minmax_bnd;
	  }
	}
      }
    } else assert(0);
  }

  double getMinmaxEst(int nodeID) {
    return nodes[nodeID].minmax_bnd;
  }

  void evalSimulation(int nodeID, int decisionLevel, bool break_from_outside) {
    int f = nodes[nodeID].fatherID;
    double value = currentSimVal;
    assert(value != dont_know);
    bool WinnerIsExist = currentSimWinnerIsExist;
    if (currentSimValid == false) {
      //std::cerr << "Warning: Mcts has no valid simulation evaluation." << std::endl;
      nodes[nodeID].innerNode = true;
      //assert(0);
      return;
    }
    if (!WinnerIsExist) {
      nodes[nodeID].winnerIsExist = false;
      nodes[nodeID].AVscore = value;
    } else {
      nodes[nodeID].winnerIsExist = true;
      nodes[nodeID].AVscore = value;
    }
    if (f >= 0) {
      nodes[f].AVscore = nodes[nodeID].AVscore;
      nodes[f].gotEvalFromSucc = true;
      nodes[f].winnerIsExist = nodes[nodeID].winnerIsExist;
      if (nodes[f].who2move == EXIST) {
	if (nodes[f].winnerIsExist) {
	  if (value != dont_know && value > nodes[f].minmax_bnd && blocks[nodes[nodeID].entryVar] ==maxBlock) {
	    nodes[f].minmax_bnd = value;
	  }
	} else {
	  if (nodes[nodeID].isClosed && nodes[nodeID].minmax_bnd > nodes[f].minmax_bnd) {
	    nodes[f].minmax_bnd = nodes[nodeID].minmax_bnd;
	  }
	}
      }
    }
    nodes[nodeID].gotEvalFromSucc = false;
    if (!nodes[nodeID].winnerIsExist) {
      if (deb_o) std::cerr << "add univ value:" << nodes[nodeID].AVscore << std::endl;
      if(fabs(nodes[nodeID].AVscore) < fabs(dont_know)/2) {
	nodes[nodeID].UnivScoreSum = (2.0*nodes[nodeID].UnivScoreSum + 3.0*nodes[nodeID].AVscore)*0.2;
	nodes[nodeID].UnivScoreCnt = nodes[nodeID].UnivScoreCnt + 1;
	nodes[nodeID].gotEvalFromSucc = true;
      }
    } else {
      if (deb_o) std::cerr << "add exist value:" << nodes[nodeID].AVscore << std::endl;
      if(fabs(nodes[nodeID].AVscore) < fabs(dont_know)/2) {
	nodes[nodeID].ExistScoreSum = nodes[nodeID].ExistScoreSum + nodes[nodeID].AVscore;
	nodes[nodeID].ExistScoreCnt = nodes[nodeID].ExistScoreCnt + 1;
	nodes[nodeID].gotEvalFromSucc = true;
      }
    } 
    nodes[nodeID].innerNode = true;
  }

  void setSimulationValue(double value,int D) {
    //if (value == dont_know) assert(0);
    if (value == dont_know) return;
    if (value < dont_know) {
      currentSimWinnerIsExist = false;
      currentSimVal = (double)D;
      //std::cerr << "set univ value:" << (double)D << std::endl;
    } else {
      currentSimWinnerIsExist = true;
      currentSimVal = value;
      if (deb_o) std::cerr << "set exist value:" << value << std::endl;
    }
    currentSimValid = true;
  }

  void write_nodeinfo(int nodeID){
    std::cout << "node " << nodeID << ": ID=" << nodes[nodeID].ID << " x" << nodes[nodeID].entryVar<< " = " << nodes[nodeID].entryVal << " lb=" << nodes[nodeID].lowerBound << " ub=" << nodes[nodeID].upperBound << " fID=" << nodes[nodeID].fatherID << " mm=" << nodes[nodeID].minmax_bnd << std::endl;
  }
  void write_successors(int nodeID) {
    int cnt =0;
    int j;
    j = GM.FirstAdjacentInGraph(nodeID);
    if (j >= 0) {
      int w = GM.getAdjacent(j);
      //w first succ
      std::cout << "successor " << cnt++ << ": " << nodes[w].ID << " mmx:" << nodes[w].minmax_bnd << " sim:" << nodes[w].ExistScoreSum / nodes[w].ExistScoreCnt << " closed:" << nodes[w].isClosed << std::endl;
      int kk = GM.NextAdjacentInGraph(j);
      while (kk >= 0) {
	int w = GM.getAdjacent(kk);
	// w next succ
	std::cout << "successor " << cnt++ << ": " << nodes[w].ID << " mmx:" << nodes[w].minmax_bnd << " sim:" << nodes[w].ExistScoreSum / nodes[w].ExistScoreCnt << " closed:" << nodes[w].isClosed << std::endl;
	kk = GM.NextAdjacentInGraph(kk);
      }
    }
  }
  
  void partialExpandOrUpdateNode(int id, std::vector<int> &vars, std::vector<int> &vals, int n,
			 double* pseudo0score,
			 double* pseudo1score,
			 uint64_t* pseudo0scoreCnt,
			 uint64_t* pseudo1scoreCnt,
			 double* pActivity,
			 double* nActivity,
			 bool enforceMemory) {
    if (!nodes[id].innerNode && !enforceMemory) return;
    if (vars.size() == 0) return;
    std::vector<int> newLabels;
    for (int i=0; i < n;i++)
      newLabels.push_back(-1);
    assert(vars.size() == vals.size());
    for (int i = 0; i < vars.size();i++) {
      assert(vals[i] == 0 || vals[i] == 1);
      newLabels[vars[i]] = vals[i];
    }

    {
      /*
	for all succesors of id {
	if (label exists in successors)
	newLabels[...] = -1;
	}
      */
      int j;
      j = GM.FirstAdjacentInGraph(id);
      if (j >= 0) {
	int w = GM.getAdjacent(j);
	//w first succ
	if ( newLabels[ nodes[w].entryVar ] == nodes[w].entryVal) {
	  if (nodes[w].entryVal == 0) {
	    nodes[w].updatePseudoCost(&pseudo0score[nodes[w].entryVar], &pseudo0scoreCnt[nodes[w].entryVar]);
	    nodes[w].updateActivity(&nActivity[nodes[w].entryVar]);
	  } else {
	    nodes[w].updatePseudoCost(&pseudo1score[nodes[w].entryVar], &pseudo1scoreCnt[nodes[w].entryVar]);
	    nodes[w].updateActivity(&pActivity[nodes[w].entryVar]);
	  } 
	  newLabels[ nodes[w].entryVar ] = -1;
	}
	int kk = GM.NextAdjacentInGraph(j);
	while (kk >= 0) {
	  int w = GM.getAdjacent(kk);
	  // w next succ
	  if ( newLabels[ nodes[w].entryVar ] == nodes[w].entryVal) {
	    if (nodes[w].entryVal == 0) {
	      nodes[w].updatePseudoCost(&pseudo0score[nodes[w].entryVar], &pseudo0scoreCnt[nodes[w].entryVar]);
	      nodes[w].updateActivity(&nActivity[nodes[w].entryVar]);
	    } else {
	      nodes[w].updatePseudoCost(&pseudo1score[nodes[w].entryVar], &pseudo1scoreCnt[nodes[w].entryVar]);
	      nodes[w].updateActivity(&pActivity[nodes[w].entryVar]);
	    } 
	    newLabels[ nodes[w].entryVar ] = -1;
	  }
	  kk = GM.NextAdjacentInGraph(kk);
	}
      }
    }

    assert(vars.size() > 0);
    if (nodes[id].who2move != eas[vars[0]]) {
      nodes[id].who2move = eas[vars[0]];
      assert(vars.size() < 3);
    }
    int w2m = nodes[id].who2move;
    for (int i = 0; i < vars.size();i++) {
      if (eas[vars[i]] != w2m)
	cerr << "w2m=" << w2m << " eas..:" << eas[vars[i]] << " var:" <<vars[i] << endl;
      assert(eas[vars[i]] == w2m);
      if (newLabels[vars[i]] != -1) {
	int newID = nodes.size();
	Node nd(newID,id,-1,-1,-1,dont_know);
	nodes.push_back(nd);
	nodes[newID].updateLabel(vars[i], vals[i]);
	if (vals[i] == 0) {
	  nodes[newID].updatePseudoCost(&pseudo0score[vars[i]], &pseudo0scoreCnt[vars[i]]);
	  nodes[newID].updateActivity(&nActivity[vars[i]]);
	} else {
	  nodes[newID].updatePseudoCost(&pseudo1score[vars[i]], &pseudo1scoreCnt[vars[i]]);
	  nodes[newID].updateActivity(&pActivity[vars[i]]);
	} 
	GM.AddEdge(id,newID);
	//(id,vars[i],vals[i]) -> (newID) with searching, first of all
      }
    }

  }
};

#endif

#ifdef DONT_FORGET_THE_OLD
  void findBestSucc(int id, int &best_succ, int &best_dir, 
		    double* pseudo0score,
		    double* pseudo1score,
		    uint64_t* pseudo0scoreCnt,
		    uint64_t* pseudo1scoreCnt,
		    double* pActivity,
		    double* nActivity,
		    int8_t *assigns,
		    double upperBnd,
		    int n) {
    int j;
    int best_w=-1;
    double best_score = -std::numeric_limits<double>::max();
    best_succ = -1;

    if (1||id < 0) return;
    if (nodes[id].isClosed) {
      best_succ = -2;
      return;
    }

    best_dir = 0;
    std::vector<bool> pseudo0Exists;
    std::vector<bool> pseudo1Exists;
    std::vector<double> loss0;
    std::vector<double> loss1;
    std::vector<int> visits0;
    std::vector<int> visits1;
    for (int i = 0; i <n; i++) {
      pseudo0Exists.push_back(false);
      pseudo1Exists.push_back(false);
      loss0.push_back(0.0);
      loss1.push_back(0.0);
      visits0.push_back(0);
      visits1.push_back(0);
    }
    j = GM.FirstAdjacentInGraph(id);
    if (j >= 0) {
      int w = GM.getAdjacent(j);
      //w first succ

      //#define MCTS_LIKE
//#define WEIGHTED_PSEUDO
#ifdef WEIGHTED_PSEUDO
	if ( nodes[w].entryVar >= 0 && assigns[nodes[w].entryVar] == 2 ) {
	  if (nodes[w].entryVal == 0) {
	    visits0[nodes[w].entryVar] = nodes[w].visits;
	    if (pseudo0scoreCnt[nodes[w].entryVar] >= 1) {
	      pseudo0Exists[nodes[w].entryVar] = true;
	      loss0[nodes[w].entryVar] = //nodes[w].visits+1;
		(1.0+0.1*(double)pseudo0scoreCnt[nodes[w].entryVar]/*nodes[w].visits*/)* (pseudo0score[nodes[w].entryVar] / (double)pseudo0scoreCnt[nodes[w].entryVar]); 
	    }
	  } else if (nodes[w].entryVal == 1) {
	    visits1[nodes[w].entryVar] = nodes[w].visits;
	    if (pseudo1scoreCnt[nodes[w].entryVar] >= 1) {
	      pseudo1Exists[nodes[w].entryVar] = true;
	      loss1[nodes[w].entryVar] = //nodes[w].visits+1;
		(1.0+0.1*(double)pseudo1scoreCnt[nodes[w].entryVar]/*nodes[w].visits*/)* (pseudo1score[nodes[w].entryVar] / (double)pseudo1scoreCnt[nodes[w].entryVar]); 
	    }
	  } else assert(0);
	}
#endif
#ifdef MCTS_LIKE
	if ( nodes[w].entryVar >= 0 && assigns[nodes[w].entryVar] == 2 ) {
	  if (nodes[w].entryVal == 0) {
	    visits0[nodes[w].entryVar] = nodes[w].visits;
	    if (1||pseudo0scoreCnt[nodes[w].entryVar] >= 1) {
	      pseudo0Exists[nodes[w].entryVar] = true;
	      loss0[nodes[w].entryVar] = 
		(10.0 / (10.0+(double)nodes[w].visits) ) * (pseudo0score[nodes[w].entryVar] / (double)pseudo0scoreCnt[nodes[w].entryVar]); 
	    }
	  } else if (nodes[w].entryVal == 1) {
	    visits1[nodes[w].entryVar] = nodes[w].visits;
	    if (1||pseudo1scoreCnt[nodes[w].entryVar] >= 1) {
	      pseudo1Exists[nodes[w].entryVar] = true;
	      loss1[nodes[w].entryVar] = 
		(10.0 / (10.0+(double)nodes[w].visits) )* (pseudo1score[nodes[w].entryVar] / (double)pseudo1scoreCnt[nodes[w].entryVar]); 
	    }
	  } else assert(0);
	}
#endif
      int kk = GM.NextAdjacentInGraph(j);
      while (kk >= 0) {
	int w = GM.getAdjacent(kk);
	// w next succ
#ifdef WEIGHTED_PSEUDO
	if ( nodes[w].entryVar >= 0 && assigns[nodes[w].entryVar] == 2 ) {
	  if (nodes[w].entryVal == 0) {
	    visits0[nodes[w].entryVar] = nodes[w].visits;
	    if (pseudo0scoreCnt[nodes[w].entryVar] >= 1) {
	      pseudo0Exists[nodes[w].entryVar] = true;
	      loss0[nodes[w].entryVar] = //nodes[w].visits+1;
		(1.0+0.1*(double)pseudo0scoreCnt[nodes[w].entryVar]/*nodes[w].visits*/)* (pseudo0score[nodes[w].entryVar] / (double)pseudo0scoreCnt[nodes[w].entryVar]); 
	    }
	  } else if (nodes[w].entryVal == 1) {
	    visits1[nodes[w].entryVar] = nodes[w].visits;
	    if (pseudo1scoreCnt[nodes[w].entryVar] >= 1) {
	      pseudo1Exists[nodes[w].entryVar] = true;
	      loss1[nodes[w].entryVar] = //nodes[w].visits+1;
		(1.0+0.1*(double)pseudo1scoreCnt[nodes[w].entryVar]/*nodes[w].visits*/)* (pseudo1score[nodes[w].entryVar] / (double)pseudo1scoreCnt[nodes[w].entryVar]); 
	    }
	  } else assert(0);
	}
#endif
#ifdef MCTS_LIKE
	if ( nodes[w].entryVar >= 0 && assigns[nodes[w].entryVar] == 2 ) {
	  if (nodes[w].entryVal == 0) {
	    visits0[nodes[w].entryVar] = nodes[w].visits;
	    if (1||pseudo0scoreCnt[nodes[w].entryVar] >= 1) {
	      pseudo0Exists[nodes[w].entryVar] = true;
	      loss0[nodes[w].entryVar] = 
		(10.0 / (10.0+(double)nodes[w].visits) ) * (pseudo0score[nodes[w].entryVar] / (double)pseudo0scoreCnt[nodes[w].entryVar]); 
	    }
	  } else if (nodes[w].entryVal == 1) {
	    visits1[nodes[w].entryVar] = nodes[w].visits;
	    if (1||pseudo1scoreCnt[nodes[w].entryVar] >= 1) {
	      pseudo1Exists[nodes[w].entryVar] = true;
	      loss1[nodes[w].entryVar] = 
		(10.0 / (10.0+(double)nodes[w].visits) )* (pseudo1score[nodes[w].entryVar] / (double)pseudo1scoreCnt[nodes[w].entryVar]); 
	    }
	  } else assert(0);
	}
#endif
	kk = GM.NextAdjacentInGraph(kk);
      }
    }

    //best_succ = -1;
    for (int i=0;i<n;i++) {
      if ((pseudo0Exists[i] && pseudo1Exists[i])) {
	double loss = loss0[i] * loss1[i];
	if (loss > best_score) {
	  best_score = loss;
	  best_succ = i;
	  if (loss0 > loss1) {//pseudo0score[i] / (double)pseudo0scoreCnt[i] < pseudo1score[i] / (double)pseudo1scoreCnt[i]) {
	    best_dir = 0;
	  } else {
	    best_dir = 1;
	  }
	}
      }
    }
    return;
    //if (best_w < 0) return;
    //best_succ = nodes[best_w].entryVar;
    //best_dir = nodes[best_w].entryVal;
  }
#endif
