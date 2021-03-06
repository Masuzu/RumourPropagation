#pragma once

//#define _DEBUG_INFO	1

#include <vector>
#include <map>

template <class T> class TPt;
template <class T, class U> class TNodeEDatNet;
class TFlt;

#ifdef _USE_libDAI
void ExactBP_Marginalization(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void ExactBP_Marginalization(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);
#endif

void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);
void PropagateFromNodeDAG(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void PropagateFromNodeDAG(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);

void ParallelBPFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void ParallelBPFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);
void ParallelBPFromNodeDAG(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void ParallelBPFromNodeDAG(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);

void ParallelBPFromNode_LevelSynchronous(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void ParallelBPFromNode_LevelSynchronous(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);
void ParallelBPFromNodeDAG_LevelSynchronous(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void ParallelBPFromNodeDAG_LevelSynchronous(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);

//! Instead of updating all the children of the node processed, we update each node only if all its parents have been updated
//! Advatage : no synchronization needed
//! Drawback : maximum distance from the seed nodes required
void ParallelBPFromNode_SingleNodeUpdate(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int>& pRankMap, int sourceNodeID);
void ParallelBPFromNode_SingleNodeUpdate(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int>& pRankMap, const std::vector<int> &vSeedIDs);
void ParallelBPFromNode_SingleNodeUpdate(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const TPt<TNodeEDatNet<TFlt, TFlt>>& pRankGraph, int sourceNodeID);
void ParallelBPFromNode_SingleNodeUpdate(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const TPt<TNodeEDatNet<TFlt, TFlt>>& pRankGraph, const std::vector<int> &vSeedIDs);

double InfluenceSpreadFromSeedNodes (const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph);

void GreedyCELF(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int setSize,std::vector<int> &vSeedSet);
void ParallelGreedyCELF(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int setSize, std::vector<int> &vSeedSet);
void NewGreedIC(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int setSize, int numRounds, std::vector<int> &vSeedSet);
void ParallelNewGreedIC(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int setSize, int numRounds, std::vector<int> &vSeedSet);
void ParallelNewGreedIC_Nested(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int setSize, int numRounds, std::vector<int> &vSeedSet);