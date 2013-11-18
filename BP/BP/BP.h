#pragma once

#include <vector>

template <class T> class TPt;
template <class T, class U> class TNodeEDatNet;
class TFlt;

#ifdef _USE_libDAI
void ExactBP_Marginalization(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void ExactBP_Marginalization(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);
#endif

void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);

void ParallelBPFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void ParallelBPFromNode_LevelSynchronous(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);

void ParallelBPFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);
void ParallelBPFromNode_LevelSynchronous(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs);

void BP2(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const TPt<TNodeEDatNet<TFlt, TFlt>>& pRankGraph, int sourceNodeID);
void ParallelBP2(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const TPt<TNodeEDatNet<TFlt, TFlt>>& pRankGraph, int sourceNodeID);