#pragma once

template <class T> class TPt;
template <class T, class U> class TNodeEDatNet;
class TFlt;

#ifdef _DEBUG
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID, bool bDisplayInfo = false);
#else
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
#endif


void ParallelBPFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
void ParallelBPFromNode_1DPartitioning(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID);
