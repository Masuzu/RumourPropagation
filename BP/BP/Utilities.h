#pragma once

template <class T> class TPt;
template <class T, class U> class TNodeEDatNet;
class TFlt;

#ifdef __linux__
void SystemPause();
#endif

int GetPhysicalProcessorCount();

//! @param minNumNodesPerRank How 'fat' the DAG should be.
//! @param minRanks How 'tall' the DAG should be.
//! @param probabilityEdge Chance of having an Edge in percent.
TPt<TNodeEDatNet<TFlt, TFlt>> GenerateRandomBayesianNetwork(unsigned int minNumNodesPerRank, unsigned int maxNumNodesPerRank, unsigned int minRanks, unsigned int maxRanks, unsigned int probabilityEdge);