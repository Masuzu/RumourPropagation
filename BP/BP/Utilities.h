#pragma once

#include <vector>
#include <functional>

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

TPt<TNodeEDatNet<TFlt, TFlt>> CopyGraph(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph);

void AddSuperRootNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedNodes, int superRootNodeID = INT_MAX);

void RandomGraphInitialization(TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph);

//! Set the values for each node of pGraph to 0
void ResetGraphBelief(TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph);

//! Input : Directed graph with initialized weight edges.
//! Edges with a propagation probability strictly greater than dThreshold are ignored.
//! @note Initial nodes values are initialized to FLT_MAX.
TPt<TNodeEDatNet<TFlt, TFlt>> MIOA(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNode, double dThreshold=0.0);

TPt<TNodeEDatNet<TFlt, TFlt>> GenerateDAG1(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, const std::vector<int>& seedNodes, double threshold = 0.0);

TPt<TNodeEDatNet<TFlt, TFlt>> GenerateDAG1(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, int sourceNode, double threshold = 0.0);

TPt<TNodeEDatNet<TFlt, TFlt>> GenerateDAG2(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs, double dThreshold=0.0);

TPt<TNodeEDatNet<TFlt, TFlt>> GenerateDAG2(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNode, double dThreshold=0.0);

//! Compute the rank from sourceNode using BFS.
//! Nodes already explored may be added back.
void CalculateRankFromSource(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, int sourceNode, std::vector<int> &vResult);
void CalculateRankFromSource(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, const std::vector<int> vSeedNodes, std::vector<int> &vResult);

//! Compute the rank from sourceNode using Bellman Ford with negative unitary weights. sourceNode has rank 0.
//! Not viable for large graphs.
void CalculateRankFromSource_BellmanFord(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, int sourceNode, std::vector<int> &vResult);
void CalculateRankFromSource_BellmanFord(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, const std::vector<int> vSeedNodes, std::vector<int> &vResult);

//! @note pGraph1 and pGraph2 must have the same nodes
//! Nodes with a null belief are not taken into account.
double BPError(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph1, const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph2, const std::function<double(double, double)> &);

void SaveEdgeWeightsToFile(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::string &fileName);

void LoadEdgeWeightsFromFile(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::string &fileName);