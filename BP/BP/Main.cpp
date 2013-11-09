#include "BP.h"
#include "Utilities.h"

#define __TBB_NO_IMPLICIT_LINKAGE	1
#include <tbb/tick_count.h>

#include <iostream>

#include <Snap.h>

#include <queue>
#include <map>

using namespace std;

// Discussion about variable declaration and loops
// http://stackoverflow.com/questions/982963/is-there-any-overhead-to-declaring-a-variable-within-a-loop-c

//#define _TEST_soc_pokec_relationships
//#define _TEST_Email_EuAll
//#define _TEST_p2p_Gnutella09

//#define _LOAD_FROM_FILE
//#define _SAVE_TO_FILE
//#define _TEST_GRAPH
#define _TEST_DAG2

/*
function Dijkstra(Graph, source):
      for each vertex v in Graph:                           // Initializations
          dist[v]      := infinity;                         // Mark distances from source to v as not yet computed
          visited[v]   := false;                            // Mark all nodes as unvisited
          previous[v]  := undefined;                        // Previous node in optimal path from source
      end for
      
      dist[source]  := 0;                                   // Distance from source to itself is zero
      insert source into Q;                                 // Start off with the source node
                                                                
      while Q is not empty:                                 // The main loop
          u := vertex in Q with smallest distance in dist[] and has not been visited;  // Source node in first case
          remove u from Q;
          visited[u] := true                                // mark this node as visited
          
          for each neighbor v of u:   
              alt := dist[u] + dist_between(u, v);          // accumulate shortest dist from source
              if alt < dist[v] && !visited[v]:                                 
                  dist[v]  := alt;                          // keep the shortest dist from src to v
                  previous[v]  := u;
                  insert v into Q;                          // Add unvisited v into the Q to be processed
              end if
          end for
      end while
      return dist;
  endfunction
*/

//! Input : Directed graph with initialized weight edges.
//! Edges with a propagation probability strictly greater than dThreshold are ignored
TPt<TNodeEDatNet<TFlt, TFlt>> MIOA(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNode, double dThreshold)
{
	double logThreshold = log(dThreshold);
	if(dThreshold==0)
		logThreshold=-DBL_MAX;

	//////////////////////////////////////////////////////////////
	// Compte the Maximum Influence Out-Arborescence with Dijkstra
	
	// Create a copy of pGraph
	auto pDAG2Graph = TNodeEDatNet<TFlt, TFlt>::New();
	int NodeID;
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
	{
		NodeID = NI.GetId();
		pDAG2Graph->AddNode(NodeID);
		pDAG2Graph->SetNDat(NodeID, FLT_MAX);
	}

	// List of visited nodes
	std::map<int, bool> visitedNodes;
	// Stores the edge vertices to build the final DAG
	std::map<int, int> mapPrevious;
	std::vector<int> vNonVisitedNodes;

	// Distance from source node to itself is 0
	pDAG2Graph->SetNDat(sourceNode, 0);
	vNonVisitedNodes.push_back(sourceNode);

	// Beginning of the loop of Dijkstra algorithm

	while(!vNonVisitedNodes.empty())
	{
		// Find the vertex in queue with the smallest distance
		auto itVertexWithSmallestDistance = std::min_element(vNonVisitedNodes.begin(), vNonVisitedNodes.end());
		int iParentID = *itVertexWithSmallestDistance;
		vNonVisitedNodes.erase(itVertexWithSmallestDistance);
		visitedNodes.insert(std::make_pair(iParentID, true));
		auto parent = pGraph->GetNI(iParentID);

		int numChildren = parent.GetOutDeg();
		// Update the belief of the children of parent
		for(int i = 0; i < numChildren; ++i)
		{
			int iChildID = parent.GetOutNId(i);
			// Accumulate shortest dist from source
			double alt = pDAG2Graph->GetNDat(iParentID) + log(parent.GetOutEDat(i).Val);
			if(alt >= logThreshold)
			{
				auto it = visitedNodes.find(iChildID);
				if (alt < pDAG2Graph->GetNDat(iChildID) && it == visitedNodes.end())
				{
					pDAG2Graph->SetNDat(iChildID, alt);
					auto itPrevious = mapPrevious.find(iChildID);
					if(itPrevious==mapPrevious.end())
						mapPrevious.insert(std::make_pair(iChildID, iParentID));
					else
						itPrevious->second = iParentID;
					// Add unvisited iChildID into queue to be processed
					vNonVisitedNodes.push_back(iChildID);                          
				}
			}
		}

	}

	// End of Dijkstra

	// Copy the initial belief of each node into pDAG2Graph
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
		pDAG2Graph->SetNDat( NI.GetId(), FLT_MAX);

	for(auto it=mapPrevious.begin(); it!= mapPrevious.end(); ++it)
	{
		pDAG2Graph->AddEdge(it->second, it->first);
		pDAG2Graph->SetEDat(it->second, it->first, pGraph->GetEDat(it->second, it->first));
	}

	// pDAG2Graph is the MIOA starting from sourceNode

	return pDAG2Graph;
}

void TestGraph(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNode, int numIterations)
{
	// Start traversing the graph
	cout << "Starting BP from nodeID " << sourceNode << ". The input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";
	tbb::tick_count tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNode(pGraph, sourceNode);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNode_1DPartitioning(pGraph, sourceNode);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel 1D partitioning BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		PropagateFromNode(pGraph, sourceNode);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for serial BP: " << dElapsedTime/numIterations << " seconds\n";
}

int main(int argc, char* argv[])
{
#ifdef _TEST_GRAPH
	auto pGraph = TNodeEDatNet<TFlt, TFlt>::New();
	pGraph->AddNode(0);
	pGraph->SetNDat(0, 1);
	pGraph->AddNode(1);
	pGraph->SetNDat(1, 0);
	pGraph->AddNode(2);
	pGraph->SetNDat(2, 0);
	pGraph->AddEdge(0,1);
	pGraph->SetEDat(0,1, 0.5);
	pGraph->AddEdge(0,2);
	pGraph->SetEDat(0,2, 0.5);
	pGraph->AddEdge(1,2);
	pGraph->SetEDat(1,2, 0.5);
	tbb::tick_count tic = tbb::tick_count::now();
	//PropagateFromNode(pGraph, 0);
	//ParallelBPFromNode(pGraph, 0);
	ParallelBPFromNode_1DPartitioning(pGraph, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	std:: cout << pGraph->GetNDat(0).Val << " " << pGraph->GetNDat(1).Val << " " << pGraph->GetNDat(2).Val << std::endl;
#endif

#ifdef _TEST_DAG2
	auto pGraph = TNodeEDatNet<TFlt, TFlt>::New();
	pGraph->AddNode(0);
	pGraph->SetNDat(0, 1);
	pGraph->AddNode(1);
	pGraph->SetNDat(1, 0);
	pGraph->AddNode(2);
	pGraph->SetNDat(2, 0);
	pGraph->AddEdge(0,1);
	pGraph->SetEDat(0,1, 0.5);
	pGraph->AddEdge(0,2);
	pGraph->SetEDat(0,2, 0.5);
	pGraph->AddEdge(1,2);
	pGraph->SetEDat(1,2, 0.5);
	pGraph->AddNode(3);
	pGraph->AddEdge(2,3);
	pGraph->SetEDat(2,3, 0.25);
	pGraph->AddEdge(3,1);
	pGraph->SetEDat(3,1, 0.35);
	pGraph = MIOA(pGraph, 0, 0.5);
	TSnap::SaveGViz(pGraph, "test.gv", "Test DAG", true);
#endif

#ifdef _LOAD_FROM_FILE
	TFIn FIn("test.graph");
	auto pGraph = TNodeEDatNet<TFlt, TFlt>::Load(FIn);
	TestGraph(pGraph, 0, 100);
#endif

#ifdef _SAVE_TO_FILE
	auto pGraph = GenerateRandomBayesianNetwork(1, 5, 5, 7, 30);
	TFOut FOut("test.graph");
	pGraph->Save(FOut);
	TSnap::SaveGViz(pGraph, "test.gv", "Test DAG", true);
#endif

#ifdef _TEST_p2p_Gnutella09
	auto pGraph = TSnap::LoadEdgeList<TPt<TNodeEDatNet<TFlt, TFlt>>>("p2p-Gnutella09.txt", 0, 1);
	TestGraph(pGraph, 0, 100);
#endif

#ifdef _TEST_Email_EuAll
	auto pGraph = TSnap::LoadEdgeList<TPt<TNodeEDatNet<TFlt, TFlt>>>("Email-EuAll.txt", 0, 1);
	TestGraph(pGraph, 0, 100);
#endif

#ifdef _TEST_soc_pokec_relationships
	auto pGraph = TSnap::LoadConnList<TPt<TNodeEDatNet<TFlt, TFlt>>>("soc-pokec-relationships.txt");
	TestGraph(pGraph, 1, 1);
#endif

#ifdef _WIN32
	system("pause");
#endif
#ifdef __linux__
	SystemPause();
#endif

	return 0;
}
