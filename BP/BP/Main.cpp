#include <queue>
#include <iostream>
#include <vector>
#include <map>
#include "Utilities.h"
#include "Snap.h"

#define __TBB_NO_IMPLICIT_LINKAGE	1
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_queue.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/tick_count.h>
#include <tbb/spin_mutex.h>

#ifdef _WIN64
#ifndef _DEBUG
#pragma comment(lib, "intel64/vc11/tbb.lib")
#else
#pragma comment(lib, "intel64/vc11/tbb_debug.lib")
#endif
#else
#ifndef _DEBUG
#pragma comment(lib, "ia32/vc11/tbb.lib")
#else
#pragma comment(lib, "ia32/vc11/tbb_debug.lib")
#endif
#endif

using namespace std;

// http://snap.stanford.edu/snap/quick.html
#ifdef _DEBUG
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>> pGraph, int sourceNodeID, bool bDisplayInfo = false)
#else
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>> pGraph, int sourceNodeID)
#endif
{
	// Used to store the nodes which have been traversed during the breadth-first search traversal
	std::map<int, bool> visitedNodes;
	std::queue<int> queue;
	queue.push(sourceNodeID);
	visitedNodes[sourceNodeID] = true;
#ifdef _DEBUG
	int numEdge = 0;
#endif
	while(!queue.empty())
	{
		int nodeID = queue.front();
		// For debug only
#ifdef _DEBUG
		if(bDisplayInfo)
			std::cout << "Parsing node " << nodeID << " - Number of edges traversed: " << ++numEdge << std::endl;
#endif

		auto parent = pGraph->GetNI(nodeID);	
		int numChildren = parent.GetOutDeg();
		// Update the belief of the children of parent
		for(int i = 0; i < numChildren; ++i)
		{
			int iChildID = parent.GetOutNId(i);
			// Do p(u) = 1-(1-p(u))*(1-p(u,v)*p(v))
			pGraph->SetNDat(iChildID,
				1.0 - (1-parent.GetOutNDat(i).Val) * (1-parent.GetOutEDat(i).Val*parent.GetDat().Val));
		
			if(visitedNodes.count(iChildID) == 0)
			{
				visitedNodes[iChildID] = true;	// Mark the child
				queue.push(iChildID);
			}
		}

		queue.pop();
	}
}

void ParallelBPFromNode(TPt<TNodeEDatNet<TFlt, TFlt>> pGraph, int sourceNodeID)
{
	// Like std::list, insertion of new items does not invalidate any iterators, nor change the order of items already in the map. Insertion and traversal may be concurrent.
	tbb::concurrent_unordered_map<int, bool> visitedNodes;
	visitedNodes[sourceNodeID] = true;

	tbb::concurrent_queue<int> queue;
	queue.push(sourceNodeID);

	int nodeID;
	while(queue.try_pop(nodeID))
	{
		auto parent = pGraph->GetNI(nodeID);
		int numChildren = parent.GetOutDeg();
		if(numChildren > 50)
			tbb::parallel_for(tbb::blocked_range<int>(0,numChildren, 50), 
			[&](const tbb::blocked_range<int>& r)
		{
			for (int i=r.begin();i!=r.end();++i)
			{
				int iChildID = parent.GetOutNId(i);
				// Do p(u) = 1-(1-p(u))*(1-p(u,v)*p(v))
				pGraph->SetNDat(iChildID,
					1.0 - (1-parent.GetOutNDat(i).Val) * (1-parent.GetOutEDat(i).Val*parent.GetDat().Val));

				if(visitedNodes.count(iChildID) == 0)
				{
					visitedNodes[iChildID] = true;	// Mark the child
					queue.push(iChildID);
				}
			}
		}
		);
		else
			// Not enough elements to do the work concurrently
		{
			for(int i = 0; i < numChildren; ++i)
			{
				int iChildID = parent.GetOutNId(i);
				// Do p(u) = 1-(1-p(u))*(1-p(u,v)*p(v))
				pGraph->SetNDat(iChildID,
					1.0 - (1-parent.GetOutNDat(i).Val) * (1-parent.GetOutEDat(i).Val*parent.GetDat().Val));

				if(visitedNodes.count(iChildID) == 0)
				{
					visitedNodes[iChildID] = true;	// Mark the child
					queue.push(iChildID);
				}
			}
		}
	}
}

// http://sc05.supercomputing.org/schedule/pdf/pap346.pdf
// This version performs concurrently the BP on the nodes at the same level.
void ParallelBPFromNode_1DPartitioning_Precise(TPt<TNodeEDatNet<TFlt, TFlt>> pGraph, int sourceNodeID)
{
	// Like std::list, insertion of new items does not invalidate any iterators, nor change the order of items already in the map. Insertion and traversal may be concurrent.
	// Stores the number of time a node has been visited.
	tbb::concurrent_unordered_map<int, tbb::atomic<int>> visitedNodes;

	// If l is the level of the nodes being updated, then N is the list of node IDs at the level l+1 
	tbb::concurrent_queue<int> N;
	auto sourceNode = pGraph->GetNI(sourceNodeID);
	int numChildren = sourceNode.GetOutDeg();
	for(int i = 0; i < numChildren; ++i)
	{
		int childID = sourceNode.GetOutNId(i);
		N.push(childID);
		visitedNodes[childID] = 1;
	}

	while(true)
	{
		int numNodesToProcess = N.unsafe_size();
		// Parse level l
		if(numNodesToProcess > 50)
		{
			tbb::parallel_for(tbb::blocked_range<int>(0, numNodesToProcess, 50), 
				[&](const tbb::blocked_range<int>& r)
			{
				for (int i=r.begin();i!=r.end();++i)
				{
					int childID;
					N.try_pop(childID);
					auto child = pGraph->GetNI(childID);
					int numTimeVisited = visitedNodes[childID];
					if(numTimeVisited >= child.GetInDeg() && numTimeVisited < INT_MAX)
					{
						// Update child
						double dBelief = 1.0;
						int numParents = child.GetInDeg();

						for(int k = 0; k < numParents; ++k)
							dBelief = dBelief*(1-child.GetInEDat(k)*child.GetInNDat(k));
						pGraph->SetNDat(childID, 1.0-dBelief);
						visitedNodes[childID] = INT_MAX;	
					}

					int numNodesAtNextLevel = child.GetOutDeg();
					// Enqueue the layer of nodes at level l+1
					for(int j = 0; j < numNodesAtNextLevel; ++j)
					{
						int nextLevelNodeID = child.GetOutNId(j);
						auto it = visitedNodes.find(nextLevelNodeID);
						if(it == visitedNodes.end())
						{
							N.push(nextLevelNodeID);
							visitedNodes[nextLevelNodeID] = 1;
						}
						else
						{
							if(it->second < pGraph->GetNI(nextLevelNodeID).GetInDeg())
							{
								N.push(nextLevelNodeID);
								++(it->second);
							}
						}
					}
				}
			}
			);
		}
		else
			// Not enough elements to do the work concurrently
		{
			for (int i=0;i<numNodesToProcess;++i)
			{
				int childID;
				N.try_pop(childID);
				auto child = pGraph->GetNI(childID);
				int numTimeVisited = visitedNodes[childID];
				if(numTimeVisited >= child.GetInDeg() && numTimeVisited < INT_MAX)
				{
					// Update child
					double dBelief = 1.0;
					int numParents = child.GetInDeg();

					for(int k = 0; k < numParents; ++k)
						dBelief = dBelief*(1-child.GetInEDat(k)*child.GetInNDat(k));
					pGraph->SetNDat(childID, 1.0-dBelief);
					visitedNodes[childID] = INT_MAX;
				}

				int numNodesAtNextLevel = child.GetOutDeg();
				// Enqueue the layer of nodes at level l+1
				for(int j = 0; j < numNodesAtNextLevel; ++j)
				{
					int nextLevelNodeID = child.GetOutNId(j);
					auto it = visitedNodes.find(nextLevelNodeID);
					if(it == visitedNodes.end())
					{
						N.push(nextLevelNodeID);
						visitedNodes[nextLevelNodeID] = 1;
					}
					else
					{
						if(it->second < pGraph->GetNI(nextLevelNodeID).GetInDeg())
						{
							N.push(nextLevelNodeID);
							++(it->second);
						}
					}
				}
			}
		}

		if(N.empty())
			return;
	}
}

// http://sc05.supercomputing.org/schedule/pdf/pap346.pdf
// This version performs concurrently the BP on the nodes at the same level.
void ParallelBPFromNode_1DPartitioning(TPt<TNodeEDatNet<TFlt, TFlt>> pGraph, int sourceNodeID)
{
	// Like std::list, insertion of new items does not invalidate any iterators, nor change the order of items already in the map. Insertion and traversal may be concurrent.
	tbb::concurrent_unordered_map<int, bool> visitedNodes;
	visitedNodes[sourceNodeID] = true;

	// If l is the level of the nodes being updated, then N is the list of node IDs at the level l+1 
	tbb::concurrent_queue<int> N;
	auto sourceNode = pGraph->GetNI(sourceNodeID);
	int numChildren = sourceNode.GetOutDeg();
	for(int i = 0; i < numChildren; ++i)
		N.push(sourceNode.GetOutNId(i));

	while(true)
	{
		int numNodesToProcess = N.unsafe_size();
		// Parse level l
		if(numNodesToProcess > 50)
		{
			tbb::parallel_for(tbb::blocked_range<int>(0, numNodesToProcess, 50), 
			[&](const tbb::blocked_range<int>& r)
			{
				for (int i=r.begin();i!=r.end();++i)
				{
					int childID;
					N.try_pop(childID);
					auto child = pGraph->GetNI(childID);
					// Update child
					double dBelief = 1.0;
					int numParents = child.GetInDeg();

					for(int k = 0; k < numParents; ++k)
						dBelief = dBelief*(1-child.GetInEDat(k)*child.GetInNDat(k));
					pGraph->SetNDat(childID, 1.0-dBelief);

					int numNodesAtNextLevel = child.GetOutDeg();
					// Enqueue the layer of nodes at level l+1
					for(int j = 0; j < numNodesAtNextLevel; ++j)
					{
						int nextLevelNodeID = child.GetOutNId(j);
						// Note that this part is not completely threadsafe and that a same node maybe be pushed twice.
						// To make it completely threadsafe we would need to make the following if statement a critical section, which would kill
						// the purpose of using a concurrent hash map for visitedNodes
						if(visitedNodes.count(nextLevelNodeID) == 0)
						{
							N.push(nextLevelNodeID);
							visitedNodes[nextLevelNodeID] = true;
						}
					}
				}
			}
			);
		}
		else
		// Not enough elements to do the work concurrently
		{
			for (int i=0;i<numNodesToProcess;++i)
			{
				int childID;
				N.try_pop(childID);
				auto child = pGraph->GetNI(childID);
				// Update child
				double dBelief = 1.0;
				int numParents = child.GetInDeg();

				for(int k = 0; k < numParents; ++k)
					dBelief = dBelief*(1-child.GetInEDat(k)*child.GetInNDat(k));
				pGraph->SetNDat(childID, 1.0-dBelief);

				int numNodesAtNextLevel = child.GetOutDeg();
				for(int j = 0; j < numNodesAtNextLevel; ++j)
				{
					int nextLevelNodeID = child.GetOutNId(j);
					if(visitedNodes.count(nextLevelNodeID) == 0)
					{
						N.push(nextLevelNodeID);
						visitedNodes[nextLevelNodeID] = true;
					}
				}
			}
		}

		if(N.empty())
			return;
	}
}

//#define _TEST_Email_EuAll
#define _TEST_p2p_Gnutella09

//#define _LOAD_FROM_FILE
//#define _SAVE_TO_FILE
//#define _TEST_GRAPH

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

#ifdef _LOAD_FROM_FILE
	TFIn FIn("test.graph");
	auto pGraph = TNodeEDatNet<TFlt, TFlt>::Load(FIn);

	// Start traversing the graph
	cout << "Starting BP from nodeID 0. The input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";
	tbb::tick_count tic = tbb::tick_count::now();
	//PropagateFromNode(pGraph, 0);
	ParallelBPFromNode_1DPartitioning(pGraph, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed: " << dElapsedTime << " seconds\n";
#endif

#ifdef _SAVE_TO_FILE
	auto pGraph = GenerateRandomBayesianNetwork(1, 5, 5, 7, 30);
	TFOut FOut("test.graph");
	pGraph->Save(FOut);
	TSnap::SaveGViz(pGraph, "test.gv", "Test DAG", true);
#endif

#ifdef _TEST_p2p_Gnutella09
	auto pGraph = TSnap::LoadEdgeList<TPt<TNodeEDatNet<TFlt, TFlt>>>("p2p-Gnutella09.txt", 0, 1);
	// Start traversing the graph
	cout << "Starting BP from nodeID 0. The input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";
	tbb::tick_count tic = tbb::tick_count::now();
	for(int i = 0; i<100; ++i)
		ParallelBPFromNode(pGraph, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel BP: " << dElapsedTime/100 << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<100; ++i)
		ParallelBPFromNode_1DPartitioning(pGraph, 0);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel 1D partitioning BP: " << dElapsedTime/100 << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<100; ++i)
		PropagateFromNode(pGraph, 0);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for serial BP: " << dElapsedTime/100 << " seconds\n";
#endif

#ifdef _TEST_Email_EuAll
	auto pGraph = TSnap::LoadEdgeList<TPt<TNodeEDatNet<TFlt, TFlt>>>("Email-EuAll.txt", 0, 1);
	// Start traversing the graph
	cout << "Starting BP from nodeID 0. The input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";
	tbb::tick_count tic = tbb::tick_count::now();
	ParallelBPFromNode(pGraph, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel BP: " << dElapsedTime << " seconds\n";

	tic = tbb::tick_count::now();
	ParallelBPFromNode_1DPartitioning(pGraph, 0);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel 1D partitioning BP: " << dElapsedTime << " seconds\n";

	tic = tbb::tick_count::now();
	PropagateFromNode(pGraph, 0);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for serial BP: " << dElapsedTime << " seconds\n";
#endif

#ifdef _WIN32
	system("pause");
#endif
#ifdef __linux__
	SystemPause();
#endif

	return 0;
}
