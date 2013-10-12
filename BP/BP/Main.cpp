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

class WorkerThread
{
protected:
	const int m_iThreadID;

public:
	WorkerThread(int id) : m_iThreadID(id)	{ }

	void operator()()
	{

	}
};

// http://snap.stanford.edu/snap/quick.html
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>> pGraph, int sourceNodeID, bool bDisplayInfo = false)
{
	// Used to store the nodes which have been traversed during the breadth-first search traversal
	std::map<int, bool> visitedNodes;
	std::queue<int> queue;
	queue.push(sourceNodeID);
	visitedNodes[sourceNodeID] = true;
	int numEdge = 0;
	while(!queue.empty())
	{
		int nodeID = queue.front();
		if(bDisplayInfo)
			std::cout << "Parsing node " << nodeID << " - Number of edges traversed: " << ++numEdge << std::endl;

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

struct ParallelBPFromNodeBody
{
	TPt<TNodeEDatNet<TFlt, TFlt>> *m_pGraph;
	TNodeEDatNet<TFlt, TFlt>::TNodeI *m_pParent;
	tbb::concurrent_queue<int> *m_pQueue;
	tbb::concurrent_unordered_map<int, bool> *m_pVisitedNodes;

	ParallelBPFromNodeBody(TPt<TNodeEDatNet<TFlt, TFlt>> *pGraph, tbb::concurrent_queue<int> *pQueue, tbb::concurrent_unordered_map<int, bool> *pVisitedNodes, TNodeEDatNet<TFlt, TFlt>::TNodeI *pParent)
		: m_pGraph(pGraph), m_pQueue(pQueue), m_pVisitedNodes(pVisitedNodes), m_pParent(pParent)
	{

	}

	ParallelBPFromNodeBody()	{}

	void operator ()(const tbb::blocked_range<int>& r)	const
	{
		for (int i=r.begin();i!=r.end();++i)
		{
			int iChildID = m_pParent->GetOutNId(i);
			// Do p(u) = 1-(1-p(u))*(1-p(u,v)*p(v))
			(*m_pGraph)->SetNDat(iChildID,
				1.0 - (1-m_pParent->GetOutNDat(i).Val) * (1-m_pParent->GetOutEDat(i).Val*m_pParent->GetDat().Val));

			if(m_pVisitedNodes->count(iChildID) == 0)
			{
				(*m_pVisitedNodes)[iChildID] = true;	// Mark the child
				m_pQueue->push(iChildID);
			}

		}

	}
};

// http://sc05.supercomputing.org/schedule/pdf/pap346.pdf
void ParallelBPFromNode(TPt<TNodeEDatNet<TFlt, TFlt>> pGraph, int sourceNodeID)
{
	// Like std::list, insertion of new items does not invalidate any iterators, nor change the order of items already in the map. Insertion and traversal may be concurrent.
	tbb::concurrent_unordered_map<int, bool> visitedNodes;

	tbb::concurrent_queue<int> queue;
	queue.push(sourceNodeID);
	visitedNodes[sourceNodeID] = true;

	int numEdge = 0;
	int nodeID;
	while(queue.try_pop(nodeID))
	{
		auto parent = pGraph->GetNI(nodeID);
		int numChildren = parent.GetOutDeg();
		if(numChildren > 10)
			tbb::parallel_for(tbb::blocked_range<int>(0,numChildren), ParallelBPFromNodeBody(&pGraph, &queue, &visitedNodes, &parent));
		else
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

//#define _TEST_Email_EuAll
#define _TEST_p2p_Gnutella09

// #define _SAVE_TO_FILE
// #define _TEST_GRAPH

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
	ParallelBPFromNode(pGraph, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	std:: cout << pGraph->GetNDat(0).Val << " " << pGraph->GetNDat(1).Val << " " << pGraph->GetNDat(2).Val << std::endl;
#endif

#ifdef _LOAD_FROM_FILE
	TFIn FIn("test.graph");
	auto pGraph = TNodeEDatNet<TFlt, TFlt>::Load(FIn);

	// Start traversing the graph
	cout << "Starting BP from nodeID 0. The input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";
	tbb::tick_count tic = tbb::tick_count::now();
	PropagateFromNode(pGraph, 0);
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
	ParallelBPFromNode(pGraph, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel BP: " << dElapsedTime << " seconds\n";

	tic = tbb::tick_count::now();
	PropagateFromNode(pGraph, 0, false);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for serial BP: " << dElapsedTime << " seconds\n";
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
	PropagateFromNode(pGraph, 0, false);
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
