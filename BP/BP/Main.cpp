#include <queue>
#include <iostream>
#include <vector>
#include <boost/thread/thread.hpp>
#include <boost/timer.hpp>
#include "Utilities.h"
#include "Snap.h"

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
	TSnap::SaveEdgeList(pGraph, "test.txt", "Save as tab-separated list of edges");
	auto pGraphTraversed =  TSnap::LoadEdgeList<TPt<TNodeEDatNet<TInt, TInt>>>("test.txt", 0, 1);
	for (auto NI = pGraphTraversed->BegNI(); NI < pGraphTraversed->EndNI(); NI++)
		pGraphTraversed->SetNDat(NI.GetId(), 0);

	std::queue<int> queue;
	queue.push(sourceNodeID);
	pGraphTraversed->SetNDat(sourceNodeID, 1);
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
				1.0 - (1-parent.GetOutNDat(i).Val) * (1-parent.GetOutEDat(i).Val)*parent.GetDat().Val);
			
			if(pGraphTraversed->GetNDat(iChildID).Val == 0)
			{
				pGraphTraversed->SetNDat(iChildID, 1);	// Mark the child
				queue.push(iChildID);
			}
		}

		queue.pop();
	}
}

#define _TEST_Email_EuAll
// #define _TEST_p2p_Gnutella09

// #define _LOAD_FROM_FILE
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
	boost::timer t; // start timing
	PropagateFromNode(pGraph, 0);
	double dElapsedTime = t.elapsed();
	std:: cout << pGraph->GetNDat(0).Val << " " << pGraph->GetNDat(1).Val << " " << pGraph->GetNDat(2).Val << std::endl;
#endif

#ifdef _LOAD_FROM_FILE
	TFIn FIn("test.graph");
	auto pGraph = TNodeEDatNet<TFlt, TFlt>::Load(FIn);
	
	// Start traversing the graph
	cout << "Starting BP from nodeID 0. The input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";
	boost::timer t; // start timing
	PropagateFromNode(pGraph, 0);
	double dElapsedTime = t.elapsed();
	cout << "Time elapsed: " << dElapsedTime << " seconds\n";
#endif

#ifdef _SAVE_TO_FILE
	auto pGraph = GenerateRandomBayesianNetwork(5, 30, 50, 60, 30);
	TFOut FOut("test.graph");
	pGraph->Save(FOut);
	TSnap::SaveGViz(pGraph, "test.gv", "Test DAG");
#endif

#ifdef _TEST_p2p_Gnutella09
	auto pGraph = TSnap::LoadEdgeList<TPt<TNodeEDatNet<TFlt, TFlt>>>("p2p-Gnutella09.txt", 0, 1);
	// Start traversing the graph
	cout << "Starting BP from nodeID 0. The input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";
	boost::timer t; // start timing
	PropagateFromNode(pGraph, 0, true);
	double dElapsedTime = t.elapsed();
	cout << "Time elapsed: " << dElapsedTime << " seconds\n";
#endif
	 
	#ifdef _TEST_Email_EuAll
	auto pGraph = TSnap::LoadEdgeList<TPt<TNodeEDatNet<TFlt, TFlt>>>("Email-EuAll.txt", 0, 1);
	// Start traversing the graph
	cout << "Starting BP from nodeID 0. The input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";
	boost::timer t; // start timing
	PropagateFromNode(pGraph, 0, false);
	double dElapsedTime = t.elapsed();
	cout << "Time elapsed: " << dElapsedTime << " seconds\n";
#endif

	/*
	vector<boost::thread> vWorkerThreads;
	for(auto it = vWorkerThreads.begin(); it != vWorkerThreads.end(); ++it)
		it->join();
		*/

#ifdef _WIN32
	system("pause");
#endif
#ifdef __linux__
	SystemPause();
#endif

	return 0;
}
