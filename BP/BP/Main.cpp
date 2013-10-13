#include "BP.h"
#include "Utilities.h"

#define __TBB_NO_IMPLICIT_LINKAGE	1
#include <tbb/tick_count.h>

#include <iostream>

#include <Snap.h>

using namespace std;

//#define _TEST_Email_EuAll
//#define _TEST_p2p_Gnutella09

#define _LOAD_FROM_FILE
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
