#include "BP.h"
#include "Utilities.h"

#define __TBB_NO_IMPLICIT_LINKAGE	1
#include <tbb/tick_count.h>

#include <iostream>

#include <Snap.h>


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

void TestGraph(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, std::vector<int> vSeedIDs, int numIterations)
{

	cout << "Starting BP from nodeID {";
	for(int sourceNode : vSeedIDs)
		cout << sourceNode << " ";
	cout << "}\nThe input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";

	// Start traversing the graph
	tbb::tick_count tic = tbb::tick_count::now();
	for(int sourceNode : vSeedIDs)
		for(int i = 0; i<numIterations; ++i)
			ParallelBPFromNode(pGraph, sourceNode);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int sourceNode : vSeedIDs)
		for(int i = 0; i<numIterations; ++i)
			ParallelBPFromNode_1DPartitioning(pGraph, sourceNode);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel 1D partitioning BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int sourceNode : vSeedIDs)
		for(int i = 0; i<numIterations; ++i)
			PropagateFromNode(pGraph, sourceNode);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for serial BP: " << dElapsedTime/numIterations << " seconds\n";
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
	for(int i=0;i<5;++i)
		pGraph->AddNode(i);

	pGraph->AddEdge(0,2);
	pGraph->AddEdge(0,3);
	pGraph->AddEdge(1,2);
	pGraph->AddEdge(1,3);
	pGraph->AddEdge(2,3);
	pGraph->AddEdge(2,4);
	pGraph->AddEdge(3,4);

	pGraph->SetEDat(0,2,0.5);
	pGraph->SetEDat(0,3,0.4);
	pGraph->SetEDat(1,2,0.3);
	pGraph->SetEDat(1,3,0.4);
	pGraph->SetEDat(2,3,0.5);
	pGraph->SetEDat(2,4,0.4);
	pGraph->SetEDat(3,4,0.5);

	std::vector<int> vSeedNodeIDs;
	vSeedNodeIDs.push_back(0);
	vSeedNodeIDs.push_back(1);
	pGraph = DAG2(pGraph, vSeedNodeIDs, 0);
	//pGraph = MIOA(pGraph, 0, 0);
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
	
	RandomGraphInitialization(pGraph);

	std::vector<int> vSeedNodeIDs;
	vSeedNodeIDs.push_back(0);
	vSeedNodeIDs.push_back(11);
	vSeedNodeIDs.push_back(21);
	pGraph = DAG2(pGraph, vSeedNodeIDs, 0);
	//pGraph = DAG2(pGraph, 0, 0);
	TestGraph(pGraph, vSeedNodeIDs, 100);
#endif

#ifdef _TEST_Email_EuAll
	auto pGraph = TSnap::LoadEdgeList<TPt<TNodeEDatNet<TFlt, TFlt>>>("Email-EuAll.txt", 0, 1);

	std::vector<int> vSeedNodeIDs;
	vSeedNodeIDs.push_back(0);
	vSeedNodeIDs.push_back(2);
	vSeedNodeIDs.push_back(6);

	RandomGraphInitialization(pGraph);

	tbb::tick_count tic = tbb::tick_count::now();
	pGraph = DAG2(pGraph, vSeedNodeIDs, 0);
	//pGraph = DAG2(pGraph, 0, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for DAG2 computation: " << dElapsedTime << " seconds\n";
	
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
