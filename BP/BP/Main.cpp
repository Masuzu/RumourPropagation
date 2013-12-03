#include "BP.h"
#include "Utilities.h"

#define __TBB_NO_IMPLICIT_LINKAGE	1
#include <tbb/tick_count.h>

#include <iostream>

#include <Snap.h>

using namespace std;

// Discussion about variable declaration and loops
// http://stackoverflow.com/questions/982963/is-there-any-overhead-to-declaring-a-variable-within-a-loop-c

//#define _TEST_webStanford
//x#define _TEST_Amazon0302
//#define _TEST_Slashdot0902
//#define _TEST_Email_EuAll
#define _TEST_p2p_Gnutella04

//#define _LOAD_FROM_FILE
//#define _SAVE_TO_FILE
//#define _TEST_GRAPH
//#define _TEST_DAG2

void TestGraph(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, std::vector<int> vSeedIDs, int numIterations)
{
	cout << "Starting BP from nodeID {";
	for(int sourceNode : vSeedIDs)
		cout << sourceNode << " ";
	cout << "}\nThe input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";

	// Start traversing the graph
	tbb::tick_count tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNode(pGraph, vSeedIDs);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNode_LevelSynchronous(pGraph, vSeedIDs);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel level synchronous BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		PropagateFromNode(pGraph, vSeedIDs);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for serial BP: " << dElapsedTime/numIterations << " seconds\n";

	std::vector<int> vResult;
	CalculateRankFromSource_BellmanFord(pGraph, vSeedIDs, vResult);
	tic = tbb::tick_count::now();	
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNode_SingleNodeUpdate(pGraph, vResult, vSeedIDs);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for ParallelBPFromNode_SingleNodeUpdate: " << dElapsedTime/numIterations << " seconds\n";
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
		ParallelBPFromNode_LevelSynchronous(pGraph, sourceNode);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel level synchronous BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		PropagateFromNode(pGraph, sourceNode);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for serial BP: " << dElapsedTime/numIterations << " seconds\n";

	std::vector<int> vResult;
	CalculateRankFromSource(pGraph, sourceNode, vResult);
	tic = tbb::tick_count::now();	
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNode_SingleNodeUpdate(pGraph, vResult, sourceNode);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for ParallelBPFromNode_SingleNodeUpdate: " << dElapsedTime/numIterations << " seconds\n";
}

void TestGraphDAG(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, std::vector<int> vSeedIDs, int numIterations)
{
	cout << "Starting BP from nodeID {";
	for(int sourceNode : vSeedIDs)
		cout << sourceNode << " ";
	cout << "}\nThe input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";

	// Start traversing the graph
	tbb::tick_count tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNodeDAG(pGraph, vSeedIDs);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNodeDAG_LevelSynchronous(pGraph, vSeedIDs);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel level synchronous BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		PropagateFromNodeDAG(pGraph, vSeedIDs);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for serial BP: " << dElapsedTime/numIterations << " seconds\n";

	std::vector<int> vResult;
	CalculateRankFromSource(pGraph, vSeedIDs, vResult);
	tic = tbb::tick_count::now();	
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNode_SingleNodeUpdate(pGraph, vResult, vSeedIDs);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for ParallelBPFromNode_SingleNodeUpdate: " << dElapsedTime/numIterations << " seconds\n";
}

void TestGraphDAG(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNode, int numIterations)
{
	// Start traversing the graph
	cout << "Starting BP from nodeID " << sourceNode << ". The input graph has " << pGraph->GetNodes() << " nodes and " << pGraph->GetEdges() << " edges.\n";
	tbb::tick_count tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNodeDAG(pGraph, sourceNode);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNodeDAG_LevelSynchronous(pGraph, sourceNode);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for parallel level synchronous BP: " << dElapsedTime/numIterations << " seconds\n";

	tic = tbb::tick_count::now();
	for(int i = 0; i<numIterations; ++i)
		PropagateFromNodeDAG(pGraph, sourceNode);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for serial BP: " << dElapsedTime/numIterations << " seconds\n";

	//auto rankGraph = CalculateRankFromSource(pGraph, sourceNode);
	std::vector<int> vResult;
	CalculateRankFromSource(pGraph, sourceNode, vResult);
	tic = tbb::tick_count::now();	
	for(int i = 0; i<numIterations; ++i)
		ParallelBPFromNode_SingleNodeUpdate(pGraph, vResult, sourceNode);
	dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for ParallelBPFromNode_SingleNodeUpdate: " << dElapsedTime/numIterations << " seconds\n";
}

#ifdef _USE_CUDA
extern void MatrixMulOnDevice(float * M, float * N, float * P, int Width);
#endif

int main(int argc, char* argv[])
{

#ifdef _USE_CUDA
	float M[4] = {1, 1, 1, 1};
	float N[4] = {1, 1, 1, 1};
	float P[4];
	MatrixMulOnDevice(M, N, P, 2);
#endif

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
	ParallelBPFromNode_LevelSynchronous(pGraph, 0);
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
	pGraph = GenerateDAG2(pGraph, vSeedNodeIDs, 0);
	//pGraph = MIOA(pGraph, 0, 0);
	TSnap::SaveGViz(pGraph, "test.gv", "Test DAG", true);
#endif

#ifdef _LOAD_FROM_FILE
	TFIn FIn("test.graph");
	auto pGraph = TNodeEDatNet<TFlt, TFlt>::Load(FIn);

	//for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
		//std::cout << EI.GetSrcNId() << " " << EI.GetDstNId() << " " << EI.GetDat() << std::endl;

	/*
	std::cout << "ExactBP_Marginalization(pGraph, 3)\n";
	ExactBP_Marginalization(pGraph, 3);
	auto ExactBPGraph = CopyGraph(pGraph);
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
		std::cout << NI.GetId() << " " << NI.GetDat() << std::endl;

	ResetGraphBelief(pGraph);
	std::cout << "PropagateFromNode(pGraph, 3)\n";
	*/

	std::vector<int> mapRanks;
	CalculateRankFromSource(pGraph, 3, mapRanks);

	ParallelBPFromNode_SingleNodeUpdate(pGraph, mapRanks, 3);
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
		std::cout << NI.GetId() << " " << NI.GetDat() << std::endl;

	//std::cout << "Error : " << BPError(ExactBPGraph, pGraph, [] (double a, double b) -> double { return abs(a-b); })  << std::endl;

	/*
	std::vector<int> vSeedNodes;
	vSeedNodes.push_back(0);
	vSeedNodes.push_back(3);
	pGraph = GenerateDAG2(pGraph, vSeedNodes);
	TSnap::SaveGViz(pGraph, "testDAG.gv", "Test DAG", true);
	*/
#endif

#ifdef _SAVE_TO_FILE
	auto pGraph = GenerateRandomBayesianNetwork(1, 5, 5, 7, 30);
	RandomGraphInitialization(pGraph);
	TFOut FOut("test.graph");
	pGraph->Save(FOut);
	TSnap::SaveGViz(pGraph, "test.gv", "Test DAG", true);
#endif

#ifdef _TEST_p2p_Gnutella04
	auto pGraph = TSnap::LoadEdgeList<TPt<TNodeEDatNet<TFlt, TFlt>>>("p2p-Gnutella04.txt", 0, 1);
	LoadEdgeWeightsFromFile(pGraph, "p2p-Gnutella04-Network.txt");

	std::vector<int> vSeedNodeIDs;
	vSeedNodeIDs.push_back(0);
	vSeedNodeIDs.push_back(11);
	pGraph = GenerateDAG1(pGraph, vSeedNodeIDs, 0);
	//pGraph = GenerateDAG2(pGraph, vSeedNodeIDs, 0);

	/*
	std::cout << "ExactBP_Marginalization(pGraph, vSeedNodeIDs)\n";	// Crash out of memory
	ExactBP_Marginalization(pGraph, vSeedNodeIDs);
	auto ExactBPGraph = CopyGraph(pGraph);

	ResetGraphBelief(pGraph);
	std::cout << "PropagateFromNode(pGraph, vSeedNodeIDs)\n";
	PropagateFromNode(pGraph, vSeedNodeIDs);

	std::cout << "Error : " << BPError(ExactBPGraph, pGraph, [] (double a, double b) -> double { return abs(a-b); })  << std::endl;
	*/

	//pGraph = GenerateDAG2(pGraph, 0, 0);
	//ParallelBPFromNode(pGraph, vSeedNodeIDs);
	TestGraph(pGraph, vSeedNodeIDs, 100);
#endif

#ifdef _TEST_Email_EuAll

	auto pGraph = TSnap::LoadEdgeList<TPt<TNodeEDatNet<TFlt, TFlt>>>("Email-EuAll.txt", 0, 1);
	LoadEdgeWeightsFromFile(pGraph, "Email-EuAll-Network.txt");

#if TEST_EDGE_WEIGHT_SAVING
	SaveEdgeWeightsToFile(pGraph, "BLABLA.txt");

	// ...

	LoadEdgeWeightsFromFile(pGraph, "BLABLA.txt");
#endif

	// ...

	std::vector<int> vSeedNodeIDs;
	vSeedNodeIDs.push_back(0);
	vSeedNodeIDs.push_back(3);

	tbb::tick_count tic = tbb::tick_count::now();
	pGraph = GenerateDAG1(pGraph, vSeedNodeIDs, 0);
	//pGraph = GenerateDAG2(pGraph, 0, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for DAG2 computation: " << dElapsedTime << " seconds\n";
	
	TestGraph(pGraph, vSeedNodeIDs, 1);

	//TestGraphDAG(pGraph, 0, 100);	// Not viable (memory overflow)

#endif

#ifdef _TEST_Slashdot0902
	auto pGraph = TSnap::LoadConnList<TPt<TNodeEDatNet<TFlt, TFlt>>>("Slashdot0902.txt");
	LoadEdgeWeightsFromFile(pGraph, "Slashdot0902-Network.txt");

	tbb::tick_count tic = tbb::tick_count::now();
	pGraph = GenerateDAG1(pGraph, 0, 0);
	//pGraph = GenerateDAG2(pGraph, 0, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for DAG2 computation: " << dElapsedTime << " seconds\n";

	TestGraph(pGraph, 0, 1);
#endif

#ifdef _TEST_webStanford
	auto pGraph = TSnap::LoadConnList<TPt<TNodeEDatNet<TFlt, TFlt>>>("web-Stanford.txt");
	LoadEdgeWeightsFromFile(pGraph, "web-Stanford-Network.txt");

	tbb::tick_count tic = tbb::tick_count::now();
	pGraph = GenerateDAG1(pGraph, 2, 0);
	//pGraph = GenerateDAG2(pGraph, 0, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for DAG2 computation: " << dElapsedTime << " seconds\n";

	//TestGraph(pGraph, 2, 1);
	ParallelBPFromNode_LevelSynchronous(pGraph, 2);
#endif

#ifdef _TEST_Amazon0302
	auto pGraph = TSnap::LoadConnList<TPt<TNodeEDatNet<TFlt, TFlt>>>("Amazon0302.txt");
	LoadEdgeWeightsFromFile(pGraph, "Amazon0302-Network.txt");

	tbb::tick_count tic = tbb::tick_count::now();
	pGraph = GenerateDAG1(pGraph, 0, 0);
	//pGraph = GenerateDAG2(pGraph, 0, 0);
	double dElapsedTime = (tbb::tick_count::now() - tic).seconds();
	cout << "Time elapsed for DAG2 computation: " << dElapsedTime << " seconds\n";

	TestGraph(pGraph, 0, 1);
#endif

#ifdef _WIN32
	system("pause");
#endif
#ifdef __linux__
	SystemPause();
#endif

	return 0;
}
