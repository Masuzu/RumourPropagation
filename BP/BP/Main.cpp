#include "Snap.h"
#include <queue>
#include <iostream>

//! @param minNumNodesPerRank How 'fat' the DAG should be.
//! @param minRanks How 'tall' the DAG should be.
//! @param probabilityEdge Chance of having an Edge in percent.
TPt<TNodeEDatNet<TFlt, TFlt>> GenerateRandomBayesianNetwork(UINT minNumNodesPerRank, UINT maxNumNodesPerRank, UINT minRanks, UINT maxRanks, UINT probabilityEdge)
{
	srand (time (NULL));
	auto pGraph = TNodeEDatNet<TFlt, TFlt>::New();
	int nodes = 0;
	int ranks = minRanks + (rand () % (maxRanks - minRanks + 1));

	for (int i = 0; i < ranks; i++)
	{
		/* New nodes of 'higher' rank than all nodes generated till now.  */
		int new_nodes = minNumNodesPerRank + (rand () % (maxNumNodesPerRank - minNumNodesPerRank + 1));

		/* Edges from old nodes ('nodes') to new ones ('new_nodes').  */
		for (int j = 0; j < nodes; ++j)
			for (int k = 0; k < new_nodes; ++k)
				if ( (rand () % 100) < probabilityEdge)
				{
					if(!pGraph->IsNode(j))
					{
						pGraph->AddNode(j);
						pGraph->SetNDat(j, (double)rand() / RAND_MAX);
					}
					if(!pGraph->IsNode(k+nodes))
					{
						pGraph->AddNode(k+nodes);
						pGraph->SetNDat(k+nodes, (double)rand() / RAND_MAX);
					}
					pGraph->AddEdge(j,k+nodes);
					pGraph->SetEDat(j,k+nodes, (double)rand() / RAND_MAX);
				}

		nodes += new_nodes; /* Accumulate into old node set.  */
	}
	return pGraph;
}

// http://snap.stanford.edu/snap/quick.html
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>> pGraph, int sourceNodeID)
{
	std::queue<int> queue;
	queue.push(sourceNodeID);
	while(!queue.empty())
	{
		int nodeID = queue.front();
		std::cout << "Parsing node " << nodeID << std::endl;

		auto parent = pGraph->GetNI(nodeID);	
		int numChildren = parent.GetOutDeg();
		// Update the belief of the children of parent
		for(int i = 0; i < numChildren; ++i)
		{
			// Do p(u) = 1-(1-p(u))*(1-p(u,v)*p(v))
			pGraph->SetNDat(parent.GetOutNId(i),
				1.0 - (1-parent.GetOutNDat(i).Val) * (1-parent.GetOutEDat(i).Val)*parent.GetDat().Val);
			queue.push(parent.GetOutNId(i));
		}

		queue.pop();
	}
}

#define _TEST_GRAPH
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
	PropagateFromNode(pGraph, 0);
	std:: cout << pGraph->GetNDat(0).Val << " " << pGraph->GetNDat(1).Val << " " << pGraph->GetNDat(2).Val << std::endl;
#endif
#ifdef _LOAD_FROM_FILE
	TFIn FIn("test.graph");
	auto pGraph = TNodeEDatNet<TFlt, TFlt>::Load(FIn);
	PropagateFromNode(pGraph, 0);
#endif
#ifdef _SAVE_TO_FILE
	auto pGraph = GenerateRandomBayesianNetwork(1, 5, 3, 5, 30);
	TFOut FOut("test.graph");
	pGraph->Save(FOut);
	TSnap::SaveGViz(pGraph, "test.gv", "Test DAG", "20");
#endif

	return 0;
}


