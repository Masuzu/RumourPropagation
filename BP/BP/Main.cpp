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

		auto it = pGraph->GetNI(nodeID);	
		int numChildren = it.GetOutDeg();
		for(int i = 0; i < numChildren; ++i)
			queue.push(it.GetOutNId(i));

		// Get the parents of the child
		int numParents = it.GetInDeg();
		for(int i = 0; i < numParents; ++i)
		{
			it.GetInNDat(i);
		}
		queue.pop();
	}
}

#define _SAVE_TO_FILE
int main(int argc, char* argv[])
{
#ifdef _LOAD_FROM_FILE
	auto pGraph = TSnap::LoadEdgeList<TPt<TNodeEDatNet<TFlt, TFlt>>>("test.txt", 0, 1);
	PropagateFromNode(pGraph, 0);
#endif
#ifdef _SAVE_TO_FILE
	auto pGraph = GenerateRandomBayesianNetwork(1, 5, 3, 5, 30);
	/*<
	TFOut FOut("test.graph");
	pGraph->Save(FOut);
	*/
	TSnap::SaveEdgeList(pGraph, "test.txt", "Save as tab-separated list of edges");
	TSnap::SaveGViz(pGraph, "test.gv", "Test DAG", "20");
#endif

	return 0;
}


