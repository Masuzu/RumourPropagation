#include "Utilities.h"
#include <Snap.h>
#include <map>
#include <algorithm>

#ifdef WIN32

#include <Windows.h>

//--------------------------------------------------------------------------------------
// Helper functions for querying information about the processors in the current
// system.  ( Copied from the doc page for GetLogicalProcessorInformation() )
//--------------------------------------------------------------------------------------
typedef BOOL (WINAPI *LPFN_GLPI)(
	PSYSTEM_LOGICAL_PROCESSOR_INFORMATION, 
	PDWORD);


//  Helper function to count bits in the processor mask
static DWORD CountBits(ULONG_PTR bitMask)
{
	DWORD LSHIFT = sizeof(ULONG_PTR)*8 - 1;
	DWORD bitSetCount = 0;
	DWORD bitTest = 1 << LSHIFT;
	DWORD i;

	for( i = 0; i <= LSHIFT; ++i)
	{
		bitSetCount += ((bitMask & bitTest)?1:0);
		bitTest/=2;
	}

	return bitSetCount;
}

int GetPhysicalProcessorCount()
{
	DWORD procCoreCount = 0;    // Return 0 on any failure.  That'll show them.

	LPFN_GLPI Glpi;

	Glpi = (LPFN_GLPI) GetProcAddress(
		GetModuleHandle(TEXT("kernel32")),
		"GetLogicalProcessorInformation");
	if (NULL == Glpi) 
	{
		// GetLogicalProcessorInformation is not supported
		return procCoreCount;
	}

	BOOL done = FALSE;
	PSYSTEM_LOGICAL_PROCESSOR_INFORMATION buffer = NULL;
	DWORD returnLength = 0;

	while (!done) 
	{
		DWORD rc = Glpi(buffer, &returnLength);

		if (FALSE == rc) 
		{
			if (GetLastError() == ERROR_INSUFFICIENT_BUFFER) 
			{
				if (buffer) 
					free(buffer);

				buffer=(PSYSTEM_LOGICAL_PROCESSOR_INFORMATION)malloc(
					returnLength);

				if (NULL == buffer) 
				{
					// Allocation failure\n
					return procCoreCount;
				}
			} 
			else 
			{
				// Unanticipated error
				return procCoreCount;
			}
		} 
		else done = TRUE;
	}

	DWORD byteOffset = 0;
	PSYSTEM_LOGICAL_PROCESSOR_INFORMATION ptr = buffer;
	while (byteOffset < returnLength) 
	{
		if (ptr->Relationship == RelationProcessorCore) 
		{
			if(ptr->ProcessorCore.Flags)
			{
				//  Hyperthreading or SMT is enabled.
				//  Logical processors are on the same core.
				procCoreCount += 1;
			}
			else
			{
				//  Logical processors are on different cores.
				procCoreCount += CountBits(ptr->ProcessorMask);
			}
		}
		byteOffset += sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
		ptr++;
	}

	free (buffer);

	return procCoreCount;
}

#endif

#ifdef __linux__
#include "stdio.h"

int GetPhysicalProcessorCount()
{
	FILE * fp;
	char res[128];
	fp = popen("/bin/cat /proc/cpuinfo |grep -c '^processor'","r");
	fread(res, 1, sizeof(res)-1, fp);
	fclose(fp);
	return res[0];
}

#ifdef __linux__
void SystemPause()
{
	char magickey;
	fflush(stdin);
	printf("Appuyez sur une touche pour continuer...");
	scanf("%c", &magickey);
	magickey = getc(stdin);
}
#endif
#endif

TPt<TNodeEDatNet<TFlt, TFlt>> GenerateRandomBayesianNetwork(unsigned int minNumNodesPerRank, unsigned int maxNumNodesPerRank, unsigned int minRanks, unsigned int maxRanks, unsigned int probabilityEdge)
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

void RandomGraphInitialization(TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph)
{
	srand(time(NULL));
	for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
		pGraph->SetEDat(EI.GetSrcNId(), EI.GetDstNId(), (double) rand() / RAND_MAX);
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
		pGraph->SetNDat(NI.GetId(), (double) rand() / RAND_MAX);
}

///////////
// DAG 2 //
///////////

//! Given pGraph with data about edge weights, computes the distance of the shortest paths from sourceNode
//! and returns the result in the nodes of pDAGGraph.
//! Updates the edges if bUpdateEdges is set to true. Default is false.
//! @note Requires initial values for the nodes of pDAGGraph (edges are not needed)
void Dijkstra(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNode, double dThreshold, TPt<TNodeEDatNet<TFlt, TFlt>>& pDAGGraph, bool bUpdateEdges = false)
{
	double logThreshold = log(dThreshold);
	if(dThreshold==0)
		logThreshold=-DBL_MAX;

	// List of visited nodes
	std::map<int, bool> visitedNodes;
	// Stores the edge vertices to build the final DAG
	std::map<int, int> mapPrevious;
	std::vector<int> vNonVisitedNodes;

	// Distance from source node to itself is 0
	pDAGGraph->SetNDat(sourceNode, 0);
	vNonVisitedNodes.push_back(sourceNode);

	// Beginning of the loop of Dijkstra algorithm

	while(!vNonVisitedNodes.empty())
	{
		// Find the vertex in queue with the smallest distance and remove it
		auto itVertexWithSmallestDistance = std::min_element(vNonVisitedNodes.begin(), vNonVisitedNodes.end());
		int iParentID = *itVertexWithSmallestDistance;
		vNonVisitedNodes.erase(itVertexWithSmallestDistance);
		visitedNodes.insert(std::make_pair(iParentID, true));
		auto parent = pGraph->GetNI(iParentID);

		int numChildren = parent.GetOutDeg();
		
		for(int i = 0; i < numChildren; ++i)
		{
			int iChildID = parent.GetOutNId(i);
			// Accumulate the shortest distance from source
			double alt = pDAGGraph->GetNDat(iParentID) - log(parent.GetOutEDat(i).Val);
			if(alt >= logThreshold)
			{
				auto it = visitedNodes.find(iChildID);
				if (alt < pDAGGraph->GetNDat(iChildID) && it == visitedNodes.end())
				{
					pDAGGraph->SetNDat(iChildID, alt);
					if(bUpdateEdges)
					{
						auto itPrevious = mapPrevious.find(iChildID);
						if(itPrevious==mapPrevious.end())
							mapPrevious.insert(std::make_pair(iChildID, iParentID));
						else
							itPrevious->second = iParentID;
					}
					// Add non-visited iChildID into vNonVisitedNodes to be processed
					vNonVisitedNodes.push_back(iChildID);                          
				}
			}
		}

	}

	if(bUpdateEdges)
		for(auto it=mapPrevious.begin(); it!= mapPrevious.end(); ++it)
			pDAGGraph->AddEdge(it->second, it->first);
}

//! Input : Directed graph with initialized weight edges.
//! Edges with a propagation probability strictly greater than dThreshold are ignored
TPt<TNodeEDatNet<TFlt, TFlt>> MIOA(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNode, double dThreshold)
{
	//////////////////////////////////////////////////////////////
	// Compte the Maximum Influence Out-Arborescence with Dijkstra

	// Copy the nodes of pGraph
	auto pDAGGraph = TNodeEDatNet<TFlt, TFlt>::New();
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
	{
		int NodeID = NI.GetId();
		pDAGGraph->AddNode(NodeID);
		pDAGGraph->SetNDat(NodeID, FLT_MAX);
	}

	Dijkstra(pGraph, sourceNode, dThreshold, pDAGGraph, true);

	// pDAGGraph is the MIOA starting from sourceNode

	return pDAGGraph;
}

TPt<TNodeEDatNet<TFlt, TFlt>> GraphUnion(const std::vector<TPt<TNodeEDatNet<TFlt, TFlt>>> &vGraphs)
{
	auto pOut = TNodeEDatNet<TFlt, TFlt>::New();
	for(auto it=vGraphs.begin(); it!=vGraphs.end(); ++it)
	{
		for (auto NI = (*it)->BegNI(); NI < (*it)->EndNI(); NI++)
		{
			int iParentID = NI.GetId();
			if(!pOut->IsNode(iParentID))
				pOut->AddNode(iParentID);
			for (int e = 0; e < NI.GetOutDeg(); ++e)
			{
				int iChildID = NI.GetOutNId(e);
				if(!pOut->IsNode(iChildID))
					pOut->AddNode(iChildID);
				pOut->AddEdge(iParentID,iChildID);
			}
		}
	}

	return pOut;
}

// Output: pGraph with the updated sum weight in each node
// Requires initial values for the nodes of pDAGGraph (edges are not needed)
void UpdateSumWeightOfShortestPath(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, TPt<TNodeEDatNet<TFlt, TFlt>>& pDAGGraph, const std::vector<int> &vSeedIDs, double dThreshold)
{
	for(auto it=vSeedIDs.begin(); it!=vSeedIDs.end(); ++it)
		Dijkstra(pGraph, *it, dThreshold, pDAGGraph);
}

TPt<TNodeEDatNet<TFlt, TFlt>> DAG2(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs, double dThreshold)
{
	// Vector of MIOA graphs per seed node
	std::vector<TPt<TNodeEDatNet<TFlt, TFlt>>> vMIOAGraphs;
	for(auto it=vSeedIDs.begin(); it!=vSeedIDs.end(); ++it)
		vMIOAGraphs.push_back(MIOA(pGraph, *it, dThreshold));
	auto pOut = GraphUnion(vMIOAGraphs);
	for (auto NI = pOut->BegNI(); NI < pOut->EndNI(); NI++)
		pOut->SetNDat(NI.GetId(), FLT_MAX);

	UpdateSumWeightOfShortestPath(pGraph, pOut, vSeedIDs, dThreshold);

	// Traverse the edges and prune the graph
	for (auto EI = pOut->BegEI(); EI < pOut->EndEI(); EI++)
	{
		if(EI.GetDstNDat().Val < EI.GetSrcNDat().Val)
			pOut->DelEdge(EI.GetSrcNId(), EI.GetDstNId());
	}

	return pOut;
}

TPt<TNodeEDatNet<TFlt, TFlt>> DAG2(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNode, double dThreshold)
{
	std::vector<int> vSeedIDs;
	vSeedIDs.push_back(sourceNode);
	return DAG2(pGraph, vSeedIDs, dThreshold);
}