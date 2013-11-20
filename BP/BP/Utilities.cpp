#include <iostream>
#include <map>
#include <algorithm>
#include <queue>
#include <float.h>
#include "Utilities.h"
#include <Snap.h>


using namespace std;

#ifdef WIN32

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

TPt<TNodeEDatNet<TFlt, TFlt>> CopyGraph(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph)
{
	// Copy the nodes of pGraph
	auto pCopy = TNodeEDatNet<TFlt, TFlt>::New();
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
	{
		int NodeID = NI.GetId();
		pCopy->AddNode(NodeID);
		pCopy->SetNDat(NodeID, NI.GetDat().Val);
	}
	for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
	{
		pCopy->AddEdge(EI.GetSrcNId(), EI.GetDstNId());
		pCopy->SetEDat(EI.GetSrcNId(), EI.GetDstNId(), EI.GetDat().Val);
	}
	return pCopy;
}

void AddSuperRootNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedNodes, int superRootNodeID)
{
	pGraph->AddNode(superRootNodeID);
	pGraph->SetNDat(superRootNodeID, 1.0);
	for(int srcNode: vSeedNodes)
	{
		pGraph->AddEdge(superRootNodeID, srcNode);
		pGraph->SetEDat(superRootNodeID, srcNode, 1.0);
	}
}

void RandomGraphInitialization(TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph)
{
	srand(time(NULL));
	for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
		pGraph->SetEDat(EI.GetSrcNId(), EI.GetDstNId(), (double) rand() / RAND_MAX);
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
		pGraph->SetNDat(NI.GetId(), (double) rand() / RAND_MAX);
}

void ResetGraphBelief(TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph)
{
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
		pGraph->SetNDat(NI.GetId(), 0.0);
}

//! Given pGraph with data about edge weights, computes the distance of the shortest paths from sourceNode
//! and returns the result in the nodes of pDAGGraph.
//! Updates the edges if bUpdateEdges is set to true. Default is false. In that case only the node data is updated with the shortest distance to sourceNode.
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

	struct Order	{inline bool operator()(std::pair<int,double> const& a, std::pair<int,double> const& b) const	{return a.second > b.second;}};
	std::priority_queue<std::pair<int,double>, std::vector<std::pair<int,double>>, Order> nodesToVisit;

	// Distance from source node to itself is 0
	pDAGGraph->SetNDat(sourceNode, 0);
	nodesToVisit.push(std::make_pair(sourceNode,0));

	// Beginning of the loop of Dijkstra algorithm

	while(!nodesToVisit.empty())
	{
		// Find the vertex in queue with the smallest distance and remove it
		int iParentID = -1;
		while (!nodesToVisit.empty() && visitedNodes[iParentID = nodesToVisit.top().first])
			nodesToVisit.pop();
		if (iParentID == -1) break;

		// mark the vertex with the shortest distance
		visitedNodes[iParentID]=true;

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
					//1. update distance
					//2. update the predecessor
					//3. push new shortest rank of chidren nodes
					pDAGGraph->SetNDat(iChildID, alt);
					mapPrevious[iChildID]= iParentID;
					nodesToVisit.push(std::make_pair(iChildID,alt));
				}
			}
		}

	}

	if(bUpdateEdges)
		for(auto it=mapPrevious.begin(); it!= mapPrevious.end(); ++it)
		{
			pDAGGraph->AddEdge(it->second, it->first);
			pDAGGraph->SetEDat(it->second,it->first, pGraph->GetEDat(it->second,it->first));
		}
}

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


///////////
// DAG 1 //
///////////

TPt<TNodeEDatNet<TFlt, TFlt>> GenerateDAG1(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, const std::vector<int>& seedNodes, double threshold)
{
	// Copy pGraph into pGraph_DAG1
	auto pGraph_DAG1 = TNodeEDatNet<TFlt, TFlt>::New();

	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
		pGraph_DAG1->AddNode(NI.GetId());
	
	for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
	{
		pGraph_DAG1->AddEdge(EI.GetSrcNId(),EI.GetDstNId());
		pGraph_DAG1->SetEDat(EI.GetSrcNId(),EI.GetDstNId(), pGraph->GetEDat(EI.GetSrcNId(),EI.GetDstNId()));
	}

	// Create a super root in order to update in one pass all the shortest paths from vSeedIDs nodes
	AddSuperRootNode(pGraph_DAG1, seedNodes);
	pGraph_DAG1 = MIOA(pGraph_DAG1, INT_MAX, threshold);
	// Remove the artificial super root node
	pGraph_DAG1->DelNode(INT_MAX);


	// Add back other edges with the condition r(u)<r(v)
	for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
	{
		int u = EI.GetSrcNId(), v = EI.GetDstNId();
		if(pGraph_DAG1->GetNDat(u)< pGraph_DAG1->GetNDat(v))
		{
			if (!pGraph_DAG1->IsEdge(u,v))
			{
				pGraph_DAG1->AddEdge(u,v);
				pGraph_DAG1->SetEDat(u,v,EI.GetDat());
			}
		}
	}

	return pGraph_DAG1;
}

TPt<TNodeEDatNet<TFlt, TFlt>> GenerateDAG1(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, int sourceNode, double threshold)
{
	std::vector<int> vSeedIDs;
	vSeedIDs.push_back(sourceNode);
	return GenerateDAG1(pGraph, vSeedIDs, threshold);
}

///////////
// DAG 2 //
///////////

//! Graph union between the graphs of vGraphs
//! @note Does not copy edge or node data.
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

TPt<TNodeEDatNet<TFlt, TFlt>> GenerateDAG2(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs, double dThreshold)
{
	// Vector of MIOA graphs per seed node
	std::vector<TPt<TNodeEDatNet<TFlt, TFlt>>> vMIOAGraphs;

	// Compute the union of MIOA for each node of vSeedIDs
	for(auto it=vSeedIDs.begin(); it!=vSeedIDs.end(); ++it)
		vMIOAGraphs.push_back(MIOA(pGraph, *it, dThreshold));
	auto pOut = GraphUnion(vMIOAGraphs);

	// Set node data
	for (auto NI = pOut->BegNI(); NI < pOut->EndNI(); NI++)
		pOut->SetNDat(NI.GetId(), FLT_MAX);

	// Copy the edge weights from pGraph
	for (auto EI = pOut->BegEI(); EI < pOut->EndEI(); EI++)
		pOut->SetEDat(EI.GetSrcNId(), EI.GetDstNId(), pGraph->GetEDat(EI.GetSrcNId(), EI.GetDstNId()));

	// Create a super root in order to update in one pass all the shortest paths from vSeedIDs nodes
	AddSuperRootNode(pOut, vSeedIDs);
	Dijkstra(pOut, INT_MAX, dThreshold, pOut);
	// Remove the artificial super root node
	pOut->DelNode(INT_MAX);

	// Traverse the edges and prune the graph
	for (auto EI = pOut->BegEI(); EI < pOut->EndEI(); EI++)
	{
		if(EI.GetDstNDat().Val < EI.GetSrcNDat().Val)
			pOut->DelEdge(EI.GetSrcNId(), EI.GetDstNId());
	}

	return pOut;
}

TPt<TNodeEDatNet<TFlt, TFlt>> GenerateDAG2(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNode, double dThreshold)
{
	std::vector<int> vSeedIDs;
	vSeedIDs.push_back(sourceNode);
	return GenerateDAG2(pGraph, vSeedIDs, dThreshold);
}

// End of DAG2
 
void CalculateRankFromSource_BellmanFord(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, int sourceNode, std::vector<int> &vResult)
{
	int numNodes = pGraph->GetNodes();
	vResult.reserve(numNodes);
	for(int i = 0; i<numNodes; ++i)
		vResult.push_back(INT_MAX);
	vResult[sourceNode]=0;

	for (size_t i = 0; i < (vResult.size() - 1); ++i)
	{
		for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
		{
			double alt = vResult[EI.GetSrcNId()] - 1;
			if (alt < vResult[EI.GetDstNId()])
				vResult[EI.GetDstNId()] = alt;
		}
	}

	for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
	{
		if ((vResult[EI.GetSrcNId()] - 1) < vResult[EI.GetDstNId()])
		{
			OutputDebugStringA("Bellman Ford error : Negative cycle detected");
			exit(-1);
		}
	}

	// Revert the rank to positive values
	for(int i = 0; i<numNodes; ++i)
		vResult[i] *= -1;
}

void CalculateRankFromSource_BellmanFord(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, const std::vector<int> vSeedNodes, std::vector<int> &vResult)
{
	auto pTemp = CopyGraph(pGraph);
	int superRootNodeID = pGraph->GetMxNId();
	AddSuperRootNode(pTemp, vSeedNodes, superRootNodeID);
	CalculateRankFromSource_BellmanFord(pTemp, superRootNodeID, vResult);
}

void CalculateRankFromSource(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, int sourceNode, std::vector<int> &vResult)
{
	int currentRank=0;
	int numNodes = pGraph->GetNodes();
	vResult.reserve(pGraph->GetNodes());
	for(int i = 0; i<numNodes; ++i)
		vResult.push_back(INT_MAX);
	vResult[sourceNode]=0;
	++currentRank;

	// Used to store the nodes which have been traversed during the breadth-first search traversal
	std::map<int, bool> visitedNodes;
	std::queue<int> queue;
	queue.push(sourceNode);
	visitedNodes[sourceNode] = true;

	while(!queue.empty())
	{
		int numNodesToProcess = queue.size();
		for(int i=0; i< numNodesToProcess; ++i)
		{
			int nodeID = queue.front();
			auto parent = pGraph->GetNI(nodeID);        
			int numChildren = parent.GetOutDeg();
			for(int i = 0; i < numChildren; ++i)
			{
				int iChildID = parent.GetOutNId(i);
				// If we update a node which has already been updated previously, we need to add it back to the queue to update the other nodes.
				if(vResult[iChildID]<currentRank)
					visitedNodes.erase(iChildID);
				vResult[iChildID] = currentRank;
				// Mark the child
				if(visitedNodes.insert(std::pair<int,bool>(iChildID,true)).second)
					queue.push(iChildID);
			}

			queue.pop();
		}
		++currentRank;
	}
}

void CalculateRankFromSource(const TPt<TNodeEDatNet<TFlt, TFlt>> &pGraph, const std::vector<int> vSeedNodes, std::vector<int> &vResult)
{
	auto pTemp = CopyGraph(pGraph);
	int superRootNodeID = pGraph->GetMxNId();
	AddSuperRootNode(pTemp, vSeedNodes, superRootNodeID);
	CalculateRankFromSource(pTemp, superRootNodeID, vResult);
}

double BPError(const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph1, const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph2, const std::function<double(double, double)> &fn)
{
	double error=0.0;
	unsigned int count=0;
	for(auto NI=pGraph1->BegNI(); NI < pGraph1->EndNI(); NI++)
	{
		
		if(NI.GetDat().Val != 0.0 || pGraph2->GetNDat(NI.GetId()).Val != 0.0)
		{
			++count;
			error+=fn(NI.GetDat().Val, pGraph2->GetNDat(NI.GetId()).Val);
		}
	}
	return error/count;
}