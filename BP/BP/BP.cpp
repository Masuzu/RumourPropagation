#ifdef _USE_libDAI
#include <dai/alldai.h> // Include main libDAI header file
#include <dai/jtree.h>
#include <dai/bp.h>
#include <dai/decmap.h>

#pragma comment(lib, "libdai.lib")
#endif

#ifdef _USE_CUDA
#pragma comment(lib, "cudart_static.lib")
#endif

#define __TBB_NO_IMPLICIT_LINKAGE	1
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_queue.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/spin_mutex.h>

#include <queue>
#include <map>
#include <iostream>

#include <Snap.h>

#include "AEWorkerThread.h"
#include "Utilities.h"
#include "BP.h"

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

#ifdef _USE_libDAI
void ExactBP_Marginalization(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID)
{
	///////////////////////////////////////////
	// Convert the graph into a factor graph //

	// Stores the variables corresponding to each node of pGraph
	std::map<int, dai::Var> VarMap;
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
	{
		// Definition of a binary variable for the node NI
		VarMap[NI.GetId()] = dai::Var(NI.GetId(), 2);
	}

	// Traverse all the nodes of the graph and build the factors
	std::vector<dai::Factor> GraphFactors;
	for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
	{
		int numParents = NI.GetInDeg();

		// Initialize the variable currentNodeID|parentsOf_currentNodeID
		dai::VarSet VarParents;	// empty
		for(int i = 0; i < numParents; ++i)
			VarParents |= VarMap[NI.GetInNId(i)];

		int currentNodeID=NI.GetId();
		dai::Factor VGivenParents(VarParents | VarMap[currentNodeID]);

		// ^__^

		if(numParents+1 > 32)
		{
			OutputDebugStringA("Runtime error : Max int limit reached");
			exit(-1);
		}

		/////////////////////////////////////////////////////////////////
		// Compute the conditional probability table for the current node
		if(numParents==0)
		{
			// The source node has a probability of activation of 1
			if(currentNodeID==sourceNodeID)
			{
				VGivenParents.set(0, 0);
				VGivenParents.set(1, 1);
			}
			else
			{
				// Default state for nodes: non activated
				VGivenParents.set(0, 1);
				VGivenParents.set(1, 0);
			}
		}
		else
			// Define the conditional probabilities for 'NI=0'
		{
			// Case in which the parents are not passing any message
			if(currentNodeID==sourceNodeID)				
				VGivenParents.set(0, 0.0);
			else
				VGivenParents.set(0, 1.0);

			// The other 2^numParents-1 cases for which 'NI=0'
			for (size_t i = 1; i<pow(2.0, numParents); ++i)
			{
				// The source node has a probability of activation of 1
				if(currentNodeID==sourceNodeID)				
					VGivenParents.set(i, 0.0);
				else
				{
					double p = 1.0;
					for(int j = 0; j < numParents; ++j)
					{
						if(i & (1<<j))
							p *= pGraph->GetEDat(NI.GetInNId(j), currentNodeID);
					}
					VGivenParents.set(i, 1.0-p);
				}
			}

			// Define the conditional probabilities for 'NI=1' (complementary event)
			int idxStart=pow(2.0, numParents);
			for (size_t i = idxStart; i<pow(2.0, numParents+1); ++i)
				VGivenParents.set(i, 1.0-VGivenParents.get(i-idxStart));
		}

		GraphFactors.push_back(VGivenParents);
	}

	// End of the factor computation
	// ^__^

	dai::FactorGraph Network( GraphFactors );

	// Write factorgraph to a file
	// Network.WriteToFile( "sprinkler.fg" );

	dai::Factor P;
	for( size_t I = 0; I < Network.nrFactors(); ++I )
		P *= Network.factor(I);

	// P.normalize(); // Not necessary: a Bayesian network is already normalized by definition
	// Calculate some probabilities (example for reference)
	/*
	dai::Real denom = P.marginal( VarMap[sourceNodeID] )[1];
	std::cout << "P(source=1) = " << denom << std::endl;
	std::cout << "P(node_4=1 | source=1) = " << P.marginal( dai::VarSet( VarMap[4], VarMap[sourceNodeID] ) )[3] / denom << std::endl;
	std::cout << "P(R=1 | W=1) = " << P.marginal( VarSet( R, W ) )[3] / denom << endl;
	*/

	// Update the belief of the nodes of pGraph
	pGraph->SetNDat(sourceNodeID, 1.0);

	std::map<int, bool> visitedNodes;
	std::queue<int> queue;
	queue.push(sourceNodeID);
	visitedNodes[sourceNodeID] = true;
	auto VarSourceNode = VarMap[sourceNodeID];
	while(!queue.empty())
	{
		int nodeID = queue.front();

		auto parent = pGraph->GetNI(nodeID);	
		int numChildren = parent.GetOutDeg();
		// Update the belief of the children of parent
		for(int i = 0; i < numChildren; ++i)
		{
			int iChildID = parent.GetOutNId(i);
			// Compute P(node_iChildID=1 | source=1)
			pGraph->SetNDat(iChildID, P.marginal( dai::VarSet( VarMap[iChildID], VarSourceNode ) )[3]);

			// Mark the child
			if(visitedNodes.insert(std::pair<int,bool>(iChildID,true)).second)
				queue.push(iChildID);
		}

		queue.pop();
	}
}

void ExactBP_Marginalization(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs)
{
	AddSuperRootNode(pGraph, vSeedIDs);

	ExactBP_Marginalization(pGraph, INT_MAX);

	// Remove the artificial super root node
	pGraph->DelNode(INT_MAX);
}
#endif

/////////////////
// One-pass BP //
/////////////////

void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID)
{
#ifdef _DEBUG_INFO
	int meter = 0;
#endif

	pGraph->SetNDat(sourceNodeID, 1.0);

	// Used to store the nodes which have been traversed during the breadth-first search traversal
	std::map<int, bool> visitedNodes;
	std::queue<int> queue;
	queue.push(sourceNodeID);
	visitedNodes[sourceNodeID] = true;

	while(!queue.empty())
	{
		int nodeID = queue.front();
#ifdef _DEBUG_INFO
		++meter;
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

			// Mark the child
			if(visitedNodes.insert(std::pair<int,bool>(iChildID,true)).second)
				queue.push(iChildID);
		}

		queue.pop();
	}

#ifdef _DEBUG_INFO
	std::cout << "Number of nodes visited: " <<  meter << std::endl;
#endif
}

void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs)
{
	AddSuperRootNode(pGraph, vSeedIDs);

	PropagateFromNode(pGraph, INT_MAX);

	// Remove the artificial super root node
	pGraph->DelNode(INT_MAX);
}

// DAG version

void PropagateFromNodeDAG(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID)
{
	pGraph->SetNDat(sourceNodeID, 1.0);

	std::queue<int> queue;
	queue.push(sourceNodeID);

	while(!queue.empty())
	{
		int nodeID = queue.front();

		auto parent = pGraph->GetNI(nodeID);	
		int numChildren = parent.GetOutDeg();
		// Update the belief of the children of parent
		for(int i = 0; i < numChildren; ++i)
		{
			int iChildID = parent.GetOutNId(i);
			// Do p(u) = 1-(1-p(u))*(1-p(u,v)*p(v))
			pGraph->SetNDat(iChildID,
				1.0 - (1-parent.GetOutNDat(i).Val) * (1-parent.GetOutEDat(i).Val*parent.GetDat().Val));

			queue.push(iChildID);
		}

		queue.pop();
	}
}

void PropagateFromNodeDAG(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs)
{
	AddSuperRootNode(pGraph, vSeedIDs);

	PropagateFromNodeDAG(pGraph, INT_MAX);

	// Remove the artificial super root node
	pGraph->DelNode(INT_MAX);
}


/////////////////
// Parallel BP //
/////////////////

static void ParallelBPFromNode_UpdateChild(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, TNodeEDatNet<TFlt, TFlt>::TNodeI &parent,
	tbb::concurrent_hash_map<int, bool> &visitedNodes, tbb::concurrent_queue<int> &queue, int i)
{
	int iChildID = parent.GetOutNId(i);
	// Do p(u) = 1-(1-p(u))*(1-p(u,v)*p(v))
	pGraph->SetNDat(iChildID,
		1.0 - (1-parent.GetOutNDat(i).Val) * (1-parent.GetOutEDat(i).Val*parent.GetDat().Val));

	// http://www.threadingbuildingblocks.org/docs/help/reference/containers_overview/concurrent_hash_map_cls/concurrent_operations.htm
	// To guarantee that only one instance of a resource exists simultaneously for a given key,
	// to construct the resource: Obtain an accessor to the key in the map before constructing the resource.
	tbb::concurrent_hash_map<int, bool>::accessor acc;
	// Mark the child
	if(visitedNodes.insert(acc,iChildID))
	{
		// Current thread inserted key and has exclusive access.
		queue.push(iChildID);
	}
}

void ParallelBPFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID)
{
#ifdef _DEBUG_INFO
	int meter = 0;
#endif
	pGraph->SetNDat(sourceNodeID, 1.0);

	// Like std::list, insertion of new items does not invalidate any iterators, nor change the order of items already in the map. Insertion and traversal may be concurrent.
	tbb::concurrent_hash_map<int, bool> visitedNodes;

	tbb::concurrent_queue<int> queue;
	queue.push(sourceNodeID);

	int nodeID;
	while(queue.try_pop(nodeID))
	{
		auto parent = pGraph->GetNI(nodeID);
		int numChildren = parent.GetOutDeg();
#ifdef _DEBUG_INFO
		meter+=numChildren;
#endif
		if(numChildren > 50)
		{
#ifdef _DEBUG_INFO
			std::cout << "Parallel update of " <<  numChildren << " nodes" << std::endl;
#endif
			tbb::parallel_for(tbb::blocked_range<int>(0,numChildren, 50), 
				[&](const tbb::blocked_range<int>& r)
			{
				for (int i=r.begin();i!=r.end();++i)
					ParallelBPFromNode_UpdateChild(pGraph, parent, visitedNodes, queue, i);
			}
			);
		}
		else
			// Not enough elements to do the work concurrently
		{
			for(int i = 0; i < numChildren; ++i)
				ParallelBPFromNode_UpdateChild(pGraph, parent, visitedNodes, queue, i);
		}
	}
#ifdef _DEBUG_INFO
	std::cout << "Number of nodes visited: " <<  meter << std::endl;
#endif
}

void ParallelBPFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs)
{
	AddSuperRootNode(pGraph, vSeedIDs);

	ParallelBPFromNode(pGraph, INT_MAX);

	// Remove the artificial super root node
	pGraph->DelNode(INT_MAX);
}

// DAG version

static void ParallelBPFromNodeDAG_UpdateChild(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, TNodeEDatNet<TFlt, TFlt>::TNodeI &parent,
	tbb::concurrent_queue<int> &queue, int i)
{
	int iChildID = parent.GetOutNId(i);
	// Do p(u) = 1-(1-p(u))*(1-p(u,v)*p(v))
	pGraph->SetNDat(iChildID,
		1.0 - (1-parent.GetOutNDat(i).Val) * (1-parent.GetOutEDat(i).Val*parent.GetDat().Val));

	queue.push(iChildID);
}

void ParallelBPFromNodeDAG(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID)
{
#ifdef _DEBUG_INFO
	int meter = 0;
#endif
	pGraph->SetNDat(sourceNodeID, 1.0);

	tbb::concurrent_queue<int> queue;
	queue.push(sourceNodeID);

	int nodeID;
	while(queue.try_pop(nodeID))
	{
		auto parent = pGraph->GetNI(nodeID);
		int numChildren = parent.GetOutDeg();
#ifdef _DEBUG_INFO
		meter+=numChildren;
#endif
		if(numChildren > 50)
		{
			tbb::parallel_for(tbb::blocked_range<int>(0,numChildren, 50), 
				[&](const tbb::blocked_range<int>& r)
			{
				for (int i=r.begin();i!=r.end();++i)
					ParallelBPFromNodeDAG_UpdateChild(pGraph, parent, queue, i);
			}
			);
		}
		else
			// Not enough elements to do the work concurrently
		{
			for(int i = 0; i < numChildren; ++i)
				ParallelBPFromNodeDAG_UpdateChild(pGraph, parent, queue, i);
		}
	}

#ifdef _DEBUG_INFO
	std::cout << "Number of nodes visited: " << meter << std::endl;
#endif
}

void ParallelBPFromNodeDAG(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs)
{
	AddSuperRootNode(pGraph, vSeedIDs);

	ParallelBPFromNodeDAG(pGraph, INT_MAX);

	// Remove the artificial super root node
	pGraph->DelNode(INT_MAX);
}

/////////////////////////////////
// Parallel BP 1D partitioning //
/////////////////////////////////

static void ParallelBPFromNode_LevelSynchronous_UpdateNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph,
	tbb::concurrent_hash_map<int, bool> &visitedNodes, tbb::concurrent_queue<int> &queue, bool MT = true)
{
	static tbb::spin_mutex sMutex;

	int iParentID;
	queue.try_pop(iParentID);
	auto parent = pGraph->GetNI(iParentID);
	int numChildren = parent.GetOutDeg();
	for(int j=0; j<numChildren; ++j)
	{
		int iChildID = parent.GetOutNId(j);
		// Do p(u) = 1-(1-p(u))*(1-p(u,v)*p(v))
		if(MT)
		{
			tbb::spin_mutex::scoped_lock lock(sMutex);
			pGraph->SetNDat(iChildID,
				1.0 - (1-parent.GetOutNDat(j).Val) * (1-parent.GetOutEDat(j).Val*parent.GetDat().Val));
		}
		else
			pGraph->SetNDat(iChildID,
			1.0 - (1-parent.GetOutNDat(j).Val) * (1-parent.GetOutEDat(j).Val*parent.GetDat().Val));

		// http://www.threadingbuildingblocks.org/docs/help/reference/containers_overview/concurrent_hash_map_cls/concurrent_operations.htm
		// To guarantee that only one instance of a resource exists simultaneously for a given key,
		// to construct the resource: Obtain an accessor to the key in the map before constructing the resource.
		tbb::concurrent_hash_map<int, bool>::accessor acc;
		// Mark the child
		if(visitedNodes.insert(acc,iChildID))
		{
			// Current thread inserted key and has exclusive access.
			queue.push(iChildID);
		}
	}
}

void ParallelBPFromNode_LevelSynchronous(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID)
{
#ifdef _DEBUG_INFO
	int meter = 0;
	int currentLevel = 0;
#endif
	pGraph->SetNDat(sourceNodeID, 1.0);

	// Like std::list, insertion of new items does not invalidate any iterators, nor change the order of items already in the map. Insertion and traversal may be concurrent.
	tbb::concurrent_hash_map<int, bool> visitedNodes;

	tbb::concurrent_queue<int> queue;
	queue.push(sourceNodeID);

	while(!queue.empty())
	{
		int numNodesToProcess = queue.unsafe_size();
#ifdef _DEBUG_INFO
		meter+=numNodesToProcess;
		++currentLevel;
#endif
		// Parse level l
		if(numNodesToProcess > 50)
		{
			tbb::parallel_for(tbb::blocked_range<int>(0, numNodesToProcess, 50), 
				[&](const tbb::blocked_range<int>& r)
			{
				for (int i=r.begin();i!=r.end();++i)
					ParallelBPFromNode_LevelSynchronous_UpdateNode(pGraph, visitedNodes, queue, true); 
			}
			);
		}
		else
			// Not enough elements to do the work concurrently
		{
			for (int i=0;i<numNodesToProcess;++i)
				ParallelBPFromNode_LevelSynchronous_UpdateNode(pGraph, visitedNodes, queue, false); 
		}
#ifdef _DEBUG_INFO
		std::cout << "Number of nodes at level " << currentLevel << " : " << numNodesToProcess << std::endl;
#endif
	}

#ifdef _DEBUG_INFO
	std::cout << "Number of nodes visited: " << meter << std::endl;
#endif
}

void ParallelBPFromNode_LevelSynchronous(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs)
{
	AddSuperRootNode(pGraph, vSeedIDs);

	ParallelBPFromNode_LevelSynchronous(pGraph, INT_MAX);

	// Remove the artificial super root node
	pGraph->DelNode(INT_MAX);
}

// DAG version

static void ParallelBPFromNodeDAG_LevelSynchronous_UpdateNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph,
	tbb::concurrent_queue<int> &queue, bool MT = true)
{
	static tbb::spin_mutex sMutex;

	int iParentID;
	queue.try_pop(iParentID);
	auto parent = pGraph->GetNI(iParentID);
	int numChildren = parent.GetOutDeg();
	for(int j=0; j<numChildren; ++j)
	{
		int iChildID = parent.GetOutNId(j);
		// Do p(u) = 1-(1-p(u))*(1-p(u,v)*p(v))
		if(MT)
		{
			tbb::spin_mutex::scoped_lock lock(sMutex);
			pGraph->SetNDat(iChildID,
				1.0 - (1-parent.GetOutNDat(j).Val) * (1-parent.GetOutEDat(j).Val*parent.GetDat().Val));
		}
		else
			pGraph->SetNDat(iChildID,
			1.0 - (1-parent.GetOutNDat(j).Val) * (1-parent.GetOutEDat(j).Val*parent.GetDat().Val));

		queue.push(iChildID);
	}
}

void ParallelBPFromNodeDAG_LevelSynchronous(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID)
{
#ifdef _DEBUG_INFO
	int meter = 0;
#endif
	pGraph->SetNDat(sourceNodeID, 1.0);

	tbb::concurrent_queue<int> queue;
	queue.push(sourceNodeID);

	while(!queue.empty())
	{
		int numNodesToProcess = queue.unsafe_size();
#ifdef _DEBUG_INFO
		meter+=numNodesToProcess;
#endif
		// Parse level l
		if(numNodesToProcess > 50)
		{
			tbb::parallel_for(tbb::blocked_range<int>(0, numNodesToProcess, 50), 
				[&](const tbb::blocked_range<int>& r)
			{
				for (int i=r.begin();i!=r.end();++i)
					ParallelBPFromNodeDAG_LevelSynchronous_UpdateNode(pGraph, queue, true); 
			}
			);
		}
		else
			// Not enough elements to do the work concurrently
		{
			for (int i=0;i<numNodesToProcess;++i)
				ParallelBPFromNodeDAG_LevelSynchronous_UpdateNode(pGraph, queue, false); 
		}
	}

#ifdef _DEBUG_INFO
	std::cout << meter << std::endl;
#endif
}

void ParallelBPFromNodeDAG_LevelSynchronous(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs)
{
	AddSuperRootNode(pGraph, vSeedIDs);

	ParallelBPFromNodeDAG_LevelSynchronous(pGraph, INT_MAX);

	// Remove the artificial super root node
	pGraph->DelNode(INT_MAX);
}


/////////////////////////////////////////////////
// ParallelBPFromNode_SingleNodeUpdate_UpdateNode

static void ParallelBPFromNode_SingleNodeUpdate_UpdateNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int>& pRankMap,
	tbb::concurrent_hash_map<int, bool> &visitedNodes, tbb::concurrent_queue<int> &queue,
	unsigned int refRank, unsigned int currentRank)
{
	unsigned int nextRank = currentRank+1;
	int nodeID;
	queue.try_pop(nodeID);
	auto node = pGraph->GetNI(nodeID);
	if((pRankMap[nodeID]-refRank) == currentRank)
	{
		int numParents = node.GetInDeg();
		double p=1.0;
		for(int j=0; j<numParents; ++j)
			p *= (1.0 - node.GetInNDat(j)*node.GetInEDat(j));
		node.GetDat().Val = 1-p;
		int numChildren = node.GetOutDeg();
		for(int j=0; j<numChildren; ++j)
		{
			int candidateNodeID = node.GetOutNId(j);

			if((pRankMap[candidateNodeID]-refRank)==nextRank)
			{
				tbb::concurrent_hash_map<int, bool>::accessor acc;
				// Mark the child
				if(visitedNodes.insert(acc,candidateNodeID))
				{
					// Current thread inserted key and has exclusive access.
					queue.push(candidateNodeID);
				}
			}
		}
	}
}

void ParallelBPFromNode_SingleNodeUpdate(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int>& pRankMap, int sourceNodeID)
{
#ifdef _DEBUG_INFO
	int meter = 0;
#endif
	pGraph->SetNDat(sourceNodeID, 1.0);

	// Like std::list, insertion of new items does not invalidate any iterators, nor change the order of items already in the map. Insertion and traversal may be concurrent.
	tbb::concurrent_hash_map<int, bool> visitedNodes;

	tbb::concurrent_queue<int> queue;
	unsigned int currentRank = 1;
	unsigned int refRank = pRankMap[sourceNodeID];

	// Initialize queue
	auto parent = pGraph->GetNI(sourceNodeID);
	int numChildren = parent.GetOutDeg();
	for(int j=0; j<numChildren; ++j)
	{
		int candidateNodeID = parent.GetOutNId(j);
		if((pRankMap[candidateNodeID]-refRank)==1)	
			queue.push(candidateNodeID);
	}
	while(!queue.empty())
	{
		int numNodesToProcess = queue.unsafe_size();
#ifdef _DEBUG_INFO
		meter += numNodesToProcess;
#endif
		if(numNodesToProcess > 50)
		{
			tbb::parallel_for(tbb::blocked_range<int>(0, numNodesToProcess, 50), 
				[&](const tbb::blocked_range<int>& r)
			{
				for (int i=r.begin();i!=r.end();++i)
					ParallelBPFromNode_SingleNodeUpdate_UpdateNode(pGraph, pRankMap, visitedNodes, queue, refRank, currentRank); 
			}
			);
		}
		else
			// Not enough elements to do the work concurrently
		{
			for (int i=0;i<numNodesToProcess;++i)
				ParallelBPFromNode_SingleNodeUpdate_UpdateNode(pGraph, pRankMap, visitedNodes, queue, refRank, currentRank); 
		}
		++currentRank;
	}

#ifdef _DEBUG_INFO
	std::cout << "Number of nodes visited: " << meter << std::endl;
#endif
}

void ParallelBPFromNode_SingleNodeUpdate(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int>& pRankMap, const std::vector<int> &vSeedIDs)
{
	int superRootNodeID = pGraph->GetMxNId();
	AddSuperRootNode(pGraph, vSeedIDs, superRootNodeID);

	ParallelBPFromNode_SingleNodeUpdate(pGraph, pRankMap, superRootNodeID);

	// Remove the artificial super root node
	pGraph->DelNode(superRootNodeID);
}

double InfluenceSpreadFromSeedNodes (const TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph)
{
	double influenceSpread = 0.0;
	for (auto NI = pGraph->BegNI() ; NI < pGraph->EndNI(); NI++)
		influenceSpread += NI.GetDat().Val;

	return influenceSpread;
}

static void GetPeerSeeds(std::map<int,TPt<TNodeEDatNet<TFlt, TFlt>> > mMIOAs, int nodeID, std::vector<int> &vPeerSeeds)
{
	for(auto it = mMIOAs.begin(); it!= mMIOAs.end(); it++)
	{
		if(it->first!=nodeID)
		{
			bool isOverlapped = false;
			auto pMIOA = it->second;
			for(auto v = mMIOAs[nodeID]->BegNI(); v < mMIOAs[nodeID]->EndNI();v++)
			{
				if (isOverlapped) break;
				for(auto u = pMIOA->BegNI();u < pMIOA->EndNI(); u++)
				{
					if(v.GetId()== u.GetId())
					{
						isOverlapped = true;
						break;
					}
				}
			}
			if(isOverlapped) vPeerSeeds.push_back(it->first);
		}
		else continue;
	}
}

void GreedyCELF(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int setSize,std::vector<int> &vSeedSet)
{

	std::map<int,double> mSpreadIncrement;
	auto pGraph_main = CopyGraph(pGraph);
	auto pGraph_temp = TNodeEDatNet<TFlt, TFlt>::New();
	double influence = 0.0;

	//Failure of using PeerSeeds due to insufficient memory
	//std::map<int,std::vector<int> > mPeerSeeds;
	//std::map<int,TPt<TNodeEDatNet<TFlt, TFlt>> > mMIOAs;

	/* Initialization*/
	for(auto NI = pGraph_main->BegNI();NI<pGraph_main->EndNI();NI++)
	{
		if(pGraph_main->IsNode(NI.GetId()))
		{
			int nodeID = NI.GetId();
			pGraph_temp = MIOA(pGraph_main, nodeID, 0);
			ResetGraphBelief(pGraph_temp);
			ParallelBPFromNode_LevelSynchronous(pGraph_temp,nodeID);
			mSpreadIncrement[nodeID]=InfluenceSpreadFromSeedNodes(pGraph_temp);
			//mMIOAs.insert(std::make_pair(i,pGraph_v));
		}
	}

	/*
	//build PeerSeeds
	//Failure due to the insufficient memory
	for (int v =0; v<pGraph->GetNodes();++v)
	if(pGraph->IsNode(v))
	mPeerSeeds[v]=GetPeerSeeds(mMIOAs,v);
	*/

	//cout<<"Precomputation finished"<<endl;

	for (int k=0;k<setSize;++k)
	{
		/* select the i'th seed by finding u = argmax(mSpreadIncrement)*/
		auto it = std::max_element(mSpreadIncrement.begin(),mSpreadIncrement.end(),
			[&](std::pair<int,double> const& a, std::pair<int,double> const& b) {
				return a.second < b.second;
		}
		);
		int SeedID = it->first;
		cout << "Round "<<k<<": the most influential node is: "<<SeedID <<endl;

		/* calculate the current influence spread */

		vSeedSet.push_back(SeedID);
		pGraph_main =CopyGraph(pGraph);
		pGraph_main = GenerateDAG1(pGraph_main, vSeedSet, 0.0);
		ParallelBPFromNode_LevelSynchronous(pGraph_main, vSeedSet);
		influence = InfluenceSpreadFromSeedNodes(pGraph_main);
		cout << "Round "<<k<<": the influential spread is: "<<influence<<endl;

		/*remove the newly selected node*/
		mSpreadIncrement.erase(SeedID);
		// update incremental influence spread for each round
		double Delta_MAX = 0.0;
		for(auto NI = pGraph_main->BegNI();NI<pGraph_main->EndNI();NI++)
		{
			int nodeID = NI.GetId();
			// exclude the nodes in seed set
			auto result = std::find(vSeedSet.begin(),vSeedSet.end(), nodeID);
			if (result != vSeedSet.end()) continue;
			//CELF
			if(pGraph_main->IsNode(nodeID) && mSpreadIncrement[nodeID] > Delta_MAX)
			{
				vSeedSet.push_back(nodeID);
				pGraph_main =CopyGraph(pGraph);
				pGraph_temp = GenerateDAG1(pGraph_main, vSeedSet, 0);
				ParallelBPFromNode_LevelSynchronous(pGraph_temp, vSeedSet);
				mSpreadIncrement[nodeID]=InfluenceSpreadFromSeedNodes(pGraph_temp)-influence;

				if (mSpreadIncrement[nodeID]> Delta_MAX)
					Delta_MAX = mSpreadIncrement[nodeID];

				vSeedSet.pop_back();
			}
		}
	}

}

void ParallelGreedyCELF(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int setSize,std::vector<int> &vSeedSet)
{

	tbb::concurrent_unordered_map<int,double> mSpreadIncrement;
	auto pGraph_main = CopyGraph(pGraph);
	auto pGraph_temp = TNodeEDatNet<TFlt, TFlt>::New();
	double influence = 0.0; int i,chunk = 70;
	static tbb::spin_mutex sMutex;

	/* Initialization*/
	int numNodes = pGraph_main->GetMxNId();
#pragma omp parallel shared(pGraph_main,chunk,mSpreadIncrement) private(pGraph_temp,i)
	{
#pragma omp for schedule(dynamic,chunk) nowait
		for (i =0;i<numNodes;++i)
		{
			if(pGraph_main->IsNode(i))
			{
				pGraph_temp = MIOA(pGraph_main, i, 0);
				ResetGraphBelief(pGraph_temp);
				ParallelBPFromNode_LevelSynchronous(pGraph_temp, i);
				mSpreadIncrement[i]=InfluenceSpreadFromSeedNodes(pGraph_temp);

			}

		}
	}
	//cout<<"Precomputation finished"<<endl;

	for (int k=0;k<setSize;++k)
	{
		/* select the i'th seed by finding u = argmax(mSpreadIncrement)*/
		auto it = std::max_element(mSpreadIncrement.begin(),mSpreadIncrement.end(),
			[&](std::pair<int,double> const& a, std::pair<int,double> const& b) {
				return a.second < b.second;
		}
		);
		int SeedID = it->first;
		cout << "Round "<<k<<": the most influential node is: "<<SeedID <<endl;

		/* calculate the current influence spread */
		vSeedSet.push_back(SeedID);
		pGraph_main = CopyGraph(pGraph);
		pGraph_main = GenerateDAG1(pGraph_main, vSeedSet, 0.0);
		ParallelBPFromNode_LevelSynchronous(pGraph_main, vSeedSet);
		influence = InfluenceSpreadFromSeedNodes(pGraph_main);
		cout << "Round "<<k<<": the influential spread is: "<<influence<<endl;
		/*remove the newly selected node*/
		mSpreadIncrement.unsafe_erase(SeedID);

		// update incremental influence spread for each round
		double Delta_MAX = 0.0;
		for(auto NI = pGraph_main->BegNI();NI<pGraph_main->EndNI();NI++)
		{
			int nodeID = NI.GetId();
			// exclude the nodes in seed set
			auto result = std::find(vSeedSet.begin(),vSeedSet.end(), nodeID);
			if (result != vSeedSet.end()) continue;
			//CELF
			if(pGraph_main->IsNode(nodeID) && mSpreadIncrement[nodeID] > Delta_MAX)
			{
				vSeedSet.push_back(nodeID);
				pGraph_main =CopyGraph(pGraph);
				pGraph_temp = GenerateDAG1(pGraph_main, vSeedSet, 0);
				ParallelBPFromNode_LevelSynchronous(pGraph_temp, vSeedSet);
				mSpreadIncrement[nodeID]=InfluenceSpreadFromSeedNodes(pGraph_temp)-influence;

				if (mSpreadIncrement[nodeID]> Delta_MAX)
					Delta_MAX = mSpreadIncrement[nodeID];

				vSeedSet.pop_back();
			}
		}


	}

}

void NewGreedIC(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int setSize, int numRounds, std::vector<int> &vSeedSet)
{

	auto pGraph_temp = TNodeEDatNet<TFlt, TFlt>::New();
	tbb::concurrent_unordered_map<int,unsigned int> mInfluenceSpread;
	std::map<int, unsigned int> mNumOfNodesFromSource;
	std::vector<int> vResult;
	int i,chunk = 70;
	static tbb::spin_mutex sMutex;

	for (int k=0; k<setSize ; ++k)
	{
		mInfluenceSpread.clear();
		mNumOfNodesFromSource.clear();
		vResult.clear();

		for (int r =0;r<numRounds;++r)
		{
			pGraph_temp = TNodeEDatNet<TFlt, TFlt>::New();
			for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
				pGraph_temp->AddNode(NI.GetId());

			//1.Compute G' by removing each edge from G with probability 1 -p
			//2.Compute the set of vertices reachable from S
			//3.Compute as well the sets of vertices from v
			for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
			{
				double coinflip = (double) rand()/RAND_MAX;
				if( coinflip > 1.0-EI.GetDat().Val)
				{
					pGraph_temp->AddEdge(EI.GetSrcNId(),EI.GetDstNId());
					pGraph_temp->SetEDat(EI.GetSrcNId(),EI.GetDstNId(),EI.GetDat().Val);
				}

			}

			if(k!=0)//skip the first time
				GetNumOfReachableNodesFromSource(pGraph_temp,vSeedSet,vResult);

			for(auto NI = pGraph_temp->BegNI();NI<pGraph_temp->EndNI();NI++)
			{
				// make sure that v is not in Rg(S)
				if (k!=0) {
					auto result = std::find(vResult.begin(),vResult.end(), NI.GetId());
					if  (result!=vResult.end()) continue;
				}
				if(pGraph_temp->IsNode(NI.GetId()))
				{
					mNumOfNodesFromSource[NI.GetId()] = GetNumOfReachableNodesFromSource(pGraph_temp,NI.GetId());
					//cout<<"node "<<i<<", num of reachable nodes: "<<mNumOfNodesFromSource[i]<<endl;
				}

				// accumulate the influence spread
				mInfluenceSpread[NI.GetId()]+=mNumOfNodesFromSource[NI.GetId()];

			}
		}

		// select the i'th seed by finding u = argmax(mInfluenceSpread)
		auto it = std::max_element(mInfluenceSpread.begin(),mInfluenceSpread.end(),
			[&](std::pair<int,unsigned int> const& a, std::pair<int,unsigned int> const& b) {
				return a.second < b.second;
		}
		);
		int SeedID = it->first;
		//cout<<" max: "<<mInfluenceSpread[SeedID]<<endl;
		mInfluenceSpread.unsafe_erase(SeedID);
		vSeedSet.push_back(SeedID);
		cout << "Round "<< k <<": the most influential node is: "<<SeedID <<endl;

		pGraph_temp = CopyGraph(pGraph);
		pGraph_temp = GenerateDAG1(pGraph_temp, vSeedSet, 0.0);
		ParallelBPFromNode_LevelSynchronous(pGraph_temp, vSeedSet);
		double influence = InfluenceSpreadFromSeedNodes(pGraph_temp);
		cout << "Round "<<k<<": the influential spread is: "<<influence<<endl;

	}


}

void ParallelNewGreedIC(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int setSize, int numRounds, std::vector<int> &vSeedSet)
{

	auto pGraph_temp = TNodeEDatNet<TFlt, TFlt>::New();
	tbb::concurrent_unordered_map<int,unsigned int> mInfluenceSpread;
	std::map<int, unsigned int> mNumOfNodesFromSource;
	std::vector<int> vResult;
	int i,chunk = 70;
	static tbb::spin_mutex sMutex;

	for (int k=0; k<setSize ; ++k)
	{
		mInfluenceSpread.clear();
		mNumOfNodesFromSource.clear();
		vResult.clear();

#pragma omp parallel shared(chunk,mInfluenceSpread,pGraph) private(i,pGraph_temp,mNumOfNodesFromSource)
		{
#pragma omp for schedule(dynamic,chunk) nowait
			for (i =0;i<numRounds;++i)
			{

				//Copy all nodes of original graph to temp one
				pGraph_temp = TNodeEDatNet<TFlt, TFlt>::New();
				for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
					pGraph_temp->AddNode(NI.GetId());

				//1.Compute G' by removing each edge from G with probability 1 -p
				//2.Compute the set of vertices reachable from S
				//3.Compute as well the sets of vertices from v
				for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
				{
					double coinflip = (double) rand()/RAND_MAX;
					if( coinflip > 1.0-EI.GetDat().Val)
					{
						pGraph_temp->AddEdge(EI.GetSrcNId(),EI.GetDstNId());
						pGraph_temp->SetEDat(EI.GetSrcNId(),EI.GetDstNId(),EI.GetDat().Val);
					}

				}

				//cout<<"graph edges: "<<pGraph_temp->GetEdges()<<endl;

				if(k!=0)//skip the first time
					GetNumOfReachableNodesFromSource(pGraph_temp,vSeedSet,vResult);

				for(auto NI = pGraph_temp->BegNI();NI<pGraph_temp->EndNI();NI++)
				{
					// make sure that v is not in Rg(S)
					if (k!=0) {
						auto result = std::find(vResult.begin(),vResult.end(), NI.GetId());
						if  (result!=vResult.end()) continue;
					}

					if(pGraph_temp->IsNode(NI.GetId()))
					{
						mNumOfNodesFromSource[NI.GetId()] = GetNumOfReachableNodesFromSource(pGraph_temp,NI.GetId());
						//cout<<"node "<<i<<", num of reachable nodes: "<<mNumOfNodesFromSource[i]<<endl;
					}

					// accumulate the influence spread
#pragma omp critical
					mInfluenceSpread[NI.GetId()]+=mNumOfNodesFromSource[NI.GetId()];

				}

			}
		}

		// select the i'th seed by finding u = argmax(mInfluenceSpread)
		auto it = std::max_element(mInfluenceSpread.begin(),mInfluenceSpread.end(),
			[&](std::pair<int,unsigned int> const& a, std::pair<int,unsigned int> const& b) {
				return a.second < b.second;
		}
		);
		int SeedID = it->first;
		//cout<<" max: "<<mInfluenceSpread[SeedID]<<endl;
		mInfluenceSpread.unsafe_erase(SeedID);
		vSeedSet.push_back(SeedID);
		cout << "Round "<< k <<": the most influential node is: "<<SeedID <<endl;

		pGraph_temp = CopyGraph(pGraph);
		pGraph_temp = GenerateDAG1(pGraph_temp, vSeedSet, 0.0);
		ParallelBPFromNode_LevelSynchronous(pGraph_temp, vSeedSet);
		double influence = InfluenceSpreadFromSeedNodes(pGraph_temp);
		cout << "Round "<<k<<": the influential spread is: "<<influence<<endl;

	}


}

void ParallelNewGreedIC_Nested(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int setSize, int numRounds, std::vector<int> &vSeedSet)
{

	auto pGraph_temp = TNodeEDatNet<TFlt, TFlt>::New();
	tbb::concurrent_unordered_map<int,unsigned int> mInfluenceSpread;
	std::map<int, unsigned int> mNumOfNodesFromSource;
	std::vector<int> vResult;
	int i,j,chunk = 50,numNodes = pGraph->GetMxNId();
	static tbb::spin_mutex sMutex;

	for (int k=0; k<setSize ; ++k)
	{
		mInfluenceSpread.clear();
		mNumOfNodesFromSource.clear();
		vResult.clear();

#pragma omp parallel shared(chunk,mInfluenceSpread,pGraph) private(i,pGraph_temp)
		{
#pragma omp for schedule(dynamic) nowait
			for (i =0;i<numRounds;++i)
			{

				//Copy all nodes of original graph to temp one
				pGraph_temp = TNodeEDatNet<TFlt, TFlt>::New();
				for (auto NI = pGraph->BegNI(); NI < pGraph->EndNI(); NI++)
					pGraph_temp->AddNode(NI.GetId());

				//1.Compute G' by removing each edge from G with probability 1 -p
				//2.Compute the set of vertices reachable from S
				//3.Compute as well the sets of vertices from v
				for (auto EI = pGraph->BegEI(); EI < pGraph->EndEI(); EI++)
				{
					double coinflip = (double) rand()/RAND_MAX;
					if( coinflip > 1.0-EI.GetDat().Val)
					{
						pGraph_temp->AddEdge(EI.GetSrcNId(),EI.GetDstNId());
						pGraph_temp->SetEDat(EI.GetSrcNId(),EI.GetDstNId(),EI.GetDat().Val);
					}

				}

				//cout<<"graph edges: "<<pGraph_temp->GetEdges()<<endl;

				if(k!=0)//skip the first time
					GetNumOfReachableNodesFromSource(pGraph_temp,vSeedSet,vResult);

#pragma omp parallel for schedule(dynamic) shared(pGraph_temp,mInfluenceSpread) private(j,mNumOfNodesFromSource)
				for (j =0;j<numNodes;++j)
				{
					if (k!=0) {
						auto result = std::find(vResult.begin(),vResult.end(), j);
						if  (result!=vResult.end()) continue;
					}
					if(pGraph_temp->IsNode(j))
					{
						mNumOfNodesFromSource[j] = GetNumOfReachableNodesFromSource(pGraph_temp,j);
						//cout<<"node "<<i<<", num of reachable nodes: "<<mNumOfNodesFromSource[i]<<endl;
					}
#pragma omp critical
					mInfluenceSpread[j]+=mNumOfNodesFromSource[j];

				}

			}
		}

		// select the i'th seed by finding u = argmax(mInfluenceSpread)
		auto it = std::max_element(mInfluenceSpread.begin(),mInfluenceSpread.end(),
			[&](std::pair<int,unsigned int> const& a, std::pair<int,unsigned int> const& b) {
				return a.second < b.second;
		}
		);

		int SeedID = it->first;
		//cout<<" max: "<<mInfluenceSpread[SeedID]<<endl;
		mInfluenceSpread.unsafe_erase(SeedID);
		vSeedSet.push_back(SeedID);
		cout << "Round "<< k <<": the most influential node is: "<<SeedID <<endl;
		//Add the new seed and remove it from the map InfluenceSpread

	}


}