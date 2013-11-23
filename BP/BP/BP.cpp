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
	pGraph->SetNDat(sourceNodeID, 1.0);

	// Used to store the nodes which have been traversed during the breadth-first search traversal
	std::map<int, bool> visitedNodes;
	std::queue<int> queue;
	queue.push(sourceNodeID);
	visitedNodes[sourceNodeID] = true;

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
		
			// Mark the child
			if(visitedNodes.insert(std::pair<int,bool>(iChildID,true)).second)
				queue.push(iChildID);
		}

		queue.pop();
	}
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

void MaxIncrementalInfluence(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int numRounds, std::vector<int> &vSeedSet)
{
	tbb::concurrent_unordered_map<int,double> mSpreadIncrement;;
	double influence = 0.0;
	static tbb::spin_mutex sMutex;

	//Failure of using PeerSeeds due to insufficient memory
	//std::map<int,std::vector<int> > mPeerSeeds;
	//std::map<int,TPt<TNodeEDatNet<TFlt, TFlt>> > mMIOAs;

	/* Initialization*/
	int numNodes = pGraph->GetMxNId();
	tbb::parallel_for(tbb::blocked_range<int>(0, numNodes, 100),
		[&](const tbb::blocked_range<int>& r)
		{
			for (int i=r.begin();i!=r.end();++i){
				if(pGraph->IsNode(i))
				{
					auto pGraph_v = MIOA(pGraph, i, 0);
					ResetGraphBelief(pGraph_v);
					ParallelBPFromNode_LevelSynchronous(pGraph_v, i);
					mSpreadIncrement[i]=InfluenceSpreadFromSeedNodes(pGraph_v);
					//mMIOAs.insert(std::make_pair(i,pGraph_v));
					//cout<< i <<endl;
				}
			}
		}
	);
/*
	//build PeerSeeds
	//Failure due to the insufficient memory
	for (int v =0; v<pGraph->GetNodes();++v)
		if(pGraph->IsNode(v))
			GetPeerSeeds(mMIOAs,v, mPeerSeeds[v]);
*/

	std::cout<<"--------------------------finish precomputation---------------------"<< std::endl;
	for (int i=0;i<numRounds;++i)
	{
		/* select the i'th seed by finding u = argmax(mSpreadIncrement)*/
		auto it = std::max_element(mSpreadIncrement.begin(),mSpreadIncrement.end(),
							[&](std::pair<int,double> const& a, std::pair<int,double> const& b) {
								 return a.second < b.second;
								 }
							);
		int SeedID = it->first;
		std::cout << SeedID << std::endl;

		/* calculate the current influence spread */
		vSeedSet.push_back(SeedID);
		pGraph = GenerateDAG1(pGraph, vSeedSet, 0.0);
		ParallelBPFromNode_LevelSynchronous(pGraph, vSeedSet);
		influence = InfluenceSpreadFromSeedNodes(pGraph);

		/*remove the newly selected node*/
		mSpreadIncrement.unsafe_erase(SeedID);

		/* update incremental influence spread for each round */
		double Delta_MAX = 0.0;
		tbb::parallel_for(tbb::blocked_range<int>(0, numNodes, 100),
			[&](const tbb::blocked_range<int>& r)
			{
				for (int i=r.begin();i!=r.end();++i)
				{
					/* exclude the nodes in seed set */
					auto result = std::find(vSeedSet.begin(),vSeedSet.end(), i);
					if (result != vSeedSet.end()) continue;
					//cout<<pGraph->IsNode(v)<<" , "<<mSpreadIncrement[v]<<endl;

					if(pGraph->IsNode(i) && mSpreadIncrement[i] > Delta_MAX)
					{
						/*different processors use different copied vSeedSet*/
						std::vector<int> vSeedSet_v = vSeedSet;
						vSeedSet_v.push_back(i);

						auto pGraph_v = GenerateDAG1(pGraph, vSeedSet_v, 0);
						ParallelBPFromNode_LevelSynchronous(pGraph_v, vSeedSet_v);
						mSpreadIncrement[i]=InfluenceSpreadFromSeedNodes(pGraph_v)-influence;
						if (mSpreadIncrement[i]> Delta_MAX)
						{
							tbb::spin_mutex::scoped_lock lock(sMutex);
							Delta_MAX = mSpreadIncrement[i];
						}
						vSeedSet_v.pop_back();
					}
				}

			}
		);


	}
}
