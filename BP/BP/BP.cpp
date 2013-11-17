#ifdef _USE_libDAI
#include <dai/alldai.h> // Include main libDAI header file
#include <dai/jtree.h>
#include <dai/bp.h>
#include <dai/decmap.h>

#pragma comment(lib, "libdai.lib")
#endif

#define __TBB_NO_IMPLICIT_LINKAGE	1
#include <tbb/concurrent_vector.h>
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
#endif

// http://snap.stanford.edu/snap/quick.html
#ifdef _DEBUG
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID, bool bDisplayInfo = false)
#else
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID)
#endif
{
	pGraph->SetNDat(sourceNodeID, 1.0);

	// Used to store the nodes which have been traversed during the breadth-first search traversal
	std::map<int, bool> visitedNodes;
	std::queue<int> queue;
	queue.push(sourceNodeID);
	visitedNodes[sourceNodeID] = true;
#ifdef _DEBUG
	int numEdge = 0;
#endif
	while(!queue.empty())
	{
		int nodeID = queue.front();
		// For debug only
#ifdef _DEBUG
		if(bDisplayInfo)
			std::cout << "Parsing node " << nodeID << " - Number of edges traversed: " << ++numEdge << std::endl;
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
}

#ifdef _DEBUG
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs, bool bDisplayInfo = false)
#else
void PropagateFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs)
#endif
{
	AddSuperRootNode(pGraph, vSeedIDs);

#ifdef _DEBUG
	PropagateFromNode(pGraph, INT_MAX, bDisplayInfo);
#else
	PropagateFromNode(pGraph, INT_MAX);
#endif

	// Remove the artificial super root node
	pGraph->DelNode(INT_MAX);
}

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
}

void ParallelBPFromNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs)
{
	AddSuperRootNode(pGraph, vSeedIDs);

	ParallelBPFromNode(pGraph, INT_MAX);

	// Remove the artificial super root node
	pGraph->DelNode(INT_MAX);
}

static void ParallelBPFromNode_1DPartitioning_UpdateNode(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph,
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

void ParallelBPFromNode_1DPartitioning(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID)
{
	pGraph->SetNDat(sourceNodeID, 1.0);

	// Like std::list, insertion of new items does not invalidate any iterators, nor change the order of items already in the map. Insertion and traversal may be concurrent.
	tbb::concurrent_hash_map<int, bool> visitedNodes;

	tbb::concurrent_queue<int> queue;
	queue.push(sourceNodeID);

	while(true)
	{
		int numNodesToProcess = queue.unsafe_size();
		// Parse level l
		if(numNodesToProcess > 50)
		{
			tbb::parallel_for(tbb::blocked_range<int>(0, numNodesToProcess, 50), 
				[&](const tbb::blocked_range<int>& r)
				{
					for (int i=r.begin();i!=r.end();++i)
						ParallelBPFromNode_1DPartitioning_UpdateNode(pGraph, visitedNodes, queue, true); 
				}
			);
		}
		else
		// Not enough elements to do the work concurrently
		{
			for (int i=0;i<numNodesToProcess;++i)
				ParallelBPFromNode_1DPartitioning_UpdateNode(pGraph, visitedNodes, queue, false); 
		}

		if(queue.empty())
			return;
	}
}

void ParallelBPFromNode_1DPartitioning(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, const std::vector<int> &vSeedIDs)
{
	AddSuperRootNode(pGraph, vSeedIDs);

	ParallelBPFromNode_1DPartitioning(pGraph, INT_MAX);

	// Remove the artificial super root node
	pGraph->DelNode(INT_MAX);
}

#if 0
void __stdcall PropagateFromNode_ThreadSafe(TPt<TNodeEDatNet<TFlt, TFlt>>& pGraph, int sourceNodeID)
{
	static tbb::spin_mutex sMutex;

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
#endif