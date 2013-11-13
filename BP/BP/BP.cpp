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
	int superRootID = pGraph->GetNodes();
	pGraph->AddNode(superRootID);
	pGraph->SetNDat(superRootID, 1.0);
	for(auto it=vSeedIDs.begin(); it!=vSeedIDs.end(); ++it)
	{
		pGraph->AddEdge(superRootID, *it);
		pGraph->SetEDat(superRootID, *it, 1.0);
	}

#ifdef _DEBUG
	PropagateFromNode(pGraph, superRootID, bDisplayInfo);
#else
	PropagateFromNode(pGraph, superRootID);
#endif

	// Remove the artificial super root node
	pGraph->DelNode(superRootID);
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
	int superRootID = pGraph->GetNodes();
	pGraph->AddNode(superRootID);
	pGraph->SetNDat(superRootID, 1.0);
	for(auto it=vSeedIDs.begin(); it!=vSeedIDs.end(); ++it)
	{
		pGraph->AddEdge(superRootID, *it);
		pGraph->SetEDat(superRootID, *it, 1.0);
	}

	ParallelBPFromNode(pGraph, superRootID);

	// Remove the artificial super root node
	pGraph->DelNode(superRootID);
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
	int superRootID = pGraph->GetNodes();
	pGraph->AddNode(superRootID);
	pGraph->SetNDat(superRootID, 1.0);
	for(auto it=vSeedIDs.begin(); it!=vSeedIDs.end(); ++it)
	{
		pGraph->AddEdge(superRootID, *it);
		pGraph->SetEDat(superRootID, *it, 1.0);
	}

	ParallelBPFromNode_1DPartitioning(pGraph, superRootID);

	// Remove the artificial super root node
	pGraph->DelNode(superRootID);
}


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
