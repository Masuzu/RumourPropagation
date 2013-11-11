#pragma once

#include <Windows.h>
#include <vector>
#include "AEFunctorSlot.h"

#define __TBB_NO_IMPLICIT_LINKAGE	1
#include <tbb/concurrent_queue.h>

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

class WorkerThreadPool
{
protected:
	// Only used for initialization
	std::vector<unsigned int> m_vThreadIndex;
	unsigned int m_uiNumThreads;

public:
	// For each thread we have an event to kick off work or to signal that all the job is done
	std::vector<HANDLE> m_vBeginWorkEvent;
	std::vector<HANDLE> m_vEndWorkEvent;

	//! Used to add new tasks. Tasks are popped up when KickOffWOrkerThreads is called, until m_pTasks is empty.
	//! Queue up all the job at the beginning and then call KickOffWorkerThreads.
	tbb::concurrent_queue<AEFunctor> * m_pTasks;
	std::vector<HANDLE> m_vWorkerThreads;

	WorkerThreadPool(unsigned int numThreads);

	~WorkerThreadPool();

	void KickOffWorkerThreads();
	void WaitForCompletion();
};