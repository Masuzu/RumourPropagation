#include "AEWorkerThread.h"

struct WorkerThreadInfo
{
	unsigned int iThreadIndex;
	WorkerThreadPool *pThreadPool;
};

static DWORD WINAPI DoWork( LPVOID lpParam )
{
	WorkerThreadInfo *pThreadInfo = (WorkerThreadInfo*)(lpParam);
	unsigned int iThreadIndex = pThreadInfo->iThreadIndex;
	WorkerThreadPool *pThreadPool = pThreadInfo->pThreadPool;

	for(;;)
	{
		WaitForSingleObject(pThreadPool->m_vBeginWorkEvent.at(iThreadIndex), INFINITE);

		// Pop up tasks

		UINT meter = 0;
		AEFunctor FunctorTask;
		while(pThreadPool->m_pTasks->try_pop(FunctorTask))
		{
			//FunctorTask(... /*params*/);
			meter++;
		}

		SetEvent(pThreadPool->m_vEndWorkEvent.at(iThreadIndex));
	}

	return 0;
}

WorkerThreadPool::WorkerThreadPool(unsigned int numThreads)
	: m_uiNumThreads(numThreads)
{
	m_pTasks = new tbb::concurrent_queue<AEFunctor>();

	for(UINT i = 0; i < numThreads; i++)
	{
		m_vThreadIndex.push_back(i);

		HANDLE hEndDeferredRenderingEvent = CreateEvent( NULL, FALSE, FALSE, NULL );
		m_vEndWorkEvent.push_back(hEndDeferredRenderingEvent);

		HANDLE hBeginDeferredRenderingEvent = CreateEvent( NULL, FALSE, FALSE, NULL );
		m_vBeginWorkEvent.push_back(hBeginDeferredRenderingEvent);
	}

	// Create the worker threads.
	for(UINT i = 0; i < numThreads; i++)
	{
		DWORD ThreadID;
		HANDLE hThread = CreateThread( 
			NULL,       // default security attributes
			0,          // default stack size
			(LPTHREAD_START_ROUTINE) DoWork, 
			&m_vThreadIndex[i],       // no thread function arguments
			0,          // default creation flags
			&ThreadID); // receive thread identifier

		m_vWorkerThreads.push_back(hThread);
	}
}

WorkerThreadPool::~WorkerThreadPool()
{
	delete m_pTasks;
}

void WorkerThreadPool::KickOffWorkerThreads()
{

	for(auto it = m_vBeginWorkEvent.begin(); it != m_vBeginWorkEvent.end(); ++it)
		SetEvent(*it);
}

void WorkerThreadPool::WaitForCompletion()
{
	for(auto itEvent = m_vEndWorkEvent.begin(); itEvent != m_vEndWorkEvent.end(); ++itEvent)
		WaitForSingleObject(*itEvent, INFINITE);
}