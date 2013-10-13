#include "Utilities.h"
#include <Snap.h>

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