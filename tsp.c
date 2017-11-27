#include <mpi.h>

int myRank = 0;
int nProcs = 0;

int citiesNum = 0;
int* xCoord = 0;
int* yCoord = 0;

int shortestRes = INT_MAX;
int* shortestPath = 0;


// comptes the distants value (wight) between cities i,j
int dist(int i, int j)
{
	return abs(xCoord[i] - xCoord[j]) + abs(yCoord[i] - yCoord[j]);
}

// finds the first and second minimum edge costs
// having an end at the vertex i
int findMinEdges(int* first, int* second, int i)
{
	*first = INT_MAX;
	*second = INT_MAX;
	for (int j=0; j< citiesNum; j++)
	{
		if (i == j)
			continue;

		int distance = dist(i, j);
		// printf("Distance between %d and %d is %d\n",i,j,distance);
		if (distance <= *first)
		{
			*second = *first;
			*first = distance;
		}
		else if (distance <= *second)
		{
			*second = distance;
		}
	}
	return *first;
}


// Master workflow
int doMaster()
{

}

// Worker workflow
int doWorker()
{
	int terminate = 0;
	while(!terminate)
	{

	}
	return 0;
}



// The dynamic parellel algorithm main function.
int tsp_main(int _citiesNum, int _xCoord[], int _yCoord[], int _shortestPath[])
{
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	if(nProcs < 2)
	{
		return -1;
	}

	citiesNum = _citiesNum;
	xCoord = _xCoord;
	yCoord = _yCoord;
	shortestPath = _shortestPath;

	int currBound = 0;

	if(myRank == 0)
	{
		return doMaster();
	}
	else
	{
		return doWorker();
	}
}
