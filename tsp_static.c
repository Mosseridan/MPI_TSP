#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

// comptes the distants value (wight) between cities i,j
int dist(int xi, int yi, int xj, int yj)
{
		return abs(xi - xj) + abs(yi - yj);
}

// finds the first and second minimum edge costs
// having an end at the vertex i
int findMinEdges(int* first, int* second, int i, int citiesNum, int xCoord[], int yCoord[])
{
	*first = INT_MAX;
	*second = INT_MAX;
	for (int j=0; j< citiesNum; j++)
	{
			if (i == j)
					continue;

			int distance = dist(xCoord[i], yCoord[i], xCoord[j], yCoord[j]);
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


// The static parellel algorithm main function.
int tsp_main(int citiesNum, int xCoord[], int yCoord[], int shortestPath[])
{
	int myRank, procCount;
	// Stores the final minimum weight of shortest tour.
	int shortestRes = INT_MAX;
	// visited[] keeps track of the already visited nodes
	// in a particular path
	int visited[citiesNum];

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &procCount);

	if(myRank == 0)
	{
		int currPath[citiesNum+1];

		// Calculate initial lower bound for the root node
		// using the formula 1/2 * (sum of first min +
		// second min) for all edges.
		// Also initialize the currPath and visited array
		int currBound = 0;
		memset(currPath, -1, sizeof(currPath));
		memset(visited, 0, sizeof(currPath));

		// Compute initial bound
		for (int i = 0; i < citiesNum; i++)
		{
			int first, second;
			findMinEdges(&first, &second, i, citiesNum, xCoord, yCoord);
			// printf("First and second of %d are: %d, %d\n",i,first,second);
			currBound += (first + second);
		}

		// We start at vertex 1 so the first vertex
		// in curr_path[] is 0
		visited[0] = 1;
		currPath[0] = 0;

		// // Call to TSPRec for curr_weight equal to
		// // 0 and level 1
		// TSPRec(adj, curr_bound, 0, 1, curr_path);
	}

	return -1;	//TODO
}
