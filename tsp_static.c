#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#define CITIES_NUM_TAG 1
#define X_COORD_TAG 2
#define Y_COORD_TAG 3
#define CURR_BOUND_TAG 4

typedef struct {
    double real,imag;
} Path;


int nCalls = 0;
int citiesNum = 0;
int* xCoord = 0;
int* yCoord = 0;
int shortestRes = INT_MAX;
int* shortestPath = 0;
int* visited = 0;
int myRank = 0;
int nProcs = 0;
int nIters = 0;

void TSPRec(int currBound, int currWeight, int level, int currPath[]);
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


void printCurrPath(int level, int currPath[]){
	fflush(stdout);
	printf("\n##[%d] currPath: ",myRank);
	for (int i = 0; i < level-1; i++)
	{
		// print the city (and its distance from the next city in the path)
		printf("%d (%d) ", currPath[i], dist(currPath[i],currPath[(i + 1) % citiesNum]));
	}
	printf("%d\n",(currPath[level-1]% citiesNum));
	fflush(stdout);
}

void printShortestPath(){
	int minPathLen = 0;
	fflush(stdout);
	printf("\n##[%d] shortestPath: ",myRank);
	for (int i = 0; i < citiesNum; i++)
	{
		// printf("\n\n#### [%d] %d %d %d\n\n",myRank,i, shortestPath, shortestPath[i]);
		fflush(stdout);
		minPathLen += dist(shortestPath[i], shortestPath[(i + 1) % citiesNum]);
		// print the city (and its distance from the next city in the path)

		printf("%d (%d) ", shortestPath[i], dist(shortestPath[i],shortestPath[(i + 1) % citiesNum]));
	}
	printf("%d  [%d]\n", shortestPath[0],minPathLen);
	fflush(stdout);
}


void printVisited(){
	fflush(stdout);
	printf("\n##[%d] visited: ",myRank);
	for (int i = 0; i < citiesNum; i++)
	{
		printf("%d (%d) ",i, visited[i]);
	}
	printf("\n");
	fflush(stdout);
}


void computeBranch(int i, int currBound, int currWeight, int level, int currPath[])
{

	int distance = dist(currPath[level-1], i);

	// Consider next vertex if it is not same (diagonal
	// entry in adjacency matrix and not visited
	// already)
	if(distance != 0 && !visited[i])
	{

		// printf("!![%d] Trying: ",myRank);
		// printCurrPath(citiesNum,xCoord,yCoord,level,currPath);
		// printf(" -> (%d) %d \n",distance,i);

		int newBound = currBound;
		int newWeight = currWeight + distance;
		// printf("@@@@ from: %d, to: %d, distance: %d, currWeight: %d, newWeight: %d\n",currPath[level-1], i,distance, currWeight, newWeight);
		int firstMin, secondMin, prevFirstMin, prevSecondMin;
		findMinEdges(&firstMin, &secondMin, i);
		findMinEdges(&prevFirstMin, &prevSecondMin, currPath[level-1]);

		// different computation of curr_bound for
		// level 2 from the other levels
		if (level == 1){
			newBound -= ((prevFirstMin + firstMin)/2);
		}
		else{
			newBound -= ((prevSecondMin + firstMin)/2);
		}

		// newBound + currWeight is the actual lower bound
		// for the node that we have arrived on
		// If current lower bound < shortestRes, we need to explore
		// the node further
		if (newBound + newWeight < shortestRes)
		{
			// printf("!![%d] %d (newBound) + %d (newWeight) = %d < %d (shortestRes) *** Exploring branch ***\n\n",myRank,newBound,newWeight,newBound+newWeight,*pShortestRes);
			currPath[level] = i;
			visited[i] = 1;

			// call TSPRec for the next level
			TSPRec(newBound, newWeight, level+1, currPath); //TODO: currPath and visited might need deep copy
		}
		else{
			// printf("!![%d] %d (newBound) + %d (newWeight) = %d >= %d (shortestRes) $$ Cuting branch $$\n\n",myRank,newBound,newWeight,newBound+newWeight,*pShortestRes);
		}

		// Also reset the visited array //TODO might want to deep copy at call instead
		memset(visited, 0, citiesNum*sizeof(int));
		for (int j=0; j<=level-1; j++)
		{
			visited[currPath[j]] = 1;
		}
	}
}


// function that takes as arguments:
// curr_bound -> lower bound of the root node
// curr_weight-> stores the weight of the path so far
// level-> current level while moving in the search
//         space tree
// curr_path[] -> where the solution is being stored which
//                would later be copied to final_path[]
void TSPRec(int currBound, int currWeight, int level, int currPath[])
{

	nCalls++;
	// fflush(stdout);
	// printf("\n@@[%d] TSPRec level: %d, currWeight: %d, currBound: %d, shortestRes: %d\n",myRank,level,currWeight, currBound, shortestRes);
	// printCurrPath(level,currPath);
	// printShortestPath(shortestPath);
	// printVisited();
	// printf("\n\n");
	// fflush(stdout);

  // base case is when we have reached level N which
  // means we have covered all the nodes once
  if (level == citiesNum)
  {

		int distance = dist(currPath[level-1], currPath[0]);

    // check if there is an edge from last vertex in
    // path back to the first vertex
    if (distance != 0)
    {
      // curr_res has the total weight of the
      // solution we got
      int currRes = currWeight + distance;


      // Update final result and final path if
      // current result is better.
      if (currRes < shortestRes)
      {
				// printf("\n!![%d] Found new shortest path that is %d long\n",myRank,currRes);
				for (int i= 0; i < citiesNum; i++)
						shortestPath[i] = currPath[i];
        shortestRes = currRes;
      }
    }
    return;
	}

	// for any other level iterate for all vertices to
	// build the search space tree recursively
	for (int i = 0; i < citiesNum; i++)
	{
		computeBranch(i, currBound, currWeight, level, currPath);
	}
}

// The static parellel algorithm main function.
int tsp_main(int _citiesNum, int _xCoord[], int _yCoord[], int _shortestPath[])
{
	MPI_Request request;
	MPI_Status status;

	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	struct {
		 int val;
		 int rank;
	} finalRes, myRes;
	finalRes.val =  INT_MAX;
	finalRes.rank = myRank;
	myRes.val = INT_MAX;
	myRes.rank = myRank;
	int currBound = 0;

	if(myRank == 0)
	{
		citiesNum = _citiesNum;
		xCoord = malloc(citiesNum*sizeof(int));
		yCoord = malloc(citiesNum*sizeof(int));

		memcpy(xCoord, _xCoord, citiesNum*sizeof(int));
		memcpy(yCoord, _yCoord, citiesNum*sizeof(int));

		for(int i = 1; i < nProcs; i++)
		{
			MPI_Isend(&citiesNum, 1, MPI_INT, i, CITIES_NUM_TAG, MPI_COMM_WORLD, &request);
			MPI_Isend(xCoord, citiesNum, MPI_INT, i, X_COORD_TAG, MPI_COMM_WORLD, &request);
			MPI_Isend(yCoord, citiesNum, MPI_INT, i, Y_COORD_TAG, MPI_COMM_WORLD, &request);
		}

		// Compute initial bound
		for (int i = 0; i < citiesNum; i++)
		{
			int first, second;
			findMinEdges(&first, &second, i);
			currBound += (first + second);
		}

		// Rounding off the lower bound to an integer
		currBound = (currBound&1)? currBound/2 + 1 : currBound/2;

		for(int i = 1; i < nProcs; i++)
		{
			MPI_Isend(&currBound, 1, MPI_INT, i, CURR_BOUND_TAG, MPI_COMM_WORLD, &request);
		}
	}
	else {
		MPI_Recv(&citiesNum, 1, MPI_INT, 0, CITIES_NUM_TAG, MPI_COMM_WORLD, &status);
		xCoord = malloc(citiesNum*sizeof(int));
		yCoord = malloc(citiesNum*sizeof(int));
		MPI_Irecv(xCoord, citiesNum, MPI_INT, 0, X_COORD_TAG, MPI_COMM_WORLD, &request);
		MPI_Irecv(yCoord, citiesNum, MPI_INT, 0, Y_COORD_TAG, MPI_COMM_WORLD, &request);
		MPI_Irecv(&currBound, 1, MPI_INT, 0, CURR_BOUND_TAG, MPI_COMM_WORLD, &request);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	int currPath[citiesNum+1];
	shortestPath = _shortestPath;
	// shortestPath = malloc(citiesNum*sizeof(int));
	visited = malloc(citiesNum*sizeof(int));
	memset(currPath, -1, citiesNum*sizeof(int));
	// memset(shortestPath, 0, citiesNum*sizeof(int));
	memset(visited, 0, citiesNum*sizeof(int));
	// We start at vertex 1 so the first vertex
	// in curr_path[] is 0
	visited[0] = 1;
	currPath[0] = 0;
	nIters = citiesNum/nProcs;


	// fflush(stdout);
	// printf("\n\n\n$@#@#@#@$$ [%d] %d\n",myRank,citiesNum);
	// for(int i = 0; i < citiesNum; i++){
	// 	printf("%d:%d ",xCoord[i], yCoord[i]);
	// }
	// printf("\n");
	// for(int i = 0; i < citiesNum; i++){
	// 	printf("%d",shortestPath[i]);
	// }
	// printf("\n\n");
	// fflush(stdout);

	int start = myRank * nIters;
	int end = (myRank == nProcs-1 ? citiesNum : start + nIters);
	for(int i = start; i < end; i++)
	{
		computeBranch(i, currBound, 0, 1, currPath);
	}

	myRes.val = shortestRes;
	MPI_Allreduce(&myRes, &finalRes, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
	MPI_Bcast(shortestPath, citiesNum, MPI_INT,finalRes.rank, MPI_COMM_WORLD);

	free(xCoord);
	free(yCoord);
	// free(shortestPath);
	free(visited);
	MPI_Barrier(MPI_COMM_WORLD);
	return finalRes.val;
}
