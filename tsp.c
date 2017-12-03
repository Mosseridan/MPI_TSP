#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#define TERMINATE_TAG 1
#define FETCH_JOB_TAG 2
#define ASSIGN_JOB_TAG 3
#define NEW_MIN_TAG 4

int myRank = 0;
int nProcs = 0;

MPI_Datatype jobType;
MPI_Datatype resultType;

void TSPRec(
	int level,
	int currentBound,
	int currentWeight,
	int* currentPath,
	int* minWeight,
	int* minPath,
	int* visited,
	int citiesNum,
	int* xCoord,
	int* yCoord);

// comptes the distants value (Weight) between cities i,j
int dist(int i, int j, int* xCoord, int* yCoord)
{
	return abs(xCoord[i] - xCoord[j]) + abs(yCoord[i] - yCoord[j]);
}

// finds the first and second minimum edge costs
// having an end at the vertex i
int findMinEdges(int* first, int* second, int i, int citiesNum, int* xCoord, int* yCoord)
{
	*first = INT_MAX;
	*second = INT_MAX;
	for (int j=0; j< citiesNum; j++)
	{
		if (i == j)
			continue;

		int distance = dist(i, j, xCoord, yCoord);
		// // printf("Distance between %d and %d is %d\n",i,j,distance);
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

void printPath(int length, int* path, int citiesNum, int* xCoord, int* yCoord){
	fflush(stdout);
	printf("!!![%d] ",myRank);
	for (int i = 0; i < length-1; i++)
	{
		// print the city (and its distance from the next city in the path)
		printf("%d (%d) ", path[i], dist(path[i],path[(i + 1) % citiesNum], xCoord, yCoord));
	}
	printf("%d \n",(path[length-1] % citiesNum));
	fflush(stdout);
}


void computeBranch(
	int i,
	int level,
	int currentBound,
	int currentWeight,
	int* currentPath,
	int* minWeight,
	int* minPath,
	int* visited,
	int citiesNum,
	int* xCoord,
	int* yCoord)
{

	int distance = dist(currentPath[level-1], i, xCoord, yCoord);

	// Consider next vertex if it is not same (diagonal
	// entry in adjacency matrix and not visited
	// already)
	if(distance != 0 && !visited[i])
	{

		// // printf("!![%d] Trying: ",myRank);
		// printcurrentPath(citiesNum,xCoord,yCoord,level,currentPath);
		// // printf(" -> (%d) %d \n",distance,i);

		currentWeight += distance;
		// // printf("@@@@ from: %d, to: %d, distance: %d, currentWeight: %d, newWeight: %d\n",currentPath[level-1], i,distance, currentWeight, newWeight);
		int firstMin, secondMin, prevFirstMin, prevSecondMin;
		findMinEdges(&firstMin, &secondMin, i, citiesNum, xCoord, yCoord);
		findMinEdges(&prevFirstMin, &prevSecondMin, currentPath[level-1], citiesNum, xCoord, yCoord);

		// different computation of curr_bound for
		// level 2 from the other levels
		if (level == 1){
			currentBound -= ((prevFirstMin + firstMin)/2);
		}
		else{
			currentBound -= ((prevSecondMin + firstMin)/2);
		}

		// newBound + currentWeight is the actual lower bound
		// for the node that we have arrived on
		// If current lower bound < minWeight, we need to explore
		// the node further
		if (currentBound + currentWeight < *minWeight)
		{
			// // printf("!![%d] %d (newBound) + %d (newWeight) = %d < %d (minWeight) *** Exploring branch ***\n\n",myRank,newBound,newWeight,newBound+newWeight,*pminWeight);
			currentPath[level] = i;
			visited[i] = 1;

			// call TSPRec for the next level
			TSPRec(level+1, currentBound, currentWeight, currentPath, minWeight, minPath, visited, citiesNum, xCoord, yCoord);
		}
		else{
			// // printf("!![%d] %d (newBound) + %d (newWeight) = %d >= %d (minWeight) $$ Cuting branch $$\n\n",myRank,newBound,newWeight,newBound+newWeight,*pminWeight);
		}

		// Also reset the visited array //TODO might want to deep copy at call instead
		memset(visited, 0, citiesNum*sizeof(int));
		for (int j=0; j<=level-1; j++)
		{
			visited[currentPath[j]] = 1;
		}
	}
}


void TSPRec(
	int level,
	int currentBound,
	int currentWeight,
	int* currentPath,
	int* minWeight,
	int* minPath,
	int* visited,
	int citiesNum,
	int* xCoord,
	int* yCoord)
{

	MPI_Status probeStatus, recvStatus;
	int flag = 0;
	// printf("!!![%d] probing\n\n",myRank);
	while (!MPI_Iprobe(MPI_ANY_SOURCE, NEW_MIN_TAG, MPI_COMM_WORLD, &flag, &probeStatus) && flag)
	{
		// fflush(stdout);
		// printf("###[%d] probed from %d\n\n",myRank, probeStatus.MPI_SOURCE);
		// fflush(stdout);
		struct{
			int minWeight;
			int minPath[citiesNum];
		} result;

		MPI_Recv(&result, 1, resultType, probeStatus.MPI_SOURCE, NEW_MIN_TAG, MPI_COMM_WORLD, &recvStatus);
		// fflush(stdout);
		// // printf("!!![%d] got path:\n",myRank);
		// // printPath(citiesNum, result.minPath, citiesNum, xCoord, yCoord);
		// printf("$$$[%d] got path with weight (%d) from [%d]\n\n",myRank,result.minWeight,recvStatus.MPI_SOURCE);
		// fflush(stdout);
		if (result.minWeight < *minWeight)
		{
			memcpy(minPath, result.minPath, citiesNum*sizeof(int));
			*minWeight = result.minWeight;

			// fflush(stdout);
			// // printf("@@@[%d] updated minPath to:\n",myRank);
			// // printPath(citiesNum, minPath, citiesNum, xCoord, yCoord);
			// printf("@@@[%d] updated minPath to weight (%d)\n\n",myRank,*minWeight);
			// fflush(stdout);
		}
	}

	// fflush(stdout);
	// // printf("\n@@[%d] TSPRec level: %d, currentWeight: %d, currentBound: %d, minWeight: %d\n",myRank,level,currentWeight, currentBound, minWeight);
	// printcurrentPath(level,currentPath);
	// printminPath(minPath);
	// printVisited();
	// // printf("\n\n");
	// fflush(stdout);

  // base case is when we have reached level N which
  // means we have covered all the nodes once
  if (level == citiesNum)
  {

		int distance = dist(currentPath[level-1], currentPath[0], xCoord, yCoord);

    // check if there is an edge from last vertex in
    // path back to the first vertex
    if (distance != 0)
    {
      // currentWeight has the total weight of the
    	currentWeight += distance;

      // Update final result and final path if
      // current result is better.
      if (currentWeight < *minWeight)
      {
				// // printf("\n!![%d] Found new shortest path that is %d long\n",myRank,currRes);
				memcpy(minPath, currentPath, citiesNum*sizeof(int));
				*minWeight = currentWeight;
				// printf("@@@[%d] broadcasting\n\n",myRank);
				MPI_Request requests[nProcs-2];
				MPI_Status statuses[nProcs-2];
				// MPI_Ibcast(minPath, 1, resultType, myRank, MPI_COMM_WORLD,&request);
				// MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm, MPI_Request *request)

				// fflush(stdout);
				// // printf("!!![%d] sending path:\n",myRank);
				// // printPath(citiesNum, minPath, citiesNum, xCoord, yCoord);
				// printf("!!![%d] sending path with weight (%d) \n\n",myRank,*minWeight);
				// fflush(stdout);

				for (int proc = 1; proc < myRank; proc++)
				{
					MPI_Isend(minWeight, 1, resultType, proc, NEW_MIN_TAG, MPI_COMM_WORLD, requests+proc-1);
				}
				for (int proc = myRank+1; proc < nProcs; proc++)
				{
					MPI_Isend(minWeight, 1, resultType, proc, NEW_MIN_TAG, MPI_COMM_WORLD, requests+proc-2);
				}
				MPI_Waitall(nProcs-2, requests, statuses);
      }
    }
    return;
	}

	// for any other level iterate for all vertices to
	// build the search space tree recursively
	for (int i = 0; i < citiesNum; i++)
	{
		computeBranch(i, level, currentBound, currentWeight, currentPath, minWeight, minPath, visited, citiesNum, xCoord, yCoord);
	}
}


// Master workflow
int doMaster(int citiesNum, int* minPath)
{
	MPI_Status status;

	int minWeight = INT_MAX;

	struct {
	 int startingCity;
	 int minWeight;
 } job;
 job.startingCity = 0;
 job.minWeight = INT_MAX;

	struct{
		int minWeight;
		int minPath[citiesNum];
	} result;
	result.minWeight = INT_MAX;
	memset(result.minPath, -1, citiesNum*sizeof(int));

	for(int city = 1; city < citiesNum; city++)
	{
		// printf("$$$[%d] handling routes that start with %d\n",myRank,city);
		MPI_Recv(&result, 1, resultType, MPI_ANY_SOURCE, FETCH_JOB_TAG, MPI_COMM_WORLD, &status);
		// printf("$$$[%d] received fetch request form %d\n",myRank,status.MPI_SOURCE);
		if(result.minWeight < minWeight)
		{
			// printf("$$$[%d] new minWeight: %d, old minWeight: %d\n",myRank,result.minWeight,minWeight);
			minWeight = result.minWeight;
			memcpy(minPath, result.minPath, citiesNum*sizeof(int));
		}
		// printf("$$$[%d] assigning routes starting with %d to %d\n",myRank,city,status.MPI_SOURCE);
		job.startingCity = city;
		job.minWeight = minWeight;
		MPI_Rsend(&job, 1, jobType, status.MPI_SOURCE, ASSIGN_JOB_TAG, MPI_COMM_WORLD);
	}

	for(int proc = 1; proc < nProcs; proc++)
	{
		MPI_Recv(&result, 1, resultType, MPI_ANY_SOURCE, FETCH_JOB_TAG, MPI_COMM_WORLD, &status);
		// printf("$$$[%d] received fetch request form %d\n",myRank,status.MPI_SOURCE);
		if(result.minWeight < minWeight)
		{
			// printf("$$$[%d] new minWeight: %d, old minWeight: %d\n",myRank,result.minWeight,minWeight);
			minWeight = result.minWeight;
			memcpy(minPath, result.minPath, citiesNum*sizeof(int));
		}
		// printf("$$$[%d] sending termination message to %d\n",myRank,status.MPI_SOURCE);
		MPI_Rsend(&job, 1, jobType, status.MPI_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD);
	}

	return minWeight;
}

// Worker workflow
int doWorker(int citiesNum, int* xCoord, int* yCoord)
{
	MPI_Request request;
	MPI_Status status;

	int terminate = 0;
	int minWeight = INT_MAX;
	int currentBound = 0;

	int visited[citiesNum];
	memset(visited, 0, citiesNum*sizeof(int));
	visited[0] = 1;

	int currentPath[citiesNum];
	memset(currentPath, -1, citiesNum*sizeof(int));
	currentPath[0] = 0;

	struct {
	 int startingCity;
	 int minWeight;
 } job;
 job.startingCity = 0;
 job.minWeight = INT_MAX;

	struct{
		int minWeight;
		int minPath[citiesNum];
	} result;
	result.minWeight = INT_MAX;
	memset(result.minPath, -1, citiesNum*sizeof(int));

	// Compute initial bound
	for (int i = 0; i < citiesNum; i++)
	{
		int first, second;
		findMinEdges(&first, &second, i, citiesNum, xCoord, yCoord);
		currentBound += (first + second);
	}
	// Rounding off the lower bound to an integer
	currentBound = (currentBound&1)? currentBound/2 + 1 : currentBound/2;

	while(!terminate)
	{
		MPI_Irecv(&job, 1, jobType, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
		MPI_Send(&result, 1, resultType, 0, FETCH_JOB_TAG, MPI_COMM_WORLD);
		MPI_Wait(&request, &status);
		if(status.MPI_TAG == ASSIGN_JOB_TAG)
		{
			// printf("## [%d] got job: %d %d\n",myRank,job.startingCity,job.minWeight);
			if(job.minWeight < minWeight)
			{
				minWeight = job.minWeight;
			}

			computeBranch(job.startingCity, 1, currentBound, 0, currentPath, &(result.minWeight), result.minPath, visited, citiesNum, xCoord, yCoord);

			if(result.minWeight < minWeight)
			{
				minWeight = result.minWeight;
			}
		}
		else if (status.MPI_TAG == TERMINATE_TAG)
		{
			// printf("## [%d] got termination message\n",myRank);
			terminate = 1;
		}
		else
		{
			// printf("## WTF?! got %d tag\n",status.MPI_TAG);
		}
	}
	return minWeight;
}


// The dynamic parellel algorithm main function.
int tsp_main(int citiesNum, int xCoord[], int yCoord[], int shortestPath[])
{
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcs);

	if(nProcs < 2)
	{
		return -1;
	}

	int blockLengths[2] = { 1, 1 };
	MPI_Aint displacements[2] = { 0, 4 };
	MPI_Datatype typelist[2] = { MPI_INT, MPI_INT };

	// Create the derived type for jobType
	MPI_Type_create_struct(2, blockLengths, displacements, typelist, &jobType);

	// Commit it so that it can be used for jobType
	MPI_Type_commit(&jobType);

	// change block length to fit resultType (the rest of the fields are identical to jobType)
	blockLengths[1] = citiesNum;

	// Create the derived type for resultType
	MPI_Type_create_struct(2, blockLengths, displacements, typelist, &resultType);

	// Commit it so that it can be used for resultType
	MPI_Type_commit(&resultType);


	if(myRank == 0)
	{
		return doMaster(citiesNum, shortestPath);
	}
	else
	{
		return doWorker(citiesNum, xCoord, yCoord);
	}
}
