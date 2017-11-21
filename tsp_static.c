#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>


int nCalls = 0;

// comptes the distants value (wight) between cities i,j
int dist(int xCoord[], int yCoord[], int i, int j)
{
	return abs(xCoord[i] - xCoord[j]) + abs(yCoord[i] - yCoord[j]);
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

		int distance = dist(xCoord, yCoord, i, j);
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


void printCurrPath(int citiesNum, int xCoord[], int yCoord[], int level, int currPath[]){
	for (int i = 0; i < level-1; i++)
	{
		// print the city (and its distance from the next city in the path)
		printf("%d (%d) ", currPath[i],
				abs(xCoord[currPath[i]] - xCoord[currPath[(i + 1) % citiesNum]]) +
				abs(yCoord[currPath[i]] - yCoord[currPath[(i + 1) % citiesNum]]) );
	}
	printf("%d",(currPath[level-1]% citiesNum));
}

void printShortestPath(int citiesNum, int xCoord[], int yCoord[], int shortestPath[]){
	int minPathLen = 0;
	for (int i = 0; i < citiesNum; i++)
	{
		minPathLen += abs(xCoord[shortestPath[i]] - xCoord[shortestPath[(i + 1) % citiesNum]]);
		minPathLen += abs(yCoord[shortestPath[i]] - yCoord[shortestPath[(i + 1) % citiesNum]]);
		// print the city (and its distance from the next city in the path)
		printf("%d (%d) ", shortestPath[i],
				abs(xCoord[shortestPath[i]] - xCoord[shortestPath[(i + 1) % citiesNum]]) +
				abs(yCoord[shortestPath[i]] - yCoord[shortestPath[(i + 1) % citiesNum]]) );
	}
	printf("%d  [%d]", shortestPath[0],minPathLen);
}

void printVisited(int citiesNum, int visited[]){
	for (int i = 0; i < citiesNum; i++)
	{
		printf("%d (%d) ",i, visited[i]);
	}
}

// function that takes as arguments:
// curr_bound -> lower bound of the root node
// curr_weight-> stores the weight of the path so far
// level-> current level while moving in the search
//         space tree
// curr_path[] -> where the solution is being stored which
//                would later be copied to final_path[]
void TSPRec(
	int citiesNum,
	int xCoord[],
	int yCoord[],
	int currBound,
	int currWeight,
	int level,
	int currPath[],
	int* pShortestRes,
	int shortestPath[],
	int visited[])
{

	nCalls++;
	printf("@@ TSPRec level: %d, currWeight: %d, currBound: %d",level,currWeight, currBound);
	printf("\n## currPath: ");
	printCurrPath(citiesNum,xCoord,yCoord,level,currPath);

	printf("\n## shortestPath: ");
	printShortestPath(citiesNum,xCoord,yCoord,shortestPath);
	printf(" [%d]",*pShortestRes);

	printf("\n## visited: ");
	printVisited(citiesNum, visited);
	printf("\n\n");


  // base case is when we have reached level N which
  // means we have covered all the nodes once
  if (level == citiesNum)
  {

		int distance = dist(xCoord, yCoord, currPath[level-1], currPath[0]);

    // check if there is an edge from last vertex in
    // path back to the first vertex
    if (distance != 0)
    {
      // curr_res has the total weight of the
      // solution we got
      int currRes = currWeight + distance;


      // Update final result and final path if
      // current result is better.
      if (currRes < *pShortestRes)
      {
				printf("\n!! Found new shortest path that is %d long\n",currRes);
        // copyToFinal(currPath); // TODO: make parallel
				for (int i= 0; i < citiesNum; i++)
						shortestPath[i] = currPath[i];
				shortestPath[citiesNum] = currPath[0];

        *pShortestRes = currRes;  // TODO: make parallel
      }
    }
    return;
	}

	// for any other level iterate for all vertices to
	// build the search space tree recursively
	for (int i = 0; i < citiesNum; i++)
	{
		int distance = dist(xCoord, yCoord, currPath[level-1], i);

	  // Consider next vertex if it is not same (diagonal
    // entry in adjacency matrix and not visited
    // already)
		if(distance != 0 && !visited[i])
    {

			printf("!! Trying: ");
			printCurrPath(citiesNum,xCoord,yCoord,level,currPath);
			printf(" -> (%d) %d \n",distance,i);

			int newBound = currBound;
			int newWeight = currWeight + distance;
			// printf("@@@@ from: %d, to: %d, distance: %d, currWeight: %d, newWeight: %d\n",currPath[level-1], i,distance, currWeight, newWeight);
			int firstMin, secondMin, prevFirstMin, prevSecondMin;
			findMinEdges(&firstMin, &secondMin, i, citiesNum, xCoord, yCoord);
			findMinEdges(&prevFirstMin, &prevSecondMin, currPath[level-1], citiesNum, xCoord, yCoord);

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
      if (newBound + newWeight < *pShortestRes)
      {
				printf("!! %d (newBound) + %d (newWeight) = %d < %d (shortestRes) *** Exploring branch ***\n\n",newBound,newWeight,newBound+newWeight,*pShortestRes);
				currPath[level] = i;
        visited[i] = 1;

        // call TSPRec for the next level
        TSPRec(citiesNum, xCoord, yCoord, newBound, newWeight, level+1, currPath, pShortestRes, shortestPath, visited); //TODO: currPath and visited might need deep copy
      }
			else{
				printf("!! %d (newBound) + %d (newWeight) = %d >= %d (shortestRes) $$ Cuting branch $$\n\n",newBound,newWeight,newBound+newWeight,*pShortestRes);
			}

      // Also reset the visited array //TODO might want to deep copy at call instead
      memset(visited, 0, citiesNum*sizeof(int));
      for (int j=0; j<=level-1; j++)
			{
				visited[currPath[j]] = 1;
			}
		}
	}
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
		memset(currPath, -1, citiesNum*sizeof(int));
		memset(visited, 0, citiesNum*sizeof(int));

		// Compute initial bound
		for (int i = 0; i < citiesNum; i++)
		{
			int first, second;
			findMinEdges(&first, &second, i, citiesNum, xCoord, yCoord);
			currBound += (first + second);
		}

		// Rounding off the lower bound to an integer
		currBound = (currBound&1)? currBound/2 + 1 : currBound/2;


		// We start at vertex 1 so the first vertex
		// in curr_path[] is 0
		visited[0] = 1;
		currPath[0] = 0;

		// Call to TSPRec for curr_weight equal to
		// 0 and level 1
		TSPRec(citiesNum, xCoord, yCoord, currBound, 0, 1, currPath, &shortestRes, shortestPath, visited);

		int factorial = 1;
		for(int i = 1; i <= citiesNum; i++){
			factorial = factorial * i;
		}
		int speedup = factorial/nCalls;
		printf("\n\n@@ Total number of calls: %d instead of %d. ~ %d Times faseter\n\n",nCalls,factorial,speedup,speedup*100);

	}

	return shortestRes;
}
