#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int m_row=20, m_column=20, zed=30, matrix[20][20], visited[20], best_path[20];
int best_cost=9999999, size=20;

void dfs(int city, int visited_in[], int path_in[], int path_i_in, int cost_in) {
  if (cost_in < best_cost) {
    int *visited=calloc(sizeof(int),size+1);
    int *path=calloc(sizeof(int),size+1);
    int path_i=path_i_in, cost=cost_in, i;

    for (i=0; i<size; i++) {
      visited[i]=visited_in[i];
      path[i]=path_in[i];
    }

    visited[city]=1;
    path[path_i]=city;
    path_i++;
    int leaf=0;

    for(i=0; i<size; i++) {
      if(visited[i]==0) {
        leaf++;
        dfs(i, visited, path, path_i, cost+matrix[city][i]);
      }
    }

    if (leaf == 0) {
      cost+=matrix[city][0];
      path[path_i]=0;
      path_i++;
      if(cost < best_cost) {
        //printf(“Found new best cost: %i\n”, cost);
        best_cost=cost;
        for(i=0; i<size; i++)
        best_path[i]=path[i];
      }
    }

    free(visited);
    free(path);
  }
}


int main(int argc, char* argv[]) {
  int rank, p, source, dest;
  int tag=0;
  MPI_Status status;
  MPI_Init(0, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  srand( time(NULL) );

  if (rank == 0) {
    int i, j;
    for(i=0; i<m_row; i++)
    for(j=0; j<m_column; j++)
    matrix[i][j]=0;

    for ( i=0; i<m_row; i++ ) {
      for( j=0; j<i; j++ ) {
        if (i != j) {
          int temp=(rand()%zed)+1;
          matrix[i][j]=temp;
          matrix[j][i]=temp;
        }
      }
    }

    for(i=1; i<p; i++)
      MPI_Send(&matrix[0][0], size*size, MPI_LONG, i, 0, MPI_COMM_WORLD);

    printf(“Matrix, %lix%li, Max Int: %li\n”, m_row,m_column,zed);

    for(i=0; i<m_row; i++) {
      for(j=0; j<m_column; j++)
      printf(“%li\t”, matrix[i][j]);
      printf(“\n”);
      fflush(NULL);
    }

    printf(“\n”);
    int winner;
    int node_array[p-1];
    int node_array_i=0;
    for(i=0; i<p-1; i++)
    node_array[i]=i+1;

    for(i=1; i<size; i++) {
      int temp_best_cost, node;
      node=node_array[node_array_i];
      if (node_array_i < p-2)
      node_array_i++;
      else
      node_array_i=0;
      int *temp_best_path=calloc(sizeof(int),size+1);
      MPI_Recv(&temp_best_cost, 1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&temp_best_path[0], size+1, MPI_INT, node, 0, MPI_COMM_WORLD, &status);
      if(temp_best_cost < best_cost) {
        winner=node;
        best_cost=temp_best_cost;
        for(j=0; j<size+1; j++)
        best_path[j]=temp_best_path[j];
      }
      MPI_Send(&best_cost, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
    }

    printf(“Best Path Found by node %i:\n”, winner);
    printf(“%i”,best_path[0]);
    for(i=1; i<size+1; i++)
    printf(” –> %i”,best_path[i]);
    printf(“\nBest Cost Found: %i\n”, best_cost);
  }
  else {
    MPI_Recv(&(matrix[0][0]), m_row*m_column, MPI_LONG, 0, 0, MPI_COMM_WORLD, &status);
    int i;

    for(i=rank; i<size; i+=(p-1)) {
      int *visited=calloc(sizeof(int),size+1);
      int *path=calloc(sizeof(int),size+1);
      int cost=matrix[0][i], path_i=1;
      path[0]=0;
      visited[0]=1;
      dfs(i, visited, path, path_i, cost);
      MPI_Send(&best_cost, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Send(&best_path[0], size+1, MPI_INT, 0, 0, MPI_COMM_WORLD);
      MPI_Recv(&best_cost, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
      free(visited);
      free(path);
    }
  }
  MPI_Finalize();
  return 0;
}
