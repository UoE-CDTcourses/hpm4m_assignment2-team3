#include <stdio.h>
#include "mpi.h"
#include <iostream>
using namespace std;
#include <cmath>

int main(){
    int N = 20, M = 5; 
    float T = 0.1, X = 1;

    float dt = T/N, dx = X/M;
    //initialize MPI    
    int nproc, rank;
    MPI_Status status;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    
    int size = (M+1-2)/nproc + 2;
    float U[size] = {0}, U_temp[size] = {0}; // Initialize two vectors, U and U_temp, U for store results, U_temp for updating.
    float U_result[nproc*size];     // Initialize U_result to store final results.

    // Give initial value to U.
    for(int i = 0; i < nproc ; i++){
        if(rank==i){            
            for(int m = 0; m < size ; m++){
                U[m] = sin(2*M_PI*(i*(size-2)+m)*dx) + 2*sin(5*M_PI*(i*(size-2)+m)*dx) + 3*sin(20*M_PI*(i*(size-2)+m)*dx);            
            }        
        }
    }

    // Update value of U.
    for(int j = 0; j < N; j++){
        for(int i = 0; i<nproc ; i++){
            if(rank == i){
                for(int m = 1; m<size-1 ; m++){
                    U_temp[m] = dt/(dx*dx) * (U[m-1] - 2*U[m] + U[m+1]) + U[m];
                }
                if(rank != 0){ // Except rank 0, send first value to current_rank-1  and receive value from current_rank-1
                    MPI_Send(&U_temp[1],1,MPI_FLOAT,(nproc+i-1)%nproc, i, MPI_COMM_WORLD); 
                    MPI_Recv(&U_temp[0],1,MPI_FLOAT,(nproc+i-1)%nproc, (nproc+i-1)%nproc, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                }
                if(rank != nproc-1){// Except last rank, send first value to current_rank+1  and receive value from current_rank+1
                    MPI_Send(&U_temp[size-2],1,MPI_FLOAT,(nproc+i+1)%nproc, i, MPI_COMM_WORLD);
                    MPI_Recv(&U_temp[size-1],1,MPI_FLOAT,(nproc+i+1)%nproc, (nproc+i+1)%nproc, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                }                

                MPI_Barrier(MPI_COMM_WORLD); // wait all process finish sending and receiving, update U.
                for(int m = 0; m<size ; m++){
                    U[m] = U_temp[m];
                    U_temp[m] = 0;
                }
            }
        } 
    }

    // Gather solution from all process.
    MPI_Gather(U, size, MPI_FLOAT, U_result, size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD); 
    MPI_Finalize();

    // Print the results
    if(rank==0){
        for(int i = 0; i < size ; i++){
            cout << U_result[i] << ' ';
        }
        for(int j = 1; j < nproc; j++){
            for(int i = 2; i < size ; i++){
                cout << U_result[j*size+i] << ' ';
            } 
        }
    }

    return 0;


}