#include <stdio.h>

#include <iostream>
using namespace std;
#include <cmath>

int main(){
    int N = 20, M = 5; 
    double T = 0.1, X = 1;

    float dt = T/N, dx = X/M;

    float U[M+1] = {0}, U_temp[M+1] = {0}; 

    for(int m = 1; m < M ; m++){
        U[m] = sin(2*M_PI*m*dx) + 2*sin(5*M_PI*m*dx) + 3*sin(20*M_PI*m*dx);
    }

    for(int n = 0; n < N ; n++){
        for(int m = 1; m < M ; m++){
            U_temp[m] = dt/(dx*dx) * (U[m-1] - 2*U[m] + U[m+1]) + U[m];
        }
        for(int m = 0; m < M ; m++){
            U[m] = U_temp[m];
            U_temp[m] = 0;
        }
    } 
     
    for(int i = 0; i < M+1 ; i++){
        cout << U[i] << ' ';
    }    

    return 0;


}

