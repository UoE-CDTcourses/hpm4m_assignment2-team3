#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

static const double PI = 3.1415926536; 

int main() { 

int M = 101;       // M number of intervals --> M+1 space points
int N = 10000;     // N time intervals
double T = 0.1;    //final time 
int J = 12;        // number of space points in each subinterval (process) 
double U[2][J];    // stores the numerical values of function U; two rows to also store values of previous time step 
double Usol[M+1];  // stores true solution 
double Unum[M+1];  // stores numerical solutiondouble Unum[M+1];
double User[2][M+1];  // stores serial numerical solution
double c[J-2];     // stores non-overlapping values in the numerical solution
double IC[M+1];    // stores initial condition vector
double dt = T/N;   // stepsize
double dx = 1./(M); 
double dtdx = dt/(dx*dx);
int rank, size;
double t_start, t_end, run_time;
int m; // for the error message
MPI_Comm comm;
comm  = MPI_COMM_WORLD;
MPI_Init(NULL,NULL);
MPI_Comm_rank(comm, &rank);            
MPI_Comm_size(comm, &size);
m = size * (J-2) + 1;
// Error message to check whether we have a correct choice of M, J and size
if(M != m and rank == 0) {
cout<<"Chosen "<< size <<" processes and "<< J << " grid points out of " << M << " for each process" << endl;
cout<<"This case cannot be handled, choose M, size and J such that M = size*(J-2)+2";
MPI_Abort(comm, 911);
MPI_Finalize();
return 0;
}

// ------Serial solution -----
User[0][0]=0, User[0][M]=0, User[1][0]=0, User[1][M]=0;
for(int m=1; m<M; ++m){
	User[0][m] = sin(2*PI*m*dx) + 2*sin(5*PI*m*dx) + 3*sin(20*PI*m*dx);
}
for(int i=1; i<=N; ++i){
	for (int m=1; m<M; ++m){
		User[1][m] = User[0][m] + dtdx * (User[0][m-1] - 2*User[0][m] + User[0][m+1]);	
	}
	for(int m=1; m<M; ++m){
		User[0][m] = User[1][m];
	}
}

// The root process generates initial conditions 
if (rank == 0) {
	cout<< "\ndx="<<dx<<", dt="<<dt<<", dt/dx²="<< dtdx<<endl;
	cout << "Employed " << size << " processes" << endl;
	for(int m=1; m<M; ++m){
		IC[m] = sin(2*PI*m*dx) + 2*sin(5*PI*m*dx) + 3*sin(20*PI*m*dx);
	}
	
	IC[0]=0, IC[M] = 0; // Dirichlet BC at time t_n
}

// Broadcast the initial condition to all processes
MPI_Bcast(IC, M+1, MPI_DOUBLE, 0, comm);

// Save the initial condition into the subinterval for the numerical solution for time tn
for (int i=0; i<J; i++) {
	U[0][i] = IC[rank*(J-2) + i];
}
// Dirichlet BC at time t_n+1 - only for the first and last process 
if (rank == 0) {
	U[1][0]=0; 
}
if (rank == size-1) {
	U[1][J-1]=0;
}
MPI_Barrier(comm);
t_start = MPI_Wtime();

// ----- NUMERICAL SCHEME ------------------------------------------------
// Use numerical scheme to obtain the future values of U on the M+1 space points
for(int i=1; i<=N; i++){	
	for (int m=1; m<J-1; m++){
		U[1][m] = U[0][m] + dtdx * (U[0][m-1] - 2*U[0][m] + U[0][m+1]);	
	}

	MPI_Barrier(comm);

	// ----- HALO SWAPPING -------------------------------------------
	if (rank > 0) {
		MPI_Ssend(&U[1][1], 1, MPI_LONG_DOUBLE, rank-1, 0, comm);
		MPI_Recv(&U[1][0], 1, MPI_LONG_DOUBLE, rank-1, 0, comm, MPI_STATUS_IGNORE);
	}

	if (rank < size-1) {
	    MPI_Recv(&U[1][J-1], 1, MPI_LONG_DOUBLE, rank+1, 0, comm, MPI_STATUS_IGNORE);
		MPI_Ssend(&U[1][J-2], 1, MPI_LONG_DOUBLE, rank+1, 0, comm);
	
	}
	// update "old" values	
	for(int m=0; m<J; m++){
		U[0][m] = U[1][m];
	}
}

MPI_Barrier(comm);
// Make shorter vectors for gather
for (int i=0; i<J-2; i++) {
	c[i] = U[1][i];
}
if (rank == size-1) {
	MPI_Ssend(&U[1][J-2], 1, MPI_DOUBLE, 0, 0, comm);
	MPI_Ssend(&U[1][J-1], 1, MPI_DOUBLE, 0, 1, comm);
}
if (rank == 0) {
	MPI_Recv(&Unum[M-1], 1, MPI_DOUBLE, size-1, 0, comm, MPI_STATUS_IGNORE);
	MPI_Recv(&Unum[M], 1, MPI_DOUBLE, size-1, 1, comm, MPI_STATUS_IGNORE);
}
MPI_Gather(&c, J-2, MPI_DOUBLE, Unum, J-2, MPI_DOUBLE, 0, comm);
MPI_Barrier(comm);
t_end = MPI_Wtime();
run_time = t_end - t_start;
cout << "Runtime: " << run_time << endl;
MPI_Barrier(comm);
if (rank == 0){
	// print out array entries of numerical solution next to true solution
	cout << "\nTrue and numerical values at M="<<M<<" space points at time T="<<T<<":"<<endl;
	cout << "\nTrue values           Parallel num.sol.     Serial num.sol.     Abs.Error\n"<<endl;
	for(int m=0; m<=M; ++m){
		Usol[m] = exp(-4*PI*PI*T)*sin(2*PI*m*dx) + 2*exp(-25*PI*PI*T)*sin(5*PI*m*dx) + 3*exp(-400*PI*PI*T)*sin(20*PI*M*dx);		
		cout << m << "  " << Usol[m] << "            " << Unum[m]<< "            " << User[1][m] <<  "            " << abs(Unum[m] - Usol[m]) << endl;
	}
}


MPI_Finalize();

}