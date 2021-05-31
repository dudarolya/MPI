/*   1.  Each process calculates "its" interval of
         integration.
     2.  Each process estimates the integral of f(x)
         over its interval.
     3a. Each process != 0 sends its integral to 0.
     3b. Process 0 sums the calculations received from
         the individual processes and prints the result.*/

#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
using namespace std;

double Rectangles(double left_endpt, double right_endpt, int trap_count, double base_len);
double f(double x);

int main(int argc, char** argv) 
{

    int procRank, procNum, n = 400000000, _n;
    double a = 0.0, b = 1.0, h, local_a, local_b;
    double local_sum = 0.0, sum = 0.0, Start, Finish, Duration;
    int source;

    MPI_Init(&argc, &argv);

    Start = MPI_Wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);

    h = (b - a) / n;
    _n = n / procNum;  //the number of rectangles

    local_a = a + procRank * _n * h;
    local_b = local_a + _n * h;
    local_sum = Rectangles(local_a, local_b, _n, h);

    std::cout << "Process number: " << procRank << " Local sum: " << local_sum << endl;// "Left bound: " << local_a << " Right bound: " << local_b << std::endl;

    //integrals calculated by each process
    if (procRank != 0) 
    {
        MPI_Send(&local_sum, 1, MPI_DOUBLE, 0, 0,
            MPI_COMM_WORLD);
    }
    else 
    {
        sum = local_sum;
        for (source = 1; source < procNum; source++) 
        {
            MPI_Recv(&local_sum, 1, MPI_DOUBLE, source, 0,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += local_sum;
        }
    }

    Finish = MPI_Wtime();
    Duration = Finish - Start;

    //TestResult(pMatrix, pVector, pResult, Size);
    if (procRank == 0)
    {
        printf("Time of execution = %f\n", Duration);
        printf("With n = %d rectangles, our estimate\n", n);
        printf("of the integral from 0 to 1 = %.15e\n", sum);
    }

    MPI_Finalize();

    return 0;
}
//Estimate of integral from a to b using rectangle_count rectangles
double Rectangles(double a, double b, int rectangle_count, double base_len) 
{

    double estimate, x;
    int i;

    estimate = f(a);

    for (i = 1; i <= rectangle_count - 1; i++) 
    {
        x = a + i * base_len;
        estimate += f(x);
    }
    estimate = estimate * base_len;

    return estimate;
}

double f(double x) 
{
    return sin(x + x * x * x);
}