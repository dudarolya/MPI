#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <iostream>
#include <iomanip>
using namespace std;

int main(int argc, char* argv[])
{
    double h, x, lBound = 0.0, rBound = 1.0, IntervalLength = 0.0, currArea = 0.0, sum;

    char processor_name[MPI_MAX_PROCESSOR_NAME], (*procNames)[MPI_MAX_PROCESSOR_NAME];
    int intervalNum = 0, procNum, procRank, namelen, proc = 0, c = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    MPI_Get_processor_name(processor_name, &namelen);

    procNames = (char(*)[128]) malloc(procNum * MPI_MAX_PROCESSOR_NAME);

    MPI_Gather(processor_name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR,
        procNames, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (procRank == 0)
    {
         intervalNum = 10000000;

        cout << "Number of rectangles = " << intervalNum << endl;

        for (proc = 0; proc < procNum; ++proc)
            printf("Process %d\n", proc);
    }

    MPI_Bcast(&intervalNum, 1, MPI_INT, 0, MPI_COMM_WORLD);

    h = (rBound - lBound) / intervalNum;

    for (int i = procRank + 1; i <= intervalNum; i += procNum) 
    {
        x = h * i;
        IntervalLength += sin(x + x * x * x);
        //c++;
    }
    //cout << c;

    currArea = h * IntervalLength;

    MPI_Reduce(&currArea, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (procRank == 0) 
        cout << std::setprecision(7) << "Result: " << sum << endl;

    MPI_Finalize();

    return 0;
}