// FinalProject.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include<stdio.h>
#include <iostream>
#include<stdlib.h>
#include "mpi.h"
#include<vector>
#include<math.h>
#define HOUSE 100000
using namespace std;
typedef pair<double, double> pairs;
vector<pairs>houses;
pairs road;
double distance(pairs& house, pairs& rail) {
    
return pow( pow(house.first - rail.first, 2)+pow(house.second-rail.second,2),0.5); 

}
void initializevector(vector<pairs>&v,int NumHouse)
{
    for (int i = 0; i < NumHouse; i++) {
        int x = rand() % 100000 -50000; //random number between 1 and 10
        int y = rand() % 100000 + -50000;
        pairs temp = make_pair((double)x, (double)y);
        houses.push_back(temp);
    }
}
void displayvector(vector<pairs>&v,int NumHouse) {
    for (int i = 0; i < NumHouse; i++) {
        cout << i+1<<"th House->[" << v[i].first << "," << v[i].second << "] ";
    }
    cout << endl;
}
int main(int argc, char*argv[])
{
    int nprocessors,    //number of the processors
        rank,           // identifier
        mystart,
        minhouse,
        finalminhouse;   
    double start_time, end_time;
    double mindistance;   //cloesest distance from the station
    double finalminDistance;

    road.first = 0;
    road.second = 0;
    
    initializevector(houses, HOUSE);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocessors);    
    mystart = rank;
    mindistance = numeric_limits<double>::max();
    minhouse = -1;
    
    if (nprocessors == 1) {
        start_time = MPI_Wtime();       
        for (int i = mystart; i < HOUSE; i += nprocessors) {
            double temp = distance(houses[i], road);
            if (temp < mindistance) {
                minhouse = i;
                mindistance = temp;
            }
        }

        end_time = MPI_Wtime();
        cout << "The " << minhouse + 1 << "th house has minimum distance " << mindistance << endl;
        cout << "TIME -> " << end_time - start_time << endl;
        MPI_Finalize();
    }
    
    start_time = MPI_Wtime();
    for (int i = mystart; i < HOUSE; i += nprocessors) {
        double temp = distance(houses[i], road);
        if (temp < mindistance) {
            minhouse = i;
            mindistance = temp;
        }
    }

    if (rank > 0) {
        MPI_Recv(&finalminhouse, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&finalminDistance, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
       
        if (finalminDistance >= mindistance) {
            finalminhouse = minhouse;
            finalminDistance = mindistance;
        }

    }
    else { finalminhouse = minhouse; finalminDistance = mindistance; }

    MPI_Send(&finalminhouse, 1, MPI_INT, (rank + 1) % nprocessors, 0, MPI_COMM_WORLD);
    MPI_Send(&finalminDistance,1,MPI_DOUBLE, (rank + 1) % nprocessors, 0, MPI_COMM_WORLD);
    
    if(rank==0){
        MPI_Recv(&finalminhouse, 1, MPI_INT, nprocessors - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&finalminDistance, 1, MPI_DOUBLE, nprocessors - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        end_time = MPI_Wtime();

        //displayvector(houses,HOUSE);
        cout << "The " << finalminhouse+1 << "th house has minimum distance " << finalminDistance << endl;
        cout << "TIME -> " << end_time - start_time << endl;
    }
   
    MPI_Finalize();

}

