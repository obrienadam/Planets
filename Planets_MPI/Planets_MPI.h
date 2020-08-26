#pragma once

#include "../Body/Body.h"
#include "mpi.h"

void Planets_MPI_Init();
  
void Planets_MPI_Finalize();

int N_Processors();

int This_Processor();

bool Main_Processor();

void Planets_MPI_Broadcast_Positions(Body *bodies, const int& NCl, const int& NCu);

void Planets_MPI_Broadcast_Positions2(Body *bodies, const int& NCl, const int& NCu);
