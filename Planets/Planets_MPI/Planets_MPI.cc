#include "Planets_MPI.h"
#include <iostream>
#include "mpi.h"

void Planets_MPI_Init(){

  MPI::Init();

}

  
void Planets_MPI_Finalize(){

  MPI::Finalize();

}

int N_Processors(){

  return MPI::COMM_WORLD.Get_size();

}

int This_Processor(){

  return MPI::COMM_WORLD.Get_rank();

}

bool Main_Processor(){

  if(MPI::COMM_WORLD.Get_rank() == 0){
    return true;
  } else {
    return false;
  }

}

void Planets_MPI_Broadcast_Positions(Body *bodies, const int& NCl, const int& NCu){

  int NCl_root, NCu_root, i, j;
  double buffer[6];

  for(i = 0; i < N_Processors(); ++i){

    if(i == This_Processor()){
      NCl_root = NCl;
      NCu_root = NCu;
    }

    MPI::COMM_WORLD.Bcast(&(NCl_root), 1, MPI::INT, i);
    MPI::COMM_WORLD.Bcast(&(NCu_root), 1, MPI::INT, i);

    for(j = NCl_root; j <= NCu_root; ++j){

      if(This_Processor() == i){

	buffer[0] = bodies[j].Xc.x;
	buffer[1] = bodies[j].Xc.y;
	buffer[2] = bodies[j].Xc.z;
	buffer[3] = bodies[j].v.x;
	buffer[4] = bodies[j].v.y;
	buffer[5] = bodies[j].v.z;

      }

      MPI::COMM_WORLD.Bcast(buffer, 6, MPI::DOUBLE, i);

      if(This_Processor() != i){

	bodies[j].Xc.x = buffer[0];
	bodies[j].Xc.y = buffer[1];
	bodies[j].Xc.z = buffer[2];
	bodies[j].v.x = buffer[3];
	bodies[j].v.y = buffer[4];
	bodies[j].v.z = buffer[5];

      }
    }

  }

}

void Planets_MPI_Broadcast_Positions2(Body *bodies, const int& NCl, const int& NCu){

  int NCl_root, NCu_root, i, j;
  double buffer[3];

  for(i = 0; i < N_Processors(); ++i){

    if(i == This_Processor()){
      NCl_root = NCl;
      NCu_root = NCu;
    }

    MPI::COMM_WORLD.Bcast(&(NCl_root), 1, MPI::INT, i);
    MPI::COMM_WORLD.Bcast(&(NCu_root), 1, MPI::INT, i);

    for(j = NCl_root; j <= NCu_root; ++j){

      if(This_Processor() == i){
	buffer[0] = bodies[j].Xc.x;
	buffer[1] = bodies[j].Xc.y;
	buffer[2] = bodies[j].Xc.z;
      }

      MPI::COMM_WORLD.Bcast(buffer, 3, MPI::DOUBLE, i);

      if(This_Processor() != i){
	bodies[j].Xc.x = buffer[0];
	bodies[j].Xc.y = buffer[1];
	bodies[j].Xc.z = buffer[2];
      }
      
    }

  }

}
