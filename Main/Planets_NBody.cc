#include <iostream>
#include "../Input/Input.h"
#include "../Body/Multi_Body_System.h"
#include "../Planets_MPI/Planets_MPI.h"
#include <cstring>

using namespace std;

int main(int argc, char *argv[]){
  
  Input NBody_IPs;
  Multi_Body_System Particles;

  Planets_MPI_Init();

  if( strcmp(argv[1], "-h") == 0){

    cout << "Help not yet available." << endl;

  } else if (strcmp(argv[1], "-f") == 0 && argc >= 2) {

    NBody_IPs.Read_Input_File(argv[2]);
    
    Particles.Simulate(NBody_IPs);

  }

  Planets_MPI_Finalize();

  return 0;
}
