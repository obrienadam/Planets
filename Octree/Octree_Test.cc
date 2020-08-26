#include <iostream>
#include <cstdlib>
#include "Octree.h"
#include "../Geometry/MathVector.h"
#include "../Body/Multi_Body_System.h"
#include "../Planets_MPI/Planets_MPI.h"


using namespace std;

int main(){
  
  Multi_Body_System test1, test2;
  double t, dt = 1e2, tmax = 365*24*60*60*3.0;
  int output_freq = int(tmax/100);

  Planets_MPI_Init();

  int kmax = 20000, k = 0;

  //  test.Generate_Random_Star_Cluster_Cold_Collapse(200, 90, 499, 4);

  test1.Generate_Solar_System();
  test2.Generate_Solar_System();

  cout << "Executing simulation..." << endl;

  for(t = 0; t <= tmax; t += dt){

    test1.Time_Integration_Leap_Frog(dt);

    if(k%10000 == 0){

      cout << "Making frame..." << endl;

      test1.Space_Partition.Write_to_Tec360(t);
      test1.Write_Particles_to_Tec360(t);

    }
    ++k;
  }

  for(int i = 0; i < 10; ++i){
    test1.bodies[i].Xc.Print();
  }

  Planets_MPI_Finalize();

  return 0;
}
