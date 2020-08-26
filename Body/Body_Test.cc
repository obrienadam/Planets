#include "Multi_Body_System.h"
#include "../Planets_MPI/Planets_MPI.h"
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int main  (){

  Planets_MPI_Init();

  Multi_Body_System test;
  Multi_Body_System test2;

  int i = 0;
  double tmax, dt = 1e7;
  tmax = dt*20000000000;

  double light_yr_conversion = 9.4605284e15;

  test2.Generate_Random_Star_Cluster_Cold_Collapse(640, 4*light_yr_conversion, 1.981e30, 0.1*1.981e30);

  Vector3D Upper_L = Vector3D(40*light_yr_conversion, 40*light_yr_conversion, 40*light_yr_conversion);
  Vector3D Lower_L = Vector3D(-40*light_yr_conversion, -40*light_yr_conversion, -40*light_yr_conversion);

  test2.Set_Domain_Limits(Lower_L, Upper_L);

  for(double t = 0; t < tmax; t += dt, ++i){

    test2.Time_Integration_Leap_Frog(dt);

    if(i%1000 == 0 && Main_Processor()){

      cout << "Day: " << t/60/60/24 << endl;
      cout << "Star Mass: " << endl;
      cout << test2.bodies[3].mass << endl;
      cout << "Star Position: " << endl;
      test2.bodies[3].Print_Position_LY();
      cout << "Star Velocity: " << endl;
      test2.bodies[3].v.Print();
      cout << "Star Acceleration: " << endl;
      test2.bodies[3].dv_dt.Print();

      if(i%2000 == 0){
	test2.Output_Particle_Map();
      }
    }
  }

  Planets_MPI_Finalize();

  return 0;
}
