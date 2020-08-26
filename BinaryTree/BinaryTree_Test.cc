#include "BinaryTree.h"
#include "../Geometry/MathVector.h"
#include "../Body/Multi_Body_System.h"
#include "../Planets_MPI/Planets_MPI.h"
#include <time.h>

using namespace std;

int main(){

  Planets_MPI_Init();

  double radius_LY = 9.4605284e15*1.0; 
  double LY_conversion = 9.4605284e15;

  BinaryTree data_tree(Vector3D(-radius_LY, -radius_LY, -radius_LY), Vector3D(radius_LY, radius_LY, radius_LY));
  Multi_Body_System particles;
  double epsilon_sqr = 100;
  int Np = 1500;
  double t, tmax, dt = 1e8;
  int k = 0, output_freq = 1e6;

  clock_t sim_time;

  tmax = 600000.0*365.0*24.0*3600.0;

  particles.Generate_Random_Star_Cluster_Cold_Collapse(Np, radius_LY, 1.981e30, 1.981*0.1e30);

  sim_time = clock();

  for(t = 0; t <= tmax; t += dt){
    
    particles.Time_Integration_Leap_Frog(dt);

    if(k%500 == 0)
      {
	if(Main_Processor()){
	  

	std::cout << "Time: " << t << std::endl
		  << "Dt: " << dt << std::endl
		  << "Simulation time: " << (clock() - sim_time)/CLOCKS_PER_SEC << std::endl
		  << "Star 1 Position: " << "(" << particles.bodies[0].Xc.x/LY_conversion << ", " << particles.bodies[0].Xc.y/LY_conversion << ", " << particles.bodies[0].Xc.z/LY_conversion << ")" << std::endl
		  << "Star 1 Velocity: " << "(" << particles.bodies[0].v.x << ", " << particles.bodies[0].v.y << ", " << particles.bodies[0].v.z << ")" << std::endl;
	
	particles.Space_Partition2.Write_to_Tec360(t);
	particles.Write_Particles_to_Tec360(t);
	}
      }

    ++k;
  }

  std::cout << "Dt: " << dt << std::endl
	    << "k: " << k << std::endl;

  // data_tree.Write_to_Tec360(0);
  //particles.Write_Particles_to_Tec360(0);

  data_tree.Clear();

  Planets_MPI_Finalize();

  return 0;
}
