#include "Multi_Body_System.h"
#include "../Planets_MPI/Planets_MPI.h"
#include "Body.h"
#include "../Octree/Octree.h"
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <sstream>
#include <algorithm>

Multi_Body_System::Multi_Body_System(){
  allocated = false;
  restrict_solution_domain = false;
  N = NCl = NCu = 0;
  N_time_steps = 0;

  Time_Integration = LEAP_FROG;
  Force_Evaluation = BINARYTREE_FML;
}

Multi_Body_System::Multi_Body_System(const int& N_init){
  allocate(N_init);
}

Multi_Body_System::~Multi_Body_System(){
  
}

void Multi_Body_System::allocate(const int& N_init){

  if(allocated == true){
    deallocate();
  }

  N = N_init;
  N_time_steps = 0;
  bodies.resize(N);
  bodies_previous_step.resize(N);
  allocated = true;
  restrict_solution_domain = false;

  Time_Integration = LEAP_FROG;
  Force_Evaluation = BINARYTREE_FML;

  Set_Iteration_Bounds();
}

void Multi_Body_System::deallocate(){

    N = NCl = NCu = 0;
    N_time_steps = 0;
    bodies.clear();
    bodies_previous_step.clear();
    allocated = false;
    restrict_solution_domain = false;
}

void Multi_Body_System::Set_Iteration_Bounds(){

  int quotient = N/N_Processors();
  int modulus = N%N_Processors();
  int NCu_local, NCl_local, N_local;
  int i;

  for(NCu_local = N-1, i = N_Processors()-1; i >= 0; --i){

    N_local = quotient;

    if(modulus > 0){
      ++N_local;
      --modulus;
    }

    NCl_local = NCu_local - N_local + 1;

    if(i == This_Processor()){

      NCu = NCu_local;
      NCl = NCl_local;

    }

    NCu_local = NCl_local - 1;

  }

  std::cout << "On process : " << This_Processor() << std::endl
	    << "NCl: " << NCl << std::endl
	    << "NCu: " << NCu << std::endl;


}

void Multi_Body_System::Generate_Earth_and_Ball_Test(){

  allocate(3);
  
  bodies[0].mass = 5.97219e24;
  bodies[1].mass = 7.34767309e22; // baseball mass
  bodies[1].Xc = Vector3D(384403e3, 0, 0);
  bodies[1].v = Vector3D(0, -800, 0);

  bodies[2].mass = 5.34767309e22;
  bodies[2].Xc = Vector3D(-384403e3, 0, 0);
  bodies[2].v = Vector3D(0, 800, 0);
}

void Multi_Body_System::Generate_Solar_System(){

  allocate(10);
  
  double eov = 67000*1.72*1000/3600.0;


  bodies[0].body_name = "Sun";
  bodies[1].body_name = "Mercury";
  bodies[2].body_name = "Venus";
  bodies[3].body_name = "Earth";
  bodies[4].body_name = "Mars";
  bodies[5].body_name = "Jupiter";
  bodies[6].body_name = "Saturn";
  bodies[7].body_name = "Uranus";
  bodies[8].body_name = "Neptune";
  bodies[9].body_name = "Pluto";

  bodies[0].mass = 1.9891e30; // The sun
  bodies[1].mass = 328.5e21; // Mercury
  bodies[2].mass = 4.867e24; // Venus
  bodies[3].mass = 5.97219e24; // Earth
  bodies[4].mass = 639e21; // Mars
  bodies[5].mass = 1.89813e27; // Jupiter
  bodies[6].mass = 568.3e24; // Saturn
  bodies[7].mass = 8.68103e25; // Uranus
  bodies[8].mass = 102.4e24; // Neptune
  bodies[9].mass = 1.30900e22; // Pluto

  
  bodies[0].Xc = Vector3D(0, 0, 0);
  bodies[1].Xc = Vector3D(57.9e9, 0, 0);
  bodies[2].Xc = Vector3D(108.2e9, 0, 0);
  bodies[3].Xc = Vector3D(149.6e9, 0, 0);
  bodies[4].Xc = Vector3D(227.9e9, 0, 0);
  bodies[5].Xc = Vector3D(778.4e9, 0, 0);
  bodies[6].Xc = Vector3D(1426.7e9, 0, 0);
  bodies[7].Xc = Vector3D(2871e9, 0, 0);
  bodies[8].Xc = Vector3D(4498.3e9, 0, 0);
  bodies[9].Xc = Vector3D(5906376272e3, 0, 0);

  bodies[0].v = Vector3D(0, 0, 0);
  bodies[1].v = Vector3D(0, 0, 1.607*eov);
  bodies[2].v = Vector3D(0, 0, 1.174*eov);
  bodies[3].v = Vector3D(0, 0, eov);
  bodies[4].v = Vector3D(0, 0, 0.802*eov);
  bodies[5].v = Vector3D(0, 0, 0.434*eov);
  bodies[6].v = Vector3D(0, 0, 0.323*eov);
  bodies[7].v = Vector3D(0, 0, 0.228*eov);
  bodies[8].v = Vector3D(0, 0, 0.182*eov);
  bodies[9].v = Vector3D(0, 0, 0.159*eov);

  Space_Partition.Set_Limits(Vector3D(-6e12, -6e12, -6e12), Vector3D(6e12, 6e12, 6e12));
  Space_Partition2.Set_Limits(Vector3D(-6e12, -6e12, -6e12), Vector3D(6e12, 6e12, 6e12));

  epsilon_sqr = 0;

}

void Multi_Body_System::Generate_Random_Star_Cluster_Cold_Collapse(const int& N_stars, const double& sphere_radius, const double& mass_avg, const double& mass_std_dev){

  allocate(N_stars);

  srand(time(NULL));

  double P, mass, radius, theta, phi, sign;
  const double PI = 4*atan(1);

  double norm_mass_avg = 1;
  double norm_mass_std_dev = mass_std_dev/mass_avg;

  for(int i = 0; i < N_stars; ++i){

    if(Main_Processor()){

      P = double(rand())/double(RAND_MAX);

      if(P <= 0.5){
	sign = -1.0;
      } else {
	sign = 1.0;
      }

      P = double(rand())/double(RAND_MAX);

      mass = sign*sqrt(log(pow(P*norm_mass_std_dev*sqrt(2.0*PI), -2.0*pow(norm_mass_std_dev, 2)))) + 1.0;
      mass *= mass_avg;

      if(mass < 0){
	std::cout << "Warning: Particle " << i << " was assigned a negative mass." << std::endl;
      }

      P = double(rand())/double(RAND_MAX);

      radius = sphere_radius*P;

      P = double(rand())/double(RAND_MAX);

      theta = 2.0*PI*P;

      P = double(rand())/double(RAND_MAX);

      phi = 2.0*PI*P;
    }

    MPI::COMM_WORLD.Bcast(&(mass), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(radius), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(theta), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(phi), 1, MPI::DOUBLE, 0);

    bodies[i].mass = mass;
    bodies[i].Xc.x = radius*sin(theta)*cos(phi);
    bodies[i].Xc.y = radius*sin(theta)*sin(phi);
    bodies[i].Xc.z = radius*cos(theta);

    bodies[i].v.x = 5000*(radius*radius)/(sphere_radius*sphere_radius)*sin(theta)*cos(phi);
    bodies[i].v.y = 5000*(radius*radius)/(sphere_radius*sphere_radius)*sin(theta)*sin(phi);
    bodies[i].v.z = 5000*(radius*radius)/(sphere_radius*sphere_radius)*cos(theta);
  }

  Space_Partition.Set_Limits(Vector3D(-1.2*sphere_radius, -1.2*sphere_radius, -1.2*sphere_radius), Vector3D(1.2*sphere_radius, 1.2*sphere_radius, 1.2*sphere_radius));
  Space_Partition2.Set_Limits(Vector3D(-1.2*sphere_radius, -1.2*sphere_radius, -1.2*sphere_radius), Vector3D(1.2*sphere_radius, 1.2*sphere_radius, 1.2*sphere_radius));
  
  Set_Domain_Limits(Vector3D(-1.2*sphere_radius, -1.2*sphere_radius, -1.2*sphere_radius), Vector3D(1.2*sphere_radius, 1.2*sphere_radius, 1.2*sphere_radius));

  if(Main_Processor()){
    std::cout << "Allocation of Star Cluster complete... \n";
  }
}

void Multi_Body_System::Generate_Random_Structure_Formation(const int& N, const double& cube_dimension, const double& mass_avg, const double& mass_std_dev){

  allocate(N);

  srand(time(NULL));

  double P, mass, x, y, z, sign;
  const double PI = 4*atan(1);

  double norm_mass_avg = 1;
  double norm_mass_std_dev = mass_std_dev/mass_avg;

  for(int i = 0; i < N; ++i){

    if(Main_Processor()){

      P = double(rand())/double(RAND_MAX);

      if(P <= 0.5){
	sign = -1.0;
      } else {
	sign = 1.0;
      }

      P = double(rand())/double(RAND_MAX);

      mass = sign*sqrt(log(pow(P*norm_mass_std_dev*sqrt(2.0*PI), -2.0*pow(norm_mass_std_dev, 2)))) + 1.0;
      mass *= mass_avg;

      if(mass < 0){
	std::cout << "Warning: Particle " << i << " was assigned a negative mass." << std::endl;
      }

      P = double(rand())/double(RAND_MAX);

      x = P*cube_dimension;

      P = double(rand())/double(RAND_MAX);

      y = P*cube_dimension;

      P = double(rand())/double(RAND_MAX);

      z = P*cube_dimension;
    }

    MPI::COMM_WORLD.Bcast(&(mass), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(x), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(y), 1, MPI::DOUBLE, 0);
    MPI::COMM_WORLD.Bcast(&(z), 1, MPI::DOUBLE, 0);

    bodies[i].mass = mass;
    bodies[i].Xc.x = x;
    bodies[i].Xc.y = y;
    bodies[i].Xc.z = z;

  }

  Set_Domain_Limits(Vector3D(0,0,0), Vector3D(cube_dimension, cube_dimension, cube_dimension));

  Space_Partition.Set_Limits(Vector3D(0,0,0), Vector3D(cube_dimension, cube_dimension, cube_dimension));
  Space_Partition2.Set_Limits(Vector3D(0,0,0), Vector3D(cube_dimension, cube_dimension, cube_dimension));

  if(Main_Processor()){
    std::cout << "Allocation of structure formation complete... \n";
  }
}

void Multi_Body_System::Time_Integration_Euler(const double& dt){

  copy_previous();

  Body::Compute_Softened_NBody_Forces(bodies.data(), N, NCl, NCu, epsilon_sqr);

  for(int i = NCl; i <= NCu; ++i){
  
    bodies[i].v += bodies[i].dv_dt*dt;
    bodies[i].Xc += bodies_previous_step[i].v*dt;

  }

  ++N_time_steps;

  Planets_MPI_Broadcast_Positions(bodies.data(), NCl, NCu);
  Remove_Out_of_Bound_Bodies();

}

void Multi_Body_System::Time_Integration_Leap_Frog(const double& dt){

  if(N_time_steps == 0){

    double h = 0.5*dt; //  we only want a half-time step for the first velocity calculation

    copy_previous();

    Time_Integration_RK4(h); // start the leapfrog integration to get vn+1/2
    
    for(int i = NCl; i <= NCu; ++i)
      {
	bodies[i].Xc = bodies_previous_step[i].Xc + bodies[i].v*dt;
      }

  } else {

    //Body::Compute_Softened_NBody_Forces(bodies.data(), N, NCl, NCu, epsilon_sqr);

    Space_Partition2.Clear();

    for(int i = 0; i < N; ++i)
      Space_Partition2.Insert(bodies[i]);

    for(int i = NCl; i <= NCu; ++i)
      {
	Space_Partition2.Compute_Softened_Force(bodies[i], 0.5, epsilon_sqr);
	bodies[i].v += bodies[i].F*dt/bodies[i].mass;
	bodies[i].Xc += bodies[i].v*dt;
      }
  }

  ++N_time_steps;

  Planets_MPI_Broadcast_Positions(bodies.data(), NCl, NCu);
   
  Remove_Out_of_Bound_Bodies();
  
}

void Multi_Body_System::Simulate(const Input& Multi_Body_IPs){

  double t, dt = Multi_Body_IPs.Dt, tmax = Multi_Body_IPs.Time_Max;
  int k = 0, Output_Freq = int(tmax/dt/100.0);

  epsilon_sqr = Multi_Body_IPs.Epsilon*Multi_Body_IPs.Epsilon;
  Time_Integration = Multi_Body_IPs.Time_Integration;

  if(Multi_Body_IPs.IC_Type == STAR_COLD_COLLAPSE){

    Generate_Random_Star_Cluster_Cold_Collapse(Multi_Body_IPs.N_Particles, Multi_Body_IPs.Sphere_Radius, Multi_Body_IPs.Mass_Average, Multi_Body_IPs.Mass_Std_Dev);

  } else if (Multi_Body_IPs.IC_Type == SOLAR_SYSTEM){

    Generate_Solar_System();

  }

  if(Main_Processor()){
    std::cout << "[";
    std::cout.flush();
  }

  for(t = 0; t < tmax; t += dt){

    switch (Time_Integration){

    case LEAP_FROG:

      Time_Integration_Leap_Frog(dt);

      break;

    case RUNGE_KUTTA_4:

      Time_Integration_RK4(dt);

      break;

    case RUNGE_KUTTA_2:

      Time_Integration_RK2(dt);

      break;

    case EULER:

      Time_Integration_Euler(dt);

      break;

    default:

      Time_Integration_Leap_Frog(dt);
    };

    ++k;

    if(Main_Processor()){
      if(k%Output_Freq == 0){
	std::cout << "|";
	std::cout.flush();
      }
    }
  }

  if(Main_Processor()){
    std::cout << "]" << std::endl
	      << "Simulation Complete." << std::endl;
  }
}

void Multi_Body_System::copy_previous(){

  for(int i = NCl; i <= NCu; ++i){
    bodies_previous_step[i] = bodies[i];
  }

}

void Multi_Body_System::Time_Integration_RK2(const double& dt){

}

void Multi_Body_System::Time_Integration_RK4(const double& dt){

  int i;

  Vector3D *k1 = new Vector3D[N];
  Vector3D *k2 = new Vector3D[N];
  Vector3D *k3 = new Vector3D[N];
  Vector3D *k4 = new Vector3D[N];

  copy_previous();
  
  Body::Compute_Softened_NBody_Forces(bodies.data(), N, NCl, NCu, epsilon_sqr); // basically f(t, x)

  for(i = NCl; i <= NCu; ++i){
    k1[i] = dt*bodies[i].dv_dt;
    bodies[i].Xc = bodies_previous_step[i].Xc + 0.5*k1[i];
  }

  Body::Compute_Softened_NBody_Forces(bodies.data(), N, NCl, NCu, epsilon_sqr);

  for(i = NCl; i <= NCu; ++i){
    k2[i] = dt*bodies[i].dv_dt;
    bodies[i].Xc = bodies_previous_step[i].Xc + 0.5*k2[i];
  }

  Body::Compute_Softened_NBody_Forces(bodies.data(), N, NCl, NCu, epsilon_sqr);

  for(i = NCl; i <= NCu; ++i){
    k3[i] = dt*bodies[i].dv_dt;
    bodies[i].Xc = bodies_previous_step[i].Xc + k3[i];
  }

  Body::Compute_Softened_NBody_Forces(bodies.data(), N, NCl, NCu, epsilon_sqr);

  for(i = NCl; i <= NCu; ++i){
    k4[i] = dt*bodies[i].dv_dt;
  }

  for(i = NCl; i <= NCu; ++i){
    bodies[i].v += (1.0/6.0)*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
    bodies[i].Xc = bodies_previous_step[i].Xc + dt*bodies[i].v;
  }

  delete[] k1;
  delete[] k2;
  delete[] k3;
  delete[] k4;

  ++N_time_steps;

  Planets_MPI_Broadcast_Positions(bodies.data(), NCl, NCu);
  
  Remove_Out_of_Bound_Bodies();
}

void Multi_Body_System::Output_File(){

  for(int i = NCl; i <= NCu; ++i){
    bodies[i].Output_File();
  }

}

void Multi_Body_System::Set_Domain_Limits(const Vector3D& lower_limits_init, const Vector3D& upper_limits_init){

  domain_lower_limits = lower_limits_init;
  domain_upper_limits = upper_limits_init;
  restrict_solution_domain = true;

}

void Multi_Body_System::Remove_Out_of_Bound_Bodies(){

  bool particles_need_removal = false;
  int Nnew = N;

  std::vector<int> remove_index;

  for(int i = 0; i < N; ++i){

    if(bodies[i].Xc.x > domain_upper_limits.x || bodies[i].Xc.y > domain_upper_limits.y || bodies[i].Xc.z > domain_upper_limits.z ||
       bodies[i].Xc.x < domain_lower_limits.x || bodies[i].Xc.y < domain_lower_limits.y || bodies[i].Xc.z < domain_lower_limits.z){
      
      remove_index.push_back(i);

      particles_need_removal = true;
    }
  }       

  if(particles_need_removal == true){

    std::vector<Body> temp = bodies;
    std::vector<Body> temp2 = bodies_previous_step;
    bodies.clear();
    bodies_previous_step.clear();

    bool keep_particle = true;

    int Ni = remove_index.size();

    for(int i = 0; i < N; ++i){

      for(int j = 0; j < Ni; ++j){
	
	if(i == remove_index[j]){

	  --Nnew;

	  if(Main_Processor()){
	    std::cout << "Removing particle " << i << "." << std::endl
		      << "Particle Location: " << bodies[i].Xc.x << ", " << bodies[i].Xc.y << ", " << bodies[i].Xc.z << std::endl
		      << "Particles Remaining: " << Nnew << std::endl;
	  }

	  keep_particle = false;
	  break;
	}
	
      }

      if(keep_particle == true){
	bodies.push_back(temp[i]);
	bodies_previous_step.push_back(temp2[i]);
      }
      
      keep_particle = true;
    }

    N = bodies.size();

    MPI::COMM_WORLD.Barrier();
    
    Set_Iteration_Bounds();
  }
}

void Multi_Body_System::Output_Particle_Map(){

  if(Main_Processor()){

    std::ofstream fout;

    fout.open("Particle_Map.csv");

    for(int i = 0; i < N; ++i){

      fout << bodies[i].Xc.x/9.4605284e15 << ", " << bodies[i].Xc.y/9.4605284e15 << ", " << bodies[i].Xc.z/9.4605284e15 << std::endl;

    }
  
    fout.close();
  }
}

void Multi_Body_System::Write_Particles_to_Tec360(const double& t){

  if(Main_Processor()){

    std::stringstream filename;
  std::ofstream fout;

  filename << "Particles_t_" << t << ".dat";

  fout.open(filename.str().c_str());

  fout << "TITLE = \"Planets Particles\"" << std::endl
       << "VARIABLES = \"X\", \"Y\", \"Z\", \"M\"" << std::endl
       << "ZONE T = \"Particles\", I=" << N << ", J=1, K=1, F=POINT, STRANDID=1, SOLUTIONTIME=" << t << std::endl;

  for(int i = 0; i < N; ++i){

    fout << bodies[i].Xc.x << " " << bodies[i].Xc.y << " " << bodies[i].Xc.z << " " << bodies[i].mass << std::endl;

  }

  fout.close();

  }

}

double Multi_Body_System::Determine_dt(const double& max_integration_distance){

  double dtx, dty, dtz, dt;

  for(int i = 0; i < N; ++i){
    dtx = max_integration_distance/bodies[i].v.x;
    dty = max_integration_distance/bodies[i].v.y;
    dtz = max_integration_distance/bodies[i].v.z;
  }

  dt = std::min(dtx, std::min(dty, dtz));

  MPI::COMM_WORLD.Allreduce(&dt, &dt, 1, MPI::DOUBLE, MPI::MIN);

  return fabs(dt);
}
