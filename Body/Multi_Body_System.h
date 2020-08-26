#pragma once

#include "../Input/Input.h"
#include "Body.h"
#include "../Geometry/MathVector.h"
#include "../Octree/Octree.h"
#include "../BinaryTree/BinaryTree.h"
#include <vector>
#include <fstream>

class Multi_Body_System{

 private:

  bool allocated, restrict_solution_domain;
  int N_time_steps;

  int N, NCl, NCu;

  void allocate(const int& N_init);
  void deallocate();
  void copy_previous();
  void Set_Iteration_Bounds();

  Vector3D domain_lower_limits, domain_upper_limits;

  int Time_Integration;
  int Force_Evaluation;
  double epsilon_sqr;

 public:
  Multi_Body_System();
  ~Multi_Body_System();
  Multi_Body_System(const int& N_init);

  std::vector<Body> bodies_previous_step;
  std::vector<Body> bodies;

  Octree Space_Partition;
  BinaryTree Space_Partition2;

  void Set_Domain_Limits(const Vector3D& lower_limits_init, const Vector3D& upper_limits_init);
  void Remove_Out_of_Bound_Bodies();

  void Generate_Earth_and_Ball_Test();
  void Generate_Solar_System();
  void Generate_Random_Star_Cluster_Cold_Collapse(const int& N_stars, const double& sphere_radius, const double& mass_avg, const double& mass_std_dev);
  void Generate_Random_Structure_Formation(const int& N, const double& cube_dimension, const double& mass_avg, const double& mass_std_dev);

  void Time_Integration_Euler(const double& dt);
  void Time_Integration_RK2(const double& dt);
  void Time_Integration_RK4(const double& dt);
  void Time_Integration_Leap_Frog(const double& dt);

  void Simulate(const Input& Multi_Body_IPs);

  double Determine_dt(const double& max_integration_distance);

  void Output_File();
  void Output_Particle_Map();
  void Write_Particles_to_Tec360(const double& t);
};

