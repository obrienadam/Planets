#pragma once

#include <string>

enum{STAR_COLD_COLLAPSE, DISK_GALAXY, SOLAR_SYSTEM, BINARYTREE_FML, OCTREE_FML, DIRECT_SUM, EULER, RUNGE_KUTTA_4, RUNGE_KUTTA_2, LEAP_FROG};

class Input{

 private:

  // Helper function to remove whitespace

  void Remove_White_Space_and_Comments(std::string& buffer);

 public:

  Input();

  int N_Particles, IC_Type, Force_Evaluation, Time_Integration;

  double Dt, Time_Max;
  double Mass_Average, Mass_Std_Dev, Epsilon, Domain_Width, Domain_Height, Domain_Length, Sphere_Radius;

  double Output_Length_Scale_Modifier, Output_Time_Scale_Modifier, Output_Mass_Scale_Modifier; // calculations are done in meters and seconds, this scales the output

  void Read_Input_File(const char* filename);
};
