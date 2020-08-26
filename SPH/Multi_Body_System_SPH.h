#pragma once

#include <vector>
#include "Multi_Body_System_SPH.h"
#include "../BinaryTree/BinaryTree.h"

class Multi_Body_System_SPH{

 public:

  Multi_Body_System_SPH();

  int N, NCl, NCu;

  BinaryTree Gravity_FMM_Tree;
  HashGrid Neighbour_Grid;

  std::vector<BodySPH> bodies;
  std::vector<BodySPH> bodies_previous_step;
  std::vector<BodySPH> boundaries;

  void Init_Liquid_Cube(const int& N_Particles);
};
