#pragma once

#include "Block.h"
#include "Flux_Functions.h"

class Solver{

 private:
 public:

  Solver();

  Block Mesh;
  void Set_ICs_Uniform(const Pstate& state_init);
  void Set_ICs_Riemann(const Pstate& state_L, const Pstate& state_R);
};
