#pragma once

#include "Cstate.h"
#include "../Geometry/MathVector.h"

class Cstate;

class Pstate{

 private:
 public:

  Pstate();
  Pstate(const Cstate& U_init);
  Pstate(const double& rho_init, const Vector3D& v_init, const double& P_init); 

  double rho, htot, a, P, T;
  Vector3D v;

  Cstate U();
};
