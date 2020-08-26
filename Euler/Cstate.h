#pragma once

#include "../Geometry/MathVector.h"
#include "Pstate.h"

class Pstate;

class Cstate{

 private:
 public:

  Cstate();
  Cstate(const Pstate& W_init);

  double rho, e;
  Vector3D mom;

  Pstate W();

};

class Cflux{

 private:
 public:

  Cflux();
  Cflux(const double& e_init, const Vector3D& mom_init, const double& rho_init);
  Cflux(const Cstate& q, const Vector3D& unit_normal);

  double rho, e;
  Vector3D mom;

  Pstate W();
  Cstate U();

};
