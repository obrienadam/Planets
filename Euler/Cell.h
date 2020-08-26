#pragma once

#include "../Geometry/MathVector.h"
#include "Pstate.h"
#include "Cstate.h"

class Cell{

 private:
  
 public:

  Cell(const Vector3D& Delta_init = Vector3D(0, 0, 0), const Vector3D& Center_init = Vector3D(0, 0, 0));

  Vector3D Delta;
  Vector3D Center;
  void Set_Areas();
  double Area_T, Area_B, Area_N, Area_S, Area_E, Area_W;

  Cstate U;
  Pstate W;

};
