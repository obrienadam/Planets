#include "Cell.h"
#include "../Geometry/MathVector.h"

Cell::Cell(const Vector3D& Delta_init, const Vector3D& Center_init){

  Delta = Delta_init;
  Center = Center_init;

  Set_Areas();
}

void Cell::Set_Areas(){
  Area_T = Delta.x*Delta.y;
  Area_B = Delta.x*Delta.y;
  Area_N = Delta.x*Delta.z;
  Area_S = Delta.x*Delta.z;
  Area_E = Delta.y*Delta.z;
  Area_W = Delta.y*Delta.z;
}

