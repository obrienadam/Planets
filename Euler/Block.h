#pragma once

#include "../Geometry/MathVector.h"
#include "Cell.h"
#include <string>

class Block{

 private:
  
  bool Allocated;

  void Deallocate();
  void Init_Geometry(const Vector3D& Delta_init, const Vector3D& Center_init);

 public:
  Block();
  ~Block();

  std::string BLOCK_ID;

  void Allocate(const int& Ni, const int& Nj, const int& Nk, const int& Nghost_init);

  int ICl, ICu, JCl, JCu, KCl, KCu, Nghost, NCi, NCj, NCk;
  Vector3D Delta, Center;

  Cell*** cells;


  Write_to_Tec360(const double& time_strand_ID, const double& time);
};
