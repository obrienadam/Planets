#pragma once

#include "../Geometry/MathVector.h"

class HashGrid_Cell{

 public:
  
  Vector3D Center;
  Vector3D Delta;

  int N_particles;
  int* particle_index;

};

class Cell_Index{

 public:

  void Set_Index(const int& I_init, const int& J_init, const int& K_init){
    
    I = I_init;
    J = J_init;
    K = K_init;

  }

  int I, J, K;

};

class HashGrid{

 private:

  int Ni, Nj, Nk, N;
  Vector3D Center;
  Vector3D Delta;

  bool allocated;

  void Deallocate();

 public:

  HashGrid(){
    Ni = Nj = Nk = N = 0;
    allocated = false;
  }
  ~HashGrid(){Deallocate();}

  void Init(const int& Ni_init, const int& Nj_init, const int& Nk_init, const Vector3D& Center_init, const Vector3D Delta_init);

  void Sort_SPH_Particles(const BodySPH *bodies, const int& N);

  HashGrid_Cell ***Cell;
  Cell_Index ***Cell_Contains_Particle;
};
