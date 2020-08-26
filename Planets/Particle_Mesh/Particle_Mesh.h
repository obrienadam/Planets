#ifndef _PARTICLE_MESH_H_
#define _PARTICLE_MESH_H_

#include "../Geometry/MathVector.h"
#include "../Body/Body.h"
#include <complex>

class Particle_Mesh{

 private:
  
  int NNi, NNj, NNk;
 
  bool allocated;
  void allocate(const int& Nx, const int& Ny, const int& Nz, const Vector3D& lower_bounds, const Vector3D& upper_bounds);
  void deallocate();

 public:

  Particle_Mesh();
  ~Particle_Mesh();
  
  double ***mass;
  std::complex<double> ***rho_ft;
  Vector3D ***Xc;

  Vector3D Delta;
  double Vol;

  void init(const int& Nx, const int& Ny, const int& Nz, const Vector3D& lower_bounds, const Vector3D& upper_bounds);

  void Reset_Mass();

  void Find_Nearest_Node(int& i, int& j, int& k, const Vector3D& point);
  void Find_Nearest_Node2(int& i, int& j, int& k, const Vector3D& point);
  void Restrict_Particles_to_Grid(const Body *particles, const int& N);
  void Compute_DFT();
  };

#endif
