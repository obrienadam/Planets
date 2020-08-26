#include "Particle_Mesh.h"
#include "../Body/Body.h"
#include <iostream>

using namespace std;

int main(){
  Multi_Body_System system;
  Particle_Mesh test_mesh;

  system.Generate_Solar_System();
  test_mesh.init(15, 2, 15, Vector3D(-6e12, -6e12, -6e12), Vector3D(6e12, 6e12, 6e12));

  test_mesh.Restrict_Particles_to_Grid(system.bodies, system.N);
  test_mesh.Compute_DFT();

  return 0;
}
