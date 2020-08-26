#include "Solver.h"

Solver::Solver(){

  Mesh.Allocate(100, 100, 100, 2);
  Mesh.Init_Geometry(Vector3D(10, 10, 10), Vector3D(0, 0, 0));

}

void Solver::Set_ICs_Uniform(const Pstate& state_init){

  

}

void Solver::Set_ICs_Riemann(const Pstate& state_L, const Pstate& state_R){

  for(int i = Mesh.ICl; i <= Mesh.ICu; ++i){
    for(int j = Mesh.JCl; i <= Mesh.JCu; ++j){
      for(int k = Mesh.KCl; k <= Mesh.KCu; ++k){

	if(Mesh.cells[i][j][k].Center.x < Mesh.Center.x &&
	   Mesh.cells[i][j][k].Center.y < Mesh.Center.y &&
	   Mesh.cells[i][j][k].Center.z < Mesh.Center.z){

	  Mesh.cells[i][j][k].W = state_L;

	} else {
	  
	  Mesh.cells[i][j][k].W = state_R;

	}

	}

      }
    }
  }
}
