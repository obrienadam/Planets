#include "HashGrid.h"
#include "../Geometry/MathVector.h"

void Hashgrid::Deallocate(){

  if(allocated == false)
    return;

  int i, j;

  for(i = 0; i < Ni; ++i){
    for(j = 0; j < Nj; ++j){
      delete[] Cell[i][j];
      delete[] Cell_Contains_Particle[i][j];
    }
    delete[] Cell[i];
    delete[] Cell_Contains_Particle[i];
  }

  delete[] Cell;
  delete[] Cell_Contains_Particle;

  allocated = false;
}

void HashGrid::Init(const int& Ni_init, const int& Nj_init, const int& Nk_init, const Vector3D& Center_init, const Vector3D Delta_init){

  int i, j, k;
  Vector3D LB, UB, Cell_Delta;

  Ni = Ni_init;
  Nj = Nj_init;
  Nk = Nk_init;
  N = Ni*Nj*Nk;

  Center = Center_init;
  Delta = Delta_init;
  Cell_Delta.x = Delta.x/double(Ni);
  Cell_Delta.y = Delta.y/double(Nj);
  Cell_Delta.z = Delta.z/double(Nk);

  LB = -0.5*Delta + Center_init;
  UB = 0.5*Delta + Center_init;

  Cell = new HashGrid_Cell**[Ni];
  Cell_Contains_Particle = new Cell_Index**[Ni];

  for(i = 0; i < Ni; ++i){

    Cell[i] = new HashGrid_Cell*[Nj];
    Cell_Contains_Particle[i] = new Cell_Index*[Nj];

    for(j = 0; j < Nj; ++j){

      Cell[i][j] = new HashGrid_Cell[Nk];
      Cell_Contains_Particle[i][j] = new Cell_Index[Nk];

      for(k = 0; i < Nk; ++k){

	Cell[i][j][k].Center.x = (double(i) + 0.5)*Cell_Delta.x + LB.x;
	Cell[i][j][k].Center.y = (double(j) + 0.5)*Cell_Delta.y + LB.y;
	Cell[i][j][k].Center.z = (double(k) + 0.5)*Cell_Delta.z + LB.z;

	Cell[i][j][k].Delta = Cell_Delta;

	Cell_Contains_Particle[i][j][k].Set_Index(0,0,0);
      }
    }
  }

  allocated = true;

}
