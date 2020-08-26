#include "Block.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>

Block::Block(){

  cells = NULL;

}

Block::~Block(){

  Deallocate();

}

void Block::Allocate(const int& Ni, const int& Nj, const int& Nk, const int& Nghost_init){

  if(Allocated){
    Deallocate();
  }

  int i, j;

  Nghost = Nghost_init;

  ICl = Nghost;
  ICu = ICl + Ni - 1;
  NCi = ICu + Nghost + 1;

  JCl = Nghost;
  JCu = JCl + Nj - 1;
  NCj = JCu + Nghost + 1;

  KCl = Nghost;
  KCu = KCl + Nk - 1;
  NCk = KCu + Nghost + 1;

  cells = new Cell**[NCi];

  for(i = 0; i < NCi; ++i){

    cells[i] = new Cell*[NCj];

    for(j = 0; j < NCj; ++j){

      cells[i][j] = new Cell[NCk];

    }
  }

  Allocated = true;
}

void Block::Deallocate(){

  if(cells == NULL)
    return;

  int i, j, k;

  for(i = 0; i < NCi; ++i){
    for(j = 0; j < NCj; ++j){
      delete[] cells[i][j];
    }
    delete[] cells[i];
  }

  delete[] cells;

  cells = NULL;
  Allocated = false;
}

void Block::Init_Geometry(const Vector3D& Delta_init, const Vector3D& Center_init){

  if(!Allocated){
    std::cout << "Error, block geometry construction attempted before initialization." << std::endl;
    return;
  }

  int i, j, k;

  Delta = Delta_init;
  Center = Center_init;

  Vector3D Cell_Delta(Delta.x/double((ICu+1-ICl)), Delta.y/double((JCu+1-JCl)), Delta.z/double((KCu+1-KCl)));

  for(i = 0; i < NCi; ++i){
    for(j = 0; j < NCj; ++j){
      for(k = 0; k < NCk; ++k){

	cells[i][j][k].Center = Vector3D(Cell_Delta.x*(0.5 + double(i)), Cell_Delta.y*(0.5 + double(j)), Cell_Delta.z*(0.5 + double(k)));
	cells[i][j][k].Delta = Cell_Delta;
	cells[i][j][k].Set_Areas();
      }
    }
  }
}

void Block::Write_to_Tec360(const double& time_strand_ID, const double& time){

  std::stringstream filename;
  std::ofstream fout;
  int i, j, k;

  filename << BLOCK_ID << "_t_" << time << ".dat";

  fout.open(filename.str().c_str());

  fout << "TITLE = \"" << BLOCK_ID << " " << time << "s\"" << std::endl
       << "VARIABLES = \"X\", \"Y\", \"Z\", \"u\", \"v\", \"w\", \"rho\", \"P\"" << std::endl
       << "ZONE T = \"Block\", I=" << Mesh.NCi - 2*Mesh.Nghost << ", J=" << Mesh.NCj - 2*Mesh.Nghost << ", K=" << Mesh.NCk - 2*Mesh.Nghost << ", F=POINT, STRANDID=" << time_strand_ID << ", SOLUTIONTIME=" << time << std::endl;

  for(i = Mesh.ICl; i <= Mesh.ICu; ++i){
    for(j = Mesh.JCl; j <= Mesh.JCu; ++j){
      for(k = Mesh.KCl; k <= Mesh.KCu; ++k){
	
	fout << Mesh.cells[i][j][k].Center.x << " " << Mesh.cells[i][j][k].Center.y << " " << Mesh.cells[i][j][k].Center.z << " " <<  Mesh.cells[i][j][k].W.v.x	<< " " <<  Mesh.cells[i][j][k].W.v.y << " " <<  Mesh.cells[i][j][k].W.v.z << " " <<  Mesh.cells[i][j][k].W.rho << " " <<  Mesh.cells[i][j][k].W.P << std::endl; 
	
      }
    }
  }

  fout.close();
}
