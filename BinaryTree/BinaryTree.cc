#include "BinaryTree.h"
#include <sstream>
#include <fstream>
#include "../Geometry/MathVector.h"
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <algorithm>

BinaryTree::BinaryTree(const Vector3D& DLL, const Vector3D& DUL){

  Set_Limits(DLL, DUL);

  parent = NULL;
  Left = NULL;
  Right = NULL;

  split_direction = NONE;
  Level = 0;

  contains_children = contains_particle = false;

  particle_p = NULL;
}

void BinaryTree::Set_Limits(const Vector3D& DLL, const Vector3D& DUL){

  HLL = DLL;
  HUL = DUL;

  MP = 0.5*(HLL + HUL);
  
  Delta = HUL - HLL;

  s = std::max(Delta.x, std::max(Delta.y, Delta.z));
}

int BinaryTree::Determine_Split_Direction(const Vector3D& point1, const Vector3D& point2){

  double SPx = Delta.x, SPy = Delta.y, SPz = Delta.z;

  if(point1.x > MP.x && MP.x > point2.x || 
     point1.x < MP.x && MP.x < point2.x){

    SPx += SPx;

  }

  if(point1.y > MP.y && MP.y > point2.y || 
     point1.y < MP.y && MP.y < point2.y){

    SPy += SPy;

  }

  if(point1.z > MP.z && MP.z > point2.z || 
     point1.z < MP.z && MP.z < point2.z){

    SPz += SPz;

  }

  if(SPx > SPy && SPx > SPz){
    return EW;
  } else if (SPx > SPy) {
    return TB;
  } else {
    return NS;
  }
     
}

int BinaryTree::Determine_Half(const Vector3D& point){

  if(split_direction == EW){

    if(point.x > MP.x){
      return RIGHT;
    } else {
      return LEFT;
    }

  } else if (split_direction == NS){

    if(point.y > MP.y){
      return RIGHT;
    } else {
      return LEFT;
    }
    
  } else {

    if(point.z > MP.z){
      return RIGHT;
    } else {
      return LEFT;
    }

  }
}

void BinaryTree::Split(const int& direction){

  parent = this;
  Left = new BinaryTree;
  Right = new BinaryTree;

  contains_children = true;

  split_direction = direction;

  if(direction == EW){

    Left->Set_Limits(HLL, Vector3D(0.5*(HLL.x + HUL.x), HUL.y, HUL.z));
    Right->Set_Limits(Vector3D(0.5*(HLL.x + HUL.x), HLL.y, HLL.z), HUL);

  } else if (direction == NS){

    Left->Set_Limits(HLL, Vector3D(HUL.x, 0.5*(HLL.y + HUL.y), HUL.z));
    Right->Set_Limits(Vector3D(HLL.x, 0.5*(HLL.y + HUL.y), HLL.z), HUL);

  } else {

    Left->Set_Limits(HLL, Vector3D(HUL.x, HUL.y, 0.5*(HLL.z + HUL.z)));
    Right->Set_Limits(Vector3D(HLL.x, HLL.y, 0.5*(HLL.z + HUL.z)), HUL);

  }

  Left->Level = Level+1;
  Right->Level = Level+1;

}

void BinaryTree::Clear(){

  if(contains_children){

    Left->Clear();
    Right->Clear();

    delete Left;
    delete Right;

    contains_children = false;
    contains_particle = false;
    mass = 0;
    Mass_Center = Vector3D(0,0,0);

  }
}

void BinaryTree::Partition_Test(){

  Split(EW);

  Left->Split(NS);
  Right->Split(TB);
  
  Left->Left->Split(NS);
  Right->Right->Split(EW);

  Write_to_Tec360(0);
}

void BinaryTree::Write_to_Tec360(const double& t){

  if(Level == 0){

    std::stringstream filename;
    std::ofstream fout;

    filename << "Binary_Tree_Mesh_t_" << t << ".dat";

    fout.open(filename.str().c_str());

    fout << "TITLE = \"Planets Binary Tree\"" << std::endl
         << "VARIABLES = \"X\", \"Y\", \"Z\", \"P\"" << std::endl;

    fout.close();
  }

  if(contains_children){

    Left->Write_to_Tec360(t);
    Right->Write_to_Tec360(t);

  } else {

    std::stringstream filename;
    std::ofstream fout;

    filename << "Binary_Tree_Mesh_t_" << t << ".dat";

    fout.open(filename.str().c_str(), std::ofstream::app);

    int physical_particle;

    if(particle_p != NULL){
	physical_particle = 1;
      }
    else{
        physical_particle = 0;
      }

    fout << "ZONE T = \"Level" << Level << "\", I=2, J=2, K=2, F=POINT, STRANDID=2, SOLUTIONTIME=" << t << std::endl
         << HLL.x << " " << HLL.y << " " << HLL.z << " " << physical_particle << std::endl
         << HUL.x << " " << HLL.y << " " << HLL.z << " " << physical_particle << std::endl
         << HLL.x << " " << HUL.y << " " << HLL.z << " " << physical_particle << std::endl
         << HUL.x << " " << HUL.y << " " << HLL.z << " " << physical_particle << std::endl
         << HLL.x << " " << HLL.y << " " << HUL.z << " " << physical_particle << std::endl
         << HUL.x << " " << HLL.y << " " << HUL.z << " " << physical_particle << std::endl
         << HLL.x << " " << HUL.y << " " << HUL.z << " " << physical_particle << std::endl
         << HUL.x << " " << HUL.y << " " << HUL.z << " " << physical_particle << std::endl;

    fout.close();

  }
}

void BinaryTree::Insert(Body& particle){

  if(particle.Xc.x > HUL.x || particle.Xc.z > HUL.z || particle.Xc.z > HUL.z ||
     particle.Xc.x < HLL.x || particle.Xc.y < HLL.y || particle.Xc.z < HLL.z){
    
    std::cout << "Error, particle is out of bounds." << std::endl;
    return;
  }

  if(contains_particle == false){ // there is room

    particle_p = &particle;

    mass = particle.mass;
    Mass_Center = particle.Xc;

    contains_particle = true;

  } else { // particle(s) is present

    int half;
    double old_mass;

    old_mass = mass;
    mass += particle.mass;
    Mass_Center *= (old_mass/mass);
    Mass_Center += (particle.mass/mass)*particle.Xc;
    

    if(contains_children){ // Node is already split
      
      half = Determine_Half(particle.Xc);
      
      if (half == LEFT){ 
	Left->Insert(particle);
      } else {
	Right->Insert(particle);
      }

    } else { // Node is not split, must be split

      Split(Determine_Split_Direction((*particle_p).Xc, particle.Xc));
   
      // Move the original particle pointer

      half = Determine_Half((*particle_p).Xc);

      if (half == LEFT){
	Left->Insert(*particle_p);
      } else {
	Right->Insert(*particle_p);
      }

      
      
      particle_p = NULL;

      // Move the current particle

      half = Determine_Half(particle.Xc);

      if (half == LEFT){
        Left->Insert(particle);
      } else {
        Right->Insert(particle);
      }

    }
  }
}

void BinaryTree::Compute_Softened_Force(Body& particle, const double& theta_MAX, const double& epsilon_sqr){

  if(contains_particle == false){
    return;
  } else if (Level == 0){
    particle.F = Vector3D(0, 0, 0);
  }
  
  if(contains_children == false){

    if(&particle == particle_p)
      return;

    double r_2 = Vector3D::Distance_Between_Squared(particle.Xc, Mass_Center);

    particle.F += (GRAVITATION_CONSTANT*particle.mass*mass/pow(epsilon_sqr + r_2, 1.5))*(Mass_Center - particle.Xc);

    return;
  }

  double theta = s/Vector3D::Distance_Between(particle.Xc, Mass_Center);

  if(theta < theta_MAX){

    double r_2 = Vector3D::Distance_Between_Squared(particle.Xc, Mass_Center);

    
    particle.F += (GRAVITATION_CONSTANT*particle.mass*mass/pow(epsilon_sqr + r_2, 1.5))*(Mass_Center - particle.Xc);

  } else {

    Left->Compute_Softened_Force(particle, theta_MAX, epsilon_sqr);
    Right->Compute_Softened_Force(particle, theta_MAX, epsilon_sqr);

  }
}  
