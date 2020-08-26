#pragma once

#include "../Geometry/MathVector.h"
#include "../Body/Body.h"

enum{TNEO, TNWO, TSEO, TSWO, BNEO, BNWO, BSEO, BSWO};

class Octree{

 private:
  
  void Split();
  int Determine_Octant(const Vector3D& point);
  bool contains_children, contains_particle;
  int Level;
  Vector3D MP;
  Vector3D OUL;
  Vector3D OLL;
  double s;

 public:

  Octree();
  Octree(const Vector3D& DUL, const Vector3D& DLL);
  ~Octree();

  void Insert(Body& particle);  // Insert a single particle
  Octree* Search(Body& particle); // Searches if a particle is present in the tree, returns pointer to that node
  void Clear();
  int Node_Level(){return Level;}
  void Set_Limits(const Vector3D& DLL, const Vector3D& DUL){
    OLL = DLL;
    OUL = DUL;
  }

  void Write_to_CSV();
  void Write_to_Tec360(const double& t);

  void Compute_Force(Body& particle, const double& Theta_MAX);
  double Theta(const Vector3D& pos); // Return the theta value for this node

  void Partition_Test(const int& N_levels);

  Body* particle_p; // Pointer to a particle, NULL if Octant is empty or has children
  double mass; // stores the mass of all sub-domains
  Vector3D Mass_Center;

  Octree *Parent;

  Octree *TNE;
  Octree *TNW;
  Octree *TSW;
  Octree *TSE;
  Octree *BNE;
  Octree *BNW;
  Octree *BSW;
  Octree *BSE;

};
