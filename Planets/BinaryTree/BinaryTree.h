#pragma once

#include "../Geometry/MathVector.h"
#include "../Body/Body.h"

enum{LEFT, RIGHT, EW, NS, TB, NONE};

class BinaryTree{

 private:

  void Split(const int& direction);
  int Determine_Split_Direction(const Vector3D& point1, const Vector3D& point2);
  int Determine_Half(const Vector3D& point);

  bool contains_children, contains_particle;
  int split_direction, Level;

  Vector3D MP;
  Vector3D HLL, HUL, Delta;
  double s;

 public:

  BinaryTree(const Vector3D& DLL = Vector3D(0,0,0), const Vector3D& DUL = Vector3D(0,0,0));
  ~BinaryTree(){Clear();}

  void Set_Limits(const Vector3D& DLL, const Vector3D& DUL);

  void Insert(Body& particle);
  void Clear();

  void Partition_Test();

  void Write_to_Tec360(const double& t);

  void Compute_Softened_Force(Body& particle, const double& theta_MAX, const double& epsilon_sqr);

  double mass;
  Vector3D Mass_Center;

  Body* particle_p;

  BinaryTree* parent;

  BinaryTree* Left;
  BinaryTree* Right;

};
