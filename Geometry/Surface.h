#ifndef _SURFACE_H_
#define _SURFACE_H_

#include "Line.h"
#include <cstdlib>
#include <iostream>

class Surface{
 private:

  int N_vertices;
  Line3D* edges;
  void allocate(const int& N_vertices_init);
  void deallocate();

  double surface_area;
  double perimeter;
  Vector3D center;
  
  void Find_Geometric_Parameters();

 public:
  Surface();
  ~Surface();
  
  double edge_length(const int& edge);
  double Surface_Area(){return surface_area;}

  void Create_Rectangle(const double& dx, const double& dy);
  void Create_Arbitrary_Polygon(const Vector3D* Vertex, const int& N_vertices_init);

  void Output_Info();
};

#endif
