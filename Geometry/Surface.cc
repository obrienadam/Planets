#include "Surface.h"
#include <iostream>

Surface::Surface()
{
  edges = NULL;
}

Surface::~Surface()
{
  deallocate();
}

void Surface::allocate(const int& N_vertices_init)
{
  deallocate();

  N_vertices = N_vertices_init;
  edges = new Line3D[N_vertices_init];
}

void Surface::deallocate()
{
  if(edges != NULL)
    {
      delete[] edges;
      edges = NULL;
    }
  surface_area = perimeter = 0;
  center = Vector3D(0, 0);
}

void Surface::Create_Rectangle(const double& dx, const double& dy){

  allocate(4);

  edges[0] = Line3D(Vector3D(0, 0), Vector3D(dx, 0));
  edges[1] = Line3D(Vector3D(dx, 0), Vector3D(dx, dy));
  edges[2] = Line3D(Vector3D(dx, dy), Vector3D(0, dy));
  edges[3] = Line3D(Vector3D(0, dy), Vector3D(0, 0));

  Find_Geometric_Parameters();
}

void Surface::Find_Geometric_Parameters()
{
  perimeter = 0;
  surface_area = 0;

  for(int i = 0; i < N_vertices; ++i)
    {
      perimeter += edges[i].Length();
      
      if(i != N_vertices-1)
	{
	  surface_area += Line3D::Area_of_Triangle(edges[0], edges[i+1]);
	}
    }
}

void Surface::Output_Info()
{
  std::cout << "Surface Info: " << std::endl
	    << "N Edges: " << N_vertices << std::endl
	    << "Perimeter: " << perimeter << std::endl
	    << "Surface Area: " << surface_area << std::endl;
}

void Surface::Create_Arbitrary_Polygon(const Vector3D* Vertex, const int& N_vertices_init){
  
  if(N_vertices_init <= 2){
    std::cerr << "Error in Surface.cc, cannot create a polygon with fewer than 3 vertices." << std::endl;
  }
  
  allocate(N_vertices_init);
  
  for(int i = 0; i < N_vertices; ++i){
    
    if (i == N_vertices-1) {

      edges[i] = Line3D(Vertex[i], Vertex[0]);

    } else {

      edges[i] = Line3D(Vertex[i], Vertex[i+1]);

    }
  }
  
  Find_Geometric_Parameters();
}
