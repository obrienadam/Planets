#include "Line.h"

Line3D::Line3D(){
  length = dx = dy = dz = 0;
}

Line3D::Line3D(const Vector3D& point1_init, const Vector3D& point2_init){
  
  Create(point1_init, point2_init);
  
}

void Line3D::Create(const Vector3D& point1_init, const Vector3D& point2_init){
  point1 = point1_init;
  point2 = point2_init;

  length = Vector3D::Distance_Between(point1, point2);

  dx = point2.x - point1.x;
  dy = point2.y - point1.y;
  dz = point2.z - point1.z;
}

double Line3D::Area_of_Triangle(const Line3D& line1, const Line3D& line2){
  return 0.5*(Vector3D::Cross_Product(line1.point2 - line1.point1, line2.point2 - line2.point1).Magnitude());
}

double Line3D::Area_of_Quadrilateral(const Line3D& line1, const Line3D& line2){
  return Vector3D::Cross_Product(line1.point2 - line1.point1, line2.point2 - line2.point1).Magnitude();
}
