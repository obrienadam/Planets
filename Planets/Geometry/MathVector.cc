#include <math.h>
#include "MathVector.h"
#include <iostream>

Vector3D::Vector3D(const double& x_init, const double& y_init, const double& z_init){

  x = x_init;
  y = y_init;
  z = z_init;

}

void Vector3D::Scale(const double& a){

  x *= a;
  y *= a;
  z *= a;

}


double Vector3D::Magnitude(){

  return sqrt(x*x + y*y + z*z);

}

Vector3D Vector3D::Unit_Vector(){

  double R = 1.0/Magnitude();

  return Vector3D(R*x, R*y, R*z);

}

void Vector3D::Print(){
 
  std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl;

}

double Vector3D::Distance_Between(const Vector3D& vec1, const Vector3D& vec2){

  return Vector3D(vec2.x - vec1.x, vec2.y - vec1.y, vec2.z - vec1.z).Magnitude();

}

double Vector3D::Distance_Between_Squared(const Vector3D& vec1, const Vector3D& vec2){

  double dx = vec2.x - vec1.x;
  double dy = vec2.y - vec1.y;
  double dz = vec2.z - vec1.z;

  return dx*dx + dy*dy + dz*dz;

}

double Vector3D::Dot_Product(const Vector3D& vec1, const Vector3D& vec2){

  return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;

}

Vector3D Vector3D::Relative_Unit_Vector(const Vector3D& vec1, const Vector3D& vec2){

  Vector3D temp(vec2.x - vec1.x, vec2.y - vec1.y, vec2.z - vec1.z);

  return temp.Unit_Vector();

}




