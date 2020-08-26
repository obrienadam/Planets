#ifndef _LINE_H_
#define _LINE_H_

#include "MathVector.h"

class Line3D{
 private:
  Vector3D point1, point2;
  double length, dx, dy, dz;

 public:

  Line3D();
  Line3D(const Vector3D& point1_init, const Vector3D& point2_init);
  void Create(const Vector3D& point1_init, const Vector3D& point2_init);
  
  double Length(){return length;}
  double Dx(){return dx;}
  double Dy(){return dy;}
  double Dz(){return dz;}
  Vector3D Line_Equation();
  
  static double Area_of_Triangle(const Line3D& line1, const Line3D& line2);
  static double Area_of_Quadrilateral(const Line3D& line1, const Line3D& line2);

  Line3D& operator+=(const Line3D& rhs)
    {
      point1 += rhs.point1;
      point2 += rhs.point2;
      return *this;
    }
  
  Line3D& operator-=(const Line3D& rhs)
    {
      point1 -= rhs.point1;
      point2 -= rhs.point2;
      return *this;
    }
};

inline Line3D operator+(Line3D lhs, const Line3D& rhs)
{
  lhs += rhs;
  return lhs;
}

inline Line3D operator-(Line3D lhs, const Line3D& rhs)
{
  lhs -= rhs;
  return lhs;
}

#endif
