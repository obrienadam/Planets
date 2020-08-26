#ifndef _MATH_VECTOR_H_
#define _MATH_VECTOR_H_

class Vector3D{
 public:
  Vector3D(const double& x_init = 0, const double& y_init = 0, const double& z_init = 0);
  double x, y, z;
  
  void Scale(const double& a);
  double Magnitude();
  Vector3D Unit_Vector();
  void Print();

  static double Distance_Between(const Vector3D& Vec1, const Vector3D& Vec2);
  static double Distance_Between_Squared(const Vector3D& vec1, const Vector3D& vec2);
  static double Dot_Product(const Vector3D& Vec1, const Vector3D& Vec2);

  static Vector3D Cross_Product(const Vector3D& Vec1, const Vector3D& Vec2);
  static Vector3D Relative_Unit_Vector(const Vector3D& Vec1, const Vector3D& Vec2);

  Vector3D& operator+=(const Vector3D& rhs)
    {
      x += rhs.x;
      y += rhs.y;
      z += rhs.z;

      return *this;
    }
  
  Vector3D& operator-=(const Vector3D& rhs)
    {
      x -= rhs.x;
      y -= rhs.y;
      z -= rhs.z;

      return *this;
    }

  Vector3D& operator*=(const double& rhs){
    x *= rhs;
    y *= rhs;
    z *= rhs;

    return *this;
  }

  Vector3D& operator/=(const double& rhs){
    x /= rhs;
    y /= rhs;
    z /= rhs;

    return *this;
  }
};

inline Vector3D operator+(Vector3D lhs, const Vector3D& rhs)
{
  lhs += rhs;
  return lhs;
}

inline Vector3D operator-(Vector3D lhs, const Vector3D& rhs)
{
  lhs -= rhs;
  return lhs;
}

inline Vector3D operator*(const double& lhs, Vector3D rhs){

  rhs *= lhs;
  return rhs;
}

inline Vector3D operator*(Vector3D lhs, const double& rhs){

  lhs *= rhs;
  return lhs;

}

inline Vector3D operator/(Vector3D lhs, const double& rhs){

  lhs /= rhs;
  return lhs;

}
#endif
