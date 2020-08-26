#ifndef _BODY_H_
#define _BODY_H_

#include "../Geometry/MathVector.h"
#include <fstream>
#include <string>
#include "../Constants/Constants.h"

enum {VELOCITY, POSITION, VELOCITY_AND_POSITION};

class Body {

 public:
  Body();
  ~Body();
  Body(const double& mass_init);
  Body(const double& mass_init, const Vector3D& Xc_init);
  Body(const double& mass_init, const Vector3D& Xc_init, const Vector3D& v_init);
  Body(const double& mass_init, const Vector3D& Xc_init, const Vector3D& v_init, const Vector3D& F_init);

  double mass;
  Vector3D Xc, v, F;
  Vector3D dXc_dt, dv_dt;

  std::string body_name;
  
  void Output_File();

  Body& operator=(const Body &rhs);

  static Vector3D Compute_Force(const Body& Body1, const Body& Body2);
  static Vector3D Compute_Softened_Force(const Body& Body1, const Body& Body2, const double& epsilon_sqr);

  static void Compute_NBody_Forces(Body *bodies, const int& N);
  static void Compute_NBody_Forces(Body *bodies, const int& N, const int& NCl, const int& NCu);
  static void Compute_Softened_NBody_Forces(Body *bodies, const int& N, const int& NCl, const int& NCu, const double& epsilon_sqr);

  void Print_Position();
  void Print_Position_LY();

};

#endif
