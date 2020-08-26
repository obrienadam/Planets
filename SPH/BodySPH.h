#pragma once

#include "../Body/Body.h"
#include "../Geometry/MathVector.h"
#include "../Constants/Constants.h"

enum{GAUSS, CUBIC_SPLINE};

class BodySPH : public Body{

 private:

  

 public:

  double h; // Local filter length scale

  double rho_0, P_0, rho, P, g;


  void Set_Reference_State(const double& rho_0_init, const double& P_0_init);

  double W_Gauss(const double& r, const double& h);
  Vector3D Grad_W_Gauss(const double& r, const double& h);
  double W_Cubic_Spline(const double& r, const double& h);
  Vector3D Grad_W_Cubic_Spline(const double& r, const double& h);

  static void Compute_SPH_Accelerations(BodySPH *bodies, const int& N, const int& NCl, const int& NCu);
};
