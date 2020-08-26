#include "Pstate.h"
#include "Cstate.h"
#include <math.h>

Pstate::Pstate(const double& rho_init, const Vector3D& v_init, const double& P_init){

  double GAMMA = 1.41, GAS_CONSTANT = 286.9;

  rho = rho_init;
  v = v_init;
  P = P_init;

  htot = P/((GAMMA-1.0)*rho) + 0.5*(v.x*v.x + v.y*v.y + v.z*v.z) + P/rho;
  a = sqrt(GAMMA*P/rho);
  T = P/(rho*GAS_CONSTANT);

}

Pstate::Pstate(const Cstate& U_init){

double GAMMA = 1.41, GAS_CONSTANT = 286.9;

  rho = U_init.rho;
  v = (1.0/rho)*U_init.mom;
  P = (U_init.e - 0.5*(v.x*v.x + v.y*v.y + v.z*v.z))*(GAMMA-1.0)*rho;

  htot = U_init.e + P/rho;
 a = sqrt(GAMMA*P/rho);
  T = P/(rho*GAS_CONSTANT);
}

Cstate Pstate::U(){

  return Cstate(*this);

}
