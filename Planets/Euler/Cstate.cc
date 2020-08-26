#include "Cstate.h"
#include "Pstate.h"

Cstate::Cstate(const Pstate& W_init){

  rho = W_init.rho;
  mom = rho*W_init.v;
  e = W_init.htot - W_init.P/rho;

}

Cflux::Cflux(const double& e_init, const Vector3D& mom_init, const double& rho_init){

  e = e_init;
  rho = rho_init;
  mom = mom_init;

}

Cflux::Cflux(const Cstate& q, const Vector3D& unit_normal){

  

}

Pstate Cstate::W(){

  return Pstate(*this);

}

