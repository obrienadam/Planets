#include "BodySPH.h"
#include <math.h>

void BodySPH::Set_Reference_State(const double& rho_0_init, const double& P_0_init){

  rho_0 = rho_0_init;
  P_0 = P_0_init;

}

double BodySPH::W_Gauss(const double& r, const double& h){

  double r_tilda = r/h;

  return 1.0/(h*h*h*pow(PI, 1.5))*exp(-r_tilda*r_tilda);

}

double BodySPH::W_Cubic_Spline(const double& r, const double& h){

  double r_tilda = r/h;

  if(r_tilda >= 0.0 && r_tilda <= 1.0){

    return (Cn/(h*h*h))*(1 - 1.5*r_tilda*r_tilda + 0.75*r_tilda*r_tilda*r_tilda); // make this more efficient

  } else if (r_tilda > 1.0 && r_tilda <= 2.0){

    return (Cn/(h*h*h))*0.25*(2.0 - r_tilda )*(2.0 - r_tilda)*(2.0 - r_tilda); // make this more efficient

  } else {

    return 0.0;

  }

}

void BodySPH::Compute_SPH_Accelerations(BodySPH *bodies, const int& N, const int& NCl, const int& NCu){

  int i, j;
  double rr;

  // Compute the densities

  for(i = NCl; i < NCu; ++i){
    for(j = 0; j < N; ++j){
    
      rr = Vector3D::Distance_Between_Squared(bodies[i].Xc, bodies[j].Xc);

      if(rr > FOUR*bodies[i].h*bodies[i].h)
	continue;

      bodies[i].rho += bodies[j].mass*bodies[i].W_Cubic_Spline(sqrt(rr), bodies[i].h);
  
    }
  }

  // Compute all other particle properties

  

}
