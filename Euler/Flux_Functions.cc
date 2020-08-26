#include "Flux_Functions.h"
#include "Pstate.h"
#include "Cstate.h"
#include <math.h>

Cflux Roe(const Pstate& L, const Pstate& R, const Vector3D& unit_normal){

  double RT, rhobar, ubar, vbar, wbar, hbar, abar, Vnbar, VnbarL, VnbarR, drho, dP, dVnbar, du, dv, dw;
  double lambda[4], LdU[4], RdU[5][4], diss[5];
  Cflux num_flux;

  RT = sqrt(R.rho/L.rho);
  VnbarL = unit_normal.x*L.v.x + unit_normal.y*L.v.y + unit_normal.z*L.v.z;
  VnbarR = unit_normal.x*R.v.x + unit_normal.y*R.v.y + unit_normal.z*R.v.z;
  
  // Calculate rho averages

  rhobar = RT*L.rho;

  ubar = (L.v.x + RT*R.v.x)/(ONE + RT);
  vbar = (L.v.y + RT*R.v.y)/(ONE + RT);
  wbar = (L.v.z + RT*R.v.z)/(ONE + RT);
  hbar = (L.htot + RT*R.htot)/(ONE + RT);
  abar = sqrt((GAMMA-1)*(hbar - HALF*(ubar*ubar + vbar*vbar + wbar*wbar)));
  Vnbar = unit_normal.x*ubar + unit_normal.y*vbar + unit_normal.z*wbar;

  // Calculate wave strengths

  drho = R.rho - L.rho;
  dP = R.P - L.P;
  dVnbar = VnbarR - VnbarL;

  LdU[0] = (dP - rhobar*abar*dVnbar )/(TWO*abar*abar); //Left-moving acoustic wave strength
  LdU[1] =  drho - dP/(abar*abar);            //Entropy wave strength
  LdU[2] = (dP + rhobar*abar*Vnbar )/(TWO*abar*abar); //Right-moving acoustic wave strength
  LdU[3] = rhobar;                         //Shear wave strength (not really, just a factor)

  // Absolute values of the wave Speeds

  lambda[0] = fabs(Vnbar-abar); //Left-moving acoustic wave
  lambda[1] = fabs(Vnbar);   //Entropy wave
  lambda[2] = fabs(Vnbar+abar); //Right-moving acoustic wave
  lambda[3] = fabs(Vnbar);   //Shear waves

  //Entropy fix

  if ( lambda[0] < FIFTH ){ 
    lambda[0] = HALF * ( lambda[0]*lambda[0]/FIFTH + FIFTH );
  }

  if ( lambda[2] < FIFTH ){ 
    lambda[2] = HALF * ( lambda[2]*lambda[2]/FIFTH + FIFTH );
  }

  /* Right Eigenvectors
     !Note: Two shear wave components are combined into one, so that tangent vectors
     !      are not required. And that's why there are only 4 vectors here.
     !      See "I do like CFD, VOL.1, Second Edition" about how tangent vectors are eliminated.
  */

  // Left-moving acoustic wave

  RdU[0][0] = ONE;    
  RdU[1][0] = ubar - abar*unit_normal.x;
  RdU[2][0] = vbar - abar*unit_normal.y;

  RdU[3][0] = wbar - abar*unit_normal.z;
  RdU[4][0] = hbar - abar*Vnbar;

  RdU[0][1] = ONE;
  RdU[1][1] = ubar;
  RdU[2][1] = vbar; 
  RdU[3][1] = wbar;
  RdU[4][1] = HALF*(ubar*ubar + vbar*vbar + wbar*wbar);

  RdU[0][2] = ONE;
  RdU[1][2] = ubar + abar*unit_normal.x;
  RdU[2][2] = vbar + abar*unit_normal.y;
  RdU[3][2] = wbar + abar*unit_normal.z;
  RdU[4][2] = hbar + abar*Vnbar;

  du = R.v.x - L.v.x;
  dv = R.v.y - L.v.y;
  dw = R.v.z - L.v.z;

  RdU[0][3] = ZERO;
  RdU[1][3] = du - dVnbar*unit_normal.x;
  RdU[2][3] = dv - dVnbar*unit_normal.y;
  RdU[3][3] = dw - dVnbar*unit_normal.z;
  RdU[4][3] = ubar*du + vbar*dv + wbar*dw - Vnbar*dVnbar;

  for(int i = 0; i < 5; ++i){

    diss[0] = lambda[0]*LdU[0]*RdU[i][0] + lambda[1]*LdU[1]*RdU[i][1] + lambda[2]*LdU[2]*RdU[i][2] + lambda[3]*LdU[3]*RdU[i][3];

  }

  num_flux.rho = HALF*(L.rho*VnbarL + R.rho*VnbarR - diss[0]);
  num_flux.mom.x = HALF*(L.rho*VnbarL*L.v.x + L.P*unit_normal.x + R.rho*VnbarR*R.v.x + R.P*unit_normal.x - diss[1]); 
  num_flux.mom.y = HALF*(L.rho*VnbarL*L.v.y + L.P*unit_normal.y + R.rho*VnbarR*R.v.y + R.P*unit_normal.y - diss[2]); 
  num_flux.mom.z = HALF*(L.rho*VnbarL*L.v.z + L.P*unit_normal.z + R.rho*VnbarR*R.v.z + R.P*unit_normal.z - diss[3]);
  num_flux.e = HALF*(L.rho*VnbarL*L.htot + R.rho*VnbarR*R.htot);

    return num_flux;
}


