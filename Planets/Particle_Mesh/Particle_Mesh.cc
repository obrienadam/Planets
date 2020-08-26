#include "Particle_Mesh.h"
#include <cstdlib>
#include <complex>
#include <iostream>
#include <math.h>

const double PI = 4*atan(1);
const std::complex<double> TWO_PI_i = std::complex<double>(0, 2.0*PI);
const std::complex<double> EIGHT_PI_3_i = std::complex<double>(0, 8.0*PI*PI*PI);

Particle_Mesh::Particle_Mesh(){

  mass = NULL;
  rho_ft = NULL;
  Xc = NULL;
  NNi = NNj = NNk = 0;
  allocated = false;

}

void Particle_Mesh::allocate(const int& Nx, const int& Ny, const int& Nz, const Vector3D& lower_bounds, const Vector3D& upper_bounds){

  if(allocated == true)
    deallocate();

  int i, j, k;

  NNi = Nx;
  NNj = Ny;
  NNk = Nz;

  Delta = Vector3D((upper_bounds.x - lower_bounds.x)/double(NNi-1),
		   (upper_bounds.y - lower_bounds.y)/double(NNj-1),
		   (upper_bounds.z - lower_bounds.z)/double(NNk-1));
  
  Vol = Delta.x*Delta.y*Delta.z;

  mass = new double**[NNi];
  rho_ft = new std::complex<double>**[NNi];
  Xc = new Vector3D**[NNi];

  for(i = 0; i < NNi; ++i){

    mass[i] = new double*[NNj];
    rho_ft[i] = new std::complex<double>*[NNj];
    Xc[i] = new Vector3D*[NNj];

    for(j = 0; j < NNj; ++j){

      mass[i][j] = new double[NNk];
      rho_ft[i][j] = new std::complex<double>[NNk];
      Xc[i][j] = new Vector3D[NNk];

      for(k = 0; k < NNk; ++k){
	
	mass[i][j][k] = 0;
	Xc[i][j][k] = lower_bounds + Vector3D(double(i)*Delta.x, double(j)*Delta.y, double(k)*Delta.z);

      }
    }
  }
}

Particle_Mesh::~Particle_Mesh(){
  deallocate();
}

void Particle_Mesh::deallocate(){

  if(allocated == true){

    int i, j;

    for(i = 0; i < NNi; ++i){
      for(j = 0; j < NNj; ++j){
	delete[] mass[i][j];
	delete[] rho_ft[i][j];
	delete[] Xc[i][j];
      }
      
      delete[] mass[i];
      delete[] rho_ft[i];
      delete[] Xc[i];
    }

    delete[] mass;
    delete[] rho_ft;
    delete[] Xc;

    NNi = NNj = NNk = 0;
    Vol = 0;
    Delta = Vector3D(0, 0, 0);
    
    allocated = false;

  }
}

void Particle_Mesh::init(const int& Nx, const int& Ny, const int& Nz, const Vector3D& lower_bounds, const Vector3D& upper_bounds){

  allocate(Nx, Ny, Nz, lower_bounds, upper_bounds);

}

void Particle_Mesh::Find_Nearest_Node(int& i, int& j, int& k, const Vector3D& point){

  if(point.x > Xc[NNi-1][NNj-1][NNk-1].x || point.y > Xc[NNi-1][NNj-1][NNk-1].y || point.z > Xc[NNi-1][NNj-1][NNk-1].z || point.x < Xc[0][0][0].x || point.y < Xc[0][0][0].y || point.z < Xc[0][0][0].z){
       
    std::cout << "Cannot find the nearest node at specified location, point is out of bounds." << std::endl;
    return;
  }

  i = (point.x - Xc[0][0][0].x)/Delta.x;
  j = (point.y - Xc[0][0][0].y)/Delta.y;
  k = (point.z - Xc[0][0][0].z)/Delta.z;
}

void Particle_Mesh::Find_Nearest_Node2(int& i, int& j, int& k, const Vector3D& point){

  if(point.x > Xc[NNi-1][NNj-1][NNk-1].x || point.y > Xc[NNi-1][NNj-1][NNk-1].y || point.z > Xc[NNi-1][NNj-1][NNk-1].z || point.x < Xc[0][0][0].x || point.y < Xc[0][0][0].y || point.z < Xc[0][0][0].z){
       
    std::cout << "Cannot find the nearest node at specified location, point is out of bounds." << std::endl;
    return;
  }
  int l, m, n;
  double R, minR;

  for(l = 0; l < NNi; ++l){
    for(m = 0; m < NNj; ++m){
      for(n = 0; n < NNk; ++n){

	R = Vector3D::Distance_Between(point, Xc[l][m][n]);

	if(R < minR || l+m+n == 0){
	  minR = R;
	  i = l;
	  j = m;
	  k = n;
	}

      }
    }
  }
}

void Particle_Mesh::Reset_Mass(){
  int i, j, k;

  for(i = 0; i < NNi; ++i){
    for(j = 0; j < NNj; ++j){
      for(k = 0; k < NNk; ++k){
	mass[i][j][k] = 0;
      }
    }
  }
  std::cout << "RESET \n";
}

void Particle_Mesh::Restrict_Particles_to_Grid(const Body *particles, const int& N){

  int i, j, k;
  int iPar;
  double ONE_OVER_CELL_VOLUME = 1.0/(Delta.x*Delta.y*Delta.z);

  
  Reset_Mass();
      
  for(iPar = 0; iPar < N; ++iPar){
    Find_Nearest_Node(i, j, k, particles[iPar].Xc);
    mass[i][j][k] += particles[iPar].mass*ONE_OVER_CELL_VOLUME;

    std::cout << "Node (" << i << ", " << j << ", " << k << ") Mass: " << mass[i][j][k] << std::endl
	      << "Particle Location: " << particles[iPar].Xc.x << ", " << particles[iPar].Xc.y << ", " << particles[iPar].Xc.z << std::endl
	      << "Node Location: " << Xc[i][j][k].x << ", " << Xc[i][j][k].y << ", " << Xc[i][j][k].z << std::endl;
  }

}


void Particle_Mesh::Compute_DFT(){

  int k1, k2, k3, n1, n2, n3;
  const double ONE_OVER_NTOT = 1.0/(double(NNi*NNj*NNk));

  for(n1 = 0; n1 < NNi; ++n1){
    for(n2 = 0; n2 < NNj; ++n2){
      for(n3 = 0; n3 < NNk; ++n3){

	rho_ft[n1][n2][n3] = std::complex<double>(0, 0);

	for(k1 = 0; k1 < NNi; ++k1){
	  for(k2 = 0; k2 < NNj; ++k2){
	    for(k3 = 0; k3 < NNk; ++k3){

	      if(mass[k1][k2][k3] == 0)
		continue;

	      rho_ft[n1][n2][n3] += mass[k1][k2][k3]*exp(-EIGHT_PI_3_i*double(k3)*double(n3)*double(k2)*double(n2)*double(k1)*double(n1)*ONE_OVER_NTOT);

	    }
	  }
	}

	rho_ft[n1][n2][n3] *= ONE_OVER_NTOT;
    
	std::cout << rho_ft[n1][n2][n3] << std::endl;
      }
    }
  }
}
