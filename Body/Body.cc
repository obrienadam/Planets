#include "Body.h"
#include <cstdlib>
#include <iostream>
#include <math.h>

Body::Body(){
  mass = 0;
  Xc = v = F = Vector3D(0, 0, 0);
}

Body::Body(const double& mass_init){
  mass = mass_init;
  Xc = v = F = Vector3D(0, 0, 0);
}

Body::Body(const double& mass_init, const Vector3D& Xc_init){
  mass = mass_init;
  Xc = Xc_init; 
  v = F = Vector3D(0, 0, 0);
}

Body::Body(const double& mass_init, const Vector3D& Xc_init, const Vector3D& v_init){
  mass = mass_init;
  Xc = Xc_init; 
  v = v_init;
  F = Vector3D(0, 0, 0);
}

Body::Body(const double& mass_init, const Vector3D& Xc_init, const Vector3D& v_init, const Vector3D& F_init){
mass = mass_init;
  Xc = Xc_init; 
  v = v_init;
  F = F_init;
}

Body::~Body(){
  
}

Body& Body::operator=(const Body &rhs){

  if (this == &rhs)      // Same object?
    return *this; 

  mass = rhs.mass;
  Xc = rhs.Xc;
  v = rhs.v;
  F = rhs.F;
  dXc_dt = rhs.dXc_dt;
  dv_dt = rhs.dv_dt;
  body_name = rhs.body_name;

  return *this;
}

void Body::Output_File(){


}

Vector3D Body::Compute_Force(const Body& Body1, const Body& Body2){

  double r_2 = Vector3D::Distance_Between_Squared(Body1.Xc, Body2.Xc);
  Vector3D force_vector = Vector3D::Relative_Unit_Vector(Body1.Xc, Body2.Xc);
  force_vector.Scale(GRAVITATION_CONSTANT*Body1.mass*Body2.mass/(r_2 + EPSILON_2));

  return force_vector;

}

Vector3D Body::Compute_Softened_Force(const Body& Body1, const Body& Body2, const double& epsilon_sqr){

  double r_2 = Vector3D::Distance_Between_Squared(Body1.Xc, Body2.Xc);

  return (GRAVITATION_CONSTANT*Body1.mass*Body2.mass/pow(epsilon_sqr + r_2, 1.5))*(Body2.Xc - Body1.Xc);

}

void Body::Compute_NBody_Forces(Body *bodies, const int& N){

  int i, j;

  for(i = 0; i < N; ++i){

    bodies[i].F = Vector3D(0, 0, 0);

    for(j = 0; j < N; ++j){
    
      if(i != j){
	bodies[i].F += Compute_Force(bodies[i], bodies[j]);
      }
    
    }

    bodies[i].dv_dt = (1.0/bodies[i].mass)*bodies[i].F;

  }

}



void Body::Compute_NBody_Forces(Body *bodies, const int& N, const int& NCl, const int& NCu){

int i, j;

  for(i = NCl; i <= NCu; ++i){

    bodies[i].F = Vector3D(0, 0, 0);

    for(j = 0; j < N; ++j){
    
      if(i != j){
	bodies[i].F += Compute_Force(bodies[i], bodies[j]);
      }
    
    }

    bodies[i].dv_dt = (1.0/bodies[i].mass)*bodies[i].F;

  }

}

void Body::Compute_Softened_NBody_Forces(Body *bodies, const int& N, const int& NCl, const int& NCu, const double& epsilon_sqr){

int i, j;

  for(i = NCl; i <= NCu; ++i){

    bodies[i].F = Vector3D(0, 0, 0);

    for(j = 0; j < N; ++j){
    
      if(i != j){
	bodies[i].F += Compute_Softened_Force(bodies[i], bodies[j], epsilon_sqr);
      }
    
    }

    bodies[i].dv_dt = (1.0/bodies[i].mass)*bodies[i].F;

  }

}

void Body::Print_Position(){

  std::cout << "(" << Xc.x << ", " << Xc.y << ", " << Xc.z << ")" << std::endl;

}

void Body::Print_Position_LY(){

std::cout << "(" << Xc.x/9.4605284e15 << ", " << Xc.y/9.4605284e15 << ", " << Xc.z/9.4605284e15 << ") LY" << std::endl;

}
