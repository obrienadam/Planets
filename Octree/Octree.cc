#include "Octree.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

Octree::Octree(){

  Parent = TNE = TNW = TSW = TSE = BNE = BNW = BSW = BSE = NULL;
  mass = 0;
  Level = 0;

  contains_children = false;
  contains_particle = false;
}

Octree::Octree(const Vector3D& DLL, const Vector3D& DUL){

  Parent = TNE = TNW = TSW = TSE = BNE = BNW = BSW = BSE = NULL;
  mass = 0;
  Level = 0;

  OUL = DUL;
  OLL = DLL;

  MP = 0.5*(OUL + OLL);

  contains_children = false;
  contains_particle = false;

}

Octree::~Octree(){

  Clear();

}

void Octree::Clear(){

  if(contains_children){
    
    TNE->Clear();
    TNW->Clear();
    TSW->Clear();
    TSE->Clear();
    BNE->Clear();
    BNW->Clear();
    BSW->Clear();
    BSE->Clear();

    delete TNE;
    delete TNW;
    delete TSW;
    delete TSE;
    delete BNE;
    delete BNW;
    delete BSW;
    delete BSE;

    contains_children = false;
    mass = 0;
    Mass_Center = Vector3D(0, 0, 0);
  }

}

void Octree::Split(){

  TNE = new Octree;
  TNW = new Octree;
  TSE = new Octree;
  TSW = new Octree;
  BNE = new Octree;
  BNW = new Octree;
  BSE = new Octree;
  BSW = new Octree;

  TNE->OUL = OUL;
  TNE->OLL = 0.5*(OUL + OLL);
  TNE->MP = 0.5*(TNE->OUL + TNE->OLL);

  TNW->OUL = Vector3D(0.5*(OUL.x + OLL.x), OUL.y, OUL.z);
  TNW->OLL = Vector3D(OLL.x, 0.5*(OUL.y + OLL.y), 0.5*(OUL.z + OLL.z));
  TNW->MP = 0.5*(TNW->OUL + TNW->OLL);

  TSW->OUL = Vector3D(0.5*(OUL.x + OLL.x), 0.5*(OUL.y + OLL.y), OUL.z);
  TSW->OLL = Vector3D(OLL.x, OLL.y, 0.5*(OUL.z + OLL.z));
  TSW->MP = 0.5*(TSW->OUL + TSW->OLL);

  TSE->OUL = Vector3D(OUL.x, 0.5*(OUL.y + OLL.y), OUL.z);
  TSE->OLL = Vector3D(0.5*(OUL.x + OLL.x), OLL.y, 0.5*(OUL.z + OLL.z));
  TSE->MP = 0.5*(TSE->OUL + TSE->OLL);

  BNE->OUL = Vector3D(OUL.x, OUL.y, 0.5*(OUL.z + OLL.z));
  BNE->OLL = Vector3D(0.5*(OUL.x + OLL.x), 0.5*(OUL.y + OLL.y), OLL.z);
  BNE->MP = 0.5*(BNE->OUL + BNE->OLL);

  BNW->OUL = Vector3D(0.5*(OUL.x + OLL.x), OUL.y, 0.5*(OUL.z + OLL.z));
  BNW->OLL = Vector3D(OLL.x, 0.5*(OUL.y + OLL.y), OLL.z);
  BNW->MP = 0.5*(BNW->OUL + BNW->OLL);

  BSW->OUL = Vector3D(0.5*(OUL.x + OLL.x), 0.5*(OUL.y + OLL.y), 0.5*(OUL.z + OLL.z));
  BSW->OLL = OLL;
  BSW->MP = 0.5*(BSW->OUL + BSW->OLL);

  BSE->OUL = Vector3D(OUL.x, 0.5*(OUL.y + OLL.y), 0.5*(OUL.z + OLL.z));
  BSE->OLL = Vector3D(0.5*(OUL.x + OLL.x), OLL.y, OLL.z);
  BSE->MP = 0.5*(BSE->OUL + BSE->OLL);

  TNE->Level = TNW->Level = TSE->Level = TSW->Level = BNE->Level = BNW->Level = BSE->Level = BSW->Level = Level+1;
  TNE->Parent = TNW->Parent = TSE->Parent = TSW->Parent = BNE->Parent = BNW->Parent = BSE->Parent = BSW->Parent = this;

  contains_children = true;
}

int Octree::Determine_Octant(const Vector3D& point){

  if(point.z > MP.z){ // Top

    if(point.y > MP.y){ // North

      if(point.x > MP.x){ // East

	return TNEO;

      } else { // West

	return TNWO;

      }

    } else { // South

      if(point.x > MP.x){ // East

	return TSEO;

      } else { // West

	return TSWO;

      }

    }

  }else{ // Bottom

    if(point.y > MP.y){ // North

      if(point.x > MP.x){ // East

	return BNEO;

      } else { // West

	return BNWO;

      }

    } else { // South

      if(point.x > MP.x){ // East

	return BSEO;

      } else { // West

	return BSWO;

      }

    }

  }

}

void Octree::Write_to_CSV(){

  if(contains_children){
    TNE->Write_to_CSV();
    TNW->Write_to_CSV();
    TSW->Write_to_CSV();
    TSE->Write_to_CSV();
    BNE->Write_to_CSV();
    BNW->Write_to_CSV();
    BSW->Write_to_CSV();
    BSE->Write_to_CSV();
  } else {
    std::ofstream fout;

    fout.open("Octree.csv", std::ofstream::app);

    fout << OLL.x << "," << OLL.y << "," << OLL.z << std::endl
	 << OUL.x << "," << OLL.y << "," << OLL.z << std::endl
	 << OUL.x << "," << OUL.y << "," << OLL.z << std::endl
	 << OLL.x << "," << OUL.y << "," << OLL.z << std::endl
	 << OLL.x << "," << OLL.y << "," << OLL.z << std::endl
	 << OLL.x << "," << OLL.y << "," << OUL.z << std::endl
	 << OUL.x << "," << OLL.y << "," << OUL.z << std::endl
	 << OUL.x << "," << OUL.y << "," << OUL.z << std::endl
	 << OLL.x << "," << OUL.y << "," << OUL.z << std::endl
	 << OLL.x << "," << OLL.y << "," << OUL.z << std::endl;

      fout.close();
      }

}

void Octree::Write_to_Tec360(const double& t){

  std::stringstream filename;

  filename << "Octree_t_" << t << ".dat";

  if(Level == 0){
    std::ofstream fout;
    
    fout.open(filename.str().c_str());

    fout << "TITLE = \"Planets Octree\"" << std::endl
	 << "VARIABLES = \"X\", \"Y\", \"Z\", \"P\"" << std::endl;
  
    fout.close(); 
  }

  if(contains_children){
    TNE->Write_to_Tec360(t);
    TNW->Write_to_Tec360(t);
    TSW->Write_to_Tec360(t);
    TSE->Write_to_Tec360(t);
    BNE->Write_to_Tec360(t);
    BNW->Write_to_Tec360(t);
    BSW->Write_to_Tec360(t);
    BSE->Write_to_Tec360(t);
  } else {
    std::ofstream fout;

    fout.open(filename.str().c_str(), std::ofstream::app);

    int physical_particle;
    
    if(particle_p != NULL)
      {
	physical_particle = 1;
      }
    else
      {
	physical_particle = 0;
	  }

    fout << "ZONE T = \"Level" << Level << "\", I=2, J=2, K=2, F=POINT, STRANDID=2, SOLUTIONTIME=" << t << std::endl
	 << OLL.x << " " << OLL.y << " " << OLL.z << " " << physical_particle << std::endl
	 << OUL.x << " " << OLL.y << " " << OLL.z << " " << physical_particle << std::endl
	 << OLL.x << " " << OUL.y << " " << OLL.z << " " << physical_particle << std::endl
	 << OUL.x << " " << OUL.y << " " << OLL.z << " " << physical_particle << std::endl
	 << OLL.x << " " << OLL.y << " " << OUL.z << " " << physical_particle << std::endl
	 << OUL.x << " " << OLL.y << " " << OUL.z << " " << physical_particle << std::endl
	 << OLL.x << " " << OUL.y << " " << OUL.z << " " << physical_particle << std::endl
	 << OUL.x << " " << OUL.y << " " << OUL.z << " " << physical_particle << std::endl;

      fout.close();
  }

}

void Octree::Partition_Test(const int& N_levels){

  int test_level = N_levels;

  Split();

  TNE->Split();
  TNE->TNE->Split();
  TNE->TNE->TNE->Split();
  TNE->TNE->BSW->Split();

  BSW->Split();
  BSW->BSE->Split();

  std::cout << "Partition testing complete." << std::endl;
}

void Octree::Insert(Body& particle){

  if(particle.Xc.x > OUL.x || particle.Xc.y > OUL.y || particle.Xc.z > OUL.z ||
     particle.Xc.x < OLL.x || particle.Xc.y < OLL.y || particle.Xc.z < OLL.z ){

    std::cout << "Error: Cannot insert particle. Physical boundary exceeded..." << std::endl;
    return;

  }


  if(contains_particle == false){ // The particle can go here!

    particle_p = &(particle);
    mass = particle.mass;
    Mass_Center = particle.Xc;
    contains_particle = true;

  } else { // No room for particle, ie particles exist in this octant

    int Octant;
    double mass_old;

    // Adjust the mass center

    mass_old = mass;
    mass += particle.mass;
    Mass_Center *= mass_old;
    Mass_Center += particle.mass*particle.Xc;
    Mass_Center /= mass;

    if(contains_children == false)
      Split();

    Octant = Determine_Octant(particle.Xc);

    switch (Octant){
    case TNEO:
      TNE->Insert(particle);
      break;
    case TNWO:
      TNW->Insert(particle);
      break;
    case TSEO:
      TSE->Insert(particle);
      break;
    case TSWO:
      TSW->Insert(particle);
      break;
    case BNEO:
      BNE->Insert(particle);
      break;
    case BNWO:
      BNW->Insert(particle);
      break;
    case BSEO:
      BSE->Insert(particle);
      break;
    case BSWO:
      BSW->Insert(particle);
      break;
    };

    if(particle_p != NULL){
    
    Octant = Determine_Octant((*particle_p).Xc);

    switch (Octant){
    case TNEO:
      TNE->Insert(*particle_p);
      break;
    case TNWO:
      TNW->Insert(*particle_p);
      break;
    case TSEO:
      TSE->Insert(*particle_p);
      break;
    case TSWO:
      TSW->Insert(*particle_p);
      break;
    case BNEO:
      BNE->Insert(*particle_p);
      break;
    case BNWO:
      BNW->Insert(*particle_p);
      break;
    case BSEO:
      BSE->Insert(*particle_p);
      break;
    case BSWO:
      BSW->Insert(*particle_p);
      break;
    };

    particle_p = NULL;

    }

  }

}

Octree* Octree::Search(Body& particle){

  if(contains_particle && &particle == particle_p){ // particle found!

    return this;

  } else if (contains_particle == false) {

    std::cout << "Error, particle not found in Octree." << std::endl;
    return NULL;

  } else {

    int Octant = Determine_Octant(particle.Xc);

    switch (Octant){

    case TNEO:
      return TNE->Search(particle);
      break;
    case TNWO:
      return TNW->Search(particle);
      break;
    case TSWO:
      return TSW->Search(particle);
      break;
    case TSEO:
      return TSE->Search(particle);
      break;
    case BNEO:
      return BNE->Search(particle);
      break;
    case BNWO:
      return BNW->Search(particle);
      break;
    case BSWO:
      return BSW->Search(particle);
      break;
    case BSEO:
      return BSE->Search(particle);
      break;
    };

  }

  return NULL;
}

double Octree::Theta(const Vector3D& pos){

  double d = Vector3D::Distance_Between(pos, MP);
  s = Vector3D::Distance_Between(OLL, OUL);

  return s/d;

}

void Octree::Compute_Force(Body& particle, const double& Theta_MAX){

  if(contains_particle == false){ // empty octant
    
    return;

  } else if (particle_p == &particle){

    return;

  } else if(Level == 0) { // root node, reset force

    particle.F = Vector3D(0, 0, 0);
  
  }

  if(contains_children == false){ // cannot go finer

    particle.F += Body::Compute_Force(particle, Body(mass, Mass_Center)); 

  } else if(Theta(particle.Xc) > Theta_MAX){ // Must go finer

    TNE->Compute_Force(particle, Theta_MAX);
    TNW->Compute_Force(particle, Theta_MAX);
    TSE->Compute_Force(particle, Theta_MAX);
    TSW->Compute_Force(particle, Theta_MAX);
    BNE->Compute_Force(particle, Theta_MAX);
    BNW->Compute_Force(particle, Theta_MAX);
    BSE->Compute_Force(particle, Theta_MAX);
    BSW->Compute_Force(particle, Theta_MAX);

  } else {

    particle.F += Body::Compute_Force(particle, Body(mass, Mass_Center));

  }

}
