#include "Input.h"
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>

Input::Input(){

  // Set default values

  N_Particles = 1000;
  IC_Type = STAR_COLD_COLLAPSE;
  Force_Evaluation = BINARYTREE_FML;
  Time_Integration = LEAP_FROG;

  Mass_Average = 1.981e30;
  Mass_Std_Dev = 0.1*Mass_Average;
  Epsilon = 6.955e8;
  Dt = 1e7;
  Time_Max = Dt*1e6;
  Sphere_Radius = 5.0*9.4605284e15;

  Output_Length_Scale_Modifier = 1.0/9.4605284e15; // Default set to LY
  Output_Time_Scale_Modifier = 1.0/(365*24*3600); // Default set to Year
  Output_Mass_Scale_Modifier = 1.0/1.981e30; // Default set to Solar Mass

  
}

void Input::Remove_White_Space_and_Comments(std::string& buffer){

  if(buffer[0] == '#'){
    buffer = "";
    return;
  }

  buffer.erase(std::remove_if(buffer.begin(), buffer.end(), isspace),buffer.end());
  buffer = buffer.substr(0, buffer.find("#"));
  buffer = buffer.substr(0, buffer.find("="));

}

void Input::Read_Input_File(const char* filename){

  enum{INVALID_PARAMETER_NAME, INVALID_PARAMETER_VALUE, VALID_PARAMETER_NAME, VALID_PARAMETER_VALUE, NO_PARAMETER};

  std::ifstream fin;
  std::string buffer;
  int error_flag = NO_PARAMETER;

  fin.open(filename);

  if(fin.fail()){

    std::cout << "Error: File \"" << filename << "\" was not found." << std::endl;
    return;

  }

  while(!fin.eof()){

    getline(fin, buffer);
    
    Remove_White_Space_and_Comments(buffer);

    if(buffer == ""){

      continue;

    } else if (buffer == "N_Particles"){

      if(fin >> N_Particles){

	error_flag = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;

      }

    } else if (buffer == "IC_Type"){

      fin >> buffer;

      if(buffer == "Star_Cold_Collapse"){

	IC_Type = STAR_COLD_COLLAPSE;
	error_flag  = VALID_PARAMETER_VALUE;

      } else if (buffer == "Disk_Galaxy") {

	IC_Type = DISK_GALAXY;
	error_flag  = VALID_PARAMETER_VALUE;

      } else if (buffer == "Solar_System") {

	IC_Type = SOLAR_SYSTEM;
	error_flag  = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;

      }

    } else if (buffer == "Force_Evaluation"){

      fin >> buffer;

      if(buffer == "BinaryTree_FML"){

	Force_Evaluation = BINARYTREE_FML;
	error_flag  = VALID_PARAMETER_VALUE;

      } else if (buffer == "Octree_FML") {

	Force_Evaluation = OCTREE_FML;
	error_flag  = VALID_PARAMETER_VALUE;

      } else if (buffer == "Direct_Summation") {

	Force_Evaluation = DIRECT_SUM;
	error_flag  = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;

      }

    } else if (buffer == "Mass_Average"){

      if(fin >> Mass_Average){

	error_flag = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;	

      }

    } else if (buffer == "Mass_Std_Dev"){

      if(fin >> Mass_Std_Dev){

	error_flag = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;	

      }

    } else if (buffer == "Domain_Width"){

      if(fin >> Domain_Width){

	error_flag = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;	

      }

    } else if (buffer == "Domain_Height"){

      if(fin >> Domain_Height){

	error_flag = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;	

      }

    } else if (buffer == "Domain_Length"){

      if(fin >> Domain_Length){

	error_flag = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;	

      }

    } else if (buffer == "Time_Integration"){

      fin >> buffer;

      if(buffer == "Runge_Kutta_4"){

	Time_Integration = RUNGE_KUTTA_4;
	error_flag  = VALID_PARAMETER_VALUE;

      } else if (buffer == "Runge_Kutta_2"){

	Time_Integration = RUNGE_KUTTA_2;
	error_flag  = VALID_PARAMETER_VALUE;

      } else if (buffer == "Leap_Frog"){

	Time_Integration = LEAP_FROG;
	error_flag  = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;

      }

    } else if (buffer == "Fixed_Time_Step"){

      if(fin >> Dt){

	error_flag = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;	

      } 

    } else if (buffer == "Time_Max"){

      if(fin >> Time_Max){

	error_flag = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;	

      } 

    } else if (buffer == "Epsilon"){

      if(fin >> Epsilon){

	error_flag = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;	

      }   

    } else if (buffer == "Sphere_Radius"){

      if(fin >> Sphere_Radius){

	error_flag = VALID_PARAMETER_VALUE;

      } else {

	error_flag = INVALID_PARAMETER_VALUE;	

      }   

    } else {

      error_flag = INVALID_PARAMETER_NAME;

    }
    
    
    switch (error_flag){

    case INVALID_PARAMETER_NAME:
      std::cout << "Error: Invalid input parameter name \"" << buffer << "\"." <<std::endl;
      break;
    case INVALID_PARAMETER_VALUE:
      std::cout << "Error: Invalid input parameter value for \"" << buffer << "\"." << std::endl;
      break;

    };

  }

  std::cout << "Reading of solution data complete." << std::endl;

  fin.close();

}
