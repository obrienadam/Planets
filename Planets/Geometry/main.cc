#include "Surface.h"
#include <iostream>

using namespace std;

int main()
{
  Surface Test_Surface;
  Vector3D vec1(1, 2, 3), vec2(2,1,4);
  
  Test_Surface.Create_Rectangle(1.2, 1.4);

  Test_Surface.Output_Info();

  return 0;
}
