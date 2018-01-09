#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defmol.h"
#include "centerxyz.h"

int centerxyz(const int natoms, const double x, const double y, const double z, Txyz *XYZ)
{
  //  printf("CENTER XYZ\n");
  //  printf("natoms=%d, x=%lf, y=%lf, z=%lf\n", natoms, x, y, z);
  int i;

  for(i = 0; i < natoms; i++ ){
    XYZ[i].x = XYZ[i].x - x;
    XYZ[i].y = XYZ[i].y - y;
    XYZ[i].z = XYZ[i].z - z;
    
    //    printf("\t %lf  %lf  %lf\n", XYZ[i].x,  XYZ[i].y,  XYZ[i].z);

  }
  return(0);
}
