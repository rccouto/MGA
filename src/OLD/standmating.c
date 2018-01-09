#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defmol.h"
#include "centerxyz.h"
#include "placeatomonzaxisforstand.h"
#include "standmating.h"

int StandMating(int natoms, int centeratom, int atomonaxis, Txyz *XYZ)
{
  int i;
  double xc, yc, zc;

    xc = XYZ[centeratom].x;
    yc = XYZ[centeratom].y;
    zc = XYZ[centeratom].z;
    /*printf("ENTER\n");    
    for (i = 0; i < natoms; i++){
      printf("\t SYMB    %lf    %lf    %lf\n", XYZ[i].x, XYZ[i].y, XYZ[i].z);
      } */   
    centerxyz(natoms, xc, yc, zc, XYZ);
    /* printf("CENTER\n");
    for (i = 0; i < natoms; i++){
      printf("\t SYMB    %lf    %lf    %lf\n", XYZ[i].x, XYZ[i].y, XYZ[i].z);
      } */       
    PlaceAtmonZAxisForStand(atomonaxis, natoms, XYZ);
    /*printf("PLACE\n");
    for (i = 0; i < natoms; i++){
      printf("\t SYMB    %lf    %lf    %lf\n", XYZ[i].x, XYZ[i].y, XYZ[i].z);
      } */   
}
