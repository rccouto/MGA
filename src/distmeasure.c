#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "defmol.h"
#include "distmeasure.h"

#define ZERO 0.0E+0
#define MAX  9999999999

int DistMeasure(char *prog, int natoms, Txyz *XYZ, double *energy, Tcnct *cnct)
{
  int    i, act;
  double R1, R;
  double deltax, deltay, deltaz;
  
  act = cnct[1].cdist;

  deltax = XYZ[1].x - XYZ[act - 1].x;
  deltay = XYZ[1].y - XYZ[act - 1].y;
  deltaz = XYZ[1].z - XYZ[act - 1].z;
  R1= pow(deltax,2) +  pow(deltay,2) + pow(deltaz,2);
  R = sqrt(R1);

  if ( R > 1.8 ){ 
    if ( strncasecmp(prog,"OBabel", 6) ) *energy = ZERO;
    else *energy = MAX;
    return(0);
  }

  act = cnct[2].cdist;

  deltax = XYZ[2].x - XYZ[act - 1].x;
  deltay = XYZ[2].y - XYZ[act - 1].y;
  deltaz = XYZ[2].z - XYZ[act - 1].z;
  R1= pow(deltax,2) +  pow(deltay,2) + pow(deltaz,2);
  R=sqrt(R1);

  if ( R > 1.8 ){ 
    if ( strncasecmp(prog,"OBabel", 6) ) *energy = ZERO;
    else *energy = MAX;
    return(0);
  }
 
  for( i = 3; i < natoms; i++){

    act = cnct[i].cdist;

    deltax = XYZ[i].x - XYZ[act - 1].x;
    deltay = XYZ[i].y - XYZ[act - 1].y;
    deltaz = XYZ[i].z - XYZ[act - 1].z;
    R1= pow(deltax,2) +  pow(deltay,2) + pow(deltaz,2);
    R=sqrt(R1);

    if ( R > 1.8 ){ 
      if ( strncasecmp(prog,"OBabel", 6) ) *energy = ZERO;
      else *energy = MAX;
      return(0);
    } 
  }
}
