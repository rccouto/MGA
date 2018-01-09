#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defmol.h"
#include "centerxyz.h"
#include "placeatomonplan.h"
#include "placeatomonzaxisforstand.h"
#include "standgeoorien.h"

int StandardGeoOrientation(int natoms, int nindiv, Txyz *XYZ)
{
  int    i, atomonaxis, atomonplane;  
  double xc, yc, zc;
  double abscoord = 1.0E-03;
 
  if ( natoms < 4 ){
    printf("Molecule is too short (Minimum of 4 atoms).");
    exit(EXIT_FAILURE);
  }
  else{
    
    if(  ( ( fabs(XYZ[0].x - XYZ[1].x) < abscoord ) && ( fabs(XYZ[1].x - XYZ[2].x) < abscoord ) )  
      || ( ( fabs(XYZ[0].y - XYZ[1].y) < abscoord ) && ( fabs(XYZ[1].y - XYZ[2].y) < abscoord ) )
      || ( ( fabs(XYZ[0].z - XYZ[1].z) < abscoord ) && ( fabs(XYZ[1].z - XYZ[2].z) < abscoord ) ) ){
      atomonplane = 3;
    }
    else if ( ( ( fabs(XYZ[0].x - XYZ[1].x) < abscoord ) && ( fabs (XYZ[1].x - XYZ[3].x) < abscoord ) )
	   || ( ( fabs(XYZ[0].y - XYZ[1].y) < abscoord ) && ( fabs (XYZ[1].y - XYZ[3].y) < abscoord ) ) 
	   || ( ( fabs(XYZ[0].z - XYZ[1].z) < abscoord ) && ( fabs (XYZ[1].z - XYZ[3].z) < abscoord ) ) ){
      atomonplane = 4;
    }
    else {
      atomonplane = 2;
    }    
    
    xc = XYZ[0].x;
    yc = XYZ[0].y;
    zc = XYZ[0].z;
    atomonaxis = 1;
    
    centerxyz(natoms, xc, yc, zc, XYZ);
    
    PlaceAtmonZAxisForStand(atomonaxis, natoms, XYZ);
    
    PlaceAtmonXZPlan(atomonplane, natoms, XYZ);
    
  }   
}
