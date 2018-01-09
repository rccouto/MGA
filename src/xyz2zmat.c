#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "defmol.h"
#include "xyz2zmat.h"
#include "placeatomonzaxis.h"
#include "placeatomonplan.h"
#include "centerxyz.h"

#define ZERO          0.0E+0

extern double pi;

int xyz2zmat(int const natoms, Txyz *XYZ, Tcnct *cnct, Tzmat *zmat)
{
  int    j, act, ang, die;
  double X, Y, Z, xc, yc, zc, cv=180/pi; 
  double deltax, deltay, deltaz, R, R2;
  double angle, dihedral;

  /* First Atom*/

  zmat[0].dist = ZERO;
  zmat[0].angle = ZERO;
  zmat[0].dihedral = ZERO;
  
  /* Second Atom*/

  act = cnct[1].cdist;

  deltax = XYZ[1].x - XYZ[act - 1].x;
  deltay = XYZ[1].y - XYZ[act - 1].y;
  deltaz = XYZ[1].z - XYZ[act - 1].z;
  R2 = pow(deltax,2) +  pow(deltay,2) + pow(deltaz,2);
  R = sqrt( R2 );

  zmat[1].dist = R;
  zmat[1].angle = ZERO;
  zmat[1].dihedral = ZERO;
  
  /*Third Atom*/

  act = cnct[2].cdist;
  ang = cnct[2].cangle;
 
  xc = XYZ[act - 1].x;
  yc = XYZ[act - 1].y;
  zc = XYZ[act - 1].z;
  
  centerxyz(natoms, xc, yc, zc, XYZ);
  
  PlaceAtmonZAxis(ang - 1, natoms - 1, XYZ);
  
  R2 =  XYZ[2].x * XYZ[2].x + XYZ[2].y * XYZ[2].y + XYZ[2].z * XYZ[2].z;
  R = sqrt( R2 );
  
  angle = acos( XYZ[ 2].z / R );

  if ( XYZ[ang - 1].z < ZERO ) angle = pi - angle;
  
  zmat[2].dist = R;
  zmat[2].angle = angle*cv;
  zmat[2].dihedral = ZERO;

  for ( j = 3; j < natoms; j++ ){

   act = cnct[j].cdist;
   ang = cnct[j].cangle;
   die = cnct[j].cdihedral;
    
    xc = XYZ[act - 1].x;
    yc = XYZ[act - 1].y;
    zc = XYZ[act - 1].z;
    
    centerxyz(natoms, xc, yc, zc, XYZ);
    
    PlaceAtmonZAxis(ang - 1, natoms - 1, XYZ);
    
    PlaceAtmonXZPlan(die - 1, natoms - 1, XYZ);
    
    R2 =  XYZ[j].x * XYZ[j].x + XYZ[j].y * XYZ[j].y + XYZ[j].z * XYZ[j].z;
    R = sqrt( R2 );
    
    angle = acos( XYZ[j].z / R );
    
    dihedral = atan2( XYZ[j].y , XYZ[j].x ); 
    
    
    if ( XYZ[ang - 1].z < ZERO ) angle = pi - angle;
    
    if ( XYZ[die - 1].x < ZERO && XYZ[die - 1].z > ZERO )  
      dihedral = - pi - dihedral;
    
    else if ( XYZ[die - 1].x < ZERO && XYZ[die - 1].z < ZERO )  
      dihedral = - pi + dihedral;
    
    else if ( XYZ[die - 1].x > ZERO && XYZ[die - 1].z > ZERO )  
      dihedral = - dihedral;
    
    else if ( XYZ[die - 1].x > ZERO && XYZ[die - 1].z < ZERO )  
      dihedral =  dihedral;
    
    /* To keep the angles between -pi and +pi */
    if ( dihedral > pi )  dihedral = 2*pi + dihedral;
    if ( dihedral < -pi )  dihedral = 2*pi + dihedral;
    
    zmat[j].dist = R;
    zmat[j].angle = angle*cv;
    zmat[j].dihedral = dihedral*cv;
  }
}
