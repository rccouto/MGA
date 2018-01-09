#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "defmol.h"
#include "zmat2xyz.h"
#include "placeatomonzaxis.h"
#include "placeatomonplan.h"
#include "centerxyz.h"

#define ZERO          0.0E+0

extern double pi;

int zmat2xyz(int const natoms, Txyz *XYZ, Tcnct *cnct, Tzmat *zmat)
{
  int    j, act, ang, die;
  double X, Y, Z, xc, yc, zc, zap, zp, modvec, cv=pi/180; 
  double ang1, angang, angdie, sdie, cang, sang, cdie;
  
  /* First Atom */
  XYZ[0].x = ZERO;
  XYZ[0].y = ZERO;
  XYZ[0].z = ZERO;
  
  /* Second Atom */
  XYZ[1].x = ZERO;
  XYZ[1].y = ZERO;
  XYZ[1].z = zmat[1].dist;

  /* Third Atom */
  
  // Compute Connections and Distance
  act    = cnct[2].cdist;
  ang    = cnct[2].cangle;
  modvec = zmat[2].dist;

  Y = ZERO + XYZ[act - 1].y;

  zap = XYZ[ang - 1].z - XYZ[act - 1].z;
  if ( zap < ZERO ) ang1 = pi - zmat[2].angle * cv ;
  else              ang1 =      zmat[2].angle * cv ;

  zp = modvec * cos( ang1 ); 

  Z = zp + XYZ[act - 1].z; 
  X = sqrt( modvec * modvec - zp * zp ) + XYZ[act - 1].x;
  
  XYZ[2].x = X;
  XYZ[2].y = Y;
  XYZ[2].z = Z;
  
  for ( j = 3; j < natoms; j++ ) {
    
    // Compute Connections
    act = cnct[j].cdist;
    ang = cnct[j].cangle;
    die = cnct[j].cdihedral;
  
    xc = XYZ[act - 1].x; 
    yc = XYZ[act - 1].y;
    zc = XYZ[act - 1].z;
    
    centerxyz(j, xc, yc, zc, XYZ);
    
    PlaceAtmonZAxis(ang - 1, j-1, XYZ);
    
    PlaceAtmonXZPlan(die - 1, j-1, XYZ);
    
    /* Compute The Vector Module */
    modvec = zmat[j].dist;
    
    /* Compute The Theta Angle */

    if ( XYZ[ang - 1].z < ZERO ) angang = pi + zmat[j].angle * cv;
    else                         angang =      zmat[j].angle * cv;
    
    /* Compute The Phi Angle */
    if ( XYZ[die - 1].x < ZERO && XYZ[die - 1].z > ZERO )  
      angdie = pi - zmat[j].dihedral * cv;
    
    else if ( XYZ[die - 1].x < ZERO && XYZ[die - 1].z < ZERO )  
      angdie = 2*pi + zmat[j].dihedral * cv;
    
    else if ( XYZ[die - 1].x > ZERO && XYZ[die - 1].z > ZERO )  
      angdie = - zmat[j].dihedral * cv;
    
    else if ( XYZ[die - 1].x > ZERO && XYZ[die - 1].z < ZERO )  
      angdie = pi + zmat[j].dihedral * cv;
    
    cang = cos(angang);
    sang = sin(angang);
    cdie = cos(angdie);
    sdie = sin(angdie);
    
    X = modvec * sang * cdie;
    Y = modvec * sang * sdie;
    Z = modvec * cang;
    
    XYZ[j].x = X;
    XYZ[j].y = Y;
    XYZ[j].z = Z;
  }
}
