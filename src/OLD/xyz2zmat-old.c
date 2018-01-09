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

int xyz2zmat(char const *Infile, int const natoms, Txyz *XYZ, Tzmat *zmat)
{
  int    j, act, ang, die;
  double distmin, distmax, angmax, angmin, diemin, diemax;
  double X, Y, Z, xc, yc, zc, cv=180/pi; 
  double deltax, deltay, deltaz, R, R2;
  double angle, dihedral;

  FILE   *fl0;
  char   rdaux[4], symb[2];
  
  if ( ! (fl0=fopen(Infile,"r")) ){
    printf("|\n\nError: The file \"%s\" is not available!\n");
    exit(EXIT_FAILURE);
  }
  
  fscanf(fl0, "%s", &rdaux);
  
  while ( strncasecmp(rdaux,"*ZMAT", 5) ) fscanf(fl0, "%s \n", &rdaux);
  
  /* First Atom*/
  fscanf(fl0,"%s \n", &symb);
  
  zmat[0].dist = ZERO;
  zmat[0].angle = ZERO;
  zmat[0].dihedral = ZERO;
  
  /* Second Atom*/
  fscanf(fl0,"%s %d %lf %lf\n", &symb, &act, &distmin, &distmax);
  
  deltax = XYZ[1].x - XYZ[act - 1].x;
  deltay = XYZ[1].y - XYZ[act - 1].y;
  deltaz = XYZ[1].z - XYZ[act - 1].z;
  R2 = pow(deltax,2) +  pow(deltay,2) + pow(deltaz,2);
  R = sqrt( R2 );

  zmat[1].dist = R;
  zmat[1].angle = ZERO;
  zmat[1].dihedral = ZERO;
  
  /*Third Atom*/
  fscanf(fl0,"%s %d %lf %lf %d %lf %lf\n", &symb, &act, &distmin, &distmax, &ang, &angmin, &angmax);
 
  xc = XYZ[act - 1].x;
  yc = XYZ[act - 1].y;
  zc = XYZ[act - 1].z;
  
  centerxyz(natoms, xc, yc, zc, XYZ);
  
  PlaceAtmonZAxis(ang - 1, natoms, XYZ);
  
  R2 =  XYZ[2].x * XYZ[2].x + XYZ[2].y * XYZ[2].y + XYZ[2].z * XYZ[2].z;
  R = sqrt( R2 );
  
  angle = acos( XYZ[ 2].z / R );

  if ( XYZ[ang - 1].z < ZERO ) angle = pi - angle;
  
  zmat[2].dist = R;
  zmat[2].angle = angle*cv;
  zmat[2].dihedral = ZERO;
  
  for ( j = 3; j < natoms; j++ ){
    fscanf(fl0,"%s %d %lf %lf %d %lf %lf %d %lf %lf \n", 
	   &symb, &act, &distmin, &distmax, &ang, &angmin, &angmax, &die, &diemin, &diemax); 
    
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
  rewind(fl0);
  fclose(fl0);
}
