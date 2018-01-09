#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "defmol.h"
#include "genindv.h"
#include "placeatomonzaxis.h"
#include "placeatomonplan.h"
#include "centerxyz.h"

#include "ran0.h"

#define ZERO          0.0E+0

extern double pi;

int GenIndv(const char flag, char const *Infile, const double distconvf, const int natoms, long *nseed, Txyz *XYZ, defmol *indv)
{
  int    j, n, i, act, ang, die;
  double distmin, distmax, angmax, angmin, diemin, diemax, modvec;
  double X, Y, Z, xc, yc, zc, ran, zap, zp; 
  double status, status2, cf = pi/180E+0;
  double ang1, angang, angdie, sdie, cang, sang, cdie;

  FILE   *fl0;
  char   rdaux[4], symb[2];

  
  if ( ! (fl0=fopen(Infile,"r")) ){
    printf("|\n\nError: The file \"%s\" is not available!\n");
    exit(EXIT_FAILURE);
  }
  
  fscanf(fl0, "%s", &rdaux);
  
  while ( strncasecmp(rdaux,"*ZMAT", 5) ) fscanf(fl0, "%s \n", &rdaux);
  
  fscanf(fl0,"%s \n", &symb);      
  
  /* First Atom*/
  indv[0].atmn = 1;
  strcpy(indv[0].symb, symb);
  XYZ[0].x = ZERO;
  XYZ[0].y = ZERO;
  XYZ[0].z = ZERO;
  indv[0].xyz = &XYZ[0];
  indv[0].np = tableZ(symb);
  indv[0].mass = tableM(symb);
  
  fscanf(fl0,"%s %d %lf %lf\n", &symb, &act, &distmin, &distmax);
  
  if ( flag == 'A' ){
    distmin = (distmin + distmax)/2;
    distmax = distmin;
  }
  
  /* Second Atom*/
  X = ZERO;
  Y = ZERO;
  
  ran = ran0(nseed);
  Z = (distmax - distmin) * ran + distmin;  
  
  strcpy(indv[1].symb, symb);
  XYZ[1].x = X * distconvf;
  XYZ[1].y = Y * distconvf;
  XYZ[1].z = Z * distconvf;
  indv[1].xyz = &XYZ[1];
  indv[1].np = tableZ(symb);
  indv[1].mass = tableM(symb);
  
  fscanf(fl0,"%s %d %lf %lf %d %lf %lf \n", &symb, &act, &distmin, &distmax, &ang, &angmin, &angmax);
  
  if ( flag == 'A' ){
    distmin = (distmin + distmax) / 2;
    distmax =  distmin;
    angmin  = (angmin + angmax) / 2;
    angmax  =  angmin;
  }
  
  /* Compute The Vector Module */
  ran = ran0(nseed);
  modvec = (distmax - distmin) * ran + distmin;
  
  /*Third Atom*/
  Y = ZERO + XYZ[act - 1].y;
  ran = ran0(nseed);
  zap = XYZ[ang - 1].z - XYZ[act - 1].z;
  
  if ( zap < ZERO ) ang1 = pi - ((angmax - angmin) * ran + angmin) * cf ;
  else              ang1 =      ((angmax - angmin) * ran + angmin) * cf ;
  zp = modvec * cos( ang1 ); 
  Z = zp + XYZ[act - 1].z; 
  X = sqrt( modvec * modvec - zp * zp ) + XYZ[act - 1].x;
  
  strcpy(indv[2].symb, symb);
  XYZ[2].x = X * distconvf;
  XYZ[2].y = Y * distconvf;
  XYZ[2].z = Z * distconvf;
  indv[2].xyz = &XYZ[2];
  indv[2].np = tableZ(symb);
  indv[2].mass = tableM(symb);
  
  for ( j = 3; j < natoms; j++ ) {
    fscanf(fl0,"%s %d %lf %lf %d %lf %lf %d %lf %lf \n", 
	   &symb, &act, &distmin, &distmax, &ang, &angmin, &angmax, &die, &diemin, &diemax);
    
    if ( flag == 'A' ){
      distmin = (distmin + distmax) / 2;
      distmax =  distmin;
      angmin  = (angmin + angmax) / 2;
      angmax  =  angmin;
      diemin  = (diemin + diemax) / 2;
      diemax  =  diemin;
    }
    
    /* degrees to radians*/
    angmin = angmin * cf;
    angmax = angmax * cf;
    diemin = diemin * cf;
    diemax = diemax * cf;
    
    xc = XYZ[act - 1].x; 
    yc = XYZ[act - 1].y;
    zc = XYZ[act - 1].z;
    
    centerxyz(j, xc, yc, zc, XYZ);
    
    PlaceAtmonZAxis(ang - 1, j-1, XYZ);
    
    PlaceAtmonXZPlan(die - 1, j-1, XYZ);
    
    /* Compute The Vector Module */
    ran = ran0(nseed);
    modvec = (distmax - distmin) * ran + distmin;
    
    /* Compute The Theta Angle */
    ran = ran0(nseed);
    if ( XYZ[ang - 1].z < ZERO ) angang = pi + ((angmax - angmin) * ran + angmin);
    else                         angang =       (angmax - angmin) * ran + angmin;
    
    /* Compute The Phi Angle */
    ran = ran0(nseed);
    if ( XYZ[die - 1].x < ZERO && XYZ[die - 1].z > ZERO )  
      angdie = pi - ((diemax - diemin) * ran + diemin);
    
    else if ( XYZ[die - 1].x < ZERO && XYZ[die - 1].z < ZERO )  
      angdie = 2*pi + ((diemax - diemin) * ran + diemin);
    
    else if ( XYZ[die - 1].x > ZERO && XYZ[die - 1].z > ZERO )  
      angdie = - ((diemax - diemin) * ran + diemin);
    
    else if ( XYZ[die - 1].x > ZERO && XYZ[die - 1].z < ZERO )  
      angdie = pi + ((diemax - diemin) * ran + diemin);
    
    cang = cos(angang);
    sang = sin(angang);
    cdie = cos(angdie);
    sdie = sin(angdie);
    
    X = modvec * sang * cdie;
    Y = modvec * sang * sdie;
    Z = modvec * cang;
    
    strcpy(indv[j].symb, symb);
    XYZ[j].x = X;
    XYZ[j].y = Y;
    XYZ[j].z = Z;
    indv[j].xyz = &XYZ[j];
    indv[j].np = tableZ(symb);
    indv[j].mass = tableM(symb);
  }
  rewind(fl0);
  fclose(fl0);
}

