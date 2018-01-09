#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "ran0.h"
#include "expdev.h"
#include "centerxyz.h"
#include "placeatomonzaxis.h"
#include "rotate.h"
#include "mutation.h"

extern double pi;
extern defmol *MPop;

int Mutation(int natoms, int nindiv, int nrotgp, int mutnum, int mutmin, int mutmax, long *nseed, Trtt *rtt, Tatom *rttatm, Txyz *XYZ)
{
  int     j, i, ip;
  int     selecmut, mutpoint, count;
  int     findiv, indiv, rot, plane, center;
  double  ran, ranlog;
  double  A, XT, YT, xc, yc, zc;
  Txyz    *fXYZ;
  
  
  /* Selecting molecule */ 
  ranlog = expdev(nseed);
  selecmut = ceil( nindiv*ranlog )-1;
  //printf("selecmut = %d\n", selecmut);
  
  /* Selecting Mutation Point */
  ran = ran0(nseed);
  mutpoint = ( ceil(nrotgp*ran) )-1;
  //printf("mutpoint = %d\n", mutpoint);
  

  // Molecule
  ip = selecmut*natoms;
  
  /*printf("17\n");
  printf("Molécula Selecionada \n");
  for (j = ip; j < (ip + natoms); j++){
    printf("\t %s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
  }/**/
  
  
  /* Generating New Molecule */
  findiv = (nindiv*natoms) + mutnum*natoms;

  for( i = 0; i < natoms; i++){
    XYZ[findiv + i].x = XYZ[ip + i].x;
    XYZ[findiv + i].y = XYZ[ip + i].y;
    XYZ[findiv + i].z = XYZ[ip + i].z;
  }
 
  /* Mutating New Molecule */
  fXYZ = &XYZ[findiv];

  count = 0;
  for( i = 0; i < nrotgp; i++ ){
    if( i == mutpoint){ 
      center = findiv + rtt[i].center;

      xc = XYZ[center - 1].x;
      yc = XYZ[center - 1].y;
      zc = XYZ[center - 1].z;
      
      centerxyz(natoms, xc, yc, zc, fXYZ);

      PlaceAtmonZAxis(rtt[i].plane - 1, natoms - 1, fXYZ);

      ran = ran0(nseed);
      A = (mutmax - mutmin)*ran + mutmin;
      A = (A*pi)/180;
    }
    for( j = 0; j < rtt[i].natmrot; j++ ){
      if( i == mutpoint){ 
	rot = (rttatm[count].atom - 1);
	rotate_z( XYZ[findiv + rot].x, XYZ[findiv + rot].y, A, &XT, &YT);
	XYZ[findiv + rot].x = XT;
	XYZ[findiv + rot].y = YT;
      }
      count=count+1;
    }
  }
  
  
  /* Cartesian Coordinates of the New Molecule */ 

  /* Printing New Molecule */
  /*printf("17\n");
  printf("Mutation[%d] - Final\n", mutnum);
  for (i = 0; i < natoms; i++){
    printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[findiv + i].symb, XYZ[findiv + i].x, XYZ[findiv + i].y, XYZ[findiv + i].z);
  }/**/
}
