#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ran0.h"
#include "expdev.h"
#include "defmol.h"
#include "mating.h"
#include "xyz2zmat.h"
#include "zmat2xyz.h"
#include "matsc.h"
#include "matdc.h"

#define ZERO          0.0E+0

int Mating(int typmat, int natoms, int nindiv, int matnum, int nonhydro, long *nseed, Txyz *XYZ, defmol *MPop, Tcnct *cnct, Tzmat *zmat)

{
  int     selecmat1=0, selecmat2=0, cutheight, centeratom1, centeratom2;
  int     j, i, m1, m2, atomonaxis1, atomonaxis2, atoms, indiv, ip;
  int     act, findiv, half, cut1, cut2;
  double  ran, ranlog, ranorm, xc, yc, zc, xc2, yc2, zc2, nonhydro2;
  Txyz    *p1XYZ, *p2XYZ, *p3XYZ;
  Tzmat   *p1zmat, *p2zmat, *p3zmat;


  /* Selecting molecules and the height of the cut */ 
  ranlog = expdev(nseed);
  selecmat1 = ceil( nindiv*ranlog ) - 1;
  //printf("%lf selecmat1 = %d\n", ranlog, selecmat1);
  
  do{
    ranlog = expdev(nseed);
    selecmat2 = ceil( nindiv*ranlog ) - 1;
    
    //printf("ranlog=%lf, s=%d\n", ranlog, selecmat2);
  }while(selecmat2 == selecmat1);
  //printf("selecmat2 = %d\n", selecmat2);
  
  ran = ran0(nseed);
  cutheight = ceil((nonhydro-3)*ran)+2;
  //printf("cutheight = %d\n", cutheight);
  
  // Molecule 1
  m1 = selecmat1*natoms;
  p1XYZ = &XYZ[m1];
  p1zmat = &zmat[0];
  
  /*printf("17\n");
  printf("MOLECULA 1 \n");
  for (j = m1; j < (m1 + natoms); j++){
    printf("\t %s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
  }/**/
  
  // Molecule 2
  m2 = selecmat2*natoms;
  p2XYZ = &XYZ[m2];
  p2zmat = &zmat[natoms];
  
  /*printf("17\n");
  printf("MOLECULA 2 \n");
  for (j = m2; j < (m2 + natoms); j++){
    printf("\t %s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
  }/**/
  
  xyz2zmat(natoms, p1XYZ, cnct, p1zmat);
  
  xyz2zmat(natoms, p2XYZ, cnct, p2zmat);
  
  /*printf("\n\t\tMolecule 1\n");
  printf("\t%lf\n", zmat[1].dist);
  printf("\t%lf    %lf\n", zmat[2].dist, zmat[2].angle);
  for(i=3; i < natoms; i++){
    printf("\t%lf    %lf  %12.6lf \n", zmat[i].dist, zmat[i].angle, zmat[i].dihedral);
   }/**/

  /* printf("\n\t\tMolecule 2\n");
  printf("\t%lf\n", zmat[natoms+1].dist);
  printf("\t%lf    %lf\n", zmat[natoms+2].dist, zmat[natoms+2].angle);
  for(i=3; i < natoms; i++){
    printf("\t%lf    %lf  %12.6lf \n", zmat[natoms+i].dist, zmat[natoms+i].angle, zmat[natoms+i].dihedral);
   }/**/

  switch(typmat){
    
  case 1: /* Single Cut */
    MatSC(natoms, cutheight, zmat);
    break;
    
  case 2: /* Double Cut */
    half = nonhydro/2;
    ran = ran0(nseed);
    cut1 = ceil( half * ran ) - 1;
    ran = ran0(nseed);
    cut2 = ceil( (nonhydro - half - 1) * ran ) + half;
    
    MatDC(natoms, cut1, cut2, zmat);
    break;
    
  case 3: /* Mixed Cut */
    if (matnum%2==0) MatSC(natoms, cutheight, zmat);
    else{
      half = nonhydro/2;
      ran = ran0(nseed);
      cut1 = ceil( half * ran ) - 1;
      ran = ran0(nseed);
      cut2 = ceil( (nonhydro - half - 1) * ran ) + half;
      MatDC(natoms, cut1, cut2, zmat);
    }
    break;
  }
   
  /* Cartesian Coordinates of the New Molecule */ 
  findiv = (nindiv*natoms) + matnum*natoms;
  //  printf("[mating] new indiv = %d\n", findiv);

  p3XYZ = &XYZ[findiv];
  p3zmat = &zmat[2*natoms];
  
  zmat2xyz(natoms, p3XYZ, cnct, p3zmat);
  
  /*printf("\n\t\tMolécula Final\n");
  printf("\t%lf\n", zmat[2*natoms+1].dist);
  printf("\t%lf    %lf\n", zmat[2*natoms+2].dist, zmat[2*natoms+2].angle);
  for(i=3; i < natoms; i++){
    printf("\t%lf    %lf  %12.6lf \n", zmat[2*natoms+i].dist, zmat[2*natoms+i].angle, zmat[2*natoms+i].dihedral);
  }/**/

  /* Printing New Molecule */
  /*printf("17\n");
  printf("Mating[%d]\n", matnum);
  for (i = 0; i < natoms; i++){
     printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[findiv + i].symb, XYZ[findiv + i].x, XYZ[findiv + i].y, XYZ[findiv + i].z);
  }/**/
}
