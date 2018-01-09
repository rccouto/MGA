#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gasdev.h"
#include "expdev.h"
#include "ran0.h"
#include "defmol.h"

#define ZERO          0.0E+0

int Mating(char const *Infile, int natoms, int nindiv, int matnum, int nonhydro, long *nseed, Txyz *XYZ, defmol *MPop, Tcnct *cnct, Tzmat *zmat)
{
  int     selecmat1, selecmat2, cutheight, centeratom1, centeratom2;
  int     j, i, m1, m2, atomonaxis1, atomonaxis2, atoms, indiv, ip;
  int     act, findiv;
  double  ran, ranlog, ranorm, xc, yc, zc, xc2, yc2, zc2, nonhydro2;
  Txyz    *p1XYZ, *p2XYZ, *p3XYZ;
  Tzmat   *p1zmat, *p2zmat, *p3zmat;

  /* Selecting molecules and the height of the cut */ 
  // indiv = nindiv/2;
 
  ran = ran0(nseed);
  selecmat1 = ceil(nindiv*ran) - 1;
  printf("selecmat1 = %d\n", selecmat1);
  
  ran = ran0(nseed);
  selecmat2 = ceil(nindiv*ran) - 1;
  printf("selecmat2 = %d\n", selecmat2);

  ran = ran0(nseed);
  cutheight = ceil(nonhydro*ran) - 1;
  printf("cutheight = %d\n", cutheight);
  
  // Molecule 1
  m1 = selecmat1*natoms;
  p1XYZ = &XYZ[m1];
  p1zmat = &zmat[m1];
  /*
  printf("17\n");
  printf("MOLECULA 1 \n");
  for (j = m1; j < (m1 + natoms); j++){
    printf("\t %s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
    }*/
  
  // Molecule 2
  m2 = selecmat2*natoms;
  p2XYZ = &XYZ[m2];
  p2zmat = &zmat[m2];

  /*
  printf("17\n");
  printf("MOLECULA 2 \n");
  for (j = m2; j < (m2 + natoms); j++){
    printf("\t %s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
    }*/
  
  findiv = (nindiv*natoms) + matnum*natoms;

  xyz2zmat(Infile, natoms, p1XYZ, cnct, p1zmat);

  xyz2zmat(Infile, natoms, p2XYZ, cnct, p2zmat);

  /* printf("17\n");
  printf("MOLECULA 2.2 \n");
  for (j = m2; j < (m2 + natoms); j++){
    printf("\t %s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
    }*/

  /* Z-Matrix of the New Molecule*/
  if ( cutheight == 1 ){
    zmat[findiv].dist = ZERO;
    zmat[findiv].angle = ZERO;
    zmat[findiv].dihedral = ZERO;
    
    zmat[findiv + 1].dist = zmat[m1 + 1].dist;
    zmat[findiv + 1].angle = ZERO;
    zmat[findiv + 1].dihedral = ZERO;
    
    zmat[findiv + 2].dist = zmat[m2 + 2].dist;
    zmat[findiv + 2].angle = zmat[m2 + 2].angle;
    zmat[findiv + 2].dihedral = ZERO;
    
    for ( j = 3; j < natoms; j++ ){
      zmat[findiv + j].dist = zmat[m2 + j].dist;
      zmat[findiv + j].angle = zmat[m2 + j].angle;
      zmat[findiv + j].dihedral = zmat[m2 + j].dihedral;
    }
  }

  else if ( cutheight == 2 ){
    zmat[findiv].dist = ZERO;
    zmat[findiv].angle = ZERO;
    zmat[findiv].dihedral = ZERO;
    
    zmat[findiv + 1].dist = zmat[m1 + 1].dist;
    zmat[findiv + 1].angle = ZERO;
    zmat[findiv + 1].dihedral = ZERO;
    
    zmat[findiv + 2].dist = zmat[m1 + 2].dist;
    zmat[findiv + 2].angle = zmat[m1 + 2].angle;
    zmat[findiv + 2].dihedral = ZERO;
    
    for ( j = 3; j < natoms; j++ ){
      zmat[findiv + j].dist = zmat[m2 + j].dist;
      zmat[findiv + j].angle = zmat[m2 + j].angle;
      zmat[findiv + j].dihedral = zmat[m2 + j].dihedral;
    }
  }
  
  else if ( cutheight > 2 ){
    zmat[findiv].dist = ZERO;
    zmat[findiv].angle = ZERO;
    zmat[findiv].dihedral = ZERO;
    
    zmat[findiv + 1].dist = zmat[m1 + 1].dist;
    zmat[findiv + 1].angle = ZERO;
    zmat[findiv + 1].dihedral = ZERO;
    
    zmat[findiv + 2].dist = zmat[m1 + 2].dist;
    zmat[findiv + 2].angle = zmat[m1 + 2].angle;
    zmat[findiv + 2].dihedral = ZERO;
    
    for ( j = 3; j < cutheight; j++ ){
      zmat[findiv + j].dist = zmat[m1 + j].dist;
      zmat[findiv + j].angle = zmat[m1 + j].angle;
      zmat[findiv + j].dihedral = zmat[m1 + j].dihedral;
    }
    
    for ( j = cutheight; j < natoms; j++ ){
      zmat[findiv + j].dist = zmat[m2 + j].dist;
      zmat[findiv + j].angle = zmat[m2 + j].angle;
      zmat[findiv + j].dihedral = zmat[m2 + j].dihedral;
    }
  }
  
  /* Cartesian Coordinates of the New Molecule */ 

  p3XYZ = &XYZ[findiv];
  p3zmat = &zmat[findiv];
  
  zmat2xyz(Infile, natoms, p3XYZ, cnct, p3zmat);

  /* Printing New Molecule */
  printf("17\n");
  printf("FINAL MOLECULE\n");
  for (i = 0; i < natoms; i++){
     printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[findiv + i].symb, XYZ[findiv + i].x, XYZ[findiv + i].y, XYZ[findiv + i].z);
  }
}
