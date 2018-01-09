#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ran0.h"
#include "expdev.h"
#include "defmol.h"
#include "xyz2zmat.h"
#include "zmat2xyz.h"
#include "mutation.h"


#define ZERO          0.0E+0

int Mutation(int natoms, int nindiv, int mutnum, int nonhydro, int mutmin, int mutmax, long *nseed, Txyz *XYZ, defmol *MPop, Tcnct *cnct, Tzmat *zmat)
{
  int     j, i, ip;
  int     selecmut, mutpoint, mutfact;
  int     act, findiv, half, cut1, cut2, indiv;
  double  ran, ranlog, ranorm, xc, yc, zc, xc2, yc2, zc2, nonhydro2;
  Txyz    *iXYZ, *fXYZ;
  Tzmat   *izmat, *fzmat;


  /* Selecting molecules and the height of the cut */ 
  indiv = nindiv/2;
  ranlog = expdev(nseed);
  selecmut = ( ceil( indiv*ranlog ) )/3;
  //  printf("selecmut = %d\n", selecmut);
  
  ran = ran0(nseed);
  mutpoint = ceil((nonhydro-3)*ran)+2;
  //  printf("mutpoint = %d\n", mutpoint);
  
  // Molecule 1
  ip = selecmut*natoms;
  iXYZ = &XYZ[ip];
  izmat = &zmat[0];
  
  /*printf("17\n");
  printf("MOLECULE \n");
  for (j = ip; j < (ip + natoms); j++){
    printf("\t %s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
  }/**/

  xyz2zmat(natoms, iXYZ, cnct, izmat);

  /*  printf("\n\t\tMolecule Z Matrix\n");
  printf("\t%lf\n", zmat[ip+1].dist);
  printf("\t%lf    %lf\n", zmat[ip+2].dist, zmat[ip+2].angle);
  for(i=3; i < natoms; i++){
    printf("\t%lf    %lf  %12.6lf \n", zmat[ip+i].dist, zmat[ip+i].angle, zmat[ip+i].dihedral);
    }*/

  /* Generating New Molecule */
  zmat[natoms].dist = ZERO;
  zmat[natoms].angle = ZERO;
  zmat[natoms].dihedral = ZERO;
  
  zmat[natoms + 1].dist = zmat[ip + 1].dist;
  zmat[natoms + 1].angle = ZERO;
  zmat[natoms + 1].dihedral = ZERO;
  
  zmat[natoms + 2].dist = zmat[ip + 2].dist;
  zmat[natoms + 2].angle = zmat[ip + 2].angle;
  zmat[natoms + 2].dihedral = ZERO;
  
  for ( i = 3; i < natoms; i++ ){
    zmat[natoms + i].dist = zmat[ip + i].dist;
    zmat[natoms + i].angle = zmat[ip + i].angle;
    zmat[natoms + i].dihedral = zmat[ip + i].dihedral;
  }

  /* Mutating New Molecule */
  ran = ran0(nseed);
  mutfact = (mutmax - mutmin) * ran + mutmin;

  zmat[ natoms + mutpoint ].dihedral = zmat[ natoms + mutpoint ].dihedral + mutfact;

  /* Cartesian Coordinates of the New Molecule */ 
  findiv = (nindiv*natoms) + mutnum*natoms;
  //printf("findiv = %d\n", findiv);

  fXYZ = &XYZ[findiv];
  fzmat = &zmat[0];
  
  zmat2xyz(natoms, fXYZ, cnct, fzmat);
  
  /* Printing New Molecule */
  /* printf("17\n");
  printf("Mutation[%d]\n", mutnum);
  for (i = 0; i < natoms; i++){
    printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[findiv + i].symb, XYZ[findiv + i].x, XYZ[findiv + i].y, XYZ[findiv + i].z);
  }/**/

}
