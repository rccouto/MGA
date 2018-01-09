#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gasdev.h"
#include "expdev.h"
#include "ran0.h"
#include "defmol.h"
#include "standmating.h"
#include "ordermolx.h"

int Mating(int natoms, int nindiv, long *nseed, Txyz *XYZ, defmol *MPop)
{
  int     selecmat1, selecmat2, cutheight, centeratom1, centeratom2;
  int     j, i, i1, i2, atomonaxis1, atomonaxis2, atoms, indiv, ip;
  int     nonhydro;
  double  ran, ranlog, ranorm, xc, yc, zc, xc2, yc2, zc2, nonhydro2;
  Txyz    *p1XYZ, *p2XYZ;
  defmol  *pMPop, *p1MPop, *p2MPop;
  
  /* Counting Non-hydrogen atoms */
  nonhydro = 0;
  for ( i = 0; i < natoms; i++){
    //   printf("symb = %s\n", MPop[i].symb); 
    if ( *MPop[i].symb != 'H' ){
      nonhydro = nonhydro + 1;
    }
  }
  printf("nonhydro=%d\n", nonhydro);
  /* Selecting molecules and the height of the cut */ 

  indiv = nindiv/2;
  ranlog = expdev(nseed);
  selecmat1 = (round( indiv*ranlog))/3;
  printf("nindiv = %d\n", nindiv);
  printf("selecmat1 = %d\n", selecmat1);
  
  ranlog = expdev(nseed);
  selecmat2 = (round( indiv*ranlog))/3;
  printf("selecmat2 = %d\n", selecmat2);
  
  nonhydro2=nonhydro/2;
  ranorm = gasdev(nseed)/3;
  cutheight = (nonhydro2*ranorm)+nonhydro2;
  
  printf("Cutheight = %d\n", cutheight);
  
  /* Standardizing molecules */
 
  // Molecule 1
   i1 = selecmat1*natoms;
   p1XYZ = &XYZ[i1];
   
   printf("\nMOLECULA 1 NÃO PADRONIZADA\n");
   for (j = i1; j < (i1 + natoms); j++){
     printf("\t%d %s %15.8lf  %15.8lf  %15.8lf\n", j, &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
   }
   
   centeratom1 = cutheight;
   atomonaxis1 = cutheight + 1;
   printf("i1 = %d\n", i1);
   StandMating(natoms, centeratom1, atomonaxis1, p1XYZ); 
 
   printf("\nMOLECULA 1 PADRONIZADA\n");
  for (j = i1; j < (i1 + natoms); j++){
    printf("\t%d %s %15.8lf  %15.8lf  %15.8lf\n", j, &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
    }
  
  // Molecule 2
  i2 = selecmat2*natoms;
  p2XYZ = &XYZ[i2];
  
  printf("\nMOLECULA 2 NÃO PADRONIZADA\n");
  for (j = i2; j < (i2 + natoms); j++){
    printf("\t%d %s %15.8lf  %15.8lf  %15.8lf\n", j, &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
  }
  
  centeratom2 = cutheight;
  atomonaxis2 = cutheight + 1;
  StandMating(natoms, centeratom2, atomonaxis2, p2XYZ); 
  

  printf("\nMOLECULA 2 PADRONIZADA\n");
  for (j = i2; j < (i2 + natoms); j++){
    printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
    }
  
  /* Ordering Atoms according to X axis */

  OrderMolX(natoms, i1, i2, MPop);

  printf("\nMOLECULA 1 REORDENADA\n");
  for (j = i1; j < (i1 + natoms); j++){
    printf("\t %s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
  }
  printf("\nMOLECULA 2 REORDENADA\n");
  for (j = i2; j < (i2 + natoms); j++){
    printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
  }

  /* Generating new molecule */
  printf("FIRST CUT\n");
  for (i = 0; i < cutheight; i++){
   *MPop[(nindiv*natoms) + 1 + i].symb = *MPop[i1 + i].symb;
    XYZ[(nindiv*natoms) + 1 + i].x = XYZ[i1 + i].x;
    XYZ[(nindiv*natoms) + 1 + i].y = XYZ[i1 + i].y; 
    XYZ[(nindiv*natoms) + 1 + i].z = XYZ[i1 + i].z;

    printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[(nindiv*natoms) +1 + i].symb, XYZ[(nindiv*natoms) + 1 + i].x, XYZ[(nindiv*natoms) + 1 + i].y, XYZ[(nindiv*natoms) + 1 + i].z);
  }
  printf("SECOND CUT\n");
  for (i = 0; i < (natoms-cutheight); i++){
    *MPop[(nindiv*natoms) + 1 + i + cutheight].symb = *MPop[i2 + i + cutheight].symb;
    XYZ[(nindiv*natoms) + 1 + i + cutheight].x = XYZ[i2 + i + cutheight].x;
    XYZ[(nindiv*natoms) + 1 + i + cutheight].y = XYZ[i2 + i + cutheight].y; 
    XYZ[(nindiv*natoms) + 1 + i + cutheight].z = XYZ[i2 + i + cutheight].z;
 
    printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[(nindiv*natoms) + 1 + i + cutheight].symb, XYZ[(nindiv*natoms) + 1 + i + cutheight].x, XYZ[(nindiv*natoms) + 1 + i + cutheight].y, XYZ[(nindiv*natoms) + 1 + i + cutheight].z);
  }
  
  /* Printing New Molecule */
  printf("\nFINAL MOLECULE\n");
  for (i = 0; i < natoms; i++){
    printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[(nindiv*natoms) + 1 + i].symb, XYZ[(nindiv*natoms) + 1 + i].x, XYZ[(nindiv*natoms) + 1 + i].y, XYZ[(nindiv*natoms) + 1 + i].z);
  }

}  


  
