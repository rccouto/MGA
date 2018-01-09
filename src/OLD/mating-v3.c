#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gasdev.h"
#include "expdev.h"
#include "ran0.h"
#include "defmol.h"
#include "standmating.h"
#include "ordermolx.h"

int Mating(char const *Infile, int natoms, int nindiv, long *nseed, Txyz *XYZ, defmol *MPop)
{
  int     selecmat1, selecmat2, cutheight, centeratom1, centeratom2;
  int     j, i, i1, i2, atomonaxis1, atomonaxis2, atoms, indiv, ip;
  int     nonhydro, act;
  double  ran, ranlog, ranorm, xc, yc, zc, xc2, yc2, zc2, nonhydro2;
  Txyz    *p1XYZ, *p2XYZ;
  defmol  *pMPop, *p1MPop, *p2MPop;
  
  FILE   *fl0;
  char   rdaux[4], symb[2];   

  /* Counting non-hydrogen atoms */
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
  
  /* printf("\nMOLECULA 1 NÃO PADRONIZADA\n");
  for (j = i1; j < (i1 + natoms); j++){
    printf("\t%d %s %15.8lf  %15.8lf  %15.8lf\n", j, &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
    }*/
  
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
  
  /*  printf("\nMOLECULA 2 NÃO PADRONIZADA\n");
  for (j = i2; j < (i2 + natoms); j++){
    printf("\t%d %s %15.8lf  %15.8lf  %15.8lf\n", j, &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
    }*/
  
  centeratom2 = cutheight;
  atomonaxis2 = cutheight + 1;
  StandMating(natoms, centeratom2, atomonaxis2, p2XYZ); 
  
  printf("\nMOLECULA 2 PADRONIZADA\n");
  for (j = i2; j < (i2 + natoms); j++){
    printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[j].x, XYZ[j].y, XYZ[j].z);
  }
 
  
  /* Generating new molecule */
  cutheight=cutheight+1;
// The Molecular Skeleton
  printf("FIRST CUT\n");
  for (i = 0; i < cutheight; i++){
    *MPop[(nindiv*natoms) + 1 + i].symb = *MPop[i1 + i].symb;
      XYZ[(nindiv*natoms) + 1 + i].x = XYZ[i1 + i].x;
      XYZ[(nindiv*natoms) + 1 + i].y = XYZ[i1 + i].y; 
      XYZ[(nindiv*natoms) + 1 + i].z = XYZ[i1 + i].z;
    
      printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[(nindiv*natoms) +1 + i].symb, XYZ[(nindiv*natoms) + 1 + i].x, XYZ[(nindiv*natoms) + 1 + i].y, XYZ[(nindiv*natoms) + 1 + i].z);
  }
  printf("SECOND CUT\n");
  for (i = 0; i < (nonhydro-cutheight); i++){
    *MPop[(nindiv*natoms) + 1 + cutheight + i].symb = *MPop[i2 + i + cutheight].symb;
      XYZ[(nindiv*natoms) + 1 + cutheight + i].x = XYZ[i2 + cutheight + i].x;
      XYZ[(nindiv*natoms) + 1 + cutheight + i].y = XYZ[i2 + cutheight + i].y; 
      XYZ[(nindiv*natoms) + 1 + cutheight + i].z = XYZ[i2 + cutheight + i].z;
    
      printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[(nindiv*natoms) + 1 + i + cutheight].symb, XYZ[(nindiv*natoms) + 1 + i + cutheight].x, XYZ[(nindiv*natoms) + 1 + i + cutheight].y, XYZ[(nindiv*natoms) + 1 + i + cutheight].z);
  }

// The Hydrogens

  // Reading Input File
  if ( ! (fl0=fopen(Infile,"r")) ){
    printf("|\n\nError: The file \"%s\" is not available!\n");
    exit(EXIT_FAILURE);
  }
  fscanf(fl0, "%s", &rdaux);
  while ( strncasecmp(rdaux,"*ZMAT", 5) ) fscanf(fl0, "%s \n", &rdaux);
  fscanf(fl0,"%s \n", &symb);
  fscanf(fl0,"%s %d %lf %lf\n", &symb, &act, &rdaux, &rdaux);
  fscanf(fl0,"%s %d %lf %lf %d %lf %lf \n", &symb, &act, &rdaux, &rdaux, &rdaux, &rdaux, &rdaux);
  for( i = 3; i < natoms; i++){  
    fscanf(fl0,"%s %d %lf %lf %d %lf %lf %d %lf %lf\n", &symb, &act, &rdaux, &rdaux, &rdaux, &rdaux, &rdaux, &rdaux, &rdaux, &rdaux);   

    if ( i >= nonhydro){
      if (act < cutheight){
	*MPop[(nindiv*natoms) + 1 + i].symb = *MPop[i1 + i].symb;
	 XYZ[(nindiv*natoms) + 1 + i].x = XYZ[i1 + i].x;
	 XYZ[(nindiv*natoms) + 1 + i].y = XYZ[i1 + i].y; 
	 XYZ[(nindiv*natoms) + 1 + i].z = XYZ[i1 + i].z;
	 
	 //	 printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[(nindiv*natoms) +1 + i].symb, XYZ[(nindiv*natoms) + 1 + i].x, XYZ[(nindiv*natoms) + 1 + i].y, XYZ[(nindiv*natoms) + 1 + i].z);
      }
      else{
	*MPop[(nindiv*natoms) + 1 + i].symb = *MPop[i2 + i ].symb;
	  XYZ[(nindiv*natoms) + 1 + i].x = XYZ[i2 + i].x;
	  XYZ[(nindiv*natoms) + 1 + cutheight + i].y = XYZ[i2 + cutheight + i].y; 
	  XYZ[(nindiv*natoms) + 1 + cutheight + i].z = XYZ[i2 + cutheight + i].z;
      }
    }
  }
  
  


  /* Printing New Molecule */
  printf("\nFINAL MOLECULE\n");
  for (i = 0; i < natoms; i++){
     printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[(nindiv*natoms) + 1 + i].symb, XYZ[(nindiv*natoms) + 1 + i].x, XYZ[(nindiv*natoms) + 1 + i].y, XYZ[(nindiv*natoms) + 1 + i].z);
  }
}
