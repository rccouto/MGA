#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "readbackup.h"

#define ZERO  0.0E+0

extern int gen;

int ReadBackup(int natoms, int nindiv, int typga, double *energy, Txyz *XYZ, defmol *indv)
{
  int    i, j, ip, ngen;
  char   symb[2], fname[14];
  double coordx, coordy, coordz, nenergy;
  
  FILE *fl;

  nindiv=nindiv/2;

  strcpy(fname, "popbackup.dat");

  if ( ! (fl=fopen(fname,"r")) ){
    printf("|\n\nError: The file \"%s\" is not available!\n", fname);
    exit(EXIT_FAILURE);
  }
  
  switch(typga){
    
  case 1: /* Non-Inclusive Genetic Algorithm (NGA) */
    nindiv = 2*nindiv;
    break;
    
  case 2: /* Inclusive Genetic Algorithm (IGA) */   
    nindiv = 3*nindiv;
    break;
  }
   
  //printf("\n\nn=%d\n\n",nindiv);

  fscanf(fl,"%d", &ngen);
  gen = ngen;
  
  for( i = 0; i < nindiv; i++){
    
    ip = i*natoms;
    
    fscanf(fl,"%lf", &nenergy);
    energy[i] = nenergy;
    printf("%lf\n", energy[i]);
	  
    for(j = ip; j < (i + 1)*natoms; j++ ){
      fscanf(fl,"%s %lf %lf %lf", &symb, &coordx, &coordy, &coordz);

      strcpy(indv[j].symb, symb);
      indv[j].np = tableZ(symb);
      indv[j].mass = tableM(symb);

      XYZ[j].x = coordx;
      XYZ[j].y = coordy;
      XYZ[j].z = coordz;
    }
  }
  rewind(fl);
  fclose(fl);
}

