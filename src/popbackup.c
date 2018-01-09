#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defmol.h"
#include "popbackup.h" 

int PopBackup (int natoms, int nindiv, int gen, int typga, double *energy, defmol *MPop, Txyz *XYZ){

  char  fname[15];
  int   i, j, ip;
  
  FILE *fl;
      
  switch(typga){
    
  case 1: /* Non-Inclusive Genetic Algorithm (NGA) */
    nindiv = nindiv;
    break;
    
  case 2: /* Inclusive Genetic Algorithm (IGA) */   
    nindiv = 3*nindiv;
    break;
  }

  strcpy(fname, "popbackup.dat");
  
  fl = fopen(fname,"w");

  fprintf(fl, "%d\n", gen);

  for ( i = 0; i < nindiv; i++ ) {
    ip = i*natoms;
    
    fprintf(fl, "%10.10lf\n", energy[i]);
    
    for(j = 0; j < natoms; j++){
      fprintf(fl, "\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[ip+j].x, XYZ[ip+j].y, XYZ[ip+j].z);
    }
  }
  fclose(fl);
}
