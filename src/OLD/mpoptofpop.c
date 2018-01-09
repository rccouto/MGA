#include <stdio.h>
#include <stdlib.h>

#include "defmol.h"
#include "mpoptofpop.h"

void MPopToFPop(int natoms, int nindiv, double *energy, double *fenergy,defmol *indiv, defmol *Findiv)
{
  int    i, j;
  int    nj;

  for ( i = 0 ; i < nindiv ; i++ ){
    fenergy[i]    = energy[i];
    Findiv[i].xyz = indiv[i].xyz;
  }
  printf("\n\n");
  for (j=0; j<100;j++){
    printf("\t%15.8lf  %15.8lf  %15.8lf\n", Findiv[j].xyz->x, Findiv[j].xyz->y, Findiv[j].xyz->z);
  }
  
}
