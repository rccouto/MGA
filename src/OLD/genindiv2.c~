#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "genindiv.h"

#define ZERO          0.0E+0

extern distconvf;

int GenIndiv(char const *Infile, const char *potential, double const cf, const int natoms, defmol *indv)
{
  
  printf("\n  entrei 1 \n");

  int    ii, j, n, i, act, nseed;
  double distmin, distmax, ang, angmax, angmin, die, diemin, diemax, modvec, X, Y, Z, xc, yc, zc; 

  FILE   *fl0;
  char   rdaux[4], symb[2];
  Txyz   *pXYZ;
  
  printf("\n  entrei \n");
  scanf(yc);
 
  
  if( !strncasecmp(potential,"Ab-Initio", 9) ){

    if ( ! (fl0=fopen(Infile,"r")) )
      {
	printf("|\n\nError: The file \"%s\" is not available!\n");
	exit(EXIT_FAILURE);
      }

    fscanf(fl0, "%s", &rdaux);
  }
  else {
    
    if( !strncasecmp(potential,"tip", 3) ) n = 6;
    else n = 3;
    
    printf("#");
  }
	
}
