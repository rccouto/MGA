#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "defmol.h"
#include "readconnect.h"


#define ZERO  0.0E+0

int ReadConnections(char const *Infile, int const natoms, Tcnct *cnct)
{
  int    i, act, ang, die;
  
  FILE   *fl0;
  char   rdaux[4], symb[2];   
  
  if ( ! (fl0=fopen(Infile,"r")) ){
    printf("|\n\nError: The file \"%s\" is not available!\n");
    exit(EXIT_FAILURE);
  }

  fscanf(fl0, "%s", &rdaux);
  printf("%s\n", rdaux);

  while ( strncasecmp(rdaux,"*ZMAT", 5) ) fscanf(fl0, "%s \n", &rdaux);
  
  fscanf(fl0,"%s \n", &symb);      

  printf("OI 2\n");


  printf("OI 3\n");
  cnct[0].cdist     = ZERO;
  cnct[0].cangle    = ZERO;
  cnct[0].cdihedral = ZERO;

  fscanf(fl0,"%s %d %lf %lf\n", &symb, &act, &rdaux, &rdaux);
  cnct[1].cdist     = act;
  cnct[1].cangle    = ZERO;
  cnct[1].cdihedral = ZERO;
  
  fscanf(fl0,"%s %d %lf %lf %d %lf %lf \n", &symb, &act, &rdaux, &rdaux, &ang, &rdaux, &rdaux);
  cnct[2].cdist     = act;
  cnct[2].cangle    = ang;
  cnct[2].cdihedral = ZERO;
  
  for( i = 3; i < natoms; i++){
    fscanf(fl0,"%s %d %lf %lf %d %lf %lf %d %lf %lf\n", &symb, &act, &rdaux, &rdaux, &ang, &rdaux, &rdaux, &die, &rdaux, &rdaux);
    cnct[3].cdist     = act;
    cnct[3].cangle    = ang;
    cnct[3].cdihedral = die;
    }
  rewind(fl0);
  fclose(fl0);
  //  for (i = 1; i < (3*natoms) - 5; i++) printf("cnct = %d\n", cnct[i]);
}
