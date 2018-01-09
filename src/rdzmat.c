#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "defmol.h"
#include "rdzmat.h"

#define ZERO          0.0E+0

extern double pi;

int RdZMat(char const *Infile, int const natoms, Tcnct *cnct)
{
  int    i, act, ang, die;
  double work;

  FILE   *fl0;
  char   rdaux[200], symb[2];
  double dlixo;
  
  if ( ! (fl0=fopen(Infile,"r")) ){
    printf("|\n\nError: The file \"%s\" is not available!\n");
    exit(EXIT_FAILURE);
  }
  
  fscanf(fl0, "%s", &rdaux);
  
  while ( strncasecmp(rdaux,"*ZMAT", 5) ) fscanf(fl0, "%s \n", &rdaux);
  
  fscanf(fl0,"%s \n", &symb);
  
  /* First Atom*/
  cnct[0].cdist     = ZERO;
  cnct[0].cangle    = ZERO;
  cnct[0].cdihedral = ZERO;
  
  /* Second Atom*/
  fscanf(fl0,"%s %d %lf %lf\n", &symb, &act, &dlixo, &dlixo);
  
  cnct[1].cdist     = act;
  cnct[1].cangle    = ZERO;
  cnct[1].cdihedral = ZERO;
  
  /*Third Atom*/
  fscanf(fl0,"%s %d %lf %lf %d %lf %lf\n", &symb, &act, &dlixo, &dlixo, &ang, &dlixo, &dlixo);
  
  cnct[2].cdist     = act;
  cnct[2].cangle    = ang;
  cnct[2].cdihedral = ZERO;
  
  for ( i = 3; i < natoms; i++ ){
    fscanf(fl0,"%s %d %lf %lf %d %lf %lf %d %lf %lf\n", 
	   &symb, &act, &dlixo, &dlixo, &ang, &dlixo, &dlixo, &die, &dlixo, &dlixo); 
    cnct[i].cdist     = act;
    cnct[i].cangle    = ang;
    cnct[i].cdihedral = die;
  }
  rewind(fl0);
  fclose(fl0);
}
  
