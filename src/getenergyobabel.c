#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "defmol.h"
#include "getenergy.h"

#define ZERO = 0.0E+0;

int GetEnergyOBabel(const char *fout, int nat, double *energy, defmol *indv)
{
  char    lixo[20];
  char    ident[60];
  char    first[12] = ("TOTAL ENERGY");
  double  energ;
  
  FILE   *arq;
  
  if ( ! (arq=fopen(fout,"r")) ) {
    printf("|\n\nError [GetEnergyOBabel]: The file \"%s\" is not available for use!\n",fout);
    exit(EXIT_FAILURE);
  }
  
  do { 
    fgets(ident, 13, arq);
  } while( !feof(arq) && strncasecmp(ident,"TOTAL ENERGY", 12) );
  
  //printf("\n*** Debug [GetXYZ]: ident = %s\n", ident);
  
  if( strncasecmp(ident,first,12) == 0 ){ 
    fscanf( arq, "%s %lf %s", &lixo, &energ, &lixo);
    *energy = energ;
    //printf("\n\t Energy: %10.10lf \n", *energy);/**/
  }
  else{
    // Energy is ZERO.
     energ = 9999999999;
     *energy = energ;
  }
  fclose(arq);
  return 0;
}
