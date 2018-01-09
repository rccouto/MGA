#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "defmol.h"
#include "getenergy.h"

int GetEnergy(const char *fout, int nat, double *energy, defmol *indv)
{
  char    lixo[20];
  char    ident[60];
  char    first[20] = ("FINAL");
  double  energ;
  
  FILE   *arq;
  
  if ( ! (arq=fopen(fout,"r")) ) {
    printf("|\n\nError [GetXYZ]: The file \"%s\" is not available for use!\n",fout);
    exit(EXIT_FAILURE);
  }
  
  do 
    { 
      fscanf( arq, "%s", &ident );
    } 
  while( !feof(arq) && strncasecmp(ident,"FINAL", 6) );

  //  printf("\n*** Debug [GetXYZ]: ident = %s#1\n", ident);
  
  if( strcmp(ident,first) == 0 ){
    
    fscanf( arq, "%s %s %s %lf", &lixo, &lixo, &lixo, &energ);
    *energy = energ;
    //    printf("\n\tElectronic Energy: %10.10lf \n", *energy);/**/
  }
  else{
    printf("\n\tElse!");
  }
  
  fclose(arq);
  return 0;
}
