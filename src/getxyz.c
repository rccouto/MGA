#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "defmol.h"
#include "getxyz.h"

int GetXYZ(const char *fdat, int nat, double *energy, defmol *indv)
{
  char    lixo[20], lixo1[20];
  char    ident[20];
  char    first[29] = ("----- RESULTS FROM SUCCESSFUL");
  int     j;
  double  energ;
  
  FILE   *arq;

  if ( ! (arq=fopen(fdat,"r")) ) {
    printf("|\n\nError [GetXYZ]: The file \"%s\" is not available for use!\n",fdat);
    exit(EXIT_FAILURE);
  }
 
  do { 
    fgets(ident, 30, arq);
  } while( !feof(arq) && strncasecmp(ident,"----- RESULTS FROM SUCCESSFUL", 29) );
  
  if( strncasecmp(ident,first,29) == 0 ){
    fseek ( arq , 250 , SEEK_CUR );
    for(j = 0; j < nat; j++){
      fscanf(arq,"%s %s %lf %lf %lf",  &lixo1, &lixo, &indv[j].xyz->x, &indv[j].xyz->y, &indv[j].xyz->z);
    } 
    fseek ( arq , 74, SEEK_CUR );
    fscanf( arq, "%lf", &energ );
    *energy = energ;
  }
  else{
    printf("@");
  }
  
  fclose(arq);
  return 0;
}
