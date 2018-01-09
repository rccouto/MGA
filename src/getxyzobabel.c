#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "defmol.h"
#include "getxyz.h"

int GetXYZOBabel(const char *fxyz, int nat, defmol *indv)
{
  char    lixo[20], lixo1[20];
  char    ident[20];
  char    first[37] = ("** Molecular cartesian coordinates **");
  int     j;
  double  energ;
  
  FILE   *arq;

  if ( ! (arq=fopen(fxyz,"r")) ) {
    printf("|\n\nError [GetXYZOBabel]: The file \"%s\" is not available for use!\n",fxyz);
    exit(EXIT_FAILURE);
  }
  
  do { 
    fgets(ident, 38, arq);
  } while( !feof(arq) && strncasecmp(ident,"** Molecular cartesian coordinates **", 37) );
  
  if( strncasecmp(ident,first,37) == 0 ){
    for(j = 0; j < nat; j++){
      fscanf(arq,"%s %lf %lf %lf",  &lixo1, &indv[j].xyz->x, &indv[j].xyz->y, &indv[j].xyz->z);
    } 
  }
  else{
    printf("@");
  }
  
  fclose(arq);
  return 0;
}
