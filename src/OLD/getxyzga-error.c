#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "defmol.h"
#include "getxyz.h"

int GetXYZGA(const char *fdat, int nat, double *energy, defmol *indv)
{
  char    lixo[20];
  char    ident[20];
  char    first[29] = ("GEOMETRY");
  int     j;
  double  energ;
  
   FILE   *arq;
   
   printf("OLA\n");
   
  if ( ! (arq=fopen(fdat,"r")) ) {
    printf("|\n\nError [GetXYZ]: The file \"%s\" is not available for use!\n",fdat);
    exit(EXIT_FAILURE);
  }

  do { 
    fscanf( arq, "%s", &ident );
  } while( !feof(arq) && strncasecmp(ident,"GEOMETRY", 8) );
  //    printf("\n%s\n", ident);
    //   printf("\n%s\n", first);
  if( strcmp(ident,first) == 0 ){
    fseek ( arq , 231 , SEEK_CUR );
    //   fscanf(arq,"%s", ident);
    //    printf("\n%s\n", ident);
    
    for(j = 0; j < nat; j++){
      fscanf(arq,"%s %s %lf %lf %lf",  &indv[j].symb, &lixo, &indv[j].xyz->x, &indv[j].xyz->y, &indv[j].xyz->z);
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
