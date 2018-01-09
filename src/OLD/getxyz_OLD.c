#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "defmol.h"
#include "getxyz.h"

int GetXYZ(const char *fdat, int nat, defmol *indv, int i, double *energy)
{
  char    lixo[20], lixo1[20];
  char    ident[20];
  char    linha[61];
  char    first[29] = ("GEOMETRY");
  char    second[11] = ("ELECTRONIC");
  char    symb[3];
  int     j, n, num, sum, en;
  double  x, y, z;

  FILE   *arq;



  arq = fopen(fdat,"r");
  printf("** en=%d sum=%d\n", en, sum);

  //  printf("Debug [GetXYZ]: #1\n");
  while(fscanf(arq,"%s", &ident) != EOF){ 
    fscanf(arq,"%s", &ident);
    if( strcmp(ident,first) == 0 ){
      //      printf("Debug [GetXYZ]: #2\n");
      for (j = 0;  j == 3; j++) fscanf(arq, "%s\n" , &lixo);
      printf("*2* en=%d sum=%d\n", en, sum);
      for(j = 0; j < nat - 1; j++){
	printf("OPA\n");
	fscanf(arq,"%s %s %lf %lf %lf",  &lixo1, &lixo, &indv[j].xyz->x, &indv[j].xyz->y, &indv[j].xyz->z);
	printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", lixo1, indv[j].xyz->x, indv[j].xyz->y,  indv[j].xyz->z);
      }    
    } 
    if(strcmp(ident,second)==0){
      //      printf("Debug [GetEnergy]: #2\n");
      fscanf(arq, "%s %s %lf" , &lixo, &lixo, &energy[i]);
      //      printf("\n\tElectronic Energy: %15.8lf \n", energy[i]);
    } 		
  }
  fclose(arq);
  return 0;
}
