#include <stdio.h>
#include <stdlib.h>

#include "defmol.h"
#include "genxyz.h" 

void GenXYZ (const int nat, const char *prog, const char *fname, defmol *indv)
{
  int    j;
  
  FILE   *fl;
  
  fl = fopen(fname,"w");
  
  if ( !strncasecmp(prog,"DALTON", 6) ){
    fprintf(fl,"\noi");
    for (j = 0; j < nat; j++){
      printf("\n%s  %3.0d  %17.12lf  %17.12lf  %17.12lf",
	     indv[j].symb, indv[j].np, &indv[j].xyz->x, &indv[j].xyz->y, &indv[j].xyz->z);
      fprintf(fl,"\n%s  %3.0d  %17.12lf  %17.12lf  %17.12lf",
	      indv[j].symb, indv[j].np, &indv[j].xyz->x, &indv[j].xyz->y, &indv[j].xyz->z);
    }
  }
  else if ( !strncasecmp(prog,"GAMESS", 6) ){
    fprintf(fl," $DATA");
    fprintf(fl,"\n ** Molecular cartesian coordinates **");
    fprintf(fl,"\n C1");
    
    /*printf("\nInside genXYZ:\n");
      for (j = 0; j < nat; j++){
      printf("\n %s  %3d   %17.12lf  %17.12lf  %17.12lf",
      indv[j].symb, indv[j].np, indv[j].xyz->x, indv[j].xyz->y, indv[j].xyz->z);
      }/**/
    
    for (j = 0; j < nat; j++){
      fprintf(fl,"\n %s  %3d   %17.12lf  %17.12lf  %17.12lf",
	      indv[j].symb, indv[j].np, indv[j].xyz->x, indv[j].xyz->y, indv[j].xyz->z); 
    }
    fprintf(fl,"\n $END\n");
  }
  else if ( !strncasecmp(prog,"OBabel", 6) ){
    fprintf(fl," %d", nat);
    fprintf(fl,"\n ** Molecular cartesian coordinates **");
    for (j = 0; j < nat; j++){
      fprintf(fl,"\n %s  %17.12lf  %17.12lf  %17.12lf",
	      indv[j].symb, indv[j].xyz->x, indv[j].xyz->y, indv[j].xyz->z); 
    }
  }
  
  else{
    printf("\n\nError: eeecp \"%s\" is not available!", prog);
    exit(1);
  }
  fclose(fl);
}
