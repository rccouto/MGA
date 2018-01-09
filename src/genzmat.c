#include <stdio.h>
#include <stdlib.h>

#include "defmol.h"
#include "genzmat.h" 
#include "xyz2zmat.h"

void GenZMat (const int nat, const char *prog, const char *fname, defmol *indv, Txyz *XYZ, Tcnct *cnct, Tzmat *zmat)
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
    xyz2zmat(nat, XYZ, cnct, zmat);
    
    fprintf(fl," $DATA");
    fprintf(fl,"\n ** Molecular Z-Matrix coordinates **");
    fprintf(fl,"\n C1");
    
    fprintf(fl,"\n %s", indv[0].symb);
    fprintf(fl,"\n %s %3d %2.2lf", indv[1].symb, cnct[1].cdist, zmat[1].dist);
    fprintf(fl,"\n %s %3d %2.2lf %3d %2.2lf", indv[2].symb, cnct[2].cdist, zmat[2].dist, cnct[2].cangle, zmat[2].angle);
    for (j = 3; j < nat; j++){
      fprintf(fl,"\n %s %3d %2.2lf %3d %2.2lf %3d %2.2lf", indv[j].symb, cnct[j].cdist, zmat[j].dist, cnct[j].cangle, zmat[j].angle, cnct[j].cdihedral, zmat[j].dihedral);
    }
    fprintf(fl,"\n $END\n");
  }

  else{
    printf("\n\nError: eeecp \"%s\" is not available!", prog);
    exit(1);
  }
  fclose(fl);
}
