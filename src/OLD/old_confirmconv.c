#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defmol.h"
#include "geninp.h"
#include "genxyz.h"
#include "runprog.h"
#include "getxyzga.h"
#include "reoptpop.h"


int ConfirmConv(char *prog, char *gamesspath, int pos, int natoms, int nindiv, int nstep, int *conv, double *energy, defmol *MPop)
{
  char   xyzname[15], datname[15], outname[15], rstname[15], vecname[15], newname[15];
  double *penergy, tempenergy;
  int    *pconv, nd;
  defmol *pMPop;
  
  FILE *fltest, *fltest2;
  
  nd = natoms*pos;

  pconv =   &conv[pos];
  penergy = &energy[pos];
  pMPop =   &MPop[nd];

  //printf("[CONFIRM]: #1 Energy=%lf\n", energy[pos]);

  /* Generate  Input for Eletronic Structure Calculations */
  tempenergy=sqrt(pow(*penergy,2));
  gcvt(tempenergy, 8, newname);

  sprintf(xyzname,"%s.xyz", newname);
  
  printf("#");fflush(stdout);
  GenXYZ (natoms, prog, xyzname, pMPop);
  printf("#");fflush(stdout);

  printf("#");fflush(stdout);
  RunProg (1, 1, "I", prog, gamesspath, newname);
  printf("#");fflush(stdout);

  sprintf(datname,"%s.dat", newname);
  sprintf(outname,"%s.out", newname);
  sprintf(rstname,"%s.rst", newname);
  sprintf(vecname,"%s.vec", newname);  

  printf("#");fflush(stdout);
  printf("[CONFIRM]");
  GetXYZGA (datname, rstname, natoms, nstep, penergy, pconv, pMPop);
  printf("#");fflush(stdout);
  
  /* Saving Wave-Function */ 
  if (conv[pos] != 1){
    tempenergy=sqrt(pow(*penergy,2));
    gcvt(tempenergy, 8, newname);
    sprintf(newname, "%s.vec", newname);
    rename(vecname, newname);
  }/**/
  
  /* Remove .dat Files */ 
  if ( ! (fltest=fopen(datname,"r")) ) {
    printf("|\n\nError [Main]: The file \"%s\" is not available for deleting!\n",datname);
    exit(EXIT_FAILURE);
  }
  fclose(fltest);
  unlink(datname);
  
  /* Remove .out Files */
  if ( ! (fltest2=fopen(outname,"r")) ) {
    printf("|\n\nError [Main]: The file \"%s\" is not available for deleting!\n",outname);
    exit(EXIT_FAILURE);
  }
  fclose(fltest2);
  unlink(outname);
  
  /* Remove .rst Files */ 
  unlink(rstname);
  
  /* Remove .vec Files */ 
  unlink(vecname);

  //printf("[CONFIRM]: #2 Energy=%lf\n", energy[pos]);

}
