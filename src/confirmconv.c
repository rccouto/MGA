#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "geninp.h"
#include "genxyz.h"
#include "genzmat.h"
#include "runprog.h"
#include "getxyzga.h"
#include "getxyzobabel.h"
#include "getenergyobabel.h"
#include "reoptpop.h"

extern char *eff;

int ConfirmConv(char *prog, char *gamesspath, int pos, int natoms, int nindiv, int nstep, int *conv, int gmsver, int ncore, double *energy, defmol *MPop, Txyz *XYZ, Tcnct *cnct, Tzmat *zmat)
{
  char   xyzname[15], datname[15], outname[15], rstname[15], vecname[15], newname[15];
  double *penergy, tempenergy;
  int    *pconv, nd, ip;
  defmol *pMPop;
  Txyz   *pXYZ;
  Tzmat  *pzmat;

  FILE *fltest, *fltest2;
  
  nd = natoms*pos;
  pconv =   &conv[pos];
  penergy = &energy[pos];
  pMPop =   &MPop[nd];
  pXYZ = &XYZ[nd];
  pzmat = &zmat[0];

  //printf("[CONFIRM]: #1 Energy=%lf\n", energy[pos]);

  if ( strncasecmp(prog,"OBabel", 6) ){
    /* Generate  Input for Eletronic Structure Calculations */
    tempenergy=sqrt(pow(*penergy,2));
    gcvt(tempenergy, 8, newname);
    
    sprintf(xyzname,"%s.xyz", newname);
    
    printf("#");fflush(stdout);
    GenZMat (natoms, prog, xyzname, pMPop, pXYZ, cnct, pzmat);
    printf("#");fflush(stdout);
    
    printf("#");fflush(stdout);
    RunProg (1, 1, gmsver, ncore, "I", prog, gamesspath, newname, eff);
    printf("#");fflush(stdout);
    
    sprintf(datname,"%s.dat", newname);
    sprintf(outname,"%s.out", newname);
    sprintf(rstname,"%s.rst", newname);
    sprintf(vecname,"%s.vec", newname);  
    
    printf("#");fflush(stdout);
    GetXYZGA (datname, rstname, natoms, nstep, penergy, pconv, pMPop);
    printf("#");fflush(stdout);
    
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
  }
  
  else{
    
    if      ( pos < 9   ) sprintf(xyzname, "OBabelM000%d.xyz", pos + 1);
    else if ( pos < 99  ) sprintf(xyzname, "OBabelM00%d.xyz", pos + 1);
    else if ( pos < 999 ) sprintf(xyzname, "OBabelM0%d.xyz", pos + 1);
    else if ( pos < 9999) sprintf(xyzname, "OBabelM%d.xyz", pos + 1);
    else               exit(EXIT_FAILURE);
    
    ip = pos*natoms;
    pMPop = &MPop[ip];
    
    GenXYZ (natoms, prog, xyzname, pMPop);
    
    RunProg (pos + 1, 0, gmsver, ncore, "I", prog, gamesspath, "Gamess", eff);
    
    if      ( pos < 9   ) sprintf(xyzname, "OBabelE000%d.xyz", pos + 1);
    else if ( pos < 99  ) sprintf(xyzname, "OBabelE00%d.xyz", pos + 1);
    else if ( pos < 999 ) sprintf(xyzname, "OBabelE0%d.xyz", pos + 1);
    else if ( pos < 9999) sprintf(xyzname, "OBabelE%d.xyz", pos + 1);
    else               exit(EXIT_FAILURE);
    
    penergy = &energy[pos];
    
    GetXYZOBabel (xyzname, natoms, pMPop);  
    
    if      ( pos < 9   ) sprintf(outname, "OBabelE000%d.out", pos + 1);
    else if ( pos < 99  ) sprintf(outname, "OBabelE00%d.out", pos + 1);
    else if ( pos < 999 ) sprintf(outname, "OBabelE0%d.out", pos + 1);
    else if ( pos < 9999) sprintf(outname, "OBabelE%d.out", pos + 1);
    else               exit(EXIT_FAILURE);
    
    GetEnergyOBabel (outname, natoms, penergy, pMPop);  
    
    /* Remove .xyz Files */
    if ( ! (fltest2=fopen(xyzname,"r")) ) {
      printf("|\n\nError [Main]: The file \"%s\" is not available for deleting!\n", xyzname);
      exit(EXIT_FAILURE);
    }
    fclose(fltest2);
    unlink(xyzname);
    
    /* Remove .out Files */
    if ( ! (fltest2=fopen(outname,"r")) ) {
      printf("|\n\nError [Main]: The file \"%s\" is not available for deleting!\n",outname);
      exit(EXIT_FAILURE);
    }
    fclose(fltest2);
    unlink(outname);
  }
}

