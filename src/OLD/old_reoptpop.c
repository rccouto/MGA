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

extern int  *conv;

int ReoptPop(char *prog, char *gamesspath, char *xyzname, char *datname, char *outname, char *rstname, int natoms, int nindiv, int nstep, int typga, double *energy, defmol *MPop)
{
  char   vecname[15], newname[15], nbsname[15];
  int    i, ip;
  int    *pconv;
  double *penergy, tempenergy;
  defmol *pMPop;
  
  FILE *fltest, *fltest2;
  
  /* Generate  Input for Eletronic Structure Calculations */
  printf("\n\n**  Re-Optimizing Population: |");
  
  for ( i = 0; i < nindiv; i++ ) {
    if      ( i < 9   ) sprintf(xyzname, "GamessP000%d.xyz", i + 1);
    else if ( i < 99  ) sprintf(xyzname, "GamessP00%d.xyz", i + 1);
    else if ( i < 999 ) sprintf(xyzname, "GamessP0%d.xyz", i + 1);
    else if ( i < 9999) sprintf(xyzname, "GamessP%d.xyz", i + 1);
    else                exit(EXIT_FAILURE);
    
    printf("#");fflush(stdout);
    ip = i*natoms;
    pMPop = &MPop[ip];
    
    printf("#");fflush(stdout);
    GenXYZ (natoms, prog, xyzname, pMPop);
    printf("#");fflush(stdout);
  }
  
  for ( i = 0; i < nindiv; i++ ) {
    printf("#");fflush(stdout);
    RunProg (i + 1, 0, "I", prog, gamesspath, "Gamess");
    printf("#");fflush(stdout);
    
    if      ( i < 9   ) sprintf(datname, "GamessP000%d.dat", i + 1);
    else if ( i < 99  ) sprintf(datname, "GamessP00%d.dat", i + 1);
    else if ( i < 999 ) sprintf(datname, "GamessP0%d.dat", i + 1);
    else if ( i < 9999) sprintf(datname, "GamessP%d.dat", i + 1);
    else               exit(EXIT_FAILURE);
    
    if      ( i < 9   ) sprintf(outname, "GamessP000%d.out", i + 1);
    else if ( i < 99  ) sprintf(outname, "GamessP00%d.out", i + 1);
    else if ( i < 999 ) sprintf(outname, "GamessP0%d.out", i + 1);
    else if ( i < 9999) sprintf(outname, "GamessP%d.out", i + 1);
    else               exit(EXIT_FAILURE);
    
    if      ( i < 9   ) sprintf(rstname, "GamessP000%d.rst", i + 1);
    else if ( i < 99  ) sprintf(rstname, "GamessP00%d.rst", i + 1);
    else if ( i < 999 ) sprintf(rstname, "GamessP0%d.rst", i + 1);
    else if ( i < 9999) sprintf(rstname, "GamessP%d.rst", i + 1);
    else               exit(EXIT_FAILURE);
    
    if      ( i < 9   ) sprintf(vecname, "GamessP000%d.vec", i + 1);
    else if ( i < 99  ) sprintf(vecname, "GamessP00%d.vec", i + 1);
    else if ( i < 999 ) sprintf(vecname, "GamessP0%d.vec", i + 1);
    else if ( i < 9999) sprintf(vecname, "GamessP%d.vec", i + 1);
    else               exit(EXIT_FAILURE);

    if      ( i < 9   ) sprintf(nbsname, "GamessP000%d.nbs", i + 1);
    else if ( i < 99  ) sprintf(nbsname, "GamessP00%d.nbs", i + 1);
    else if ( i < 999 ) sprintf(nbsname, "GamessP0%d.nbs", i + 1);
    else if ( i < 9999) sprintf(nbsname, "GamessP%d.nbs", i + 1);
    else               exit(EXIT_FAILURE);

    ip = i*natoms;
    pMPop = &MPop[ip];  
    penergy = &energy[i];
    pconv = &conv[i];
    
    printf("#");fflush(stdout);
    GetXYZGA (datname, rstname, natoms, nstep, penergy, pconv, pMPop);
    printf("#");fflush(stdout);
    

    /* Saving Wave-Function */ 
    if ( typga == 2 ){
	tempenergy=sqrt(energy[i]*energy[i]);
	gcvt(tempenergy, 8, newname);
	sprintf(newname, "%s.vec", newname);
	rename(vecname, newname);
    }

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
  printf("|                                                ");fflush(stdout);
}
