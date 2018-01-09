#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "genindv.h"
#include "geninp.h"
#include "genxyz.h"
#include "runprog.h"
#include "getxyzga.h"
#include "getenergy.h"
#include "order.h"
#include "mating.h"
#include "rdzmat.h"
#include "mutation.h"
#include "rtgroup.h"
#include "checkconv.h"
#include "reoptpop.h"


int NGA(const char *Infile, char *prog, char *gamesspath, char *optstruct, char *rmoutfile, int natoms, int nindiv, int nstep, int percmut, int typmat, int nrotgp, int mutmin, int mutmax, int genafconv, int ngen, int *conv, double *energy, double absenrgy, long *nseed, defmol *MPop, Txyz *XYZ, Tcnct *cnct, Tzmat *zmat, Trtt *rtt, Tatom  *rttatm)
{
  char     xyzname[15], datname[15], outname[15], rstname[15], fname[15];
  int      i, ip;
  int      gen = 0, *pconv, tindiv;
  int      nonhydro, mutpop, matpop;
  int      genconv=0, nconv, halfindiv, halfgen, gentest;
  double  *penergy;
  defmol  *pMPop;
  Txyz    *pXYZ;
  
  FILE   *fltest, *fltest2;

  nindiv = nindiv/2;

  /* Print Generation's Energy */
  printf("\n\n**  Generation's [0] Energy:\n");
  for ( i = 0 ; i < nindiv; i++){
    if ( conv[i] == 1 ) printf("Energy[%d] = %lf *\n", i+1, energy[i]);
    else printf("Energy[%d] = %lf\n", i+1, energy[i]);
  }

  RdZMat(Infile, natoms, cnct);
  
  /* Counting Non-Hydrogen Atoms */
  nonhydro = 0;
  for ( i = 0; i < natoms; i++){
    if ( *MPop[i].symb != 'H' ){
      nonhydro = nonhydro + 1;
    }
  }

  /* Seting the number of molecules for Mating and Mutation */
  mutpop = ceil( ( nindiv*percmut ) /100 ); 
  matpop = nindiv - mutpop;
  
  /* GA Loop! */
  do{  
    printf("\n\n\n** GENERATION [%d] **", gen+1);
    
    /* Mating */
    printf("\n\n**  Generating New Molecules with Mating Operator: |");
    for (i = 0; i < matpop; i++){
      printf("#");fflush(stdout);
      Mating(typmat, natoms, nindiv, i, nonhydro, nseed, XYZ, MPop, cnct, zmat);
      printf("#");fflush(stdout);
    }
    printf("|");fflush(stdout);
    
    /* Mutation */
    printf("\n\n**  Generating New Molecules with Mutation Operator: |");
    for (i = matpop; i < nindiv; i++){
      printf("#");fflush(stdout);
      Mutation(natoms, nindiv, nrotgp, i, mutmin, mutmax, nseed, rtt, rttatm, XYZ);
      printf("#");fflush(stdout);
    }
    printf("|");fflush(stdout);
    
    /* Generate  Input for Eletronic Structure Calculations */
    printf("\n\n**  Generating Input for Eletronic Structure Calculations in GA: |");
    
    for ( i = nindiv; i < 2*nindiv; i++ ) {
      
      if      ( i < 9   ) sprintf(xyzname, "GamessP000%d.xyz", i + 1);
      else if ( i < 99  ) sprintf(xyzname, "GamessP00%d.xyz", i + 1);
      else if ( i < 999 ) sprintf(xyzname, "GamessP0%d.xyz", i + 1);
      else if ( i < 9999) sprintf(xyzname, "GamessP%d.xyz", i + 1);
      else               exit(EXIT_FAILURE);
      
      printf("#");fflush(stdout);
      ip = i*natoms;
      pMPop = &MPop[ip];

      printf("#");fflush(stdout);
      GenXYZ (natoms, prog, xyzname, pMPop);
      printf("#");fflush(stdout);
    }
    printf("|");fflush(stdout);
    
    printf("\n\n**  Running Eletronic Structure Calculations in GA: |");
    
    for ( i = nindiv; i < 2*nindiv; i++ ) {
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
      
      ip = i*natoms;
      pMPop = &MPop[ip];  
      penergy = &energy[i];
      pconv = &conv[i];
      
      if ( !strncasecmp(optstruct,"on", 2) ){
	printf("#");fflush(stdout);
	GetXYZGA (datname, rstname, natoms, nstep, penergy, pconv, pMPop);  
	printf("#");fflush(stdout);
      }
      else {
	printf("#");fflush(stdout);
	GetEnergy (outname, natoms, penergy, pMPop);  
	printf("#");fflush(stdout);
      }
      
      /* Remove .dat Files */ 
      if ( !strncasecmp(rmoutfile,"on", 2 ) ){
	if ( ! (fltest=fopen(datname,"r")) ) {
	  printf("|\n\nError [Main]: The file \"%s\" is not available for deleting!\n",datname);
	  exit(EXIT_FAILURE);
	}
	fclose(fltest);
	unlink(datname);/**/
      }
      
      /* Remove .rst Files */ 
      unlink(rstname);
      
      /* Remove .out Files */
      if ( !strncasecmp(rmoutfile,"on", 2 ) ){
	if ( ! (fltest2=fopen(outname,"r")) ) {
	  printf("|\n\nError [Main]: The file \"%s\" is not available for deleting!\n",outname);
	  exit(EXIT_FAILURE);
	}
	fclose(fltest2);
	unlink(outname);
      }
    }
    printf("|                                                ");fflush(stdout);
    
    tindiv = 2*nindiv;
    
    /* Checking Structures Convergence */
    nconv = 0;
    for ( i = 0; i < nindiv; i++ ){
      nconv = nconv + conv[i];
    }
    halfindiv = ceil(nindiv/2); 
    if (nconv == halfindiv ){
      ReoptPop(prog, gamesspath, xyzname, datname, outname, rstname, natoms, tindiv, nstep, 1, energy, MPop);
      genconv = 0;
    }
    
    /* Confirming Convergence */
    halfgen = ceil(genafconv/2);
    if( genconv == halfgen ){ 
      ReoptPop(prog, gamesspath, xyzname, datname, outname, rstname, natoms, tindiv, nstep, 1, energy, MPop);
    }
    gentest = genafconv - 5;
    if( genconv == gentest ) {
      ReoptPop(prog, gamesspath, xyzname, datname, outname, rstname, natoms, tindiv, nstep, 1, energy, MPop);
    }
    
    gentest = genafconv - 2;
    if( genconv == gentest ) {
      ReoptPop(prog, gamesspath, xyzname, datname, outname, rstname, natoms, tindiv, nstep, 1, energy, MPop);
    }
    
    gentest = ngen - 10;
    if( gen == gentest ) ReoptPop(prog, gamesspath, xyzname, datname, outname, rstname, natoms, tindiv, nstep, 1, energy, MPop);
    
    /* Counting Generations before convergence */
    if ( fabs(energy[0] - energy[nindiv-1]) < absenrgy ) genconv = genconv + 1;
    
    tindiv = 2*nindiv;
    Order(natoms, tindiv, energy, conv, MPop);
    
    /* Print the Generation's Energy */
    printf("\n\n**  Generation's [%d] Energy:\n", gen+1);
    for ( i = 0 ; i < nindiv; i++){
      if ( conv[i] == 1) printf("Energy[%d] = %lf *\n", i+1, energy[i]);
      else printf("Energy[%d] = %lf\n", i+1, energy[i]);
    }
    
    /* Counting the number of Generations */
    gen = gen + 1;
    if( gen == ngen ) break;
    
  }while( genconv != genafconv );
  unlink(fname);
}

