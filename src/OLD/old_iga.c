#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "alloc.h"
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
#include "predator.h"
#include "confirmconv.h"
#include "history.h"
#include "mpoptofpop.h"

#define ZERO          0.0E+0

int IGA(const char *Infile, char *prog, char *gamesspath, char *optstruct, char *rmoutfile, int natoms, int nindiv, int nstep, int percmut, int typmat, int nrotgp, int mutmin, int mutmax, int genafconv, int ngen, int *conv, double *energy, double *Fenergy, double absenrgy, long *nseed, defmol *MPop, defmol *FPop, Txyz *XYZ, Tcnct *cnct, Tzmat *zmat, Trtt *rtt, Tatom  *rttatm)
{
  char     vecname[15], newname[15];
  char     xyzname[15], datname[15], outname[15], rstname[15], fname[15];
  int      i, j, ip, ni, pos;
  int      gen = 0, *pconv, tindiv;
  int      nonhydro, mutpop, matpop, nenerg;
  int      genconv = 0, halfindiv, halfgen, gentest;
  double  *penergy, Eavrg, EavrgPrev = ZERO, tempenergy, cac;
  double   Emin, Emax, enerc;
  defmol  *pMPop;
  Txyz    *pXYZ, *XYZ_temp;
  
  FILE   *fltest, *fltest2;
  
  
  /* Applying Predator Operator */
  printf("\n\n** Applying Predator Operator: |");
  printf("#");fflush(stdout);
  Predator(natoms, nindiv, conv, absenrgy, energy, MPop);
  printf("#");fflush(stdout);
  printf("|");fflush(stdout);
  
  nindiv = nindiv/2;

  /* Print Generation's Energy */
  printf("\n\n**  Generation's [0] Energy:\n");
  for ( i = 0 ; i < nindiv; i++){
    if ( conv[i] == 1 ) printf("Energy[%d] = %lf *\n", i+1, energy[i]);
    else printf("Energy[%d] = %lf\n", i+1, energy[i]);
  }
  
  /* Passing First Generation's Structures to the Final Population */
  for ( i = 0; i < nindiv; i++ ){
    energy[2*nindiv + i] = energy[i];
    pos=(2*nindiv + i)*natoms;
    ip=i*natoms;
    for ( j = 0; j < natoms; j++ ){
      XYZ[pos + j].x = XYZ[ip + j].x;
      XYZ[pos + j].y = XYZ[ip + j].y;
      XYZ[pos + j].z = XYZ[ip + j].z;
    }
  }
  
  /* Print Final Population's Energy */
  printf("\n\n**  Final Population's Energy - GEN[0]\n");
  for ( i = 2*nindiv ; i < 3*nindiv; i++){
    printf("Energy[%d] = %lf\n", i+1, energy[i]);
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
  //printf("Structures Generated by: Mutation = %d, Mating=%d\n", mutpop, matpop); 

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
      
      if      ( i < 9   ) sprintf(vecname, "GamessP000%d.vec", i + 1);
      else if ( i < 99  ) sprintf(vecname, "GamessP00%d.vec", i + 1);
      else if ( i < 999 ) sprintf(vecname, "GamessP0%d.vec", i + 1);
      else if ( i < 9999) sprintf(vecname, "GamessP%d.vec", i + 1);
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
      
      /* Saving Wave-Function */ 
      tempenergy=sqrt(energy[i]*energy[i]);
      gcvt(tempenergy, 8, newname);
      sprintf(newname, "%s.vec", newname);
      rename(vecname, newname);
      
      /* Using History Operator */
      printf("#");fflush(stdout);
      History(prog, gamesspath, i, nindiv, natoms, nstep, conv, energy, absenrgy, MPop, XYZ);
      printf("#");fflush(stdout);
      
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
      /* Remove .vec Files */ 
      unlink(vecname);
    }
    printf("|                                                ");fflush(stdout);
 
    tindiv=2*nindiv;

    printf("\n\n** Applying Predator Operator: |");
    printf("#");fflush(stdout);
    Predator(natoms, tindiv, conv, absenrgy, energy, MPop);
    printf("#");fflush(stdout);
    printf("|");fflush(stdout);

    /* Counting Generations before convergence */
    Eavrg = ZERO;
    for ( i = 2*nindiv; i < 3*nindiv; i++ ){
      Eavrg = Eavrg + energy[i];
    }
    Eavrg = Eavrg/nindiv;
 
    if( Eavrg == EavrgPrev ){ 
      genconv = genconv + 1;
    }
    else genconv = ZERO;
    EavrgPrev = Eavrg;
    
    /* Print the Generation's Energy */
    printf("\n\n**  Generation's [%d] Energy:\n", gen+1);
    for ( i = 0 ; i < nindiv; i++){
      if ( conv[i] == 1) printf("Energy[%d] = %lf *\n", i+1, energy[i]);
      else printf("Energy[%d] = %lf\n", i+1, energy[i]);
    }

    /* Print Final Population's Energy */
    printf("\n\n***  Final Population's Energy - GEN[%d]:\n =======================================\n", gen+1);
    for ( i = 2*nindiv ; i < 3*nindiv; i++){
      printf("    Energy[%d] = %lf\n", i+1, energy[i]);
    } /**/

    /* Final Population's Energy Analysis */
    Eavrg = 0.0E+0;
    for ( i = 2*nindiv; i < 3*nindiv; i++ ){
      Eavrg = Eavrg + energy[i];
      if ( energy[i] < Emin || i == 2*nindiv ) Emin = energy[i];
      if ( energy[i] != ZERO && energy[i] > Emax || i == 2*nindiv ) Emax = energy[i];
    }
    Eavrg = Eavrg/nindiv;
    printf("\n ** Analysis: \n");
    printf("   Minimum Electronic Energy = %10.10lf\n   Average Electronic Energy = %10.10lf\n   Maximum Electronic Energy = %10.10lf\n", Emin, Eavrg, Emax); 
    printf("\n =============================================\n");

    /* Counting the number of Generations */
    gen = gen + 1;
    if( gen == ngen ) break;
    
  }while( genconv != genafconv );
  unlink(fname);


    /*---------------------------------FINAL ANALYSIS--------------------------------------*/
    Eavrg = 0.0E+0;
    for ( i = 2*nindiv; i < 3*nindiv; i++ ){
      Eavrg = Eavrg + energy[i];
      if ( energy[i] < Emin || i == 2*nindiv ) Emin = energy[i];
      if ( energy[i] != ZERO && energy[i] > Emax || i == 2*nindiv ) Emax = energy[i];
    }
    Eavrg = Eavrg/nindiv;
    printf("\n\n ** Energy analysis ** \n =====================\n");
    printf("   Minimum Electronic Energy = %10.10lf\n   Average Electronic Energy = %10.10lf\n   Maximum Electronic Energy = %10.10lf\n", Emin, Eavrg, Emax); 
    
    printf("\n       *Energy*\t\t*No. of Structures*\n");
    for ( i = 2*nindiv; i < 3*nindiv; i++){
      printf("    %12.6lf\t         1\n", energy[i]);
    }
    
    printf("\n ** Molecular Structures in Cartesian Coordinates **\n ===================================================\n");
    for ( i = 2*nindiv; i < 3*nindiv; i++ ) {
      ip = i*natoms;
      printf("%d\n", natoms);/**/
      printf("Electronic Energy (Hartree) = %10.10lf\n", energy[i]);/**/
      for(j = 0; j < natoms; j++ ){
	printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[ip+j].x, XYZ[ip+j].y, XYZ[ip+j].z);
      }
    }
    printf("\nDone!\n");

}
