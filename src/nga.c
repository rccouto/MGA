#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "genindv.h"
#include "geninp.h"
#include "genxyz.h"
#include "genzmat.h"
#include "runprog.h"
#include "getxyzga.h"
#include "getenergy.h"
#include "order.h"
#include "mating.h"
#include "rdzmat.h"
#include "mutation.h"
#include "rtgroup.h"
#include "reoptpop.h"
#include "distmeasure.h"
#include "readbackup.h"

#define ZERO          0.0E+0

extern char *eff;
extern int  gen;

int NGA(const char *Infile, char *prog, char *gamesspath, char *optstruct, char *rmoutfile, char *backup, int natoms, int nindiv, int nstep, int percmut, int typmat, int nrotgp, int mutmin, int mutmax, int genafconv, int ngen, int *conv, int gmsver, int ncore, double *energy, double absenrgy, long *nseed, defmol *MPop, Txyz *XYZ, Tcnct *cnct, Tzmat *zmat, Trtt *rtt, Tatom  *rttatm)
{
  char     xyzname[15], datname[15], outname[15], rstname[15], fname[15], vecname[15];
  int      i, j, ip;
  int     *pconv, tindiv, typga;
  int      nonhydro, mutpop, matpop, popstg = ZERO;
  int      genconv = ZERO, nconv, halfindiv, halfgen, gentest;
  double  *penergy, Emin, Emax, Eavrg, EavrgPrev;
  defmol  *pMPop;
  Txyz    *pXYZ;
  Tzmat   *pzmat;

  FILE   *fltest, *fltest2;
  
  nindiv = nindiv/2;

  /* Print Generation's Energy */
  printf("\n\n**  Generation's [%d] Energy (Hartree):\n ==============================================\n\n", gen);
  for ( i = 0 ; i < nindiv; i++){
    if ( conv[i] == 1 ) printf("   E. Indiv.[%d] = %lf *\n", i+1, energy[i]);
    else printf("   E. Indiv.[%d] = %lf\n", i+1, energy[i]);
  }

  /* Population's Energy Analysis */ 
  Eavrg = 0.0E+0;
  for ( i = 0; i < nindiv; i++ ){
    Eavrg = Eavrg + energy[i];
    if ( energy[i] < Emin || i == 0 ) Emin = energy[i];
    if ( energy[i] != ZERO && energy[i] > Emax || i == 0 ) Emax = energy[i];
  }
  Eavrg = Eavrg/nindiv;
  
  printf("\n  --------------------------- \n\n");
  printf("   Minimum Electronic Energy = %10.6lf\n   Average Electronic Energy = %10.6lf\n   Maximum Electronic Energy = %10.6lf\n", Emin, Eavrg, Emax); 
  printf("\n ==============================================\n"); 

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
    
    /* If it's converged, do only Mutation */
    if ( genconv != 0 ){
      mutpop = nindiv;
      matpop = 0;	
    }

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
    
    
    // FOR GAMESS
    if ( strncasecmp(prog,"OBabel", 6) ){
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
	pXYZ = &XYZ[ip];
	pzmat = &zmat[0];
	
	printf("#");fflush(stdout);
	GenZMat (natoms, prog, xyzname, pMPop, pXYZ, cnct, pzmat);
	printf("#");fflush(stdout);
      }
      printf("|");fflush(stdout);
      
      printf("\n\n**  Running Eletronic Structure Calculations in GA: |");
      
      for ( i = nindiv; i < 2*nindiv; i++ ) {
	printf("#");fflush(stdout);
	RunProg (i + 1, 0, gmsver, ncore, "I", prog, gamesspath, "Gamess", eff);
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

	/* Remove .dat Files */ 
	if ( !strncasecmp(rmoutfile,"on", 2 ) ){
	  if ( ! (fltest=fopen(datname,"r")) ) {
	    printf("|\n\nError [Main]: The file \"%s\" is not available for deleting!\n",datname);
	    exit(EXIT_FAILURE);
	  }
	  fclose(fltest);
	  //unlink(datname);/**/
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
	  //unlink(outname);
	}
	/* Remove .vec Files */ 
	unlink(vecname);
      }
      printf("|                                                ");fflush(stdout);
    }
    
    // FOR OPEN BABEL
    else{
      /* Generate  Input for Eletronic Structure Calculations */
      printf("\n\n**  Generating Input for Eletronic Structure Calculations in GA: |");
     
      for ( i = nindiv; i < 2*nindiv; i++ ) {

	if      ( i < 9   ) sprintf(xyzname, "OBabelM000%d.xyz", i + 1);
	else if ( i < 99  ) sprintf(xyzname, "OBabelM00%d.xyz", i + 1);
	else if ( i < 999 ) sprintf(xyzname, "OBabelM0%d.xyz", i + 1);
	else if ( i < 9999) sprintf(xyzname, "OBabelM%d.xyz", i + 1);
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
	RunProg (i + 1, 0, gmsver, ncore, "I", prog, gamesspath, "Gamess", eff);
	printf("#");fflush(stdout);
	
	if      ( i < 9   ) sprintf(xyzname, "OBabelE000%d.xyz", i + 1);
	else if ( i < 99  ) sprintf(xyzname, "OBabelE00%d.xyz", i + 1);
	else if ( i < 999 ) sprintf(xyzname, "OBabelE0%d.xyz", i + 1);
	else if ( i < 9999) sprintf(xyzname, "OBabelE%d.xyz", i + 1);
	else               exit(EXIT_FAILURE);
	
	ip = i*natoms;
	pMPop = &MPop[ip];  
	penergy = &energy[i];
	
	printf("#");fflush(stdout);
	GetXYZOBabel (xyzname, natoms, pMPop);  
	printf("#");fflush(stdout);
	
	if      ( i < 9   ) sprintf(outname, "OBabelE000%d.out", i + 1);
	else if ( i < 99  ) sprintf(outname, "OBabelE00%d.out", i + 1);
	else if ( i < 999 ) sprintf(outname, "OBabelE0%d.out", i + 1);
	else if ( i < 9999) sprintf(outname, "OBabelE%d.out", i + 1);
	else               exit(EXIT_FAILURE);
	
	printf("#");fflush(stdout);
	//printf("[%d]", i);
	GetEnergyOBabel (outname, natoms, penergy, pMPop);  
	printf("#");fflush(stdout);
	
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
	
	/* Confirm if structure corresponds to molecule in study */
	ip = i*natoms;
	pXYZ  = &XYZ[ip];
	penergy = &energy[i];
	DistMeasure(prog, natoms, pXYZ, penergy, cnct);
      }
      printf("|                                                ");fflush(stdout);
    }
    
    tindiv = 2*nindiv;

    /* Checking Structures Convergence */
    nconv = 0;
    for ( i = 0; i < nindiv; i++ ){
      nconv = nconv + conv[i];
    }
    halfindiv = ceil(nindiv/2); 
    if (nconv == halfindiv ){
      ReoptPop(prog, gamesspath, xyzname, datname, outname, rstname, natoms, nindiv, nstep, 1, gmsver, ncore, energy, XYZ, cnct, zmat, MPop);
      genconv = 0;
    }

    /* Eliminate Spurious Structures */
    Eavrg = ZERO;
    for ( i = 0; i < nindiv; i++ ){
      Eavrg = Eavrg + energy[i];
    }
    Eavrg = Eavrg/nindiv;
    if( Eavrg == EavrgPrev ){ 
      popstg = popstg + 1;
    }
    else popstg = ZERO;
    EavrgPrev = Eavrg;

    if ( fabs(energy[0] - energy[nindiv-1]) > absenrgy &&  popstg == 5 ){
      ReoptPop(prog, gamesspath, xyzname, datname, outname, rstname, natoms, nindiv, nstep, 1, gmsver, ncore, energy, XYZ, cnct, zmat, MPop);
      genconv = 0;
    }
    
    tindiv = 2*nindiv;
    Order(natoms, tindiv, energy, conv, MPop);
    
    /* Counting Generations before convergence */
    if ( fabs(energy[0] - energy[nindiv-1]) < absenrgy ) genconv = genconv + 1;
    else genconv = 0;

    /* Print the Generation's Energy */
    printf("\n\n***  Generation's [%d] Energy (Hartree):\n ==============================================\n\n", gen+1);
    for ( i = 0 ; i < nindiv; i++){
      if ( conv[i] == 1) printf("   E. Indiv.[%d] = %lf *\n", i+1, energy[i]);
      else printf("   E. Indiv.[%d] = %lf\n", i+1, energy[i]);
    }
    
    /* Population's Energy Analysis */ 
    Eavrg = ZERO;
    for ( i = 0; i < nindiv; i++ ){
      Eavrg = Eavrg + energy[i];
      if ( energy[i] < Emin || i == 0 ) Emin = energy[i];
      if ( energy[i] != ZERO && energy[i] > Emax || i == 0 ) Emax = energy[i];
    }
    Eavrg = Eavrg/nindiv;
    
    printf("\n  --------------------------- \n\n");
    printf("   Minimum Electronic Energy = %10.6lf\n   Average Electronic Energy = %10.6lf\n   Maximum Electronic Energy = %10.6lf\n", Emin, Eavrg, Emax); 
    printf("\n ==============================================\n"); 

    /* Counting the number of Generations */
    
    gen = gen + 1;
    //printf("gen=%d\n",gen);
    /* Backup of Final Population */
    if ( strncasecmp(backup,"off", 6) ){
      typga = 1;
      PopBackup(natoms, nindiv, typga, gen, energy, MPop, XYZ);
    }

    /* Stoping Algorithm Test */
    if( gen == ngen ) break;
    
  }while( genconv != genafconv );
  unlink(fname);
  
  
  /*---------------------------------FINAL ANALYSIS--------------------------------------*/
  printf("\n\n\n **************************************************************");
  printf("\n\n *  *  *  *  *  *  *  ALGORITHM IS CONVERGED *  *  *  *  *  *  \n");
  printf("\n **************************************************************\n");
  printf("\n\n\n ----------------------- FINAL ANALYSIS -----------------------\n\n");

  printf(" ** Algorithm Analysis: \n");
  printf(" =========================\n");
  printf(" N° Generations Needed: %d\n", gen);
  printf(" Total Population:      %d\n\n", nindiv);
  printf(" --------------------------------------------------\n");
  printf("  Structure's Energy Founded:  %10.6lf  Hartree\n", energy[0]);
  printf(" --------------------------------------------------\n");

  //printf("\n\n ** Structure Founded ** \n =====================\n");
  //printf("\n       *Energy*\t\t*No. of Structures*\n");
  //printf("    %12.6lf\t        %d\n", energy[0], nindiv);
  
  printf("\n ** Molecular Structures in Cartesian Coordinates **\n ===========================================================\n");
  for ( i = 0; i < nindiv; i++ ) {
    ip = i*natoms;
    
    printf(" %d\n", natoms);/**/
    printf(" Electronic Energy (Hartree) = %10.10lf\n", energy[i]);/**/
    
    for(j = ip; j < (i + 1)*natoms; j++ ){
      printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, MPop[j].xyz->x, MPop[j].xyz->y, MPop[j].xyz->z);
    }
  }
  printf("\n ===========================================================\n");
  printf("\n  Done!\n");
}

