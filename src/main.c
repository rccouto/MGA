#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "defmol.h"
#include "genindv.h"
#include "geninp.h"
#include "genxyz.h"
#include "runprog.h"
#include "getxyz.h"
#include "getxyzga.h"
#include "getxyzobabel.h"
#include "getenergy.h"
#include "getenergyobabel.h"
#include "order.h"
#include "rdzmat.h"
#include "mutation.h"
#include "rtgroup.h"
#include "distmeasure.h"
#include "nga.h"
#include "iga.h"

#define ZERO          0.0E+0
#define ONE           1.0E+0
#define TIMLIM        15000000
#define MEMORY        33554432
#define NCLTSTEP      1.0E-2
#define MAXINT        200
#define CHARMAX       256

#define HARTREE2KCAL  27.2114E+0 * 26.03E+0 
#define BOHR2ANGSTRON 1.88972613392E+0

#define DALTONPATH "/usr/local/bin"
#define GAMESSPATH "/usr/local/bin/rungms"

double pi = 3.1415926535897932384626433832795028E+0;
double hartree2kcal = HARTREE2KCAL;
double bohr2angstron = BOHR2ANGSTRON;

int    natoms = ZERO, nos = ZERO, nindiv=15, ndiv, percmut = 30;
int    timl = TIMLIM, mem = MEMORY, ncv = 8, nstep = 500, *conv = 0;
int    mutmin = 180, mutmax = 360, typmat = 1, typga = 1, ngen = 20, gen = 0;
int    genafconv, nrotgp, natmrot, ncore, gmsver;
double *energy, *Fenergy;

char   *daltonpath = NULL, *gamesspath = NULL, *obabelpath = NULL;
char   *toc = NULL, *title = NULL, *potential = NULL, *backup = "off";
char   *print_pop = "off", *rmdatfile= "on", *prog = NULL, *rstbackup = "off";
char   *optstruct = "on", *rmoutfile = "on", *typcalc = "O", *flag = "R";
char   *loc = NULL, *bs = NULL, *eeecp = NULL, *psdpot = NULL, *eff = NULL;
int    *noas = NULL, *mltpl = NULL, *icharg = NULL;
double absenrgy = 1.0E-06;
long   nseed=1;


Txyz   *XYZ  = NULL;
defmol *MPop = NULL;
defmol *FPop = NULL;
Tzmat  *zmat = NULL;
Tcnct  *cnct = NULL;
Trtt   *rtt  = NULL;
Tatom  *rttatm  = NULL;

main() /* Main routine of MGA program*/
{
  char    path, fname[10], rstname[10], xyzname[15];
  char    datname[15], outname[15], vecname[15], nbsname[15];
  int     i = ZERO, j = ZERO, ij = ZERO, ip, tindiv, *pconv;
  int     nat, nbs, intee, k, nenerg, nonhydro; 
  int     halfindiv, halfgen, gentest;
  double  cf = 1.0E+0, Emin, Emax, Eavrg, EavrgPrev = ZERO, enerc, diff, *penergy;
  time_t  start, end, rawtime;

  Txyz   *pXYZ;
  Tzmat  *pzmat;
  defmol *pMPop;
  defmol *indv;

  FILE   *fltest, *fltest2, *fltest3;
  char   Infile [] = "input.mga";

  time (&start);
  printf("Initialized at %s\n", ctime (&start) );

  printf("\t***********************************************************\n");
  printf("\t**                      Program MGA                      **\n");
  printf("\t**             (Molecular Genetic Algorithm)             **\n");
  printf("\t**                                                       **\n");
  printf("\t**                                                       **\n");
  printf("\t** Principal Author:                                     **\n");
  printf("\t**        Freddy F. Guimaraes                            **\n");
  printf("\t**        email: freddy@quimica.ufg.br                   **\n");
  printf("\t**                                                       **\n");
  printf("\t** Authors:                                              **\n");
  printf("\t**        Rafael C. Couto, Fabrício S. Paranhos.         **\n");
  printf("\t**                                                       **\n");
  printf("\t***********************************************************\n\v");

  printf("Reading input file: |");
  printf("#");fflush(stdout);
  RdInput(Infile, CHARMAX);
  printf("|\n");
 
  /* Printing Input Description */
  printf("\n **  Input description  **  \n =========================\n");
  printf("  * Main *\n");
  if (title != NULL ) printf("    Title: %s\n", title);
  printf("    ABS Energy Comparison: %lf;\n", absenrgy);
  printf("    NSTEP: %d;\n", nstep);
  printf("    Optimize Structures: %s;\n", optstruct);
  printf("    Remove .OUT Files:   %s;\n", rmoutfile);
  printf("    Remove .DAT Files:   %s;\n", rmdatfile);
  printf("    Print Init. Pop.:    %s;\n", print_pop);
  if (gamesspath != NULL )
    printf("    Gamess path: %s;\n", gamesspath);
  else if (daltonpath != NULL )
    printf("    Dalton path: %s;\n", daltonpath);
  else if (obabelpath != NULL )
    printf("    Open Babel path: %s;\n", obabelpath);
  switch(typga){
  case 1: /* Non-Inclusive Genetic Algorithm (NGA) */
    printf("  * Genetic Algorithm *\n");
    printf("    Type of GA: Non-Inclusive Genetic Algorithm;\n");
    break;
  case 2: /* Inclusive Genetic Algorithm (IGA) */
    printf("  * Genetic Algorithm *\n");
    printf("    Type of GA: Inclusive Genetic Algorithm;\n");
    break;
  case 4: /* Random Search */
    printf("  * Random Search *\n");
  }
  if (typga == 4){
    printf("    Population Size: %d;\n", nindiv);
    printf("    Seed: %d;\n", nseed);
  }
  else{
    printf("    Population Size: %d;\n", nindiv);
    printf("    Seed: %d;\n", nseed);
    printf("    Generations: %d;\n", ngen);
    printf("    Generations After Convergence: %d;\n", genafconv);
    switch(typmat){
    case 1: /* Single Cut */
      printf("    Type of Mating: Single Cut;\n");
      break;
    case 2: /* Double Cut */
      printf("    Type of Mating: Double Cut;\n");
      break;
    case 3: /* Mixed Cut */
      printf("    Type of Mating: Mixed Cut;\n");
      break;
    }
    printf("    Percentage of Mutation: %d;\n", percmut);
    printf("    Rotation Mutation Range: %d - %d;\n", mutmin, mutmax);
    printf("    Number of Rotation Groups : %d;\n", nrotgp);
  }
  printf("  * Molecule *\n");
  printf("    Total number of atoms: %d;\n",natoms);
  printf("  * Force Field *\n");
  printf("    Interatomic Potential: %s;\n", potential);
  if ( strncasecmp(potential,"Empiric", 7) ){
    printf("    Levels Of Calculation: %s;\n", loc);
    printf("    Basis Set: %s;\n", bs);
    printf("    Multiplicity: %d;\n", mltpl);
    printf("    Molecular Charge: %d;\n", icharg);
  }
  if ( !strncasecmp(potential,"Empiric", 7) ){
    printf("    Empiric Force Field: %s;\n", eff);
  }
  printf("    Electronic Structure Package: %s;\n", prog);
  if(ncore > 1) printf("    GAMESS Parallel Computing Cores: %d;\n", ncore);
  ij = 0;

  if (typga != 4 ){ 
    nindiv = 2*nindiv;
  }
  /* Generate Average Molecule */ 
  flag = "A";
  GenIndv(*flag, Infile, cf, natoms, &nseed, XYZ, MPop);
  
  printf("  * Cartesian Coordinates of Average Initial Molecule *\n");
  for(i = 0; i < natoms; i++ ) {
    printf("    %s  %14.8lf  %14.8lf  %14.8lf \n", MPop[i].symb, XYZ[i].x, XYZ[i].y, XYZ[i].z);
  }
  
  // Used in Backup
  if ( strncasecmp(rstbackup,"on", 6) ){
    
    /* Generate Initial Population */
    printf("\n\n**  Generating Initial Population: |");
    cf = 1.0E+0;
    
    for ( i = 0; i < nindiv; i++ ) {
      printf("#");fflush(stdout);
      ip = i*natoms;
      pXYZ  = &XYZ[ ip ];
      pMPop = &MPop[ ip ];
      
      flag = "R";
      printf("#");fflush(stdout);
      GenIndv(*flag, Infile, cf, natoms, &nseed, pXYZ, pMPop);
      printf("#");fflush(stdout);
      
      /* Applying Rotation Operator */
      RtGroup(natoms, nrotgp, &nseed, rtt, rttatm, pXYZ);
    }
    printf("|"); ;fflush(stdout);/**/
    
    /* Print the Non-optimized Structures for the Initial Population */
    if ( !strncasecmp(print_pop,"on", 2) ){
      printf("\n\n**  Molecular Non-optimized Structures for the Initial Population:\n");
      for ( i = 0; i < nindiv; i++ ) {
	ip = i*natoms;
	pMPop = &MPop[ip];
	printf("%d\n\n", natoms);
	for(j = ip; j < (i + 1)*natoms; j++ ){
	  printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, MPop[j].xyz->x, MPop[j].xyz->y, MPop[j].xyz->z);
	}
      }
    }
    
    /* Generate  Input for Eletronic Structure Calculations */
    
    if ( !strncasecmp(optstruct,"on", 2) )printf("\n\n**  Generating Input for Eletronic Structure Calculations: |");
    else printf("\n\n**  Generating Input for Energy Calculations: |");
    
    // FOR GAMESS
    if ( strncasecmp(prog,"OBabel", 6) ){
      
      pMPop = &MPop[0];
      nbs = 0;
      nat = 1;
      intee = 200;
      //strcpy(prog,    "GAMESS");
      strcpy(fname,   "Gamess.opt");
      strcpy(nbsname, "Gamess.nbs");
      strcpy(rstname, "Gamess.dat");
      
      if ( !strncasecmp(optstruct,"off", 3) ) typcalc = "E";
      
      printf("#");fflush(stdout);
      GenInp (*typcalc, natoms, nbs, mem, timl, intee, ncv, icharg, mltpl, fname, prog, 
	      loc, bs, psdpot, rstname, pMPop, nstep );
      printf("#");fflush(stdout);
      
      /* Read Z-Matrix Connections */
      printf("#");fflush(stdout);
      RdZMat(Infile, natoms, cnct);
      printf("#");fflush(stdout);
      
      for ( i = 0; i < nindiv; i++ ) {
	
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
    }
    
    // FOR OPEN BABEL
    else{
      for ( i = 0; i < nindiv; i++ ) {
	
	/* Read Z-Matrix Connections */
	printf("#");fflush(stdout);
	RdZMat(Infile, natoms, cnct);
	printf("#");fflush(stdout);
	
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
    }
    
    if ( !strncasecmp(optstruct,"on", 2) )printf("\n\n**  Running Eletronic Structure Calculations: |");
    else printf("\n\n**  Running Energy Calculations: |");
    
    // FOR GAMESS
    if ( strncasecmp(prog,"OBabel", 6) ){
      
      for ( i = 0; i < nindiv; i++ ) {
	
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
	  if (typga != 4 ){       
	    printf("#");fflush(stdout);
	    GetXYZGA (datname, rstname, natoms, nstep, penergy, pconv, pMPop);  
	    printf("#");fflush(stdout);
	  }
	  else{
	    printf("#");fflush(stdout);
	    GetXYZ (datname, natoms, penergy, pMPop);  
	    printf("#");fflush(stdout);
	  }
	}
	else{
	  printf("#");fflush(stdout);
	  GetEnergy (outname, natoms, penergy, pMPop);  
	  printf("#");fflush(stdout);
	}
	
	/* Remove .dat Files */
	if ( !strncasecmp(rmdatfile,"on", 2 ) ){
	  if ( ! (fltest=fopen(datname,"r")) ) {
	    printf("|\n\nError [Main]: The file \"%s\" is not available for deleting!\n",datname);
	    exit(EXIT_FAILURE);
	  }
	  fclose(fltest);
	  unlink(datname);
	}/**/
	
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
    }
    
    // FOR OPEN BABEL
    else{
      for ( i = 0; i < nindiv; i++ ) {
	
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
      }
      printf("|                                                ");fflush(stdout);
    }
    
    
    if (typga == 2 || typga == 4 && strncasecmp(optstruct,"off", 3) ){
      if ( !strncasecmp(potential,"Empiric", 7) ){
	for (i = 0; i < 4; i++) {
	  ReoptPop(prog, gamesspath, xyzname, datname, outname, rstname, natoms, nindiv, nstep, typga, gmsver, ncore, energy, XYZ, cnct, zmat, MPop);
	}
      }
      else{
	ReoptPop(prog, gamesspath, xyzname, datname, outname, rstname, natoms, nindiv, nstep, typga, gmsver, ncore, energy, XYZ, cnct, zmat, MPop);
      }
    }
    
    /* Confirm if all structures corresponds to Molecule in study */
    for ( i = 0; i < nindiv; i++){
      ip = i*natoms;
      pXYZ  = &XYZ[ip];
      penergy = &energy[i];
      DistMeasure(prog, natoms, pXYZ, penergy, cnct);
    }
    Order(natoms, nindiv, energy, conv, MPop);
  }

  else{ 

    printf("\n\n**  Reading Backup Population: |");
    printf("#");fflush(stdout);
    ReadBackup(natoms, nindiv, typga, energy, XYZ, MPop);
    printf("#");fflush(stdout);
    printf("|");fflush(stdout);
    //printf("gen=%d\n", gen);
    
    pMPop = &MPop[0];
    nbs = 0;
    nat = 1;
    intee = 200;
    //strcpy(prog,    "GAMESS");
    strcpy(fname,   "Gamess.opt");
    strcpy(nbsname, "Gamess.nbs");
    strcpy(rstname, "Gamess.dat");
    
    if ( !strncasecmp(optstruct,"off", 3) ) typcalc = "E";
    
    printf("#");fflush(stdout);
    GenInp (*typcalc, natoms, nbs, mem, timl, intee, ncv, icharg, mltpl, fname, prog, 
	    loc, bs, psdpot, rstname, pMPop, nstep );
    printf("#");fflush(stdout);
   
    /* Read Z-Matrix Connections */
    printf("#");fflush(stdout);
    RdZMat(Infile, natoms, cnct);
    printf("#");fflush(stdout);
    
  }

  /*------------------------- Initiating Genetic Algorithm ------------------------------*/
  
  switch(typga){
  case 1: /* Non-Inclusive Genetic Algorithm (NGA) */
    NGA(Infile, prog, gamesspath, optstruct, rmoutfile, backup, natoms, nindiv, nstep, percmut, typmat, nrotgp, mutmin, mutmax, genafconv, ngen, conv, gmsver, ncore, energy, absenrgy, &nseed, MPop, XYZ, cnct, zmat, rtt, rttatm);
    break;
    
    
  case 2: /* Inclusive Genetic Algorithm (IGA) */
    IGA(Infile, prog, gamesspath, optstruct, rmoutfile, backup, natoms, nindiv, nstep, percmut, typmat, nrotgp, mutmin, mutmax, genafconv, ngen, conv, gmsver, ncore, energy, Fenergy, absenrgy, &nseed, MPop, FPop, XYZ, cnct, zmat, rtt, rttatm);
    break;
  }
  unlink(fname);
  unlink(nbsname);
 
  if ( typga == 4 ){
    /*---------------------------------FINAL ANALYSIS--------------------------------------*/
    Eavrg = 0.0E+0;
    for ( i = 0; i < nindiv; i++ ){
      Eavrg = Eavrg + energy[i];
      if ( energy[i] < Emin || i == 0 ) Emin = energy[i];
      if ( energy[i] != ZERO && energy[i] > Emax || i == 0 ) Emax = energy[i];
    }
    Eavrg = Eavrg/nindiv;
    printf("\n\n ** Energy analysis ** \n =====================\n");
    printf("   Minimum Electronic Energy = %10.10lf\n   Average Electronic Energy = %10.10lf\n   Maximum Electronic Energy = %10.10lf\n", Emin, Eavrg, Emax); 
    
    nenerg = 1;
    enerc = energy[0];
    printf("\n       *Energy*\t\t*No. of Structures*\n");
    for ( i = 1; i <= nindiv; i++){
      if ( fabs(enerc - energy[i]) < absenrgy){
	nenerg = nenerg + 1;
	if ( i == nindiv ) printf("    %12.6lf\t        %d\n", enerc, nenerg);
      }
      else {
	if (enerc != ZERO) printf("    %12.6lf\t        %d\n", enerc, nenerg);
	nenerg = 1;
	enerc = energy[i];
      }
    }
    
    printf("\n ** Molecular Structures in Cartesian Coordinates **\n ===================================================\n");
    for ( i = 0; i < nindiv; i++ ) {
      if (energy[i] != 0){	  
	ip = i*natoms;
	
	printf("%d\n", natoms);
	printf("Electronic Energy (Hartree) = %10.10lf\n", energy[i]);
	
	for(j = ip; j < (i + 1)*natoms; j++ ){
	  printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, MPop[j].xyz->x, MPop[j].xyz->y, MPop[j].xyz->z);
	}
      }
    }
    printf("\nDone!\n");
  }
  
  /* CPU Time Determination */
  time(&end);
  diff = difftime (end,start);
  if( diff < 60 ) printf("\n  Total CPU Time: %.2lf seconds.\n", diff );
  else if( diff < 3600 ) printf("\n  Total CPU Time: %.2lf minutes.\n", diff/60 );
  else if( diff >= 3600 ) printf("\n  Total CPU Time: %.2lf hours.\n", diff/3600 );
  
  time( &rawtime );
  printf("  Finished at %s\n", ctime (&rawtime) );
  
  return(0);
}
