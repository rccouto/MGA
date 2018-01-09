#include <stdio.h>
#include <stdlib.h>

#include "defmol.h"
#include "geninp.h" 

int GenInp ( const char flag, const int nat, const int nbs, const int mem,
	     const int timl, const int intee, const int ncv, const int *icharg,
	     const int *mult, const char *fname, const char *prog, 
	     const char *loc, const char *bs, const char *psdpot, 
	     const char *rstname, defmol *molec, const int nstep)

{
  char   rdaux[8];
  int    j, st;
  defmol *pmolec;
  FILE   *fl, *flrst;
  
  fl = fopen(fname,"w");
  
  pmolec = &*molec;

  if ( !strncasecmp(prog,"DALTON", 6) ){
    fprintf(fl,"Not implemented, yet.");
  }
  
  else if ( !strncasecmp(prog,"GAMESS", 6) ){ 
    if ( flag == 'O' ){
      if ( !strncasecmp(loc,"DFT", 3) )
	fprintf(fl," $CONTRL  SCFTYP=%c%c%c MULT=%d ICHARG=%d RUNTYP=OPTIMIZE MAXIT=%d  COORD=ZMT $END\n", *(loc + 4), *(loc + 5), *(loc + 6), mult, icharg, intee);
      else 
	fprintf(fl," $CONTRL  SCFTYP=%s MULT=%d ICHARG=%d RUNTYP=OPTIMIZE MAXIT=%d COORD=ZMT $END\n", loc, mult, icharg, intee);
      if ( !strncasecmp(bs,"PC1", 3) || !strncasecmp(bs,"PC2", 3) ) fprintf(fl," $CONTRL  ISPHER=1 $END\n");
      fprintf(fl," $STATPT  OPTTOL=1.0E-6 NSTEP=%d $END\n", nstep);
      fprintf(fl," $GUESS   GUESS=HUCKEL   $END\n");
    }
    
    else if ( flag == 'E' ){
      if ( !strncasecmp(loc,"DFT", 3) ){
	fprintf(fl," $CONTRL  SCFTYP=%c%c%c MULT=%d ICHARG=%d RUNTYP=ENERGY MAXIT=%d UNITS=ANGS $END\n", *(loc + 4), *(loc + 5), *(loc + 6), mult, icharg, intee);
	if ( !strncasecmp(bs,"PC1", 3) || !strncasecmp(bs,"PC2", 3))fprintf(fl," $CONTRL  ISPHER=1 $END\n");
      }
      else
	fprintf(fl," $CONTRL  SCFTYP=%s MULT=%d ICHARG=%d RUNTYP=ENERGY MAXIT=%d UNITS=ANGS $END\n", loc, mult, icharg, intee);
      fprintf(fl," $STATPT  OPTTOL=1.0E-6 NSTEP=%d $END\n", nstep);
      fprintf(fl," $GUESS   GUESS=HUCKEL   $END\n");
    }
      
    else if (flag == 'G' )
      if ( !strncasecmp(loc,"DFT", 3) )
	fprintf(fl," $CONTRL SCFTYP=%c%c%c MULT=%d ICHARG=%d RUNTYP=GRADIENT MAXIT=%i UNITS=BOHR $END\n", loc[4],loc[5],loc[6], mult, icharg, intee);
      else
	fprintf(fl," $CONTRL  SCFTYP=%s MULT=%d ICHARG=%d RUNTYP=GRADIENT MAXIT=%i UNITS=BOHR $END\n", loc, mult, icharg, intee);
    else if ( flag == 'H' )
      if ( !strncasecmp(loc,"DFT", 3) )
	fprintf(fl," $CONTRL  SCFTYP=%c%c%c MULT=%d ICHARG=%d RUNTYP=OPTIMIZE MAXIT=%i $END\n", loc[5],loc[6],loc[7], mult, icharg, intee);
      else
	fprintf(fl," $CONTRL  SCFTYP=%s MULT=%d ICHARG=%d RUNTYP=OPTIMIZE MAXIT=%d $END\n", loc, mult, icharg, intee); 
    else {
      perror("Error [GenInp]: Wrong flag inside GenInp");
      exit(EXIT_FAILURE);
    }
    
    
    fprintf(fl," $SCF     DIRSCF=.T. FDIFF=.F. NCONV=%d SOSCF=.F. DAMP=.T. SHIFT=.T.   $END\n",ncv);
    fprintf(fl," $SYSTEM  TIMLIM=%d MEMORY=%d   $END\n", mem, timl);  
    fprintf(fl," $FORCE   PROJCT=.TRUE.   $END\n");
    
    if ( !strncasecmp(bs,"6-31G**", 7) )
      fprintf(fl," $BASIS   GBASIS=N31 NGAUSS=6 NPFUNC=1 NDFUNC=1  $END\n");
    else if ( !strncasecmp(bs,"6-31G*", 6) )
      fprintf(fl," $BASIS   GBASIS=N31 NGAUSS=6 NPFUNC=1 NDFUNC=0  $END\n");
    else if ( !strncasecmp(bs,"6-31G", 5) )
      fprintf(fl," $BASIS   GBASIS=N31 NGAUSS=6 $END\n");
    else if ( !strncasecmp(bs,"6-311G**", 8) )
      fprintf(fl," $BASIS   GBASIS=N311 NGAUSS=6 NPFUNC=1 NDFUNC=1  $END\n");
    else if ( !strncasecmp(bs,"6-311G*", 7) )
      fprintf(fl," $BASIS   GBASIS=N311 NGAUSS=6 NPFUNC=1 NDFUNC=0  $END\n");
    else if ( !strncasecmp(bs,"6-311G", 6) )
      fprintf(fl," $BASIS   GBASIS=N311 NGAUSS=6 $END\n");
    else if ( !strncasecmp(bs,"4-311G**", 6) )
      fprintf(fl," $BASIS   GBASIS=N31 NGAUSS=4 NPFUNC=1 NDFUNC=1  $END\n");
    else if ( !strncasecmp(bs,"4-31G*", 6) )
      fprintf(fl," $BASIS   GBASIS=N31 NGAUSS=4 NPFUNC=1 NDFUNC=0  $END\n");
    else if ( !strncasecmp(bs,"4-31G", 5) )
      fprintf(fl," $BASIS   GBASIS=N31 NGAUSS=4 $END\n");
    else if ( !strncasecmp(bs,"4-311G**", 8) )
      fprintf(fl," $BASIS   GBASIS=N311 NGAUSS=4 NPFUNC=1 NDFUNC=1  $END\n");
    else if ( !strncasecmp(bs,"4-311G*", 7) )
      fprintf(fl," $BASIS   GBASIS=N311 NGAUSS=4 NPFUNC=1 NDFUNC=0  $END\n");
    else if ( !strncasecmp(bs,"4-311G", 6) )
      fprintf(fl," $BASIS   GBASIS=N311 NGAUSS=4 $END\n");
    
    else
      fprintf(fl," $BASIS   GBASIS=%s   $END\n", bs);
    
    if ( psdpot != NULL && !strncasecmp(loc,"DFT", 3) ) 
      if ( flag == 'I' )
	fprintf(fl," $DFT     DFTTYP=%s   $END\n", psdpot);
      else 
	fprintf(fl," $DFT     DFTTYP=%s   $END\n", psdpot);
    /* fprintf(fl," $DFT     DFTTYP=%s NRAD0=96 NTHE0=12 NPHI0=24  $END\n", psdpot); */
    
    
    if ( flag == 'G' || flag == 'H' ){
      if( flrst = fopen(rstname,"r") ) {
	fscanf(flrst, "%s %d", &rdaux, &nbs);
	/* printf("\n\nrdaux = %s, nbs = %d\n\n",rdaux, nbs); */
	fclose(flrst);
      }
      
      fprintf(fl," $GUESS   GUESS=MOREAD NORB=%i   $END \n", nbs);
      /* fprintf(fl," $GUESS   GUESS=HUCKEL   $END\n"); */
      /* fprintf(fl," $STATPT  NSTEP=0 PROJCT=.TRUE. HESS=CALC   $END \n"); */
      
      st = 0;
      for (j = 0; j < nat; j++) {
	if ( pmolec->atmn == 0 ) {
	  if (st == 0 )fprintf(fl," $STATPT  IFREEZ(1)=");
	  else fprintf(fl, ","); 
	  st = 1;
	  fprintf(fl, "%i,%i,%i", j*3 + 1, j*3 + 2, j*3 + 3);
	}
	pmolec++;
      }
      if ( st == 1 )fprintf(fl,"   $END\n");
    }
  }
  
  else{
    printf("Error [GenInp]: eeecp \"%s\" is not available!", prog);
    exit(EXIT_FAILURE);
  }
  fclose(fl); 
  return 0;
}
