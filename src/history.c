#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "history.h"
#include "orderfpop.h"
#include "confirmconv.h"

#define INF   999999999

int History(char *prog, char *gamesspath, int pos, int nindiv, int natoms, int nstep, int *conv, int gmsver, int ncore, double *energy, double absenrgy, defmol *MPop, Txyz *XYZ, Tcnct *cnct, Tzmat *zmat)
{
  int     i, j, k, newstruct, nd, *pconv, ip, nwmol;
  double *penergy;
  defmol *pMPop;
  
  nd = pos*natoms;
  newstruct = 3;
  
  // FOR GAMESS
  if ( strncasecmp(prog,"OBabel", 6) ){
    
    if( energy[pos] < 0 ){
      for ( i = 2*nindiv; i < 3*nindiv; i++){
	if ( fabs(energy[pos] - energy[i]) < absenrgy ) {
	  newstruct = 1;
	}
      }
      if (newstruct != 1) ConfirmConv(prog, gamesspath, pos, natoms, nindiv, nstep, conv, gmsver, ncore, energy, MPop, XYZ, cnct, zmat);
    }
    
    newstruct = 3;
    if ( conv[pos] != 1 && energy[pos] < 0){ 
      for ( i = 2*nindiv; i < 3*nindiv; i++){
	if ( fabs(energy[pos] - energy[i]) < absenrgy ) {
	  newstruct = 1;
	}
      }
      if (newstruct != 1) newstruct = 2;
    }
    
    if ( newstruct == 2 ) {
      nwmol=3*nindiv*natoms;
      energy[3*nindiv] = energy[pos];
      for ( j = 0; j < natoms; j++ ){
	XYZ[nwmol + j].x = XYZ[nd + j].x;
	XYZ[nwmol + j].y = XYZ[nd + j].y;
	XYZ[nwmol + j].z = XYZ[nd + j].z;
      }  
      OrderFPop(natoms, nindiv, energy, conv, XYZ);
    }
  }

  // FOR OPEN BABEL
  else{
    
    // Compare and Optimize
    if( energy[pos] < INF ){
      for ( i = 2*nindiv; i < 3*nindiv; i++){
	if ( fabs(energy[pos] - energy[i]) < absenrgy ) {
	  newstruct = 1;
	}
      }
      if (newstruct != 1){
	for( k = 0; i < 4; i++ ){
	  ConfirmConv(prog, gamesspath, pos, natoms, nindiv, nstep, conv, gmsver, ncore, energy, MPop, XYZ, cnct, zmat);
	}
      }
    }
    
    newstruct = 3;
    // Compare
    if ( energy[pos] < INF ){ 
      for ( i = 2*nindiv; i < 3*nindiv; i++){
	if ( fabs(energy[pos] - energy[i]) < absenrgy ) {
	  newstruct = 1;
	}
      }
      if (newstruct != 1) newstruct = 2;
    }

    if ( newstruct == 2 ) {
      nwmol=3*nindiv*natoms;
      energy[3*nindiv] = energy[pos];
      for ( j = 0; j < natoms; j++ ){
	XYZ[nwmol + j].x = XYZ[nd + j].x;
	XYZ[nwmol + j].y = XYZ[nd + j].y;
	XYZ[nwmol + j].z = XYZ[nd + j].z;
      }  
      OrderFPop(natoms, nindiv, energy, conv, XYZ);
    }
  }
}  
