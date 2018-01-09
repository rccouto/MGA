#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "mnga_history.h"
#include "orderfpop.h"

int MNGA_History(int pos, int nindiv, int natoms, int *conv, double *energy, double absenrgy, defmol *MPop, Txyz *XYZ)
{
  int     i, j, newstruct, nd, *pconv, ip, nwmol;
  double *penergy;
  defmol *pMPop;
  
  nd = pos*natoms;
  newstruct = 3;

  /* Adding Structure to Backup Population */  
  if ( energy[pos] < 0){ 
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
