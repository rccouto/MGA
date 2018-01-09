#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "history.h"
#include "orderfpop.h"

int History(char *prog, char *gamesspath, int pos, int nindiv, int natoms, int nstep, int *conv, double *energy, double absenrgy, defmol *MPop, Txyz *XYZ)
{
  int     i, j, newstruct, nd, *pconv, ip, nwmol;
  double *penergy;
  defmol *pMPop;
  
  nd = pos*natoms;
  newstruct = 3;
  
  if( energy[pos] < 0 ) ConfirmConv(prog, gamesspath, pos, natoms, nindiv, nstep, conv, energy, MPop);

  if ( conv[pos] != 1 && energy[pos] < 0){ 
    for ( i = 2*nindiv; i < 3*nindiv; i++){
      if ( fabs(energy[pos] - energy[i]) < absenrgy ) {
	newstruct = 1;
      }
    }
    if (newstruct != 1) newstruct = 2;
  }
  
  if ( newstruct == 2 ) {
    
    /* printf("\nBEFORE PASSING TO FINAL POPULATION\n");
    for(j = 0; j < natoms; j++ ){
      printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[nd+j].x, XYZ[nd+j].y, XYZ[nd+j].z);
      }*/
    nwmol=3*nindiv*natoms;
    energy[3*nindiv] = energy[pos];
    for ( j = 0; j < natoms; j++ ){
      XYZ[nwmol + j].x = XYZ[nd + j].x;
      XYZ[nwmol + j].y = XYZ[nd + j].y;
      XYZ[nwmol + j].z = XYZ[nd + j].z;
    }  
    
    /*printf("IN FINAL POPULATION with XYZ\n");
    for(j = 0; j < natoms; j++ ){
      printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[nwmol+j].x, XYZ[nwmol+j].y, XYZ[nwmol+j].z);
      }*/



    /*printf("\n\n ** BEFORE [OrderFPop] **\n \n");
    for ( i = 2*nindiv; i < 3*nindiv; i++ ) {
      ip = i*natoms;
      printf("%d\n", natoms);
      printf("Electronic Energy (Hartree) = %10.10lf\n", energy[i]);
      for(j = 0; j < natoms; j++ ){
	printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[ip+j].x, XYZ[ip+j].y, XYZ[ip+j].z);
      }
    }*/



    OrderFPop(natoms, nindiv, energy, conv, XYZ);



    /*printf("\n\n ** AFTER [OrderFPop] **\n\n");
    for ( i = 2*nindiv; i < 3*nindiv; i++ ) {
      ip = i*natoms;
      printf("%d\n", natoms);
      printf("Electronic Energy (Hartree) = %10.10lf\n", energy[i]);
      for(j = 0; j < natoms; j++ ){
	printf("\t%s %15.8lf  %15.8lf  %15.8lf\n", &MPop[j].symb, XYZ[ip+j].x, XYZ[ip+j].y, XYZ[ip+j].z);
      }
    }*/


  }
}
