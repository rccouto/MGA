#include <stdio.h>
#include <stdlib.h>

#include "defmol.h"
#include "orderfpop.h"
#include "alloc.h"

void OrderFPop(int natoms, int nindiv, double *energy, int *conv, Txyz *XYZ)
{
  int    i, j, k, l;
  int    ni, nj, nk, nks, nkm;
  int    temp2;
  double temp;
  Txyz   *XYZ_temp;

  XYZ_temp = (Txyz *) calloc( 3*natoms, sizeof(Txyz) );

  for ( i = 2*nindiv + 1; i <= 3*nindiv; i++ ){
    for ( j = 2*nindiv ; j < i ; j++ ){
      if ( energy[j] > energy[i] ){
	temp = energy[j];
	temp2 =  conv[j];

	nj = j * natoms;
	for ( l = 0 ; l < natoms ; l++ ){
	  XYZ_temp[l].x = XYZ[nj + l].x;
	  XYZ_temp[l].y = XYZ[nj + l].y;
	  XYZ_temp[l].z = XYZ[nj + l].z;
	}

	energy[j] = energy[i];
	  conv[j] =   conv[i];

	ni = i * natoms;
	for ( l = 0 ; l < natoms ; l++ ){
	  XYZ[nj + l].x = XYZ[ni + l].x;
	  XYZ[nj + l].y = XYZ[ni + l].y;
	  XYZ[nj + l].z = XYZ[ni + l].z;
	}

	for ( k = i ; k > j ; k-- ){
	
	  energy[k] = energy[k - 1];
	    conv[k] =   conv[k - 1];

	  nk = k * natoms;
	  nkm = (k - 1) * natoms;
	  for ( l = 0 ; l < natoms ; l++ ){
	    XYZ[nk + l].x = XYZ[nkm + l].x;
	    XYZ[nk + l].y = XYZ[nkm + l].y;
	    XYZ[nk + l].z = XYZ[nkm + l].z;
	    }
	}
	
	energy[k + 1] = temp;
	  conv[k + 1] = temp2;

	nks = (k + 1) * natoms;
	for ( l = 0 ; l < natoms ; l++ ){
	  XYZ[nks + l].x = XYZ_temp[l].x; 
	  XYZ[nks + l].y = XYZ_temp[l].y;
	  XYZ[nks + l].z = XYZ_temp[l].z;
	}
      }
    }
  }  
}
