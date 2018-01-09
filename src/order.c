#include <stdio.h>
#include <stdlib.h>

#include "defmol.h"
#include "order.h"
#include "alloc.h"

void Order(int natoms, int nindiv, double *energy, int *conv, defmol *MPop)
{
  int    i, j, k, l;
  int    ni, nj, nk, nks, nkm;
  int    temp2;
  double temp;
  Txyz   *XYZ_temp;

  XYZ_temp = (Txyz *) calloc( natoms, sizeof(Txyz) );

  for ( i = 1 ; i < nindiv ; i++ ){
    for ( j = 0 ; j < i ; j++ ){
      if ( energy[j] > energy[i] ){
	temp = energy[j];
	temp2 =  conv[j];

	nj = j * natoms;
	for ( l = 0 ; l < natoms ; l++ ){
	  XYZ_temp[l].x = MPop[nj + l].xyz->x;
	  XYZ_temp[l].y = MPop[nj + l].xyz->y;
	  XYZ_temp[l].z = MPop[nj + l].xyz->z;
	}

	energy[j] = energy[i];
	  conv[j] =   conv[i];

	ni = i * natoms;
	for ( l = 0 ; l < natoms ; l++ ){
	  MPop[nj + l].xyz->x = MPop[ni + l].xyz->x;
	  MPop[nj + l].xyz->y = MPop[ni + l].xyz->y;
	  MPop[nj + l].xyz->z = MPop[ni + l].xyz->z;
	}

	for ( k = i ; k > j ; k-- ){
	
	  energy[k] = energy[k - 1];
	    conv[k] =   conv[k - 1];

	  nk = k * natoms;
	  nkm = (k - 1) * natoms;
	  for ( l = 0 ; l < natoms ; l++ ){
	    MPop[nk + l].xyz->x = MPop[nkm + l].xyz->x;
	    MPop[nk + l].xyz->y = MPop[nkm + l].xyz->y;
	    MPop[nk + l].xyz->z = MPop[nkm + l].xyz->z;
	    }
	}
	
	energy[k + 1] = temp;
	  conv[k + 1] = temp2;

	nks = (k + 1) * natoms;
	for ( l = 0 ; l < natoms ; l++ ){
	  MPop[nks + l].xyz->x = XYZ_temp[l].x; 
	  MPop[nks + l].xyz->y = XYZ_temp[l].y;
	  MPop[nks + l].xyz->z = XYZ_temp[l].z;
	}
      }
    }
  }  
}
