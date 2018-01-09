#include <stdio.h>
#include <stdlib.h>

#include "defmol.h"
#include "alloc.h"
#include "ordermolx.h"

/* Order only the molecular skeleton atoms*/

void OrderMolX(int natoms, int nonhydro, int m1, int m2, defmol *MPop)
{
  int   i, j, k, i1, j1, k1, i2, j2, k2;
  char  temp1, temp2;
  Txyz  *XYZ_temp;  
  
  XYZ_temp = (Txyz *) calloc( natoms, sizeof(Txyz) );

  for ( i = 1 ; i < (nonhydro - 1) ; i++ ){
    for ( j = 0 ; j < i ; j++ ){
      i1 = m1 + i;
      j1 = m1 + j;
      i2 = m2 + i;
      j2 = m2 + j;
       if ( MPop[j1].xyz->x > MPop[i1].xyz->x ){

	temp1 = *MPop[j1].symb;
	XYZ_temp[0].x = MPop[j1].xyz->x;
	XYZ_temp[0].y = MPop[j1].xyz->y;
	XYZ_temp[0].z = MPop[j1].xyz->z;
	
	temp2 = *MPop[j2].symb;
	XYZ_temp[1].x = MPop[j2].xyz->x;
	XYZ_temp[1].y = MPop[j2].xyz->y;
	XYZ_temp[1].z = MPop[j2].xyz->z;

	*MPop[j1].symb = *MPop[i1].symb;
	MPop[j1].xyz->x = MPop[i1].xyz->x;
	MPop[j1].xyz->y = MPop[i1].xyz->y;
	MPop[j1].xyz->z = MPop[i1].xyz->z;
	
	*MPop[j2].symb = *MPop[i2].symb;
	MPop[j2].xyz->x = MPop[i2].xyz->x;
	MPop[j2].xyz->y = MPop[i2].xyz->y;
	MPop[j2].xyz->z = MPop[i2].xyz->z;

	for ( k = i ; k > j ; k-- ){
	  k1 = m1 + k;
	  k2 = m2 + k;
	  *MPop[k1].symb = *MPop[k1 - 1].symb;
	  MPop[k1].xyz->x = MPop[k1 - 1].xyz->x;
	  MPop[k1].xyz->y = MPop[k1 - 1].xyz->y;
	  MPop[k1].xyz->z = MPop[k1 - 1].xyz->z;

	  *MPop[k2].symb = *MPop[k2 - 1].symb;
	  MPop[k2].xyz->x = MPop[k2 - 1].xyz->x;
	  MPop[k2].xyz->y = MPop[k2 - 1].xyz->y;
	  MPop[k2].xyz->z = MPop[k2 - 1].xyz->z;
	} 
	k1 = m1 + k;
	k2 = m2 + k;
	*MPop[k1 + 1].symb = temp1;
	MPop[k1 + 1].xyz->x = XYZ_temp[0].x; 
	MPop[k1 + 1].xyz->y = XYZ_temp[0].y;
	MPop[k1 + 1].xyz->z = XYZ_temp[0].z;

	*MPop[k2 + 1].symb = temp2;
	MPop[k2 + 1].xyz->x = XYZ_temp[1].x; 
	MPop[k2 + 1].xyz->y = XYZ_temp[1].y;
	MPop[k2 + 1].xyz->z = XYZ_temp[1].z;
      }
    }
  }
}

