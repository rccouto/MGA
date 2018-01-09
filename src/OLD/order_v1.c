#include <stdio.h>
#include <stdlib.h>

#include "defmol.h"
#include "order.h"

void Order(int nindiv, double *energy, defmol *MPop)
{
  int i, j, k;
  double temp, temp2;
  
  for ( i = 1 ; i < nindiv ; i++ ){
    for ( j = 0 ; j < i ; j++ ){
      if ( energy[j] > energy[i] ){
	temp = energy[j];
	temp2 = MPop[j].xyz; 
	
	energy[j] = energy[i];
	MPop[j].xyz = MPop[i].xyz;
	
	for ( k = i ; k > j ; k-- ){
	  energy[k] = energy[k - 1];
	  MPop[k].xyz = MPop[k - 1].xyz;
	}
	energy[k + 1] = temp;
	MPop[k + 1].xyz = temp2;
      }
    }
  }  
}
