#include <stdlib.h>
#include <stdio.h>

void Order(int natoms, int nindiv, double *energy, double *energy_old, defmol *MPop, defmol *MPop_old)
{
  int i, j, k;
  double temp;

  //*************************Insertion Sorting******************************
  
  for ( i = 1 ; i < nindiv; i++ ){
    for ( j = 0 ; j < i ; j++ ){
      if ( energy[j] > energy[i] ){
	
	temp = energy[j];
	
	energy[j] = energy[i];
	
	for ( k = i ; k > j ; k-- ){
	  
	  energy[k] = energy[k - 1];
	}
	energy[k + 1] = temp;
      }
    }
  }
  
  //*****************************Order**************************************
  
  for( i = 0; i < nindiv; i++ ){
    for( j = 0; j < nindiv; j++ ){
      if( energy_old[i] == energy[j] ){
	for( k = 0; k < natoms; k++ ){
	  MPop[ ( natoms * j ) + k ].xyz = MPop_old[ ( natoms * i ) + k ].xyz;
	}
      }
    }
  }
}
