#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "predator.h"
#include "order.h"

#define ZERO          0.0E+0

int MNGA_Predator(int natoms, int nindiv, int *conv, double absenrgy, double *energy, defmol *MPop)
{
  int    i, j;
  int    neq, elim;
  double enerc;
  
  /*Ordering Population*/
  Order(natoms, nindiv, energy, conv, MPop); 
  
  /* Eliminating Equal Structures */
  for ( i = 0 ; i < nindiv ; i++){

    for ( j = 0 ; j < nindiv ; j++){
      if ( fabs(energy[i+(3*nindiv+1)] - energy[j]) < absenrgy ){
	energy[j] = ZERO; 
      }
    }

  }

  /*Re-Ordering Population*/
  Order(natoms, nindiv, energy, conv, MPop); 
}


