#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "predator.h"
#include "order.h"

#define ZERO          0.0E+0

int Predator(char *prog, int natoms, int nindiv, int *conv, double absenrgy, double *energy, defmol *MPop)
{
  int    i, j;
  int    neq, elim;
  double enerc;
  
  //printf("\nABS=%lf\n", absenrgy);

  // FOR GAMESS
  if ( strncasecmp(prog,"OBabel", 6) ){
    /*Ordering Population*/
    Order(natoms, nindiv, energy, conv, MPop); 
    
    /* Eliminating Equal Structures */
    for ( i = 0 ; i < nindiv ; i++){
      for ( j = i + 1 ; j < nindiv ; j++){
	if ( fabs(energy[i] - energy[j]) < absenrgy ){
	  energy[j] = ZERO; 
	}
      }
    }
    
    /*Re-Ordering Population*/
    Order(natoms, nindiv, energy, conv, MPop); 
  }

  // FOR OPEN BABEL
  else{
    /*Ordering Population*/
    Order(natoms, nindiv, energy, conv, MPop); 
    
    /* Eliminating Equal Structures */
    for ( i = 0 ; i < nindiv ; i++){
      for ( j = i + 1 ; j < nindiv ; j++){
	if ( fabs(energy[i] - energy[j]) < absenrgy ){
	  energy[j] = 9999999999; 
	}
      }
    }
    
    /*Re-Ordering Population*/
    Order(natoms, nindiv, energy, conv, MPop); 
  }
}

