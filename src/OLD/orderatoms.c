#include <stdio.h>
#include <stdlib.h>

#include "defmol.h"
#include "alloc.h"
#include "orderatoms.h"

void OrderAtoms(int natoms, defmol *MPop)
{
  int   i, k, l, nonhydro;
  char  temp;
  Txyz  *XYZ_temp;  

  nonhydro = 0;
  for ( i = 0; i < natoms; i++){
    printf("symb = %s\n", MPop[i].symb); 
    if ( *MPop[i].symb != 'H' ){
      nonhydro = nonhydro + 1;
    }
  }
  printf("\nNumber of non-hydrogen atoms = %d\n", nonhydro);
  
  XYZ_temp = (Txyz *) calloc( natoms, sizeof(Txyz) );
  
  for ( i = 0; i < (natoms-1); i++){
    //    printf("symb = %s\n", MPop[i].symb);
    if ( *MPop[i].symb == 'H' && *MPop[i+1].symb != 'H'  ){
      printf("i = %d\n", i);
      temp = *MPop[i].symb;
      
      XYZ_temp[0].x = MPop[i].xyz->x;
      XYZ_temp[0].y = MPop[i].xyz->y;
      XYZ_temp[0].z = MPop[i].xyz->z;

      printf("MPop[%d].xyz->x = %lf\n", i, MPop[i].xyz->x);
      printf("XYZ_temp[0].x = %lf\n", XYZ_temp[0].x);
            
      *MPop[i].symb = *MPop[i+1].symb;
      
      MPop[i].xyz->x = MPop[i+1].xyz->x;
      MPop[i].xyz->y = MPop[1+1].xyz->y;
      MPop[i].xyz->z = MPop[i+1].xyz->z;
      
      for ( k = i + 1 ; k > i ; k-- ){
	
	*MPop[k].symb = *MPop[k-1].symb;
	
	MPop[k].xyz->x = MPop[k-1].xyz->x;
	MPop[k].xyz->y = MPop[k-1].xyz->y;
	MPop[k].xyz->z = MPop[k-1].xyz->z;
      } 
      
      *MPop[k+1].symb = temp;
       
      MPop[k+1].xyz->x = XYZ_temp[0].x; 
      MPop[k+1].xyz->y = XYZ_temp[0].y;
      MPop[k+1].xyz->z = XYZ_temp[0].z;
      
    }
  }
}
  
