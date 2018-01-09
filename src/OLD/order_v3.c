#include <stdlib.h>
#include <stdio.h>

#include "defmol.h"
#include "order.h"
#include "alloc.h"

void Order(int natoms, int nindiv, double *energy, defmol *MPop)
{
  int    i, j, k, l, ni, nj, nk, nks, nkm, ip;
  double E_temp;
  Txyz   *XYZ_temp;
  defmol *pMPop;

  XYZ_temp = (Txyz *) calloc( natoms, sizeof(Txyz) );

  for(i = 0; i < nindiv; i++) printf("1######## Energy[%d] = %.10lf\n", i, energy[i]);
  printf("\n");
  
  for ( i = 1; i < nindiv; i++ ){
    for ( j = 0; j < i; j++ ){
      //printf("energy[%d] = %lf, energy[%d] = %lf\n", j, energy[j], i, energy[i]);
      if ( energy[j] > energy[i] ){
	//printf("*IF energy[%d] = %lf, energy[%d] = %lf\n", j, energy[j], i, energy[i]);	
	for(i = 0; i < nindiv; i++) printf("2######## Energy[%d] = %.10lf\n", i, energy[i]);

	printf("\n");

	E_temp = energy[j];

	for(i = 0; i < nindiv; i++) printf("2.5######## Energy[%d] = %.10lf\n", i, energy[i]);

	printf("\n");
	
	//printf("*IF E_temp = %lf\n", E_temp);

	/*nj = j * natoms;	
	for ( l = 0 ; l < natoms ; l++ ){
	  XYZ_temp[l].x = MPop[nj + l].xyz->x;
	  XYZ_temp[l].y = MPop[nj + l].xyz->y;
	  XYZ_temp[l].z = MPop[nj + l].xyz->z;
	  }*/
	
	energy[j] = energy[i];

	for(i = 0; i < nindiv; i++) printf("3######## Energy[%d] = %.10lf\n", i, energy[i]);

	printf("\n");
	
	//ni = i * natoms;

	/*for ( l = 0 ; l < natoms ; l++ ){
	  MPop[nj + l].xyz->x = MPop[ni + l].xyz->x;
	  MPop[nj + l].xyz->y = MPop[ni + l].xyz->y;
	  MPop[nj + l].xyz->z = MPop[ni + l].xyz->z;
	  }*/
	
	for ( k = i ; k > j ; k-- ){

	  //printf("* FOR K: i = %d, j = %d, k = %d\n", i, j, k);

	  energy[k] = energy[k - 1];

	  for(i = 0; i < nindiv; i++) printf("4######## Energy[%d] = %.10lf\n", i, energy[i]);

	  printf("\n");
	  

	  //printf("* FOR K1: energy[k - 1] = %lf energy[1] = %lf\n", energy[k - 1], energy [1]);

	  /*nk = k * natoms;
	  nkm = (k - 1) * natoms;
	  for ( l = 0 ; l < natoms ; l++ ){

	    MPop[nk + l].xyz->x = MPop[nkm + l].xyz->x;
	    MPop[nk + l].xyz->y = MPop[nkm + l].xyz->y;
	    MPop[nk + l].xyz->z = MPop[nkm + l].xyz->z;
	    }*/

	  //printf("* FOR K2: E_temp = %lf\n", E_temp);

	  for(i = 0; i < nindiv; i++) printf("5######## Energy[%d] = %.10lf\n", i, energy[i]);

	  printf("\n");

	  energy[k + 1] = E_temp;

	  for(i = 0; i < nindiv; i++) printf("6######## Energy[%d] = %.10lf\n", i, energy[i]);

	  printf("\n");

	  /*nks = (k + 1) * natoms;
	  
	  printf("*  FOR K3: energy[k + 1] = %lf energy[1] = %lf \n", energy[k + 1], energy [1]); 

	  for ( l = 0 ; l < natoms ; l++ ){
	    MPop[nks + l].xyz->x = XYZ_temp[l].x; 
	    MPop[nks + l].xyz->y = XYZ_temp[l].y;
	    MPop[nks + l].xyz->z = XYZ_temp[l].z;
	    }*/
	  //printf("*END IF energy[%d] = %lf, energy[%d] = %lf\n", j, energy[j], i, energy[i]);  
	}
	for(i = 0; i < nindiv; i++) printf("7######## Energy[%d] = %.10lf\n", i, energy[i]);

	printf("\n");
      }
    }
  }
}
 
    /*
      printf("energy = %10.10lf\n", energy[0]);
      printf("energy OLD = %10.10lf\n", energy_old[0]);
      
      for( i = 0; i < nindiv; i++ ){
      for( j = 0; j < nindiv; j++ ){
      if( energy[j] == energy_old[i]){
      for( k = 0; k < natoms; k++ ){
      MPop[ ( natoms * j ) + k ].xyz->x = MPop_old[ ( natoms * i ) + k ].xyz->x;
      MPop[ ( natoms * j ) + k ].xyz->y = MPop_old[ ( natoms * i ) + k ].xyz->y;
      MPop[ ( natoms * j ) + k ].xyz->z = MPop_old[ ( natoms * i ) + k ].xyz->z;
      }
      }
      }
      }
    */
 
