#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "genindiv.h"

#define ZERO          0.0E+0

extern distconvf;

int GenIndiv(char const *Infile, const char *potential, double const cf, const int natoms, defmol *indv)
{
  
  printf("\n  entrei 1 \n");

  int    ii, j, n, i, act, nseed;
  double distmin, distmax, ang, angmax, angmin, die, diemin, diemax, modvec, X, Y, Z, xc, yc, zc; 

  FILE   *fl0;
  char   rdaux[4], symb[2];
  Txyz   *pXYZ;
  
  printf("\n  entrei \n");
  scanf(yc);
  
  if( !strncasecmp(potential,"Ab-Initio", 9) ){

    if ( ! (fl0=fopen(Infile,"r")) )
      {
	printf("|\n\nError: The file \"%s\" is not available!\n");
	exit(EXIT_FAILURE);
      }

    fscanf(fl0, "%s", &rdaux);

    while ( strncasecmp(rdaux,"*ZMAT", 5) ) fscanf(fl0, "%s \n", &rdaux);
      
    fscanf(fl0,"%s \n", &symb[0]);      
    /*
    indv[0].atmn = 1;
    strcpy(indv[0].symb, symb);
    indv[0].xyz.x = ZERO;
    indv[0].xyz.y = ZERO;
    indv[0].xyz.z = ZERO;
    indv[0].np = tableZ(symb);
    indv[0].mass = tableM(symb);
    */
    X = ZERO;
    Y = ZERO;
    fscanf(fl0,"%s %d %lf %lf\n", &symb, act, distmin, distmax);
    Z = (distmax - distmin) * ran3_(nseed) + distmin;  
    /*
    strcpy(indv[1].symb, symb);
    indv[1].xyz.x = X * distconvf;
    indv[1].xyz.y = Y * distconvf;
    indv[1].xyz.z = Z * distconvf;
    indv[1].np = tableZ(symb);
    indv[1].mass = tableM(symb);
    */
    for ( j = 2; j = natoms; j = j + 3){

      if ( j == 2 ){
	fscanf(fl0,"%s %d %lf %lf %d %lf %lf \n", &symb, act, distmin, distmax, ang, angmin, angmax);
      }
      else{
	fscanf(fl0,"%s %d %lf %lf %d %lf %lf %d %lf %lf \n", &symb, act, distmin, distmax, ang, angmin, angmax, die, diemin, diemax); 
      }
      
      xc = indv[act - 1].xyz.x; 
      yc = indv[act - 1].xyz.y;
      zc = indv[act - 1].xyz.z;
      
      pXYZ = &indv[0].xyz;
	
      centerxyz(j-1, xc, yc, zc, pXYZ);

      atomonaxis(ang, act, pXYZ);

      /* compute the vector  module */
      modvec = (distmax - distmin) * ran3_(nseed) + distmin;

      /* Compute Z coordinate */
      Z = modvec * cos((angmax - angmin) * ran3_(nseed) + angmin);
	/* Compute Y coordinate */
	if ( j == 2 )  Y = ZERO;
	else           Y = modvec * sin((diemax - diemin) * ran3_(nseed) + diemin) ;	       
	/* Compute X coordinate */
	X = sqrt( modvec * modvec - Z * Z - Y * Y );
    }
    /*
    strcpy(indv[j].symb, symb);
    indv[j].xyz.x = X;
    indv[j].xyz.y = Y;
    indv[j].xyz.z = Z;
    indv[j].np = tableZ(symb);
    indv[j].mass = tableM(symb);
    */
    rewind(fl0);
    fclose(fl0);
  }

  else {
    
    if( !strncasecmp(potential,"tip", 3) ) n = 6;
    else n = 3;
    
    printf("#");
    for ( j=0; natoms; j=j+3){
      
      X = ran3_(nseed);
      Y = ran3_(nseed);
      Z = ran3_(nseed); 
      /*
      strcpy(indv[j].symb, symb);
      indv[j].xyz.x = X;
      indv[j].xyz.y = Y;
      indv[j].xyz.z = Z;
      indv[j].np = tableZ(symb);
      indv[j].mass = tableM(symb);
      */
      if( !strncasecmp(potential,"tip", 3) ) {
	X = ran3_(nseed);
	Y = ran3_(nseed);
	Z = ran3_(nseed);
	/*
	strcpy(indv[j].symb, symb);
	indv[j].xyz.x = X;
	indv[j].xyz.y = Y;
	indv[j].xyz.z = Z;
	indv[j].np = tableZ(symb);
	indv[j].mass = tableM(symb);
	*/
      }
	
    }
    
  }

}
