#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "defmol.h"
#include "genindv.h"
#include "placeatomonzaxis.h"

#define ZERO          0.0E+0

extern double pi;

int GenIndv(char const *Infile, const char *potential, const double distconvf, const int natoms, Txyz *XYZ, defmol *indv)
{
  int    ii, j, n, i, act, ang, die, nseed=1;
  double distmin, distmax, angmax, angmin, diemin, diemax, modvec, modvec2;
  double X, Y, Z, xc, yc, zc, ran, angdie; 
  double status, cf = pi/180E+0, stetad, stetav, steta, teta;
  double yz, xz, sdie, a, b, c, delta, Y1, Y2, YA, delta1;
  FILE   *fl0;
  char   rdaux[4], symb[2];
  
  srand ( time(NULL) );

  if( !strncasecmp(potential,"Ab-Initio", 3) ){

    if ( ! (fl0=fopen(Infile,"r")) )
      {
	printf("|\n\nError: The file \"%s\" is not available!\n");
	exit(EXIT_FAILURE);
      }
    
    fscanf(fl0, "%s", &rdaux);

    while ( strncasecmp(rdaux,"*ZMAT", 5) ) fscanf(fl0, "%s \n", &rdaux);
      
    fscanf(fl0,"%s \n", &symb);      

    indv[0].atmn = 1;
    strcpy(indv[0].symb, symb);
    XYZ[0].x = ZERO;
    XYZ[0].y = ZERO;
    XYZ[0].z = ZERO;
    indv[0].xyz = &XYZ[0];
    indv[0].np = tableZ(symb);
    indv[0].mass = tableM(symb);

    fscanf(fl0,"%s %d %lf %lf\n", &symb, &act, &distmin, &distmax);
    printf("# 4.1 symb=%s act=%d distmin=%lf distmax=%lf  nseed=%d\n",symb, act, distmin, distmax, nseed);

    X = ZERO;
    Y = ZERO;
    ran = rand() % 101;
    ran = ran/100;
    Z = (distmax - distmin) * ran + distmin;  
 
    strcpy(indv[1].symb, symb);
    XYZ[1].x = X * distconvf;
    XYZ[1].y = Y * distconvf;
    XYZ[1].z = Z * distconvf;
    indv[1].xyz = &XYZ[1];
    indv[1].np = tableZ(symb);
    indv[1].mass = tableM(symb);

    for ( j = 2; j < natoms; j++ ) 
      {

	if ( j == 2 ){
	  fscanf(fl0,"%s %d %lf %lf %d %lf %lf \n", 
		 &symb, &act, &distmin, &distmax, &ang, &angmin, &angmax);
	  
	  printf("# 8-2 j=%d symb=%s, act=%d, distmin=%lf, distmax=%lf, ang=%d, angmin=%lf, angmax=%lf \n", j, symb, act, distmin, distmax, ang, angmin, angmax);
	}
	
	else {
	  fscanf(fl0,"%s %d %lf %lf %d %lf %lf %d %lf %lf \n", 
		 &symb, &act, &distmin, &distmax, &ang, &angmin, &angmax, 
		 &die, &diemin, &diemax); 
	}
	
	/* degrees to radians*/
	angmin = angmin * cf;
	angmax = angmax * cf;
	diemin = diemin * cf;
	diemax = diemax * cf;
	
	printf("# 8.1 j=%d symb=%s, act=%d, distmin=%lf, distmax=%lf, ang=%d, angmin=%lf, angmax=%lf die=%d, diemin=%lf, diemax=%lf \n", j, symb, act, distmin, distmax, ang, angmin, angmax, die, diemin, diemax);

	xc = XYZ[act - 1].x; 
	yc = XYZ[act - 1].y;
	zc = XYZ[act - 1].z;

	centerxyz(j-1, xc, yc, zc, XYZ);
	for(i = 0; i < j; i++ )
	  printf("##GENINDV1 i=%d, x=%lf, y=%lf, z=%lf \n", i, XYZ[i].x, XYZ[i].y, XYZ[i].z);

	status = PlaceAtmonZAxis(ang - 1, j-1, XYZ);
	printf("# status=%lf \n", status);
	
	status = PlaceAtmonXPlan(die - 1, j-1, XYZ);

	for(i = 0; i < j; i++ )
	  printf("##GENINDV2 i=%d, x=%lf, y=%lf, z=%lf \n", i, XYZ[i].x, XYZ[i].y, XYZ[i].z);
	
	if( XYZ[ang - 1].z < ZERO ) {
	  angmax = pi - angmax;
	  angmin = pi - angmin;
	  printf("distmax=%lf, distmin=%lf, cos(angmin)=%lf\n",distmax,distmin, cos(angmin));
	}
       
	/* compute the vector  module */
	ran =  rand() % 101;
	ran = ran /100;
	modvec = (distmax - distmin) * ran + distmin;
	
	/* Compute Z coordinate */
	ran = rand() % 1001;
	ran = ran/1000;
	Z = modvec * cos( (angmax - angmin) * ran + angmin );
	
	printf("# ZZZ j=%d Z=%lf \n", j, Z);
	
	/* Compute Y coordinate */
	if ( j != 2 )  {
	  ran = rand() % 1001;
	  ran=ran/1000;
	  angdie = (diemax - diemin) * ran + diemin;
	  
	  YA = modvec * sin((diemax - diemin) * ran + diemin);/**/ 
	 	  
	  sdie = sin(angdie);
	  
	  yz = - XYZ[die - 1].y * XYZ[ang - 1].z;
	  
	  xz =   XYZ[die - 1].x * XYZ[ang - 1].z;

	  YA = modvec * sdie;/**/ 

	  printf("# YYY yz=%lf xz=%lf sdie=%lf YA=%lf\n", yz, xz, sdie, YA);	  
	  
	  a =  yz*yz + xz*xz;
	  
	  b = -2.0E+0 * modvec * xz * sqrt(a) * sdie;

	  modvec2 = modvec*modvec;

	  c = yz*yz*Z*Z - yz*yz*modvec2 + modvec2 * a * sdie*sdie;

	  delta1 = -4*a*yz*yz*(Z*Z + modvec2*sdie*sdie - modvec2);

	  delta = b * b - 4.0E+0 * a * c;
	  
	  printf("# YYY a=%lf b=%lf c=%lf YA =%lf delta=%lf delta1=%lf\n", a, b, c, YA, delta, delta1);
	  printf("# YYY Z*Z =%lf modvec2*sdie*sdie=%lf modvec2=%lf \n",Z*Z, modvec2*sdie*sdie, modvec2);
	  if (delta < -1.0E-5){
	    printf ("Error [GenIndv]:Delta smaller than ZERO. Delta = %lf! \n", delta); 
	      exit(1);}
	  
	  if (delta < ZERO) delta = ZERO;
	  
	  Y1 = ( -b + sqrt(delta) ) / ( 2 * a ) ;
	  
	  Y =  ( -b - sqrt(delta) ) / ( 2 * a ) ;
	  
	  printf("# YYY angdie=%lf modvec=%lf yz=%lf xz=%lf sdie=%lf a=%lf b =%lf c=%lf delta=%17.15lf\n ", angdie, modvec, yz, xz, sdie, a, b, c, delta);
	  	  
	}
	
	else  Y = ZERO;
	
	printf("# YYY j=%d Y1=%lf Y2=%lf \n", j, Y1, Y2);

	/* Compute X coordinate */
	X = sqrt( modvec * modvec - Z * Z - Y * Y );
	
	if ( j > 2 ) {
	  stetav = sqrt(
			XYZ[ang-1].z * XYZ[die-1].y * XYZ[ang-1].z * XYZ[die-1].y +
			XYZ[ang-1].z * XYZ[die-1].x * XYZ[ang-1].z * XYZ[die-1].x  
			);
	  stetad = XYZ[ang-1].z * XYZ[die-1].x;
	  if ( sqrt(stetad * stetad) > stetav ) {
	    if      (stetad > ZERO ) stetad = stetad - 1.0E-15;
	    else if (stetad < ZERO ) stetad = stetad + 1.0E-15;
	  }


	  steta = stetad/stetav;
	 
	  /* 
	  if ( (sqrt(steta*steta) - 1.0E+0) > 1.0E-10 ) {
	    exit(1);
	  }
	  else if ( steta >  1.0E+0 ) steta =  1.0E+0;
	  else if ( steta < -1.0E+0 ) steta = -1.0E+0;
	  */
	  teta = asin ( steta );
	  
	  printf("teta =%lf,  stetad=%lf,  stetav=%lf, steta=%20.19lf, pi/2.0E+0=%lf\n",teta, stetad, stetav, steta, pi/2.0E+0);
	  printf("\nangdie =%lf,  (pi - teta)=%lf,  (2*pi - teta)=%lf",angdie,(pi - sqrt(teta*teta)),(2*pi - sqrt(teta*teta)) );
	  if ( teta > ZERO ) {
	    if ( angdie >=  (pi - teta) && angdie <=  (2*pi - teta) )  X = -X;
	    if ( angdie <= -(pi - teta) && angdie >= -(2*pi - teta) )  X = -X; 
	  }
	  else  if ( teta < ZERO ) {
	    if ( angdie <=  (pi + teta) || angdie >=  (2*pi + teta) )  X = -X;
	    if ( angdie >= -(pi + teta) || angdie <= -(2*pi + teta) )  X = -X; 
	  }
	}

	printf("# 10 j=%d modvec=%lf X=%lf Y=%lf Z=%lf \n", j, modvec, X, Y, Z);

	strcpy(indv[j].symb, symb);
	XYZ[j].x = X;
	XYZ[j].y = Y;
	XYZ[j].z = Z;
	indv[j].xyz = &XYZ[j];
	indv[j].np = tableZ(symb);
	indv[j].mass = tableM(symb);
	
	printf("#  XYZ.x=%lf, XYZ.y=%lf, XYZ.z=%lf, \n", XYZ[j].x, XYZ[j].y, XYZ[j].z);

      }
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
      
      strcpy(indv[j].symb, symb);
      XYZ[j].x = X;
      XYZ[j].y = Y;
      XYZ[j].z = Z;
      indv[j].xyz = &XYZ[j];
      indv[j].np = tableZ(symb);
      indv[j].mass = tableM(symb);
      
      if( !strncasecmp(potential,"tip", 3) ) {
	X = ran3_(nseed);
	Y = ran3_(nseed);
	Z = ran3_(nseed);
	
	strcpy(indv[j].symb, symb);
	XYZ[j].x = X;
        XYZ[j].y = Y;
	XYZ[j].z = Z;
	indv[j].xyz = &XYZ[j];
	indv[j].np = tableZ(symb);
	indv[j].mass = tableM(symb);
      }
    }
  }
}
