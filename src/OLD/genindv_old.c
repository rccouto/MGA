#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "defmol.h"
#include "genindv.h"
#include "placeatomonzaxis.h"
#include "placeatomonplan.h"

#define ZERO          0.0E+0

extern double pi;

int GenIndv(char const *Infile, const char *potential, const double distconvf, const int natoms, Txyz *XYZ, defmol *indv)
{
  int    ii, j, n, i, act, ang, die, nseed=1;
  double distmin, distmax, angmax, angmin, diemin, diemax, modvec, modvec2;
  double X, Y, Z, xc, yc, zc, ran, angdie, ang1; 
  double status, status2, cf = pi/180E+0, stetad, stetav, steta, teta;
  double xl, yl, zl, xa, ya, za, xd, yd, zd, angang, modvecang, sdie, cang;
  double xap, yap, zap, xdp, ydp, zdp, nx, ny, nz, ax, az, bx, by, c1, c2, c22, c12;
  double ax2, bx2, by2, az2, xp2, xp, a, b, c, delta, b2, xp1, yp, zp, modvecn;
  
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

    /* First Atom*/
    indv[0].atmn = 1;
    strcpy(indv[0].symb, symb);
    XYZ[0].x = ZERO;
    XYZ[0].y = ZERO;
    XYZ[0].z = ZERO;
    indv[0].xyz = &XYZ[0];
    indv[0].np = tableZ(symb);
    indv[0].mass = tableM(symb);
    
    fscanf(fl0,"%s %d %lf %lf\n", &symb, &act, &distmin, &distmax);
    printf("# FIRST symb=%s act=%d distmin=%lf distmax=%lf  nseed=%d\n",symb, act, distmin, distmax, nseed);

    /* Second Atom*/
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

    printf ("# SECOND X=%lf  Y=%lf  Z=%lf \n", X, Y, Z); 

    
    fscanf(fl0,"%s %d %lf %lf %d %lf %lf \n", &symb, &act, &distmin, &distmax, &ang, &angmin, &angmax);
    
    printf ("# 0 symb=%s, act=%d, distmin=%lf, distmax=%lf, ang=%d, angmin=%lf, angmax=%lf \n", symb, act, distmin, distmax, ang, angmin, angmax);
    
    /* Compute The Vector Module */
    ran =  rand() % 101;
    ran = ran /100;
    modvec = (distmax - distmin) * ran + distmin;

    /*Third Atom*/
    Y = ZERO + XYZ[act - 1].y;
    ran = rand() % 101;
    ran = ran/100;
    zap = XYZ[ang - 1].z - XYZ[act - 1].z;

    if ( zap < ZERO ) ang1 = pi - ((angmax - angmin) * ran + angmin) * cf ;
    else              ang1 =      ((angmax - angmin) * ran + angmin) * cf ;
    zp = modvec * cos( ang1 ); 
    Z = zp + XYZ[act - 1].z; 
    X = sqrt( modvec * modvec - zp * zp ) + XYZ[act - 1].x;
    
    strcpy(indv[2].symb, symb);
    XYZ[2].x = X * distconvf;
    XYZ[2].y = Y * distconvf;
    XYZ[2].z = Z * distconvf;
    indv[2].xyz = &XYZ[2];
    indv[2].np = tableZ(symb);
    indv[2].mass = tableM(symb);
    
    printf ("# THIRD X=%lf  Y=%lf  Z=%lf  \n", XYZ[2].x , XYZ[2].y, XYZ[2].z );
    
    for ( j = 3; j < natoms; j++ ) 
      {
	fscanf(fl0,"%s %d %lf %lf %d %lf %lf %d %lf %lf \n", 
	       &symb, &act, &distmin, &distmax, &ang, &angmin, &angmax, &die, &diemin, &diemax); 

	
	/* degrees to radians*/
	angmin = angmin * cf;
	angmax = angmax * cf;
	diemin = diemin * cf;
	diemax = diemax * cf;
	
	printf("# 1 j=%d symb=%s, act=%d, distmin=%lf, distmax=%lf, ang=%d, angmin=%lf, angmax=%lf die=%d, diemin=%lf, diemax=%lf \n", j, symb, act, distmin, distmax, ang, angmin, angmax, die, diemin, diemax);

	/* Compute The Vector Module */
	ran =  rand() % 101;
	ran = ran /100;
	modvec = (distmax - distmin) * ran + distmin;

	
	/* Identify Variables */
	xl = XYZ[act-1].x;
	yl = XYZ[act-1].y;
	zl = XYZ[act-1].z;

	printf ("# IDENTIFY 1 modvec=%lf xl=%lf  yl=%lf  zl=%lf \n", modvec, xl, yl, zl);	
    
	xa = XYZ[ang-1].x;
	ya = XYZ[ang-1].y;
	za = XYZ[ang-1].z;

	printf ("# IDENTIFY 2 xa=%lf  ya=%lf  za=%lf \n", xa, ya, za);
	
	xd = XYZ[die-1].x;
	yd = XYZ[die-1].y;
	zd = XYZ[die-1].z;

	printf ("# IDENTIFY 3 xd=%lf  yd=%lf  zd=%lf \n", xd, yd, zd);

	/* Compute Variables */
	
	xap = xa - xl;
	
	yap = ya - yl;
	
	zap = za - zl;
	
	printf ("# COMPUTE 1 xap=%lf  yap=%lf zap=%lf \n", xap, yap, zap);		
	/*
	xdp = xd - xap;
	
	ydp = yd - yap;
	
	zdp = zd - zap;
	
	printf ("# ******** COMPUTE 2.1 xdp=%lf  ydp=%lf  zdp=%lf \n", xdp, ydp, zdp);
	*/
	xdp = xd - xl;
	
	ydp = yd - yl;
	
	zdp = zd - zl;
	printf ("# ******** COMPUTE 2.2 xdp=%lf  ydp=%lf  zdp=%lf \n", xdp, ydp, zdp);
	
	ran =  rand() % 101;
	ran = ran /100;
	
	
	if ( zap < ZERO ) angang =  pi - ((angmax - angmin) * ran + angmin );
	else              angang =       ((angmax - angmin) * ran + angmin );
  

	ran =  rand() % 101;
	ran = ran/100;
	angdie = -(diemax - diemin) * ran + diemin; 
	
	modvecang = sqrt ( xap * xap + yap * yap + zap * zap );

	printf ("# COMPUTE 3 angang=%lf  angdie=%lf modvecang=%lf \n", angang, angdie, modvecang);
	
	nx = yap * zdp - zap * ydp;
	
	ny = zap * xdp - xap * zdp;
	
	nz = xap * ydp - yap * xdp;
	
	modvecn = sqrt(nx * nx + ny * ny + nz * nz);
	
	printf ("# COMPUTE 4 nx=%lf ny=%lf nz=%lf modvecn=%lf \n", nx, ny, nz, modvecn);
	
	ax = ny * xap - yap * nx; 
	
	az = ny * zap - yap * nz;
	
	bx = nz * xap - zap * nx;
	
	by = nz * yap - zap * ny;	

	cang = cos(angang);
	sdie = sin(angdie);

	printf ("# COMPUTE 5 ax=%lf az=%lf bx=%lf by=%lf sdie=%lf, cang=%lf \n", ax, az, bx, by, sdie, cang);
	
	c1 = ny * modvec * modvecang * cang - yap * modvecn * modvec * sdie;
	
	c2 = nz * modvec * modvecang * cang - zap * modvecn * modvec * sdie;
	
	printf ("# COMPUTE 6 angang=%lf angdie=%lf c1=%lf c2=%lf \n", angang, angdie, c1, c2);	
	
	c22 = c2 * c2;
	
	c12 = c1 * c1;
	
	ax2 = ax * ax;
	
	bx2 = bx * bx;
	
	by2 = by * by;
	
	az2 = az * az;
	
	modvec2 = modvec * modvec;

	printf ("# COMPUTE 7 c22=%lf c12=%lf ax2=%lf bx2=%lf by2=%lf az2=%lf xp2=%lf \n", c22, c12, ax2, bx2, by2, az2, xp2);	
	
	a = az2 * by2 + az2 * bx2 + ax2 * by2;
	
	b = c1 * ax * by2 + c2 * bx * az2;
	
	c = az2 * c22 + by2 * c12 - az2 * by2 * modvec2;
       
	if ( sqrt(a*a) < 1.0E-10 ) {
	  xp = c/(2 * b);
	  printf ("#### COMPUTE 8 a=%lf b=%lf c=%lf xp=%lf \n", a, b, c, xp);
	}
	else {

	  delta = 4.0E+0 * (b * b - a * c);
	  printf ("# COMPUTE 9 b*b=%lf a*c=%lf \n", b*b, a*c);
	  printf ("# COMPUTE 10 a=%lf b=%lf c=%lf delta=%lf \n", a, b, c, delta);
	
	  /* Compute X coordinate */
	
	  xp1 = ( 2.0E+0 * b - sqrt ( delta ) ) / (2*a);
	
	  xp2 = ( 2.0E+0 * b + sqrt ( delta ) ) / (2*a);	
	  
	  if ( zap < ZERO ) xp = xp2;
	  else              xp = xp1;
	  //	  xp = -0.513360;
	}
	
	X = xp + xl;

	/*
	stetav = sqrt(
		      XYZ[ang - 1].z * XYZ[die - 1].y * XYZ[ang - 1].z * XYZ[die - 1].y +
		      XYZ[ang - 1].z * XYZ[die - 1].x * XYZ[ang - 1].z * XYZ[die - 1].x  
		      );
	stetad = XYZ[ang - 1].z * XYZ[die - 1].x;

	if ( sqrt(stetad * stetad) > stetav ) {
	  if      (stetad > ZERO ) stetad = stetad - 1.0E-15;
	  else if (stetad < ZERO ) stetad = stetad + 1.0E-15;
	}
	steta = stetad/stetav;
	teta = asin ( steta );
	if ( teta > ZERO ) {
	  if ( angdie >=  (pi - teta) && angdie <=  (2*pi - teta) )  X = -X;
	  if ( angdie <= -(pi - teta) && angdie >= -(2*pi - teta) )  X = -X; 
	}
	else  if ( teta < ZERO ) {
	  if ( angdie <=  (pi + teta) || angdie >=  (2*pi + teta) )  X = -X;
	  if ( angdie >= -(pi + teta) || angdie <= -(2*pi + teta) )  X = -X; 
	} 
	*/
	printf("# XXX j=%d xp1=%lf xp2=%lf xp=%lf X=%lf \n", j, xp1, xp2, xp, X);
	
	/* Compute Y coordinate */
	
	yp = ( c2 - bx * xp ) / by; 
	
	Y = yp + yl;
	
	printf("# YYY j=%d yp=%lf Y=%lf \n", j, yp, Y);
	
	
	/* Compute Z coordinate */
	
	zp = ( c1 - ax * xp )/ az;
	
	Z = zp + zl; 
	
	printf("# ZZZ j=%d zp=%lf Z=%lf \n", j, zp, Z);
	
	/*if ( j > 2 )
	  {
	    stetav = sqrt(
			  XYZ[ang-1].z * XYZ[die-1].y * XYZ[ang-1].z * XYZ[die-1].y +
			  XYZ[ang-1].z * XYZ[die-1].x * XYZ[ang-1].z * XYZ[die-1].x  
			  );
	    stetad = XYZ[ang-1].z * XYZ[die-1].x;

	    /* if ( sqrt(stetad * stetad) > stetav ) {
	      if      (stetad > ZERO ) stetad = stetad - 1.0E-15;
	      else if (stetad < ZERO ) stetad = stetad + 1.0E-15;
	      }*/
	    
	    
	    /* steta = stetad/stetav;
	    
	    
	    /* 
	       if ( (sqrt(steta*steta) - 1.0E+0) > 1.0E-10 ) {
	       exit(1);
	       }
	       else if ( steta >  1.0E+0 ) steta =  1.0E+0;
	       else if ( steta < -1.0E+0 ) steta = -1.0E+0;
	    */
	    /* teta = asin ( steta );
	    
	    printf("teta =%lf,  stetad=%lf,  stetav=%lf, steta=%20.19lf, pi/2.0E+0=%lf\n",teta, stetad, stetav, steta, pi/2.0E+0);
	    printf("\nangdie =%lf,  (pi - teta)=%lf,  (2*pi - teta)=%lf",angdie,(pi - sqrt(teta*teta)),(2*pi - sqrt(teta*teta)) );
	    /* if ( teta > ZERO ) 
	      {
		if ( angdie >=  (pi - teta) && angdie <=  (2*pi - teta) )  X = -X;
		if ( angdie <= -(pi - teta) && angdie >= -(2*pi - teta) )  X = -X; 
		}*/
	    /*  else  if ( teta < ZERO )
	      {
		if ( angdie <=  (pi + teta) || angdie >=  (2*pi + teta) )  X = -X;
		if ( angdie >= -(pi + teta) || angdie <= -(2*pi + teta) )  X = -X; 
		}*//*
	    
	  }*/
	
	printf("# FINAL j=%d X=%lf Y=%lf Z=%lf \n", j, X, Y, Z);
	
	strcpy(indv[j].symb, symb);
	XYZ[j].x = X;
	XYZ[j].y = Y;
	XYZ[j].z = Z;
	indv[j].xyz = &XYZ[j];
	indv[j].np = tableZ(symb);
	indv[j].mass = tableM(symb);

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
