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
  double xl, yl, zl, xa, ya, za, xd, yd, zd, angang, modvecang, sdie, cang, sang, cdie;
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
	
	xc = XYZ[act - 1].x; 
	yc = XYZ[act - 1].y;
	zc = XYZ[act - 1].z;

	centerxyz(j-1, xc, yc, zc, XYZ);
	for(i = 0; i < j; i++ )
	  printf("#CENTER OUT i=%d, x=%lf, y=%lf, z=%lf \n", i, XYZ[i].x, XYZ[i].y, XYZ[i].z);
	
	status = PlaceAtmonZAxis(ang - 1, j-1, XYZ);
	printf("# status=%lf \n", status);
	
	status = PlaceAtmonXZPlan(die - 1, j-1, XYZ);
	
	/* Compute The Vector Module */
	ran =  rand() % 101;
	ran = ran /100;
	modvec = (distmax - distmin) * ran + distmin;
	
	/* Compute The Theta Angle */
	ran =  rand() % 101;
	ran = ran /100;
	if ( XYZ[ang - 1].z < ZERO ) angang = pi - ((angmax - angmin) * ran + angmin);
	else                         angang =        (angmax - angmin) * ran + angmin;
  
	/* Compute The Phi Angle */
	ran =  rand() % 101;
	ran = ran/100;
	if ( XYZ[die - 1].x < ZERO ) angdie = pi + ((diemax - diemin) * ran + diemin);
	else                         angdie =       -(diemax - diemin) * ran + diemin;

	//angdie = -(diemax - diemin) * ran + diemin; 
	
	cang = cos(angang);
	sang = sin(angang);
	cdie = cos(angdie);
	sdie = sin(angdie);
	
	X = modvec * sang * cdie;
	Y = modvec * sang * sdie;
	Z = modvec * cang;

	
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
