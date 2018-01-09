#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rotate.h"
/* 
   Ref: http://en.wikipedia.org/wiki/Rotation_matrix  
*/

double alpha_angle(double SX, double YS)
{
/*
  Compute angle to place an atom on a cartesian axis usig the rotation matrix.
  For rotate_x with   +Y (ou Z) place atom on XZ axis   -Y (ou Z) place atom on XY axis.
  For rotate_y with   +Z (ou X) place atom on XY axis   -Z (ou X) place atom on YZ axis. 
  For rotate_z with   +X (ou Y) place atom on YZ axis   -X (ou Y) place atom on XZ axis. 
*/  
  double alpha;
  alpha = atan(SX/YS);
  return (alpha);
}


void rotate_x(double Y, double Z, double A, double *YT, double *ZT)
{
  *YT = Y*cos(A) - Z*sin(A);
  *ZT = Y*sin(A) + Z*cos(A);
}
void rotate_y(double X, double Z, double A, double *XT, double *ZT)
{
  *XT =  X*cos(A) + Z*sin(A);
  *ZT = -X*sin(A) + Z*cos(A);
}
void rotate_z(double X, double Y, double A, double *XT, double *YT)
{
  *XT = X*cos(A) - Y*sin(A); 
  *YT = X*sin(A) + Y*cos(A);
}
