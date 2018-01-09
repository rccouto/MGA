#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defmol.h"
#include "rotate.h"
#include "placeatomonzaxis.h"

extern double pi;

int PlaceAtmonZAxis(const int N, const int natoms, Txyz *XYZ)
{
  int i;
  double A, XT, YT, ZT;
  double raiz;

  if( sqrt(XYZ[N].y * XYZ[N].y) <= 1.0E-10 &&  
      sqrt(XYZ[N].x * XYZ[N].x) <= 1.0E-10 ) return(10);
  
  if ( sqrt(XYZ[N].z * XYZ[N].z) <= 1.0E-10 ) A = pi/2.0E+0;
  else                                        A = alpha_angle( XYZ[N].y, XYZ[N].z );
  
  for( i = 0; i <= natoms; i++){
    rotate_x( XYZ[i].y, XYZ[i].z, A, &YT, &ZT);
    XYZ[i].y = YT;
    XYZ[i].z = ZT;
  }
  
  if ( sqrt(XYZ[N].z * XYZ[N].z) <= 1.0E-10 ) A = pi/2.0E+0;
  else                                        A = alpha_angle( XYZ[N].x, - XYZ[N].z );
  
  for( i = 0; i <= natoms; i++ ){
    rotate_y( XYZ[i].x, XYZ[i].z, A, &XT, &ZT );
    XYZ[i].x = XT;
    XYZ[i].z = ZT;
  }

  return(0);
}
