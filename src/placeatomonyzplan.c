#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defmol.h"
#include "rotate.h"
#include "placeatomonyzplan.h" 

extern double pi;

int PlaceAtmonYZPlan(const int N, const int natoms, Txyz *XYZ)
{
  int i;
  double A, B, XT, YT, ZT;
  double raiz;
  
   if( sqrt(XYZ[N].y * XYZ[N].y) <= 1.0E-10 &&
       sqrt(XYZ[N].z * XYZ[N].z) <= 1.0E-10 ) return(10);
   

   if ( sqrt(XYZ[N].x * XYZ[N].x) <= 1.0E-10 ) A = pi/2.0E+0;
   
   A = alpha_angle( XYZ[N].x , XYZ[N].y);
   
   for( i = 0; i < natoms; i++){
     rotate_z( XYZ[i].x, XYZ[i].y, A, &XT, &YT);
     XYZ[i].x = XT;
     XYZ[i].y = YT;
   }
 
   return(0);
}


