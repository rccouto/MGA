#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defmol.h"
#include "rtgroup.h"
#include "placeatomonzaxis.h"
#include "centerxyz.h"
#include "rotate.h"
#include "ran0.h"

extern double pi;

int RtGroup(int natoms, int nrotgp, long *nseed, Trtt *rtt, Tatom *rttatm, Txyz *XYZ)
{
  int    i, j;
  int    count, rot;
  double ran, A, XT, YT;
  double xc, yc, zc;
  
  count = 0;
  for( i = 0; i < nrotgp; i++ ){
    
    xc = XYZ[rtt[i].center - 1].x;
    yc = XYZ[rtt[i].center - 1].y;
    zc = XYZ[rtt[i].center - 1].z;

    centerxyz(natoms, xc, yc, zc, XYZ);
    
    PlaceAtmonZAxis(rtt[i].plane - 1, natoms - 1, XYZ);
 
    ran = ran0(nseed);
    A = (pi+pi)*ran + pi;

    for( j = 0; j < rtt[i].natmrot; j++ ){

      rot = (rttatm[count].atom - 1);

      rotate_z( XYZ[rot].x, XYZ[rot].y, A, &XT, &YT);
      XYZ[rot].x = XT;
      XYZ[rot].y = YT;
 
      count=count+1;
    }
  }
}
