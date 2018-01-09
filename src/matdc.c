#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defmol.h"
#include "matdc.h"

#define ZERO          0.0E+0

int MatDC(int natoms, int  cut1, int cut2, Tzmat *zmat)
{
  int i, j, k, m2, m3;

  m2 = natoms;
  m3 = 2*natoms;

  if ( cut1 == 0 ){
    zmat[m3].dist = ZERO;
    zmat[m3].angle = ZERO;
    zmat[m3].dihedral = ZERO;
    
    zmat[m3 + 1].dist = zmat[1].dist;
    zmat[m3 + 1].angle = ZERO;
    zmat[m3 + 1].dihedral = ZERO;
    
    zmat[m3 + 2].dist = zmat[m2 + 2].dist;
    zmat[m3 + 2].angle = zmat[m2 + 2].angle;
    zmat[m3 + 2].dihedral = ZERO;
    
    for ( i = 3; i < natoms; i++ ){
      zmat[m3 + i].dist = zmat[i].dist;
      zmat[m3 + i].angle = zmat[i].angle;
      zmat[m3 + i].dihedral = zmat[i].dihedral;
    }
  }
  
  else if ( cut1 == 1 ){
    zmat[m3].dist = ZERO;
    zmat[m3].angle = ZERO;
    zmat[m3].dihedral = ZERO;
    
    zmat[m3 + 1].dist = zmat[1].dist;
    zmat[m3 + 1].angle = ZERO;
    zmat[m3 + 1].dihedral = ZERO;
    
    zmat[m3 + 2].dist = zmat[2].dist;
    zmat[m3 + 2].angle = zmat[2].angle;
    zmat[m3 + 2].dihedral = ZERO;
    
    for ( i = 3; i < cut2; i++ ){
      zmat[m3 + i].dist = zmat[m2 + i].dist;
      zmat[m3 + i].angle = zmat[m2 + i].angle;
      zmat[m3 + i].dihedral = zmat[m2 + i].dihedral;
    }
    
    for ( j = cut2; j < natoms; j++ ){
      zmat[m3 + j].dist = zmat[j].dist;
      zmat[m3 + j].angle = zmat[j].angle;
      zmat[m3 + j].dihedral = zmat[j].dihedral;
    }
  }
  
  else if ( cut1 > 1 ){
    zmat[m3].dist = ZERO;
    zmat[m3].angle = ZERO;
    zmat[m3].dihedral = ZERO;
    
    zmat[m3 + 1].dist = zmat[1].dist;
    zmat[m3 + 1].angle = ZERO;
    zmat[m3 + 1].dihedral = ZERO;
    
    zmat[m3 + 2].dist = zmat[2].dist;
    zmat[m3 + 2].angle = zmat[2].angle;
    zmat[m3 + 2].dihedral = ZERO;
    
    for ( i = 3; i <= cut1; i++ ){
      zmat[m3 + i].dist = zmat[i].dist;
      zmat[m3 + i].angle = zmat[i].angle;
      zmat[m3 + i].dihedral = zmat[i].dihedral;
    }
    
    for ( j = cut1 + 1; j < cut2; j++ ){
      zmat[m3 + j].dist = zmat[m2 + j].dist;
      zmat[m3 + j].angle = zmat[m2 + j].angle;
      zmat[m3 + j].dihedral = zmat[m2 + j].dihedral;
    }
    
    for ( k = cut2; k < natoms; k++ ){
      zmat[m3 + k].dist = zmat[k].dist;
      zmat[m3 + k].angle = zmat[k].angle;
      zmat[m3 + k].dihedral = zmat[k].dihedral;
    }
  }
}
