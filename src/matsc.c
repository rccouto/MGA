#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defmol.h"
#include "matsc.h"

#define ZERO          0.0E+0

int MatSC(int natoms, int cut, Tzmat *zmat)
{
  int i, m2, m3;
  
  m2 = natoms;
  m3 = 2*natoms;

  if ( cut == 0 ){
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
      zmat[m3 + i].dist = zmat[m2 + i].dist;
      zmat[m3 + i].angle = zmat[m2 + i].angle;
      zmat[m3 + i].dihedral = zmat[m2 + i].dihedral;
    }
  }

  else if ( cut == 1 ){
    zmat[m3].dist = ZERO;
    zmat[m3].angle = ZERO;
    zmat[m3].dihedral = ZERO;
    
    zmat[m3 + 1].dist = zmat[1].dist;
    zmat[m3 + 1].angle = ZERO;
    zmat[m3 + 1].dihedral = ZERO;
    
    zmat[m3 + 2].dist = zmat[2].dist;
    zmat[m3 + 2].angle = zmat[2].angle;
    zmat[m3 + 2].dihedral = ZERO;
    
    for ( i = 3; i < natoms; i++ ){
      zmat[m3 + i].dist = zmat[m2 + i].dist;
      zmat[m3 + i].angle = zmat[m2 + i].angle;
      zmat[m3 + i].dihedral = zmat[m2 + i].dihedral;
    }
  }
  
  else if ( cut > 1 ){
    zmat[m3].dist = ZERO;
    zmat[m3].angle = ZERO;
    zmat[m3].dihedral = ZERO;
    
    zmat[m3 + 1].dist = zmat[1].dist;
    zmat[m3 + 1].angle = ZERO;
    zmat[m3 + 1].dihedral = ZERO;
    
    zmat[m3 + 2].dist = zmat[2].dist;
    zmat[m3 + 2].angle = zmat[2].angle;
    zmat[m3 + 2].dihedral = ZERO;
    
    for ( i = 3; i < cut; i++ ){
      zmat[m3 + i].dist = zmat[i].dist;
      zmat[m3 + i].angle = zmat[i].angle;
      zmat[m3 + i].dihedral = zmat[i].dihedral;
    }
    
    for ( i = cut; i < natoms; i++ ){
      zmat[m3 + i].dist = zmat[m2 + i].dist;
      zmat[m3 + i].angle = zmat[m2 + i].angle;
      zmat[m3 + i].dihedral = zmat[m2 + i].dihedral;
    }
  }
}
