#include <math.h>
#include "ran1.h"

float expdev(long *idum)
/* Returns an exponentially distributed, positive, random deviate of unit mean, using ran1(idum) as the source of uniform deviates. */
{
  float ran1(long *idum);
  float dum;
  
  do
    dum = ran1(idum);
  while (dum < 1.0E-1 );
  return -log10(dum);
}
