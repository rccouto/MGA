#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

extern double pi;

int atantwo(double y, double x)
{
  double a, reslt;
 
  a = y/x;
  
  if ( x > 0 ) {

    reslt = atan(a);

  }
  
  else if ( x < 0 ){
    
    if ( y >= 0 ) reslt = atan(a) + pi;
    
    if ( y < 0 )  reslt = atan(a) - pi;    
  }
  
  else if ( x == 0 ) {

    if ( y == 0 ) reslt = 0.0;
      
    else if ( y > 0 ) reslt = pi/2;
      
    else if ( y < 0 ) reslt = - pi/2;
    
  }
  
  return(reslt);
}
