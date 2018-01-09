#define    IA 16807
#define    IM 2147483647
#define    AM (1.0/IM)
#define    IQ 127773
#define    IR 2836
#define    MASK 123459876

float ran0(long *idum)

/* “Minimal” random number generator of Park and Miller. Returns a uniform random deviate
between 0.0 and 1.0. Set or reset idum to any integer value (except the unlikely value 
MASK) to initialize the sequence; idum must not be altered between calls for successive 
deviates in a sequence. */

{
  long k;
  float ans;
  *idum ^= MASK;                 // XORing with MASK allows use of zero and other
  k=(*idum)/IQ;                  // simple bit patterns for idum.
  *idum=IA*(*idum-k*IQ)-IR*k;    // Compute idum=(IA*idum) % IM without over-
  if (*idum < 0) *idum += IM;    // flows by Schrage’s method.
  ans=AM*(*idum);                // Convert idum to a floating result.
  
  *idum ^= MASK;                 //  Unmask before return.
  
  return ans;
  
}
