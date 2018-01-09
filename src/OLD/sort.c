#include "sort.h"
#include "nrutil.h"
#include "defmol.h"

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

/* Here M is the size of subarrays sorted by straight insertion and NSTACK is the required auxiliary storage.*/

void sort(unsigned long n, defmol *indv)
  
/* Sorts an array arr[1..n] into ascending numerical order using the Quicksort algorithm. n is input; arr is replaced on output by its sorted rearrangement.*/
  
{
  unsigned long i, ir = n, j, k, l = 1, *istack;
  int jstack = 0;
  float a, temp;
  float arr[] = indv[].ergy
    
  istack = lvector(1,NSTACK);
  for (;;) {                              /* Insertion sort when subarray small enough.*/
    if ( ir - l < M ) {
      for ( j = l + 1 ; j <= ir ; j++ ) {
	a = arr[j];
	for ( i = j - 1 ; i >= l ; i-- ) {
	  if ( arr[i] <= a ) break;
	  arr[ i + 1 ] = arr[i];
	}
	arr[i+1]=a;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];               /* Pop stack and begin a new round of partitioning. */
      l = istack[jstack--];                   
    } else {
      k=( l + ir ) >> 1;                   /* Choose median of left, center, and right elements as */
      SWAP( arr[k],arr[l+1] )              /* partitioning element a.*/         
	if ( arr[l] > arr[ir] ) {          /* Also rearrange so that a[l] ≤ a[l+1] ≤ a[ir]. */
	  SWAP( arr[l],arr[ir] )
	    } 
      if ( arr[l+1] > arr[ir] ) {
	SWAP( arr[l+1],arr[ir] )
	  }
      if ( arr[l] > arr[l+1] ) {
	SWAP( arr[l],arr[l+1] )
	  }
      i = l + 1;                            /* Initialize pointers for partitioning. */
      j=ir;
      a=arr[l+1];                           /* Partitioning element. */
      for (;;) {                            /* Beginning of innermost loop. */
	do i++; while (arr[i] < a);                /* Scan up to find element > a. */
	do j--; while (arr[j] > a);                /* Scan down to find element < a. */
	if (j < i) break;                   /* Pointers crossed. Partitioning complete. */
	SWAP(arr[i],arr[j]);                /* Exchange elements. */
      }                                     /*  End of innermost loop. */
      arr[l+1]=arr[j];                      /*  Insert partitioning element. */
      arr[j]=a;
      jstack += 2;
      /* Push pointers to larger subarray on stack, process smaller subarray immediately. */
      if (jstack > NSTACK) nrerror("NSTACK too small in sort.");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free_lvector(istack,1,NSTACK);
  printf("%lf\n", arr[1]);
}

