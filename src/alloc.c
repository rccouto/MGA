#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "alloc.h"

char *Char_alloc (char *name2)
{
  char *name1;
  
  name1 = (char *) malloc ((strlen(name2) + 1) * sizeof(char)); /*
  name1 = (char *) calloc ((strlen(name2) + 1), sizeof(char)); /**/
  if (name1 == NULL) 
    {
      printf ("\n\n ** Error: Insufficient Memory **"); 
      exit(1);
      return (NULL);
    }
  strcpy(name1, name2);
  return (name1);
}

char *Char_free (char *name)
{
  if (name == NULL) return (NULL);
  free(name);        
  return (NULL); 
}



int *Int_alloca (int n)
{
  int *v;

  if (n < 1) 
    { 
      printf ("\n\n** [Int_alloca] **\n"); 
      printf ("** Error: Invalid parameter **\n");
      printf ("** n must be bigger than 0 **\n"); 
      printf ("** n = %d **\n",n);
      exit(1);
      return (NULL);
    }
  /* aloca o vetor /**/
  v = (int *) calloc (n, sizeof(int));  /**/
  /* v = (int *) malloc (n * sizeof(int)); /**/ 
  if (v == NULL) 
    {
      printf ("\n\n** Error: Insufficient Memory **");
      exit(1);
      return (NULL);
    }
  return (v);    
}

int *Int_freea (int *v)
{
  if (v == NULL) return (NULL);
  free(v);        
  return (NULL); 
}



double *Real_alloca (int n)
{
  double *v;        
  if (n < 1) 
    { 
      printf ("\n\n** Error: Invalid parameter **\n");
      printf ("** n must be bigger than 0 **\n"); 
      printf ("** n = %d **\n",n);
      exit(1);
      return (NULL);
    }
  /* aloca o vetor */
  v = (double *) calloc (n, sizeof(double));
  if (v == NULL) 
    {
      printf ("\n** Error: Insufficient Memory **");
      exit(1);
      return (NULL);
    }
  return (v);    
}

double *Real_freea (double *v)
{
  if (v == NULL) return (NULL);
  free(v);        
  return (NULL); 
}


int **Int_allocm (int m, int n)
{
  int  **v; 
  int    i;    
  if (m < 1 || n < 1) 
    { 
      printf ("\n\n** Error: Invalid parameter **\n");
      printf ("** m and n must be bigger than 0 **\n");
      printf ("** m = %d, n = %d **\n",m,n);
      exit(1);
      return (NULL);
    }
  /* aloca as linhas da matriz */
  v = (int **) calloc (m, sizeof(int *)); 
  if (v == NULL) 
    {
      printf ("**  Error: Insufficient Memory **");
      exit(1);
      return (NULL);
    }
  /* aloca as colunas da matriz */
  for ( i = 0; i < m; i++ ) 
    {
      v[i] = (int*) calloc (n, sizeof(int));      
      if (v[i] == NULL) 
	{
	  printf ("\n**  Error: Insufficient Memory  **");
	  exit(1);
	  return (NULL);
	}
    }
  return (v); /* retorna o ponteiro para a matriz */
}

int **Int_freem (int m, int n, int **v)
{
  int i; 
  if (v == NULL) return (NULL);
  if (m < 1 || n < 1) 
    { 
      printf ("\n\n** Error: Invalid parameter **\n");
      printf ("** m and n must be bigger than 0 **\n");
      printf ("** m = %d, n = %d **\n",m,n);
      exit(1);
      return (NULL);
    }
  for (i=0; i<m; i++) free (v[i]); /* libera as linhas da
				      matriz */
  free (v);       /* libera a matriz (vetor de ponteiros) */
  return (NULL);  /* retorna um ponteiro nulo */
}



double **Real_allocm (int m, int n)
{
  double **v; 
  int    i;    
  if (m < 1 || n < 1) 
    { 
      printf ("\n\n** Error: Invalid parameter **\n");
      printf ("** m and n must be bigger than 0 **\n");
      printf ("** m = %d, n = %d **\n",m,n);
      exit(1);
      return (NULL);
    }
  /* aloca as linhas da matriz */
  v = (double **) calloc (m, sizeof(double *)); 
  if (v == NULL) 
    {
      printf ("\n**  Error: Insufficient Memory **");
      exit(1);
      return (NULL);
    }
  /* aloca as colunas da matriz */
  for ( i = 0; i < m; i++ ) 
    {
      v[i] = (double*) calloc (n, sizeof(double));      
      if (v[i] == NULL) 
	{
	  printf ("\n**  Error: Insufficient Memory  **");
	  exit(1);
	  return (NULL);
	}
    }
  return (v); /* retorna o ponteiro para a matriz */
}

double **Real_freem (int m, int n, double **v)
{
  int i; 
  if (v == NULL) return (NULL);
  if (m < 1 || n < 1) 
    { 
      printf ("\n\n** Error: Invalid parameter **\n");
      printf ("** m and n must be bigger than 0 **\n");
      printf ("** m = %d, n = %d **\n",m,n);
      exit(1);
      return (NULL);
    }
  for (i=0; i<m; i++) free (v[i]); /* libera as linhas da
				      matriz */
  free (v);       /* libera a matriz (vetor de ponteiros) */
  return (NULL);  /* retorna um ponteiro nulo */
}
