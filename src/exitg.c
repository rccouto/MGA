#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "exitg.h"

void exitg(char *name, char *keyword)
{
  int i, length;

  length = strlen(keyword);
  printf("|\nError: %s SubGroup \"",name);
  for (i=0; i<=length-1; i++) printf("%c",keyword[i]);
  printf("\" is not available!\n");
  exit (EXIT_FAILURE);
}
