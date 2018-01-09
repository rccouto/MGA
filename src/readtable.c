#include <stdio.h>
#include <stdlib.h>

#define LNCHAR    61

int tableZ (char *atom) {
  
  char symb[3], rdaux[LNCHAR]; 
  int np, scanout;
  double mass;

  FILE *fp;
  fp = fopen("data/atominfo.dat",  "r");
  if (! fp) {
    printf("\n\nError: Atom data file not found!\n");
    exit(EXIT_FAILURE);
  }
  /* Pre-check the file atominfo.dat */  
  fgets(rdaux, LNCHAR, fp );
  if( strncmp(rdaux, "**Symbol**      **Atomic Number (Z)**      **Mass (u.m.a.)**", LNCHAR) != 0) {
    printf("\n\nError:Atom data file seens to be crashed!\n");
    fclose(fp); 
    exit(EXIT_FAILURE);
    return 1; 
  }

  do{ 
    fscanf(fp, "%s %d %lf \n", symb, &np, &mass);
    //printf("scanout = %d %s\n", scanout, atom);
    if (strcmp(atom, symb) == 0) {
      fclose(fp);
      // printf(" is %d\n", np);
      return np;
    }
  } while (strcmp(symb, "EOF") != 0);
    
  printf("\n\nError: %s's atomic number is not included in the data table!\n", atom); 
  fclose(fp); 
  exit(EXIT_FAILURE);
  return 1;
}

double tableM (char *atom) {

  char symb[3], rdaux[LNCHAR];
  int np;
  double mass;
  
  FILE *fp; 
  fp = fopen("data/atominfo.dat",  "r");
  if (! fp) {
    printf("Atom data file not found!\n");
    exit(EXIT_FAILURE);
  }

  /* Pre-check the file atominfo.dat */  
  fgets(rdaux, LNCHAR, fp );
  if( strncmp(rdaux, "**Symbol**      **Atomic Number (Z)**      **Mass (u.m.a.)**", LNCHAR) != 0) {
    printf("\n\nError:Atom data file seens to be crashed!\n");
    fclose(fp); 
    exit(EXIT_FAILURE);
    return 1; 
  }

  //printf("[DEBUG] %s's mass", atom);
  
  do{
    fscanf(fp, "%s %d %lf \n", symb, &np, &mass);
    if (strcmp(atom, symb) == 0) {
      fclose(fp);
      // printf(" is %lf\n", mass);
      mass = mass * (1.672621637/9.10938215) * 1.0e4;           // Mass conversion to Atomic Units
      return mass;
    }
  } while (strcmp(symb, "EOF") != 0);

  printf("\n\nError: %s's mass is not included in the data table!\n", atom);  
  fclose(fp); 
  exit(EXIT_FAILURE);
  return 1;
}
