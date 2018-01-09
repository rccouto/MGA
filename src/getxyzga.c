#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "defmol.h"
#include "getxyz.h"

int GetXYZGA(const char *fdat, const char *frst, int nat, int nstep, double *energy, int *conv, defmol *indv)
{
  char    lixo[20], lixo1[20];
  char    ident[65], str[42];
  char    first[29] = ("GEOMETRY");
  char    second[29] = ("*-*-*INCOMPLETE*-*-*");
  char    third[29] = ("$GRAD");
  char    nserch[100];
  int     j, pconv;
  double  energ;
  
  FILE   *arqdat, *arqrst;



  if ( ! (arqdat=fopen(fdat,"r")) ) {
    printf("|\n\nError [GetXYZGA]: The file \"%s\" is not available for use!\n",fdat);
    exit(EXIT_FAILURE);
  }
  
  pconv = 1;

  printf("\nGETXYZ\n");

  /* Used when geometry couldn't be optimized */

  if ( arqrst=fopen(frst,"r") ) {
    printf("#1\n");
    fgets(nserch, 65, arqrst);
    
    do { 
      fgets(ident, 65, arqdat);
    } while( !feof(arqdat) && strcmp(ident,nserch) );
    printf("#2\n");
    if ( strncasecmp(ident,nserch,65) == 0 ){
      printf("#3 nat=%d\n", nat); 
      fseek ( arqdat , 161  , SEEK_CUR );
      for(j = 0; j < nat; j++){
	printf("#4\n");
	fscanf(arqdat,"%s %s %lf %lf %lf",  &lixo1, &lixo, &indv[j].xyz->x, &indv[j].xyz->y, &indv[j].xyz->z);
      }
      fseek ( arqdat , 158, SEEK_CUR );
      fscanf( arqdat, "%lf", &energ );
      *energy = energ;
      pconv = 1;
      //printf("[GetXYZGA]: (NO) conv = %d\n", pconv);
    }
    fclose(arqrst);
  }
  
  else{
    /* Used when wasn't enought NSTEP's to optimize the geometry */
    printf("\n[Second Case]\n");
    do { 
      fscanf( arqdat, "%s", &ident );
    } while( !feof(arqdat) && strncasecmp(ident,"*-*-*INCOMPLETE*-*-*", 20) );
    
    if( strcmp(ident,second) == 0 ) {
      rewind(arqdat);
      
      if ( nstep < 100 ) sprintf(str,"-------------------- DATA FROM NSERCH=  %d", nstep);
      else if ( nstep < 1000 ) sprintf(str,"-------------------- DATA FROM NSERCH= %d", nstep);
      else if ( nstep >= 1000 ) sprintf(str,"-------------------- DATA FROM NSERCH=%d", nstep);
      
      do { 
	fgets(ident, 43, arqdat);
      } while( !feof(arqdat) && strcmp(ident,str) );
      
      if ( strncasecmp(ident,str,43) == 0 ){
	fseek ( arqdat , 183  , SEEK_CUR );
	for(j = 0; j < nat; j++){
	  fscanf(arqdat,"%s %s %lf %lf %lf",  &lixo1, &lixo, &indv[j].xyz->x, &indv[j].xyz->y, &indv[j].xyz->z);
	}
	fseek ( arqdat , 158, SEEK_CUR );
	fscanf( arqdat, "%lf", &energ );
	*energy = energ;
	pconv = 1;
	//printf("[GetXYZGA]: (ES) conv = %d\n", pconv);
      }
    }
    
    /* Used when geometry could be optimized */
    else{
      printf("\n[Third Case]\n");
      rewind(arqdat);
      do { 
	fscanf( arqdat, "%s", &ident );
      } while( !feof(arqdat) && strncasecmp(ident,"GEOMETRY", 12) );

      if( strcmp(ident,first) == 0 ){

	fseek ( arqdat , 231 , SEEK_CUR );

	for(j = 0; j < nat; j++){
	  fscanf(arqdat,"%s %s %lf %lf %lf",  &indv[j].symb, &lixo, &indv[j].xyz->x, &indv[j].xyz->y, &indv[j].xyz->z);
	}
	
	do { 
	  fscanf( arqdat, "%s", &ident );
	  //printf("%s\n", ident);
	} while( !feof(arqdat) && strncasecmp(ident,"$GRAD", 12) );	
	
	if( strcmp(ident,third) == 0 ){
	  fseek ( arqdat , 9 , SEEK_CUR );
	  fscanf( arqdat, "%lf", &energ );
	}
      	//fseek ( arqdat , 74, SEEK_CUR );
	//fscanf( arqdat, "%lf", &energ );
	*energy = energ;
	pconv = 0;
	//printf("[GetXYZGA]: (FO) conv = %d\n", pconv);
      }
    }
  }

  *conv=pconv;
  if( pconv == 1 ) printf("[GetXYZGA]: Energy(%s) = %lf *\n", fdat, energ);
  else printf("[GetXYZGA]: Energy(%s) = %lf\n", fdat, energ);
  fclose(arqdat);
  return 0;
}
