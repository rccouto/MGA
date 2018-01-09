#include <stdio.h>
#include <stddef.h> 
#include <stdlib.h> 
#include <unistd.h> 
#include <string.h>
#include <sys/types.h> 
#include <sys/wait.h> 

#include "runprog.h"
 
/* Execute the command using this shell program.  */
#define SHELL "/bin/bash"
     
int RunProg (const int k, int typgs, int gmsver, int ncore, const char *flag, const char *prog, 
	     const char *path, char *name, char *eff)
{
  char   ck[3], command[50], file[15], gv[3], nc[3];
  int    status;
  pid_t  pid;
  
  // Dalton
  if ( !strncasecmp(prog,"DALTON", 5) )
    strcpy(command,"dalton.sh");

  // GAMESS
  else if ( !strncasecmp(prog,"GAMESS", 5) ){
    if (typgs==1) {
      strcpy(command,"gamessiga.sh"); 
      sprintf(file,"%s", name);
    }
      else {
	strcpy(command,"gamess.sh"); 
	strcpy(file, "Gamess");
      }
    }
  
  // Open Babel
  else if ( !strncasecmp(prog,"OBabel", 6) ){
    strcpy(command,"obabel.sh");
    flag = eff;
  }
  else 
    exit (EXIT_FAILURE);

  
  if      ( k <= 9  )  sprintf(ck, "000%d", k);
  else if ( k <= 99 )  sprintf(ck, "00%d", k);
  else if ( k <= 999)  sprintf(ck, "0%d", k);
  else if ( k <= 9999) sprintf(ck, "%d", k);
  else               exit(EXIT_FAILURE);

  if ( gmsver < 9 ) sprintf(gv, "0%d", gmsver);
  else              sprintf(gv, "%d", gmsver);

  sprintf(nc, "%d", ncore);

  pid = fork ();
  if (pid == 0)
    {
      /* This is the child process.  Execute the shell command. */
      execl (SHELL, SHELL, command, path, ck, flag, file, gv, nc, NULL); 
      _exit (EXIT_FAILURE);
    }
  else if (pid < 0)
    /* The fork failed.  Report failure.  */
    status = -1;
  else
    /* This is the parent process.  Wait for the child to complete.  */
    if (waitpid (pid, &status, 0) != pid)
      status = -1;

  return status;
}

