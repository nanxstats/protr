#include <Rinternals.h>
#include <R.h>
#include <Rdefines.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>

/*##################################################*/
/*# Converts a String into a vector of characters  #*/
/*##################################################*/

SEXP s2c(SEXP seq){
  char *string;
  int lseq, i;
  char mot[2] = {'\0', '\0'};
    
  SEXP chaine;

  string = (char *) CHAR(STRING_ELT(seq, 0));
  
  lseq = strlen(string);
  
  PROTECT(chaine = NEW_CHARACTER(lseq));

  for(i = 0 ; i < lseq ; i++){  
    mot[0] = string[i];
    SET_STRING_ELT(chaine, i, mkChar(mot));
    }

  UNPROTECT(1);
  return(chaine);
}

