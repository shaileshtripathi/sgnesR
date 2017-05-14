
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
//#include <Rinternals.h>
//#include <Rinterface.h>
#include<stdio.h>
using namespace std;
#include "sgns.h"
//////////////////////////////////////////
//int main (int argc, char *argv[]){     
//	int p = mainxx(argc-1, &argv[1]);
//}
///////////////////////////////////////////
extern "C" {

	SEXP getArgs(SEXP myint, SEXP mychar){
		int n;
		n = INTEGER_VALUE(myint);
		char *Pmychar[n];
		PROTECT(mychar = AS_CHARACTER(mychar));
		for(int i=0;i<n;i++){
			Pmychar[i] = R_alloc(strlen(CHAR(STRING_ELT(mychar, i))), sizeof(char));
		}	
		for(int i=0;i<n;i++){
			strcpy(Pmychar[i], CHAR(STRING_ELT(mychar, i)));
		}
		//for(int i=0;i<n;i++){
                	//printf(" %s \n",Pmychar[i]);
        	//}

		int p = mainxx(n, &Pmychar[0]);
	UNPROTECT(1);
	return(R_NilValue);

	}
R_CallMethodDef callMethods[]  = {
  {"getArgs", (DL_FUNC) &getArgs, 2},
  {NULL, NULL, 0}
};
void R_init_sgnesR(DllInfo *info)
{
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}



}
