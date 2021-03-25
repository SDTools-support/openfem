/*
 * MAT-file diagnose program
 *
 * See the MATLAB API Guide for compiling information.
 *
 * Calling syntax:
 *
 *   matdgns <matfile>
 *
 * It will diagnose the MAT-file named <matfile>.
 *
 * This program demonstrates the use of the following functions:
 *
 *  matClose
 *  matGetDir
 *  matGetNextVariable
 *  matGetNextVariableInfo
 *  matOpen
 *
 * Copyright 1984-2003 The MathWorks, Inc.
 */
/* $Revision: 1.1 $ */
#include <stdio.h>
#include <stdlib.h>
#include "mat.h"
#include <vector>

extern "C" void mexFunction (int nlhs, mxArray *plhs[],
                          int nrhs, mxArray *prhs[]); 

int main(int argc, char **argv)
{
	if (argc <= 1)return 0; 
 
	const char *file=argv[1];
	MATFile *pmat;
	const char **dir;
	const char *name;
	int	  ndir;
	int	  i;
	mxArray *pa;
	std::vector<mxArray *> vec_pa;
  
	printf("Reading file %s...\n\n", file);

	/* Open file to get directory */
	pmat = matOpen(file, "r");
	if (pmat == NULL)
	{	printf("Error opening file %s\n", file);return(1); }
  
	/* get directory of MAT-file */
	dir = (const char **)matGetDir(pmat, &ndir);
	if (dir == NULL)
	{
		printf("Error reading directory of file %s\n", file);
		return(1);
	}
	else
	{
		printf("Directory of %s:\n", file);
		for (i=0; i < ndir; i++)  printf("%s\n",dir[i]);
	}
	mxFree(dir);

	/* In order to use matGetNextXXX correctly, reopen file to read in headers. */
	if (matClose(pmat) != 0)  {printf("Error closing file %s\n",file); return(1);}
	pmat = matOpen(file, "r");
	if (pmat == NULL)
	{ printf("Error reopening file %s\n", file); return(1); }

	/* Read in each array. */

	printf("\nReading in the actual array contents:\n");
	for (i=0; i<ndir; i++)
	{
		pa = matGetNextVariable(pmat, &name);
		if (pa == NULL)
		{
			printf("Error reading in file %s\n", file); return(1);
		} 
		/*
		* Diagnose array pa
		*/
		printf("According to its contents, array %s has %d dimensions\n",name, mxGetNumberOfDimensions(pa));
		if (mxIsFromGlobalWS(pa))printf("  and was a global variable when saved\n");
		else printf("  and was a local variable when saved\n");
		//mxDestroyArray(pa);
	}
	for (i=0;i<mxGetNumberOfElements(pa);i++) vec_pa.push_back(mxGetCell(pa,i));


	if (matClose(pmat) != 0) {printf("Error closing file %s\n",file); return(1);}

	//---------------------------------------------------------------------------------------------

	int nlhs=10;
	mxArray *plhs[10];

	mexFunction (nlhs,plhs,vec_pa.size(),&(vec_pa[0])); 

	printf("Creating file %s...\n\n", file);
	pmat = matOpen("..\\..\\K_of_mk.mat", "w");
	if (pmat == NULL)
	{
		printf("Error creating file %s\n", file);
		printf("(Do you have write permission in this directory?)\n");
		return(EXIT_FAILURE);
	}

	//matPutVariableAsGlobal
	int status = matPutVariable(pmat,"K",vec_pa[12]);
	if (status != 0)
	{
		printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
		return(EXIT_FAILURE);
	}  
  
	if (matClose(pmat) != 0)
	{
		printf("Error closing file %s\n",file);
		return(EXIT_FAILURE);
	}

	//---------------------------------------------------------------------------------------------

	for (i=0; i<ndir; i++) {mxDestroyArray(vec_pa[i]); }
	/*for (i=0; i<10; i++) { if (plhs[i]!=NULL) mxDestroyArray(plhs[i]); }*/
	
	printf("Done\n");
	
	return(0);
}