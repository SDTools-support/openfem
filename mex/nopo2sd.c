/************************/
/* NOPO2SD.C            */
/* Auteur : Frank Genot */
/* Date : 13/09/01      */
/* Copyright INRIA      */
/************************/
#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef min
#define min(x,y) (((x) < (y))?(x):(y))
#endif

/**********************************/
/*                                */
/* Le type descripteur de tableau */
/*                                */
/**********************************/
typedef struct {
	char name[5];
	int address;
	int size_int;
	int size; 
	int type;
	char comment[73];
	int *donnees_int;
	float *donnees_float;
	char *donnees_char;
	double *donnees_double;
} descripteur_tableau;

/************************************/
/*                                  */
/* Les differents types de probleme */
/*                                  */
/************************************/
#define MAX_PROBLEM 7
char *problemname[MAX_PROBLEM] = {
  "2D",
  "3D", 
  "AXI", 
  "FOURIER",
  "INCOMPRESSIBLE",
  "PLAQUE", 
  "COQUE"
};

#define D2             0 
#define D3             1
#define AXI            2
#define FOURIER        3
#define INCOMPRESSIBLE 4
#define PLAQUE         5
#define COQUE          6

/******************************/
/*                            */
/* Les types d'elements finis */
/*                            */
/******************************/
#define UNKNOWN       -1

/* 2D */
#define TRIA_2P1D      0
#define TRIA_2P2D      1
#define QUAD_2Q1D      2
#define QUAD_2Q2R      3

/* 3D */
#define TETR_3P1D      4
#define TETR_3P2D      5
#define PENT_3R1D      6
#define PENT_3R2D      7
#define HEXA_3Q1D      8
#define HEXA_3Q2R      9

/* AXI */
#define TRIA_AP1D      10
#define TRIA_AP2D      11
#define QUAD_AQ1D      12
#define QUAD_AQ2R      13

/* FOURIER */
#define QUAD_AQ2F      14

/* INCOMPRESSIBLE */
#define QUAD_5NOE      15

/* PLAQUE */
#define TRIA_DKTP      16

/* COQUE */
#define QUAD_MITC4     17


/*******************************/
/*                             */
/* Les noms des elements (sdt) */
/*                             */
/*******************************/

#define MAX_FEM 18
char *femname[MAX_FEM] = {
  "t3p",          /* 2D */
  "t6p",
  "q4p",
  "q8p",
  "tetra4",       /* 3D */
  "tetra10",
  "penta6",
  "penta15",
  "hexa8",
  "hexa20",
  "t3a",          /* AXI */
  "t6a",
  "q4a",
  "q8a",
  "q9a",          /* FOURIER */
  "q5p",          /* INCOMPRESSIBLE */
  "tria3",        /* PLAQUE */
  "mitc4"         /* COQUE */
};

typedef struct{
  int ncge; /* Code de la geometrie */
  int nmae; /* Nombre de mots necessaires au stockage des numeros de references des faces, aretes et sommets */ 
  int ndsde; /* Numero de sous-domaine */
  int nno; /* Nombre de noeuds */
  int *nono; /* numero de chaque noeud */
  int npo; /* si ncopnp=0, nombre de points de l'element */
  int *nopo; /* numero du point */
  int ining; /* si nmae != 0, indicateur de niveau d'information */
  int *ref; /* les refs */
  int fem;
} ELEMENT;

/* La taille du champ 'otherinfo' dans la 'finite element model description matrice' (sdt) */
int femotherinfo[MAX_FEM] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

/* Les differents types d'elements du maillage (modulef) */
#define NOEUD      0
#define SEGMENT    1
#define TRIANGLE   2
#define QUADRANGLE 3
#define TETRAEDRE  4
#define PENTAEDRE  5
#define HEXAEDRE   6

#define max(a,b) ((a)>(b)?(a):(b))

typedef struct{
  int som[2]; /* les 2 sommets */
} ARETE;

typedef struct{
  int nsom;   /* nombre de sommets numerotes de 1 a nsom */
  int *som;   /* liste des sommets */ 
} FACE;

typedef struct{
  int nsom;
  int narete;
  int nface;
  ARETE *arete;
  FACE *face;
} ELEM;

typedef enum{small,big} ENDIAN;


/********************************************/
/*                                          */
/* Devine le type d'element fini a utiliser */
/*                                          */
/********************************************/
void guess_FEM(ELEMENT *elem, int problem, int n1, int iset, int iseq, int isete, int isepe, int isehe)
{
  elem->fem = UNKNOWN;
  switch(problem){
  case D2:
    switch(elem->ncge - 1){
    case NOEUD:
      break;
    case SEGMENT:
      break;
    case TRIANGLE:
      if(iset == 0)
        switch(n1){
        case 0: elem->fem = TRIA_2P1D; break;
        case 1: elem->fem = TRIA_2P2D; break;
        }
      break;
    case QUADRANGLE:
      if(iseq == 0)
        switch(n1){
        case 0: elem->fem = QUAD_2Q1D; break;
        case 1: elem->fem = QUAD_2Q2R; break;
        }
      break;
    }
    break;
  case D3:
    switch(elem->ncge - 1){
    case TETRAEDRE:
      if(isete == 0 && iset == 0)
        switch(n1){
        case 0: elem->fem = TETR_3P1D; break;
        case 1: elem->fem = TETR_3P2D; break;
        }
      break;
    case PENTAEDRE:
      if(isepe == 0 && iset == 0 && iseq == 0)
        switch(n1){
        case 0: elem->fem = PENT_3R1D; break;
        case 1: elem->fem = PENT_3R2D; break;
        }
      break;
    case HEXAEDRE:
      if(isehe == 0 && iseq == 0)
        switch(n1){
        case 0: elem->fem = HEXA_3Q1D; break;
        case 1: elem->fem = HEXA_3Q2R; break;
        }
      break;
    }
    break;
  case AXI:
    switch(elem->ncge - 1){
    case TRIANGLE:
      if(iset == 0)
        switch(n1){
        case 0: elem->fem = TRIA_AP1D; break;
        case 1: elem->fem = TRIA_AP2D; break;
        }
      break;
    case QUADRANGLE:
      if(iseq == 0)
        switch(n1){
        case 0: elem->fem = QUAD_AQ1D; break;
        case 1: elem->fem = QUAD_AQ2R; break;
        }
      break;
    }
    break;
  case FOURIER:    
    switch(elem->ncge - 1){
    case QUADRANGLE:
      switch(iseq){
      case 1:
        switch(n1){
        case 1: elem->fem = QUAD_AQ2F; break;
        }
        break;
      }
      break;
    }
    break;
  case INCOMPRESSIBLE:
    switch(elem->ncge - 1){
    case QUADRANGLE:
      if(n1 == 0 && iseq == 1)
        elem->fem = QUAD_5NOE;
      break;
    }
    break;
  case PLAQUE:
    switch(elem->ncge - 1){
    case TRIANGLE: 
      if(iset == 0)
        switch(n1){
        case 0: elem->fem = TRIA_DKTP; break;
        }
      break;
    }
    break;
  case COQUE:
    switch(elem->ncge - 1){
    case QUADRANGLE: 
      if(iseq == 0)
        switch(n1){
        case 0: elem->fem = QUAD_MITC4; break;
        }
      break;
    }
    break; 
  }
}


void add_end_of_string(char *s, int n)
{
  int i;
  
  for(i = n - 1; (i >= 0) && (s[i] == ' '); i--);
  s[i + 1] = '\0';
}

void small2big(void *i, int nbyte)
{
  unsigned char c[8], *ii;
  int j;

  ii = (unsigned char*)i;
  for(j = 0; j < nbyte; j++)
    c[nbyte - j - 1] = *(ii + j);

  memcpy(i, c, nbyte);
}

void interpol(int ndim, int nsomi, float *somi[], float *res[], int n)
{
  int i, j;
  float dir[3];
  
  if(n != 1 && nsomi != 2){
    fprintf(stderr, "Computing more than one internal point in face not implemented...\n");
    exit(5);
  }
  
  if(n == 1){
    for(i = 0; i < ndim; i++)
      res[0][i] = 0;
    for(i = 0; i < ndim; i++){
      for(j = 0; j < nsomi; j++)
        res[0][i] += somi[j][i];
      
      res[0][i] /= nsomi;
    }
  }
  else{ /* Cutting edge in n */
    for(i = 0; i < ndim; i++)
      dir[i] = (somi[1][i] - somi[0][i])/(n + 1);

    for(j = 1; j <=n; j++)
      for(i = 0; i < ndim; i++)
        res[j][i] = somi[0][i] + j*dir[i];
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char *fname, *pname, *oname, *verb, errormessage[1024];
  int lenstr[3], status;
  
  int VERBOSE;

  FILE *fnopo, *fnode, *fmd;

  int problem, elementfini;
  int dummy_int, i, j, k;
  double dummy_double;
  ENDIAN endian;

  /* Variable for nopofile */
  char titre[81];
  char date[9];
  char nomcre[25];
  char typesd[5];
  int niveau, etat, ntacm;
	descripteur_tableau *tab;

  int ndim, ndsr, ndsd, ncopnp, ne, nepo, nseg, ntri, nqua;
  int ntet, npen, nhex, nsup, nef, noe, n1, iset, iseq, isete;
  int isepe, isehe, np, ntycoo, lpgdn, nbegm, lnop5, ntacoo;

  float r, theta, z, phi;
  float *som;
  int *refsom;
  float *node;
  int *refnode;
  float *somi[4]; /* Tableau de pointeurs sur les sommets a interpoler */
  float **resi;   /* Tableau de pointeurs sur les noeuds interpoles */
  int offsetnode, offsetref;
  
  ELEMENT *elem;
  int nref;

  double *node_matrix; /* Node matrix */

  double *elm_matrix, *redge_matrix, *rface_matrix;   /* Model description matrix */

  int nbfem[MAX_FEM]; /* Nombre de FEMs de chaque type */

  int nbrow_FEelt, nbcol_FEelt, nbcol_rface, nbcol_redge; /* Nombre de lignes, colonnes de la 'finite element model description matrix' + matrices de references des edges et faces*/
  int nfacecur, nedgecur;
  int row, col;

  ELEM element[7];

  if (nrhs==1){ 
   plhs[0]=mxCreateString("$Revision: 1.10 $  $Date: 2005/07/04 08:29:18 $");
   return;
  }

   
  /* Check for proper number of input arguments */
  if(nrhs != 2 && nrhs != 3)
    mexErrMsgTxt("There should be at least 2 input arguments.");       

  /* First input argument must be a string */
  if(mxIsChar(prhs[0]) != 1)
    mexErrMsgTxt("First argument <nopofile basename> must be a string.");
   
  /* First input argument must be a row vector */
  if(mxGetM(prhs[0]) != 1)
    mexErrMsgTxt("First argument <nopofile basename> must be a row vector.");
   
  /* Get length of first input argument */
  lenstr[0] = mxGetN(prhs[0]) + 1;
   
  /* Allocate memory for nopofilename */
  fname = mxCalloc(lenstr[0] + strlen(".nopo"), sizeof(char));
   
  /* Copy the string data from prhs[0] into fname */
  status = mxGetString(prhs[0], fname, lenstr[0]);
  if(status != 0)
    mexErrMsgTxt("Not enough space.");
   
  /* Add extension .nopo to fname */
  strcat(fname, ".nopo");

  if((fnopo = fopen(fname, "rb")) == NULL){
    sprintf(errormessage, "%s: No such file.", fname);
    mexErrMsgTxt(errormessage);
  }

  /* Second input argument must be a string */
  if(mxIsChar(prhs[1]) != 1)
    mexErrMsgTxt("Second argument <problem> must be a string.");
      
  /* Second input argument must be a row vector */
  if(mxGetM(prhs[1]) != 1)
    mexErrMsgTxt("Second argument <problem> must be a row vector.");
      
  /* Get length of second input argument */
  lenstr[1] = mxGetN(prhs[1]) + 1;
   
  /* Allocate memory for problemname */
  pname = mxCalloc(lenstr[1], sizeof(char));
   
  /* Copy the string data from prhs[1] into pname */
  status = mxGetString(prhs[1], pname, lenstr[1]);
  if(status != 0)
    mexErrMsgTxt("Not enough space.");

  for(problem = 0; problem < MAX_PROBLEM; problem++)
    if(!strcmp(problemname[problem], pname))
      break;
  if(problem == MAX_PROBLEM){
    sprintf(errormessage, "Second argument <problem> should be in {");
    for(j = 0; j < MAX_PROBLEM - 1; j++)
      sprintf(errormessage + strlen(errormessage), "%s, ", problemname[j]);
    sprintf(errormessage + strlen(errormessage)," %s}.", problemname[MAX_PROBLEM-1]);
    mexErrMsgTxt(errormessage);
  }   
   
  VERBOSE = 0;
  if(nrhs == 3){
    /* Third input argument must be a string */
    if(mxIsChar(prhs[2]) != 1)
      mexErrMsgTxt("Optional third argument should be a string.");
         
    /* Third input argument must be a row vector */
    if(mxGetM(prhs[2]) != 1)
      mexErrMsgTxt("Optional third argument should be a row vector.");
         
    /* Get length of third input argument */
    lenstr[2] = mxGetN(prhs[2]) + 1;
      
    /* Allocate memory for verb */
    verb = mxCalloc(lenstr[2], sizeof(char));
      
    /* Copy the string data from prhs[2] into verb */
    status = mxGetString(prhs[2], verb, lenstr[2]);
    if(status != 0)
      mexErrMsgTxt("Not enough space.");
      
    if (verb[0]=='v') {VERBOSE = 1;}

  }

  else VERBOSE = 0;

  /* Check for proper number of output arguments */
  if(nlhs > 4)
    mexErrMsgTxt("Too many output arguments.");


if(VERBOSE){   
  mexPrintf("Reading nopo-file...\n");
}
  
  /*****************************************/
  /* TABLEAU NOP0 : informations generales */
  /*****************************************/
    
  /* Check for endian */
  fread(&dummy_int, sizeof(int), 1, fnopo);
  if(dummy_int > 10000){
    endian = small; 
    small2big(&dummy_int, sizeof(int));
  }
  else endian = big;
if(VERBOSE){
  mexPrintf("Offset de depart = %d\n", dummy_int);
}
  fseek(fnopo, dummy_int*sizeof(char) + 3*sizeof(int), SEEK_CUR);

  /* Read titre : 80 char */
  fread(titre, sizeof(char), 80, fnopo);
  add_end_of_string(titre, 80);
if(VERBOSE){
  mexPrintf("titre = %s\n", titre);
}

  /* Read date : 8 char */
  fread(date, sizeof(char), 8, fnopo);
  add_end_of_string(date, 8);
if(VERBOSE){
  mexPrintf("date = %s\n", date);
}

  /* Read nomcre (nom du createur) : 24 char */
  fread(nomcre, sizeof(char), 24, fnopo);
  add_end_of_string(nomcre, 24);
if(VERBOSE){
  mexPrintf("nomcre = %s\n", nomcre);
}
  
  /* Verify sd type (NOPO) */
  fread(typesd, sizeof(char), 4, fnopo);
  add_end_of_string(typesd, 4);
if(VERBOSE){
  mexPrintf("typesd = %s\n", typesd);
}
  if(strcmp(typesd, "NOPO") != 0){
    sprintf(errormessage, "%s is not a 'NOPO' file.\n", fname);
    fclose(fnopo);
    mexErrMsgTxt(errormessage);
  }

  /* Read niveau : 1 int */
  fread(&niveau, sizeof(int), 1, fnopo); if(endian == small) small2big(&niveau, sizeof(int));
if(VERBOSE){
  mexPrintf("niveau = %d\n", niveau);
}
  
  /* Read etat : 1 int */
  fread(&etat, sizeof(int), 1, fnopo); if(endian == small) small2big(&etat, sizeof(int));
if(VERBOSE){
  mexPrintf("etat = %d\n", etat);
}

  /* Read ntacm (nombre de tableaux supplementaires associes) : 1 int */
  fread(&ntacm, sizeof(int), 1, fnopo); if(endian == small) small2big(&ntacm, sizeof(int));
if(VERBOSE){
  mexPrintf("ntacm = %d\n", ntacm);
}

	/* On saute l'entier final (longueur de NOP0) */
	fseek(fnopo, 1*sizeof(int), SEEK_CUR);
	
  /**********************************************************************/
  /* TABLEAU NOP1 : Descripteur des eventuels tableaux supplementaires */
  /**********************************************************************/
  
  /* On ne peut exploiter les donnees de ce tableau car internes a Modulef */
  if(ntacm != 0){
    mexPrintf("Don't know what to do with table NOP1 : ignored.\n");
		
		/* On saute l'entier initial (longueur de NOP1) */
		/* On saute l'entier suivant (22*ntacm) */
		fseek(fnopo, 2*sizeof(int), SEEK_CUR);
		
		/* On alloue la memoire necessaire pour stocker les infos */
		tab = mxCalloc(ntacm, sizeof(descripteur_tableau));
		
		/* Lecture de ntacm tableaux */
    for(i = 0; i < ntacm; i++){
if(VERBOSE){
	mexPrintf("-->Lecture de la description du tableau %d\n", i+1);
}
			/* On lit le nom du tableau */
			fread(tab[i].name, sizeof(char), 4, fnopo);
			add_end_of_string(tab[i].name, 4);
if(VERBOSE){
	mexPrintf("---->Nom du tableau : %s\n",tab[i].name);
}
	
			/* On lit l'adresse dans le super tableau M du tableau */
			fread(&tab[i].address, sizeof(int), 1, fnopo);if(endian == small) small2big(&tab[i].address, sizeof(int));
if(VERBOSE){
	mexPrintf("---->Adresse dans le super tableau : %d\n",tab[i].address);
}
			/* On lit le nombre de mots de ce tableau */
			fread(&tab[i].size_int, sizeof(int), 1, fnopo);if(endian == small) small2big(&tab[i].size_int, sizeof(int));
if(VERBOSE){
	mexPrintf("---->Nombre de mots du tableau : %d\n",tab[i].size_int);
}
			
			/* On lit le type de ce tableau */
			fread(&tab[i].type, sizeof(int), 1, fnopo);if(endian == small) small2big(&tab[i].type, sizeof(int));
if(VERBOSE){
	mexPrintf("---->Type du tableau : %d\n",tab[i].type);
}

			/* On alloue la memoire necessaire */
			tab[i].donnees_int = NULL;
			tab[i].donnees_float = NULL;
			tab[i].donnees_char = NULL;
			tab[i].donnees_double = NULL;
			switch(tab[i].type){
				case 1 : /* entier */
					tab[i].size = tab[i].size_int*sizeof(int)/sizeof(int);
					tab[i].donnees_int = mxCalloc(tab[i].size, sizeof(int));
					break;
				case 2 : /* reel simple */
					tab[i].size = tab[i].size_int*sizeof(int)/sizeof(float);
					tab[i].donnees_float = mxCalloc(tab[i].size, sizeof(float));					
					break;
				case 4 : /* caractere */
					tab[i].size = tab[i].size_int*sizeof(int)/sizeof(char);
					tab[i].donnees_char = mxCalloc(tab[i].size+1, sizeof(char));
					break;
				case 5 : /* reel double */
					tab[i].size = tab[i].size_int*sizeof(int)/sizeof(double);
					tab[i].donnees_double = mxCalloc(tab[i].size, sizeof(double));
					break;
			}
			
			/* On lit le commentaire associe a ce tableau */
			fread(tab[i].comment, sizeof(char), 72, fnopo);
			add_end_of_string(tab[i].comment, 72);
if(VERBOSE){
	mexPrintf("---->Commentaire du tableau : %s\n",tab[i].comment);
}

		} 
		/* On saute le mot de fin */
		fseek(fnopo, 1*sizeof(int), SEEK_CUR);
				
		/* On lit maintenant le contenu des tableaux */
		for(i = 0; i < ntacm; i++){
			
			/* On saute le debut de l'enregistrement */
			fseek(fnopo, 1*sizeof(int), SEEK_CUR);
			
if(VERBOSE){
	mexPrintf("-->Les donnees du tableau %d\n", i+1);
} 
			/* On saute l'entier de la longueur, on l'a deja */
			fseek(fnopo, 1*sizeof(int), SEEK_CUR);
			
			for(j = 0; j < tab[i].size; j++)
				switch(tab[i].type){
					case 1: /* Entier */
						fread(&tab[i].donnees_int[j], sizeof(int), 1, fnopo); if(endian == small) small2big(&tab[i].donnees_int[j], sizeof(int));		
if(VERBOSE){
	mexPrintf("%d ", tab[i].donnees_int[j]);
} 
						break;
					case 2: /* Reel simple */
						fread(&tab[i].donnees_float[j], sizeof(float), 1, fnopo); if(endian == small) small2big(&tab[i].donnees_float[j], sizeof(float));								
if(VERBOSE){
	mexPrintf("%f ", tab[i].donnees_float[j]);
} 						
						break;
					case 4: /* Caractere */
						fread(&tab[i].donnees_char[j], sizeof(char), 1, fnopo);
if(VERBOSE){
	mexPrintf("%c", tab[i].donnees_char[j]);
}						
						break;
					case 5: /* Reel double */
						fread(&tab[i].donnees_double[j], sizeof(double), 1, fnopo); if(endian == small) small2big(&tab[i].donnees_double[j], sizeof(double));		
if(VERBOSE){
	mexPrintf("%lf ", tab[i].donnees_double[j]);
}						
						break; 
				}
			
			if(tab[i].type == 4)
				add_end_of_string(tab[i].donnees_char, tab[i].size);
			
if(VERBOSE){
	mexPrintf("\n");
}				
			/* On saute la fin de l'enregistrement */
			fseek(fnopo, 1*sizeof(int), SEEK_CUR);
		}
		
		/* On saute sur la premiere entree du tableau suivant */
 		fseek(fnopo, 2*sizeof(int), SEEK_CUR);
  }
	else
		/* On saute l'entier de debut et la longueur du tableau */
  	fseek(fnopo, 2*sizeof(int), SEEK_CUR);

  /***************************************************/
  /* TABLEAU NOP2 : Description generale du maillage */
  /***************************************************/ 
  
  /* Read ndim (dimension de l'espace (2,3)) */
  fread(&ndim, sizeof(int), 1, fnopo); if(endian == small) small2big(&ndim, sizeof(int));
if(VERBOSE){
  mexPrintf("ndim = %d\n", ndim);
}

  /* Read ndsr (maximum des numeros de reference) */
  fread(&ndsr, sizeof(int), 1, fnopo); if(endian == small) small2big(&ndsr, sizeof(int));
if(VERBOSE){
  mexPrintf("ndsr = %d\n", ndsr);
}

  /* Read ndsd (maximum des numeros de sous-domaine */
  fread(&ndsd, sizeof(int), 1, fnopo); if(endian == small) small2big(&ndsd, sizeof(int));
if(VERBOSE){
  mexPrintf("ndsd = %d\n", ndsd);
}
  
  /* Read ncopnp (code de coincidence des noeuds et des sommets) */
  fread(&ncopnp, sizeof(int), 1, fnopo); if(endian == small) small2big(&ncopnp, sizeof(int));
if(VERBOSE){
  mexPrintf("ncopnp = %d\n", ncopnp);
}

  /* Read ne (nombre d'elements du maillage) */
  fread(&ne, sizeof(int), 1, fnopo); if(endian == small) small2big(&ne, sizeof(int));
if(VERBOSE){
  mexPrintf("ne = %d\n", ne);
}

  /* Read nepo (nombre d'elements reduits a un point) */
  fread(&nepo, sizeof(int), 1, fnopo); if(endian == small) small2big(&nepo, sizeof(int));
if(VERBOSE){
  mexPrintf("nepo = %d\n", nepo);    
}
  
  /* Read nseg (nombre de segments) */
  fread(&nseg, sizeof(int), 1, fnopo); if(endian == small) small2big(&nseg, sizeof(int));
if(VERBOSE){
  mexPrintf("nseg = %d\n", nseg);
}
  
  /* Read ntri (nombre de triangles) */
  fread(&ntri, sizeof(int), 1, fnopo); if(endian == small) small2big(&ntri, sizeof(int));
if(VERBOSE){
  mexPrintf("ntri = %d\n", ntri);
}

  /* Read nqua (nombre de quadrangles) */
  fread(&nqua, sizeof(int), 1, fnopo); if(endian == small) small2big(&nqua, sizeof(int));
if(VERBOSE){
  mexPrintf("nqua = %d\n", nqua);
}
  
  /* Read ntet (nombre de tetraedres) */
  fread(&ntet, sizeof(int), 1, fnopo); if(endian == small) small2big(&ntet, sizeof(int));
if(VERBOSE){
  mexPrintf("ntet = %d\n", ntet);
}

  /* Read npen (nombre de pentaedres) */
  fread(&npen, sizeof(int), 1, fnopo); if(endian == small) small2big(&npen, sizeof(int));
if(VERBOSE){
  mexPrintf("npen = %d\n", npen);
}

  /* Read nhex (nombre d'hexaedres) */
  fread(&nhex, sizeof(int), 1, fnopo); if(endian == small) small2big(&nhex, sizeof(int));
if(VERBOSE){
  mexPrintf("nhex = %d\n", nhex);
}

  /* Read nsup (nombre de super-elements) */
  fread(&nsup, sizeof(int), 1, fnopo); if(endian == small) small2big(&nsup, sizeof(int));
if(VERBOSE){
  mexPrintf("nsup = %d\n", nsup);
}

  /* Read nef (nombre d'elements frontaliers) */
  fread(&nef, sizeof(int), 1, fnopo); if(endian == small) small2big(&nef, sizeof(int));
if(VERBOSE){
  mexPrintf("nef = %d\n", nef);
}

  /* Read noe (nombre de noeuds) */
  fread(&noe, sizeof(int), 1, fnopo); if(endian == small) small2big(&noe, sizeof(int));
if(VERBOSE){
  mexPrintf("noe = %d\n", noe);
}

  /* Read n1 (nombre de noeuds sur un segment ou une arete (hors extremites)) */
  fread(&n1, sizeof(int), 1, fnopo); if(endian == small) small2big(&n1, sizeof(int));
if(VERBOSE){
  mexPrintf("n1 = %d\n", n1);
}

  /* Read iset (le nombre de noeuds internes a un triangle ou une face triangulaire)) */
  fread(&iset, sizeof(int), 1, fnopo); if(endian == small) small2big(&iset, sizeof(int));
if(VERBOSE){
  mexPrintf("iset = %d\n", iset);
}

  /* Read iseq (nombre de noeuds internes a un quadrangle ou une face quadrangulaire) */
  fread(&iseq, sizeof(int), 1, fnopo); if(endian == small) small2big(&iseq, sizeof(int));
if(VERBOSE){
  mexPrintf("iseq = %d\n", iseq);
}

  /* Read isete (nombre de noeuds internes a un tetraedre) */
  fread(&isete, sizeof(int), 1, fnopo); if(endian == small) small2big(&isete, sizeof(int));
if(VERBOSE){
  mexPrintf("isete = %d\n", isete);
}

  /* Read isepe (nombre de noeuds internes a un pentaedre) */
  fread(&isepe, sizeof(int), 1, fnopo); if(endian == small) small2big(&isepe, sizeof(int));
if(VERBOSE){
  mexPrintf("isepe = %d\n", isepe);
}

  /* Read isehe (nombre de noeuds internes a hexaedre) */
  fread(&isehe, sizeof(int), 1, fnopo); if(endian == small) small2big(&isehe, sizeof(int));
if(VERBOSE){
  mexPrintf("isehe = %d\n", isehe);
}

  /* Read np (nombre de points) */
  fread(&np, sizeof(int), 1, fnopo); if(endian == small) small2big(&np, sizeof(int));
if(VERBOSE){
  mexPrintf("np = %d\n", np);
}
  
  /* Read ntycoo (type des valeurs des coordonnees) */
  fread(&ntycoo, sizeof(int), 1, fnopo); if(endian == small) small2big(&ntycoo, sizeof(int));
if(VERBOSE){
  mexPrintf("ntycoo = %d\n", ntycoo);
}

  /* Read lpgdn (la plus grande difference entre les numeros des noeuds d'une element + 1) */
  fread(&lpgdn, sizeof(int), 1, fnopo); if(endian == small) small2big(&lpgdn, sizeof(int));
if(VERBOSE){
  mexPrintf("lpgdn = %d\n", lpgdn);
}

  /* Read nbegm (nombre de super-elements ou de descriptions dans NOP3) */
  fread(&nbegm, sizeof(int), 1, fnopo); if(endian == small) small2big(&nbegm, sizeof(int));
if(VERBOSE){
  mexPrintf("nbegm = %d\n", nbegm);
}

  /* Read lnop5 (nombre de mots du tableau NOP5) */
  fread(&lnop5, sizeof(int), 1, fnopo); if(endian == small) small2big(&lnop5, sizeof(int));
if(VERBOSE){
  mexPrintf("lnop5 = %d\n", lnop5);
}

  /* Read ntacoo (type des axes de coordonnees) */
  fread(&ntacoo, sizeof(int), 1, fnopo); if(endian == small) small2big(&ntacoo, sizeof(int)); 
if(VERBOSE){
  mexPrintf("ntacoo = %d\n", ntacoo);
}

  /************************************/
  /* TABLEAU NOP3 : Pointeur eventuel */
  /************************************/

  if(nbegm != 0){
    mexWarnMsgTxt("Don't know what to do with table NOP3 : ignored.");
    fseek(fnopo, (6 + nbegm)*sizeof(int), SEEK_CUR);
  }
  else fseek(fnopo, 3*sizeof(int), SEEK_CUR);

  /******************************************/
  /* TABLEAU NOP4 : Coordonnees des sommets */
  /******************************************/

  /* Allocate memory for 'som' */
  som = (float *)mxCalloc(ndim*np, sizeof(float));
  refsom = (int *)mxCalloc(np, sizeof(int));

if(VERBOSE){
  mexPrintf("Coordonnees des sommets :\n");
}
  for(i = 0; i < np; i++){
    refsom[i] = 0; /* Default : no ref */
    fread(som + ndim*i, sizeof(float), ndim, fnopo);
    if(endian == small)
      for(j = 0; j < ndim; j++)
        small2big(som + ndim*i + j, sizeof(float));
if(VERBOSE){
    mexPrintf("sommet %d : ", i + 1);
    for(j = 0; j < ndim; j++)
      mexPrintf("%f ", som[i*ndim + j]);
    mexPrintf("\n");
}
    /* Bring back to Cartesian coordinates if needed */
    if(ntacoo != 1){
      r = som[ndim*i];
      theta = som[ndim*i + 1];
      if(ntacoo == 2){
        som[ndim*i] = r*cos(theta);
        som[ndim*i + 1] = r*sin(theta);
        /* z if exists not modified */
      }
      else /* ntacoo == 3 */
        if(ndim == 3){
          phi = som[ndim*i + 2];
          som[ndim*i] = r*cos(phi)*sin(theta);
          som[ndim*i + 1] = r*sin(phi)*sin(theta);
          som[ndim*i + 2] = r*cos(theta);
        }
    } 
  }


  /*----------------*/
  /* define element */
  /*----------------*/

    element[NOEUD].nsom = 1;
    element[NOEUD].narete = 0;
    element[NOEUD].nface = 0;
    
    element[SEGMENT].nsom = 2;
    element[SEGMENT].narete = 1;
    element[SEGMENT].nface = 0;
    element[SEGMENT].arete = (ARETE *)mxCalloc(1, sizeof(ARETE));
    element[SEGMENT].arete[0].som[0] = 1;
    element[SEGMENT].arete[0].som[1] = 2;
    
    element[TRIANGLE].nsom = 3;
    element[TRIANGLE].narete = 3;
    element[TRIANGLE].nface = 1;
    element[TRIANGLE].arete = (ARETE *)mxCalloc(3, sizeof(ARETE));
    element[TRIANGLE].arete[0].som[0] = 1;
    element[TRIANGLE].arete[0].som[1] = 2;
    element[TRIANGLE].arete[1].som[0] = 2;
    element[TRIANGLE].arete[1].som[1] = 3;
    element[TRIANGLE].arete[2].som[0] = 3;
    element[TRIANGLE].arete[2].som[1] = 1;
    element[TRIANGLE].face = (FACE *)mxCalloc(1, sizeof(FACE));
    element[TRIANGLE].face[0].nsom = 3;
    element[TRIANGLE].face[0].som = (int *)mxCalloc(3, sizeof(int));
    element[TRIANGLE].face[0].som[0] = 1;
    element[TRIANGLE].face[0].som[1] = 2;
    element[TRIANGLE].face[0].som[2] = 3;
    
    element[QUADRANGLE].nsom = 4;
    element[QUADRANGLE].narete = 4;
    element[QUADRANGLE].nface = 1;
    element[QUADRANGLE].arete = (ARETE *)mxCalloc(4, sizeof(ARETE));
    element[QUADRANGLE].arete[0].som[0] = 1;
    element[QUADRANGLE].arete[0].som[1] = 2;
    element[QUADRANGLE].arete[1].som[0] = 2;
    element[QUADRANGLE].arete[1].som[1] = 3;
    element[QUADRANGLE].arete[2].som[0] = 3;
    element[QUADRANGLE].arete[2].som[1] = 4;
    element[QUADRANGLE].arete[3].som[0] = 4;
    element[QUADRANGLE].arete[3].som[1] = 1;
    element[QUADRANGLE].face = (FACE *)mxCalloc(1, sizeof(FACE));
    element[QUADRANGLE].face[0].nsom = 4;
    element[QUADRANGLE].face[0].som = (int *)mxCalloc(4, sizeof(int));
    element[QUADRANGLE].face[0].som[0] = 1;
    element[QUADRANGLE].face[0].som[1] = 2;
    element[QUADRANGLE].face[0].som[2] = 3;
    element[QUADRANGLE].face[0].som[3] = 4;
    
    element[TETRAEDRE].nsom = 4;
    element[TETRAEDRE].narete = 6;
    element[TETRAEDRE].nface = 4;
    element[TETRAEDRE].arete = (ARETE *)mxCalloc(6, sizeof(ARETE));
    element[TETRAEDRE].arete[0].som[0] = 1;
    element[TETRAEDRE].arete[0].som[1] = 2;
    element[TETRAEDRE].arete[1].som[0] = 2;
    element[TETRAEDRE].arete[1].som[1] = 3;
    element[TETRAEDRE].arete[2].som[0] = 3;
    element[TETRAEDRE].arete[2].som[1] = 1;
    element[TETRAEDRE].arete[3].som[0] = 1;
    element[TETRAEDRE].arete[3].som[1] = 4;
    element[TETRAEDRE].arete[4].som[0] = 2;
    element[TETRAEDRE].arete[4].som[1] = 4;
    element[TETRAEDRE].arete[5].som[0] = 3;
    element[TETRAEDRE].arete[5].som[1] = 4;
    element[TETRAEDRE].face = (FACE *)mxCalloc(4, sizeof(FACE));
    element[TETRAEDRE].face[0].nsom = 3;
    element[TETRAEDRE].face[0].som = (int *)mxCalloc(3, sizeof(int));
    element[TETRAEDRE].face[0].som[0] = 1;
    element[TETRAEDRE].face[0].som[1] = 3;
    element[TETRAEDRE].face[0].som[2] = 2;
    element[TETRAEDRE].face[1].nsom = 3;
    element[TETRAEDRE].face[1].som = (int *)mxCalloc(3, sizeof(int));
    element[TETRAEDRE].face[1].som[0] = 1;
    element[TETRAEDRE].face[1].som[1] = 4;
    element[TETRAEDRE].face[1].som[2] = 3;
    element[TETRAEDRE].face[2].nsom = 3;
    element[TETRAEDRE].face[2].som = (int *)mxCalloc(3, sizeof(int));
    element[TETRAEDRE].face[2].som[0] = 1;
    element[TETRAEDRE].face[2].som[1] = 2;
    element[TETRAEDRE].face[2].som[2] = 4;
    element[TETRAEDRE].face[3].nsom = 3;
    element[TETRAEDRE].face[3].som = (int *)mxCalloc(3, sizeof(int));
    element[TETRAEDRE].face[3].som[0] = 2;
    element[TETRAEDRE].face[3].som[1] = 3;
    element[TETRAEDRE].face[3].som[2] = 4;
    
    element[PENTAEDRE].nsom = 6;
    element[PENTAEDRE].narete = 9;
    element[PENTAEDRE].nface = 5;
    element[PENTAEDRE].arete = (ARETE *)mxCalloc(9, sizeof(ARETE));
    element[PENTAEDRE].arete[0].som[0] = 1;
    element[PENTAEDRE].arete[0].som[1] = 2;
    element[PENTAEDRE].arete[1].som[0] = 2;
    element[PENTAEDRE].arete[1].som[1] = 3;
    element[PENTAEDRE].arete[2].som[0] = 3;
    element[PENTAEDRE].arete[2].som[1] = 1;
    element[PENTAEDRE].arete[3].som[0] = 1;
    element[PENTAEDRE].arete[3].som[1] = 4;
    element[PENTAEDRE].arete[4].som[0] = 2;
    element[PENTAEDRE].arete[4].som[1] = 5;
    element[PENTAEDRE].arete[5].som[0] = 3;
    element[PENTAEDRE].arete[5].som[1] = 6;
    element[PENTAEDRE].arete[6].som[0] = 4;
    element[PENTAEDRE].arete[6].som[1] = 5;
    element[PENTAEDRE].arete[7].som[0] = 5;
    element[PENTAEDRE].arete[7].som[1] = 6;
    element[PENTAEDRE].arete[8].som[0] = 6;
    element[PENTAEDRE].arete[8].som[1] = 4;
    element[PENTAEDRE].face = (FACE *)mxCalloc(5, sizeof(FACE));
    element[PENTAEDRE].face[0].nsom = 3;
    element[PENTAEDRE].face[0].som = (int *)mxCalloc(3, sizeof(int));
    element[PENTAEDRE].face[0].som[0] = 1;
    element[PENTAEDRE].face[0].som[1] = 3;
    element[PENTAEDRE].face[0].som[2] = 2;
    element[PENTAEDRE].face[1].nsom = 4;
    element[PENTAEDRE].face[1].som = (int *)mxCalloc(4, sizeof(int));
    element[PENTAEDRE].face[1].som[0] = 1;
    element[PENTAEDRE].face[1].som[1] = 4;
    element[PENTAEDRE].face[1].som[2] = 6;
    element[PENTAEDRE].face[1].som[3] = 3;
    element[PENTAEDRE].face[2].nsom = 4;
    element[PENTAEDRE].face[2].som = (int *)mxCalloc(4, sizeof(int));
    element[PENTAEDRE].face[2].som[0] = 1;
    element[PENTAEDRE].face[2].som[1] = 2;
    element[PENTAEDRE].face[2].som[2] = 5;
    element[PENTAEDRE].face[2].som[3] = 4;
    element[PENTAEDRE].face[3].nsom = 3;
    element[PENTAEDRE].face[3].som = (int *)mxCalloc(3, sizeof(int));
    element[PENTAEDRE].face[3].som[0] = 4;
    element[PENTAEDRE].face[3].som[1] = 5;
    element[PENTAEDRE].face[3].som[2] = 6;
    element[PENTAEDRE].face[4].nsom = 4;
    element[PENTAEDRE].face[4].som = (int *)mxCalloc(4, sizeof(int));
    element[PENTAEDRE].face[4].som[0] = 2;
    element[PENTAEDRE].face[4].som[1] = 3;
    element[PENTAEDRE].face[4].som[2] = 5; 
    element[PENTAEDRE].face[4].som[3] = 6; 
    
    element[HEXAEDRE].nsom = 8;
    element[HEXAEDRE].narete = 12;
    element[HEXAEDRE].nface = 6;
    element[HEXAEDRE].arete = (ARETE *)mxCalloc(12, sizeof(ARETE));
    element[HEXAEDRE].arete[0].som[0] = 1;
    element[HEXAEDRE].arete[0].som[1] = 2;
    element[HEXAEDRE].arete[1].som[0] = 2;
    element[HEXAEDRE].arete[1].som[1] = 3;
    element[HEXAEDRE].arete[2].som[0] = 3;
    element[HEXAEDRE].arete[2].som[1] = 4;
    element[HEXAEDRE].arete[3].som[0] = 4;
    element[HEXAEDRE].arete[3].som[1] = 1;
    element[HEXAEDRE].arete[4].som[0] = 1;
    element[HEXAEDRE].arete[4].som[1] = 5;
    element[HEXAEDRE].arete[5].som[0] = 2;
    element[HEXAEDRE].arete[5].som[1] = 6;
    element[HEXAEDRE].arete[6].som[0] = 3;
    element[HEXAEDRE].arete[6].som[1] = 7;
    element[HEXAEDRE].arete[7].som[0] = 4;
    element[HEXAEDRE].arete[7].som[1] = 8;
    element[HEXAEDRE].arete[8].som[0] = 5;
    element[HEXAEDRE].arete[8].som[1] = 6;
    element[HEXAEDRE].arete[9].som[0] = 6;
    element[HEXAEDRE].arete[9].som[1] = 7;
    element[HEXAEDRE].arete[10].som[0] = 7;
    element[HEXAEDRE].arete[10].som[1] = 8;
    element[HEXAEDRE].arete[11].som[0] = 8;
    element[HEXAEDRE].arete[11].som[1] = 5;
    element[HEXAEDRE].face = (FACE *)mxCalloc(6, sizeof(FACE));
    element[HEXAEDRE].face[0].nsom = 4;
    element[HEXAEDRE].face[0].som = (int *)mxCalloc(4, sizeof(int));
    element[HEXAEDRE].face[0].som[0] = 1;
    element[HEXAEDRE].face[0].som[1] = 4;
    element[HEXAEDRE].face[0].som[2] = 3;
    element[HEXAEDRE].face[0].som[3] = 2;
    element[HEXAEDRE].face[1].nsom = 4;
    element[HEXAEDRE].face[1].som = (int *)mxCalloc(4, sizeof(int));
    element[HEXAEDRE].face[1].som[0] = 1;
    element[HEXAEDRE].face[1].som[1] = 5;
    element[HEXAEDRE].face[1].som[2] = 8;
    element[HEXAEDRE].face[1].som[3] = 4;
    element[HEXAEDRE].face[2].nsom = 4;
    element[HEXAEDRE].face[2].som = (int *)mxCalloc(4, sizeof(int));
    element[HEXAEDRE].face[2].som[0] = 1;
    element[HEXAEDRE].face[2].som[1] = 2;
    element[HEXAEDRE].face[2].som[2] = 6;
    element[HEXAEDRE].face[2].som[3] = 5;
    element[HEXAEDRE].face[3].nsom = 4;
    element[HEXAEDRE].face[3].som = (int *)mxCalloc(4, sizeof(int));
    element[HEXAEDRE].face[3].som[0] = 5;
    element[HEXAEDRE].face[3].som[1] = 6;
    element[HEXAEDRE].face[3].som[2] = 7;
    element[HEXAEDRE].face[3].som[3] = 8;
    element[HEXAEDRE].face[4].nsom = 4;
    element[HEXAEDRE].face[4].som = (int *)mxCalloc(4, sizeof(int));
    element[HEXAEDRE].face[4].som[0] = 2;
    element[HEXAEDRE].face[4].som[1] = 3;
    element[HEXAEDRE].face[4].som[2] = 7; 
    element[HEXAEDRE].face[4].som[3] = 6;
    element[HEXAEDRE].face[5].nsom = 4;
    element[HEXAEDRE].face[5].som = (int *)mxCalloc(4, sizeof(int));
    element[HEXAEDRE].face[5].som[0] = 3;
    element[HEXAEDRE].face[5].som[1] = 4;
    element[HEXAEDRE].face[5].som[2] = 8; 
    element[HEXAEDRE].face[5].som[3] = 7;


  /********************************************************/
  /* TABLEAU NOP5 : Description sequentielle des elements */
  /********************************************************/
  fseek(fnopo, 12, SEEK_CUR);

  for(i = 0; i < MAX_FEM; i++)
    nbfem[i] = 0;

  /* Allocate memory for 'elem' */
  elem = (ELEMENT *)mxCalloc(ne, sizeof(ELEMENT));

if(VERBOSE){
  mexPrintf("Les elements :\n");
}
  for(i = 0; i < ne; i++){
if(VERBOSE){
    mexPrintf("Element %d :\n", i + 1);
}

    /* Read ncge */
    fread(&(elem[i].ncge), sizeof(int), 1, fnopo); if(endian == small) small2big(&(elem[i].ncge), sizeof(int));
if(VERBOSE){
    mexPrintf("   ncge = %d\n", elem[i].ncge);
}

    /* Read nmae */
    fread(&(elem[i].nmae), sizeof(int), 1, fnopo); if(endian == small) small2big(&(elem[i].nmae), sizeof(int));
if(VERBOSE){
    mexPrintf("   nmae = %d\n", elem[i].nmae);
}

    /* Read ndsde */
    fread(&(elem[i].ndsde), sizeof(int), 1, fnopo); if(endian == small) small2big(&(elem[i].ndsde), sizeof(int));
if(VERBOSE){
    mexPrintf("   ndsde = %d\n", elem[i].ndsde);
}
    
     /* Read nno */
    fread(&(elem[i].nno), sizeof(int), 1, fnopo); if(endian == small) small2big(&(elem[i].nno), sizeof(int));
if(VERBOSE){
    mexPrintf("   nno = %d\n", elem[i].nno);
}
    
    if(elem[i].nno != 0){
      /* Allocate memory for 'nono' */
      elem[i].nono = (int *)mxCalloc(elem[i].nno, sizeof(int));

      /* Read nono */
      fread(elem[i].nono, sizeof(int), elem[i].nno, fnopo);
      if(endian == small)
        for(j = 0; j < elem[i].nno; j++)
          small2big(elem[i].nono + j, sizeof(int));

if(VERBOSE){
      mexPrintf("      nono = ");
      for(j = 0; j < elem[i].nno; j++)
        mexPrintf("%d ", *(elem[i].nono + j));
      mexPrintf("\n");
}
    }

    if(ncopnp == 0){
      /* Read npo */
      fread(&(elem[i].npo), sizeof(int), 1, fnopo); if(endian == small) small2big(&(elem[i].npo), sizeof(int));
if(VERBOSE){
      mexPrintf("   npo = %d\n", elem[i].npo);
}
    
      if(elem[i].npo != 0){
        /* Allocate memory for 'nopo' */
        elem[i].nopo = (int *)mxCalloc(elem[i].npo, sizeof(int));

        /* Read nopo */
        fread(elem[i].nopo, sizeof(int), elem[i].npo, fnopo);
        if(endian == small)
          for(j = 0; j < elem[i].npo; j++)
            small2big(elem[i].nopo + j, sizeof(int));

if(VERBOSE){
        mexPrintf("      nopo = ");
        for(j = 0; j < elem[i].npo; j++)
          mexPrintf("%d ", *(elem[i].nopo + j));
        mexPrintf("\n");
}
      }
    }

    if(elem[i].nmae != 0){
      /* Read ining */
      fread(&(elem[i].ining), sizeof(int), 1, fnopo); if(endian == small) small2big(&(elem[i].ining), sizeof(int));
if(VERBOSE){
      mexPrintf("   ining = %d\n", elem[i].ining);
}      

      if(elem[i].nmae > 1){
        /* Allocate memory for 'ref' */
        elem[i].ref = (int *)mxCalloc(elem[i].nmae - 1, sizeof(int));
        
        /* Read ref */
        fread(elem[i].ref, sizeof(int), (elem[i].nmae - 1), fnopo);
        if(endian == small)
          for(j = 0; j < elem[i].nmae - 1; j++)
            small2big(elem[i].ref + j, sizeof(int));

if(VERBOSE){
        mexPrintf("      ref = ");
        for(j = 0; j < elem[i].nmae - 1; j++)
          mexPrintf("%d ", *(elem[i].ref + j));
        mexPrintf("\n");
}

        /* Store 'refsom' */
        if(ncopnp == 0){
          offsetref = elem[i].nmae - 1 - elem[i].npo;
          for(j = 0; j < elem[i].npo; j++)
            refsom[elem[i].nopo[j] - 1] = elem[i].ref[offsetref + j];
        }
        else{
	  offsetnode = element[elem[i].ncge-1].nsom;
	  offsetref = element[elem[i].ncge-1].narete;
	  if (elem[i].ining==1)
	    offsetref += element[elem[i].ncge-1].nface;
	  else if (elem[i].ining==3)
	    offsetref = 0;

	  if ((n1!=0)&&(elem[i].ining<3))
	    {
	      for (j=0;j<element[elem[i].ncge-1].narete;j++)
		{
		  refsom[elem[i].nono[j+offsetnode]-1] = elem[i].ref[j];
		}
	    }
	  for (j=0;j<element[elem[i].ncge-1].nsom;j++)
	    {
	      refsom[elem[i].nono[j]-1] = elem[i].ref[offsetref+j];
	    }
        }
      }
    }
    else elem[i].ining = 0;

    guess_FEM(elem+i, problem, n1, iset, iseq, isete, isepe, isehe);
    
    if(elem[i].fem == UNKNOWN){
      sprintf(errormessage, "Couldn't guess FEM for element %d\n", i+1);
      mexErrMsgTxt(errormessage);
    }
    nbfem[elem[i].fem]++;

if(VERBOSE){
    mexPrintf("FEM for element %d : %s\n", i, femname[elem[i].fem]);
}
  }
  
  /********************/
  /* READING FINISHED */
  /********************/
if(VERBOSE){
  mexPrintf("Reading finished.\n");
}

  fclose(fnopo);
if(VERBOSE){
  mexPrintf("Exporting Node and Model Description matrices...\n");
}

  if(ncopnp == 0){
    /* Il faut interpoler les noeuds qui ne sont pas sommet */
    /* Allocate memory for the node coordinates */
    node = (float *)mxCalloc(ndim*noe, sizeof(float));
    
    /* Allocate memory for 'refnode' */
    refnode = (int *)mxCalloc(noe, sizeof(int));
    for(i = 0; i < noe; i++)
      refnode[i] = 0; /* Default : no ref */

    /* Allocate memory for resi */
    resi = (float **)mxCalloc(max(n1,max(iset,max(iseq,max(isete,max(isepe,isehe))))), sizeof(float *));

    /* Allocate memory */
    for(i = 0; i < ne; i++){ 
      offsetref = (elem[i].nmae > 1)?(elem[i].nmae - 1 - elem[i].npo):0;

      /* nno >= npo */
      /* Process the first npo nodes */
      for(j = 0; j < elem[i].npo; j++){
        memcpy(node + ndim*(elem[i].nono[j] - 1), som + ndim*(elem[i].nopo[j] - 1), ndim*sizeof(float));
        if(elem[i].ining != 0)
          refnode[elem[i].nono[j] - 1] = elem[i].ref[offsetref + j];
      }

      /* Now process the rest */
      offsetnode = elem[i].npo;

      if(n1 != 0){ 
        /* noeuds supplementaires par arete */
        offsetref = (elem[i].ining == 1)?(element[elem[i].ncge - 1].nface):0;
        
        for(j = 0; j < element[elem[i].ncge - 1].narete; j++){
          for(k = 0; k < 2; k++)
            somi[k] = som + ndim*(elem[i].nopo[element[elem[i].ncge - 1].arete[j].som[k] - 1] - 1);
          
          for(k = 0; k < n1; k++){
            resi[k] = node + ndim*(elem[i].nono[offsetnode + k] - 1);
            if(elem[i].ining == 1 || elem[i].ining == 2)
              refnode[elem[i].nono[offsetnode + k] - 1] = elem[i].ref[offsetref + j];
          }
            
          interpol(ndim, 2, somi, resi, n1);
          offsetnode += n1;
        }
      }
      
      /* noeud supplementaire par face */
      for(j = 0; j < element[elem[i].ncge - 1].nface; j++)
        if(iset != 0 && element[elem[i].ncge - 1].face[j].nsom == 3){
          for(k = 0; k < 3; k++)
            somi[k] = som + ndim*(elem[i].nopo[element[elem[i].ncge - 1].face[j].som[k] - 1] - 1);
          
          for(k = 0; k < iset; k++){
            resi[k] = node + ndim*(elem[i].nono[offsetnode + k] - 1);
            if(elem[i].ining == 1)
              refnode[elem[i].nono[offsetnode + k] - 1] = elem[i].ref[j];
          }
          
          interpol(ndim, 3, somi, resi, iset);
          offsetnode += iset;
        }
        else if(iseq != 0 && element[elem[i].ncge - 1].face[j].nsom == 4){
          for(k = 0; k < 4; k++)
            somi[k] = som + ndim*(elem[i].nopo[element[elem[i].ncge - 1].face[j].som[k] - 1] - 1);
          
          for(k = 0; k < iseq; k++){
            resi[k] = node + ndim*(elem[i].nono[offsetnode + k] - 1);
            if(elem[i].ining == 1)
              refnode[elem[i].nono[offsetnode + k] - 1] = elem[i].ref[j];
          }
          
          interpol(ndim, 4, somi, resi, iseq);
          offsetnode += iseq;
        }
      
      if(isete != 0 && elem[i].ncge == TETRAEDRE - 1){ /* noeud supplementaire dans un tetraedre */
        for(k = 0; k < 4; k++)
          somi[k] = som + ndim*(elem[i].nono[k] - 1);
        
        for(k = 0; k < isete; k++)
          resi[k] = node + ndim*(elem[i].nono[offsetnode + k] - 1);
        
        interpol(ndim, 4, somi, resi, isete);
        offsetnode += isete;
      }
      else if(isepe != 0 && elem[i].ncge == PENTAEDRE - 1){ /* noeud supplementaire dans un tetraedre */
        for(k = 0; k < 6; k++)
          somi[k] = som + ndim*(elem[i].nono[k] - 1);
        
        for(k = 0; k < isepe; k++)
          resi[k] = node + ndim*(elem[i].nono[offsetnode + k] - 1);
        
        interpol(ndim, 6, somi, resi, isete);
        offsetnode += isepe;
      }
      else if(isehe != 0 && elem[i].ncge == HEXAEDRE - 1){ /* noeud supplementaire dans un tetraedre */
        for(k = 0; k < 8; k++)
          somi[k] = som + ndim*(elem[i].nono[k] - 1);
        
        for(k = 0; k < isehe; k++)
          resi[k] = node + ndim*(elem[i].nono[offsetnode + k] - 1);
        
        interpol(ndim, 8, somi, resi, isete);
        offsetnode += isehe;
      }      
    }
if(VERBOSE){
    for(i = 0; i < noe; i++){
      mexPrintf("Noeud %d : ", i + 1);
      for(j = 0; j < ndim; j++)
        mexPrintf("%f ", *(node + ndim*i + j));
      mexPrintf("\n");
    }
}
  }
  else{
    node = som;
    refnode = refsom;
  }
  
  
  /************************/
  /* Building Node matrix */
  /************************/

  /* A node matrix has seven columns. Each row of gives 
                                                      
       NodeID PID DID GID x y z                         
                                                      
    where :                                            
    - NodeId are node numbers (int > 0, no continuity),                         
    - PID and DID are coordinate systems numbers for position and displacement respectively (not currently used),
    - GID is a node group number (0 or int > 0),
    - x,y,z are the coordinates in the global coordinate system (no local, since no NASTRAN) */
  
  /* Allocate memory for node matrix */
  plhs[0] = mxCreateDoubleMatrix(noe, 7, mxREAL);
  node_matrix = (double *)mxGetPr(plhs[0]); 

  for(i = 0; i < noe; i++){
    /* Store 'numero des sommets' as NodeId */
    node_matrix[0*noe + i] = (double)(i + 1);
    
    /* Store 0 for PID and DID */
    node_matrix[1*noe + i] = (double)0;
    node_matrix[2*noe + i] = (double)0;
    
    /* Store nodes references */
    node_matrix[3*noe + i] = refnode[i];

    /* Store nodes coordinates */
    node_matrix[4*noe + i] = (double)(*(node + ndim*i));
    node_matrix[5*noe + i] = (double)(*(node + ndim*i + 1));
    node_matrix[6*noe + i] = (ndim == 2)?((double)0):((double)(*(node + ndim*i + 2)));

if(VERBOSE){
    mexPrintf("Node %d : ", i + 1);
    for(j = 0; j < 7; j++)
      mexPrintf("%lf ", node_matrix[j*noe + i]);
    mexPrintf("\n");
}
  }

  /**********************************/
  /* Build Model Description matrix */
  /**********************************/

  /* Calcul du nombre de colonnes (et de ligne) de la matrice */
  /* col = max_{j,i} (c1[j], c2[i]) */
  /* c1[j] = 1 /-1/ + strlen(femname[j]) + 1 /0/ + 1 /opt/ */
  /* c2[i] = elem[i].nno + 1 /MatID/(defaut : 0) + 1 /ProID/(defaut : 0) + sizeof(/OTHERINFO/) */
  nbrow_FEelt = 0;
  nbcol_FEelt = 0;
  for(j = 0; j < MAX_FEM; j++)
    if(nbfem[j]){
      /* Calcul de c1[j] */
      nbcol_FEelt = max(nbcol_FEelt, 1 + strlen(femname[j]) + 1 + 1);
      
      nbrow_FEelt += 1 + nbfem[j];
    }
     

  nbcol_redge = 0;
  nbcol_rface = 0;
  
  /* Calculs du nombre de colonnes pour elm_matrix, redge_matrix, rface_matrix */
  for(i = 0; i < ne; i++)
    {
      /* Calcul de c2[i] */
      nbcol_FEelt = max(nbcol_FEelt, elem[i].nno + 1 + 1 + femotherinfo[elem[i].fem]);
      /* nombre de colonnes pour redge_matrix */
      nbcol_redge = max(nbcol_redge, element[elem[i].ncge-1].narete);
      /* nombre de colonnes pour rface_matrix */
      nbcol_rface = max(nbcol_rface, element[elem[i].ncge-1].nface);
    }
  
  /* Memoire pour elm_matrix, redge_matrix, rface_matrix */
  plhs[1] = mxCreateDoubleMatrix(nbrow_FEelt, nbcol_FEelt, mxREAL);
  elm_matrix = mxGetPr(plhs[1]);
  plhs[2] = mxCreateDoubleMatrix(nbrow_FEelt, nbcol_redge, mxREAL);
  redge_matrix = mxGetPr(plhs[2]);
  plhs[3] = mxCreateDoubleMatrix(nbrow_FEelt, nbcol_rface, mxREAL);
  rface_matrix = mxGetPr(plhs[3]);
  

  row = 0;
  for(j = 0; j < MAX_FEM; j++)
    if(nbfem[j]){
      /* Il y a des elements du type femname[i] */
      /* On commence par le 'header row' */
      
      /* Inf */
      elm_matrix[0*nbrow_FEelt + row] = -1.0;
      col = 1;
      
      /* abs('femname') */
      for(k = 0; k < strlen(femname[j]); k++)
        elm_matrix[(col + k)*nbrow_FEelt + row] = (double) femname[j][k];
      col += strlen(femname[j]);

      /* 0 */
      elm_matrix[(col++)*nbrow_FEelt + row] = 0;

      /* Opt(1) = 0 */
      elm_matrix[(col++)*nbrow_FEelt + row] = 0;

      /* On finit la ligne par des zeros, si necessaire */
      for(k = col; k < nbcol_FEelt; k++)
        elm_matrix[k*nbrow_FEelt + row] = 0;

      /* pour redge_matrix */
      for (k = 0; k < nbcol_redge; k++)
	redge_matrix[row+k] = 0;
      
      /* pour rface_matrix */
      for (k = 0; k < nbcol_rface; k++)
	rface_matrix[row+k] = 0;

      row++;

      /* For each element of this type : [NodeNumbers MatId ProId OtherInfo] 
        http://www-rocq.inria.fr/modulef/Doc/FR/Guide2-14/node49.html#SECTION051140000000000000000
      */
      for(i = 0; i < ne; i++)
        if(elem[i].fem == j){
          /* NodeNumbers */
          for(k = 0; k < elem[i].nno; k++)
            elm_matrix[k*nbrow_FEelt + row] = (double) elem[i].nono[k]; 
           
          col = elem[i].nno;

          /* MatId <=> ndsde */
          elm_matrix[(col++)*nbrow_FEelt + row] = (double) elem[i].ndsde;
          
          /* ProId <=> 0 */
          elm_matrix[(col++)*nbrow_FEelt + row] = 0;

          /* OtherInfo <=> 0 */
          for(k = 0; k < femotherinfo[j]; k++)
            elm_matrix[(col + k)*nbrow_FEelt + row] = 0;
          col += femotherinfo[j];
          
          /* On finit la ligne par des zeros, si necessaire */
          for(k = col; k < nbcol_FEelt; k++)
            elm_matrix[k*nbrow_FEelt + row] = 0;

	  /* redge_matrix */
	  nfacecur = element[elem[i].ncge-1].nface;
	  nedgecur = element[elem[i].ncge-1].narete;
	
	  for(k = 0; k < nedgecur; k++)
	    if (elem[i].ining == 3)
	      redge_matrix[k*nbrow_FEelt+row] = 0;
	    else if (elem[i].ining == 2)
	      redge_matrix[k*nbrow_FEelt+row] = elem[i].ref[k];
	    else if (elem[i].ining == 1)
	      redge_matrix[k*nbrow_FEelt+row] = elem[i].ref[nfacecur + k];
	  for(k = nedgecur; k < nbcol_redge; k++)
	    redge_matrix[k*nbrow_FEelt+row] = 0;

	  /* rface_matrix */
	  for(k = 0; k < nfacecur; k++)
	    if (elem[i].ining > 1)
	      rface_matrix[k*nbrow_FEelt+row] = 0;
	    else if (elem[i].ining == 1)
	      rface_matrix[k*nbrow_FEelt+row] = elem[i].ref[k];
	  for (k = nfacecur; k < nbcol_rface; k++)
	    rface_matrix[k*nbrow_FEelt+row] = 0;
	  
          row++;
        }
    }
if(VERBOSE){
  for(i = 0; i < nbrow_FEelt; i++){
    mexPrintf("Elm %d : ", i + 1);
    for(j = 0; j < nbcol_FEelt; j++)
      mexPrintf("%lf ", elm_matrix[j*nbrow_FEelt + i]);
    mexPrintf("\n");
  }
}

if(VERBOSE){
  mexPrintf("Exporting finished.\n");
}
}
