/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   GOBNILP Copyright (C) 2012 James Cussens, Mark Bartlett             */
/*                                                                       */
/*   This program is free software; you can redistribute it and/or       */
/*   modify it under the terms of the GNU General Public License as      */
/*   published by the Free Software Foundation; either version 3 of the  */
/*   License, or (at your option) any later version.                     */
/*                                                                       */
/*   This program is distributed in the hope that it will be useful,     */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of      */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    */
/*   General Public License for more details.                            */
/*                                                                       */
/*   You should have received a copy of the GNU General Public License   */
/*   along with this program; if not, see                                */
/*   <http://www.gnu.org/licenses>.                                      */
/*                                                                       */
/*   Additional permission under GNU GPL version 3 section 7             */
/*                                                                       */
/*   If you modify this Program, or any covered work, by linking or      */
/*   combining it with SCIP (or a modified version of that library),     */
/*   containing parts covered by the terms of the ZIB Academic License,  */
/*   the licensors of this Program grant you additional permission to    */
/*   convey the resulting work.                                          */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
This file was created by editing the file probdata_lop.c that comes with the linear ordering example
in SCIP
*/

/*#define SCIP_DEBUG*/
#include <string.h>
#include <stdlib.h>

#include "probdata_bn.h"

#include "cons_dagcluster.h"
#include "scip/pub_misc.h"
#include "scip/scipdefplugins.h"

/* ----------------- SCIP interface functions ------------------------ */

/** delete problem data */

static
SCIP_DECL_PROBDELORIG(probdelorigBN)
{
   int i, k;

   assert( probdata != NULL );
   assert( *probdata != NULL );

   assert( (*probdata)->Scores != NULL );
   assert( (*probdata)->PaVars != NULL );
   assert( (*probdata)->nParents != NULL );
   assert( (*probdata)->ParentSets != NULL );
   assert( (*probdata)->nParentSets != NULL );

   for (i = 0; i < (*probdata)->n; ++i)
   {
      for (k = 0; k < (*probdata)->nParentSets[i]; ++k)
      {
       SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->PaVars[i][k]) );
       SCIPfreeMemoryArray(scip, &((*probdata)->ParentSets[i][k]));
      }
      SCIPfreeMemoryArray(scip, &(*probdata)->PaVars[i]);
      SCIPfreeMemoryArray(scip, &((*probdata)->Scores[i]));
      SCIPfreeMemoryArray(scip, &(*probdata)->nParents[i]);
      SCIPfreeMemoryArray(scip, &(*probdata)->ParentSets[i]);
   }
   SCIPfreeMemoryArray(scip, &(*probdata)->PaVars);
   SCIPfreeMemoryArray(scip, &((*probdata)->Scores));
   SCIPfreeMemoryArray(scip, &((*probdata)->nParents));
   SCIPfreeMemoryArray(scip, &((*probdata)->nParentSets));
   SCIPfreeMemoryArray(scip, &((*probdata)->ParentSets));

   // Free the sexual consistency variables if we have used them
   if ((*probdata)->SexVars != NULL ) {
      for (i = 0; i < (*probdata)->n; ++i) {
         SCIPreleaseVar(scip, &(*probdata)->SexVars[i]);
      }
      SCIPfreeMemoryArray(scip, &(*probdata)->SexVars);
   }

   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}



#define probtransBN NULL
#define probdeltransBN NULL
#define probinitsolBN NULL
#define probexitsolBN NULL
#define probcopyBN NULL

/* ----------------- auxiliary functions ------------------------ */

/** read local score BN file (in Jaakkola format)
 *
 *  Format:
 *  first line is: number of variables
 *  then local scores for each variable
 *  first line of section of local scores is "variable number_of_parent_sets"
 *  other lines are like "score 3 parent1 parent2 parent3" (e.g. when there are 3 parents)
 *  NB Variables are 0, 1 ,2 , n-1
 */

static
SCIP_RETCODE BNreadFile(
   SCIP*        scip,          /**< SCIP data structure */
   const char*  filename,      /**< name of file to read */
   SCIP_PROBDATA* probdata     /**< problem data to be filled */
   )
{
  int i, k, l;
  int i0;
  FILE *file;
  int status;
  int n;            /* number of variables */
  SCIP_Real** Scores;
  int* nParentSets;
  int** nParents;
  int*** ParentSets;
  char* pedOutput;
  SCIP_Bool usingPedigrees;

   SCIPgetStringParam(scip,"gobnilp/outputfile/pedigree", &pedOutput);
   usingPedigrees = (strcmp(pedOutput,"") != 0);

  /* open file */
  if (strcmp(filename,"-") == 0)
       file = stdin;
  else
     file = fopen(filename, "r");
  if ( file == NULL )
    {
      SCIPerrorMessage("Could not open file %s.\n", filename);
      return SCIP_NOFILE;
    }

  /* read number of elements */
  status = fscanf(file, "%d", &n);
  if ( ! status )
    {
      SCIPerrorMessage("Reading failed: first line did not state number of variables.\n");
      return SCIP_READERROR;
    }
  assert( 0 < n );
  printf("Number of variables: %d\n\n", n);
  probdata->n = n;

  SCIP_CALL( SCIPallocMemoryArray(scip, &Scores, n) );
  SCIP_CALL( SCIPallocMemoryArray(scip, &nParentSets, n) );
  SCIP_CALL( SCIPallocMemoryArray(scip, &nParents, n) );
  SCIP_CALL( SCIPallocMemoryArray(scip, &ParentSets, n) );

  probdata->Scores = Scores;
  probdata->nParentSets = nParentSets;
  probdata->nParents = nParents;
  probdata->ParentSets = ParentSets;

  for (i = 0; i < n; ++i)
    {
      status = fscanf(file, "%d %d", &i0, &(nParentSets[i]));
      if ( ! status )
     {
       SCIPerrorMessage("Reading failed: did not get number of parents for variable %d.\n", i);
       return SCIP_READERROR;
     }
      assert(i0 == i);
      /*SCIPmessagePrintInfo("%d %d\n\n", i, nParentSets[i]);*/

      SCIP_CALL( SCIPallocMemoryArray(scip, &(nParents[i]), nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(Scores[i]), nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(ParentSets[i]), nParentSets[i]) );

      for (k = 0; k < nParentSets[i]; ++k) {
   status = fscanf(file, "%lf %d", &(Scores[i][k]), &(nParents[i][k]));
   if ( ! status )
     {
       SCIPerrorMessage("Reading failed: did not get size of parent set %d for variable %d.\n", k, i);
       return SCIP_READERROR;
     }

      if (nParents[i][k] > 2 && usingPedigrees) {
         printf("Can only use pedigree output format if no node can have more than 2 parents\nDefaulting to normal output format\n");
         SCIPsetParam(scip, "gobnilp/solutionformat", (void *)"normal");
      }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(ParentSets[i][k]), nParents[i][k]) );

   for (l = 0; l < nParents[i][k]; ++l)
   {
      status = fscanf(file, "%d", &(ParentSets[i][k][l]));
      if ( ! status )
      {
         SCIPerrorMessage("Reading failed: did not get parent %d for parent set %d for variable %d.\n", l, k, i);
         return SCIP_READERROR;
      }
   }
   SCIPsortInt(ParentSets[i][k],nParents[i][k]);
      }
    }

   fclose( file );

   printf("File reading successful \n");

   return SCIP_OKAY;
}

static
SCIP_RETCODE immorality_constraint(
   SCIP* scip,
   int i,
   int j,
   int child,
   SCIP_Bool truthvalue
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS* cons;

   int k,l;

   int found;

   char s[SCIP_MAXSTRLEN];

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

  /* scip will upgrade from linear to specialised */
   if ( truthvalue )
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "user_immorality_yes#%d#%d#%d", child, i, j);
   else
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "user_immorality_no#%d#%d#%d",  child, i, j);
   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, s, 0, NULL, NULL,
               truthvalue ? 1 : -SCIPinfinity(scip),
               truthvalue ? SCIPinfinity(scip) : 0,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   for (k = 0; k < probdata->nParentSets[child]; ++k)
   {
      found = 0;
      for ( l = 0; l < probdata->nParents[child][k]; ++l )
      {
    if ( probdata->ParentSets[child][k][l] == i || probdata->ParentSets[child][k][l] == j )
       found++;

    if ( found == 2 )
    {
       SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[child][k], 1) );
       break;
    }
      }
   }

   for (k = 0; k < probdata->nParentSets[i]; ++k)
      for ( l = 0; l < probdata->nParents[i][k]; ++l )
    if ( probdata->ParentSets[i][k][l] == j )
    {
       SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][k], -1) );
       break;
    }

   for (k = 0; k < probdata->nParentSets[j]; ++k)
      for ( l = 0; l < probdata->nParents[j][k]; ++l )
    if ( probdata->ParentSets[j][k][l] == i )
    {
       SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[j][k], -1) );
       break;
    }

   SCIP_CALL( SCIPaddCons(scip, cons) );
   /*SCIP_CALL( SCIPprintCons(scip, cons, NULL) );*/
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}

static
SCIP_Bool differ(
   SCIP* scip,
   int i,
   int ki1,
   int ki2,
   int j
   )
{
   int big,small,l_big,l_small;
   SCIP_PROBDATA* probdata;

   int nParents_diff;

  /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   nParents_diff = probdata->nParents[i][ki1] - probdata->nParents[i][ki2];

   if ( nParents_diff == 1 )
   {
      small = ki2;
      big = ki1;
   }
   else if ( nParents_diff == -1 )
   {
      small = ki1;
      big = ki2;
   }
   else
      return FALSE;

   l_small = 0;
   for (l_big = 0;  l_big < probdata->nParents[i][big]; ++l_big )
   {
      if ( probdata->ParentSets[i][big][l_big] == j )
    continue;

      if ( l_small < probdata->nParents[i][small] && probdata->ParentSets[i][big][l_big] == probdata->ParentSets[i][small][l_small] )
    l_small++;
      else
    return FALSE;
   }
   return TRUE;
}

static
SCIP_RETCODE edge_constraint(
   SCIP* scip,
   int i,
   int j,
   SCIP_Bool undirected,
   SCIP_Bool truthvalue
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS* cons;

   int k,l;

   char s[SCIP_MAXSTRLEN];

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   /* scip will upgrade from linear to specialised */
   if ( undirected )
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "user_edge#%d#%d", i, j);
   else
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "user_arrow#%d#%d", i, j);
   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, s, 0, NULL, NULL,
               truthvalue,
               truthvalue,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   for (k = 0; k < probdata->nParentSets[i]; ++k)
      for ( l = 0; l < probdata->nParents[i][k]; ++l )
    if ( probdata->ParentSets[i][k][l] == j )
    {
       SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][k], 1) );
       break;
    }

   if ( undirected )
      for (k = 0; k < probdata->nParentSets[j]; ++k)
    for ( l = 0; l < probdata->nParents[j][k]; ++l )
       if ( probdata->ParentSets[j][k][l] == i )
       {
          SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[j][k], 1) );
          break;
       }

   SCIP_CALL( SCIPaddCons(scip, cons) );
   /*SCIP_CALL( SCIPprintCons(scip, cons, NULL) );*/
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   return SCIP_OKAY;
}



static
SCIP_RETCODE process_constraint(
   SCIP* scip,
   const char* line
   )
{

   int i,j,child;

   if ( line[0] == '#' )
      return SCIP_OKAY;

   if ( sscanf(line,"%d-%d",&i,&j) == 2 )
      edge_constraint(scip,i,j,TRUE,TRUE);
   else if ( sscanf(line,"~%d-%d",&i,&j) == 2 )
      edge_constraint(scip,i,j,TRUE,FALSE);
   else if ( sscanf(line,"%d<-%d",&i,&j) == 2 )
      edge_constraint(scip,i,j,FALSE,TRUE);
   else if ( sscanf(line,"~%d<-%d",&i,&j) == 2 )
      edge_constraint(scip,i,j,FALSE,FALSE);
   else if ( sscanf(line,"%d->%d<-%d",&i,&child,&j) == 3 )
      immorality_constraint(scip,i,j,child,TRUE);
   else if ( sscanf(line,"~%d->%d<-%d",&i,&child,&j) == 3 )
      immorality_constraint(scip,i,j,child,FALSE);
   else
   {
      SCIPerrorMessage("Not recognised as a DAG constraint: %s\n",line);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/*        this function copied from LOP example written by March Pfetsch */
/** get problem name
 *
 *  Returns NULL on error
 */
static
SCIP_RETCODE getProblemName(
   const char* filename,         /**< input filename */
   char*       probname,         /**< output problemname */
   int         maxSize           /**< maximum size of probname */
   )
{
   int i = 0;
   int j = 0;
   int l;

   /* first find end of string */
   while ( filename[i] != 0)
      ++i;
   l = i;

   /* go back until '.' or '/' or '\' appears */
   while ((i > 0) && (filename[i] != '.') && (filename[i] != '/') && (filename[i] != '\\'))
      --i;

   /* if we found '.', search for '/' or '\\' */
   if (filename[i] == '.')
   {
      l = i;
      while ((i > 0) && (filename[i] != '/') && (filename[i] != '\\'))
    --i;
   }

   /* correct counter */
   if ((filename[i] == '/') || (filename[i] == '\\'))
      ++i;

   /* copy name */
   while ( (i < l) && (filename[i] != 0) )
   {
      probname[j++] = filename[i++];
      if (j > maxSize-1)
    return SCIP_ERROR;
   }
   probname[j] = 0;

   return SCIP_OKAY;
}



/* ----------------- public interface functions ------------------------ */


/** create BN learning problem instance */
SCIP_RETCODE BNcreateProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of file to read */
   )
{
   SCIP_PROBDATA* probdata = NULL;
   char probname[SCIP_MAXSTRLEN];

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, &probdata) );

   /* take filename as problem name */
   SCIP_CALL( getProblemName(filename, probname, SCIP_MAXSTRLEN) );

   printf("File name:\t\t%s\n", filename);
   printf("Problem name:\t\t%s\n", probname);

   /* read file */
   SCIP_CALL( BNreadFile(scip, filename, probdata) );
   probdata->PaVars = NULL;
   probdata->SexVars = NULL;

   SCIP_CALL( SCIPcreateProb(scip, probname, probdelorigBN, probtransBN, probdeltransBN,
    probinitsolBN, probexitsolBN, probcopyBN, probdata) );

   return SCIP_OKAY;
}

SCIP_RETCODE BNgenerateModel(
   SCIP*                 scip               /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS* cons;
   int i;        /* indexes variables */
   int k;        /* indexes parentsets for a particular variable */
   int l;

   char s[SCIP_MAXSTRLEN];
   char tmp[SCIP_MAXSTRLEN];

   int n;

   SCIP_Bool noimmoralities;
   SCIP_Bool orderedcoveredarcs;
   SCIP_Bool implicitfounders;
   /* SCIP_Bool sosscores; */
   /* SCIP_Bool sosnparents; */
   SCIP_Bool sexuallyconsistent;

   int small_i, big_i, small_j, big_j;



   int ki1, ki2, kj1, kj2;
   int j, jj, l2;
   SCIP_Bool ok2;
   SCIP_VAR* arc_tmp[2];

   int minfounders;
   int maxfounders;

   int minedges;
   int maxedges;

   char* dagconstraintsfile;
   /*char* adhocconstraintsfile;*/

   FILE* dagconstraints;
   /*FILE* adhocconstraints;*/

   int status;

   SCIPgetBoolParam(scip,"gobnilp/noimmoralities",&noimmoralities);
   SCIPgetBoolParam(scip,"gobnilp/orderedcoveredarcs",&orderedcoveredarcs);
   SCIPgetBoolParam(scip,"gobnilp/implicitfounders",&implicitfounders);
   /* SCIPgetBoolParam(scip,"gobnilp/sosscores",&sosscores); */
   /* SCIPgetBoolParam(scip,"gobnilp/sosnparents",&sosnparents); */
   SCIPgetBoolParam(scip,"gobnilp/sexconsistent",&sexuallyconsistent);

   SCIPgetIntParam(scip,"gobnilp/minfounders",&minfounders);
   SCIPgetIntParam(scip,"gobnilp/maxfounders",&maxfounders);

   SCIPgetIntParam(scip,"gobnilp/minedges",&minedges);
   SCIPgetIntParam(scip,"gobnilp/maxedges",&maxedges);

   SCIPgetStringParam(scip,"gobnilp/dagconstraintsfile",&dagconstraintsfile);
   /*SCIPgetStringParam(scip,"gobnilp/adhocconstraintsfile",&adhocconstraintsfile);*/


   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   n = probdata->n;

   /* generate variables */

   // Create sexual consistency variables if necessary
   if (sexuallyconsistent) {
      SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->SexVars, n) );
      for (i = 0; i < n; i++) {
         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "isFemale(%d)", i);
         SCIP_CALL( SCIPcreateVar(scip, &(probdata->SexVars[i]), s, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, probdata->SexVars[i]) );
      }
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &probdata->PaVars, n) );
   for (i = 0; i < n; ++i)
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->PaVars[i]), probdata->nParentSets[i]) );

      /* sort, best parent set first */
      /* for primal heuristics */
      for (k=0; k < probdata->nParentSets[i]; ++k)
    probdata->Scores[i][k] = -probdata->Scores[i][k];

      SCIPsortRealPtrPtrInt(probdata->Scores[i],(void**)probdata->PaVars[i],
             (void**)probdata->ParentSets[i],probdata->nParents[i],probdata->nParentSets[i]);

      for (k=0; k < probdata->nParentSets[i]; ++k)
    probdata->Scores[i][k] = -probdata->Scores[i][k];


      for (k = 0; k < probdata->nParentSets[i]; ++k)
   {
     (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "I(%d<-{", i);
     for (l = 0; l < (probdata->nParents[i][k]); ++l)
       {
         (void) SCIPsnprintf(tmp, SCIP_MAXSTRLEN, "%d,", probdata->ParentSets[i][k][l]);
         (void) strcat(s, tmp);
       }
     (void) SCIPsnprintf(tmp, SCIP_MAXSTRLEN, "})");
     (void) strcat(s, tmp);

     /* still create variable for empty parent set */
     /* but if implicit founders it appears in no constraints and has zero objective coefficient */

     if ( implicitfounders )
        SCIP_CALL( SCIPcreateVar(scip, &(probdata->PaVars[i][k]), s, 0.0, 1.0, (probdata->Scores[i][k])-(probdata->Scores[i][(probdata->nParentSets[i])-1]),
                  SCIP_VARTYPE_BINARY,
                  TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));
     else
        SCIP_CALL( SCIPcreateVar(scip, &(probdata->PaVars[i][k]), s, 0.0, 1.0, probdata->Scores[i][k],
                  SCIP_VARTYPE_BINARY,
                  TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));

     SCIP_CALL( SCIPaddVar(scip, probdata->PaVars[i][k]) );
     SCIPdebugMessage("adding variable %s with obj coefficient %f\n", SCIPvarGetName(probdata->PaVars[i][k]), SCIPvarGetObj(probdata->PaVars[i][k]));

   }


      /* constraint that at most one parent set chosen for variable i */

      if ( implicitfounders )
      {
    (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "setpack#%d", i);
    SCIP_CALL( SCIPcreateConsSetpack(scip,&cons,s,probdata->nParentSets[i],probdata->PaVars[i],
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      }
      else
      {
    (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "setpart#%d", i);
    SCIP_CALL( SCIPcreateConsSetpart(scip,&cons,s,probdata->nParentSets[i],probdata->PaVars[i],
                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      }
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* constraint that at most one parent set chosen for variable i */
      /* SOS versions */
      /* if ( sosscores ) */
      /* { */
      /*     (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sos#%d", i); */
      /*     SCIP_CALL( SCIPcreateConsSOS1(scip,&cons,s,probdata->nParentSets[i],probdata->PaVars[i], probdata->Scores[i], */
      /*                    TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) ); */
      /*     SCIP_CALL( SCIPaddCons(scip, cons) ); */
      /*     SCIP_CALL( SCIPreleaseCons(scip, &cons) ); */
      /* } */

      /* if ( sosnparents ) */
      /* { */
      /*     (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sos#%d", i); */
      /*     SCIP_CALL( SCIPcreateConsSOS1(scip,&cons,s,probdata->nParentSets[i],probdata->PaVars[i],(double *) probdata->nParents[i], */
      /*                    TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) ); */
      /*     SCIP_CALL( SCIPaddCons(scip, cons) ); */
      /*     SCIP_CALL( SCIPreleaseCons(scip, &cons) ); */
      /* } */



   }

   if ( orderedcoveredarcs )
   {
      /*

    rule out covered arcs from higher to lower vertex

    for each i,j,C (j>i,C a set) such that the following variables exist:
    I(i<-C), I(i<-C+j), I(j<-C), I(j<-C+i)

    if I(i<-C+j)=1,I(j<-C)=1
    then the edge i<-j is 'covered' and
    can be replaced by:
    I(i<-C)=1,I(j<-C+i)
    which reverse the edge


    without creating a cycle or changing Markov equivalence
    class, so add constraint:
    I(i<-C+j) + I(j<-C) <= 1
    to rule out first case
      */

      for (i = 0; i < n; ++i)
      {
    for (ki1 = 0; ki1 < probdata->nParentSets[i]; ++ki1)
    {
       for (ki2 = ki1+1; ki2 < probdata->nParentSets[i]; ++ki2)
       {
          if ( probdata->nParents[i][ki1] - probdata->nParents[i][ki2] == 1 )
          {
        big_i = ki1;
        small_i = ki2;
          }
          else if ( probdata->nParents[i][ki2] - probdata->nParents[i][ki1] == 1 )
          {
        big_i = ki2;
        small_i = ki1;
          }
          else
        /* need to at least find two parent sets for i
           differing in size by 1 */
        continue;

          for (j = i+1; j < n; ++j)
          {
        /* ki1 and ki2 can only differ by j */
        if ( !differ(scip,i,ki1,ki2,j) )
           continue;

        for (kj1 = 0; kj1 < probdata->nParentSets[j]; ++kj1)
        {
           for (kj2 = kj1+1; kj2 < probdata->nParentSets[j]; ++kj2)
           {
         if ( probdata->nParents[j][kj1] - probdata->nParents[j][kj2] == 1 )
         {
            small_j = kj2;
            big_j = kj1;
         }
         else if ( probdata->nParents[j][kj2] - probdata->nParents[j][kj1] == 1 )
         {
            small_j = kj1;
            big_j = kj2;
         }
         else
            /* need to at least find two parent sets for j
               differing in size by 1 */
            continue;

         if ( !differ(scip,j,kj1,kj2,i) )
            continue;

         /* small_i and small_j  must be the same */
         /* if so big_i is small_j with j added */
         if (  probdata->nParents[i][small_i] !=  probdata->nParents[j][small_j] )
            continue;

         l2 = 0;
         ok2 = TRUE;
         for (l = 0; l < probdata->nParents[i][small_i]; ++l)
         {
            if ( probdata->ParentSets[i][small_i][l] != probdata->ParentSets[j][small_j][l2] )
            {
               ok2 = FALSE;
               break;
            }
            else
               l2++;
         }

         if ( !ok2 )
            continue;

         SCIPdebugMessage("Ruling out having both %s and %s since arc %d<-%d is covered and there exists %s and %s\n",
                SCIPvarGetName(probdata->PaVars[i][big_i]),SCIPvarGetName(probdata->PaVars[j][small_j]),i,j,
                SCIPvarGetName(probdata->PaVars[i][small_i]),SCIPvarGetName(probdata->PaVars[j][big_j]));
         arc_tmp[0] = probdata->PaVars[i][big_i];
         arc_tmp[1] = probdata->PaVars[j][small_j];
         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "covered_arc#%d#%d", i, j);
         SCIP_CALL( SCIPcreateConsSetpack(scip,&cons,s,2,arc_tmp,
                      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         /*SCIP_CALL( SCIPprintCons(scip, cons, NULL) );*/
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
           }
        }
          }
       }
    }
      }
   }
         /* post constraint */

         /* do ki1, ki2, kj1, kj2 lead to a constraint? */


    /* { */

    /* /\* check each as a potential 'C' *\/ */
    /* /\* small_i = 'k' *\/ */
    /* for (small_i = 0; small_i < probdata->nParentSets[i]; ++small_i) */
    /* { */

    /*    c = probdata->ParentSets[i][small_i]; */
    /*    csize = probdata->nParents[i][small_i]; /\* can be zero *\/ */
    /*    /\* is parent set k a suitable 'C'? *\/ */



    /*    { */
    /*       if ( probdata->nParents[i][kk] != csize+1 ) */
    /*     continue; */


    /*    /\* is there another parent set of i which is C+j?*\/ */

    /*    { */
    /*       big_i = -1; */




    /*     ok2 = TRUE; */
    /*     l2 = 0; */
    /*     for ( l = 0; l < probdata->nParents[i][kk]; ++l ) */
    /*     { */
    /*        /\*printf("%d,%d,%d\n",l,l2,csize);*\/ */
    /*        if ( probdata->ParentSets[i][kk][l] == j ) */
    /*      continue; */

    /*        if ( csize == 0 || probdata->ParentSets[i][kk][l] != c[l2] ) */
    /*        { */
    /*      ok2 = FALSE; */
    /*      break; */
    /*        } */
    /*        l2++; */
    /*     } */
    /*     if ( ok2 ) */
    /*     { */
    /*        big_i = kk; */
    /*        break; */
    /*     } */
    /*     /\* otherwise keep looking *\/ */
    /*       } */
    /*       if ( big_i < 0 ) */
    /*     /\* couldn't find one *\/ */
    /*       { */
    /*     /\*printf("No arc-reversal constraint for i=%d, j=%d, c=%s.\n",i,j,SCIPvarGetName(probdata->PaVars[i][small_i]));*\/ */
    /*     continue; */
    /*       } */

    /*       /\* are C and C+i parent sets of j ? *\/ */
    /*       small_j = -1; */
    /*       big_j = -1; */
    /*       for (kk = 0; kk < probdata->nParentSets[j]; ++kk) */
    /*       { */
    /*     if ( small_j < 0 && probdata->nParents[j][kk] == csize ) */
    /*     { */
    /*        /\* need to see whether it is 'c' *\/ */

    /*        ok2 = TRUE; */
    /*        l2 = 0; */
    /*        for ( l = 0; l < probdata->nParents[j][kk]; ++l ) */
    /*        { */
    /*      if ( probdata->ParentSets[j][kk][l] != c[l] ) */
    /*      { */
    /*         ok2 = FALSE; */
    /*         break; */
    /*      } */
    /*      l2++; */
    /*        } */
    /*        if ( ok2 ) */
    /*      small_j = kk; */
    /*     } */

    /*     if ( big_j < 0 && probdata->nParents[j][kk] == csize+1 ) */
    /*     { */
    /*         /\* need to see whether it is 'c'+i *\/ */

    /*        ok2 = TRUE; */
    /*        l2 = 0; */
    /*        for ( l = 0; l < probdata->nParents[j][kk]; ++l ) */
    /*        { */
    /*      if ( probdata->ParentSets[j][kk][l] == i ) */
    /*         continue; */

    /*      if ( csize == 0 || probdata->ParentSets[j][kk][l] != c[l2] ) */
    /*      { */
    /*         ok2 = FALSE; */
    /*         break; */
    /*      } */
    /*      l2++; */
    /*        } */
    /*        if ( ok2 ) */
    /*      big_j = kk; */
    /*     } */

    /*     if ( small_j>0 && big_j>0 ) */
    /*        break; */
    /*       } */

   /*           if ( small_j>0 && big_j>0 ) */
   /*           { */
   /*         /\*printf("Found arc-reversal constraint for i=%d, j=%d, c=%s.\n",i,j,SCIPvarGetName(probdata->PaVars[i][small_i]));*\/ */
   /*         /\*ok, constraint can be posted *\/ */
   /*         /\*printf("%d,%d,%d,%d\n",small_i,big_i,small_j,big_j);*\/ */
   /*         arc_tmp[0] = probdata->PaVars[i][big_i];   */
   /*         arc_tmp[1] = probdata->PaVars[j][small_j]; */
   /*         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "covered_arc#%d#%d", i, j); */
   /*         SCIP_CALL( SCIPcreateConsSetpack(scip,&cons,s,2,arc_tmp, */
   /*                      TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) ); */
   /*         SCIP_CALL( SCIPaddCons(scip, cons) ); */
   /*         /\*SCIP_CALL( SCIPprintCons(scip, cons, NULL) );*\/ */
   /*         SCIP_CALL( SCIPreleaseCons(scip, &cons) ); */
   /*           } */
   /*           /\*printf("No arc-reversal constraint for i=%d, j=%d, c=%s.\n",i,j,SCIPvarGetName(probdata->PaVars[i][small_i]));*\/ */
   /*        } */
   /*     } */
   /*    } */
   /* } */

   SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "edges", 0, NULL, NULL,
               minedges,
               maxedges > 0 ? maxedges : SCIPinfinity(scip),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
   for (i = 0; i < n; ++i)
      for (k=0; k < probdata->nParentSets[i]; ++k)
    SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][k], probdata->nParents[i][k]) );

   SCIP_CALL( SCIPaddCons(scip, cons) );
   /*SCIP_CALL( SCIPprintCons(scip, cons, NULL) );*/
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* add in constraint on number of founders */
   if ( implicitfounders )
   {
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "founders", 0, NULL, NULL,
                  minfounders-n,
                  (maxfounders > 0 ? maxfounders : SCIPinfinity(scip))-n,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      for (i = 0; i < n; ++i)
    for (k=0; k < probdata->nParentSets[i]; ++k)
       if ( probdata->nParents[i][k] > 0 )
          SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][k], -1) );
   }
   else
   {
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "founders", 0, NULL, NULL,
                  minfounders,
                  maxfounders > 0 ? maxfounders : SCIPinfinity(scip),
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      for (i = 0; i < n; ++i)
    for (k=0; k < probdata->nParentSets[i]; ++k)
       if ( probdata->nParents[i][k] == 0 )
       {
          SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][k], 1) );
          break;
       }
   }

   SCIP_CALL( SCIPaddCons(scip, cons) );
   /*SCIP_CALL( SCIPprintCons(scip, cons, NULL) );*/
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* 2 clique cuts ...*/

   /* for (i = 0; i < n; ++i) */
   /* { */
   /*    for (j = i+1; j < n; ++j) */
   /*    { */
    /* SCIP_CALL( SCIPcreateConsLinear(scip, &cons, "todo", 0, NULL, NULL, */
    /*             -SCIPinfinity(scip), */
    /*             1, */
    /*             TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) ); */

    /* for (k=0; k < probdata->nParentSets[i]; ++k) */
    /* { */
    /*    for ( l=0; l<probdata->nParents[i][k]; ++l ) */
    /*    { */
    /*       if ( probdata->ParentSets[i][k][l] == j ) */
    /*       { */
    /*     SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][k], 1) ); */
    /*     break; */
    /*       } */
    /*    } */
    /* } */
    /* for (k=0; k < probdata->nParentSets[j]; ++k) */
    /* { */
    /*    for ( l=0; l<probdata->nParents[j][k]; ++l ) */
    /*    { */
    /*       if ( probdata->ParentSets[j][k][l] == i ) */
    /*       { */
    /*     SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[j][k], 1) ); */
    /*     break; */
    /*       } */
    /*    } */
    /* } */
    /* SCIP_CALL( SCIPaddCons(scip, cons) ); */
    /* /\*SCIP_CALL( SCIPprintCons(scip, cons, NULL) );*\/ */
    /* SCIP_CALL( SCIPreleaseCons(scip, &cons) ); */


   /*     for (jj = j+1; jj < n; ++jj) */
   /*     { */
   /*        (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "clique#%d#%d#%d", i, j, jj); */
   /*        SCIP_CALL( SCIPcreateConsLinear(scip, &cons, s, 0, NULL, NULL, */
   /*                    -SCIPinfinity(scip), */
   /*                    2, */
   /*                    TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) ); */
   /*        ok2=FALSE; */
   /*        for (k=0; k < probdata->nParentSets[i]; ++k) */
   /*        { */
   /*           found = 0; */
   /*           for ( l=0; l<probdata->nParents[i][k]; ++l ) */
   /*         if ( probdata->ParentSets[i][k][l] == j || probdata->ParentSets[i][k][l] == jj ) */
   /*            found++; */
   /*           if ( found == 2 ) */
   /*         ok2 = TRUE; */
   /*           SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][k], found) ); */
   /*        } */

   /*        for (k=0; k < probdata->nParentSets[j]; ++k) */
   /*        { */
   /*           found = 0; */
   /*           for ( l=0; l<probdata->nParents[j][k]; ++l ) */
   /*         if ( probdata->ParentSets[j][k][l] == i  || probdata->ParentSets[j][k][l] == jj ) */
   /*            found++; */
   /*           if ( found == 2 ) */
   /*         ok2=TRUE; */
   /*           SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[j][k], found) ); */
   /*        } */

   /*        for (k=0; k < probdata->nParentSets[jj]; ++k) */
   /*        { */
   /*           found = 0; */
   /*           for ( l=0; l<probdata->nParents[jj][k]; ++l ) */
   /*         if ( probdata->ParentSets[jj][k][l] == i  || probdata->ParentSets[jj][k][l] == j ) */
   /*            found++; */
   /*           if ( found == 2 ) */
   /*         ok2 = TRUE;  */
   /*           SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[jj][k], found) ); */
   /*        } */
   /*        if (!ok2) */
   /*           continue; */
   /*        SCIP_CALL( SCIPaddCons(scip, cons) ); */
   /*        SCIP_CALL( SCIPprintCons(scip, cons, NULL) ); */
   /*        SCIP_CALL( SCIPreleaseCons(scip, &cons) ); */
   /*     } */
   /*    } */
   /* } */




   /* generate DAG cluster constraint */
   SCIP_CALL( SCIPcreateConsDagcluster(
                   scip,
                   &cons,
                   "DagCluster",
                   probdata->n,
                   probdata->nParentSets,
                   probdata->nParents,
                   probdata->ParentSets,
                   probdata->PaVars,
                   TRUE,
                   TRUE,
                   TRUE,
                   TRUE,
                   TRUE,
                   FALSE,
                   FALSE,
                   FALSE,
                   FALSE,
                   FALSE
                   ));
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );


   if ( strcmp(dagconstraintsfile,"") != 0 )
   {
      dagconstraints = fopen(dagconstraintsfile, "r");
      if ( dagconstraints == NULL )
      {
    SCIPerrorMessage("Could not open file %s.\n", dagconstraintsfile);
    return SCIP_NOFILE;
      }
      status = fscanf(dagconstraints,"%[^\n]%*c", s);
      while ( status == 1 )
      {
    process_constraint(scip,s);
    status = fscanf(dagconstraints,"%[^\n]%*c", s);
      }
      fclose(dagconstraints);
   }

   if ( noimmoralities )
   {
      for (i = 0; i < n; ++i)
      {
    for (j = 0; j < n; ++j)
    {
       if ( i == j )
          continue;

       for (jj = j+1; jj < n; ++jj)
       {
          if ( i == jj  )
        continue;

          SCIP_CALL( immorality_constraint(scip,j,jj,i,FALSE) );
       }
    }
      }
   }

     // Create constraints for sexual consistency if necessary
   if (sexuallyconsistent) {
      for (i = 0; i < n; i++)
         for (j = 0; j < probdata->nParentSets[i]; j++)
       if (probdata->nParents[i][j] > 2) {
          SCIPerrorMessage("Parent set %d for variable %d has more than two parents and yet gobnilp/sexconsistent is set to true.\n", j, i);
          return SCIP_ERROR;
       }
            else if (probdata->nParents[i][j] == 2) {
               SCIPdebugMessage("Creating sexual consistency constraints for %s\n", SCIPvarGetName(probdata->PaVars[i][j]));
               SCIPdebugMessage("   %s\n", SCIPvarGetName(probdata->PaVars[i][j]));
               SCIPdebugMessage("   %s\n", SCIPvarGetName(probdata->SexVars[probdata->ParentSets[i][j][1]]));
               SCIPdebugMessage("   %s\n", SCIPvarGetName(probdata->SexVars[probdata->ParentSets[i][j][0]]));
               (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sex_consistency_1_for_%s", SCIPvarGetName(probdata->PaVars[i][j]));
               SCIP_CALL( SCIPcreateConsLinear(scip, &cons, s, 0, NULL, NULL, -SCIPinfinity(scip), 2, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][j], 1) );
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->SexVars[probdata->ParentSets[i][j][1]], 1) );
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->SexVars[probdata->ParentSets[i][j][0]], 1) );
               SCIP_CALL( SCIPaddCons(scip, cons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "sex_consistency_2_for_%s", SCIPvarGetName(probdata->PaVars[i][j]));
               SCIP_CALL( SCIPcreateConsLinear(scip, &cons, s, 0, NULL, NULL, -SCIPinfinity(scip), 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][j], 1) );
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->SexVars[probdata->ParentSets[i][j][1]], -1) );
               SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->SexVars[probdata->ParentSets[i][j][0]], -1) );
               SCIP_CALL( SCIPaddCons(scip, cons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            }
   }

   /* set maximization */
   SCIP_CALL_ABORT( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );

   return SCIP_OKAY;
}

/** Prints the solution in the traditional GOBNILP format.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param probdata The problem data used by the solution.
 *  @param sol The solution to print.
 *  @param stream Where to print the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionLegacyFormat(SCIP* scip, SCIP_PROBDATA* probdata, SCIP_SOL* sol, FILE* stream) {
   int i,k,l;
   SCIP_Real val = 0.0;
   SCIP_Real bn_score = 0.0;
   SCIP_Bool no_parents;

   for (i = 0; i < probdata->n; ++i) {
      no_parents = TRUE;
      fprintf(stream, "%d<-",i);
      for (k = 0; k < probdata->nParentSets[i]; ++k) {
         val = SCIPgetSolVal(scip, sol, probdata->PaVars[i][k]);
         assert( SCIPisIntegral(scip, val) );
         if ( val > 0.5 ) {
            bn_score = bn_score + probdata->Scores[i][k];
            for ( l=0; l<probdata->nParents[i][k]; ++l )
               fprintf(stream, "%d,",probdata->ParentSets[i][k][l]);
            fprintf(stream, " %f",probdata->Scores[i][k]);
            no_parents = FALSE;
            break;
         }
      }
      if ( no_parents ) {
         fprintf(stream, " ,");
         fprintf(stream, " %f",probdata->Scores[i][probdata->nParentSets[i]-1]);
         bn_score = bn_score + probdata->Scores[i][probdata->nParentSets[i]-1];
      }
      fprintf(stream, "\n");
   }
   fprintf(stream, "BN score is %f\n",bn_score);
   return SCIP_OKAY;
}

/** Prints the solution in a Bayesian network format.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param probdata The problem data used by the solution.
 *  @param sol The solution to print.
 *  @param stream Where to print the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionBNFormat(SCIP* scip, SCIP_PROBDATA* probdata, SCIP_SOL* sol, FILE* stream) {
   int i,k,l;
   SCIP_Real val = 0.0;
   SCIP_Real bn_score = 0.0;
   SCIP_Bool no_parents;

   for (i = 0; i < probdata->n; ++i) {
      no_parents = TRUE;
      fprintf(stream, "%d <- ",i);
      for (k = 0; k < probdata->nParentSets[i]; ++k) {
         val = SCIPgetSolVal(scip, sol, probdata->PaVars[i][k]);
         assert( SCIPisIntegral(scip, val) );
         if ( val > 0.5 ) {
            fprintf(stream, "{");
            bn_score = bn_score + probdata->Scores[i][k];
            for ( l=0; l<probdata->nParents[i][k]-1; ++l ) {
               fprintf(stream, "%d,",probdata->ParentSets[i][k][l]);
            }
            if (probdata->nParents[i][k] > 0)
               fprintf(stream, "%d",probdata->ParentSets[i][k][probdata->nParents[i][k]-1]);
            fprintf(stream, "}\t%f",probdata->Scores[i][k]);
            no_parents = FALSE;
            break;
         }
      }
      if ( no_parents ) {
         fprintf(stream, "{}\t%f",probdata->Scores[i][probdata->nParentSets[i]-1]);
         bn_score = bn_score + probdata->Scores[i][probdata->nParentSets[i]-1];
      }
      fprintf(stream, "\n");
   }
   fprintf(stream, "BN score is %f\n",bn_score);
   return SCIP_OKAY;
}

/** Prints the solution as a pedigree.
 *
 *  The pedigree consists of three columns.  The first is the individual, the second is its father and the third is its mother.
 *  If either the father or the mother is not present in the sample, a '-' is printed instead.  If sex consistency is enforced,
 *  then individuals will not appear as father of one individual and mother of another.  The pedigree is sorted such that all
 *  individuals are declared on a earlier line than those in which they appear as parents.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param probdata The problem data used by the solution.
 *  @param sol The solution to print.
 *  @param stream Where to print the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionPedigreeFormat(SCIP* scip, SCIP_PROBDATA* probdata, SCIP_SOL* sol, FILE* stream) {
   int i,k;
   SCIP_Real val = 0.0;
   SCIP_Bool no_parents;
   SCIP_Bool Done[probdata->n];
   int numDone = 0;
   SCIP_Bool sexconsistent;
   char sex;

   SCIPgetBoolParam(scip,"gobnilp/sexconsistent", &sexconsistent);

   for (i = 0; i < probdata->n; i++)
      Done[i] = FALSE;

   while (numDone < probdata->n) {
      for (i = 0; i < probdata->n; ++i) {
         if (Done[i] == FALSE) {
            no_parents = TRUE;
            if (sexconsistent == FALSE)
               sex = 'U';
            else if (SCIPgetSolVal(scip, sol, probdata->SexVars[i]) > 0.5)
               sex = 'F';
            else
               sex = 'M';
            for (k = 0; k < probdata->nParentSets[i]; ++k) {
               val = SCIPgetSolVal(scip, sol, probdata->PaVars[i][k]);
               assert( SCIPisIntegral(scip, val) );
               if ( val > 0.5 ) {
                  no_parents = FALSE;
                  if (probdata->nParents[i][k] == 0) {
                     fprintf(stream, "%d\t%c\t-\t-\n",i,sex);
                     Done[i] = TRUE;
                     numDone++;
                  } else if (probdata->nParents[i][k] == 1) {
                     if (Done[probdata->ParentSets[i][k][0]] == TRUE) {
                        if (sexconsistent == FALSE || SCIPgetSolVal(scip, sol, probdata->SexVars[probdata->ParentSets[i][k][0]]) < 0.5)
                           fprintf(stream, "%d\t%c\t%d\t-\n",i,sex,probdata->ParentSets[i][k][0]);
                        else
                           fprintf(stream, "%d\t%c\t-\t%d\n",i,sex,probdata->ParentSets[i][k][0]);
                        Done[i] = TRUE;
                        numDone++;
                     }
                  } else {
                     if (Done[probdata->ParentSets[i][k][0]] == TRUE && Done[probdata->ParentSets[i][k][1]] == TRUE) {
                        if (sexconsistent == FALSE || SCIPgetSolVal(scip, sol, probdata->SexVars[probdata->ParentSets[i][k][0]]) < 0.5)
                           fprintf(stream, "%d\t%c\t%d\t%d\n",i,sex,probdata->ParentSets[i][k][0],probdata->ParentSets[i][k][1]);
                        else
                           fprintf(stream, "%d\t%c\t%d\t%d\n",i,sex,probdata->ParentSets[i][k][1],probdata->ParentSets[i][k][0]);
                        Done[i] = TRUE;
                        numDone++;
                     }
                  }
                  break;
               }
            }
            if ( no_parents ) {
               fprintf(stream, "%d\t%c\t-\t-\n",i,sex);
               Done[i] = TRUE;
               numDone++;
            }
         }
      }
   }
   return SCIP_OKAY;
}

/** Prints the solution as a file suitable for plotting using the dot command from graphviz.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param probdata The problem data used by the solution.
 *  @param sol The solution to print.
 *  @param stream Where to print the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionDotFormat(SCIP* scip, SCIP_PROBDATA* probdata, SCIP_SOL* sol, FILE* stream) {
   int i,k,l;
   SCIP_Real val = 0.0;
   SCIP_Real bn_score = 0.0;
   SCIP_Bool no_parents;

   fprintf(stream, "digraph {\n");
   for (i = 0; i < probdata->n; ++i) {
      no_parents = TRUE;
      for (k = 0; k < probdata->nParentSets[i]; ++k) {
         val = SCIPgetSolVal(scip, sol, probdata->PaVars[i][k]);
         assert( SCIPisIntegral(scip, val) );
         if ( val > 0.5 ) {
            bn_score = bn_score + probdata->Scores[i][k];
            for ( l=0; l<probdata->nParents[i][k]-1; ++l )
               fprintf(stream, "   %d -> %d;\n",probdata->ParentSets[i][k][l], i);
            if (probdata->nParents[i][k] > 0)
               fprintf(stream, "   %d -> %d;\n",probdata->ParentSets[i][k][probdata->nParents[i][k]-1], i);
            no_parents = FALSE;
            break;
         }
      }
      if ( no_parents )
         fprintf(stream, "   %d;\n",i);
   }
   fprintf(stream, "}\n");
   return SCIP_OKAY;
}

/** Prints just the objective value of the given solution and the time taken to find the solution.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param probdata The problem data used by the solution.
 *  @param sol The solution to print.
 *  @param stream Where to print the solution.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
static SCIP_RETCODE printSolutionScoreAndTimeFormat(SCIP* scip, SCIP_PROBDATA* probdata, SCIP_SOL* sol, FILE* stream) {
   int i,k;
   SCIP_Real val = 0.0;
   SCIP_Real bn_score = 0.0;
   SCIP_Bool no_parents;

   for (i = 0; i < probdata->n; ++i) {
      no_parents = TRUE;
      for (k = 0; k < probdata->nParentSets[i]; ++k) {
         val = SCIPgetSolVal(scip, sol, probdata->PaVars[i][k]);
         assert( SCIPisIntegral(scip, val) );
         if ( val > 0.5 ) {
            bn_score = bn_score + probdata->Scores[i][k];
            no_parents = FALSE;
            break;
         }
      }
      if ( no_parents )
         bn_score = bn_score + probdata->Scores[i][probdata->nParentSets[i]-1];
   }
   fprintf(stream, "%f\t%f\n",bn_score, SCIPgetSolTime(scip, sol));
   return SCIP_OKAY;
}

/** Prints the solution of the problem after solving.
 *
 *  @param scip The SCIP instance for which to print the solution.
 *  @param filename The filename to output to, "stdout" for stdout or "" for nowhere.
 *  @param format The format in which to print the solution.  Recognised values are dot, pedigree, scoreandtime, legacy and normal.
 *  @return SCIP_OKAY if the solution was printed correctly or an appropriate error message otherwise.
 */
SCIP_RETCODE printSolution(SCIP* scip, char* filename, char* format) {
   SCIP_PROBDATA* probdata;
   SCIP_SOL* sol;
   FILE* file;

   sol = SCIPgetBestSol(scip);

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   assert( probdata->PaVars != NULL );
   assert( probdata->nParentSets != NULL );
   assert( probdata->ParentSets != NULL );
   assert( probdata->nParents != NULL );
   assert( probdata->Scores != NULL );

   // Open the file for writing
   if (strcmp(filename,"") == 0)
      return SCIP_OKAY;
   else if (strcmp(filename,"stdout") == 0)
      file = stdout;
   else {
      file = fopen(filename, "w");
      printf("Writing output to %s\n", filename);
   }
   if ( file == NULL ) {
      SCIPerrorMessage("Could not open file %s for writing.\n", filename);
      return SCIP_WRITEERROR;
   }

   // Print the solution to the file
   if ( sol == NULL )
      printf("No solution found.\n");
   else if (strcmp(format,"dot") == 0)
      printSolutionDotFormat(scip, probdata, sol, file);
   else if (strcmp(format,"pedigree") == 0)
      printSolutionPedigreeFormat(scip, probdata, sol, file);
   else if (strcmp(format,"scoreandtime") == 0)
      printSolutionScoreAndTimeFormat(scip, probdata, sol, file);
   else if (strcmp(format,"legacy") == 0)
      printSolutionLegacyFormat(scip, probdata, sol, file);
   else
      printSolutionBNFormat(scip, probdata, sol, file);

   // Close the file
   if (file != stdout)
      fclose(file);

   return SCIP_OKAY;
}

/** evaluate solution */
SCIP_RETCODE BNevalSolution(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_SOL* sol;
   int i,j,k,l;

   SCIP_Real val = 0.0;
   SCIP_Bool printscipsol;
   SCIP_Bool printmecinfo;

   int** edges;
   int parent1;
   int parent2;
   int ll;

   SCIPgetBoolParam(scip,"gobnilp/printscipsol",&printscipsol);
   SCIPgetBoolParam(scip,"gobnilp/printmecinfo",&printmecinfo);

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   assert( probdata->PaVars != NULL );
   assert( probdata->nParentSets != NULL );
   assert( probdata->ParentSets != NULL );
   assert( probdata->nParents != NULL );
   assert( probdata->Scores != NULL );

   sol = SCIPgetBestSol(scip);

   if ( printscipsol )
      SCIP_CALL( SCIPprintSol(scip,sol,NULL,FALSE) );

   if ( printmecinfo )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &edges, probdata->n) );

      for (i = 0; i < probdata->n; ++i)
	 SCIP_CALL( SCIPallocClearMemoryArray(scip, &(edges[i]), probdata->n) );

      /* skeleton */

      printf("START MEC info\n");

      for (i = 0; i < probdata->n; ++i)
    for (k = 0; k < probdata->nParentSets[i]; ++k)
    {
       val = SCIPgetSolVal(scip, sol, probdata->PaVars[i][k]);
       assert( SCIPisIntegral(scip, val) );
       if ( val > 0.5 )
       {
          for ( l = 0; l < probdata->nParents[i][k]; ++l )
          {
        j = probdata->ParentSets[i][k][l];
        edges[i][j] = TRUE;
        edges[j][i] = TRUE;

          }
       }
    }

      for (i = 0; i < probdata->n; ++i)
    for (j = i+1; j < probdata->n; ++j)
       if ( edges[i][j] )
          printf("%d-%d\n",i,j);


      /* immoralities */

      for (i = 0; i < probdata->n; ++i)
    for (k = 0; k < probdata->nParentSets[i]; ++k)
    {
       val = SCIPgetSolVal(scip, sol, probdata->PaVars[i][k]);
       assert( SCIPisIntegral(scip, val) );
       if ( val > 0.5 )
          for ( l = 0; l < probdata->nParents[i][k]; ++l )
          {
        parent1 = probdata->ParentSets[i][k][l];
        for ( ll = l+1; ll < probdata->nParents[i][k]; ++ll )
        {
           parent2 = probdata->ParentSets[i][k][ll];
           if ( !edges[parent1][parent2] )
         printf("%d->%d<-%d\n",parent1,i,parent2);
        }
          }
    }

      printf("END MEC info\n");
      for ( i = 0 ; i < probdata->n ; ++i)
      {
    SCIPfreeMemoryArray(scip, &(edges[i]));
      }
      SCIPfreeMemoryArray(scip, &edges);
   }
   return SCIP_OKAY;
}


