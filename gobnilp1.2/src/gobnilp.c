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
This file was created by editing the cmain.c file that comes with the linear ordering example
in SCIP
*/

#include <string.h>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include "probdata_bn.h"
#include "cons_dagcluster.h"
#include "heur_sinks.h"

#include "versiongit.c"

#include <stdio.h>

#define DEFAULT_GOBNILP_PARAMS_FILE "gobnilp.set"

/** read parameters from a file */
static
SCIP_RETCODE readParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< parameter file name, or NULL */
   )
{
   if (filename != NULL)
   {
      if( SCIPfileExists(filename) )
      {
         printf("reading parameter file <%s>\n", filename);
         SCIP_CALL( SCIPreadParams(scip, filename) );
      }
      else
      {
         printf("parameter file <%s> not found.\n", filename);
    return SCIP_NOFILE;
      }
   }
   else
   {
      printf("parameter file empty.\n");
      return SCIP_NOFILE;
   }

   return SCIP_OKAY;
}


/** Adds a boolean parameter to those recognised by SCIP.
 *
 *  This is just a shortcut for SCIPaddBoolParam() with various options set to their most common values.
 *  Use the full function if you need any of the more advanced options.
 *
 *  @param scip The SCIP instance to add the parameter to.
 *  @param name The parameter's name.
 *  @param desc A description of the parameter.
 *  @param value The parameter's initial value.
 *  @return SCIP_OKAY if the operation suceeded.  Otherwise an appropriate error message.
 */
static SCIP_RETCODE addBoolParam(SCIP* scip, const char* name, const char* desc, SCIP_Bool value) {
   return SCIPaddBoolParam(scip, name, desc, NULL, FALSE, value, NULL, NULL);
}

/** Adds an integer parameter to those recognised by SCIP.
 *
 *  This is just a shortcut for SCIPaddIntParam() with various options set to their most common values.
 *  Use the full function if you need any of the more advanced options.  The minimum value of the
 *  parameter is the initial value; the maximum value of the parameter is INT_MAX.
 *
 *  @param scip The SCIP instance to add the parameter to.
 *  @param name The parameter's name.
 *  @param desc A description of the parameter.
 *  @param value The parameter's initial value.
 *  @return SCIP_OKAY if the operation suceeded.  Otherwise an appropriate error message.
 */
static SCIP_RETCODE addIntParam(SCIP* scip, const char* name, const char* desc, int value) {
   return SCIPaddIntParam(scip, name, desc, NULL, FALSE, value, value, INT_MAX, NULL, NULL);
}

/** Adds a string parameter to those recognised by SCIP.
 *
 *  This is just a shortcut for SCIPaddStringParam() with various options set to their most common values.
 *  Use the full function if you need any of the more advanced options.
 *
 *  @param scip The SCIP instance to add the parameter to.
 *  @param name The parameter's name.
 *  @param desc A description of the parameter.
 *  @param value The parameter's initial value.
 *  @return SCIP_OKAY if the operation suceeded.  Otherwise an appropriate error message.
 */
static SCIP_RETCODE addStringParam(SCIP* scip, const char* name, const char* desc, const char* value) {
   return SCIPaddStringParam(scip, name, desc, NULL, FALSE, value, NULL, NULL);
}

/** Adds GOBNILP specific parameters to those recognised by SCIP.
 *
 *  @param scip The SCIP instance to add to the parameters to.
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error coede otherwise.
 */
static SCIP_RETCODE addParameters(SCIP* scip) {
   SCIP_CALL(addBoolParam(scip,
      "gobnilp/noimmoralities",
      "whether to disallow immoralities",
      FALSE
   ));

   SCIP_CALL(addBoolParam(scip,
      "gobnilp/orderedcoveredarcs",
      "whether to only allow a covered arc i<-j if i<j",
      FALSE
   ));

   SCIP_CALL(addBoolParam(scip,
      "gobnilp/implicitfounders",
      "whether to represent empty parent sets implicitly",
      FALSE
   ));

   SCIP_CALL(addBoolParam(scip,
      "gobnilp/printscipsol",
      "whether to (additionally) print BNs in SCIP solution format",
      FALSE
   ));

   SCIP_CALL(addIntParam(scip,
      "gobnilp/nbns",
      "gobnilp to find the 'nbns' best BNs ( in decreasing order of score )",
      1
   ));

   SCIP_CALL(addIntParam(scip,
      "gobnilp/minfounders",
      "minimum number of founders",
      0
   ));

   SCIP_CALL(addIntParam(scip,
      "gobnilp/maxfounders",
      "maximum number of founders (-1 for no upper bound )",
      -1
   ));

   SCIP_CALL(addIntParam(scip,
      "gobnilp/minedges",
      "minimum number of edges",
      0
   ));

   SCIP_CALL(addIntParam(scip,
      "gobnilp/maxedges",
      "maximum number of edges (-1 for no upper bound )",
      -1
   ));

   SCIP_CALL(addBoolParam(scip,
      "gobnilp/printparameters",
      "whether to print parameters not at default values",
      TRUE
   ));

   SCIP_CALL(addBoolParam(scip,
      "gobnilp/printmecinfo",
      "whether to print edges in the undirected skeleton and any immoralities",
      FALSE
   ));

   SCIP_CALL(addStringParam(scip,
      "gobnilp/dagconstraintsfile",
      "file containing constraints on dag structure",
      ""
   ));

   SCIP_CALL(addStringParam(scip,
      "gobnilp/statisticsfile",
      "file for statistics",
      ""
   ));

   SCIP_CALL(addBoolParam(scip,
      "gobnilp/printstatistics",
      "whether to print solving statistics",
      FALSE
   ));

   SCIP_CALL(addBoolParam(scip,
      "gobnilp/printbranchingstatistics",
      "whether to print variable branching statistics",
      FALSE
   ));

   SCIP_CALL(addBoolParam(scip,
      "gobnilp/sexconsistent",
      "whether to enforce sexual consistency in the dag",
      FALSE
   ));

   SCIP_CALL(addStringParam(scip,
      "gobnilp/outputfile/solution",
      "file which solution should be printed to (stdout for standard out, empty string for nowhere)",
      "stdout"
   ));

   SCIP_CALL(addStringParam(scip,
      "gobnilp/outputfile/dot",
      "file which dot output should be printed to (stdout for standard out, empty string for nowhere)",
      ""
   ));

   SCIP_CALL(addStringParam(scip,
      "gobnilp/outputfile/pedigree",
      "file which pedigree output should be printed to (stdout for standard out, empty string for nowhere)",
      ""
   ));

   SCIP_CALL(addStringParam(scip,
      "gobnilp/outputfile/scoreandtime",
      "file which additional score and time data should be printed to (stdout for standard out, empty string for nowhere)",
      ""
   ));

   return SCIP_OKAY;
}



/** Finds the location in a string of the file extension.
 *
 *  @param filename The filename to inspect.
 *  @return The location of the "." marking the beginning of the file extension, or -1 if there is no extension.
 */
static int findExtension(char* filename) {
   int i;
   for (i = strlen(filename); i >= 0; i--)
      if (filename[i] == '.')
         return i;
      else if (filename[i] == '/')
         return -1;
   return -1;
}

/** Constructs a filename for output for a particular iteration of the program.
 *
 *  The filename constructed will be the input filename the iteration number appearing before the file extension
 *  or appended if there is no extension.  If there is only a single network to find, then the iteration number is
 *  not inserted.  For special values ("" and "stdout"), the value returned is the same value as given as input.
 *
 *  @param filename The filename to insert the string in to.
 *  @param nbns The number of Bayesian networks that the program is trying to find.
 *  @param iteration The current iteration of the program.
 *  @return The filename that should be used for output.
 */
static char* createFilename(char* filename, int nbns, int iteration) {
   if (strcmp(filename, "") == 0)
      // Blank string is a special string meaning nowhere
      return filename;
   else if (strcmp(filename, "stdout") == 0)
      // stdout is a special string meaning standard output
      return filename;
   else if (nbns == 1)
      // If only one BN, there is no need to add iteration numbers
      return filename;
   else {
      // Need to add iteration numbers for each BN
      int i;
      int extpos;
      char* ans;
      char insertion[SCIP_MAXSTRLEN];
      sprintf(insertion, "_%d", iteration+1);
      ans = malloc((strlen(filename)+strlen(insertion)+1) * sizeof(char));
      extpos = findExtension(filename);
      if (extpos == -1) {
         for (i = 0; i < (int)strlen(filename); i++)
            ans[i] = filename[i];
         for (i = 0; i < (int)strlen(insertion); i++)
            ans[i+strlen(filename)] = insertion[i];
         ans[strlen(filename)+strlen(insertion)] = '\0';
      } else {
         for (i = 0; i < extpos; i++)
            ans[i] = filename[i];
         for (i = 0; i < (int)strlen(insertion); i++)
            ans[i+extpos] = insertion[i];
         for (i = extpos; i < (int)strlen(filename); i++)
            ans[i+strlen(insertion)] = filename[i];
         ans[strlen(filename)+strlen(insertion)] = '\0';
      }
      return ans;
   }
}


/** main function, which starts the solution of the BN learning problem */
int main(
   int    argc,
   char** argv
   )
{
   SCIP* scip = NULL;

   SCIP_SOL* sol;
   SCIP_PROBDATA* probdata;
   SCIP_Bool no_parents;
   int i,k;
   int n_empty;
   SCIP_Real val;
   int run;

   char opt;

   char consname[SCIP_MAXSTRLEN];

   int nbns;
   const char *g = DEFAULT_GOBNILP_PARAMS_FILE;

   int* chosen;
   SCIP_CONS *cons;
   SCIP_Bool printparameters;
   SCIP_Bool sexconsistent;

   SCIP_Bool printstatistics;
   SCIP_Bool printbranchingstatistics;
   char* statisticsfile;
   FILE* statsfile;

   char* solfile;
   char* dotfile;
   char* pedfile;
   char* satfile;

   /* check paramaters */
   /* if (argc < 2 || argc > 4) */
   /* { */
   /*    printf("usage: %s <BN file> [<parameter file>] [<subscip parameter file>]\n", argv[0]); */
   /*    return 1; */
   /* } */

   /* output version information */
   printf("GOBNILP version %s [GitHash: %s ]\n", GOBNILP_VERSION, GOBNILP_GITHASH);
   printf("Solving the BN structure learning problem using SCIP.\n\n");

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   // Add the custom parameters
   SCIP_CALL( addParameters(scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include DAG cluster constraint handler */
   SCIP_CALL( SCIPincludeConshdlrDagcluster(scip) );



   /* read parameters if requested */
   /* only parameter is name of parameter file  */


   for ( i = 1; i < argc-1; i++ )
   {
      if ( argv[i][0] != '-' )
      {
    printf( "Each optional argument must be preceded by '-'.\n" );
    return 1;
      }
      opt = argv[i][1];
      switch ( opt )
      {
      case 'g':
    g = argv[i]+2;
    break;
      default:
    printf( "Unrecognised optional argument. Can only be g.\n" );
    return 1;
      }
   }

#if SCIP_VERSION == 300
   assert( scip != NULL );
   SCIPprintVersion(scip,NULL);
#else
   SCIPprintVersion(NULL);
#endif
   printf("\n");

   /* include sink heuristic */
   SCIP_CALL( SCIPincludeHeurSinks(scip) );

   SCIP_CALL( readParams(scip, g) );

   SCIPgetIntParam(scip,"gobnilp/nbns", &nbns);
   SCIPgetBoolParam(scip,"gobnilp/printstatistics", &printstatistics);
   SCIPgetBoolParam(scip,"gobnilp/printparameters", &printparameters);
   SCIPgetBoolParam(scip,"gobnilp/printbranchingstatistics", &printbranchingstatistics);
   SCIPgetBoolParam(scip,"gobnilp/sexconsistent", &sexconsistent);
   SCIPgetStringParam(scip,"gobnilp/statisticsfile", &statisticsfile);

   SCIPgetStringParam(scip,"gobnilp/outputfile/solution", &solfile);
   SCIPgetStringParam(scip,"gobnilp/outputfile/dot", &dotfile);
   SCIPgetStringParam(scip,"gobnilp/outputfile/pedigree", &pedfile);
   SCIPgetStringParam(scip,"gobnilp/outputfile/scoreandtime", &satfile);

   if ( printparameters )
   {
      printf("START Parameters not at default value\n");
      SCIP_CALL( SCIPwriteParams(scip,NULL,FALSE,TRUE) );
      printf("END Parameters not at default value\n");
      fflush(stdout);
   }

   if ( strcmp(statisticsfile,"") != 0 )
   {
      statsfile = fopen(statisticsfile, "w");
      if ( statsfile == NULL )
      {
    SCIPerrorMessage("Could not open file %s.\n", statisticsfile);
    return SCIP_NOFILE;
      }
   }
   else
      statsfile = NULL;


   /*SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) ); */

   /* SCIP_CALL( SCIPtransformProb(scip) );*/

   /*SCIP_CALL( BNgivesol(scip) );*/

   /* read problem data */
   SCIP_CALL( BNcreateProb(scip, argv[argc-1]) );

   /* generate BN learning model */
   SCIP_CALL( BNgenerateModel(scip) );

   probdata = SCIPgetProbData(scip);
   SCIP_CALL( SCIPallocMemoryArray(scip, &chosen, probdata->n) );

   for ( run = 0; run < nbns; ++run )
   {

      /* solve the model */
      SCIP_CALL( SCIPsolve(scip) );

      if ( printstatistics )
    SCIP_CALL( SCIPprintStatistics(scip, statsfile) );

      if ( printbranchingstatistics )
    SCIP_CALL( SCIPprintBranchingStatistics(scip, statsfile) );

      SCIP_CALL( printSolution(scip, createFilename(solfile, nbns, run), (char*)"legacy") );
      SCIP_CALL( printSolution(scip, createFilename(dotfile, nbns, run), (char*)"dot") );
      SCIP_CALL( printSolution(scip, createFilename(pedfile, nbns, run), (char*)"pedigree") );
      SCIP_CALL( printSolution(scip, createFilename(satfile, nbns, run), (char*)"scoreandtime") );

      SCIP_CALL( BNevalSolution(scip) );

      sol =  SCIPgetBestSol(scip);

      /* record which BN just found before doing 'free transform' */

      n_empty = 0;
      for (i = 0; i < probdata->n; ++i)
      {
    no_parents = TRUE;
          for (k = 0; k < probdata->nParentSets[i]; ++k)
          {
             val = SCIPgetSolVal(scip, sol, probdata->PaVars[i][k]);
             assert( SCIPisIntegral(scip, val) );
                if ( val > 0.5 )
          {
        chosen[i] = k;
        no_parents = FALSE;
              break;
          }
    }
    if ( no_parents )
    {
       n_empty++;
       chosen[i] = -1;
    }
      }

      SCIP_CALL( SCIPfreeTransform(scip) );
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "ruleout#%d", run);
      /* maybe change this to set covering constraint */
      /* rather than rely on upgrading */
      /* basically the same as CUTOFF_CONSTRAINT(addBinaryCons) in
    cons_countsols.c */
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, -SCIPinfinity(scip),(probdata->n)-1-n_empty,
                         TRUE,
                         TRUE,
                         TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

      for (i = 0; i < probdata->n; ++i)
      {
    if ( chosen[i] == -1 )
       for (k = 0; k < probdata->nParentSets[i]-1; ++k)
                SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][k], -1) );
    else
        SCIP_CALL( SCIPaddCoefLinear(scip, cons, probdata->PaVars[i][chosen[i]], 1) );
      }

      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }


   SCIPfreeMemoryArray(scip, &chosen);

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return 0;
}
