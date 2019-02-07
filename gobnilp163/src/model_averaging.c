/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *   GOBNILP Copyright (C) 2012-2015 James Cussens, Mark Bartlett        *
 *                                                                       *
 *   This program is free software; you can redistribute it and/or       *
 *   modify it under the terms of the GNU General Public License as      *
 *   published by the Free Software Foundation; either version 3 of the  *
 *   License, or (at your option) any later version.                     *
 *                                                                       *
 *   This program is distributed in the hope that it will be useful,     *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    *
 *   General Public License for more details.                            *
 *                                                                       *
 *   You should have received a copy of the GNU General Public License   *
 *   along with this program; if not, see                                *
 *   <http://www.gnu.org/licenses>.                                      *
 *                                                                       *
 *   Additional permission under GNU GPL version 3 section 7             *
 *                                                                       *
 *   If you modify this Program, or any covered work, by linking or      *
 *   combining it with SCIP (or a modified version of that library),     *
 *   containing parts covered by the terms of the ZIB Academic License,  *
 *   the licensors of this Program grant you additional permission to    *
 *   convey the resulting work.                                          *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/** @file
 *  Implements the functions needed to perform model averaging over the n best Bayesian networks.
 */

#include "model_averaging.h"
#include "parent_set_data.h"

/** The number of variables that this modelling averaging is performed on. */
static int num_vars;
/** The variables which the model averaging applies to. */
static SCIP_VAR** vars;
/** The average scores of each of the variables in the program.
 *
 *  average_scores[i] is the likelihood weighted average score of vars[i].
 */
static SCIP_Real* average_scores;
/** The sum of the likelihoods of all networks found so far.
 *
 *  This is used to normalise the scores in @link average_scores @endlink .
 */
static SCIP_Real total_score = 0.0;
/** The number of seconds spent finding all of the solutions included in the average. */
static SCIP_Real total_time = 0.0;
/** Find the index of @link vars @endlink and @link average_scores @endlink
 *  relating to a variable.
 *
 *  @param var The variable of interest.
 *  @return The index of the var in @link vars @endlink.  If the variable is not
 *  included in the averaging, -1 is returned.
 */
static int indexOf(SCIP_VAR* var)
{
   int i;
   for( i = 0; i < num_vars; i++ )
      if( vars[i] == var )
         return i;
   return -1;
}

/** Whether the objective function is logarithmic.*/
static SCIP_Bool is_log_score = TRUE;
/** Whether the next solution is the first solution found. */
static SCIP_Bool is_first_score = TRUE;
/** The score of the first solution found. */
static SCIP_Real first_score = 0.0;
/** The normalising constant for logarithmic scores. */
static SCIP_Real normalising_constant = 0.0;

/** Adds parameters for controlling the model averaging.
 *  @param scip The SCIP instance to which the parameter is to be added.
 *  @return SCIP_OKAY if the parameters were added successfully or an error code otherwise.
 */
SCIP_RETCODE MA_addAveragingParameters(SCIP* scip)
{
   SCIP_CALL(SCIPaddBoolParam(scip, "gobnilp/logaverage", "whether model averaging should assume logarithmic objective function",
                              &is_log_score, FALSE, TRUE, NULL, NULL));
   return SCIP_OKAY;
}

/** Allocates memory for the data structures used for model averaging.
 *
 *  @param scip The SCIP instance on which the model averaging will be performed.
 *  @return SCIP_OKAY if memory allocation was successful or an appropriate error
*   message otherwise.
 */
SCIP_RETCODE MA_createAverageDataStructure(SCIP* scip)
{
   int i, j;
   int num_all_vars;
   SCIP_VAR** all_vars;

   num_vars = SCIPgetNBinVars(scip);
   num_all_vars = SCIPgetNVars(scip);
   all_vars = SCIPgetVars(scip);

   SCIP_CALL( SCIPallocMemoryArray(scip, &vars, num_vars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &average_scores, num_vars) );

   /* Copy the binary variables */
   j = 0;
   for( i = 0; i < num_all_vars; i++ )
      if( SCIPvarIsBinary(all_vars[i]) )
      {
         vars[j] = all_vars[i];
         j++;
      }
   assert(j == num_vars);

   /* Set the averages to initial values of 0 */
   for( i = 0; i < num_vars; i++ )
      average_scores[i] = 0.0;

   return SCIP_OKAY;
}
/** Frees memory used for the data structures used for model averaging.
 *
 *  @param scip The SCIP instance on which the model averaging was performed.
 *  @return SCIP_OKAY if memory deallocation was successful or an appropriate error
 *  message otherwise.
 */
SCIP_RETCODE MA_destroyAverageDataStructure(SCIP* scip)
{
   SCIPfreeMemoryArray(scip, &vars);
   SCIPfreeMemoryArray(scip, &average_scores);
   return SCIP_OKAY;
}

/** Updates the average scores based on a newly found solution.
 *
 *  @param scip The SCIP instance to which the solution belongs.
 *  @param sol The new solution to incorporate in to the averages.
 *  @return SCIP_OKAY if the operation succeeded or an appropriate error message otherwise.
 */
SCIP_RETCODE MA_updateAverageDataStructure(SCIP* scip, SCIP_SOL* sol)
{
   int i;
   SCIP_Real overall_score = 0;
   SCIP_Real exp_score = 0;

   /* Calculate the total solution value */
   for( i = 0; i < num_vars; i++ )
      overall_score += SCIPgetSolVal(scip, sol, vars[i]) * SCIPvarGetObj(vars[i]);
   if( is_first_score )
   {
      first_score = overall_score;
      is_first_score = FALSE;
   }
   exp_score = exp(overall_score - first_score);

   /* Update the averages */
   if( is_log_score )
   {
      for( i = 0; i < num_vars; i++ )
         if( SCIPgetSolVal(scip, sol, vars[i]) > 0.5 )
            average_scores[i] += exp_score;
   }
   else
   {
      for( i = 0; i < num_vars; i++ )
         if( SCIPgetSolVal(scip, sol, vars[i]) > 0.5 )
            average_scores[i] += overall_score;
   }

   /* Update the totals */
   total_score += overall_score;
   total_time += SCIPgetSolvingTime(scip);
   normalising_constant += exp_score;

   return SCIP_OKAY;
}

/** Returns the model average value of a given variable.
 *
 *  If the variable is not part of this model averaging, the value -1 will be returned.
 *
 *  @param variable The variable to get the model averaging score for.
 *  @return The model average score of the variable.
 */
SCIP_Real MA_getAverageValue(SCIP_VAR* variable)
{
   int index = indexOf(variable);
   if( index == -1 )
      return -1;
   else if( is_log_score )
      return average_scores[index] / normalising_constant;
   else
      return average_scores[index] / total_score;
}
/** Returns the total time spent solving for all the solutions included
 *  in the model average.
 *
 *  @return The number of seconds spend on solving all of the solutions
 *  that have been used to make up the average.
 */
SCIP_Real MA_getTotalAveragesTime(void)
{
   return total_time;
}
/** Returns the total likelihood of all the solutions included
 *  in the model average.
 *
 *  @return The sum of the likelihoods of each of the solutions included.
 */
SCIP_Real MA_getTotalAveragesScore(void)
{
   return total_score;
}
