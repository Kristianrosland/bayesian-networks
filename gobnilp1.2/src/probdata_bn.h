/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   GOBNILP Copyright (C) 2012 James Cussens                            */
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
This file was created by editing the file probdata_lop.h that comes with the linear ordering example
in SCIP
*/

#ifndef __BN_PROBDATA_BN__
#define __BN_PROBDATA_BN__

#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#ifdef __cplusplus
extern "C" {
#endif

struct SCIP_ProbData
{
   int n;                   /** number of elements */
   int* nParentSets;        /** nParentSets[i] is the number of  parent sets for variable i*/
   int** nParents;          /** nParents[i][k] is the number of  parents in the kth parent set for variable i*/
   SCIP_Real** Scores;      /** Scores[i][k] is score of kth parent set of ith variable**/
   int*** ParentSets;       /** ParentSets[i][k][l] is the lth parent in the kth parent set of ith variable **/
   SCIP_VAR*** PaVars;      /** PaVars[i][k] = 1 if kth parent set of ith variable is selected */
   SCIP_VAR** SexVars;      /** SexVars[i] = 1 if the ith variable is female or = 0 if it is male */
};


extern
SCIP_RETCODE BNcreateProb(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< name of file to read */
   );

extern
SCIP_RETCODE BNgenerateModel(
   SCIP*                 scip               /**< SCIP data structure */
   );

extern
SCIP_RETCODE BNevalSolution(
   SCIP*                 scip
   );


extern
SCIP_RETCODE BNgivesol(
   SCIP*                 scip
   );

extern
SCIP_RETCODE printSolution(SCIP* scip, char* filename, char* format);

#ifdef __cplusplus
}
#endif

#endif
