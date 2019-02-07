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
 *  Function declarations for probdata_bn.c
 */

/*
This file was created by editing the file psd_lop.h that comes with the linear ordering example
in SCIP
*/

#ifndef __BN_probdata_bn__
#define __BN_probdata_bn__

#include <scip/scip.h>

extern SCIP_RETCODE BN_readCommandLineArgs(SCIP* scip, int argc, char** argv);
extern SCIP_RETCODE BN_setParamaterDefaults(SCIP* scip);
extern SCIP_RETCODE BN_suppresscols(SCIP* scip);
extern        char* BN_getParameterFile(void);
extern    SCIP_Bool BN_exitBeforeSolving(void);
extern SCIP_RETCODE BN_printScores(SCIP* scip);
extern SCIP_RETCODE BN_printProblem(SCIP* scip, int run);
extern SCIP_RETCODE BN_doIterativePrint(SCIP* scip, int run);
extern SCIP_RETCODE BN_printParameters(SCIP* scip);
extern SCIP_RETCODE BN_printHeader(SCIP* scip);
extern SCIP_RETCODE BN_includePlugins(SCIP* scip);
extern SCIP_RETCODE BN_readProblem(SCIP* scip, const char* filename);
extern SCIP_RETCODE BN_addNonRepetitionConstraint(SCIP* scip, int run);
extern SCIP_RETCODE BN_addParameters(SCIP* scip);
extern          int BN_getNumberOfRepeats(SCIP* scip);
#endif
