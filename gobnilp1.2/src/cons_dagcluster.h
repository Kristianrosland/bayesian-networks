/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*   GOBNILP Copyright (C) 2012 James Cussens 	       	       	       	         */
/*    		    	       						 */
/*   This program is free software; you can redistribute it and/or       */
/*   modify it under the terms of the GNU General Public License as	 */
/*   published by the Free Software Foundation; either version 3 of the	 */
/*   License, or (at your option) any later version.   	       	    	 */
/*   	      	     	  	      	    				 */
/*   This program is distributed in the hope that it will be useful,	 */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of	 */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU	 */
/*   General Public License for more details.	  	       	   	 */
/*    	      	     	     	      					 */ 
/*   You should have received a copy of the GNU General Public License	 */
/*   along with this program; if not, see   	 	 		 */
/*   <http://www.gnu.org/licenses>.   					 */ 
/*    									 */
/*   Additional permission under GNU GPL version 3 section 7		 */
/*    		 	    	      	  	    	    		 */
/*   If you modify this Program, or any covered work, by linking or	 */
/*   combining it with SCIP (or a modified version of that library),	 */
/*   containing parts covered by the terms of the ZIB Academic License,  */
/*   the licensors of this Program grant you additional permission to	 */
/*   convey the resulting work.    	      		 	    	 */
/*   	    		  						 */ 
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* 
This file was created by editing the constraint handler template header file
in SCIP
*/ 

#ifndef __SCIP_CONS_DAGCLUSTER_H__
#define __SCIP_CONS_DAGCLUSTER_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for dagcluster constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeConshdlrDagcluster(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a dagcluster constraint */
extern
SCIP_RETCODE SCIPcreateConsDagcluster(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int n,                   /** number of elements */
   int* nParentSets,        /** nParentSets[i] is the number of  parent sets for variable i*/
   int** nParents,          /** nParents[i][k] is the number of  parents in the kth parent set for variable i*/
   int*** ParentSets,       /** ParentSets[i][k][l] is the lth parent in the kth parent set of ith variable **/
   SCIP_VAR*** PaVars,      /** PaVars[i][k] = 1 if kth parent set of ith variable is selected */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which 
                                              *   are seperated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   );

#ifdef __cplusplus
}
#endif

#endif
