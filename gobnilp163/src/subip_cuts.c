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
 *  Code using a sub-IP to find "cluster constraint" cutting planes
 */

/*#define SCIP_DEBUG*/
#include "subip_cuts.h"
#include "scip/scipdefplugins.h"

/** Data for subIP for searching for cluster cuts */
typedef struct
{
   SCIP*           subscip;          /**< sub MIP for finding good clusters for cutting planes */
   SCIP_VAR***     family;           /**< subIP variables: family[i][k] = 1 if kth parent set of ith variable is one of those in cluster cut */
   SCIP_VAR**      incluster;        /**< subIP variables: incluster[i] if variable i in cluster */
   SCIP_VAR*       kvar;             /**< lower bound on number of parents to be in cluster for family variable to be set */
   SCIP_CONS***    clausecons;       /**< clausecons[i][k] is the constraint  "if family[i][k]=1 then incluster[i]=1" */
   SCIP_CONS***    overlapcons;      /**< overlapcons[i][k] is the constraint  "if family[i][k]=1 then \sum_{u \in W} >= kvar" */
   SCIP_CONS*      ck_cons;          /**< constraint for lb for cluster size (depends on kvar) */
   SCIP_CONS*      incumbentcons;    /**< (optional) constraint that incumbent must lie on cluster cut */
} DAGCLUSTER_AUXIPDATA;


/** Frees sub-IP data */
static SCIP_RETCODE AuxIPDataFree(
   SCIP* scip,                      /**< (Main) SCIP instance */
   DAGCLUSTER_AUXIPDATA* auxipdata, /**< Data for the sub-IP (to be freed) */
   int n                            /**< Number of BN variables in the acyclicity (dagcluster) constraint */
)
{
   int i;

   if( auxipdata->subscip != NULL )
      SCIP_CALL( SCIPfree(&(auxipdata->subscip)) );

   for( i = 0 ; i < n ; ++i )
   {
      SCIPfreeMemoryArray(scip, &(auxipdata->family[i]));
      SCIPfreeMemoryArray(scip, &(auxipdata->clausecons[i]));
      SCIPfreeMemoryArray(scip, &(auxipdata->overlapcons[i]));
   }
   SCIPfreeMemoryArray(scip, &(auxipdata->family));
   SCIPfreeMemoryArray(scip, &(auxipdata->clausecons));
   SCIPfreeMemoryArray(scip, &(auxipdata->overlapcons));
   SCIPfreeMemoryArray(scip, &(auxipdata->incluster));

   SCIPfreeMemory(scip, &auxipdata);
   auxipdata = NULL;

   return SCIP_OKAY;
}

/** Adds a found cluster cut */
static SCIP_RETCODE AddClusterCut(
   SCIP*           scip,             /**< (Main) SCIP instance */
   const ParentSetData*  psd,        /**< family variable information */
   SCIP_CONSHDLR*  conshdlr,         /**< the constraint handler responsible for adding these cuts (will be 'dagcluster') */
   SCIP_SOL*       sol,              /**< solution to be separated */
   SCIP_Bool*      incluster,        /**< the cluster itself: incluster[i] = TRUE iff i is in the cluster */
   int             cluster_size,     /**< the size of the found cluster */
   int             kval,             /**< kval = 1 for normal cluster cuts, kval = k for a 'k-cluster cut' */
   SCIP_Bool       addtopool,               /**< whether to add the cut to the global cut pool */
   SCIP_Bool       forcecuts,        /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr,   /**< for recording whether the cutting plane is efficacious */
   SCIP_Bool       ci_cut,                        /**< whether we are adding a CI cut */
   SCIP_Bool*      cutoff
)
{
   SCIP_ROW* cut;
   int rhs = cluster_size - kval;
   int i;
   int k;
   int l;
   int n_included;
   int n_excluded;
   SCIP_VAR** included;
   SCIP_VAR** excluded;
   int overlap;
   SCIP_Bool include_in_cut;
   int* parent_set;

   assert(psd != NULL);
   assert(psd->PaVars != NULL);

   /* CI cuts are tighter than cluster cuts */
   if( ci_cut )
      rhs--;

#if SCIP_VERSION >= 300
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cut, conshdlr, "clustercut", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
#else
   SCIP_CALL( SCIPcreateEmptyRow(scip, &cut, "clustercut", -SCIPinfinity(scip), rhs, FALSE, FALSE, TRUE) );
#endif


   for( i = 0 ; i < psd->n ; ++i )
   {
      if( !incluster[i] )
         continue;

      SCIP_CALL( SCIPallocMemoryArray(scip, &included, psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &excluded, psd->nParentSets[i]) );
      n_included = 0;
      n_excluded = 0;


      /* include all parents sets with at least kval parents in cluster */
      for( k = 0;  k < psd->nParentSets[i]; ++k )
      {
         include_in_cut = FALSE;
         overlap = 0;
         parent_set = psd->ParentSets[i][k];
         for( l = 0; l < psd->nParents[i][k]; ++l )
         {
            /* if ( incluster[psd->ParentSets[i][k][l]] ) */
            if( incluster[parent_set[l]] )
            {
               overlap++;
               if( SCIPisGE(scip, overlap, kval) )
               {
                  include_in_cut = TRUE;
                  break;
               }
            }
         }

         if( include_in_cut )
            included[n_included++] = psd->PaVars[i][k];
                  else
            excluded[n_excluded++] = psd->PaVars[i][k];
      }

      /* Use convexity constraint to reduce number of variables in the cut */
      if( n_excluded < psd->nParentSets[i]  / 2 )
      {
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_excluded, excluded, -1.0) );
         rhs--;
      }
      else
         SCIP_CALL( SCIPaddVarsToRowSameCoef(scip, cut, n_included, included, 1.0) );

      SCIPfreeMemoryArray(scip, &included);
      SCIPfreeMemoryArray(scip, &excluded);

      }


   SCIP_CALL( SCIPchgRowRhs(scip, cut, rhs) );
   assert(SCIPisIntegral(scip, rhs));
   SCIPdebugMessage(" -> Cluster-cut <clustercut>: act=%f, rhs=%f, norm=%f, eff=%f, min=%f, max=%f (range=%f)\n",
                    SCIPgetRowLPActivity(scip, cut), SCIProwGetRhs(cut), SCIProwGetNorm(cut),
                    SCIPgetCutEfficacy(scip, NULL, cut),
                    SCIPgetRowMinCoef(scip, cut), SCIPgetRowMaxCoef(scip, cut),
                    SCIPgetRowMaxCoef(scip, cut) / SCIPgetRowMinCoef(scip, cut));
   SCIPdebug(SCIP_CALL( SCIPprintRow(scip, cut, NULL) ));


#if SCIP_VERSION >= 310
   {
      SCIP_CALL( SCIPaddCut(scip, sol, cut, forcecuts, cutoff) );
      if( *cutoff )
      {
         SCIPdebugMessage("Cluster cut led to cutoff\n");
         SCIP_CALL( SCIPreleaseRow(scip, &cut) );
         return SCIP_OKAY;
      }
   }
#else
   SCIP_CALL( SCIPaddCut(scip, sol, cut, forcecuts) );
#endif

   if( addtopool )
      SCIP_CALL( SCIPaddPoolCut(scip, cut) );

   if( SCIPisCutEfficacious(scip, sol, cut) )
      *found_efficacious_ptr = TRUE;

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** main routine for looking for cutting planes

The number of found cutting planes is recorded in *nGen. A positive value indicates that the current solution is infeasible.
 */
extern SCIP_RETCODE IP_findCuts(
   SCIP*           scip,                    /**< SCIP data structure */
   ParentSetData*  psd,                     /**< family variable information */
   SolutionInfo*   solinfo,                 /**< information about the solution to be separated */
   SCIP_SOL*       sol,                     /**< solution to be separated */
   int*            nGen,                    /**< *nGen is number of cutting planes added ( even non-efficacious ones are added ) */
   int             k_lb,                    /**< lowerbound on 'k' values for k-cluster searching, always positive */
   int             k_ub,                    /**< upperbound on 'k' values for k-cluster searching */
   SCIP_CONSHDLR*  conshdlr,                /**< constraint handler */
   SCIP_Bool       addtopool,               /**< whether to add any found cut to the global cut pool */
   SCIP_Bool       forcecuts,               /**< whether to force cuts to be added */
   SCIP_Bool*      found_efficacious_ptr,   /**< to return whether an efficacious cutting plane was found */
   SCIP_Real       limits_time,             /**< limit on how long to spend sub-IP solving */
   SCIP_Real       limits_gap,              /**< limit on size of gap in sub-IP */
   SCIP_Real       limits_absgap,           /**< limit on size of the absolute gap in sub-IP */
   SCIP_Bool       incumbent_cons,          /**< whether to consider only cutting planes on which the incumbent lies */
   int*            must_be_included,        /**< set of nodes which must be included in any found cluster */
   int             n_must_be_included,      /**< size of the set of nodes which must be included in any found cluster */
   int*            must_be_excluded,        /**< set of nodes which must be excluded from any found cluster */
   int             n_must_be_excluded,      /**< size of the set of nodes which must be excluded from any found cluster */
   SCIP_Bool       ci_cut,                  /**< TRUE if we are looking for conditional independence cuts */
   SCIP_Bool*      cutoff                  
)
{

   int i;
   int k;
   int l;
   SCIP_STATUS status;
   int s;
   int nsols;
   SCIP_SOL** subscip_sols;
   SCIP_SOL* subscip_sol;
   SCIP_SOL* best_sol;

   char consname[SCIP_MAXSTRLEN];
   char varname[SCIP_MAXSTRLEN];

   SCIP_Real kval;
   SCIP_Real val;

   SCIP_VAR** clausevars;
   int nvars;

   DAGCLUSTER_AUXIPDATA* auxipdata;          /* data for subscip */

   int cluster_size;
   SCIP_Bool* incluster;

   int ki;

   int* parent_set;

   SCIP_Bool infeasible;
   SCIP_Bool fixed;


   /* check called with sensible 'k' values */

   assert(k_lb > 0);
   assert(k_ub >= k_lb);

   /* check called with sensible must_be_included, must_be_excluded values */

   assert(must_be_included != NULL || n_must_be_included == 0);
   assert(must_be_excluded != NULL || n_must_be_excluded == 0);

   /* check that if a CI cut then k_lb = k_ub = 1 */

   assert(!ci_cut || (k_lb == 1 && k_ub == 1));

   (*nGen) = 0;
   (*found_efficacious_ptr) = FALSE;

   /* allocate temporary memory for building clausal constraints */

   SCIP_CALL( SCIPallocMemoryArray(scip, &clausevars, (psd->n) + 1) );

   /* create and initialise auxiliary IP data structure */

   SCIP_CALL( SCIPallocMemory(scip, &auxipdata) );
   auxipdata->subscip = NULL;
   auxipdata->family = NULL;
   auxipdata->incluster = NULL;
   auxipdata->kvar = NULL;
   auxipdata->clausecons = NULL;
   auxipdata->overlapcons = NULL;
   auxipdata->ck_cons = NULL;

   /* allocate temporary memory for subscip elements */

   SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->family), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->incluster), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->clausecons), psd->n) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->overlapcons), psd->n) );

   for( i = 0 ; i < psd->n ; ++i )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->family[i]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->clausecons[i]), psd->nParentSets[i]) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(auxipdata->overlapcons[i]), psd->nParentSets[i]) );
   }

   /* initialize allocated data structures */

   BMSclearMemoryArray(auxipdata->incluster, psd->n);
   for( i = 0 ; i < psd->n ; ++i )
   {
      BMSclearMemoryArray(auxipdata->family[i], psd->nParentSets[i]);
      BMSclearMemoryArray(auxipdata->clausecons[i], psd->nParentSets[i]);
      BMSclearMemoryArray(auxipdata->overlapcons[i], psd->nParentSets[i]);
   }

   /* create and initialise subscip */

   SCIP_CALL( SCIPcreate(&(auxipdata->subscip)) );


   SCIP_CALL( SCIPincludeDefaultPlugins(auxipdata->subscip) );

   SCIP_CALL( SCIPcreateProb(auxipdata->subscip, "DAG cluster separating MIP", NULL, NULL , NULL , NULL , NULL , NULL , NULL) );

   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetCharParam(auxipdata->subscip, "nodeselection/childsel", 'd') );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "limits/maxsol", 100000) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "limits/maxorigsol", 2000) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "nodeselection/dfs/stdpriority", 536870911) );
   SCIP_CALL( SCIPsetHeuristics(auxipdata->subscip, SCIP_PARAMSETTING_OFF, TRUE) );

   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "lp/solvefreq", 1) );

   SCIP_CALL( SCIPsetRealParam(auxipdata->subscip, "limits/time", limits_time) );
   SCIP_CALL( SCIPsetRealParam(auxipdata->subscip, "limits/gap", limits_gap) );
   SCIP_CALL( SCIPsetRealParam(auxipdata->subscip, "limits/absgap", limits_absgap) );

   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/closecuts/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/cgmip/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/cmir/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/flowcover/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/impliedbounds/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/intobj/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/mcf/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/oddcycle/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/rapidlearning/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/strongcg/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/zerohalf/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "separating/clique/freq", -1) );

   /* experimental */
   /*   SCIP_CALL(  SCIPsetIntParam(auxipdata->subscip, "separating/gomory/freq", -1)  );*/

   /* forbid recursive call of heuristics solving subMIPs */
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "heuristics/rins/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "heuristics/rens/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "heuristics/localbranching/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(auxipdata->subscip, "heuristics/crossover/freq", -1) );

   SCIP_CALL( SCIPsetBoolParam(auxipdata->subscip, "constraints/logicor/negatedclique", FALSE) );


   /* create subscip family variables for each main-problem family variable that is positive in the current solution */
   /* This solution will typically be the solution to the linear relaxation */
   /* use the same name for both */

   for( i = 0 ; i < psd->n ; ++i )
   {

      assert(solinfo->nposvars[i] > -1);
      for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
      {
         k = solinfo->posvars[i][ki];
         assert(k > -1 && k < psd->nParentSets[i]);
         /* essential not to consider the empty parent set */
         if( psd->nParents[i][k] == 0 )
            continue;

         /* val = SCIPgetSolVal(scip, sol, psd->PaVars[i][k]); */
         val = solinfo->lpsolvals[i][k];
         /*if ( !SCIPisPositive(scip,val) )
           printf("foo %f, %d, %d\n",val,i,ki);*/
         assert(SCIPisPositive(scip, val));
         SCIP_CALL( SCIPcreateVar(auxipdata->subscip, &(auxipdata->family[i][k]), SCIPvarGetName(psd->PaVars[i][k]), 0.0, 1.0, val, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(auxipdata->subscip, auxipdata->family[i][k]) );

      }
      /* create variables to identify clusters */
      /* would make it a dummy variable if no non-empty parent sets have positive value, but this oddly causes a slow down (why?) */

      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "incluster#%d", i);
      SCIP_CALL( SCIPcreateVar(auxipdata->subscip, &(auxipdata->incluster[i]), varname, 0.0, 1.0, -1.0, SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(auxipdata->subscip, auxipdata->incluster[i]) );
      SCIP_CALL( SCIPchgVarBranchPriority(auxipdata->subscip, auxipdata->incluster[i], 10) );

   }


   /* create variable for lower bound */
   /* convenient to create it, even if k_lb=k_ub=1 */
   SCIP_CALL( SCIPcreateVar(auxipdata->subscip, &(auxipdata->kvar), "kvar", k_lb, k_ub, 1.0, SCIP_VARTYPE_INTEGER, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(auxipdata->subscip, auxipdata->kvar) );

   /* constraint that incumbent lies on the cutting plane */
   if( incumbent_cons )
   {
      best_sol = SCIPgetBestSol(scip);
      if( best_sol != NULL )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "incumbent");
         SCIP_CALL(SCIPcreateConsLinear(auxipdata->subscip, &(auxipdata->incumbentcons), consname, 0, NULL, NULL,
                                        -1, -1,
                                        TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
         for( i = 0 ; i < psd->n ; ++i )
         {
            SCIP_CALL( SCIPaddCoefLinear(auxipdata->subscip, auxipdata->incumbentcons, auxipdata->incluster[i], -1) );
            for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
            {
               k = solinfo->posvars[i][ki];
               if( solinfo->lpsolvals[i][k] > 0.5 )
               {
                  if( psd->nParents[i][k] > 0 )
                     SCIP_CALL( SCIPaddCoefLinear(auxipdata->subscip, auxipdata->incumbentcons, auxipdata->family[i][k], 1) );
                  break;
               }
            }
         }
         SCIP_CALL( SCIPaddCons(auxipdata->subscip, auxipdata->incumbentcons) );
      }
   }

   /* if family[i][k]=1 then incluster[i]=1 */
   /* (ie in consistent notation with constraints below: If I(W->v)=1 then incluster[v] */
   /*  ~family[i][k]=1 + incluster[i] >= 1 */
   for( i = 0 ; i < psd->n ; ++i )
      for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
      {
         k = solinfo->posvars[i][ki];
         /* essential not to consider the empty parent set */
         if( psd->nParents[i][k] == 0 )
            continue;

         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "clause#%d#%d", i, k);
         SCIP_CALL( SCIPgetNegatedVar(auxipdata->subscip, auxipdata->family[i][k], &(clausevars[0])) );
         clausevars[1] = auxipdata->incluster[i];
         SCIP_CALL(SCIPcreateConsLogicor(auxipdata->subscip, &(auxipdata->clausecons[i][k]), consname, 2, clausevars,
                                         TRUE, TRUE, TRUE, TRUE,
                                         TRUE, /*propagate*/
                                         FALSE, FALSE, FALSE, FALSE, FALSE));
         SCIP_CALL( SCIPaddCons(auxipdata->subscip, auxipdata->clausecons[i][k]) );
      }



   /* k_ub*I(W->v) <= \sum_{u \in W} - k + k_ub */
   /* if I(W->v)=1 this becomes \sum_{u \in W} >= k */
   /* if I(W->v)=0 this becomes vacuous */
   /* just post as a normal linear constraint:
      -inf <= k_ub*I(W->v) - \sum_{u \in W} + k <= k_ub

      note: in the code below 'k' has a different meaning. It indexes
      parent sets
      ' u \in W' is represented by the binary variable incluster[psd->ParentSets[i][k][l]]

      if k_ub == 1 use an equivalent logicor representation

   */

   if( k_ub == 1 )
   {
      for( i = 0 ; i < psd->n ; ++i )
         for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
         {
            k = solinfo->posvars[i][ki];
            /* essential not to consider the empty parent set */
            if( psd->nParents[i][k] == 0 )
               continue;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "overlap#%d#%d", i, k);

            SCIP_CALL( SCIPgetNegatedVar(auxipdata->subscip, auxipdata->family[i][k], &(clausevars[0])) );
            nvars = 1;
            parent_set = psd->ParentSets[i][k];
            for( l = 0; l < psd->nParents[i][k]; ++l )
            {
               /* clausevars[nvars++] = auxipdata->incluster[psd->ParentSets[i][k][l]]; */
               clausevars[nvars++] = auxipdata->incluster[parent_set[l]];
            }
            SCIP_CALL(SCIPcreateConsLogicor(auxipdata->subscip, &(auxipdata->clausecons[i][k]), consname, nvars, clausevars,
                                            TRUE, TRUE, TRUE, TRUE,
                                            TRUE, /* propagate*/
                                            FALSE, FALSE, FALSE, FALSE, FALSE));
            SCIP_CALL( SCIPaddCons(auxipdata->subscip, auxipdata->clausecons[i][k]) );
         }
   }
   else
   {
      for( i = 0 ; i < psd->n ; ++i )
         for( ki = 0; ki < solinfo->nposvars[i]; ++ki )
         {
            k = solinfo->posvars[i][ki];

            /* essential not to consider the empty parent set */
            if( psd->nParents[i][k] == 0 )
               continue;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "overlap#%d#%d", i, k);
            SCIP_CALL(SCIPcreateConsLinear(auxipdata->subscip, &(auxipdata->overlapcons[i][k]), consname, 0, NULL, NULL,
                                           -SCIPinfinity(scip), k_ub,
                                           TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
            SCIP_CALL( SCIPaddCoefLinear(auxipdata->subscip, auxipdata->overlapcons[i][k], auxipdata->family[i][k], k_ub) );
            for( l = 0; l < psd->nParents[i][k]; ++l )
            {
               SCIP_CALL( SCIPaddCoefLinear(auxipdata->subscip, auxipdata->overlapcons[i][k], auxipdata->incluster[psd->ParentSets[i][k][l]], -1) );
            }
            SCIP_CALL( SCIPaddCoefLinear(auxipdata->subscip, auxipdata->overlapcons[i][k], auxipdata->kvar, 1) );
            SCIP_CALL( SCIPaddCons(auxipdata->subscip, auxipdata->overlapcons[i][k]) );
         }

   }


   /* 2 <= |C|-k  <= inf : for the added cut to make sense */
   /* for k=1, this becomes 3 <= |C| <= inf */
   /* can ignore clusters of size two, since arrow variables already deal with them
      and no use for conditional independence cuts */
   if( k_ub == 1 )
      SCIP_CALL(SCIPcreateConsLinear(auxipdata->subscip, &(auxipdata->ck_cons), "ck_constraint", 0, NULL, NULL,
                                     2, SCIPinfinity(scip),
                                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
   else
      SCIP_CALL(SCIPcreateConsLinear(auxipdata->subscip, &(auxipdata->ck_cons), "ck_constraint", 0, NULL, NULL,
                                     1, SCIPinfinity(scip),
                                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

   for( i = 0 ; i < psd->n ; ++i )
   {
      SCIP_CALL( SCIPaddCoefLinear(auxipdata->subscip, auxipdata->ck_cons, auxipdata->incluster[i], 1) );
   }
   if( k_ub != 1 )
      SCIP_CALL( SCIPaddCoefLinear(auxipdata->subscip, auxipdata->ck_cons, auxipdata->kvar, -1) );
   SCIP_CALL( SCIPaddCons(auxipdata->subscip, auxipdata->ck_cons) );

   /* all constraints posted - free temporary memory */

   SCIPfreeMemoryArray(scip, &clausevars);

   /* let I(u) denote u is in the cluster, then
      objective function is \sum_{v,W}I(W->v) - \sum_{u}I(u) + k
      If a feasible solution has a positive objective value,
      then a cutting plane has been found,
      so maximise and rule out non-positive solutions
   */

   SCIP_CALL_ABORT(SCIPsetObjsense(auxipdata->subscip, SCIP_OBJSENSE_MAXIMIZE));

   /* for cluster cuts rule out non-positive solutions using SCIPsetObjlimit
      for CI cuts -1 is enough */
   if( ci_cut )
      SCIP_CALL( SCIPsetObjlimit(auxipdata->subscip, -1) );
   else
      SCIP_CALL( SCIPsetObjlimit(auxipdata->subscip, 0) );

   /*SCIP_CALL( SCIPwriteOrigProblem(auxipdata->subscip,NULL,NULL,FALSE)  );
     SCIP_CALL( SCIPwriteParams(auxipdata->subscip,NULL,FALSE,TRUE)  );
    */

   for( i = 0; i < n_must_be_included; ++i )
   {
      SCIP_CALL( SCIPfixVar(auxipdata->subscip, auxipdata->incluster[must_be_included[i]], 1.0, &infeasible, &fixed) );
      assert(!infeasible && fixed);
   }
   for( i = 0; i < n_must_be_excluded; ++i )
   {
      SCIP_CALL( SCIPfixVar(auxipdata->subscip, auxipdata->incluster[must_be_excluded[i]], 0.0, &infeasible, &fixed) );
      assert(!infeasible && fixed);
   }

   if( ci_cut )
      SCIPdebugMessage("Looking for a conditional independence cluster cut using a subIP...\n");
   else
      SCIPdebugMessage("Looking for a cluster cut using a subIP...\n");


   SCIP_CALL( SCIPsolve(auxipdata->subscip) );


   status = SCIPgetStatus(auxipdata->subscip);
   /*SCIP_CALL( SCIPprintStatus(auxipdata->subscip,NULL) );*/

   if( status == SCIP_STATUS_USERINTERRUPT || status == SCIP_STATUS_INFEASIBLE || status == SCIP_STATUS_INFORUNBD )
   {
      /* /\* print out LP relaxation that could not be separated *\/ */
      /*SCIP_CALL(  SCIPwriteMIP(scip,"foo",TRUE,TRUE)  );
      SCIP_CALL(  SCIPprintSol(scip,NULL,NULL,FALSE)  );
      exit(1);*/

      /*printf("infeasible.\n");
      SCIP_CALL(  SCIPprintSol(scip,NULL,NULL,FALSE)  );*/

      SCIPdebugMessage("could not find a cluster cut.\n");
      SCIP_CALL( AuxIPDataFree(scip, auxipdata, psd->n) );

      return SCIP_OKAY;
   }

   /* if there are feasible solutions but the best has objective value not better that
      0, then we have not found a cutting plane.
      This code snippet from Timo Berthold
   */
   nsols = SCIPgetNSols(auxipdata->subscip);
   if( nsols > 0 && SCIPisFeasLE(auxipdata->subscip, SCIPgetSolOrigObj(auxipdata->subscip, SCIPgetBestSol(auxipdata->subscip)), 0.0) )
   {
      /*printf("obj value too low.\n");
        SCIP_CALL(  SCIPprintSol(scip,NULL,NULL,FALSE)  );
      */
      SCIPdebugMessage("could not find a cluster cut: best cluster objective too low\n");
      SCIP_CALL( AuxIPDataFree(scip, auxipdata, psd->n) );
      return SCIP_OKAY;
   }

   if( status != SCIP_STATUS_SOLLIMIT && status != SCIP_STATUS_GAPLIMIT && status != SCIP_STATUS_OPTIMAL && status != SCIP_STATUS_NODELIMIT  && status != SCIP_STATUS_TIMELIMIT )
   {
      SCIPerrorMessage("Solution of subscip for DAG cluster separation returned with invalid status %d.\n", status);
      SCIP_CALL( AuxIPDataFree(scip, auxipdata, psd->n) );
      return SCIP_ERROR;
   }

   /* To get here a cutting plane must have been found */

   subscip_sols = SCIPgetSols(auxipdata->subscip);

   SCIP_CALL( SCIPallocMemoryArray(scip, &incluster, psd->n) );

   for( s = 0; s <  nsols; ++s )
   {

      subscip_sol = subscip_sols[s];

      /*SCIP_CALL(  SCIPprintSol(auxipdata->subscip,subscip_sol,NULL,FALSE)  );*/

#ifdef SCIP_DEBUG
      if( ci_cut )
         SCIPdebugMessage("found conditional independence cut for this cluster: ");
      else
         SCIPdebugMessage("found cut for this cluster: ");
      for( i = 0 ; i < psd->n ; ++i )
         if( SCIPisPositive(scip, SCIPgetSolVal(auxipdata->subscip, subscip_sol, auxipdata->incluster[i])) )
            SCIPdebugPrintf("%d,", i);
      SCIPdebugPrintf(" %f", SCIPgetSolOrigObj(auxipdata->subscip, subscip_sol));
      SCIPdebugPrintf("\n");
#endif

      kval = SCIPgetSolVal(auxipdata->subscip, subscip_sol, auxipdata->kvar);
      cluster_size = 0;
      for( i = 0 ; i < psd->n ; ++i )
         if( SCIPgetSolVal(auxipdata->subscip, subscip_sol, auxipdata->incluster[i]) > 0.5 )
         {
            incluster[i] = TRUE;
            cluster_size++;
         }
         else
            incluster[i] = FALSE;

      SCIP_CALL( AddClusterCut(scip, psd, conshdlr, sol, incluster, cluster_size, kval, addtopool, forcecuts, found_efficacious_ptr, ci_cut, cutoff) );
      if( *cutoff )
      {
         SCIPdebugMessage("Adding cluster cut led to cutoff\n"); 
         SCIPfreeMemoryArray(scip, &incluster);
         SCIP_CALL( AuxIPDataFree(scip, auxipdata, psd->n) );
         return SCIP_OKAY;
      }
      (*nGen)++;

   }

   SCIPfreeMemoryArray(scip, &incluster);

   SCIPdebugMessage("added %d cluster cuts.\n", *nGen);

   SCIP_CALL( AuxIPDataFree(scip, auxipdata, psd->n) );

   return SCIP_OKAY;
}


