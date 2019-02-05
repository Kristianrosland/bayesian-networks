# -*- coding: iso-8859-1 -*-

# Copyright 2014 Hossein Shahrabi Farahani, Pekka Parviainen
#
# This file is part of TWILP.
#
# TWILP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# TWILP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with TWILP.  If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import math
import random


import cplex
from cplex.exceptions import CplexSolverError
import time

from threading import Thread
import threading
from cplex.callbacks import MIPInfoCallback

import networkx as nx
from numpy.ma.bench import timer
import json


def comb_list(seq, k):
    "returns a list of all k-combinations of the elements of sequence seq"
    n = len(seq)
    if not 0 <= k <= n:
        raise Exception, "0<=k<=len(seq) is not true"
    v = []  # list of combinations

    def f(x, y, a):
        if x == k:
            # we have taken enough elements, reject all remaining elements
            v.append(a)
            return
        if y == n-k:
            # we have rejected enough elements, take all remaining elements
            a.extend(seq[x+y:])
            v.append(a)
            return
        if (x < k):
            # take element seq[x+y]
            h = a+[seq[x+y]]
            f(x+1, y, h)
        if (y < n-k):
            # don't take element seq[x+y]
            f(x, y+1, a)
    f(0, 0, [])
    return v

# Solves the sub-IP for prime cycle cuts (constraint 12)


class C_ILP_solver:

    def ILP_solver_time_limit(self, time_lim, prob):
        I_list = []

        start = time.time()

        prob.parameters.mip.limits.treememory.set(6144)
        prob.parameters.timelimit.set(int(time_lim))

        start_t = time.time()

        prob.solve()

        end_t = time.time()
        cplex_elapsed_time = end_t-start_t

        I_list = []

        sol = prob.solution
        print 'the cluster-finder sol status $$$$$$ is ' + \
            prob.solution.status[prob.solution.get_status()]
        if sol.is_primal_feasible():
            x = prob.solution.get_values(0, prob.variables.get_num()-1)
            gap = prob.solution.MIP.get_mip_relative_gap()
            sol_status = prob.solution.status[prob.solution.get_status()]

        else:
            sol_status = 'no solution'
            gap = 100

        for i in range(prob.solution.pool.get_num()):
            I_set = set()
            x = prob.solution.pool.get_values(i)
            for j in range(prob.variables.get_num()):
                if x[j] > 0.9 and prob.variables.get_names(j).find('I') != -1:
                    I_set.add(prob.variables.get_names(j))
            I_list.append(I_set)

        return I_list, sol_status, cplex_elapsed_time, gap

# Generates a sub-IP for prime cycles (constraint 12)


class C_cluster_finder():
    def __init__(self, vars_dict, waiting_time, vars_no, tw):
        self.vars_dict = vars_dict
        self.waiting_time = waiting_time
        self.vars_no = int(vars_no)
        self.tw = int(tw)

    def ilp_generator(self, var_letter):
        my_lb = []
        my_ub = []
        prob = cplex.Cplex()
        my_obj = []
        my_column_name = []
        rows = []
        my_rhs = []

        n_row_params = []
        n_row_coeff = []

        my_senses = ''
        ctype = ''
        prob.objective.set_sense(prob.objective.sense.maximize)

        edgs = comb_list(range(self.vars_no), 2)

        # adding the first part of the objective function
        for e in edgs:
            u = e[0]
            v = e[1]
            key_1 = var_letter + '_'+str(u)+'_'+str(v)
            key_2 = var_letter + '_'+str(v)+'_'+str(u)

            my_obj.append(self.vars_dict[key_1])
            my_obj.append(self.vars_dict[key_2])
            my_column_name.append("e_"+str(u)+"_"+str(v))
            my_column_name.append("e_"+str(v)+"_"+str(u))
            ctype += 'B'
            ctype += 'B'

            # adding the boundaries

            # condition 16
            # first combination
            row = [["I_"+str(u), "I_"+str(v), "e_" +
                    str(u)+"_"+str(v)], [1, 1, -2]]
            rows.append(row)
            my_senses += 'G'
            my_rhs.append(0)

            my_lb.append(0)
            my_ub.append(1)

            # second combination
            row = [["I_"+str(u), "I_"+str(v), "e_" +
                    str(v)+"_"+str(u)], [1, 1, -2]]
            rows.append(row)
            my_senses += 'G'
            my_rhs.append(0)
            my_lb.append(0)
            my_ub.append(1)

            # condition 17
            # first combination
            row2 = [["I_"+str(u), "I_"+str(v), "e_" +
                     str(u)+"_"+str(v)], [1, 1, -1]]
            my_rhs.append(1)
            my_senses += 'L'
            rows.append(row2)

            # second combination
            row2 = [["I_"+str(u), "I_"+str(v), "e_" +
                     str(v)+"_"+str(u)], [1, 1, -1]]
            my_rhs.append(1)
            my_senses += 'L'
            rows.append(row2)

            # adding efficient ILP boundary PART I
            n_row_params.append("e_"+str(u)+"_"+str(v))
            n_row_coeff.append(self.vars_dict[key_1])

            n_row_params.append("e_"+str(v)+"_"+str(u))
            n_row_coeff.append(self.vars_dict[key_2])

        # other part of the objective function
        r1 = []
        r2 = []
        for i in range(self.vars_no):
            my_obj.append(-self.tw)
            my_column_name.append('I_'+str(i))
            ctype += 'B'
            my_lb.append(0)
            my_ub.append(1)

            # my condition
            r1.append('I_'+str(i))
            r2.append(1)

            # efficient ilp condition
            n_row_coeff.append(-self.tw)
            n_row_params.append('I_'+str(i))

        # adding condition 19
        rows.append([r1, r2])
        my_rhs.append(self.tw+2)
        my_senses += 'G'

        # adding the cutting condition
        # condition 18
        rows.append([n_row_params, n_row_coeff])
        my_rhs.append((-0.5*(self.tw+1)*self.tw)+1e-5)
        my_senses += 'G'

        prob.variables.add(obj=my_obj,  types=ctype, names=my_column_name)
        prob.linear_constraints.add(
            lin_expr=rows, senses=my_senses, rhs=my_rhs)
        return prob

    def give_cutting_plane(self, var_letter):
        prob = self.ilp_generator(var_letter)
        csolver = C_ILP_solver()
        I_set, sol_status, cplex_elapsed_time, gap = csolver.ILP_solver_time_limit(
            float(self.waiting_time), prob)
        return I_set, sol_status, cplex_elapsed_time, gap


#########################################################################################################################
# ILP solver for Cussens cuts
class Cussens_ILP_solver:

    def __init__(self):
        pass

    def ILP_solver_time_limit(self, time_lim, prob):
        J_set = set([])

        start = time.time()

        prob.parameters.mip.limits.treememory.set(6144)
        prob.parameters.timelimit.set(int(time_lim))

        start_t = time.time()

        prob.solve()

        end_t = time.time()
        cplex_elapsed_time = end_t-start_t
        J = []

        sol = prob.solution
        print 'the cluster-finder sol status is ' + \
            prob.solution.status[prob.solution.get_status()]
        if sol.is_primal_feasible():
            x = prob.solution.get_values(0, prob.variables.get_num()-1)
            gap = prob.solution.MIP.get_mip_relative_gap()
            sol_status = prob.solution.status[prob.solution.get_status()]
        else:
            sol_status = 'no solution'
            gap = 100

        for i in range(prob.solution.pool.get_num()):
            J_set = set()
            x = prob.solution.pool.get_values(i)
            for j in range(prob.variables.get_num()):
                if x[j] > 0.9 and prob.variables.get_names(j).find('z') != -1:
                    J_set.add(prob.variables.get_names(j))
            J.append(J_set)

        return J, sol_status, cplex_elapsed_time, gap


#####################
# Cluster cuts by Cussens (constraint 11)
class Cussens_cluster_finder:

    def __init__(self, vars_dict, waiting_time, vars_no):
        self.vars_dict = vars_dict
        self.waiting_time = waiting_time
        self.vars_no = vars_no

    def cnd_2_generator(self, vars_dict):
        rows = []
        for k, v in vars_dict.iteritems():
            tmp = k.split('_')
            if len(tmp) > 2 and float(v) > 0:
                r1 = []
                r2 = []
                r1.append(k)
                r2.append(-1)
                child = tmp[-1]
                r1.append('z_' + str(int(child)))
                r2.append(1)
                parents = tmp[1:-1]
                for p in parents:
                    r1.append('z_' + str(int(p)))
                    r2.append(-1)
                rows.append([r1, r2])
        return rows

    def cnd_3_generator(self, vars_dict):
        rows = []
        for k, v in vars_dict.iteritems():
            tmp = k.split('_')
            if len(tmp) > 2 and float(v) > 0:
                child = tmp[-1]
                r1 = [k, 'z_' + str(int(child))]
                r2 = [1, -1]
                rows.append([r1, r2])
        return rows

    def cnd_4_generator(self, vars_dict):
        rows = []
        for k, v in vars_dict.iteritems():
            tmp = k.split('_')
            if len(tmp) > 2 and float(v) > 0:
                child = tmp[-1]
                parents = tmp[1:-1]
                for p in parents:
                    r1 = [k, 'z_' + str(int(p))]
                    r2 = [1, 1]
                    rows.append([r1, r2])
        return rows

    def cnd_5_generator(self, vars_no):
        r1 = []
        r2 = []
        for v in range(vars_no):
            r1.append('z_' + str(v))
            r2.append(1)
        return [r1, r2]

    def ILP_writer(self, vars_dict, var_no):
        my_lb = []
        my_ub = []
        prob = cplex.Cplex()
        my_obj = []
        my_column_name = []
        rows = []
        my_rhs = []

        my_senses = ''
        ctype = ''
        prob.objective.set_sense(prob.objective.sense.minimize)

        for k, v in vars_dict.iteritems():
            tmp = k.split('_')
            if float(v) > 0:
                my_obj.append(v)
                my_column_name.append(k)
                ctype += 'B'
            elif len(tmp) == 2:
                my_obj.append(0)
                my_column_name.append(k)
                ctype += 'B'

        rows.append([my_column_name, my_obj])
        my_rhs.append(0.99999)
        my_senses += 'L'

        cnd_2 = self.cnd_2_generator(vars_dict)
        for r in cnd_2:
            rows.append(r)
            my_rhs.append(0)
            my_senses += 'L'

        cnd_3 = self.cnd_3_generator(vars_dict)
        for r in cnd_3:
            rows.append(r)
            my_rhs.append(0)
            my_senses += 'L'

        cnd_4 = self.cnd_4_generator(vars_dict)
        for r in cnd_4:
            rows.append(r)
            my_rhs.append(1)
            my_senses += 'L'

        cnd_5 = self.cnd_5_generator(var_no)
        rows.append(cnd_5)
        my_rhs.append(2)
        my_senses += 'G'

        prob.variables.add(obj=my_obj,  types=ctype, names=my_column_name)
        prob.linear_constraints.add(
            lin_expr=rows, senses=my_senses, rhs=my_rhs)
        return prob

    def give_cutting_plane(self):
        prob = self.ILP_writer(self.vars_dict, self.vars_no)
        csolver = Cussens_ILP_solver()
        z_set, sol_status, time, gap = csolver.ILP_solver_time_limit(
            float(self.waiting_time), prob)
        return z_set, sol_status, time, gap
