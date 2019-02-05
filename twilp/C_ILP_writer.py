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

import math
import C_weight_file_reader_bounded_parents
import cplex

# For writing the ILP


class C_main_ILP_writer():

    def __init__(self, weight_dict, tw, var_number, mode):
        self.weight_dict = weight_dict
        self.tw = int(tw)
        self.var_number = int(var_number)
        self.bags = self.var_number-self.tw+1
        self.mode = mode

    def comb_list(self, seq, k):
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

    def z_maker(self, s):
        # converts 'a,b,c' to 'z_a_b_c'
        t = s.replace(',', '_')
        output = 'z'+'_'+t
        return output

    ################################################################################################
    # Generates the objective function                                                             #
    ################################################################################################
    def objective_function_generator(self):
        my_obj = []
        my_column_name = []
        ctype = ''
        # objective=''
        for edg, weight in self.weight_dict.iteritems():
            my_obj.append(float(weight))
            my_column_name.append(self.z_maker(edg))
            ctype += 'B'
            #objective+=weight+' '+self.z_maker(edg)+' '
        if self.mode == 1:
            for i in range(self.var_number+self.tw):
                for j in range(self.var_number+self.tw):
                    if i != j:
                        my_obj.append(0)
                        my_column_name.append('y_' + str(i) + '_' + str(j))
                        ctype += 'B'
        elif self.mode == 2:
            for i in range(self.var_number):
                for j in range(self.var_number):
                    for l in range(self.bags):
                        if i < j:
                            my_obj.append(0)
                            my_column_name.append(
                                'J_' + str(i) + '_' + str(j) + '_' + str(l))
                            ctype += 'B'
                        if i == j:
                            my_obj.append(0)
                            my_column_name.append('I_' + str(i) + '_' + str(l))
                            ctype += 'B'
        elif self.mode == 3:
            for i in range(self.var_number+self.tw):
                for j in range(self.var_number+self.tw):
                    if i < j:
                        my_obj.append(0)
                        my_column_name.append('y_' + str(i) + '_' + str(j))
                        ctype += 'B'
                    if i == j:
                        my_obj.append(0)
                        my_column_name.append('a_' + str(i))
                        ctype += 'B'
        return my_obj, my_column_name, ctype

    ################################################################################################
    # Constraints:										   #
    # (Constraint 11 and 12 will be added later as cutting planes)				   #
    ################################################################################################

    def cond_1(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        eqs_dict = {}
        eqs_coeff = {}
        # making the dictionary contains equations for each child vertex
        for i in range(self.var_number):
            eqs_dict[str(i)] = []
            eqs_coeff[str(i)] = []

        for z in self.weight_dict.keys():
            child = z.split(',')[-1]
            z_ijk = self.z_maker(z)
            eqs_dict[child].append(z_ijk)
            eqs_coeff[child].append(1)

        # adding =1 to the left side of the equations
        for k, v in eqs_dict.iteritems():
            rows1.append([eqs_dict[k], eqs_coeff[k]])
            my_senses1 += 'E'
            my_rhs1.append(1)
        return rows1, my_rhs1, my_senses1

    #####################################################################################################################
    # Treewidth constraints                                                                                             #
    #####################################################################################################################

    def cond_2(self):
        # for the non-root vertices
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            r1 = []
            c1 = []
            for j in range(int(self.var_number)+int(self.tw)):
                if i != j:
                    r1.append('y_'+str(i)+'_'+str(j))
                    c1.append(1)
            rows1.append([r1, c1])
            my_senses1 += 'E'
            my_rhs1.append(self.tw)
        return rows1, my_rhs1, my_senses1

    def cond_3(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        r1 = []
        c1 = []
        # edges between real vertices
        for i in range(self.var_number):
            for j in range(self.var_number):
                if i != j:
                    r1.append('y_'+str(i)+'_'+str(j))
                    c1.append(1)
        rows1.append([r1, c1])
        edg_number = self.var_number*self.tw-0.5*self.tw*(self.tw+1)
        my_senses1 += 'E'
        my_rhs1.append(edg_number)
        return rows1, my_rhs1, my_senses1

    def cond_4(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            for j in range(self.var_number):
                if i != j:
                    r1 = []
                    c1 = []
                    for k in self.weight_dict.keys():
                        tmp = k.split(',')
                        parents_list = tmp[0:-1]
                        child = tmp[-1]
                        if str(i) in parents_list and j == int(child):
                            z_ijk = self.z_maker(k)
                            r1.append(z_ijk)
                            c1.append(1)
                    if len(r1) > 0:
                        r1.append('y_' + str(i) + '_' + str(j))
                        r1.append('y_' + str(j) + '_' + str(i))
                        c1.append(-1)
                        c1.append(-1)
                        rows1.append([r1, c1])
                        my_senses1 += 'L'
                        my_rhs1.append(0)
        return rows1, my_rhs1, my_senses1

    def cond_5(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for k in self.weight_dict.keys():
            tmp = k.split(',')
            parents_list = tmp[0:-1]
            z_ijk = self.z_maker(k)
            if len(parents_list) >= 2:
                parents_pairs = self.comb_list(parents_list, 2)
                for pair in parents_pairs:
                    rows1.append(
                        [[z_ijk, 'y_'+pair[0]+'_'+pair[1], 'y_'+pair[1]+'_'+pair[0]], [1, -1, -1]])
                    my_senses1 += 'L'
                    my_rhs1.append(0)
        return rows1, my_rhs1, my_senses1

    def cond_6(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            pairs = self.comb_list(range(self.var_number+self.tw), 2)
            # removing outgoing edges from the roots
            for pair in pairs:
                j = pair[0]
                k = pair[1]
                if str(i) != str(j) and str(i) != str(k):
                    rows1.append([['y_'+str(i)+'_'+str(j),  'y_'+str(i)+'_'+str(k),
                                   'y_'+str(j)+'_'+str(k), 'y_'+str(k)+'_'+str(j)], [1, 1, -1, -1]])
                    my_senses1 += 'L'
                    my_rhs1.append(1)
        return rows1, my_rhs1, my_senses1

    def cond_7(self):
        # impose that no edge leaves the roots
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for r in range(self.var_number, self.var_number+self.tw):
            for v in range(self.var_number):
                rows1.append([['y_'+str(r)+'_'+str(v)], [1]])
                my_senses1 += 'E'
                my_rhs1.append(0)
        return rows1, my_rhs1, my_senses1

    def cond_8(self):
        # and root parent connectivity
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        # root parents connectivity
        for i in range(self.var_number, self.var_number+self.tw):
            for j in range(self.var_number, self.var_number+self.tw):
                if i < j:
                    rows1.append([['y_'+str(i)+'_'+str(j)], [1]])
                    my_senses1 += 'E'
                    my_rhs1.append(1)
        return rows1, my_rhs1, my_senses1

    def cond_9(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number+self.tw):
            for j in range(self.var_number+self.tw):
                if i != j:
                    rows1.append(
                        [['y_'+str(i)+'_'+str(j), 'y_'+str(j)+'_'+str(i)], [1, 1]])
                    my_senses1 += 'L'
                    my_rhs1.append(1)
        return rows1, my_rhs1, my_senses1

    def cond_10(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        triples = self.comb_list(range(self.var_number+self.tw), 3)
        for tr in triples:
            i = tr[0]
            j = tr[1]
            k = tr[2]
            rows1.append([['y_'+str(i)+'_'+str(j), 'y_'+str(j) +
                           '_'+str(k), 'y_'+str(k)+'_'+str(i)], [1, 1, 1]])
            rows1.append([['y_'+str(i)+'_'+str(k), 'y_'+str(k) +
                           '_'+str(j), 'y_'+str(j)+'_'+str(i)], [1, 1, 1]])
            my_senses1 += 'LL'
            my_rhs1.append(2)
            my_rhs1.append(2)
        return rows1, my_rhs1, my_senses1

    ###########################################################################################################
    # Pathwidth conditions                                                                                    #
    ###########################################################################################################
    def cond_11(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            for j in range(self.var_number):
                if i < j:
                    r1 = []
                    c1 = []
                    for k in self.weight_dict.keys():
                        tmp = k.split(',')
                        parents_list = tmp[0:-1]
                        child = tmp[-1]
                        if (str(i) in parents_list and j == int(child)) or (str(j) in parents_list and i == int(child)):
                            z_ijk = self.z_maker(k)
                            r1.append(z_ijk)
                            c1.append(1)
                    if len(r1) > 0:
                        for l in range(self.bags):
                            r1.append('J_' + str(i) + '_' +
                                      str(j) + '_' + str(l))
                            c1.append(-1)
                        my_senses1 += 'L'
                        my_rhs1.append(0)
                        rows1.append([r1, c1])
        return rows1, my_rhs1, my_senses1

    def cond_12(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for k in self.weight_dict.keys():
            tmp = k.split(',')
            parents_list = tmp[0:-1]
            z_ijk = self.z_maker(k)
            if len(parents_list) >= 2:
                parents_pairs = self.comb_list(parents_list, 2)
                for pair in parents_pairs:
                    r1 = []
                    c1 = []
                    if int(pair[0]) < int(pair[1]):
                        r1.append(z_ijk)
                        c1.append(1)
                        for l in range(self.bags):
                            r1.append('J_'+pair[0]+'_'+pair[1]+'_' + str(l))
                            c1.append(-1)
                    else:
                        r1.append(z_ijk)
                        c1.append(1)
                        for l in range(self.bags):
                            r1.append('J_'+pair[1]+'_'+pair[0]+'_' + str(l))
                            c1.append(-1)
                    my_senses1 += 'L'
                    my_rhs1.append(0)
                    rows1.append([r1, c1])
        return rows1, my_rhs1, my_senses1

    def cond_13(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            for j in range(self.var_number):
                for l in range(self.bags):
                    if i < j:
                        rows1.append(
                            [['J_' + str(i) + '_' + str(j) + '_' + str(l), 'I_' + str(i) + '_' + str(l)], [1, -1]])
                        rows1.append(
                            [['J_' + str(i) + '_' + str(j) + '_' + str(l), 'I_' + str(j) + '_' + str(l)], [1, -1]])
                        my_senses1 += 'LL'
                        my_rhs1.append(0)
                        my_rhs1.append(0)
        return rows1, my_rhs1, my_senses1

    def cond_14(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            for j in range(self.var_number):
                for l in range(self.bags):
                    if i < j:
                        rows1.append([['J_' + str(i) + '_' + str(j) + '_' + str(l), 'I_' + str(
                            i) + '_' + str(l), 'I_' + str(j) + '_' + str(l)], [2, -1, -1]])
                        my_senses1 += 'L'
                        my_rhs1.append(0)
        return rows1, my_rhs1, my_senses1

    def cond_15(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            for j in range(self.var_number):
                for l in range(self.bags):
                    if i < j:
                        rows1.append([['I_' + str(i) + '_' + str(l), 'I_' + str(j) + '_' + str(
                            l),  'J_' + str(i) + '_' + str(j) + '_' + str(l)], [1, 1, -1]])
                        my_senses1 += 'L'
                        my_rhs1.append(1)
        return rows1, my_rhs1, my_senses1

    def cond_16(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for l in range(self.bags):
            r1 = []
            c1 = []
            for i in range(self.var_number):
                r1.append('I_' + str(i) + '_' + str(l))
                c1.append(1)
            rows1.append([r1, c1])
            my_senses1 += 'L'
            my_rhs1.append(self.tw + 1)
        return rows1, my_rhs1, my_senses1

    def cond_17(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            r1 = []
            c1 = []
            for l in range(self.bags):
                r1.append('I_' + str(i) + '_' + str(l))
                c1.append(1)
            rows1.append([r1, c1])
            my_senses1 += 'G'
            my_rhs1.append(1)
        return rows1, my_rhs1, my_senses1

    def cond_18(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            for l in range(self.bags - 2):
                for ll in range(l + 2, self.bags):
                    for x in range(l + 1, ll):
                        rows1.append([['I_' + str(i) + '_' + str(l), 'I_' + str(i) +
                                       '_' + str(ll), 'I_' + str(i) + '_' + str(x)], [1, 1, -1]])
                        my_senses1 += 'L'
                        my_rhs1.append(1)
        return rows1, my_rhs1, my_senses1

      ################################################################################################################################
      # Vertex cover constraints                                                                                                     #
      ################################################################################################################################

    def cond_21(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            for j in range(self.var_number):
                if i < j:
                    r1 = []
                    c1 = []
                    for k in self.weight_dict.keys():
                        tmp = k.split(',')
                        parents_list = tmp[0:-1]
                        child = tmp[-1]
                        if (str(i) in parents_list and j == int(child)) or (str(j) in parents_list and i == int(child)):
                            z_ijk = self.z_maker(k)
                            r1.append(z_ijk)
                            c1.append(1)
                    if len(r1) > 0:
                        r1.append('y_' + str(i) + '_' + str(j))
                        c1.append(-1)
                        rows1.append([r1, c1])
                        my_senses1 += 'L'
                        my_rhs1.append(0)
        return rows1, my_rhs1, my_senses1

    def cond_22(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for k in self.weight_dict.keys():
            tmp = k.split(',')
            parents_list = tmp[0:-1]
            z_ijk = self.z_maker(k)
            if len(parents_list) >= 2:
                parents_pairs = self.comb_list(parents_list, 2)
                for pair in parents_pairs:
                    if int(pair[0]) < int(pair[1]):
                        rows1.append(
                            [[z_ijk, 'y_'+pair[0]+'_'+pair[1]], [1, -1]])
                    else:
                        rows1.append(
                            [[z_ijk, 'y_'+pair[1]+'_'+pair[0]], [1, -1]])
                    my_senses1 += 'L'
                    my_rhs1.append(0)
        return rows1, my_rhs1, my_senses1

    def cond_23(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        for i in range(self.var_number):
            for j in range(self.var_number):
                if i < j:
                    rows1.append(
                        [['y_' + str(i) + '_' + str(j), 'a_' + str(i), 'a_' + str(j)], [1, -1, -1]])
                    my_senses1 += 'L'
                    my_rhs1.append(0)
        return rows1, my_rhs1, my_senses1

    def cond_24(self):
        rows1 = []
        my_senses1 = ''
        my_rhs1 = []
        r1 = []
        c1 = []
        for i in range(self.var_number):
            r1.append('a_' + str(i))
            c1.append(1)
        rows1.append([r1, c1])
        my_senses1 += 'L'
        my_rhs1.append(self.tw)
        return rows1, my_rhs1, my_senses1

    #################################################################################################################################
    # Set packing inequalities a la Cussens                                                                                         #
    # Not necessary for correct results, may speed up computations though                                                           #
    #################################################################################################################################

    def cond_set_packing(self, set_size):
        eq_list = []
        sets = self.comb_list(range(self.var_number), set_size)
        for s in sets:
            eq = ''
            for i in s:
                for k in self.weight_dict.keys():
                    tmp = k.split(',')
                    parents_list = tmp[0:-1]
                    child = tmp[-1]
                    if int(i) == int(child):
                        subsets = 0
                        for p in parents_list:
                            if int(p) in s:
                                subsets = subsets + 1
                        if subsets >= set_size - 1:
                            if len(eq) == 0:
                                eq = self.z_maker(k)
                            else:
                                eq = eq + ' + ' + self.z_maker(k)
            if len(eq) > 0:
                eq = eq + ' <= 1'
                eq_list.append(eq)
        return eq_list

    #################################################################################################################################
    # Writes the (relaxed) problem as ILP                                                                                           #
    #################################################################################################################################

    def problem_writer(self):
        prob = cplex.Cplex()
        my_lb = []
        my_ub = []
        rows = []
        my_rhs = []
        my_senses = ''
        prob.objective.set_sense(prob.objective.sense.maximize)

        my_obj, my_column_name, ctype = self.objective_function_generator()

        rows1, my_rhs1, my_senses1 = self.cond_1()
        for i in range(len(rows1)):
            rows.append(rows1[i])
            my_rhs.append(my_rhs1[i])
        my_senses += my_senses1

        # Treewidth
        if self.mode == 1 and (self.tw > 2 and self.tw < self.var_number):

            rows1, my_rhs1, my_senses1 = self.cond_2()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_3()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_4()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_5()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_6()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_7()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_8()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_9()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_10()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

        # Pathwidth
        elif self.mode == 2:

            rows1, my_rhs1, my_senses1 = self.cond_11()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_12()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_13()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_14()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_15()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_16()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_17()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_18()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

        # Vertex cover
        elif self.mode == 3:

            rows1, my_rhs1, my_senses1 = self.cond_21()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_22()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_23()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

            rows1, my_rhs1, my_senses1 = self.cond_24()
            for i in range(len(rows1)):
                rows.append(rows1[i])
                my_rhs.append(my_rhs1[i])
            my_senses += my_senses1

        prob.variables.add(obj=my_obj,  types=ctype, names=my_column_name)
        prob.linear_constraints.add(
            lin_expr=rows, senses=my_senses, rhs=my_rhs)
        return prob
