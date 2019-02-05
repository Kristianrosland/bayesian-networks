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

# This class reads the weights from a file formatted according to GOBNILP package standards
# sample output
#{'10,31': '-98.0461862476', '13,2,17': '-58.0735242873', '14,32,33,12': '-25.6585567939', '12,32,3,17': '-76.2888395164', '31,28,13': '-21.776313716'}


class C_weight_reader():
    def __init__(self, weightfile_path, max_parents):
        self.weightfile_path = weightfile_path
        # 0 = unlimitied number of parents
        self.max_parent_number = int(max_parents)

    def max_parent_finder(self):
        # this fucntion parses the input and find the maximum number of parents
        infile = open(self.weightfile_path)
        l = infile.readlines()
        max_par_num = 1
        for i in range(len(l)):
            if l[i][0] == '-':
                tmp = l[i].strip().split(' ')
                parent_num = int(tmp[1])
                if parent_num > max_par_num:
                    max_par_num = parent_num
        infile.close()
        return max_par_num

    def z_param_maker(self, child_vertex, parents_ver_list):
        # takes the child and parents and makes a comma separated list
        s = ''
        for p in parents_ver_list:
            s += p+','
        s += child_vertex
        return s

    def file_reader(self):
        weight_dict = {}
        infile = open(self.weightfile_path)
        l = infile.readlines()
        self.var_number = int(l[0].strip())
        current_line_no = 1

        while current_line_no < len(l):
            tmp = l[current_line_no].strip().split(' ')
            current_child = tmp[0]
            start_loop = current_line_no
            end_loop = current_line_no+int(tmp[1])
            for line_number in range(start_loop, end_loop+1):
                tmp2 = l[line_number].strip().split(' ')

                weight = tmp2[0]
                parents_list = tmp2[2:]
                if self.max_parent_number == 0 or int(tmp2[1]) <= self.max_parent_number:
                    z_string = self.z_param_maker(current_child, parents_list)
                    if float(weight) < 0:
                        weight_dict[z_string] = weight
                    if float(weight) > 0:
                        weight_dict[z_string] = '+'+weight
                current_line_no += 1
        infile.close()
        return weight_dict

    def give_weight(self):
        weight_dict = self.file_reader()
        return weight_dict

    def give_var_number(self):
        return self.var_number
