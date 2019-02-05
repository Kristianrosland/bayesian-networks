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


import prime_pipe
import C_weight_file_reader_bounded_parents
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="input file")
parser.add_argument("-o", help="output directory")
parser.add_argument(
    "-t", help="tree-width bound (0 = unlimited tree-width)", type=int)
parser.add_argument(
    "-p", help="maximum number of parents (0 = unlimited number of parents)", type=int)
parser.add_argument("-r", help="maximum running time (in seconds)", type=int)
parser.add_argument(
    "-s", help="maximum running time for an sub-IP (in seconds)", type=int)
parser.add_argument(
    "-m", help="mode: 1 = treewidth (default), 2 = pathwidth, 3 = vertex cover", type=int)
parser.add_argument(
    "-d", help="debugging mode (0 = debugging off, 1= debugging on)", type=int)
args = parser.parse_args()

if args.f is None:
    print "No input file defined"
    print "Exiting ..."
    sys.exit()
else:
    scores_path = args.f

# Setting arguments
if args.o is None:
    output_path = os.getcwd()
else:
    output_path = args.o

if args.t is None:
    tree_width = int(var_no) - 1
else:
    tree_width = args.t

if args.r is None:
    waiting_time_main = 86400
else:
    waiting_time_main = args.r

if args.s is None:
    waiting_time_sub_IP = 3600
else:
    waiting_time_sub_IP = args.s

if args.p is None:
    max_parent_no = 0
else:
    max_parent_no = args.p

if args.m is None:
    mode = 1
elif args.m > 0 and args.m <= 3:
    mode = args.m
else:
    print 'Invalid mode'
    sys.exit()

if args.d is None:
    debug = False
elif args.d == 0:
    debug = False
else:
    debug = True


print "Reading the input ..."

if max_parent_no == 0 and tree_width == 0:
    pp = 0
elif max_parent_no == 0 and tree_width > 0:
    pp = tree_width
else:
    pp = min(max_parent_no, tree_width)

cwr = C_weight_file_reader_bounded_parents.C_weight_reader(scores_path, pp)
weights = cwr.give_weight()
var_no = cwr.give_var_number()

print 'The data contains ' + \
    str(var_no) + ' variables and in total ' + \
    str(len(weights)) + ' possible parent sets'


# Running the algorithm
obj_val = prime_pipe.pipe(scores_path, tree_width, output_path, waiting_time_sub_IP,
                          waiting_time_main, max_parent_no, weights, var_no, mode, debug)
