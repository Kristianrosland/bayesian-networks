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

import cplex.callbacks as CPX_CB
import sys
import time
import random
import math
import cplex as CPX


import C_prime_cluster_finder
import C_weight_file_reader_bounded_parents
import C_ILP_writer as WR

import networkx as nx
import json
import pickle
import subprocess


from threading import Thread
import threading


class MyBranch(CPX_CB.BranchCallback):

    def __call__(self):
        self.times_called += 1

        nid = self.get_node_ID()
        print '***********************************************************'
        print 'nid nid nid nid nid nid nid nid nid nid nid nid '+str(nid)
        print '***********************************************************'


# Cplex callback that is called at every node
class MyNode(CPX_CB.NodeCallback):

    def __call__(self):
        self.times_called += 1

        # Store the incumbent score and gap
        if self.has_incumbent():
            x = self.get_incumbent_values()
            vars_dict = {}
            for i in range(len(x)):
                if var_names[i].find('y') != -1 or var_names[i].find('z') != -1:
                    vars_dict[var_names[i]] = x[i]
            if feasible_sol_checker(vars_dict) or mode != 1:
                obj_val = self.get_incumbent_objective_value()
                gap = self.get_MIP_relative_gap()
                f = open(res_file, 'a')
                f.write(str(obj_val)+','+str(gap)+',' +
                        str(time.time()-elapsedd)+'\n')
                f.close()


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

# Computes a binomial coefficient


def nCr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

# Writes a DAG to file using gml format


def graph_writer(hyper_edg_set, file_path):
    G = nx.DiGraph()
    for z in hyper_edg_set:
        tmp = z.split('_')
        if len(tmp) >= 3:
            child = tmp[-1]
            for v in tmp[1:-1]:
                G.add_edge(v, child)
        else:
            G.add_node(tmp[1])
    nx.write_gml(G, file_path)


def cutting_plane_utility(I_list):
    # calculates and returns the cutting plane parameters
    ind_list = []
    val_list = []

    vertices = []
    for i in I_list:
        v = i.split('_')[1]
        vertices.append(v)
    edgs_list = comb_list(vertices, 2)
    for e in edgs_list:
        u = e[0]
        v = e[1]
        ind_list.append('y_'+u+'_'+v)
        ind_list.append('y_'+v+'_'+u)
        val_list.append(1)
        val_list.append(1)
    my_rhs = len(I_list)*int(tw)-0.5*(int(tw)+1)*int(tw)

    return ind_list, val_list, my_rhs

# Checks whether the z-graph is acyclic


def z_acyclicity_checker(z_edg_set):
    G = nx.DiGraph()
    for z in z_edg_set:
        tmp = z.split('_')
        if len(tmp) >= 3:
            child = tmp[-1]
            for v in tmp[0:-1]:
                G.add_edge(v, child)
    is_acyclic = nx.is_directed_acyclic_graph(G)
    return is_acyclic

# Returns (some of the) cycles in the z-graph


def z_cycles(z_edg_set):
    G = nx.DiGraph()
    for z in z_edg_set:
        tmp = z.split('_')
        if len(tmp) >= 3:
            child = tmp[-1]
            for v in tmp[0:-1]:
                G.add_edge(v, child)
    cycles = []
    for node in G.nodes():
        Gprime = G.copy()
        Gprime.add_node('X')
        for parent in Gprime.predecessors(node):
            Gprime.add_edge(parent, 'X')
        try:
            path = nx.shortest_path(Gprime, node, 'X')
        except nx.exception.NetworkXNoPath:
            path = []
        if len(path) > 1:
            cycles.append(path[0:-1])
    return cycles

# Checks whether the y-graph is acyclic


def acyclicity_checker(y_edg_set):
    G = nx.DiGraph()
    for y in y_edg_set:
        tmp = y.split('_')
        if len(tmp) >= 3:
            child = tmp[2]
            parent = tmp[1]
            G.add_edge(parent, child)
    is_acyclic = nx.is_directed_acyclic_graph(G)
    return is_acyclic

# Returns (some of the) cycles in the y-graph


def y_cycles(y_edg_set):
    G = nx.DiGraph()
    for y in y_edg_set:
        tmp = y.split('_')
        if len(tmp) >= 3:
            child = tmp[2]
            parent = tmp[1]
            G.add_edge(parent, child)
    cycles = []
    for node in G.nodes():
        Gprime = G.copy()
        Gprime.add_node('X')
        for parent in Gprime.predecessors(node):
            Gprime.add_edge(parent, 'X')
        try:
            path = nx.shortest_path(Gprime, node, 'X')
        except nx.exception.NetworkXNoPath:
            path = []
        if len(path) > 1:
            cycles.append(path[0:-1])
    return cycles

# Returns a constraint that forbids a cycle


def y_cycle_constraint(cycle):
    ind_list = []
    val_list = []
    for i in range(len(cycle)):
        if i < len(cycle) - 1:
            ind_list.append('y_' + cycle[i] + '_' + cycle[i + 1])
            val_list.append(1)
        else:
            ind_list.append('y_' + cycle[i] + '_' + cycle[0])
            val_list.append(1)
    my_rhs = len(cycle) - 1
    return ind_list, val_list, my_rhs


def feasible_sol_checker(vars_dict):
    # gets the currect solution
    # checks if z and y graphs are acyclic
    # returns 1 for acyclic and 0 for cyclic solutions
    y_set = set([])
    z_set = set([])

    for k, v in vars_dict.iteritems():
        if k.find('y') != -1 and v > 0.9:
            y_set.add(k)

    for k, v in vars_dict.iteritems():
        if k.find('z') != -1 and v > 0.9:
            z_set.add(k)

    y_cyclicity = acyclicity_checker(y_set)
    z_cyclicity = z_acyclicity_checker(z_set)

    if y_cyclicity == True and z_cyclicity == True:
        res = True
    else:
        res = False
    return res


def constraint_maker(J_set, vars_dict):
    # gets the set of constrains and make a raw for cplex according to the template in lpex1.py
    #    rows = [[["x1","x2","x3"],[-1.0, 1.0,1.0]],
    #   [["x1","x2","x3"],[ 1.0,-3.0,1.0]]]
    J_list = []
    one_list = []
    s = []
    for J in J_set:
        tmp = J.split('_')
        s.append(tmp[-1])
    for k in vars_dict.keys():
        tmp = k.split('_')
        child = tmp[-1]
        if child in s:
            ok = 1
            for p in tmp[1:-1]:
                if p in s:
                    ok = 0
                    break
            if ok == 1:
                J_list.append(k)
                one_list.append(1)
    rows = [J_list, one_list]
    return rows

# Returns the moralized graph of a directed graph


def moral_graph(directed_graph):
    mg = nx.Graph(directed_graph)
    for child in directed_graph.nodes():
        for p1 in directed_graph.predecessors(child):
            for p2 in directed_graph.predecessors(child):
                if p1 != p2:
                    mg.add_edge(p1, p2)
    return mg


# Minimum fill-in heuristic for an upper bound for tree-width
def min_fillin(input_graph):
    graph = nx.Graph(input_graph)
    chordal_graph = nx.Graph(input_graph)
    elim_order = []
    nodes_remaining = graph.nodes()
    critical_edge = None
    tw_ub = 0
    while len(nodes_remaining) > 0:
        current_min_node = None
        current_min_fillin = len(nodes_remaining)
        for node in nodes_remaining:
            fillin = 0
            for n1 in graph.neighbors(node):
                for n2 in graph.neighbors(node):
                    if int(n1) < int(n2) and not graph.has_edge(n1, n2):
                        fillin += 1
            if fillin < current_min_fillin:
                current_min_node = node
                current_min_fillin = fillin
        neighbors = graph.neighbors(current_min_node)
        for n1 in neighbors:
            for n2 in neighbors:
                if n1 != n2 and not graph.has_edge(n1, n2):
                    graph.add_edge(n1, n2)
                    chordal_graph.add_edge(n1, n2)
            graph.remove_edge(current_min_node, n1)
        elim_order.append(current_min_node)
        nodes_remaining.remove(current_min_node)
        if len(neighbors) > tw_ub:
            tw_ub = len(neighbors)
            nn = input_graph.neighbors(current_min_node)
            found = False
            for node in nn:
                if node in nodes_remaining:
                    if input_graph.has_edge(node, current_min_node):
                        critical_edge = (node, current_min_node)
                    found = True
                    break
            if not found:
                for node in nn:
                    if input_graph.has_edge(node, current_min_node):
                        critical_edge = (node, current_min_node)
                    found = True
                    break
    return tw_ub, elim_order, chordal_graph, critical_edge


# Minimum cardinality heuristic for an upper bound for tree-width
def min_cardinality(input_graph):
    graph = nx.Graph(input_graph)
    chordal_graph = nx.Graph(input_graph)
    elim_order = []
    nodes_remaining = graph.nodes()
    tw_ub = 0
    critical_edge = None
    while len(nodes_remaining) > 0:
        current_min_node = nodes_remaining[0]
        current_min_degree = graph.degree(current_min_node)
        for node in nodes_remaining:
            if graph.degree(node) < current_min_degree:
                current_min_node = node
                current_min_degree = graph.degree(node)
        neighbors = graph.neighbors(current_min_node)
        for n1 in neighbors:
            for n2 in neighbors:
                if n1 != n2 and not graph.has_edge(n1, n2):
                    graph.add_edge(n1, n2)
                    chordal_graph.add_edge(n1, n2)
            graph.remove_edge(current_min_node, n1)
        elim_order.append(current_min_node)
        nodes_remaining.remove(current_min_node)
        if current_min_degree > tw_ub:
            tw_ub = current_min_degree
            nn = input_graph.neighbors(current_min_node)
            found = False
            for node in nn:
                if node in nodes_remaining:
                    if input_graph.has_edge(node, current_min_node):
                        critical_edge = (node, current_min_node)
                    found = True
                    break
            if not found:
                for node in nn:
                    if input_graph.has_edge(node, current_min_node):
                        critical_edge = (node, current_min_node)
                    found = True
                    break

    return tw_ub, elim_order, chordal_graph, critical_edge

# A heuristic for pathwidth upper bound


def min_cardinality_pw(input_graph):
    n_nodes = input_graph.number_of_nodes()
    graph = nx.Graph(input_graph)
    elim_order = []
    path_decomposition = []
    current_bag = set([])
    nodes_remaining = graph.nodes()
    pw_ub = 0
    critical_edge = None
    while len(nodes_remaining) > int(tw) - 1:
        current_min_node = -1
        current_min_degree = 2*n_nodes
        for node in nodes_remaining:
            c = graph.degree(node)
            if node not in current_bag:
                c += len(current_bag)
            if c < current_min_degree:
                current_min_node = node
                current_min_degree = c
        neighbors = graph.neighbors(current_min_node)
        for n1 in neighbors:
            graph.remove_edge(current_min_node, n1)
        elim_order.append(current_min_node)
        current_bag = current_bag.union([current_min_node])
        current_bag = current_bag.union(neighbors)
        path_decomposition.append(current_bag)
        current_pw = len(current_bag) - 1
        if current_min_node in current_bag:
            current_bag = current_bag.difference([current_min_node])
        nodes_remaining.remove(current_min_node)
        if current_pw > int(pw_ub):
            pw_ub = current_pw
            nn = input_graph.neighbors(current_min_node)
            found = False
            for node in nn:
                if node in current_bag:
                    if input_graph.has_edge(node, current_min_node):
                        critical_edge = (node, current_min_node)
                        found = True
                        break
            if not found:
                for node in nn:
                    if input_graph.has_edge(node, current_min_node):
                        critical_edge = (node, current_min_node)
                        found = True
                        break
    for node in nodes_remaining:
        elim_order.append(node)
    final_bag = current_bag.union(nodes_remaining)
    path_decomposition.append(final_bag)
    if len(final_bag) - 1 > pw_ub:
        pw_ub = len(final_bag) - 1
        found = False
        for node1 in final_bag:
            for node2 in final_bag:
                if input_graph.has_edge(node1, node2):
                    critical_edge = (node1, node2)
                    found = True
                    break
            if not found:
                for node1 in final_bag:
                    nn = input_graph.neighbors(node1)
                    for node2 in nn:
                        if input_graph.has_edge(node1, node2):
                            critical_edge = (node1, node2)
                            found = True
                            break
    max_bag = 0
    for b in path_decomposition:
        if len(b) - 1 > max_bag:
            max_bag = len(b) - 1
    if max_bag > pw_ub:
        print "max_bag: " + str(max_bag) + ' pw_ub: ' + str(pw_ub)
        sys.exit()
    return pw_ub, elim_order, path_decomposition, critical_edge


# Own heuristic for tree-width upper bound
def tw_y_heuristic(input_graph, y_dict, root_vars):
    graph1 = nx.Graph(input_graph)
    chordal_graph1 = nx.Graph(input_graph)
    elim_order1 = []
    graph2 = nx.Graph(input_graph)
    chordal_graph2 = nx.Graph(input_graph)
    elim_order2 = []
    tw_ub = 0
    tw_ub1 = 0
    tw_ub2 = 0
    critical_edge1 = None
    critical_edge2 = None

    y_graph = nx.DiGraph()
    for node in graph1.nodes():
        y_graph.add_node(node)
    for k, v in y_dict.iteritems():
        tmp = k.split('_')
        opposite_edge = 'y_' + tmp[2] + '_' + tmp[1]
        if float(v) > float(y_dict[opposite_edge]) and tmp[1] not in root_vars and tmp[2] not in root_vars:
            y_graph.add_edge(tmp[1], tmp[2])

        if nx.is_directed_acyclic_graph(y_graph):
            elim_order1 = nx.topological_sort(y_graph)
        else:
            y_graph2 = nx.DiGraph()
            for node in graph1.nodes():
                y_graph2.add_node(node)
            for k, v in y_dict.iteritems():
                tmp = k.split('_')
                opposite_edge = 'y_' + tmp[2] + '_' + tmp[1]
                if float(v) > float(y_dict[opposite_edge]) and tmp[1] not in root_vars and tmp[2] not in root_vars and float(v) > 0.5:
                    y_graph2.add_edge(tmp[1], tmp[2])
            if nx.is_directed_acyclic_graph(y_graph2):
                elim_order1 = nx.topological_sort(y_graph2)
            else:
                elim_order1 = range(int(len(y_graph.nodes())))

    for node in elim_order1:
        neighbors = graph1.neighbors(node)
        for n1 in neighbors:
            for n2 in neighbors:
                if n1 != n2 and not graph1.has_edge(n1, n2):
                    graph1.add_edge(n1, n2)
                    chordal_graph1.add_edge(n1, n2)
            graph1.remove_edge(node, n1)
        if len(neighbors) > tw_ub1:
            tw_ub1 = len(neighbors)
            nn = input_graph.neighbors(node)
            found = False
            for n in nn:
                if input_graph.has_edge(n, node):
                    critical_edge1 = (n, node)
                found = True
                break

    elim_order2 = elim_order1[:]
    elim_order2.reverse()

    for node in elim_order2:
        neighbors = graph2.neighbors(node)
        for n1 in neighbors:
            for n2 in neighbors:
                if n1 != n2 and not graph2.has_edge(n1, n2):
                    graph2.add_edge(n1, n2)
                    chordal_graph2.add_edge(n1, n2)
            graph2.remove_edge(node, n1)
        if len(neighbors) > tw_ub2:
            tw_ub2 = len(neighbors)
            nn = input_graph.neighbors(node)
            found = False
            for n in nn:
                if input_graph.has_edge(n, node):
                    critical_edge2 = (n, node)
                found = True
                break

    if tw_ub1 < tw_ub2:
        tw_ub = tw_ub1
        elim_order = elim_order1
        chordal_graph = chordal_graph1
        critical_edge = critical_edge1
    else:
        tw_ub = tw_ub2
        elim_order = elim_order2
        chordal_graph = chordal_graph2
        critical_edge = critical_edge2

    return tw_ub, elim_order, chordal_graph, critical_edge


# Transforms a moralized graph to a k-tree
def to_k_tree(graph, elim_order, tw, y_vars):
    new_graph = nx.DiGraph()
    elim_order.reverse()
    root_vars = []
    for v in y_vars:
        if v not in elim_order:
            root_vars.append(v)
    for node in root_vars:
        new_graph.add_node(node)
    for n1 in new_graph.nodes():
        for n2 in new_graph.nodes():
            if int(n1) < int(n2):
                new_graph.add_edge(n1, n2)
    for node in elim_order:
        new_graph.add_node(node)
        for n in graph.neighbors(node):
            if n in new_graph.nodes():
                new_graph.add_edge(node, n)
        counter = 0
        while int(new_graph.degree(node)) < int(tw):
            counter += 1
            for n in new_graph.nodes():
                cur_node = None
                if n != node and n not in new_graph.neighbors(node):
                    cliq = True
                    for m in new_graph.neighbors(node):
                        if not (new_graph.has_edge(n, m) or new_graph.has_edge(m, n)):
                            cliq = False
                            break
                    if cliq == True:
                        cur_node = n
                    if cur_node != None:
                        if cur_node != node and not new_graph.has_edge(node, cur_node):
                            if cur_node not in root_vars or int(len(new_graph.predecessors(cur_node))) < int(tw):
                                new_graph.add_edge(node, cur_node)
                                if int(new_graph.degree(node)) >= int(tw):
                                    break

    elim_order.reverse()
    return new_graph

# Returns an upper bound for tree-width using several heuristics


def treewidth_heuristics_combined(directed_graph, y_dict, y_vars):
    mg = moral_graph(directed_graph)
    tw_ub1, elim_order1, chordal_graph1, critical_edge1 = min_cardinality(mg)
    tw_ub2, elim_order2, chordal_graph2, critical_edge2 = min_fillin(mg)
    tw_ub3, elim_order3, chordal_graph3, critical_edge3 = tw_y_heuristic(mg, y_dict, y_vars)
    critical_edge = None
    if int(tw_ub1) < int(tw_ub2) and int(tw_ub1) < int(tw_ub3):
        tw_ub = tw_ub1
        elim_order = elim_order1
        chordal_graph = chordal_graph1
        critical_edge = critical_edge1
    elif int(tw_ub2) < int(tw_ub3):
        tw_ub = tw_ub2
        elim_order = elim_order2
        chordal_graph = chordal_graph2
        critical_edge = critical_edge2
    else:
        tw_ub = tw_ub3
        elim_order = elim_order3
        chordal_graph = chordal_graph3
        critical_edge = critical_edge3
    return tw_ub, elim_order, chordal_graph, critical_edge

# Removes edges from G so that its tree-width will be at most tw


def make_low_treewidth(G, y_dict, y_vars, tw):
    tw_ub, elim_order, chordal_graph, critical_edge = treewidth_heuristics_combined(G, y_dict, y_vars)
    while int(tw_ub) > int(tw):
        if G.has_edge(critical_edge[0], critical_edge[1]):
            G.remove_edge(critical_edge[0], critical_edge[1])
        elif G.has_edge(critical_edge[1], critical_edge[0]):
            G.remove_edge(critical_edge[1], critical_edge[0])
        else:
            children0 = G.successors(critical_edge[0])
            children1 = G.successors(critical_edge[1])
            for node in children0:
                if node in children1:
                    G.remove_edge(critical_edge[0], node)
                    break
        tw_ub, elim_order, chordal_graph, critical_edge = treewidth_heuristics_combined(G, y_dict, y_vars)
    return tw_ub, elim_order, chordal_graph, critical_edge

# Removes edges from G so that its pathwidth will be at most tw


def make_low_pathwidth(G, tw):
    mg = moral_graph(G)
    pw_ub, elim_order, path_decomposition, critical_edge = min_cardinality_pw(mg)
    while int(pw_ub) > int(tw):
        if G.has_edge(critical_edge[0], critical_edge[1]):
            G.remove_edge(critical_edge[0], critical_edge[1])
        elif G.has_edge(critical_edge[1], critical_edge[0]):
            G.remove_edge(critical_edge[1], critical_edge[0])
        else:
            children0 = G.successors(critical_edge[0])
            children1 = G.successors(critical_edge[1])
            for node in children0:
                if node in children1:
                    G.remove_edge(critical_edge[0], node)
                    break
        mg = moral_graph(G)
        pw_ub, elim_order, path_decomposition, critical_edge = min_cardinality_pw(mg)
    return pw_ub, elim_order, path_decomposition, critical_edge

# Transforms a heuristic solution for treewidth to a form that can be inputted to cplex during the heuristic callback


def make_solution(chordal_graphx, elim_orderx, twx, y_varsx, G, tuple_map, coeffs):
    k_tree = to_k_tree(chordal_graphx, elim_orderx, twx, y_varsx)
    variables = []
    values = []
    root_vars = []
    score = 0
    for v in y_varsx:
        if v not in elim_orderx:
            root_vars.append(v)
    z_variables = construct_z_variables(G, tuple_map)
    for i in range(len(var_names)):
        tmp = var_names[i].split('_')
        if var_names[i].find('z') != -1:
            if var_names[i] in z_variables:
                variables.append(i)
                values.append(1)
                score += coeffs[i]
            else:
                variables.append(i)
                values.append(0)
        elif var_names[i].find('y') != -1:
            if k_tree.has_edge(tmp[1], tmp[2]):
                variables.append(i)
                values.append(1)
            else:
                variables.append(i)
                values.append(0)
    return variables, values, score

# Transforms a heuristic solution for pathwidth to a form that can be inputted to cplex during the heuristic callback


def make_solution_pw(path_decomposition, G, tuple_map, coeffs):
    variables = []
    values = []
    score = 0
    z_variables = construct_z_variables(G, tuple_map)
    for i in range(len(var_names)):
        tmp = var_names[i].split('_')
        if var_names[i].find('z') != -1:
            if var_names[i] in z_variables:
                variables.append(i)
                values.append(1)
                score += coeffs[i]
            else:
                variables.append(i)
                values.append(0)
        elif var_names[i].find('I') != -1:
            if tmp[1] in path_decomposition[int(tmp[2])]:
                variables.append(i)
                values.append(1)
            else:
                variables.append(i)
                values.append(0)
        elif var_names[i].find('J') != -1:
            if tmp[1] in path_decomposition[int(tmp[3])] and tmp[2] in path_decomposition[int(tmp[3])]:
                variables.append(i)
                values.append(1)
            else:
                variables.append(i)
                values.append(0)
    return variables, values, score

# A heuristic for finding an order that hopefully produces a small feedback arc set by Eades, Lin and Smyth


def feedback_arc_set_heuristic(digraph):
    G = nx.DiGraph(digraph)
    s1 = []
    s2 = []
    while G.number_of_nodes() > 0:
        while G.number_of_nodes() > 0 and min(G.out_degree().values()) == 0:
            sinks = []
            for node in G.nodes():
                if int(G.out_degree(node)) == 0:
                    sinks.append(node)
                    G.remove_node(node)
            for node in sinks:
                s2.append(node)
        while G.number_of_nodes() > 0 and min(G.in_degree().values()) == 0:
            sources = []
            for node in G.nodes():
                if int(G.in_degree(node)) == 0:
                    sources.append(node)
                    G.remove_node(node)
            for node in sources:
                s1.append(node)
        if G.number_of_nodes() > 0:
            md = max(G.degree().values())
            md_node = None
            for node in G.nodes():
                if int(G.degree(node)) == int(md):
                    md_node = node
                    break
            if md_node != None:
                s1.append(md_node)
                G.remove_node(md_node)
    s2.reverse()
    s = s1
    for node in s2:
        s.append(node)
    return s

# Returns an acyclic subgraph of digraph


def acyclisize(digraph):
    node_order = feedback_arc_set_heuristic(digraph)
    acyclic_graph = nx.DiGraph()
    predecessors = []
    removed_arcs = 0
    for node in node_order:
        acyclic_graph.add_node(node)
        parents = digraph.in_edges(node)
        for p in parents:
            if p[0] in predecessors:
                acyclic_graph.add_edge(p[0], node)
            else:
                removed_arcs += 1
        predecessors.append(node)
    return acyclic_graph

# Chack whether variables are integral


def integrality_checker(vars_dict):
    for v in vars_dict.values():
        if float(v) > 0 and float(v) < 1:
            return False
    return True

# Returns a list of the z variables names in a digraph


def construct_z_variables(digraph, tuple_map):
    z_variables = []
    for child in digraph.nodes():
        parents = digraph.predecessors(child)
        if tuple_map.get((child, frozenset(parents))) != None:
            t = tuple_map[(child, frozenset(parents))]
            z_variables.append(t[0])
        else:
            s = len(parents)
            best_score = float("-inf")
            best_set = None
            parent_set_found = False
            while parent_set_found == False:
                combs = comb_list(parents, s - 1)
                for comb in combs:
                    par = tuple_map.get((child, frozenset(comb)))
                    if par != None:
                        if float(par[1]) > best_score:
                            best_score = float(par[1])
                            best_set = comb
                            parent_set_found = True
                s = s - 1
            t = tuple_map[(child, frozenset(best_set))]
            z_variables.append(t[0])
    return z_variables


class MyIncumbentCallback(CPX_CB.IncumbentCallback):

    def __call__(self):
        print 'Incumbent'


# Cplex callback for a heuristic for finding feasible solutions
class MyHeuristicCallback(CPX_CB.HeuristicCallback):

    def __call__(self):
        print 'Heuristic callback called'

        self.times_called += 1

        incumbent_score = float("-1e+75")
        if self.has_incumbent():
            incumbent_score = self.get_incumbent_objective_value()

        x = self.get_values()
        coeffs = self.get_objective_coefficients()
        best_score = {}
        best_parents = {}
        mapping = {}
        best_score2 = {}
        best_parents2 = {}
        mapping2 = {}
        tuple_map = {}
        y_vars = set()
        y_dict = {}
        for i in range(len(x)):
            if var_names[i].find('z') != -1:
                tmp = var_names[i].split('_')
                if (tmp[-1] in best_score and float(x[i]) > best_score[tmp[-1]]) or tmp[-1] not in best_score:
                    best_score[tmp[-1]] = float(x[i])
                    best_parents[tmp[-1]] = tmp[1:-1]
                    mapping[tmp[-1]] = var_names[i]
                if float(x[i]) > 0 and ((tmp[-1] in best_score2 and float(coeffs[i]) > best_score2[tmp[-1]]) or tmp[-1] not in best_score2):
                    best_score2[tmp[-1]] = float(coeffs[i])
                    best_parents2[tmp[-1]] = tmp[1:-1]
                    mapping2[tmp[-1]] = var_names[i]
                tuple_map[(tmp[-1], frozenset(tmp[1:-1]))
                          ] = (var_names[i], self.get_objective_coefficients(i))
            elif var_names[i].find('y') != -1:
                tmp = var_names[i].split('_')
                y_vars.add(tmp[1])
                y_vars.add(tmp[2])
                y_dict[var_names[i]] = x[i]

        G = nx.DiGraph()
        for child in best_parents.keys():
            G.add_node(child)
            for parent in best_parents[child]:
                G.add_edge(parent, child)
        if nx.freeze(G) in heuristic_dict:
            G_found = True
        else:
            heuristic_dict[nx.freeze(G)] = []
            G_found = False

        H = nx.DiGraph()
        for child in best_parents2.keys():
            H.add_node(child)
            for parent in best_parents2[child]:
                H.add_edge(parent, child)
        if nx.freeze(H) in heuristic_dict:
            H_found = True
        else:
            heuristic_dict[nx.freeze(H)] = []
            H_found = False

        if (not G_found or not H_found):
            is_acyclic = nx.is_directed_acyclic_graph(G)
            if not is_acyclic:
                G = acyclisize(G)
                is_acyclic = True
            is_acyclic = nx.is_directed_acyclic_graph(H)
            if not is_acyclic:
                H = acyclisize(H)
                is_acyclic = True

        if is_acyclic:
            if mode == 1:
                if int(tw) < int(var_no) - 1 and (not G_found or not H_found):
                    G = nx.DiGraph(G)
                    H = nx.DiGraph(H)
                    tw_ub1, elim_order1, chordal_graph1, critical_edge1 = make_low_treewidth(
                        G, y_dict, y_vars, tw)
                    tw_ub2, elim_order2, chordal_graph2, critical_edge2 = make_low_treewidth(
                        H, y_dict, y_vars, tw)
                    tw_ub = min([tw_ub1, tw_ub2])
                    if int(tw_ub) <= int(tw):
                        variables1, values1, score1 = make_solution(
                            chordal_graph1, elim_order1, tw, y_vars, G, tuple_map, coeffs)
                        variables2, values2, score2 = make_solution(
                            chordal_graph2, elim_order2, tw, y_vars, H, tuple_map, coeffs)
                        if float(score1) >= float(score2):
                            variables = variables1
                            values = values1
                        else:
                            variables = variables2
                            values = values2
                        self.set_solution([variables, values])
                elif int(tw) == int(var_no) - 1:
                    variables = []
                    values = []
                    z_variables = construct_z_variables(G, tuple_map)
                    for i in range(len(var_names)):
                        tmp = var_names[i].split('_')
                        if var_names[i].find('z') != -1:
                            if var_names[i] in z_variables:
                                variables.append(i)
                                values.append(1)
                            else:
                                variables.append(i)
                                values.append(0)
                    self.set_solution([variables, values])
            if mode == 2 and (not G_found or not H_found):
                G = nx.DiGraph(G)
                H = nx.DiGraph(H)
                pw_ub1, elim_order1, path_decomposition1, critical_edge1 = make_low_pathwidth(
                    G, tw)
                pw_ub2, elim_order2, path_decomposition2, critical_edge2 = make_low_pathwidth(
                    H, tw)
                pw_ub = min([pw_ub1, pw_ub2])
                if int(pw_ub) <= int(tw):
                    variables1, values1, score1 = make_solution_pw(
                        path_decomposition1, G, tuple_map, coeffs)
                    variables2, values2, score2 = make_solution_pw(
                        path_decomposition2, H, tuple_map, coeffs)
                    if float(score1) >= float(score2):
                        variables = variables1
                        values = values1
                    else:
                        variables = variables2
                        values = values2
                    self.set_solution([variables, values])
            if self.get_incumbent_objective_value() > incumbent_score:
                print "Improved the incumbent"
                obj_val = self.get_incumbent_objective_value()
                gap = self.get_MIP_relative_gap()
                f = open(res_file, 'a')
                f.write(str(obj_val)+','+str(gap)+',' +
                        str(time.time()-elapsedd)+'\n')
                f.close()

        print "Heuristic callback finished"

# Cplex callback for adding constraints as cutting planes


class MyLazyCutCallback(CPX_CB.LazyConstraintCallback):

    def __call__(self):

        self.times_called += 1
        print 'User Lazy cut opportunity '
        print self.get_objective_value()
        print self.get_MIP_relative_gap()
        if self.has_incumbent():
            print self.get_incumbent_objective_value()
            if mode == 1 or mode == 2:
                obj_val = self.get_incumbent_objective_value()
                gap = self.get_MIP_relative_gap()
                f = open(res_file, 'a')
                f.write(str(obj_val)+','+str(gap)+',' +
                        str(time.time()-elapsedd)+'\n')
                f.close()

        #####################################
        # Cuts for y

        # making the dictionary
        if mode == 1:
            vars_dict = {}
            vars_dict2 = {}
            x = self.get_values()
            for i in range(len(x)):
                if var_names[i].find('y') != -1:
                    vars_dict[var_names[i]] = x[i]
                vars_dict2[var_names[i]] = x[i]
            y_set = set([])
            z_set = set([])

            for k, v in vars_dict.iteritems():
                if k.find('y') != -1 and float(v) > 0.9:
                    y_set.add(k)

            for k, v in vars_dict.iteritems():
                if k.find('z') != -1 and float(v) > 0.9:
                    z_set.add(k)

            counter = 0
            s = ''
            for i in range(len(x)):
                if counter == 0:
                    s = ''
                if float(x[i]) > 0:
                    s += var_names[i] + '=' + str(x[i]) + '\t'
                    if counter == 5:
                        counter = 0
                    else:
                        counter += 1

            if int(tw) < int(var_no) - 1:
                ###########################################

                # getting the cutting plane
                cluster_finder = C_prime_cluster_finder.C_cluster_finder(
                    vars_dict, ip_waiting_time, var_no, tw)
                I_list, sol_status, cplex_elapsed_time, gap = cluster_finder.give_cutting_plane(
                    'y')

                if len(I_list) > 0 and sol_status != 'no solution':
                    for sol in I_list:
                        if len(sol) > 0:
                            self.times_called += 1
                            # logging
                            log_dict[self.times_called] = len(sol)

                            # adding the cutting plane
                            ind_list, val_list, my_rhs = cutting_plane_utility(
                                sol)

                            print '************ Lazy User cut applied for y  ************'

                            lhs = CPX.SparsePair(ind=ind_list, val=val_list)
                            self.add(lhs, 'L', my_rhs)

        ########################################
        # Cussens cuts for z

        vars_dict = {}
        x = self.get_values()
        for i in range(len(x)):
            if var_names[i].find('z') != -1:
                vars_dict[var_names[i]] = x[i]
            # checking

        # getting the cutting plane
        cluster_finder2 = C_prime_cluster_finder.Cussens_cluster_finder(
            vars_dict, ip_waiting_time, var_no)
        I_list, sol_status, cplex_elapsed_time, gap = cluster_finder2.give_cutting_plane()

        if len(I_list) > 0 and sol_status != 'no solution':
            for sol in I_list:
                if len(sol) > 0:
                    self.times_called += 1
                    # logging
                    log_dict[self.times_called] = len(I_list)
                    print '************ Lazy User cut for z applied  ************'

                    rows = constraint_maker(sol, vars_dict)
                    self.add(rows, 'G', 1)

        # Just in case,forbid at least one cycle if it exists
        z_edg_set = []
        for i in range(len(x)):
            if var_names[i].find('z') != -1 and x[i] > 0.9:
                z_edg_set.append(var_names[i])
        cycles = z_cycles(z_edg_set)
        for cycle in cycles:
            r1 = []
            c1 = []
            for node in cycle:
                for z_edg in z_edg_set:
                    tmp = z_edg.split('_')
                    if node == tmp[-1]:
                        r1.append(z_edg)
                        c1.append(1)
            row = [r1, c1]
            self.add(row, 'L', len(r1) - 1)

        print 'Lazy callback completed'


def pipe(scores_path, tree_width, output_path, sub_ip_waiting_time, waiting_time_main, max_parent_no, weights, vars_no, input_mode, debug):
    global var_names
    global var_no
    global tw
    global ip_waiting_time
    global log_dict
    global c
    global sol_params
    global res_file

    global elapsedd
    global mode

    global heuristic_dict  # Network structures that have been seen before

    mode = input_mode  # 1 = tw, 2 = pw, 3 = vc
    heuristic_dict = {}
    elapsedd = time.time()

    score_filename = scores_path.split('/')[-1]
    if mode == 1:
        file_prefix = output_path+'/'+'tw_' + \
            str(tree_width)+'_mp_'+str(max_parent_no)+'_'+score_filename
    elif mode == 2:
        file_prefix = output_path+'/'+'pw_' + \
            str(tree_width)+'_mp_'+str(max_parent_no)+'_'+score_filename
    elif mode == 3:
        file_prefix = output_path+'/'+'vc_' + \
            str(tree_width)+'_mp_'+str(max_parent_no)+'_'+score_filename
    res_file = file_prefix + '_gap_scores.csv'
    f = open(res_file, 'w')
    f.close()
    sol_params = []

    time_lim = waiting_time_main

    log_dict = {}
    start = time.time()
    #global variables
    ip_waiting_time = sub_ip_waiting_time
    tw = tree_width

    var_no = vars_no

    # Write the initial ILP
    cwriter = WR.C_main_ILP_writer(weights, tree_width, var_no, mode)
    c = cwriter.problem_writer()
    if debug:
        c.write(output_path+'/'+'debug.lp')

    c.set_log_stream(sys.stdout)
    c.set_results_stream(sys.stdout)

    # Register callbacks
    lazycut_instance = c.register_callback(MyLazyCutCallback)
    lazycut_instance.times_called = 0

    branch_instance = c.register_callback(MyBranch)
    branch_instance.times_called = 0

    node_instance = c.register_callback(MyNode)
    node_instance.times_called = 0

    if mode == 1 or mode == 2:
        heuristic_instance = c.register_callback(MyHeuristicCallback)
        heuristic_instance.times_called = 0
        

    # Set cplex parameters
    # Due to the cutting planes, some of the presolve methods need to be turned off
    c.parameters.mip.strategy.search.set(
        c.parameters.mip.strategy.search.values.traditional)
    c.parameters.preprocessing.repeatpresolve.set(0)
    c.parameters.preprocessing.relax.set(0)
    c.parameters.preprocessing.presolve.set(0)
    c.parameters.mip.strategy.presolvenode.set(-1)
    c.parameters.preprocessing.reduce.set(0)

    c.parameters.mip.tolerances.mipgap.set(0.0000001)

    c.parameters.timelimit.set(int(waiting_time_main))

    # Force to branch once and awhile
    # c.parameters.mip.limits.cutpasses.set(-1)

    var_names = c.variables.get_names()

    # Solve the ILP
    c.solve()

    solution = c.solution
    objective_value = c.solution.get_objective_value()
    gap = c.solution.MIP.get_mip_relative_gap()
    sol_status = c.solution.status[c.solution.get_status()]
    end = time.time()
    elapsed_time = end-start

    if mode == 1:
        print 'The treewidth bound was ' + str(tw)
    elif mode == 2:
        print 'The pathwidth bound was ' + str(tw)
    elif mode == 3:
        print 'The vertex cover bound was ' + str(tw)
    print 'The objective value is '+str(objective_value)
    if sol_status == 'MIP_optimal':
        print 'The solution is optimal'
    else:
        print 'The gap is ' + str(gap)

    # making the y_graph
    G = nx.DiGraph()
    G_x = nx.DiGraph()

    J_set = set([])
    z_set = set([])
    I_set = set([])
    y_set = set([])
    a_set = set([])
    x = c.solution.get_values(0, c.variables.get_num()-1)
    for j in range(len(x)):
        if c.variables.get_names(j).find('J') != -1 and x[j] > 0.9:
            J = c.variables.get_names(j)
            J_set.add(J)

        if c.variables.get_names(j).find('y') != -1 and x[j] > 0.9:
            y = c.variables.get_names(j)
            y_set.add(y)
            tmp = y.split('_')
            G.add_edge(tmp[1], tmp[2])

        if c.variables.get_names(j).find('z') != -1 and x[j] > 0.9:
            z = c.variables.get_names(j)
            z_set.add(z)
            if mode == 2 or mode == 3:
                tmp = z.split('_')
                if len(tmp) > 1:
                    for i in range(len(tmp) - 2):
                        G.add_edge(tmp[i + 1], tmp[-1])
        if c.variables.get_names(j).find('I') != -1 and x[j] > 0.9:
            I = c.variables.get_names(j)
            I_set.add(I)
        if c.variables.get_names(j).find('a') != -1 and x[j] > 0.9:
            a = c.variables.get_names(j)
            a_set.add(a)
        graph_writer(z_set, file_prefix + '_z.gml')

    if mode == 1:
        print y_set
        if acyclicity_checker(y_set):
            print 'y is acyclic'

    print z_set
    if z_acyclicity_checker(z_set):
        print 'z is acyclic'
    else:
        print 'z is cyclic'

    if mode == 2 or mode == 3:
        t = min(min_fillin(moral_graph(G)), min_cardinality(moral_graph(G)))
        print 'Treewidth of the optimal graph is at most ' + str(t[0])

    if mode == 1:
        s = nx.topological_sort(G)
        nx.write_gml(G, file_prefix + '_y.gml')
        print 'elimination order is'
        print s

    if mode == 2:
        print 'path decomposition is '
        print I_set
        d = {}
        for I in I_set:
            tmp = I.split('_')
            if tmp[2] not in d:
                d[tmp[2]] = [tmp[1]]
            else:
                t = d[tmp[2]]
                t.append(tmp[1])
                d[tmp[2]] = t
        keys = d.keys()
        keys.sort()
        for key in keys:
            print "bag " + str(key)
            print d[key]

    if mode == 3:
        print 'vertex cover is '
        print a_set

    # writing the result file
    learned_score_file = open(file_prefix + '.result', 'w')
    learned_score_file.write(score_filename+'\t'+str(tree_width)+'\t'+str(
        max_parent_no)+'\t'+str(objective_value)+'\t'+str(gap)+'\t'+str(elapsed_time)+'\n')

    learned_score_file.close()

    print "Algorithm finished!"
