import sys
import argparse
import pandas as pd #pylint: disable=import-error
import networkx as nx
import numpy as np
from pgmpy.estimators import BaseEstimator
from pgmpy.readwrite import BIFReader
from pgmpy.models import BayesianModel
from pgmpy.estimators import MaximumLikelihoodEstimator
from pgmpy.factors.discrete import TabularCPD
from pgmpy_utils import import_and_translate_data
from pgmpy.utils.mathext import cartesian

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="input .bif file")
parser.add_argument("-d", help="input .dat file")
parser.add_argument("-tr", help="input translation .txt file")
args = parser.parse_args()


data = pd.DataFrame(data={
    'X': ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',     '1','1','1','1','1','1',       '0',        '1', '1', '1', '1', '1', '1','1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',       '0','0','0',       '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',      '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '1'],
    'Y': ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',     '0','0','0','0','0','0',       '1',        '1', '1', '1', '1', '1', '1','1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',       '0','0','0',       '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',      '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    # 'W': ['0', '0', '0', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',     '1','1','1','1','0','0',       '1',        '0', '0', '0', '0', '0', '0','0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',       '1','1','1',       '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '0', '0', '0', '0',      '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1'],
    'Z': ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',     '0','0','0','0','0','0',       '0',        '0', '0', '0', '0', '0', '0','0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0',       '1','1','1',       '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1',      '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1']})


# G = nx.DiGraph()
# G.add_edges_from([('X', 'Z'), ('Y', 'Z')])
# model = BayesianModel(G.edges())

reader = BIFReader(args.f)
model = reader.get_model()
model.check_model()


data = import_and_translate_data(args.d, args.tr)
model.fit(data, estimator=MaximumLikelihoodEstimator)
estimator = BaseEstimator(data)

print model.nodes()
print data
if model.nodes() != None:
    sys.exit()

# Creates a new counts-dataframe without variable 'v'
def remove_variable_from_counts(child, v, cs):
    parents = filter(lambda x: x != v and x != child, cs.index.names)
    return estimator.state_counts(child, parents=parents).unstack()

# Turns an array into a chain of .loc[]'s
def locs(dataframe, arr):
    if (len(arr) == 0):
        return dataframe
    else:
        return locs(dataframe.loc[arr[0]], arr[1:])


# Creates a new cpd with variable 'v' and counts from 'cs'
def make_cpd(v, cs):
    if cs.index.names[-1] != v:
        print "ERROR: Last variable of counts-dataframe is not {} \n\n{}".format(v, cs)
        sys.exit()

    child = cs.index.names[-1]
    child_states = cs.index.levels[-1]

    parents = cs.index.names[:-1]
    parent_states = cs.index.levels[:-1]
    parent_cardinalities = map(len, parent_states)

    cartesian_product = cartesian(parent_states)
    
    values = []
    for _ in range(len(child_states)):
        values.append([])
    
    for cart in cartesian_product:
        conditional_counts = locs(cs, cart)
        sum_conditional_counts = reduce(lambda x,y: x+y, conditional_counts)
        for i in range(len(child_states)):
            probability = float(conditional_counts[i])/sum_conditional_counts
            values[i].append(probability)
    
    return TabularCPD(child, len(child_states), values, evidence=parents, evidence_card=parent_cardinalities)

def calculate_likelihood(cpd, cs):
    values = cpd.values.reshape(1, np.prod(cpd.cardinality))[0]
    counts = cs.values
    
    if (len(values) != len(counts)):
        print "Error in calculate likelihood: List lengths are not equal"
        sys.exit()

    chunk_size = len(values) / cpd.cardinality[0]
    chunks = [values[i:i + chunk_size] for i in xrange(0, len(values), chunk_size)]

    likelihood = 1
    for entry in range(len(chunks[0])):
        for i in range(cpd.cardinality[0]):
            a = chunks[i][entry]
            b = counts[entry*2 + i]
            likelihood *= a**b


    return likelihood

#counts_without_Y = remove_variable_from_counts('Z', 'Y', counts)
#cpd_without_Y = make_cpd('Z', counts_without_Y)

#counts_without_X = remove_variable_from_counts('Z', 'X', counts)
#cpd_without_X = make_cpd('Z', counts_without_X)

# Evaluates the edges coming into v
def evaluate_edges_into(v):
    cpd = model.get_cpds(v)
    parents = cpd.variables[1:]
    counts = estimator.state_counts(v, parents=parents).unstack()

    likelihood = calculate_likelihood(cpd, counts)
    edge_dependency_strengths = []

    for p in parents:
            counts_without_p = remove_variable_from_counts(v, p, counts)
            cpd_without_p = make_cpd(v, counts_without_p)
            likelihood_without_p = calculate_likelihood(cpd_without_p, counts_without_p)
            ratio = likelihood/likelihood_without_p ## TODO: Change this ratio to log of likelihood ratio
            edge_dependency_strengths.append((p, v, ratio))

    return edge_dependency_strengths

## Try to remove all the edges
for node in model.nodes():
    print "Dependency strength for arcs coming into {}: {}".format(node, evaluate_edges_into(node))


