import pandas as pd # pylint: disable=E0401
import networkx as nx
import itertools
import argparse
import sys
from pgmpy.models import BayesianModel
from pgmpy.estimators import MaximumLikelihoodEstimator, BayesianEstimator
from pgmpy.inference import VariableElimination
from pgmpy.factors.discrete import TabularCPD
from bifWriter import write_model

parser = argparse.ArgumentParser()
parser.add_argument("-t", help="training data file")
parser.add_argument("-g", help="graph file")
parser.add_argument("-o", help="output location")
args = parser.parse_args()

if args.t is None or args.g is None or args.o is None:
    print "Training data (-t), gml file (-g) and output location (-o) must be provided. Exiting parameter estimation \n"
    sys.exit()

input_data = args.t
input_gml = args.g
output_file = args.o

G = nx.read_gml(input_gml)
input_file = open(input_data)

num_vars = int(input_file.readline())
num_values = map(int, input_file.readline().split(' '))
num_datapoints = int(input_file.readline())

datapoints = []
for line in input_file:
    datapoints.append(map(int, line.split(' ')))

dataframe = pd.DataFrame(datapoints, columns=map(str, range(num_vars)))

model = BayesianModel(G.edges())
model.fit(dataframe, estimator=BayesianEstimator)

write_model(output_file, model)