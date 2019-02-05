import sys
import argparse
import random
from pgmpy.inference import VariableElimination
from pgmpy.readwrite import BIFReader

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="input .bif file")
args = parser.parse_args()

if args.f is None:
    print "Input file (-f) or output location (-o) not provided \n"
    sys.exit()

reader = BIFReader(args.f)
model = reader.get_model()
model.check_model()

inference = VariableElimination(model)

# RANDOM QUERY
'''
for _ in range(1):
    query_variable = str(random.randint(0, model.number_of_nodes()-1))
    evidence_variable = query_variable
    while (evidence_variable == query_variable):
        evidence_variable = str(random.randint(0, model.number_of_nodes()-1))
    
    evidence_value = random.randint(0, model.get_cardinality(evidence_variable)-1)
    query = inference.query(variables=[query_variable], evidence={evidence_variable: evidence_value})
    print "Query variable: {}\nEvidence: {}={}\nq:\n{}\n".format(query_variable, evidence_variable, evidence_value, query[query_variable])
'''

# SPECIFIC QUERY
query_variable = '21'

if not model.has_node(query_variable):
    print "Error: Query variable {} not present \n".format(query_variable)
    sys.exit()

print "Running inference on {} \n".format(query_variable)
query = inference.query(variables=[query_variable], evidence={'9': 2})
print "Query result\n{}\n".format(query[query_variable])
